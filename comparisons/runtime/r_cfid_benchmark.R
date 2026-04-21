library(cfid)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  length(file_arg) > 0 || stop("Cannot determine script path")
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
}

now_ms <- function() {
  as.numeric(Sys.time()) * 1000
}

stage_summary <- function(x) {
  list(
    min_ms = min(x),
    median_ms = as.numeric(stats::median(x)),
    mean_ms = mean(x),
    max_ms = max(x)
  )
}

benchmark_stages <- function(f, repeats = 30, warmups = 5) {
  for (i in seq_len(warmups)) {
    f()
  }

  cold <- f()
  samples_setup <- numeric(repeats)
  samples_ident <- numeric(repeats)
  samples_total <- numeric(repeats)
  last_res <- cold

  for (i in seq_len(repeats)) {
    cur <- f()
    samples_setup[[i]] <- cur$setup_ms
    samples_ident[[i]] <- cur$identify_ms
    samples_total[[i]] <- cur$total_ms
    last_res <- cur
  }

  list(
    cold = cold,
    warm = list(
      setup_ms = stage_summary(samples_setup),
      identify_ms = stage_summary(samples_ident),
      total_ms = stage_summary(samples_total)
    ),
    result = last_res
  )
}

parse_hepar2_dag <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- gsub("//.*$", "", lines)
  lines <- gsub("#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  edges <- character()
  for (line in lines) {
    if (!startsWith(line, "potential")) {
      next
    }

    m <- regexec("^potential\\s*\\(\\s*([^|)]+?)\\s*(?:\\|\\s*([^)]+?)\\s*)?\\)$", line)
    caps <- regmatches(line, m)[[1]]
    if (length(caps) == 0) {
      next
    }

    child <- trimws(caps[2])
    parents_raw <- if (length(caps) >= 3) trimws(caps[3]) else ""
    if (!nzchar(parents_raw)) {
      next
    }

    parents <- strsplit(parents_raw, "\\s+")[[1]]
    parents <- parents[nzchar(parents)]
    for (parent in parents) {
      edges <- c(edges, sprintf("%s -> %s", parent, child))
    }
  }

  dag(paste(edges, collapse = "; "))
}

val <- function(x) {
  structure(0L, names = x)
}

run_drug_example <- function() {
  total_t0 <- now_ms()

  setup_t0 <- now_ms()
  g <- dag("X -> W -> Y <- Z <- D; X <-> Y")
  x <- 0L
  xt <- 1L
  d <- 0L
  z <- 0L
  y <- 0L

  gamma <- conj(cf(var = "Y", obs = y, sub = c(X = x)))
  delta <- conj(
    cf(var = "X", obs = xt),
    cf(var = "D", obs = d),
    cf(var = "Z", obs = z, sub = c(D = d))
  )
  setup_ms <- now_ms() - setup_t0

  ident_t0 <- now_ms()
  res <- identifiable(g, gamma, delta, data = "interventions")
  identify_ms <- now_ms() - ident_t0

  total_ms <- now_ms() - total_t0
  list(
    setup_ms = setup_ms,
    identify_ms = identify_ms,
    total_ms = total_ms,
    output = res
  )
}

run_party_example <- function() {
  total_t0 <- now_ms()

  setup_t0 <- now_ms()
  g <- dag("A -> B; A -> C; B -> S; C -> S")
  b <- 0L
  bp <- 1L
  s <- 0L

  gamma <- conj(cf(var = "S", obs = s, sub = c(B = b)))
  delta <- conj(cf(var = "B", obs = bp))
  setup_ms <- now_ms() - setup_t0

  ident_t0 <- now_ms()
  res <- identifiable(g, gamma, delta)
  identify_ms <- now_ms() - ident_t0

  total_ms <- now_ms() - total_t0
  list(
    setup_ms = setup_ms,
    identify_ms = identify_ms,
    total_ms = total_ms,
    output = res
  )
}

run_hepar2_conditional_example <- function() {
  total_t0 <- now_ms()

  setup_t0 <- now_ms()
  g <- parse_hepar2_dag(file.path(project_root, "net", "hepar2.net"))

  gamma <- conj(
    cf(var = "carcinoma", obs = val("carc_obs"), sub = c(age = 0L))
  )

  delta <- conj(
    cf(var = "age", obs = val("age_obs")),
    cf(var = "sex", obs = val("sex_obs")),
    cf(var = "pbc", obs = val("pbc_obs"), sub = c(age = 0L))
  )
  setup_ms <- now_ms() - setup_t0

  ident_t0 <- now_ms()
  res <- identifiable(g, gamma, delta, data = "interventions")
  identify_ms <- now_ms() - ident_t0

  total_ms <- now_ms() - total_t0
  list(
    setup_ms = setup_ms,
    identify_ms = identify_ms,
    total_ms = total_ms,
    output = res
  )
}

print_summary <- function(label, stats) {
  cold <- stats$cold
  warm <- stats$warm
  cat(
    sprintf(
      "impl=r_cfid example=%s cold_total_ms=%.3f warm_total_min_ms=%.3f warm_total_median_ms=%.3f warm_total_mean_ms=%.3f warm_total_max_ms=%.3f warm_setup_median_ms=%.3f warm_identifiable_median_ms=%.3f\n",
      label,
      cold$total_ms,
      warm$total_ms$min_ms, warm$total_ms$median_ms, warm$total_ms$mean_ms, warm$total_ms$max_ms,
      warm$setup_ms$median_ms, warm$identify_ms$median_ms
    )
  )
  cat(sprintf("impl=r_cfid alias example=%s warm_identify_median_ms=%.3f\n", label, warm$identify_ms$median_ms))
  if (isTRUE(stats$result$output$id)) {
    cat(sprintf("formula=%s\n", format(stats$result$output$formula, use_do = TRUE)))
  } else {
    cat("formula=FAIL\n")
  }
}

args <- commandArgs(trailingOnly = TRUE)
repeats <- if (length(args) == 0) 30 else as.integer(args[[1]])
warmups <- if (length(args) <= 1) 5 else as.integer(args[[2]])
project_root <- normalizePath(file.path(dirname(script_path()), "..", ".."), winslash = "/", mustWork = TRUE)

cat(sprintf("impl=r_cfid benchmark_config repeats=%d warmups=%d timer=Sys.time\n", repeats, warmups))
print_summary("drug", benchmark_stages(run_drug_example, repeats = repeats, warmups = warmups))
print_summary("party", benchmark_stages(run_party_example, repeats = repeats, warmups = warmups))
print_summary("hepar2_conditional", benchmark_stages(run_hepar2_conditional_example, repeats = repeats, warmups = warmups))
