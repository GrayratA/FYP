library(cfid)

now_ms <- function() {
  as.numeric(Sys.time()) * 1000
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

run_drug <- function() {
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
  list(setup_ms = setup_ms, identify_ms = identify_ms, res = res)
}

run_party <- function() {
  setup_t0 <- now_ms()
  g <- dag("A -> B; A -> C; B -> S; C -> S")
  b <- 0L
  bp <- 1L
  s <- 0L
  gamma <- conj(cf(var = "S", obs = s, sub = c(B = b)))
  delta <- conj(cf(var = "B", obs = bp))
  setup_ms <- now_ms() - setup_t0

  ident_t0 <- now_ms()
  res <- identifiable(g, gamma, delta, data = "interventions")
  identify_ms <- now_ms() - ident_t0
  list(setup_ms = setup_ms, identify_ms = identify_ms, res = res)
}

run_hepar2 <- function(project_root) {
  setup_t0 <- now_ms()
  g <- parse_hepar2_dag(file.path(project_root, "net", "hepar2.net"))
  gamma <- conj(cf(var = "carcinoma", obs = val("carc_obs"), sub = c(age = 0L)))
  delta <- conj(
    cf(var = "age", obs = val("age_obs")),
    cf(var = "sex", obs = val("sex_obs")),
    cf(var = "pbc", obs = val("pbc_obs"), sub = c(age = 0L))
  )
  setup_ms <- now_ms() - setup_t0

  ident_t0 <- now_ms()
  res <- identifiable(g, gamma, delta, data = "interventions")
  identify_ms <- now_ms() - ident_t0
  list(setup_ms = setup_ms, identify_ms = identify_ms, res = res)
}

top_stats <- function(summary_obj, fname) {
  tab <- summary_obj[["by.total"]]
  if (is.null(tab)) {
    return(c(total_time = 0, total_pct = 0))
  }
  rows <- rownames(tab)
  key <- gsub("^\"|\"$", "", rows)
  idx <- which(key == fname)
  if (length(idx) == 0) {
    return(c(total_time = 0, total_pct = 0))
  }
  i <- idx[[1]]
  c(total_time = tab[i, "total.time"], total_pct = tab[i, "total.pct"])
}

profile_once <- function(label, runner) {
  runner()
  runner()

  f <- tempfile(pattern = paste0("rprof_", label, "_"), fileext = ".out")
  Rprof(f, interval = 0.005)
  out <- runner()
  Rprof(NULL)
  s <- summaryRprof(f)

  keys <- c("idc_star", "id_star", "make_cg", "cg", "ancestors", "dsep", "%*%")
  kv <- lapply(keys, function(k) top_stats(s, k))
  names(kv) <- keys

  cat(sprintf("example=%s setup_ms=%.3f identify_ms=%.3f\n", label, out$setup_ms, out$identify_ms))
  for (k in keys) {
    cat(sprintf("module=%s total_time=%.3f total_pct=%.2f\n", k, kv[[k]][["total_time"]], kv[[k]][["total_pct"]]))
  }
}

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)

profile_once("drug", run_drug)
profile_once("party", run_party)
profile_once("hepar2_conditional", function() run_hepar2(project_root))
