library(cfid)

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

v <- function(i) {
  sprintf("V%d", i)
}

edges_chain <- function(n) {
  out <- character()
  for (i in seq_len(n - 1)) {
    out <- c(out, sprintf("%s -> %s", v(i), v(i + 1)))
  }
  out
}

edges_fanin <- function(n) {
  out <- character()
  for (i in 2:(n - 1)) {
    out <- c(out, sprintf("%s -> %s", v(1), v(i)))
    out <- c(out, sprintf("%s -> %s", v(i), v(n)))
  }
  out <- c(out, sprintf("%s -> %s", v(1), v(n)))
  unique(out)
}

edges_chain_bi <- function(n) {
  out <- edges_chain(n)
  for (i in seq_len(n - 1)) {
    out <- c(out, sprintf("%s <-> %s", v(i), v(i + 1)))
  }
  unique(out)
}

edges_dense <- function(n) {
  out <- edges_chain(n)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      if (((i * 37 + j * 17) %% 10) < 3) {
        out <- c(out, sprintf("%s -> %s", v(i), v(j)))
      }
    }
  }
  unique(out)
}

edge_string <- function(family, n) {
  edges <- switch(
    family,
    chain = edges_chain(n),
    fanin = edges_fanin(n),
    chain_bi = edges_chain_bi(n),
    dense = edges_dense(n),
    stop("unknown family")
  )
  paste(edges, collapse = "; ")
}

run_case <- function(edge_str, treat, out) {
  total_t0 <- now_ms()

  setup_t0 <- now_ms()
  g <- dag(edge_str)
  gamma <- conj(cf(var = out, obs = 0L, sub = setNames(0L, treat)))
  delta <- conj(cf(var = treat, obs = 1L))
  setup_ms <- now_ms() - setup_t0

  ident_t0 <- now_ms()
  ok <- TRUE
  res <- tryCatch(
    identifiable(g, gamma, delta),
    error = function(e) {
      ok <<- FALSE
      list(id = FALSE, formula = NULL, error = conditionMessage(e))
    }
  )
  identify_ms <- now_ms() - ident_t0

  total_ms <- now_ms() - total_t0
  list(
    ok = ok,
    setup_ms = setup_ms,
    identify_ms = identify_ms,
    total_ms = total_ms,
    output = res
  )
}

benchmark_case <- function(f, repeats = 7, warmups = 3) {
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

run_family <- function(family, n, repeats = 7, warmups = 3) {
  treat <- v(1)
  out <- v(n)
  es <- edge_string(family, n)
  stats <- benchmark_case(function() run_case(es, treat, out), repeats = repeats, warmups = warmups)
  w <- stats$warm
  cat(sprintf(
    "impl=r_struct family=%s n=%d ok=%s warm_total_median_ms=%.3f warm_setup_median_ms=%.3f warm_identify_median_ms=%.3f\n",
    family, n, as.character(stats$result$ok), w$total_ms$median_ms, w$setup_ms$median_ms, w$identify_ms$median_ms
  ))
}

args <- commandArgs(trailingOnly = TRUE)
repeats <- if (length(args) >= 1) as.integer(args[[1]]) else 7
warmups <- if (length(args) >= 2) as.integer(args[[2]]) else 3
ns <- if (length(args) >= 3) as.integer(strsplit(args[[3]], ",", fixed = TRUE)[[1]]) else c(8, 16, 24, 32)

cat(sprintf("impl=r_struct benchmark_config repeats=%d warmups=%d timer=Sys.time\n", repeats, warmups))
for (family in c("chain", "fanin", "dense", "chain_bi")) {
  for (n in ns) {
    run_family(family, n, repeats = repeats, warmups = warmups)
  }
}
