library(cfid)

benchmark_ms <- function(f, repeats = 30) {
  cold <- system.time(cold_res <- f())[["elapsed"]] * 1000
  samples <- numeric(repeats)
  last_res <- cold_res

  for (i in seq_len(repeats)) {
    t <- system.time(last_res <- f())[["elapsed"]] * 1000
    samples[[i]] <- t
  }

  list(
    cold_ms = cold,
    min_ms = min(samples),
    median_ms = median(samples),
    mean_ms = mean(samples),
    max_ms = max(samples),
    result = last_res
  )
}

run_drug_example <- function() {
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

  identifiable(g, gamma, delta, data = "interventions")
}

run_party_example <- function() {
  g <- dag("A -> B; A -> C; B -> S; C -> S")
  b <- 0L
  bp <- 1L
  s <- 0L

  gamma <- conj(cf(var = "S", obs = s, sub = c(B = b)))
  delta <- conj(cf(var = "B", obs = bp))

  identifiable(g, gamma, delta, data = "interventions")
}

print_summary <- function(label, stats) {
  cat(
    sprintf(
      "impl=r_cfid example=%s cold_ms=%.3f min_ms=%.3f median_ms=%.3f mean_ms=%.3f max_ms=%.3f\n",
      label, stats$cold_ms, stats$min_ms, stats$median_ms, stats$mean_ms, stats$max_ms
    )
  )
  if (isTRUE(stats$result$id)) {
    cat(sprintf("formula=%s\n", format(stats$result$formula, use_do = TRUE)))
  } else {
    cat("formula=FAIL\n")
  }
}

args <- commandArgs(trailingOnly = TRUE)
repeats <- if (length(args) == 0) 30 else as.integer(args[[1]])

print_summary("drug", benchmark_ms(run_drug_example, repeats = repeats))
print_summary("party", benchmark_ms(run_party_example, repeats = repeats))
