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

runner_drug <- function() {
  g <- dag("X -> W -> Y <- Z <- D; X <-> Y")
  gamma <- conj(cf(var = "Y", obs = 0L, sub = c(X = 0L)))
  delta <- conj(cf(var = "X", obs = 1L), cf(var = "D", obs = 0L), cf(var = "Z", obs = 0L, sub = c(D = 0L)))
  identifiable(g, gamma, delta, data = "interventions")
}

runner_party <- function() {
  g <- dag("A -> B; A -> C; B -> S; C -> S")
  gamma <- conj(cf(var = "S", obs = 0L, sub = c(B = 0L)))
  delta <- conj(cf(var = "B", obs = 1L))
  identifiable(g, gamma, delta, data = "interventions")
}

runner_hepar2 <- function(project_root) {
  g <- parse_hepar2_dag(file.path(project_root, "net", "hepar2.net"))
  gamma <- conj(cf(var = "carcinoma", obs = val("carc_obs"), sub = c(age = 0L)))
  delta <- conj(
    cf(var = "age", obs = val("age_obs")),
    cf(var = "sex", obs = val("sex_obs")),
    cf(var = "pbc", obs = val("pbc_obs"), sub = c(age = 0L))
  )
  identifiable(g, gamma, delta, data = "interventions")
}

extract <- function(tab, name, nloops) {
  if (is.null(tab)) {
    return(c(total_time_s = 0, total_pct = 0, per_call_ms = 0))
  }
  rn <- gsub("^\"|\"$", "", rownames(tab))
  i <- which(rn == name)
  if (length(i) == 0) {
    return(c(total_time_s = 0, total_pct = 0, per_call_ms = 0))
  }
  idx <- i[[1]]
  t <- tab[idx, "total.time"]
  p <- tab[idx, "total.pct"]
  c(total_time_s = t, total_pct = p, per_call_ms = (t * 1000) / nloops)
}

run_profile <- function(label, runner, nloops) {
  runner()
  f <- tempfile(pattern = paste0("rprof_", label, "_"), fileext = ".out")
  Rprof(f, interval = 0.005)
  t0 <- now_ms()
  for (i in seq_len(nloops)) {
    runner()
  }
  elapsed <- now_ms() - t0
  Rprof(NULL)
  s <- summaryRprof(f)
  tab <- s[["by.total"]]
  keys <- c("idc_star", "id_star", "make_cg", "cg", "ancestors", "dsep", "%*%")
  cat(sprintf("example=%s loops=%d avg_identify_ms=%.3f\n", label, nloops, elapsed / nloops))
  for (k in keys) {
    e <- extract(tab, k, nloops)
    cat(sprintf("module=%s total_pct=%.2f per_call_ms=%.3f\n", k, e[["total_pct"]], e[["per_call_ms"]]))
  }
}

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)

run_profile("drug", runner_drug, 300)
run_profile("party", runner_party, 400)
run_profile("hepar2_conditional", function() runner_hepar2(project_root), 3)
