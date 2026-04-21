library(cfid)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  length(file_arg) > 0 || stop("Cannot determine script path")
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
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

project_root <- normalizePath(file.path(dirname(script_path()), "..", ".."), winslash = "/", mustWork = TRUE)
hepar2_path <- file.path(project_root, "net", "hepar2.net")
g <- parse_hepar2_dag(hepar2_path)

# Fraction-producing full-HEPAR2 query corresponding to the Julia output version:
#   P(carcinoma_{do(age)} = carc_obs | age = age_obs, sex = sex_obs, pbc_{do(age)} = pbc_obs)
gamma <- conj(
  cf(var = "carcinoma", obs = val("carc_obs"), sub = c(age = 0L))
)

delta <- conj(
  cf(var = "age", obs = val("age_obs")),
  cf(var = "sex", obs = val("sex_obs")),
  cf(var = "pbc", obs = val("pbc_obs"), sub = c(age = 0L))
)

res <- identifiable(g, gamma, delta, data = "interventions")

cat("Identifiable? ", res$id, "\n", sep = "")

if (isTRUE(res$id)) {
  cat("Formula:\n")
  cat(format(res$formula, use_do = TRUE), "\n")
}
