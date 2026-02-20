# ------------------------------------------------------------
# party_example_cfid.R
# Counterfactual identification for the "party example"
# using the R package cfid
# ------------------------------------------------------------

library(cfid)

# ------------------------------------------------------------
# A -> B, A -> C, B -> S, C -> S
# ------------------------------------------------------------
g <- dag("A -> B; A -> C; B -> S; C -> S")

# ------------------------------------------------------------
# (These are just labels, not numeric values)
# ------------------------------------------------------------
b  <- 0L   # do(B=b)
bp <- 1L   # observed B=b'
s  <- 0L   # S_b = s

# ------------------------------------------------------------
#   P( S_b = s | B = b' )
# ------------------------------------------------------------

# gamma: S under intervention do(B=b)
gamma <- conj(
  cf(var = "S", obs = s, sub = c(B = b))
)

# delta: observed B = b'
delta <- conj(
  cf(var = "B", obs = bp)
)

# ------------------------------------------------------------
# 4. Identify
# ------------------------------------------------------------
res <- identifiable(g, gamma, delta, data = "interventions")

cat("Identifiable? ", res$id, "\n")

if (isTRUE(res$id)) {
  cat("\nIdentified formula:\n")
  cat(format(res$formula, use_do = TRUE), "\n")
}