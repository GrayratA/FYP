library(cfid)

# Commute ADMG:
# W -> T, W -> M, T -> L, M -> L, L -> Y, and T <-> Y
g <- dag("W -> T; W -> M; T -> L; M -> L; L -> Y; T <-> Y")

# Labels used in the Julia notebook block:
# q1: Y under do(T=t)
# q2: observed W=w, T=t_
# q3: observed M=m under do(W=w)
t_do  <- 0L
y_obs <- 0L
w_obs <- 0L
t_obs <- 1L
m_obs <- 0L

gamma <- conj(
  cf(var = "Y", obs = y_obs, sub = c(T = t_do))
)

delta <- conj(
  cf(var = "W", obs = w_obs),
  cf(var = "T", obs = t_obs),
  cf(var = "M", obs = m_obs, sub = c(W = w_obs))
)

res <- identifiable(g, gamma, delta, data = "interventions")

cat("Identifiable? ", res$id, "\n", sep = "")
if (isTRUE(res$id)) {
  cat("Formula:\n")
  cat(format(res$formula, use_do = TRUE), "\n")
}
