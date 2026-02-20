library(cfid)

# ADMG: X -> W -> Y <- Z <- D and X <-> Y
g <- dag("X -> W -> Y <- Z <- D; X <-> Y")
x  <- 0L          # world1: do(X=x)
xt <- 1L          # world2: observed X = x~
d  <- 0L          # world2 observed D=d, world3 do(D=d)
z  <- 0L          # world3: Z_d = z
y  <- 0L          

# gamma: Y(1)=y
vY <- cf(var = "Y", obs = y, sub = c(X = x))
gamma <- conj(vY)

# delta: X(2)=x~, D(2)=d, Z(3)=z with world3 do(D=d)
vX <- cf(var = "X", obs = xt)             # observed world (no sub)
vD <- cf(var = "D", obs = d)              # observed world (no sub)
vZ <- cf(var = "Z", obs = z, sub = c(D = d))  # world3 do(D=d)
delta <- conj(vX, vD, vZ)

# conditional counterfactual
res <- identifiable(g, gamma, delta, data = "interventions")
print(res$id)
if (isTRUE(res$id)) {
  cat(format(res$formula, use_do = TRUE), "\n")
}