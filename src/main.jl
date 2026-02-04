

using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics: Graphviz, to_graphviz

include("admg_compile.jl")
include("id_cf.jl")
include("simplify_cf.jl")
include("utils.jl")



@present CausalR(FreeSymmetricMonoidalCategory) begin
  # objects
  (R, X, W, Z, Y, D)::Ob

  # mechanisms (the c* boxes in your figure)
  cUR::Hom(munit(), R)            # I -> R
  cUX::Hom(R, X)                  # R -> X
  cUW::Hom(X, W)                  # X -> W
  cUY::Hom(R ⊗ W ⊗ Z, Y)          # (R,W,Z) -> Y

  cUD::Hom(munit(), D)            # I -> D
  cUZ::Hom(D, Z)                  # D -> Z

  # interventions (do boxes)
  doX::Hom(munit(), X)            # I -> X      (label: do_X=doX)
  doD_d::Hom(munit(), D)          # I -> D      (label: do_D=d)
  doZ::Hom(munit(), Z)            # I -> Z      (label: doZ)

  # sharp effects (observations)
  obsX_x::Hom(X, munit())         # X -> I      (label: obs_X=x)
  obsD_d::Hom(D, munit())         # D -> I      (label: obs_D=d)
  obsZ_z::Hom(Z, munit())         # Z -> I      (label: obs_Z=z)
end

# 2) Program = the overall wiring diagram in your picture
prog_overall = @program CausalR () begin
  # --- bottom causal core ---
  r  = cUR()

  x_obs = cUX(r)        # this X goes to obs_X=x
  obsX_x(x_obs)

  x_do  = doX()         # do_X=doX, this X goes to cUW
  w     = cUW(x_do)

  z_do  = doZ()         # doZ goes into cUY
  y     = cUY(r, w, z_do)




  # output (figure only outputs Y)
  return y
end

draw(prog_overall)
draw(add_junctions(prog_overall))
display_var = Set([:Z, :X, :W, :Y])
r_frags = id_cf_step41!(prog_overall, display_var)
# wd = id_cf_step42!(prog_overall, r_frags)