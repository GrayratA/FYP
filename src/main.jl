

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
# include("chyp_export.jl") # uncomment if you want trace_dir output



@present CausalR(FreeSymmetricMonoidalCategory) begin
  # objects
  (R, X, W, Z, Y, D)::Ob

  # mechanisms (the c* boxes in your figure)
  c_R::Hom(munit(), R)            # I -> R
  c_X::Hom(R, X)                  # R -> X
  c_W::Hom(X, W)                  # X -> W
  c_Y::Hom(R ⊗ W ⊗ Z, Y)          # (R,W,Z) -> Y

  c_D::Hom(munit(), D)            # I -> D
  c_Z::Hom(D, Z)                  # D -> Z

  # interventions (do boxes)
  doX::Hom(munit(), X)            # I -> X      (label: do_X=doX)
  doD_d::Hom(munit(), D)          # I -> D      (label: do_D=d)
  do_Z::Hom(munit(), Z)           # I -> Z      (label: do_Z)

  # sharp effects (observations)
  obsX_x::Hom(X, munit())         # X -> I      (label: obs_X=x)
  obsD_d::Hom(D, munit())         # D -> I      (label: obs_D=d)
  obsZ_z::Hom(Z, munit())         # Z -> I      (label: obs_Z=z)
end

# 2) Program = the overall wiring diagram in your picture
prog_overall = @program CausalR () begin
  # --- bottom causal core ---
  r  = c_R()

  x_obs = c_X(r)        # this X goes to obs_X=x
  obsX_x(x_obs)

  x_do  = doX()         # do_X=doX, this X goes to c_W
  w     = c_W(x_do)

  z_do  = do_Z()        # do_Z goes into c_Y
  y     = c_Y(r, w, z_do)




  # output (figure only outputs Y)
  return y
end

draw(prog_overall)
draw(add_junctions(prog_overall))

# Unified pipeline example: identify_counterfactual
commute_admg = ADMGModel(
    [
        :W => :T,
        :W => :M,
        :T => :L,
        :M => :L,
        :L => :Y,
    ],
    [
        :T => :Y,
    ],
)

q1 = CounterfactualQuery(
    :World1,
    Dict(:T => :t),
    Dict{Symbol, Symbol}(),
    [:Y],
)

q2 = CounterfactualQuery(
    :World2,
    Dict{Symbol, Symbol}(),
    Dict(:W => :w, :T => :t_),
    Symbol[],
)

q3 = CounterfactualQuery(
    :World3,
    Dict(:W => :w),
    Dict(:M => :m),
    Symbol[],
)

queries = [q1, q2, q3]

res = identify_counterfactual(
    commute_admg,
    queries;
    display_syms=[:W, :T, :M, :L, :Y],
    output_vars=["Y"],
    run_simplify=true,
    step4_verbose=true,
    rules=Step5RuleConfig(
        enabled=true,
        directed_edges=[("W", "T"), ("W", "M"), ("T", "L"), ("M", "L"), ("L", "Y")],
    ),
    display=Step5DisplayConfig(
        symbols=Dict("W" => "w", "T" => "t", "M" => "m", "L" => "l", "Y" => "y"),
        value_rename=Dict("t_" => "t'"),
    ),
    data=Step5DataConfig(
        mode=:interventions,
        mix_var="W",
        mix_sym="w^*",
        anchor_var="W",
        anchor_token="w",
    ),
    # trace_dir="trace/commute",
    # trace_prefix="commute",
)

println("identifiable: ", res.identifiable)
println("failure_stage: ", res.failure_stage)
println("error: ", res.error)
println("raw formula:")
println(res.raw_tex)
println("data formula:")
println(res.data_tex)

draw(add_junctions(res.wd))
