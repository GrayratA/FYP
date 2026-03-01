using Test
using Catlab
using Catlab.WiringDiagrams
const WD = Catlab.WiringDiagrams

include(joinpath(@__DIR__, "..", "..", "src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "simplify_cf.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "id_cf.jl"))

function build_step4_drug_wd()
    model = ConfoundedModel(
        [
            :X => :W,
            :W => :Y,
            :D => :Z,
            :Z => :Y,
        ],
        Dict(
            :R => [:X, :Y],
        ),
    )

    display_var = Set([:D, :Z, :X, :W, :Y])
    base_scm = graph_b_to_scm(model; outputs=display_var)

    q1 = CounterfactualQuery(
        :World1,
        Dict(:X => :doX_),
        Dict{Symbol, Symbol}(),
        [:Y],
    )

    q2 = CounterfactualQuery(
        :World2,
        Dict(),
        Dict(:X => :x_hat, :D => :d),
        Symbol[],
    )

    q3 = CounterfactualQuery(
        :World3,
        Dict(:D => :d),
        Dict(:Z => :z),
        Symbol[],
    )

    wd = build_multiverse(base_scm, [q1, q2, q3])
    drop_discard_branches!(wd)
    merge_identical_deterministic_boxes(wd)
    replace_fx_with_cx!(wd)
    id_cf_step4!(wd, display_var; verbose=false)
    return wd
end

function build_step4_party_wd()
    model = ConfoundedModel(
        [
            :A => :B,
            :A => :C,
            :B => :S,
            :C => :S,
        ],
        Dict(),
    )

    display_var = Set([:B, :S])
    base_scm = graph_b_to_scm(model; outputs=display_var)

    q_real = CounterfactualQuery(
        :Real,
        Dict(),
        Dict(:B => :bp),
        Symbol[],
    )

    q_cf = CounterfactualQuery(
        :CF,
        Dict(:B => :b),
        Dict{Symbol, Symbol}(),
        [:S],
    )

    wd = build_multiverse(base_scm, [q_real, q_cf])
    drop_discard_branches!(wd)
    merge_identical_deterministic_boxes(wd)
    replace_fx_with_cx!(wd)
    id_cf_step4!(wd, display_var; verbose=false)
    return wd
end

@testset "id-cf Step5: drug example" begin
    wd = build_step4_drug_wd()

    res = id_cf_step5(
        wd;
        output_vars=["Y"],
        simplify=true,
        enable_cfid_rules=true,
        directed_edges=[("X", "W"), ("W", "Y"), ("Z", "Y"), ("D", "Z")],
        var2sym=Dict("Y" => "y", "W" => "w", "D" => "d", "Z" => "z", "X" => "x"),
        obs_value_rename=Dict("x_hat" => "x'"),
        data_mode=:interventions,
        mix_var="D",
        mix_sym="d^*",
        mix_target_var="Z",
        anchor_var="X",
    )

    @test res.raw_tex ==
        "\\frac{\\sum_{w} P(z|do(d))P(x',y|do(z,w))P(d)P(w|do(doX_))}{P(z|do(d))P(x')P(d)}"
    @test res.simplified_tex == res.raw_tex
    @test res.data_tex ==
        "\\frac{\\sum_{w,d^*} P(z|do(d))P(x',y|do(z,w))P(d)P(w|do(doX_))P(d^*)}{\\sum_{d^*} P(x')P(z|do(d^*))P(d^*)}"
end

@testset "id-cf Step5: party example" begin
    wd = build_step4_party_wd()

    res = id_cf_step5(
        wd;
        output_vars=["S"],
        simplify=true,
        enable_cfid_rules=true,
        directed_edges=[("A", "B"), ("A", "C"), ("B", "S"), ("C", "S")],
        var2sym=Dict("A" => "a", "B" => "b", "C" => "c", "S" => "s"),
        obs_value_rename=Dict("bp" => "b'"),
        anchor_var="B",
        data_mode=:none,
    )

    @test res.raw_tex ==
        "\\frac{\\sum_{a,c} P(s|do(b,c))P(c|do(a))P(a)P(b'|do(a))}{\\sum_{a} P(a)P(b'|do(a))}"
    @test res.simplified_tex == res.raw_tex
    @test res.data_tex == res.raw_tex
end
