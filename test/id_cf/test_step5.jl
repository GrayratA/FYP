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
    separate_sharp_effect_copy_maps!(wd)
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
    separate_sharp_effect_copy_maps!(wd)
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
        rules=Step5RuleConfig(
            enabled=true,
            directed_edges=[("X", "W"), ("W", "Y"), ("Z", "Y"), ("D", "Z")],
        ),
        display=Step5DisplayConfig(
            symbols=Dict("Y" => "y", "W" => "w", "D" => "d", "Z" => "z", "X" => "x"),
            value_rename=Dict("x_hat" => "x'"),
        ),
        data=Step5DataConfig(
            mode=:interventions,
            mix_var="D",
            mix_sym="d^*",
            mix_target_var="Z",
            anchor_var="X",
        ),
    )

    @test res.raw_tex ==
        "\\frac{\\sum_{w} P(z|do(d))P(x',y|do(z,w))P(d)P(w|do(doX_))}{P(z|do(d))P(x')P(d)}"
    @test res.simplified_tex == res.raw_tex
    @test res.data_tex ==
        "\\frac{\\sum_{w,d^*} P(z|do(d))P(x',y|do(z,w))P(d^*)P(w|do(doX_))}{\\sum_{d^*} P(x')P(z|do(d^*))P(d^*)}"
end

@testset "id-cf Step5: party example" begin
    wd = build_step4_party_wd()

    res = id_cf_step5(
        wd;
        output_vars=["S"],
        simplify=true,
        rules=Step5RuleConfig(
            enabled=true,
            directed_edges=[("A", "B"), ("A", "C"), ("B", "S"), ("C", "S")],
        ),
        display=Step5DisplayConfig(
            symbols=Dict("A" => "a", "B" => "b", "C" => "c", "S" => "s"),
            value_rename=Dict("bp" => "b'"),
        ),
        data=Step5DataConfig(mode=:none, anchor_var="B"),
    )

    @test res.raw_tex ==
        "\\frac{\\sum_{a,c} P(s|do(b,c))P(c|do(a))P(a)P(b'|do(a))}{\\sum_{a} P(a)P(b'|do(a))}"
    @test res.simplified_tex == res.raw_tex
    @test res.data_tex == res.raw_tex
end

@testset "id-cf Step5: parse grouped do-vars and drop trivial empty-left factors" begin
    f = parse_prob_box_name(Symbol("P_7(A ; do(B,C),do(D))"))
    @test f !== nothing
    @test f.left == ["A"]
    @test f.dovars == ["B", "C", "D"]

    wd = WiringDiagram([], Any[])
    WD.add_box!(wd, Box(Symbol("P_1(∅ ; do(X,Y),do(Z))"), Any[], Any[]))
    WD.add_box!(wd, Box(Symbol("P_2(A ; do(B,C),do(D))"), Any[], Any[:A]))

    res = id_cf_step5(
        wd;
        output_vars=["A"],
        simplify=false,
    )

    @test !occursin("P(|do(", res.raw_tex)
    @test occursin("P(a|do(b,c,d))", res.raw_tex)
end

@testset "id-cf Step5: empty product renders as 1" begin
    @test cf_latex(CFProd(CFExpr[])) == "1"
    @test cf_latex(CFFrac(CFAtom("P(x)", ["X"], String[], ["x"], String[]), CFProd(CFExpr[]))) ==
        "\\frac{P(x)}{1}"
end

@testset "id-cf Step5: empty output_vars returns event kernel, not self-normalized fraction" begin
    wd = WiringDiagram([], Any[])
    pAge = WD.add_box!(wd, Box(Symbol("P_1(age ; ∅)"), Any[], Any[:age]))
    pB = WD.add_box!(wd, Box(Symbol("P_2(B ; do(age))"), Any[:age], Any[:B]))
    obsAge = WD.add_box!(wd, Box(Symbol("obs_age=age_obs"), Any[:age], Any[]))
    obsB = WD.add_box!(wd, Box(Symbol("obs_B=b_obs"), Any[:B], Any[]))

    WD.add_wire!(wd, (pAge, 1) => (pB, 1))
    WD.add_wire!(wd, (pAge, 1) => (obsAge, 1))
    WD.add_wire!(wd, (pB, 1) => (obsB, 1))

    res = id_cf_step5(
        wd;
        output_vars=String[],
        simplify=false,
    )

    @test !occursin("\\frac{", res.raw_tex)
    @test res.raw_tex == "P(age_obs)P(b_obs|do(age))"
    @test res.simplified_tex == res.raw_tex
    @test res.data_tex == res.raw_tex
end

@testset "id-cf Step5: query-fixed observed vars are not re-summed as latent" begin
    wd = WiringDiagram([], Any[])
    pA = WD.add_box!(wd, Box(Symbol("P_1(A ; ∅)"), Any[], Any[:A]))
    pB = WD.add_box!(wd, Box(Symbol("P_2(B ; do(A))"), Any[:A], Any[:B]))
    obsB = WD.add_box!(wd, Box(Symbol("obs_B=b_obs"), Any[:B], Any[]))

    WD.add_wire!(wd, (pA, 1) => (pB, 1))
    WD.add_wire!(wd, (pB, 1) => (obsB, 1))

    res = id_cf_step5(
        wd;
        output_vars=String[],
        simplify=false,
        fixed=Step5FixedConfig(obs=Dict("a" => "a_obs", "b" => "b_obs")),
    )

    @test !occursin("\\sum_{a}", res.raw_tex)
    @test occursin("P(a_obs)", res.raw_tex)
    @test occursin("P(b_obs|do(a))", res.raw_tex)
end

@testset "id-cf Step5: conditional_queries rewrites fixed conditions into starred-style fraction" begin
    wd = WiringDiagram([], Any[])
    pA = WD.add_box!(wd, Box(Symbol("P_1(A ; ∅)"), Any[], Any[:A]))
    pS = WD.add_box!(wd, Box(Symbol("P_2(S ; ∅)"), Any[], Any[:S]))
    pP = WD.add_box!(wd, Box(Symbol("P_3(P ; do(A,S))"), Any[:A, :S], Any[:P]))
    pY = WD.add_box!(wd, Box(Symbol("P_4(Y ; do(P))"), Any[:P], Any[:Y]))

    WD.add_wire!(wd, (pA, 1) => (pP, 1))
    WD.add_wire!(wd, (pS, 1) => (pP, 2))
    WD.add_wire!(wd, (pP, 1) => (pY, 1))

    queries = [
        CounterfactualQuery(:Real, Dict{Symbol, Symbol}(), Dict(:A => :a_obs, :S => :s_obs), Symbol[]),
        CounterfactualQuery(:CF_p, Dict(:A => :a_do), Dict(:P => :p_obs), Symbol[]),
        CounterfactualQuery(:CF_y, Dict(:A => :a_do), Dict(:Y => :y_obs), Symbol[]),
    ]

    res = id_cf_step5(
        wd;
        output_vars=["Y"],
        simplify=true,
        queries=queries,
        data=Step5DataConfig(mode=:conditional_queries),
        display=Step5DisplayConfig(
            symbols=Dict("A" => "a", "S" => "s", "P" => "p", "Y" => "y"),
            value_rename=Dict("a_do" => "a"),
        ),
    )

    @test occursin("\\sum_{s^*}", res.data_tex)
    @test occursin("P(s^*)", res.data_tex)
    @test occursin("P(p|do(a,s))", res.data_tex)
    @test occursin("\\sum_{s^*}", res.data_tex)
    @test occursin("P(p|do(a,s^*))", res.data_tex)
    @test occursin("P(y|do(p))", res.data_tex)
    @test !occursin("a_obs", res.data_tex)
    @test !occursin("s_obs", res.data_tex)
end

@testset "id-cf Step5: display_only keeps fixed vars in summation domain" begin
    wd = WiringDiagram([], Any[])
    WD.add_box!(wd, Box(Symbol("P_1(age ; ∅)"), Any[], Any[:age]))
    WD.add_box!(wd, Box(Symbol("P_2(B ; do(age))"), Any[:age], Any[:B]))
    WD.add_box!(wd, Box(Symbol("obs_age=age_obs"), Any[:age], Any[]))

    legacy = id_cf_step5(
        wd;
        output_vars=["B"],
        simplify=false,
    )
    display_only = id_cf_step5(
        wd;
        output_vars=["B"],
        latent=Step5LatentConfig(mode=:include_fixed),
        simplify=false,
    )

    @test !occursin("\\sum_{age}", legacy.raw_tex)
    @test occursin("\\sum_{age}", display_only.raw_tex)
    @test occursin("P(age_obs)", display_only.raw_tex)
end

@testset "id-cf Step5: starred_style display uses starred fixed roots and symbolic do-vars" begin
    wd = WiringDiagram([], Any[])
    WD.add_box!(wd, Box(Symbol("P_1(age ; ∅)"), Any[], Any[:age]))
    WD.add_box!(wd, Box(Symbol("P_2(B ; do(age,C))"), Any[:age, :C], Any[:B]))
    WD.add_box!(wd, Box(Symbol("obs_age=age_obs"), Any[:age], Any[]))
    WD.add_box!(wd, Box(Symbol("do_C=C_do"), Any[], Any[:C]))

    res = id_cf_step5(
        wd;
        output_vars=["B"],
        latent=Step5LatentConfig(mode=:include_fixed),
        display=Step5DisplayConfig(
            mode=:starred_style,
            symbols=Dict("age" => "a", "B" => "b", "C" => "c"),
        ),
        simplify=false,
    )

    @test occursin("P(a^*)", res.raw_tex)
    @test occursin("P(b|do(a,c))", res.raw_tex)
    @test !occursin("age_obs", res.raw_tex)
    @test !occursin("C_do", res.raw_tex)
end

@testset "id-cf Step5: starred_style synthesizes starred root factor for do-only variable" begin
    wd = WiringDiagram([], Any[])
    WD.add_box!(wd, Box(Symbol("P_1(B ; do(H))"), Any[:H], Any[:B]))
    WD.add_box!(wd, Box(Symbol("do_H=h_do"), Any[], Any[:H]))

    res = id_cf_step5(
        wd;
        output_vars=["B"],
        latent=Step5LatentConfig(mode=:include_fixed),
        display=Step5DisplayConfig(
            mode=:starred_style,
            symbols=Dict("B" => "b", "H" => "h"),
        ),
        simplify=false,
    )

    @test occursin("P(h^*)", res.raw_tex)
    @test occursin("P(b|do(h))", res.raw_tex)
end
