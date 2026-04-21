using Statistics

using Catlab
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Theories
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics: Graphviz, to_graphviz

include(joinpath(@__DIR__, "..", "..", "src", "utils.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "bn_import.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "simplify_cf.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "id_cf.jl"))

const STAGE_KEYS = (:build_ms, :simplify_ms, :step4_ms, :step5_ms, :setup_ms, :identify_ms, :total_ms)

stage_summary(xs::Vector{Float64}) = (
    min_ms = minimum(xs),
    median_ms = median(xs),
    mean_ms = mean(xs),
    max_ms = maximum(xs),
)

function benchmark_stages(f::Function; repeats::Int=30, warmups::Int=5)
    for _ in 1:warmups
        f()
    end

    cold = f()
    warm = Dict{Symbol, Vector{Float64}}(k => Float64[] for k in STAGE_KEYS)
    last = cold

    for _ in 1:repeats
        cur = f()
        for k in STAGE_KEYS
            push!(warm[k], cur[k])
        end
        last = cur
    end

    return (
        cold = cold,
        warm = Dict(k => stage_summary(warm[k]) for k in STAGE_KEYS),
        result = last.result,
    )
end

function run_drug_example()
    total_t0 = time_ns()

    build_t0 = time_ns()
    model = ADMGModel(
        [
            :X => :W,
            :W => :Y,
            :D => :Z,
            :Z => :Y,
        ],
        [
            :X => :Y,
        ],
    )

    display_var = Set([:D, :Z, :X, :W, :Y])
    base_scm = graph_b_to_scm(model; outputs=display_var)

    q1 = CounterfactualQuery(
        :World1,
        Dict(:X => :x),
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
    build_ms = (time_ns() - build_t0) / 1e6

    simplify_t0 = time_ns()
    drop_discard_branches!(wd)
    separate_sharp_effect_copy_maps!(wd)
    merge_identical_deterministic_boxes(wd)
    replace_fx_with_cx!(wd)
    simplify_ms = (time_ns() - simplify_t0) / 1e6

    step4_t0 = time_ns()
    id_cf_step4!(wd, display_var; verbose=false)
    step4_ms = (time_ns() - step4_t0) / 1e6

    step5_t0 = time_ns()
    result = id_cf_step5(
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
    step5_ms = (time_ns() - step5_t0) / 1e6
    setup_ms = build_ms + simplify_ms
    identify_ms = step4_ms + step5_ms
    total_ms = (time_ns() - total_t0) / 1e6

    return (
        build_ms = build_ms,
        simplify_ms = simplify_ms,
        step4_ms = step4_ms,
        step5_ms = step5_ms,
        setup_ms = setup_ms,
        identify_ms = identify_ms,
        total_ms = total_ms,
        result = result,
    )
end

function run_party_example()
    total_t0 = time_ns()

    build_t0 = time_ns()
    model = ADMGModel(
        [
            :A => :B,
            :A => :C,
            :B => :S,
            :C => :S,
        ],
        Pair{Symbol, Symbol}[],
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
    build_ms = (time_ns() - build_t0) / 1e6

    simplify_t0 = time_ns()
    drop_discard_branches!(wd)
    separate_sharp_effect_copy_maps!(wd)
    merge_identical_deterministic_boxes(wd)
    replace_fx_with_cx!(wd)
    simplify_ms = (time_ns() - simplify_t0) / 1e6

    step4_t0 = time_ns()
    id_cf_step4!(wd, display_var; verbose=false)
    step4_ms = (time_ns() - step4_t0) / 1e6

    step5_t0 = time_ns()
    result = id_cf_step5(
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
    step5_ms = (time_ns() - step5_t0) / 1e6
    setup_ms = build_ms + simplify_ms
    identify_ms = step4_ms + step5_ms
    total_ms = (time_ns() - total_t0) / 1e6

    return (
        build_ms = build_ms,
        simplify_ms = simplify_ms,
        step4_ms = step4_ms,
        step5_ms = step5_ms,
        setup_ms = setup_ms,
        identify_ms = identify_ms,
        total_ms = total_ms,
        result = result,
    )
end

function run_hepar2_conditional_example()
    total_t0 = time_ns()

    build_t0 = time_ns()
    hepar2 = read_bn_structure(joinpath(@__DIR__, "..", "..", "net", "hepar2.net"))
    model = hepar2.model

    display_var = Set([:age, :sex, :pbc, :carcinoma])
    base_scm = graph_b_to_scm(model; outputs=display_var)

    queries = [
        CounterfactualQuery(:Real, Dict{Symbol, Symbol}(), Dict(:age => :age_obs, :sex => :sex_obs), Symbol[]),
        CounterfactualQuery(:CF_pbc, Dict(:age => :age_do), Dict(:pbc => :pbc_obs), Symbol[]),
        CounterfactualQuery(:CF_carc, Dict(:age => :age_do), Dict(:carcinoma => :carc_obs), Symbol[]),
    ]

    wd = build_multiverse(base_scm, queries)
    build_ms = (time_ns() - build_t0) / 1e6

    simplify_t0 = time_ns()
    drop_discard_branches!(wd)
    separate_sharp_effect_copy_maps!(wd)
    merge_identical_deterministic_boxes(wd)
    replace_fx_with_cx!(wd)
    simplify_ms = (time_ns() - simplify_t0) / 1e6

    step4_t0 = time_ns()
    id_cf_step4!(wd, display_var; verbose=false)
    step4_ms = (time_ns() - step4_t0) / 1e6

    step5_t0 = time_ns()
    result = id_cf_step5(
        wd;
        output_vars=["carcinoma"],
        simplify=true,
        queries=queries,
        data=Step5DataConfig(mode=:conditional_queries),
        display=Step5DisplayConfig(
            symbols=Dict("age" => "age", "sex" => "sex", "pbc" => "pbc", "carcinoma" => "carcinoma"),
            value_rename=Dict("age_do" => "age"),
        ),
    )
    step5_ms = (time_ns() - step5_t0) / 1e6
    setup_ms = build_ms + simplify_ms
    identify_ms = step4_ms + step5_ms
    total_ms = (time_ns() - total_t0) / 1e6

    return (
        build_ms = build_ms,
        simplify_ms = simplify_ms,
        step4_ms = step4_ms,
        step5_ms = step5_ms,
        setup_ms = setup_ms,
        identify_ms = identify_ms,
        total_ms = total_ms,
        result = result,
    )
end

function print_summary(label::String, stats)
    c = stats.cold
    w = stats.warm
    println(
        "impl=julia example=$(label) " *
        "cold_total_ms=$(round(c.total_ms, digits=3)) " *
        "warm_total_min_ms=$(round(w[:total_ms].min_ms, digits=3)) " *
        "warm_total_median_ms=$(round(w[:total_ms].median_ms, digits=3)) " *
        "warm_total_mean_ms=$(round(w[:total_ms].mean_ms, digits=3)) " *
        "warm_total_max_ms=$(round(w[:total_ms].max_ms, digits=3)) " *
        "warm_setup_median_ms=$(round(w[:setup_ms].median_ms, digits=3)) " *
        "warm_identify_median_ms=$(round(w[:identify_ms].median_ms, digits=3)) " *
        "warm_build_median_ms=$(round(w[:build_ms].median_ms, digits=3)) " *
        "warm_simplify_median_ms=$(round(w[:simplify_ms].median_ms, digits=3)) " *
        "warm_step4_median_ms=$(round(w[:step4_ms].median_ms, digits=3)) " *
        "warm_step5_median_ms=$(round(w[:step5_ms].median_ms, digits=3))"
    )
    println("formula=$(stats.result.data_tex)")
end

repeats = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 30
warmups = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5

println("impl=julia benchmark_config repeats=$(repeats) warmups=$(warmups) timer=time_ns")
print_summary("drug", benchmark_stages(run_drug_example; repeats=repeats, warmups=warmups))
print_summary("party", benchmark_stages(run_party_example; repeats=repeats, warmups=warmups))
print_summary("hepar2_conditional", benchmark_stages(run_hepar2_conditional_example; repeats=repeats, warmups=warmups))
