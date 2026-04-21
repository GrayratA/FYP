using Statistics

using Catlab
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Theories
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics: Graphviz, to_graphviz

include(joinpath(@__DIR__, "..", "src", "utils.jl"))
include(joinpath(@__DIR__, "..", "src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "src", "simplify_cf.jl"))
include(joinpath(@__DIR__, "..", "src", "id_cf.jl"))

const STAGE_KEYS = (:setup_ms, :identify_ms, :total_ms, :build_ms, :simplify_ms, :step4_ms, :step5_ms)

stage_summary(xs::Vector{Float64}) = (
    min_ms = minimum(xs),
    median_ms = median(xs),
    mean_ms = mean(xs),
    max_ms = maximum(xs),
)

v(i::Int) = Symbol("V", i)

function model_chain(n::Int)
    directed = Pair{Symbol, Symbol}[v(i) => v(i + 1) for i in 1:(n - 1)]
    return ADMGModel(directed, Pair{Symbol, Symbol}[])
end

function model_fanin(n::Int)
    directed = Pair{Symbol, Symbol}[]
    for i in 2:(n - 1)
        push!(directed, v(1) => v(i))
        push!(directed, v(i) => v(n))
    end
    push!(directed, v(1) => v(n))
    return ADMGModel(unique(directed), Pair{Symbol, Symbol}[])
end

function model_chain_bi(n::Int)
    directed = Pair{Symbol, Symbol}[v(i) => v(i + 1) for i in 1:(n - 1)]
    bidirected = Pair{Symbol, Symbol}[v(i) => v(i + 1) for i in 1:(n - 1)]
    return ADMGModel(directed, bidirected)
end

function model_dense(n::Int)
    directed = Pair{Symbol, Symbol}[v(i) => v(i + 1) for i in 1:(n - 1)]
    for i in 1:(n - 1)
        for j in (i + 1):n
            ((i * 37 + j * 17) % 10 < 3) || continue
            push!(directed, v(i) => v(j))
        end
    end
    return ADMGModel(unique(directed), Pair{Symbol, Symbol}[])
end

function run_case(model::ADMGModel, treat::Symbol, out::Symbol)
    total_t0 = time_ns()

    build_t0 = time_ns()
    display_var = Set([treat, out])
    base_scm = graph_b_to_scm(model; outputs=display_var)

    q_real = CounterfactualQuery(
        :Real,
        Dict{Symbol, Symbol}(),
        Dict(treat => :tp),
        Symbol[],
    )
    q_cf = CounterfactualQuery(
        :CF,
        Dict(treat => :t),
        Dict{Symbol, Symbol}(),
        [out],
    )

    wd = build_multiverse(base_scm, [q_real, q_cf])
    build_ms = (time_ns() - build_t0) / 1e6

    simplify_t0 = time_ns()
    drop_discard_branches!(wd)
    separate_sharp_effect_copy_maps!(wd)
    merge_identical_deterministic_boxes(wd)
    replace_fx_with_cx!(wd)
    simplify_ms = (time_ns() - simplify_t0) / 1e6

    ok = true
    result = nothing
    step4_ms = 0.0
    step5_ms = 0.0
    try
        step4_t0 = time_ns()
        id_cf_step4!(wd, display_var; verbose=false)
        step4_ms = (time_ns() - step4_t0) / 1e6

        step5_t0 = time_ns()
        result = id_cf_step5(
            wd;
            output_vars=[string(out)],
            simplify=true,
            display=Step5DisplayConfig(value_rename=Dict("tp" => "t'")),
            data=Step5DataConfig(mode=:none, anchor_var=string(treat)),
        )
        step5_ms = (time_ns() - step5_t0) / 1e6
    catch
        ok = false
    end

    setup_ms = build_ms + simplify_ms
    identify_ms = step4_ms + step5_ms
    total_ms = (time_ns() - total_t0) / 1e6

    return (
        ok = ok,
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

function benchmark_case(f::Function; repeats::Int, warmups::Int)
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
        last = last,
    )
end

function run_family(family::String, n::Int; repeats::Int, warmups::Int)
    treat = v(1)
    out = v(n)
    model =
        family == "chain" ? model_chain(n) :
        family == "fanin" ? model_fanin(n) :
        family == "dense" ? model_dense(n) :
        family == "chain_bi" ? model_chain_bi(n) :
        error("unknown family")

    stats = benchmark_case(() -> run_case(model, treat, out); repeats=repeats, warmups=warmups)
    w = stats.warm
    println(
        "impl=julia_struct family=$(family) n=$(n) " *
        "ok=$(stats.last.ok) " *
        "warm_total_median_ms=$(round(w[:total_ms].median_ms, digits=3)) " *
        "warm_setup_median_ms=$(round(w[:setup_ms].median_ms, digits=3)) " *
        "warm_build_median_ms=$(round(w[:build_ms].median_ms, digits=3)) " *
        "warm_simplify_median_ms=$(round(w[:simplify_ms].median_ms, digits=3)) " *
        "warm_identify_median_ms=$(round(w[:identify_ms].median_ms, digits=3)) " *
        "warm_step4_median_ms=$(round(w[:step4_ms].median_ms, digits=3)) " *
        "warm_step5_median_ms=$(round(w[:step5_ms].median_ms, digits=3))"
    )
end

function main()
    repeats = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 7
    warmups = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
    ns = length(ARGS) >= 3 ? parse.(Int, split(ARGS[3], ",")) : [8, 16, 24, 32]

    println("impl=julia_struct benchmark_config repeats=$(repeats) warmups=$(warmups) timer=time_ns")
    for family in ("chain", "fanin", "dense", "chain_bi")
        for n in ns
            run_family(family, n; repeats=repeats, warmups=warmups)
        end
    end
end

main()
