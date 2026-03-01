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
include(joinpath(@__DIR__, "..", "..", "src", "simplify_cf.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "id_cf.jl"))

function benchmark_ms(f::Function; repeats::Int=30)
    cold_t0 = time_ns()
    cold_res = f()
    cold_ms = (time_ns() - cold_t0) / 1e6

    samples = Float64[]
    last_res = cold_res
    for _ in 1:repeats
        t0 = time_ns()
        last_res = f()
        push!(samples, (time_ns() - t0) / 1e6)
    end

    return (
        cold_ms = cold_ms,
        min_ms = minimum(samples),
        median_ms = median(samples),
        mean_ms = mean(samples),
        max_ms = maximum(samples),
        result = last_res,
    )
end

function run_drug_example()
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
    drop_discard_branches!(wd)
    merge_identical_deterministic_boxes(wd)
    replace_fx_with_cx!(wd)
    id_cf_step4!(wd, display_var; verbose=false)

    return id_cf_step5(
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
end

function run_party_example()
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
    drop_discard_branches!(wd)
    merge_identical_deterministic_boxes(wd)
    replace_fx_with_cx!(wd)
    id_cf_step4!(wd, display_var; verbose=false)

    return id_cf_step5(
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
end

function print_summary(label::String, stats)
    println("impl=julia example=$(label) cold_ms=$(round(stats.cold_ms, digits=3)) min_ms=$(round(stats.min_ms, digits=3)) median_ms=$(round(stats.median_ms, digits=3)) mean_ms=$(round(stats.mean_ms, digits=3)) max_ms=$(round(stats.max_ms, digits=3))")
    println("formula=$(stats.result.data_tex)")
end

repeats = isempty(ARGS) ? 30 : parse(Int, ARGS[1])

print_summary("drug", benchmark_ms(run_drug_example; repeats=repeats))
print_summary("party", benchmark_ms(run_party_example; repeats=repeats))
