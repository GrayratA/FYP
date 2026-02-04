using Test
using Catlab
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Theories
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics: Graphviz, to_graphviz

include(joinpath(@__DIR__, "..", "..","src", "utils.jl"))
include(joinpath(@__DIR__, "..", "..","src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "..","src", "simplify_cf.jl"))
include(joinpath(@__DIR__, "..", "..","src", "id_cf.jl"))

normalize_frags(frags::Vector{<:AbstractVector{<:Integer}}) =
    sort([sort(collect(f)) for f in frags], by = x -> (length(x), x))

@testset "R-fragments ids" begin
    model_2 = ConfoundedModel(
        [
            :X => :Y,
            :X => :Z,
        ],
        Dict(
            :R1 => [:Y,:W1],
            :R2 => [:W1, :W2],
            :R3 => [:W2, :Z],
        )
    )

    display_var = Set([:X, :Y, :Z, :W1, :W2])
    base_scm = graph_b_to_scm(model_2; outputs=display_var)

    q1 = CounterfactualQuery(
        :World1,
        Dict(:X => :doX),
        Dict{Symbol, Symbol}(),
        [:Y, :Z, :W1, :W2]
    )

    queries = [q1]

    full_wd = build_multiverse(base_scm, queries)
    full_wd = replace_fx_with_cx!(full_wd)
    full_wd = drop_discard_branches!(full_wd)

    got = find_r_fragments(full_wd)
    expected = [[1, 4, 5, 2, 3, 8, 7]]

    @test normalize_frags(got) == normalize_frags(expected)
end
