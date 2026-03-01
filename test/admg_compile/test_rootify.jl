using Test
using Catlab
using Catlab.WiringDiagrams
const WD = Catlab.WiringDiagrams

include(joinpath(@__DIR__, "..", "..", "src", "admg_compile.jl"))

function box_id_by_name(wd::WiringDiagram, name::Symbol)::Int
    for b in WD.box_ids(wd)
        b == WD.input_id(wd) && continue
        b == WD.output_id(wd) && continue
        if WD.box(wd, b).value == name
            return b
        end
    end
    error("Box $name not found")
end

@testset "rootify pairwise bidirected edges" begin
    admg = ADMGModel(
        [
            :X => :Y,
            :X => :Z,
        ],
        [
            :Y => :W1,
            :W1 => :W2,
            :W2 => :Z,
        ],
    )

    expected = ConfoundedModel(
        [
            :X => :Y,
            :X => :Z,
        ],
        Dict(
            :R1 => [:Y, :W1],
            :R2 => [:W1, :W2],
            :R3 => [:W2, :Z],
        ),
    )

    @test rootify(admg) == expected
end

@testset "rootify supports custom latent prefix" begin
    admg = ADMGModel(
        [
            :A => :B,
        ],
        [
            :A => :C,
            :B => :C,
        ],
    )

    @test rootify(admg; latent_prefix=:U) == ConfoundedModel(
        [
            :A => :B,
        ],
        Dict(
            :U1 => [:A, :C],
            :U2 => [:B, :C],
        ),
    )
end

@testset "graph_b_to_scm accepts non-rootified ADMG input" begin
    admg = ADMGModel(
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

    rooted = ConfoundedModel(
        [
            :X => :W,
            :W => :Y,
            :D => :Z,
            :Z => :Y,
        ],
        Dict(
            :R1 => [:X, :Y],
        ),
    )

    display_var = Set([:D, :Z, :X, :W, :Y])
    got = graph_b_to_scm(admg; outputs=display_var)
    expected = graph_b_to_scm(rooted; outputs=display_var)

    @test WD.output_ports(got) == WD.output_ports(expected)
    @test length(WD.wires(got)) == length(WD.wires(expected))

    for name in (:f_D, :PU_D, :f_W, :PU_W, :f_X, :PU_X, :f_Y, :PU_Y, :f_Z, :PU_Z, :f_R1, :PU_R1)
        bid_got = box_id_by_name(got, name)
        bid_expected = box_id_by_name(expected, name)
        @test WD.input_ports(got, bid_got) == WD.input_ports(expected, bid_expected)
        @test WD.output_ports(got, bid_got) == WD.output_ports(expected, bid_expected)
    end
end
