using Test

include(joinpath(@__DIR__, "..", "..", "src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "bn_import.jl"))

@testset "read BIF structure into ADMGModel" begin
    bif_text = """
    network sample { }

    variable A {
      type discrete [ 2 ] { yes, no };
    }
    variable B {
      type discrete [ 2 ] { yes, no };
    }
    variable C {
      type discrete [ 2 ] { yes, no };
    }

    probability ( A ) {
      table 0.5, 0.5;
    }

    probability ( B | A ) {
      table 0.1, 0.9, 0.8, 0.2;
    }

    probability ( C | A, B ) {
      table 0.1, 0.9, 0.2, 0.8, 0.3, 0.7, 0.4, 0.6;
    }
    """

    temp = tempname() * ".bif"
    write(temp, bif_text)

    structure = read_bn_structure(temp)

    @test structure.nodes == [:A, :B, :C]
    @test structure.model == ADMGModel(
        Pair{Symbol, Symbol}[
            :A => :B,
            :A => :C,
            :B => :C,
        ],
        Pair{Symbol, Symbol}[],
    )
end

@testset "read Hugin net structure into ADMGModel" begin
    net_text = """
    net { }

    node Smoking {
      states = ("yes" "no");
    }

    node "Lung Cancer" {
      states = ("yes" "no");
    }

    node XRay {
      states = ("yes" "no");
    }

    potential ( Smoking ) {
      data = (0.5 0.5);
    }

    potential ( "Lung Cancer" | Smoking ) {
      data = ((0.1 0.9) (0.8 0.2));
    }

    potential ( XRay | "Lung Cancer" Smoking ) {
      data = (((0.1 0.9) (0.2 0.8)) ((0.3 0.7) (0.4 0.6)));
    }
    """

    temp = tempname() * ".net"
    write(temp, net_text)

    structure = read_bn_structure(temp)

    @test structure.nodes == [:Smoking, :Lung_Cancer, :XRay]
    @test structure.model == ADMGModel(
        Pair{Symbol, Symbol}[
            :Smoking => :Lung_Cancer,
            :Lung_Cancer => :XRay,
            :Smoking => :XRay,
        ],
        Pair{Symbol, Symbol}[],
    )
end

@testset "emit Julia literal for imported model" begin
    structure = BNStructure(
        ADMGModel(
            Pair{Symbol, Symbol}[
                :A => :B,
                :B => :C,
            ],
            Pair{Symbol, Symbol}[],
        ),
        [:A, :B, :C],
    )

    out = sprint(io -> write_admg_model(io, structure; model_name="hepar2_model"))

    @test occursin("hepar2_model = ADMGModel(", out)
    @test occursin("nodes = [:A, :B, :C]", out)
    @test occursin(":A => :B", out)
    @test occursin(":B => :C", out)
end
