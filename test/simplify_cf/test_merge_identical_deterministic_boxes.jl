using Test
using Catlab
using Catlab.WiringDiagrams
const WD = Catlab.WiringDiagrams

include(joinpath(@__DIR__, "..", "..", "src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "id_cf.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "simplify_cf.jl"))

# ---------- helpers ----------
function box_id_by_name(wd::WiringDiagram, name::Symbol)::Int
    for b in WD.box_ids(wd)
        b == WD.input_id(wd)  && continue
        b == WD.output_id(wd) && continue
        if WD.box(wd, b).value == name
            return b
        end
    end
    error("Box $name not found")
end

function has_wire_pair(wd::WiringDiagram, src::Tuple{Int,Int}, tgt::Tuple{Int,Int})::Bool
    for w in WD.wires(wd)
        if (w.source.box, w.source.port) == src &&
           (w.target.box, w.target.port) == tgt
            return true
        end
    end
    return false
end
# --------------------------------


@testset "multiverse test" begin
    model = ConfoundedModel(
        [
            :X => :W,
            :W => :Y,
            :D => :Z,
            :Z => :Y
        ],
        Dict(:R => [:X, :Y])
    )

    base = graph_b_to_scm(model)

    queries = [
        CounterfactualQuery(:World1, Dict(:X => :x), Dict{Symbol,Symbol}(), [:Y]),
        CounterfactualQuery(:World2, Dict(), Dict(:X => :x0, :D => :d), Symbol[]),
        CounterfactualQuery(:World3, Dict(:D => :d), Dict(:Z => :z), Symbol[]),
    ]

    wd = build_multiverse(base, queries)
    drop_discard_branches!(wd)
    merge_identical_deterministic_boxes(wd)

    # 1) outputs
    @test WD.output_ports(wd) == Any[:Y]

    # 2) expected boxes exist
    expected_boxes = [
        :PU_D, :PU_R, :PU_W, :PU_X, :PU_Y, :PU_Z,
        :f_D, :f_R, :f_W, :f_Z, :f_X, :f_Y,
        Symbol("do_X=x"),
        Symbol("do_D=d"),
        :doZ,
        Symbol("obs_X=x0"),
        Symbol("obs_D=d"),
        Symbol("obs_Z=z"),
    ]

    for nm in expected_boxes
        @test box_id_by_name(wd, nm) isa Int
    end

    # 3) check some key port signatures (spot-check, not everything)
    fY = box_id_by_name(wd, :f_Y)
    @test WD.input_ports(wd, fY)  == Any[:W, :Z, :R, :UY]
    @test WD.output_ports(wd, fY) == Any[:Y]

    fX = box_id_by_name(wd, :f_X)
    @test WD.input_ports(wd, fX) == Any[:UX, :R]

    obsX = box_id_by_name(wd, Symbol("obs_X=x0"))
    @test WD.input_ports(wd, obsX) == Any[:X]

    doX = box_id_by_name(wd, Symbol("do_X=x"))
    @test WD.output_ports(wd, doX) == Any[:X]

    # 4) key wires that MUST exist
    out_id = WD.output_id(wd)

    @test has_wire_pair(wd, (box_id_by_name(wd, :f_Y), 1), (out_id, 1))

    @test has_wire_pair(wd, (box_id_by_name(wd, :PU_R), 1),
                             (box_id_by_name(wd, :f_R), 1))

    @test has_wire_pair(wd, (box_id_by_name(wd, :f_R), 1),
                             (box_id_by_name(wd, :f_X), 2))

    @test has_wire_pair(wd, (box_id_by_name(wd, Symbol("do_X=x")), 1),
                             (box_id_by_name(wd, :f_W), 1))

    @test has_wire_pair(wd, (box_id_by_name(wd, :f_Z), 1),
                             (box_id_by_name(wd, Symbol("obs_Z=z")), 1))

    # 5) total wire count
    @test length(WD.wires(wd)) == 16
end
