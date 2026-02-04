using Test
using Catlab
using Catlab.WiringDiagrams
const WD = Catlab.WiringDiagrams

include(joinpath(@__DIR__, "..", "..", "src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "simplify_cf.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "id_cf.jl"))

# ---------------- helpers ----------------
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

function has_any_box_with_prefix(wd::WiringDiagram, prefix::String)::Bool
    for b in WD.box_ids(wd)
        b == WD.input_id(wd)  && continue
        b == WD.output_id(wd) && continue
        n = WD.box(wd, b).value
        if n isa Symbol && startswith(String(n), prefix)
            return true
        end
    end
    return false
end
# ----------------------------------------


@testset "replace_fx_with_cx! on multiverse fig15b yields expected wd" begin
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
    replace_fx_with_cx!(wd)

    # 1) outputs must be [:Y]
    @test WD.output_ports(wd) == Any[:Y]

    # 2) no more f_ boxes after replace_fx_with_cx!
    @test !has_any_box_with_prefix(wd, "f_")

    # 3) expected boxes exist (names exactly as your print)
    expected_boxes = [
        :cUY,
        :cUX,
        Symbol("obs_X=x0"),
        Symbol("obs_D=d"),
        :doZ,
        Symbol("do_D=d"),
        :cUD,
        :cUR,
        :cUW,
        Symbol("obs_Z=z"),
        :cUZ,
        Symbol("do_X=x"),
    ]

    for nm in expected_boxes
        @test box_id_by_name(wd, nm) isa Int
    end

    # 4) spot-check port signatures (match your print)
    cUY = box_id_by_name(wd, :cUY)
    @test WD.input_ports(wd, cUY)  == Any[:R, :W, :Z]
    @test WD.output_ports(wd, cUY) == Any[:Y]

    cUX = box_id_by_name(wd, :cUX)
    @test WD.input_ports(wd, cUX)  == Any[:R]
    @test WD.output_ports(wd, cUX) == Any[:X]

    cUW = box_id_by_name(wd, :cUW)
    @test WD.input_ports(wd, cUW)  == Any[:X]
    @test WD.output_ports(wd, cUW) == Any[:W]

    cUZ = box_id_by_name(wd, :cUZ)
    @test WD.input_ports(wd, cUZ)  == Any[:D]
    @test WD.output_ports(wd, cUZ) == Any[:Z]

    # 5) verify ALL wires you printed exist
    out_id = WD.output_id(wd)

    doD = box_id_by_name(wd, Symbol("do_D=d"))
    doX = box_id_by_name(wd, Symbol("do_X=x"))
    doZ_id = box_id_by_name(wd, :doZ)

    cUD = box_id_by_name(wd, :cUD)
    cUR = box_id_by_name(wd, :cUR)
    cUX = box_id_by_name(wd, :cUX)
    cUW = box_id_by_name(wd, :cUW)
    cUZ = box_id_by_name(wd, :cUZ)
    obsX = box_id_by_name(wd, Symbol("obs_X=x0"))
    obsD = box_id_by_name(wd, Symbol("obs_D=d"))
    obsZ = box_id_by_name(wd, Symbol("obs_Z=z"))

    @test has_wire_pair(wd, (doD, 1), (cUZ, 1))          # (6,1) => (11,1)
    @test has_wire_pair(wd, (doZ_id, 1), (cUY, 3))       # (5,1) => (1,3)
    @test has_wire_pair(wd, (cUX, 1), (obsX, 1))         # (2,1) => (3,1)
    @test has_wire_pair(wd, (cUZ, 1), (obsZ, 1))         # (11,1)=> (10,1)
    @test has_wire_pair(wd, (cUR, 1), (cUY, 1))          # (8,1) => (1,1)
    @test has_wire_pair(wd, (cUW, 1), (cUY, 2))          # (9,1) => (1,2)
    @test has_wire_pair(wd, (cUD, 1), (obsD, 1))         # (7,1) => (4,1)
    @test has_wire_pair(wd, (doX, 1), (cUW, 1))          # (12,1)=> (9,1)
    @test has_wire_pair(wd, (cUR, 1), (cUX, 1))          # (8,1) => (2,1)
    @test has_wire_pair(wd, (cUY, 1), (out_id, 1))       # (1,1) => (-1,1)

    # 6) optional strict check: wire count matches your print (=10)
    @test length(WD.wires(wd)) == 10
end
