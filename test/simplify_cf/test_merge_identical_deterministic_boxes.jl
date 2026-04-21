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

function count_boxes_by_name(wd::WiringDiagram, name::Symbol)::Int
    count = 0
    for b in WD.box_ids(wd)
        b == WD.input_id(wd)  && continue
        b == WD.output_id(wd) && continue
        WD.box(wd, b).value == name && (count += 1)
    end
    return count
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
    separate_sharp_effect_copy_maps!(wd)
    merge_identical_deterministic_boxes(wd)

    # 1) outputs
    @test WD.output_ports(wd) == Any[:Y]

    # 2) expected boxes exist
    expected_boxes = [
        :PU_D, :PU_R, :PU_W, :PU_X, :PU_Y, :PU_Z,
        :f_D, :f_R, :f_W, :f_Z, :f_X, :f_Y,
        Symbol("do_X=x"),
        Symbol("do_D=d"),
        Symbol("do_Z=z"),
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

@testset "merge_identical_deterministic_boxes requires truly shared inputs" begin
    wd = WiringDiagram([], Any[:Y, :Y])

    src_w1 = WD.add_box!(wd, Box(:src_w1, Any[], Any[:W]))
    src_w2 = WD.add_box!(wd, Box(:src_w2, Any[], Any[:W]))
    src_z  = WD.add_box!(wd, Box(:src_z, Any[], Any[:Z]))
    src_r  = WD.add_box!(wd, Box(:src_r, Any[], Any[:R]))
    src_u  = WD.add_box!(wd, Box(:PU_Y, Any[], Any[:UY]))

    fy1 = WD.add_box!(wd, Box(:f_Y, Any[:W, :Z, :R, :UY], Any[:Y]))
    fy2 = WD.add_box!(wd, Box(:f_Y, Any[:W, :Z, :R, :UY], Any[:Y]))

    WD.add_wire!(wd, (src_w1, 1) => (fy1, 1))
    WD.add_wire!(wd, (src_w2, 1) => (fy2, 1))
    WD.add_wire!(wd, (src_z, 1) => (fy1, 2))
    WD.add_wire!(wd, (src_z, 1) => (fy2, 2))
    WD.add_wire!(wd, (src_r, 1) => (fy1, 3))
    WD.add_wire!(wd, (src_r, 1) => (fy2, 3))
    WD.add_wire!(wd, (src_u, 1) => (fy1, 4))
    WD.add_wire!(wd, (src_u, 1) => (fy2, 4))
    WD.add_wire!(wd, (fy1, 1) => (WD.output_id(wd), 1))
    WD.add_wire!(wd, (fy2, 1) => (WD.output_id(wd), 2))

    merge_identical_deterministic_boxes(wd)

    @test count_boxes_by_name(wd, :f_Y) == 2
end

@testset "merge_identical_deterministic_boxes does not stop after first non-mergeable name" begin
    wd = WiringDiagram([], Any[:A, :A, :B, :B])

    src_x1 = WD.add_box!(wd, Box(:src_x1, Any[], Any[:X]))
    src_x2 = WD.add_box!(wd, Box(:src_x2, Any[], Any[:X]))
    src_z  = WD.add_box!(wd, Box(:src_z, Any[], Any[:Z]))

    fa1 = WD.add_box!(wd, Box(:f_A, Any[:X], Any[:A]))
    fa2 = WD.add_box!(wd, Box(:f_A, Any[:X], Any[:A]))
    fb1 = WD.add_box!(wd, Box(:f_B, Any[:Z], Any[:B]))
    fb2 = WD.add_box!(wd, Box(:f_B, Any[:Z], Any[:B]))

    WD.add_wire!(wd, (src_x1, 1) => (fa1, 1))
    WD.add_wire!(wd, (src_x2, 1) => (fa2, 1))
    WD.add_wire!(wd, (src_z, 1) => (fb1, 1))
    WD.add_wire!(wd, (src_z, 1) => (fb2, 1))

    WD.add_wire!(wd, (fa1, 1) => (WD.output_id(wd), 1))
    WD.add_wire!(wd, (fa2, 1) => (WD.output_id(wd), 2))
    WD.add_wire!(wd, (fb1, 1) => (WD.output_id(wd), 3))
    WD.add_wire!(wd, (fb2, 1) => (WD.output_id(wd), 4))

    merge_identical_deterministic_boxes(wd)

    @test count_boxes_by_name(wd, :f_A) == 2
    @test count_boxes_by_name(wd, :f_B) == 1
end
