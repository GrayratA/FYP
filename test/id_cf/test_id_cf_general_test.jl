using Test
using Catlab
using Catlab.WiringDiagrams
const WD = Catlab.WiringDiagrams

include(joinpath(@__DIR__, "..", "..", "src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "simplify_cf.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "id_cf.jl"))

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

@testset "id-cf: drug example" begin
    fig_15b_model = ConfoundedModel(
        [
            :X => :W,
            :W => :Y,
            :D => :Z,
            :Z => :Y
        ],
        Dict(
            :R => [:X, :Y]
        )
    )

    display_var = Set([:D, :Z, :X, :W, :Y])
    base_scm = graph_b_to_scm(fig_15b_model; outputs=display_var)

    q1 = CounterfactualQuery(
        :World1,
        Dict(:X => :doX_),
        Dict{Symbol, Symbol}(),
        [:Y]
    )

    q2 = CounterfactualQuery(
        :World2,
        Dict(),
        Dict(:X => :x, :D => :d),
        Symbol[]
    )

    q3 = CounterfactualQuery(
        :World2,
        Dict(:D => :d),
        Dict(:Z => :z),
        Symbol[]
    )

    full_wd = build_multiverse(base_scm, [q1, q2, q3])

    simplified_wd = full_wd
    drop_discard_branches!(simplified_wd)
    merge_identical_deterministic_boxes(simplified_wd)
    replace_fx_with_cx!(simplified_wd)

    wd = id_cf_step4!(simplified_wd, display_var)

    # 1) outputs exactly [:Y]
    @test WD.output_ports(wd) == Any[:Y]

    # 2) expected boxes exist
    P4   = box_id_by_name(wd, Symbol("P_4(Z ; do(D))"))
    doX  = box_id_by_name(wd, Symbol("do_X=doX_"))
    obsX = box_id_by_name(wd, Symbol("obs_X=x"))
    obsD = box_id_by_name(wd, Symbol("obs_D=d"))
    doZ  = box_id_by_name(wd, :doZ)
    doD  = box_id_by_name(wd, Symbol("do_D=d"))
    P1   = box_id_by_name(wd, Symbol("P_1(X,Y ; do(Z),do(W))"))
    P2   = box_id_by_name(wd, Symbol("P_2(D ; ∅)"))
    P3   = box_id_by_name(wd, Symbol("P_3(W ; do(X))"))
    obsZ = box_id_by_name(wd, Symbol("obs_Z=z"))

    # 3) port signatures must match printed WD
    @test WD.input_ports(wd, P4) == Any[:D]
    @test WD.output_ports(wd, P4) == Any[:Z]

    @test WD.input_ports(wd, doX) == Any[]
    @test WD.output_ports(wd, doX) == Any[:X]

    @test WD.input_ports(wd, obsX) == Any[:X]
    @test WD.output_ports(wd, obsX) == Any[]

    @test WD.input_ports(wd, obsD) == Any[:D]
    @test WD.output_ports(wd, obsD) == Any[]

    @test WD.input_ports(wd, doZ) == Any[]
    @test WD.output_ports(wd, doZ) == Any[:Z]

    @test WD.input_ports(wd, doD) == Any[]
    @test WD.output_ports(wd, doD) == Any[:D]

    @test WD.input_ports(wd, P1) == Any[:W, :Z]
    @test WD.output_ports(wd, P1) == Any[:X, :Y]

    @test WD.input_ports(wd, P2) == Any[]
    @test WD.output_ports(wd, P2) == Any[:D]

    @test WD.input_ports(wd, P3) == Any[:X]
    @test WD.output_ports(wd, P3) == Any[:W]

    @test WD.input_ports(wd, obsZ) == Any[:Z]
    @test WD.output_ports(wd, obsZ) == Any[]

    # 4) wires must match printed list
    out_id = WD.output_id(wd)

    @test has_wire_pair(wd, (P2,1),  (obsD,1))     # (8,1) => (4,1)
    @test has_wire_pair(wd, (P3,1),  (P1,1))       # (9,1) => (7,1)
    @test has_wire_pair(wd, (P1,1),  (obsX,1))     # (7,1) => (3,1)
    @test has_wire_pair(wd, (doX,1), (P3,1))       # (2,1) => (9,1)
    @test has_wire_pair(wd, (doZ,1), (P1,2))       # (5,1) => (7,2)
    @test has_wire_pair(wd, (doD,1), (P4,1))       # (6,1) => (1,1)
    @test has_wire_pair(wd, (P4,1),  (obsZ,1))     # (1,1) => (10,1)
    @test has_wire_pair(wd, (P1,2),  (out_id,1))   # (7,2) => (-1,1)

    # 5) strict wire count
    @test length(WD.wires(wd)) == 8

    # 6) (optional) sanity: should still have no f_* boxes
    @test !has_any_box_with_prefix(wd, "f_")
end
