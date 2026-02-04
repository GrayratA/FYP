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
# ----------------------------


@testset "simplify-cf end-to-end: party example" begin
    model_party = ConfoundedModel(
      [
        :A => :B,
        :A => :C,
        :B => :S,
        :C => :S,
      ],
      Dict()
    )

    display_var = Set([:B, :S])
    base_scm = graph_b_to_scm(model_party; outputs=display_var)

    q_real = CounterfactualQuery(
        :Real,
        Dict(),                 # no interventions
        Dict(:B => :b),         # observe B=b
        Symbol[]                # outputs none
    )

    q_cf = CounterfactualQuery(
        :CF,
        Dict(:B => :doB),       # do_B=doB
        Dict{Symbol,Symbol}(),  # no observations
        [:B, :S]
    )

    queries = [q_real, q_cf]

    wd = build_multiverse(base_scm, queries)
    drop_discard_branches!(wd)
    merge_identical_deterministic_boxes(wd)
    replace_fx_with_cx!(wd)

    # 1) outputs exactly [:B,:S]
    @test WD.output_ports(wd) == Any[:B, :S]

    # 2) should have no f_* boxes after replace_fx_with_cx!
    @test !has_any_box_with_prefix(wd, "f_")

    # 3) boxes must exist (exact names from your printed simplified_wd)
    obsBb = box_id_by_name(wd, Symbol("obs_B=b"))
    cUS   = box_id_by_name(wd, :cUS)
    doB   = box_id_by_name(wd, Symbol("do_B=doB"))
    cUC   = box_id_by_name(wd, :cUC)
    cUA   = box_id_by_name(wd, :cUA)
    cUB   = box_id_by_name(wd, :cUB)

    # 4) port signatures must match your printed WD
    @test WD.input_ports(wd, obsBb) == Any[:B]
    @test WD.output_ports(wd, obsBb) == Any[]

    @test WD.input_ports(wd, cUS) == Any[:C, :B]
    @test WD.output_ports(wd, cUS) == Any[:S]

    @test WD.input_ports(wd, doB) == Any[]
    @test WD.output_ports(wd, doB) == Any[:B]

    @test WD.input_ports(wd, cUC) == Any[:A]
    @test WD.output_ports(wd, cUC) == Any[:C]

    @test WD.input_ports(wd, cUA) == Any[]
    @test WD.output_ports(wd, cUA) == Any[:A]

    @test WD.input_ports(wd, cUB) == Any[:A]
    @test WD.output_ports(wd, cUB) == Any[:B]

    # 5) wires must match your printed list
    out_id = WD.output_id(wd)

    @test has_wire_pair(wd, (cUC,1), (cUS,1))        # (4,1) => (2,1)
    @test has_wire_pair(wd, (doB,1), (cUS,2))        # (3,1) => (2,2)
    @test has_wire_pair(wd, (cUA,1), (cUC,1))        # (5,1) => (4,1)
    @test has_wire_pair(wd, (cUA,1), (cUB,1))        # (5,1) => (6,1)
    @test has_wire_pair(wd, (cUB,1), (obsBb,1))      # (6,1) => (1,1)
    @test has_wire_pair(wd, (doB,1), (out_id,1))     # (3,1) => (-1,1)
    @test has_wire_pair(wd, (cUS,1), (out_id,2))     # (2,1) => (-1,2)

    # 6) strict wire count
    @test length(WD.wires(wd)) == 7
end

@testset "simplify-cf end-to-end: drug example" begin
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
        Dict(:X => :doX_),          # interventions
        Dict{Symbol, Symbol}(),     # no observations
        [:Y]                        # output
    )

    q2 = CounterfactualQuery(
        :World2,
        Dict(),                     # no interventions
        Dict(:X => :x, :D => :d),   # observations
        Symbol[]
    )

    q3 = CounterfactualQuery(
        :World2,
        Dict(:D => :d),             # interventions
        Dict(:Z => :z),             # observations
        Symbol[]
    )

    full_wd = build_multiverse(base_scm, [q1, q2, q3])

    simplified_wd = full_wd
    drop_discard_branches!(simplified_wd)
    merge_identical_deterministic_boxes(simplified_wd)
    replace_fx_with_cx!(simplified_wd)

    # 1) outputs exactly [:Y]
    @test WD.output_ports(simplified_wd) == Any[:Y]

    # 2) should have no f_* boxes after replace_fx_with_cx!
    @test !has_any_box_with_prefix(simplified_wd, "f_")

    # 3) boxes must exist (exact names from your expected printed WD)
    cUY   = box_id_by_name(simplified_wd, :cUY)
    cUX   = box_id_by_name(simplified_wd, :cUX)
    obsXx = box_id_by_name(simplified_wd, Symbol("obs_X=x"))
    obsDd = box_id_by_name(simplified_wd, Symbol("obs_D=d"))
    doZ   = box_id_by_name(simplified_wd, :doZ)
    doDd  = box_id_by_name(simplified_wd, Symbol("do_D=d"))
    cUD   = box_id_by_name(simplified_wd, :cUD)
    cUR   = box_id_by_name(simplified_wd, :cUR)
    cUW   = box_id_by_name(simplified_wd, :cUW)
    obsZz = box_id_by_name(simplified_wd, Symbol("obs_Z=z"))
    cUZ   = box_id_by_name(simplified_wd, :cUZ)
    doX   = box_id_by_name(simplified_wd, Symbol("do_X=doX_"))

    # 4) port signatures must match your expected printed WD
    @test WD.input_ports(simplified_wd, cUY) == Any[:R, :W, :Z]
    @test WD.output_ports(simplified_wd, cUY) == Any[:Y]

    @test WD.input_ports(simplified_wd, cUX) == Any[:R]
    @test WD.output_ports(simplified_wd, cUX) == Any[:X]

    @test WD.input_ports(simplified_wd, obsXx) == Any[:X]
    @test WD.output_ports(simplified_wd, obsXx) == Any[]

    @test WD.input_ports(simplified_wd, obsDd) == Any[:D]
    @test WD.output_ports(simplified_wd, obsDd) == Any[]

    @test WD.input_ports(simplified_wd, doZ) == Any[]
    @test WD.output_ports(simplified_wd, doZ) == Any[:Z]

    @test WD.input_ports(simplified_wd, doDd) == Any[]
    @test WD.output_ports(simplified_wd, doDd) == Any[:D]

    @test WD.input_ports(simplified_wd, cUD) == Any[]
    @test WD.output_ports(simplified_wd, cUD) == Any[:D]

    @test WD.input_ports(simplified_wd, cUR) == Any[]
    @test WD.output_ports(simplified_wd, cUR) == Any[:R]

    @test WD.input_ports(simplified_wd, cUW) == Any[:X]
    @test WD.output_ports(simplified_wd, cUW) == Any[:W]

    @test WD.input_ports(simplified_wd, obsZz) == Any[:Z]
    @test WD.output_ports(simplified_wd, obsZz) == Any[]

    @test WD.input_ports(simplified_wd, cUZ) == Any[:D]
    @test WD.output_ports(simplified_wd, cUZ) == Any[:Z]

    @test WD.input_ports(simplified_wd, doX) == Any[]
    @test WD.output_ports(simplified_wd, doX) == Any[:X]

    # 5) wires must match your expected printed list
    out_id = WD.output_id(simplified_wd)

    @test has_wire_pair(simplified_wd, (doDd,1), (cUZ,1))     # (6,1) => (11,1)
    @test has_wire_pair(simplified_wd, (doZ,1),  (cUY,3))     # (5,1) => (1,3)
    @test has_wire_pair(simplified_wd, (cUX,1),  (obsXx,1))   # (2,1) => (3,1)
    @test has_wire_pair(simplified_wd, (cUZ,1),  (obsZz,1))   # (11,1) => (10,1)
    @test has_wire_pair(simplified_wd, (cUR,1),  (cUY,1))     # (8,1) => (1,1)
    @test has_wire_pair(simplified_wd, (cUW,1),  (cUY,2))     # (9,1) => (1,2)
    @test has_wire_pair(simplified_wd, (cUD,1),  (obsDd,1))   # (7,1) => (4,1)
    @test has_wire_pair(simplified_wd, (doX,1),  (cUW,1))     # (12,1) => (9,1)
    @test has_wire_pair(simplified_wd, (cUR,1),  (cUX,1))     # (8,1) => (2,1)
    @test has_wire_pair(simplified_wd, (cUY,1),  (out_id,1))  # (1,1) => (-1,1)

    # 6) strict wire count
    @test length(WD.wires(simplified_wd)) == 10
end
