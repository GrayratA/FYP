using Test
using Catlab
using Catlab.WiringDiagrams
const WD = Catlab.WiringDiagrams

include(joinpath(@__DIR__, "..", "..", "src", "admg_compile.jl"))  # 视你的项目调整

# --- helpers ---
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

function has_wire_pair(wd::WiringDiagram, src::Tuple{Int,Int}, tgt::Tuple{Int,Int})::Bool
    for w in WD.wires(wd)
        if (w.source.box, w.source.port) == src && (w.target.box, w.target.port) == tgt
            return true
        end
    end
    return false
end

@testset "graph_b_to_scm test 1" begin
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
    wd = graph_b_to_scm(model_party; outputs=display_var)

    # 1) output ports exactly [:B,:S]
    @test WD.output_ports(wd) == Any[:B, :S]

    # 2) all expected boxes exist with expected port signatures
    fA  = box_id_by_name(wd, :f_A)
    puA = box_id_by_name(wd, :PU_A)
    fB  = box_id_by_name(wd, :f_B)
    puB = box_id_by_name(wd, :PU_B)
    fC  = box_id_by_name(wd, :f_C)
    puC = box_id_by_name(wd, :PU_C)
    fS  = box_id_by_name(wd, :f_S)
    puS = box_id_by_name(wd, :PU_S)

    @test WD.input_ports(wd, fA)  == Any[:UA]
    @test WD.output_ports(wd, fA) == Any[:A]

    @test WD.input_ports(wd, puA)  == Any[]
    @test WD.output_ports(wd, puA) == Any[:UA]

    @test WD.input_ports(wd, fB)  == Any[:A, :UB]
    @test WD.output_ports(wd, fB) == Any[:B]

    @test WD.input_ports(wd, puB)  == Any[]
    @test WD.output_ports(wd, puB) == Any[:UB]

    @test WD.input_ports(wd, fC)  == Any[:A, :UC]
    @test WD.output_ports(wd, fC) == Any[:C]

    @test WD.input_ports(wd, puC)  == Any[]
    @test WD.output_ports(wd, puC) == Any[:UC]

    @test WD.input_ports(wd, fS)  == Any[:B, :C, :US]
    @test WD.output_ports(wd, fS) == Any[:S]

    @test WD.input_ports(wd, puS)  == Any[]
    @test WD.output_ports(wd, puS) == Any[:US]

    # 3) expected wires (exactly the ones you printed)
    out_id = WD.output_id(wd)

    @test has_wire_pair(wd, (puA,1), (fA,1))
    @test has_wire_pair(wd, (puB,1), (fB,2))
    @test has_wire_pair(wd, (puC,1), (fC,2))
    @test has_wire_pair(wd, (puS,1), (fS,3))

    @test has_wire_pair(wd, (fA,1), (fB,1))
    @test has_wire_pair(wd, (fA,1), (fC,1))
    @test has_wire_pair(wd, (fB,1), (fS,1))
    @test has_wire_pair(wd, (fC,1), (fS,2))

    @test has_wire_pair(wd, (fB,1), (out_id,1))   # B -> output port 1
    @test has_wire_pair(wd, (fS,1), (out_id,2))   # S -> output port 2

    # 4) optional: ensure wire count matches exactly 10 (like your print)
    @test length(WD.wires(wd)) == 10
end

@testset "graph_b_to_scm test 2" begin
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

    wd = graph_b_to_scm(fig_15b_model; outputs=display_var)

    # 1) outputs should match the printed order (sorted)
    @test WD.output_ports(wd) == Any[:D, :W, :X, :Y, :Z]

    # 2) expected boxes exist
    fD  = box_id_by_name(wd, :f_D)
    puD = box_id_by_name(wd, :PU_D)

    fR  = box_id_by_name(wd, :f_R)
    puR = box_id_by_name(wd, :PU_R)

    fW  = box_id_by_name(wd, :f_W)
    puW = box_id_by_name(wd, :PU_W)

    fX  = box_id_by_name(wd, :f_X)
    puX = box_id_by_name(wd, :PU_X)

    fY  = box_id_by_name(wd, :f_Y)
    puY = box_id_by_name(wd, :PU_Y)

    fZ  = box_id_by_name(wd, :f_Z)
    puZ = box_id_by_name(wd, :PU_Z)

    # sanity: total real boxes = 12
    real_boxes = [b for b in WD.box_ids(wd) if b != WD.input_id(wd) && b != WD.output_id(wd)]
    @test length(real_boxes) == 12

    # 3) port signatures (must match your printed WD)
    @test WD.input_ports(wd, fD)  == Any[:UD]
    @test WD.output_ports(wd, fD) == Any[:D]
    @test WD.input_ports(wd, puD) == Any[]
    @test WD.output_ports(wd, puD) == Any[:UD]

    @test WD.input_ports(wd, fR)  == Any[:UR]
    @test WD.output_ports(wd, fR) == Any[:R]
    @test WD.input_ports(wd, puR) == Any[]
    @test WD.output_ports(wd, puR) == Any[:UR]

    @test WD.input_ports(wd, fW)  == Any[:X, :UW]
    @test WD.output_ports(wd, fW) == Any[:W]
    @test WD.input_ports(wd, puW) == Any[]
    @test WD.output_ports(wd, puW) == Any[:UW]

    @test WD.input_ports(wd, fX)  == Any[:R, :UX]
    @test WD.output_ports(wd, fX) == Any[:X]
    @test WD.input_ports(wd, puX) == Any[]
    @test WD.output_ports(wd, puX) == Any[:UX]

    @test WD.input_ports(wd, fY)  == Any[:W, :Z, :R, :UY]
    @test WD.output_ports(wd, fY) == Any[:Y]
    @test WD.input_ports(wd, puY) == Any[]
    @test WD.output_ports(wd, puY) == Any[:UY]

    @test WD.input_ports(wd, fZ)  == Any[:D, :UZ]
    @test WD.output_ports(wd, fZ) == Any[:Z]
    @test WD.input_ports(wd, puZ) == Any[]
    @test WD.output_ports(wd, puZ) == Any[:UZ]

    # 4) exogenous wires PU_* -> f_* (last input port)
    @test has_wire_pair(wd, (puD,1), (fD,1))
    @test has_wire_pair(wd, (puR,1), (fR,1))
    @test has_wire_pair(wd, (puW,1), (fW,2))
    @test has_wire_pair(wd, (puX,1), (fX,2))
    @test has_wire_pair(wd, (puY,1), (fY,4))
    @test has_wire_pair(wd, (puZ,1), (fZ,2))

    # 5) causal + latent edges (match your printed wires)
    @test has_wire_pair(wd, (fX,1), (fW,1))   # X -> W
    @test has_wire_pair(wd, (fW,1), (fY,1))   # W -> Y
    @test has_wire_pair(wd, (fD,1), (fZ,1))   # D -> Z
    @test has_wire_pair(wd, (fZ,1), (fY,2))   # Z -> Y

    @test has_wire_pair(wd, (fR,1), (fX,1))   # R -> X (latent)
    @test has_wire_pair(wd, (fR,1), (fY,3))   # R -> Y (latent)

    # 6) wires to external outputs (-1, i)
    out_id = WD.output_id(wd)
    @test has_wire_pair(wd, (fD,1), (out_id,1))  # D
    @test has_wire_pair(wd, (fW,1), (out_id,2))  # W
    @test has_wire_pair(wd, (fX,1), (out_id,3))  # X
    @test has_wire_pair(wd, (fY,1), (out_id,4))  # Y
    @test has_wire_pair(wd, (fZ,1), (out_id,5))  # Z

    # optional strict check: wire count equals your print (17)
    @test length(WD.wires(wd)) == 17
end