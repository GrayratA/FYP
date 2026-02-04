using Test
const WD = Catlab.WiringDiagrams

@testset "Step4.2 collapse fragment into prob box" begin
  wd = WiringDiagram([], Any[])

  # Build a tiny fragment with: 
  # doX -> cY -> obsY   plus latent root cR feeding U
  doX = WD.add_box!(wd, Box(:doX, Any[], Any[:X]))
  cR  = WD.add_box!(wd, Box(:cR, Any[], Any[:U]))
  cY  = WD.add_box!(wd, Box(:cY, Any[:X, :U], Any[:Y]))
  obsY = WD.add_box!(wd, Box(:obsY, Any[:Y], Any[]))

  WD.add_wire!(wd, (doX,1)=>(cY,1))
  WD.add_wire!(wd, (cR,1)=>(cY,2))
  WD.add_wire!(wd, (cY,1)=>(obsY,1))

  frag = [doX, cR, cY]  # collapse this fragment

  replace_fragment_with_prob_box_grouped!(wd, frag, 1)

  names = [string(box(wd,b).value) for b in box_ids(wd) if b!=input_id(wd) && b!=output_id(wd)]
  @test any(n -> startswith(n, "P_") || startswith(n, "P("), names)
end
