using Test
using Catlab
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Theories
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics: Graphviz, to_graphviz

# ===============================
# 1. include
# ===============================
include(joinpath(@__DIR__, "..", "src", "utils.jl"))
include(joinpath(@__DIR__, "..", "src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "src", "simplify_cf.jl"))
include(joinpath(@__DIR__, "..", "src", "id_cf.jl"))


# ===============================
# 2. admg_compile tests
# ===============================
@testset "admg_compile" begin
    include(joinpath(@__DIR__, "admg_compile", "test_admg.jl"))
    include(joinpath(@__DIR__, "admg_compile", "test_build_multiverse_labels.jl"))
end

# ===============================s
# 2. simplified-cf tests
# ===============================
@testset "simplified-cf" begin
    include(joinpath(@__DIR__, "simplify_cf", "test_merge_identical_deterministic_boxes.jl"))
    include(joinpath(@__DIR__, "simplify_cf", "test_replace_fx_with_cx.jl"))
    include(joinpath(@__DIR__, "simplify_cf", "test_simplify_cf_general_test.jl"))
end


# ===============================
# 3. id-cf Step 4 tests (all cases)
# ===============================
@testset "id-cf Step4" begin
    include(joinpath(@__DIR__, "id_cf", "test_step41.jl"))
    include(joinpath(@__DIR__, "id_cf", "test_step42_collapse_prob_box.jl"))
    include(joinpath(@__DIR__, "id_cf", "test_Rfragements.jl"))
    include(joinpath(@__DIR__, "id_cf", "test_id_cf_general_test.jl"))
    
end
