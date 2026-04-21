include(joinpath(@__DIR__, "..", "src", "admg_compile.jl"))
include(joinpath(@__DIR__, "..", "src", "bn_import.jl"))

function usage()
    println("Usage:")
    println("  julia tools/import_bn_to_admg.jl <input.bif|input.net> [output.jl] [model_name]")
end

if isempty(ARGS)
    usage()
    exit(1)
end

input_path = ARGS[1]
output_path = length(ARGS) >= 2 ? ARGS[2] : ""
model_name = length(ARGS) >= 3 ? ARGS[3] : "model"

structure = read_bn_structure(input_path)

if isempty(output_path)
    write_admg_model(stdout, structure; model_name=model_name)
else
    write_admg_model(output_path, structure; model_name=model_name)
    println("Wrote $(length(structure.nodes)) nodes and $(length(structure.model.directed)) directed edges to $output_path")
end
