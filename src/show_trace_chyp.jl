trace_dir = joinpath(@__DIR__, "trace")

if !isdir(trace_dir)
    error("trace directory not found: $trace_dir")
end

files = sort(filter(f -> endswith(f, ".chyp"), readdir(trace_dir)))

if isempty(files)
    println("No .chyp files found in $trace_dir")
else
    for (i, file) in enumerate(files)
        path = joinpath(trace_dir, file)
        println("="^80)
        println(file)
        println("="^80)
        println(read(path, String))
        if i != length(files)
            println()
        end
    end
end
