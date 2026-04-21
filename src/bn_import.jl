struct BNStructure
    model::ADMGModel
    nodes::Vector{Symbol}
end

function _strip_bn_comments(text::AbstractString)
    buf = IOBuffer()
    for line in split(text, '\n')
        cleaned = replace(line, r"//.*$" => "")
        cleaned = replace(cleaned, r"#.*$" => "")
        cleaned = replace(cleaned, r"%.*$" => "")
        println(buf, cleaned)
    end
    return String(take!(buf))
end

function _normalize_bn_name(name::AbstractString)
    stripped = strip(name)
    if startswith(stripped, "\"") && endswith(stripped, "\"") && ncodeunits(stripped) >= 2
        stripped = stripped[2:end-1]
    end

    normalized = replace(stripped, r"[^A-Za-z0-9_]" => "_")
    isempty(normalized) && error("Invalid empty node name in BN file")
    isnothing(match(r"^[A-Za-z_]", normalized)) && (normalized = "_" * normalized)
    return Symbol(normalized)
end

function _split_parent_list(raw::AbstractString)
    cleaned = strip(raw)
    isempty(cleaned) && return Symbol[]

    parts = String[]
    for m in eachmatch(r"\"[^\"]+\"|[^,\s]+", cleaned)
        push!(parts, m.match)
    end

    return [_normalize_bn_name(p) for p in parts]
end

function _parse_bif_structure(text::AbstractString)
    nodes = Symbol[]
    node_set = Set{Symbol}()
    directed = Pair{Symbol, Symbol}[]

    for m in eachmatch(r"variable\s+(\"[^\"]+\"|[^\s\{]+)\s*\{"s, text)
        node = _normalize_bn_name(m.captures[1])
        node in node_set || push!(nodes, node)
        push!(node_set, node)
    end

    for m in eachmatch(r"probability\s*\(\s*(\"[^\"]+\"|[^\s\|\)]+)\s*(?:\|\s*([^\)]*?))?\s*\)"s, text)
        child = _normalize_bn_name(m.captures[1])
        child in node_set || error("BIF references undeclared variable $child")
        parents = isnothing(m.captures[2]) ? Symbol[] : _split_parent_list(m.captures[2])
        for parent in parents
            parent in node_set || error("BIF references undeclared parent $parent")
            push!(directed, parent => child)
        end
    end

    return BNStructure(ADMGModel(unique(directed), Pair{Symbol, Symbol}[]), nodes)
end

function _parse_net_structure(text::AbstractString)
    nodes = Symbol[]
    node_set = Set{Symbol}()
    directed = Pair{Symbol, Symbol}[]

    for m in eachmatch(r"node\s+(\"[^\"]+\"|[^\s\{]+)\s*\{"s, text)
        node = _normalize_bn_name(m.captures[1])
        node in node_set || push!(nodes, node)
        push!(node_set, node)
    end

    for m in eachmatch(r"potential\s*\(\s*(\"[^\"]+\"|[^\s\|\)]+)\s*(?:\|\s*([^\)]*?))?\s*\)"s, text)
        child = _normalize_bn_name(m.captures[1])
        child in node_set || error(".net references undeclared node $child")
        parents = isnothing(m.captures[2]) ? Symbol[] : _split_parent_list(m.captures[2])
        for parent in parents
            parent in node_set || error(".net references undeclared parent $parent")
            push!(directed, parent => child)
        end
    end

    return BNStructure(ADMGModel(unique(directed), Pair{Symbol, Symbol}[]), nodes)
end

function read_bn_structure(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    ext in (".bif", ".net") || error("Unsupported BN file extension '$ext'. Expected .bif or .net.")

    text = _strip_bn_comments(read(path, String))
    if ext == ".bif"
        return _parse_bif_structure(text)
    end
    return _parse_net_structure(text)
end

function read_bn_as_admg(path::AbstractString)
    return read_bn_structure(path).model
end

function admg_model_literal(model::ADMGModel)
    pairs = ["$(repr(edge.first)) => $(repr(edge.second))" for edge in model.directed]
    lines = [
        "ADMGModel(",
        "    Pair{Symbol, Symbol}[",
        isempty(pairs) ? "" : "        " * join(pairs, ",\n        "),
        "    ],",
        "    Pair{Symbol, Symbol}[],",
        ")",
    ]
    return join(lines, "\n")
end

function write_admg_model(io::IO, structure::BNStructure; model_name::AbstractString="model")
    println(io, "# Auto-generated from a Bayesian network structure file")
    println(io, "$model_name = $(admg_model_literal(structure.model))")
    println(io, "nodes = $(repr(structure.nodes))")
end

function write_admg_model(path::AbstractString, structure::BNStructure; model_name::AbstractString="model")
    open(path, "w") do io
        write_admg_model(io, structure; model_name=model_name)
    end
    return path
end
