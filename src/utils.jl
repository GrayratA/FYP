const WD = Catlab.WiringDiagrams

using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics: Graphviz, to_graphviz

# domain size  = number of boundary input wires
n_in_boundary(wd::WiringDiagram)  = length(WD.output_ports(wd))

# codomain size = number of boundary output wires
n_out_boundary(wd::WiringDiagram) = length(WD.input_ports(wd))


draw(d::WiringDiagram) = to_graphviz(d;
  orientation=LeftToRight,
  labels=true,
  label_attr=:xlabel,
  node_attrs=Graphviz.Attributes(:fontname => "Courier"),
  edge_attrs=Graphviz.Attributes(:fontname => "Courier")
)


# get box name
box_name(wd::WiringDiagram, b::Int) = begin
    bx = WD.box(wd, b)
    bx === nothing && return nothing
    v = bx.value
    if v isa Symbol
        return v
    elseif v isa NamedTuple
        return (:name ∈ keys(v)) ? v.name : nothing
    elseif v isa Pair
        return v.first
    else
        return nothing
    end
end

function safe_box_name(wd::WiringDiagram, b::Int)
    try
        return box_name(wd, b)
    catch e
        return nothing
    end
end

# count input/output ports
function nin(wd, b::Int)
    # input / output vertices (-2, -1) have no ports
    b <= 0 && return 0
    length(input_ports(wd, b))
end

function nout(wd, b::Int)
    b <= 0 && return 0
    length(output_ports(wd, b))
end


# check if a box is a sharp effect
function is_sharp_effect(wd, b::Int)
    # output / input vertices are never sharp effects
    b <= 0 && return false

    n = box_name(wd, b)
    n isa Symbol || return false
    startswith(String(n), "obs") || return false
    nin(wd, b) == 1 || return false
    nout(wd, b) == 0 || return false
    return true
end




# detect discard boxes: have inputs but no outputs
is_discard_box(wd::WiringDiagram, b::Int) = (nout(wd, b) == 0) && (nin(wd, b) >= 1) && !is_sharp_effect(wd, b)

# check if a box has no consumers (no outgoing wires)
function no_consumers(wd::WiringDiagram, b::Int)
    isempty(out_wires(wd, b))
end

#chech if the box is isolated
function is_lonely_box(wd::WiringDiagram, b::Int)
    nin  = length(input_ports(wd, b))
    nout = length(output_ports(wd, b))

    if nin == 0 && nout == 0
        return true
    end

    connected_to_other_box = false

    for j in 1:nin
        for w in in_wires(wd, b, j)
            sbox = w.source.box
            if sbox != 0 && sbox != b
                connected_to_other_box = true
            end
        end
    end

    for j in 1:nout
        for w in out_wires(wd, b, j)
            tbox = w.target.box
            if tbox != 0 && tbox != b
                connected_to_other_box = true
            end
        end
    end

    return !connected_to_other_box
end

# remove all the isolated box
function remove_lonely_boxes!(wd::WiringDiagram)
    lonely_boxes = Int[]
    for b in 1:WD.nboxes(wd)
        if is_lonely_box(wd, b)
            push!(lonely_boxes, b)
        end
    end

    # Delete in reverse order to avoid index changes
    sort!(lonely_boxes, rev = true)
    for b in lonely_boxes
        WD.rem_box!(wd, b)
    end

    return wd
end

is_state_box(wd::WiringDiagram, b::Int) =
    nin(wd,b) == 0 && nout(wd,b) == 1 && begin
        n = safe_box_name(wd,b)
        n isa Symbol && startswith(String(n), "PU")
    end

canonical_c_name(var::Symbol) = Symbol("c_", var)
canonical_do_name(var::Symbol) = Symbol("do_", var)

function state_input_ports(wd::WiringDiagram, b::Int)
    ports = Int[]
    for j in 1:nin(wd, b)
        ws = in_wires(wd, b, j)
        any(w -> is_state_box(wd, w.source.box), ws) && push!(ports, j)
    end
    return ports
end

function lambda_is_directly_fed(wd::WiringDiagram, b::Int)
    ports = state_input_ports(wd, b)
    isempty(ports) && return false

    for j in ports
        ws = in_wires(wd, b, j)
        length(ws) == 1 || return false

        w = only(ws)
        is_state_box(wd, w.source.box) || return false

        src_ws = out_wires(wd, w.source.box, w.source.port)
        length(src_ws) == 1 || return false

        src_w = only(src_ws)
        same_edge =
            src_w.source.box == w.source.box &&
            src_w.source.port == w.source.port &&
            src_w.target.box == w.target.box &&
            src_w.target.port == w.target.port
        same_edge || return false
    end

    return true
end

function unabsorbed_lambda_boxes(wd::WiringDiagram)
    boxes = Int[]
    for b in WD.box_ids(wd)
        n = safe_box_name(wd, b)
        n isa Symbol && startswith(String(n), "f_") || continue
        !isempty(state_input_ports(wd, b)) && push!(boxes, b)
    end
    return boxes
end

function topological_order(wd::WiringDiagram)
    boxes = WD.box_ids(wd)
    
    indeg = Dict(b => length(WD.in_wires(wd, b)) for b in boxes)

    for b in boxes
        indeg[b] = count(w -> w.source.box != WD.input_id(wd), WD.in_wires(wd, b))
    end

    order = Int[]
    queue = [b for b in boxes if indeg[b] == 0]

    while !isempty(queue)
        v = popfirst!(queue)
        push!(order, v)

        for w in WD.out_wires(wd, v)
            tgt = w.target.box
            tgt == WD.output_id(wd) && continue
            indeg[tgt] -= 1
            if indeg[tgt] == 0
                push!(queue, tgt)
            end
        end
    end

    return order
end

function duplicate_mechanism_groups(wd::WiringDiagram)
    boxes_all = WD.box_ids(wd)

    indeg = Dict(b => 0 for b in boxes_all)
    for b in boxes_all
        indeg[b] = count(w -> w.source.box != WD.input_id(wd), WD.in_wires(wd, b))
    end

    order = Int[]
    queue = sort([b for b in boxes_all if indeg[b] == 0])
    while !isempty(queue)
        v = popfirst!(queue)
        push!(order, v)
        for w in WD.out_wires(wd, v)
            tgt = w.target.box
            tgt == WD.output_id(wd) && continue
            indeg[tgt] -= 1
            if indeg[tgt] == 0
                push!(queue, tgt)
                sort!(queue)
            end
        end
    end

    topo_pos = Dict{Int,Int}()
    for (i, b) in enumerate(order)
        topo_pos[b] = i
    end

    groups = Dict{Symbol,Vector{Int}}()
    for b in boxes_all
        n = safe_box_name(wd, b)
        n isa Symbol && startswith(String(n), "f_") || continue
        push!(get!(groups, n, Int[]), b)
    end

    filtered = [(name, sort(boxes)) for (name, boxes) in groups if length(boxes) > 1]
    sort!(filtered; by = pair -> minimum(get(topo_pos, b, typemax(Int)) for b in pair[2]))
    return filtered
end

function find_same_name_boxes(wd::WiringDiagram)
    groups = Dict{Symbol,Vector{Int}}()
    for (name, boxes) in duplicate_mechanism_groups(wd)
        groups[name] = boxes
    end
    return groups
end


# find input source set for a box
function input_sources_set(wd::WiringDiagram, b::Int)
    S = Set{Symbol}()
    nin = length(input_ports(wd, b))
    for j in 1:nin
        ws = in_wires(wd, b, j)
        for w in ws
            var = port_value(wd, w.source)
            push!(S, var)
        end
    end
    return S
end

function merge_source_signature(wd::WiringDiagram, src_box::Int, src_port::Int)
    if is_sharp_state(wd, src_box)
        return (:sharp_state, safe_box_name(wd, src_box))
    else
        return (:wire, src_box, src_port)
    end
end

function input_source_signature(wd::WiringDiagram, b::Int)
    sig = Tuple{Int,Vector{Any}}[]
    nin = length(input_ports(wd, b))
    for j in 1:nin
        sources = Any[]
        for w in in_wires(wd, b, j)
            push!(sources, merge_source_signature(wd, w.source.box, w.source.port))
        end
        sort!(sources)
        push!(sig, (j, sources))
    end
    return sig
end

function merge_signature(wd::WiringDiagram, b::Int)
    # Paper step 3(alpha) only merges mechanisms that truly share their inputs.
    # In this representation that means matching the exact source endpoint(s)
    # feeding each input port, not just matching variable types.
    return input_source_signature(wd, b)
end

# Parsing Box Type and Variable Name
function get_node_info(box_val::Any)
    val_str = string(box_val)
    if startswith(val_str, "f_")
        return :Mechanism, Symbol(replace(val_str, "f_" => ""))
    elseif startswith(val_str, "PU")
        return :Exogenous, Symbol(val_str)
    else
        return :Other, Symbol(val_str)
    end
end

is_c_mechanism(wd::WiringDiagram, b::Int) = begin
    n = safe_box_name(wd, b)
    n isa Symbol && (startswith(String(n), "c_") || startswith(String(n), "cU"))
end

is_root_mechanism(wd::WiringDiagram, b::Int) = begin
    n = safe_box_name(wd, b)
    n isa Symbol && (
        startswith(String(n), "c_R") ||
        startswith(String(n), "cUR") ||
        startswith(String(n), "cR")
    )
end

is_cX_box(wd, b::Int, var_type::Symbol) = begin
    n = safe_box_name(wd, b)
    old_name = Symbol("cU", var_type)
    new_name = canonical_c_name(var_type)
    n isa Symbol && (n == new_name || n == old_name)
end

function is_sharp_state(wd, b::Int)
    # Support both canonical do_X labels and older doX labels.
    b <= 0 && return false

    n = box_name(wd, b)
    n isa Symbol || return false
    (startswith(String(n), "do_") || startswith(String(n), "do")) || return false
    nin(wd, b) == 0 || return false
    nout(wd, b) == 1 || return false
    return true
end

function vars_in_fragment(wd, frag_boxes::Vector{Int}, display_var)
    vars = Set{Symbol}()
    for b in frag_boxes
        for w in in_wires(wd, b)
            v = port_value(wd, w.target)
            v isa Symbol && v in display_var && push!(vars, v)
        end
        for w in out_wires(wd, b)
            v = port_value(wd, w.source)
            v isa Symbol && v in display_var && push!(vars, v)
        end
    end
    return collect(vars)
end





