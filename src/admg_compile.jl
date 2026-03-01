using Catlab
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra

include("utils.jl")


struct ConfoundedModel
    directed::Vector{Pair{Symbol, Symbol}}
    latents::Dict{Symbol, Vector{Symbol}}
end

struct ADMGModel
    directed::Vector{Pair{Symbol, Symbol}}
    bidirected::Vector{Pair{Symbol, Symbol}}
end

Base.:(==)(a::ConfoundedModel, b::ConfoundedModel) =
    a.directed == b.directed && a.latents == b.latents

Base.:(==)(a::ADMGModel, b::ADMGModel) =
    a.directed == b.directed && a.bidirected == b.bidirected

# Define query structure
struct CounterfactualQuery
    world_name::Symbol
    interventions::Dict{Symbol, Symbol}
    observations::Dict{Symbol, Symbol}
    outputs::Vector{Symbol}
end

function rootify(model::ADMGModel; latent_prefix::Symbol=:R)
    latents = Dict{Symbol, Vector{Symbol}}()

    for (i, edge) in enumerate(model.bidirected)
        a, b = edge.first, edge.second
        a == b && error("rootify: bidirected edge $edge must connect two distinct variables")
        latents[Symbol(latent_prefix, i)] = [a, b]
    end

    return ConfoundedModel(copy(model.directed), latents)
end

function graph_b_to_scm(model::ConfoundedModel; outputs=nothing)
    # collect all nodes that appear anywhere
    nodes = Set{Symbol}()
    for (src, tgt) in model.directed
        push!(nodes, src, tgt)
    end
    for (root, children) in model.latents
        push!(nodes, root)
        union!(nodes, children)
    end
    all_nodes = sort(collect(nodes))

    # choose output variables
    latent_roots = Set(keys(model.latents))
    if outputs === nothing
        output_vars = sort([v for v in all_nodes if !(v in latent_roots)])
    else
        output_vars = sort(collect(outputs))
    end

    wd = WiringDiagram([], output_vars)
    out_port_map = Dict(v => i for (i, v) in enumerate(output_vars))
    node_to_fbox = Dict{Symbol, Int}()

    sorted_latent_roots = sort(collect(keys(model.latents)))

    # create nodes + connect to external outputs if requested
    for v in all_nodes
        parents = [src for (src, tgt) in model.directed if tgt == v]
        for root in sorted_latent_roots
            if v in model.latents[root] && root != v
                push!(parents, root)
            end
        end

        u_var_name = Symbol("U", v)
        f_inputs = Vector{Any}(vcat(parents, [u_var_name]))

        f_box_id = add_box!(wd, Box(Symbol("f_", v), f_inputs, [v]))
        node_to_fbox[v] = f_box_id

        u_box_id = add_box!(wd, Box(Symbol("PU_", v), [], [u_var_name]))
        add_wire!(wd, Port(u_box_id, OutputPort, 1) => Port(f_box_id, InputPort, length(f_inputs)))

        if haskey(out_port_map, v)
            out_idx = out_port_map[v]
            add_wire!(wd, Port(f_box_id, OutputPort, 1) => Port(output_id(wd), InputPort, out_idx))
        end
    end

    # connect f -> f
    for v in all_nodes
        f_child_id = node_to_fbox[v]
        parents = [src for (src, tgt) in model.directed if tgt == v]
        for root in sorted_latent_roots
            if v in model.latents[root] && root != v
                push!(parents, root)
            end
        end
        for (p_idx, parent) in enumerate(parents)
            add_wire!(wd, Port(node_to_fbox[parent], OutputPort, 1) => Port(f_child_id, InputPort, p_idx))
        end
    end

    return wd
end

function graph_b_to_scm(model::ADMGModel; outputs=nothing, latent_prefix::Symbol=:R)
    return graph_b_to_scm(rootify(model; latent_prefix=latent_prefix); outputs=outputs)
end



function build_multiverse(base_wd::WiringDiagram, queries::Vector{CounterfactualQuery})
    
    # Calculate the total output ports
    total_outputs = Symbol[]
    for q in queries; append!(total_outputs, q.outputs); end
    combined = WiringDiagram([], total_outputs)
    
    # Scan the reference image and create an index
    base_pu_ids = Int[]
    base_f_ids = Int[]
    base_id_to_var = Dict{Int, Symbol}() # ID -> :X
    
    for bid in box_ids(base_wd)
        if bid == input_id(base_wd) || bid == output_id(base_wd); continue; end
        type, var = get_node_info(box(base_wd, bid).value)
        
        if type == :Exogenous
            push!(base_pu_ids, bid)
        elseif type == :Mechanism
            push!(base_f_ids, bid)
            base_id_to_var[bid] = var
        end
    end

    # Establish a Global Shared Map
    # Map: Base_PU_ID -> New_Shared_ID
    global_u_map = Dict{Int, Int}()
    for bid in base_pu_ids
        global_u_map[bid] = add_box!(combined, box(base_wd, bid))
    end
    
    # Local World Maps
    # Map: World_Index -> Base_f_ID -> {source: New_ID, sink: New_ID}
    world_maps = Dict{Int, Dict{Int, NamedTuple{(:source, :sink), Tuple{Int, Int}}}}()
    
    for (w_idx, query) in enumerate(queries)
        world_maps[w_idx] = Dict()
        
        for bid in base_f_ids
            b_orig = box(base_wd, bid)
            var_name = base_id_to_var[bid]
            b_outs = output_ports(b_orig)
            
            if haskey(query.interventions, var_name)
                
                # Sink
                sink_id = add_box!(combined, b_orig) 
                
                # Discard
                disc_id = add_box!(combined, Box(Symbol("discard", var_name), b_outs, []))
                for i in 1:length(b_outs)
                    add_wire!(combined, Port(sink_id, OutputPort, i) => Port(disc_id, InputPort, i))
                end
                
                # Source
                do_val = query.interventions[var_name]               # :x
                do_label = Symbol("do_", var_name, "=", do_val)      # :do_X=x
                source_id = add_box!(combined, Box(do_label, [], b_outs))
                
                world_maps[w_idx][bid] = (source = source_id, sink = sink_id)
            else
                # Source = Sink
                new_id = add_box!(combined, b_orig)
                world_maps[w_idx][bid] = (source = new_id, sink = new_id)
            end
            
            # Obs
            if haskey(query.observations, var_name)
                src_id = world_maps[w_idx][bid].source
                obs_val = query.observations[var_name]                 # e.g. :x0 or :z
                obs_label = Symbol("obs_", var_name, "=", obs_val)     # :obs_X=x0
                obs_id = add_box!(combined, Box(obs_label, b_outs, []))
                for i in 1:length(b_outs)
                    add_wire!(combined, Port(src_id, OutputPort, i) => Port(obs_id, InputPort, i))
                end
            end
        end
    end
    
    # Context-Aware Wiring
    current_out_idx = 1
    
    for (w_idx, query) in enumerate(queries)
        local_map = world_maps[w_idx]
        
        # Reproduce Connection
        for wire in wires(base_wd)
            src_base = wire.source.box
            tgt_base = wire.target.box
            
            if src_base == input_id(base_wd) || src_base == output_id(base_wd) ||
               tgt_base == input_id(base_wd) || tgt_base == output_id(base_wd); continue; end
            
            # find New Source 
            new_src_id = -1
            if haskey(global_u_map, src_base)
                new_src_id = global_u_map[src_base]
            elseif haskey(local_map, src_base)
                new_src_id = local_map[src_base].source
            end
            
            # find New Target
            new_tgt_id = -1
            if haskey(local_map, tgt_base)
                new_tgt_id = local_map[tgt_base].sink
            end
            
            # Connection
            if new_src_id != -1 && new_tgt_id != -1
                add_wire!(combined, Port(new_src_id, OutputPort, wire.source.port) => 
                                    Port(new_tgt_id, InputPort, wire.target.port))
            end
        end
        
        # connect to output
        for req_out in query.outputs
            # check ID 
            for bid in base_f_ids
                if base_id_to_var[bid] == req_out
                    real_src = local_map[bid].source
                    add_wire!(combined, Port(real_src, OutputPort, 1) => 
                                        Port(output_id(combined), InputPort, current_out_idx))
                    current_out_idx += 1
                    break
                end
            end
        end
    end
    
    return combined
end
