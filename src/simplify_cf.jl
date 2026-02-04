const WD = Catlab.WiringDiagrams

using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics: Graphviz, to_graphviz

include("utils.jl")

# remove all discard branches and prune dead boxes
function drop_discard_branches!(wd::WiringDiagram)
    changed = true
    while changed
        changed = false
        discards = [b for b in 1:WD.nboxes(wd) if is_discard_box(wd, b)]
        if !isempty(discards)
            for b in reverse(discards)
                for w in reverse(in_wires(wd, b))
                    WD.rem_wire!(wd, w)
                end
                WD.rem_box!(wd, b)
            end
            changed = true
        end
        
        pruned = true
        while pruned
            pruned = false
            for b in reverse(1:WD.nboxes(wd))
                if nout(wd, b) > 0 && no_consumers(wd, b)
                    for w in reverse(in_wires(wd, b))
                        WD.rem_wire!(wd, w)
                    end
                    WD.rem_box!(wd, b)
                    pruned = true
                    changed = true
                end
            end
        end
        
    end
    return wd
end

# merge identical deterministic boxes
function merge_identical_deterministic_boxes(wd::WiringDiagram)
    flag = true
    while flag
        flag = false
        same_name_boxes = find_same_name_boxes(wd)
        for (name, boxes) in same_name_boxes
            # group boxes by input source set
            groups = Dict{Set{Symbol}, Vector{Int}}()
            for b in boxes
                S = input_sources_set(wd, b)
                push!(get!(groups, S, Int[]), b)
            end
            filter!(kv -> length(kv[2]) > 1, groups)  # remove groups with only one box
            # println(groups)
            if isempty(groups)
                flag = false
                continue
            end
            sub_groups = Set{Int}()
            # merge boxes in each group
            for (S, bs) in groups
                sorted_bs = sort(bs)
                b_keep = first(sorted_bs)
                push!(sub_groups, b_keep)
                for dup in Iterators.reverse(sorted_bs[2:end])
                    for w in out_wires(wd, dup)
                        src_out_idx = w.source.port
                        tgt_box = w.target.box
                        tgt_in_idx = w.target.port

                        pair = (b_keep, src_out_idx) => (tgt_box, tgt_in_idx)
                        if !has_wire(wd, pair)
                            WD.add_wire!(wd, pair)
                        end
                        WD.rem_wire!(wd, w)
                    end
                    WD.rem_box!(wd, dup)
                    flag = true
                end
            end
            # println(1)
            for i in sub_groups
                sharpEffectGroup = Set{Symbol}()
                if length(output_ports(wd, i)) == 1
                    for j in 1:length(output_ports(wd, i))
                        ws = out_wires(wd, i, j)
                        if length(ws) > 1
                            for w in ws
                                tgt_box = w.target.box
                                if is_sharp_effect(wd, tgt_box)
                                    push!(sharpEffectGroup, port_value(wd, w.target))
                                end
                            end
                        end
                        if length(sharpEffectGroup) > 0
                            var_type = first(collect(sharpEffectGroup))
                            for w in ws
                                if port_value(wd, w.target) != var_type || !is_sharp_effect(wd, w.target.box)
                                    tgt_box = w.target.box
                                    tgt_in_idx = w.target.port
                                    WD.rem_wire!(wd, w)
                                    do_box = Box(Symbol("do" * String(var_type)), Any[], Any[var_type])
                                    do_box_id = WD.add_box!(wd, do_box)
                                    WD.add_wire!(wd, (do_box_id, 1) => (w.target.box, tgt_in_idx))
                                else
                                    tgt_box = w.target.box
                                    tgt_in_idx = w.target.port
                                    src_out_idx = w.source.port
                                    WD.rem_wire!(wd, w)
                                    WD.add_wire!(wd, (i, src_out_idx) => (tgt_box, tgt_in_idx))
                                end
                            end    
                        end
                    end
                end
                flag = true
            end
        end
        
    end
    
    remove_lonely_boxes!(wd)
    return wd

end

function replace_fx_with_cx!(wd::WiringDiagram)
    while true
        changed = false

        # use the topological order
        order = topological_order(wd)
 
        for b in order
            n = box_name(wd, b)
            # search the box startswith "f_"
            n isa Symbol && startswith(String(n), "f_") || continue

            nin_b = nin(wd, b)
            nin_b == 0 && continue

            kept_input_indices = Int[]
            has_state_input = false

            for j in 1:nin_b
                ws = in_wires(wd, b, j)

                #If any wire of this port comes from the state box, it is considered a "λ input"
                if any(is_state_box(wd, w.source.box) for w in ws)
                    has_state_input = true
                else
                    push!(kept_input_indices, j)
                end
            end

            # If this f_box has no λ-state input, it is not the target of step 4.
            has_state_input || continue

            old_in_types  = input_ports(wd, b)
            old_out_types = output_ports(wd, b)

            new_in_types = Any[ old_in_types[j] for j in kept_input_indices ]
            new_out_types = copy(old_out_types)
            state_syms = Symbol[]
            for j in 1:nin_b
                for w in in_wires(wd, b, j)
                    if is_state_box(wd, w.source.box)
                        v = port_value(wd, w.source)
                        v isa Symbol && push!(state_syms, v)
                    end
                end
            end
            isempty(state_syms) && continue
            v = first(state_syms)
            cname = Symbol("c" * String(v))

            # create new "c_box"
            c_box = Box(cname, new_in_types, new_out_types)
            c_id  = WD.add_box!(wd, c_box)


            # Change the wire originally connected to the f_box to connect to the c_box
            for (new_idx, j) in enumerate(kept_input_indices)
                for w in copy(in_wires(wd, b, j))
                    src = w.source
                    WD.rem_wire!(wd, w)
                    WD.add_wire!(wd, (src.box, src.port) => (c_id, new_idx))
                end
            end

            for j in 1:nin_b
                j in kept_input_indices && continue
                for w in copy(in_wires(wd, b, j))
                    WD.rem_wire!(wd, w)
                end
            end

            for ow in copy(out_wires(wd, b))
                src_port_idx = ow.source.port
                tgt_box  = ow.target.box
                tgt_port = ow.target.port
                WD.rem_wire!(wd, ow)
                WD.add_wire!(wd, (c_id, src_port_idx) => (tgt_box, tgt_port))
            end

            # remove the f_box
            WD.rem_box!(wd, b)

            changed = true
            break   # This round only processes one f_box, and then recalculates the topological order.
        end
        changed || break
    end
    remove_lonely_boxes!(wd)
    return wd
end