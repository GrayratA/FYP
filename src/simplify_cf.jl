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
function sharp_state_name_from_effects(wd::WiringDiagram, ws, var_type::Symbol)
    effect_names = Symbol[]
    for w in ws
        is_sharp_effect(wd, w.target.box) || continue
        n = safe_box_name(wd, w.target.box)
        n isa Symbol && push!(effect_names, n)
    end

    unique_names = unique(effect_names)
    if length(unique_names) == 1
        only_name = only(unique_names)
        m = match(r"^obs_(.*)=(.*)$", String(only_name))
        if m !== nothing
            V, v = m.captures
            V == String(var_type) && return Symbol("do_", V, "=", v)
        end
    end

    return canonical_do_name(var_type)
end

function separate_sharp_effect_fanout_once!(wd::WiringDiagram; source_boxes=nothing)
    changed = false

    boxes =
        source_boxes === nothing ?
        collect(1:WD.nboxes(wd)) :
        [b for b in unique(source_boxes) if 1 <= b <= WD.nboxes(wd)]

    for b in boxes
        for j in 1:nout(wd, b)
            ws = copy(out_wires(wd, b, j))
            length(ws) > 1 || continue

            sharp_ws = [w for w in ws if is_sharp_effect(wd, w.target.box)]
            isempty(sharp_ws) && continue

            var_type = port_value(wd, first(sharp_ws).target)
            nonsharp_ws = [
                w for w in ws
                if !(is_sharp_effect(wd, w.target.box) && port_value(wd, w.target) == var_type)
            ]
            isempty(nonsharp_ws) && continue

            do_name = sharp_state_name_from_effects(wd, sharp_ws, var_type)
            do_box_id = WD.add_box!(wd, Box(do_name, Any[], Any[var_type]))

            for w in nonsharp_ws
                tgt_box = w.target.box
                tgt_in_idx = w.target.port
                WD.rem_wire!(wd, w)
                WD.add_wire!(wd, (do_box_id, 1) => (tgt_box, tgt_in_idx))
            end

            changed = true
        end
    end

    return changed
end

function separate_sharp_effect_copy_maps!(wd::WiringDiagram; source_boxes=nothing)
    changed = true
    while changed
        changed = separate_sharp_effect_fanout_once!(wd; source_boxes=source_boxes)
    end
    remove_lonely_boxes!(wd)
    return wd
end

# Split do-box fanout for visualization:
# if one do_* box feeds multiple targets, clone it so each clone feeds exactly one target.
# This preserves semantics and removes junction copy nodes caused by do-box fanout.
function split_do_fanout!(wd::WiringDiagram)
    changed = true
    while changed
        changed = false
        # Snapshot current ids because we may add boxes during iteration.
        for b in collect(WD.box_ids(wd))
            b <= 0 && continue
            is_sharp_state(wd, b) || continue
            nout(wd, b) == 1 || continue

            ws = copy(out_wires(wd, b, 1))
            length(ws) > 1 || continue

            # Keep the first outgoing edge on original box; clone for the rest.
            for w in ws[2:end]
                tgt_box = w.target.box
                tgt_port = w.target.port

                do_name = safe_box_name(wd, b)
                out_types = copy(output_ports(wd, b))
                do_clone = WD.add_box!(wd, Box(do_name, Any[], out_types))

                WD.rem_wire!(wd, w)
                WD.add_wire!(wd, (do_clone, 1) => (tgt_box, tgt_port))
            end

            changed = true
        end
    end

    remove_lonely_boxes!(wd)
    return wd
end

function merge_identical_deterministic_boxes(wd::WiringDiagram)
    flag = true
    while flag
        flag = false

        # Eq.(6)-style separation step:
        # whenever a sharp effect is connected to a copy map (fanout), separate
        # the involved wires before attempting deterministic-box merging.
        before_sep = length(box_ids(wd)) + length(wires(wd))
        separate_sharp_effect_copy_maps!(wd)
        after_sep = length(box_ids(wd)) + length(wires(wd))
        flag = flag || (after_sep != before_sep)

        for (name, boxes) in duplicate_mechanism_groups(wd)
            # group boxes by input source set
            groups = Dict{Any, Vector{Int}}()
            for b in boxes
                S = merge_signature(wd, b)
                push!(get!(groups, S, Int[]), b)
            end
            filter!(kv -> length(kv[2]) > 1, groups)  # remove groups with only one box
            if isempty(groups)
                continue
            end
            # merge boxes in each group
            merged_keeps = Int[]
            for (S, bs) in groups
                sorted_bs = sort(bs)
                b_keep = first(sorted_bs)
                merged_here = false
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
                    merged_here = true
                end
                merged_here && push!(merged_keeps, b_keep)
            end

            if !isempty(merged_keeps)
                before = length(box_ids(wd)) + length(wires(wd))
                separate_sharp_effect_copy_maps!(wd; source_boxes=merged_keeps)
                after = length(box_ids(wd)) + length(wires(wd))
                flag = flag || (after != before)
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

            state_ports = state_input_ports(wd, b)
            isempty(state_ports) && continue

            # Paper step 4 only applies when λ_X is fed directly into f_X,
            # not when the same state is still shared across multiple worlds.
            lambda_is_directly_fed(wd, b) || continue

            state_port_set = Set(state_ports)
            kept_input_indices = [j for j in 1:nin_b if !(j in state_port_set)]

            old_in_types  = input_ports(wd, b)
            old_out_types = output_ports(wd, b)

            new_in_types = Any[ old_in_types[j] for j in kept_input_indices ]
            new_out_types = copy(old_out_types)
            state_syms = Symbol[]
            for j in state_ports
                for w in in_wires(wd, b, j)
                    v = port_value(wd, w.source)
                    v isa Symbol && push!(state_syms, v)
                end
            end
            isempty(state_syms) && continue
            unique_state_syms = unique(state_syms)
            length(unique_state_syms) == 1 || continue

            mech_var = Symbol(replace(String(n), "f_" => ""))
            cname = canonical_c_name(mech_var)

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
