const WD = Catlab.WiringDiagrams

using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics: Graphviz, to_graphviz
# Splitting R-fragment

function find_r_fragments(wd::WiringDiagram)
    boxes = WD.box_ids(wd)

    # All c_* mechanism boxes
    mech_boxes = [b for b in boxes if is_c_mechanism(wd, b)]

    # cX,cY,cZ...（except cR*）
    obs_boxes  = [b for b in mech_boxes if !is_root_mechanism(wd, b)]
    # cR,cR1,cR2...
    root_boxes = [b for b in mech_boxes if  is_root_mechanism(wd, b)]

    obs_set  = Set(obs_boxes)
    root_set = Set(root_boxes)

    # Only record adjacency relationships connected through latent roots
    adj = Dict{Int,Set{Int}}()
    for b in mech_boxes
        adj[b] = Set{Int}()
    end

    # Only root -> observable (a line with object type R) counts as an 'ambiguous connection'
    for r in root_boxes
        for j in 1:nout(wd, r)
            for w in out_wires(wd, r, j)
                tgt = w.target.box
                if tgt in obs_set
                    push!(adj[r], tgt)
                    push!(adj[tgt], r)
                end
            end
        end
    end

    # Run connected components on this graph containing only (root, observable)
    visited   = Set{Int}()
    fragments = Vector{Vector{Int}}()

    for b in mech_boxes
        b in visited && continue
        comp  = Int[]
        stack = [b]
        while !isempty(stack)
            v = pop!(stack)
            v in visited && continue
            push!(visited, v)
            push!(comp, v)
            for nb in adj[v]
                nb in visited || push!(stack, nb)
            end
        end
        push!(fragments, comp)
    end

    # Arrange the boxes inside each fragment in topological order
    topo = topological_order(wd)
    pos  = Dict(b => i for (i,b) in enumerate(topo))
    for frag in fragments
        sort!(frag, by = b -> get(pos, b, typemax(Int)))
    end

    return fragments
end

function show_r_fragments(wd::WiringDiagram)
    frags = find_r_fragments(wd)
    println("Found $(length(frags)) R-fragments:")
    for (i, frag) in enumerate(frags)
        names = [string(safe_box_name(wd, b)) for b in frag]
        println("  Fragment $i: ", join(names, ", "))
    end
    return frags
end


function handle_case1_for_var!(
    wd,
    frag_boxes::Vector{Int},
    var_type::Symbol
)
    # Find all mechanisms cX that generate var_type within the fragment
    cX_boxes = [b for b in frag_boxes if is_cX_box(wd, b, var_type)]
    isempty(cX_boxes) && return :no_cX

    cX = first(cX_boxes)

    # Collect all X-type wires within the fragment
    X_wires = Wire[]
    for b in frag_boxes
        # in port
        for w in in_wires(wd, b)
            if port_value(wd, w.target) == var_type
                push!(X_wires, w)
            end
        end
        # out port
        for w in out_wires(wd, b)
            if port_value(wd, w.source) == var_type
                push!(X_wires, w)
            end
        end
    end

    # Find the X-type output line going out from cX
    cX_out_X_wires = [w for w in out_wires(wd, cX) if port_value(wd, w.source) == var_type]

    # Are any of these output lines connected to the sharp effect?
    has_effect = any(cX_out_X_wires) do w
        tgt = w.target
        tgt_box = tgt.box
        # Connecting to the output (-1) does not count as a sharp effect
        tgt_box < 0 && return false
        is_sharp_effect(wd, tgt_box) && port_value(wd, tgt) == var_type
    end

    #   Case 1.a
    if !has_effect
        # Output of cX has no sharp effect
        # Check: whether there is only this one/few X-shaped lines in the fragment (that is, these outputs of cX)
        if length(X_wires) == length(cX_out_X_wires) &&
           all(w -> w in X_wires, cX_out_X_wires)
            return :case1a_do_nothing
        end
        # else FAIL
    end

    #   Case 1.b 
    if has_effect

        state_sources = Int[]  # record sharp state boxs' id

        for b in frag_boxes
            for w in in_wires(wd, b)
                if port_value(wd, w.target) == var_type
                    src_box = w.source.box

                    # If the output of cX comes from cX itself, skip
                    if src_box == cX
                        continue
                    end

                    # If the source is a sharp state, then make a note of it
                    if is_sharp_state(wd, src_box)
                        push!(state_sources, src_box)
                    else
                        # fail
                        error("FAIL (case1, var=$(var_type)): X-input not from cX or sharp state")
                    end
                end
            end
        end

        unique_states = unique(state_sources)

        if isempty(unique_states)
            # All X usage comes directly from cX (add an effect), and can be seen as a degenerate case of 1.b
            return :case1b_detected
        end

        # Require "value x consistency" for all sharp states
        names = Set{Symbol}()
        for s in unique_states
            n = box_name(wd, s)
            n isa Symbol || error("FAIL: sharp state without Symbol name")
            push!(names, n)
        end

        if length(names) == 1
            # rewiting

            state_box = unique_states[1]  # only one sharp state box

            # Find all within the fragment:
            # 'Exit from this sharp state and feed the X-type input line to the box inside the fragment'
            targets = Tuple{Int,Int}[]  # (tgt_box, tgt_in_port)
            wires_to_remove = Wire[]

            for b in frag_boxes
                for w in in_wires(wd, b)
                    if port_value(wd, w.target) == var_type &&
                       w.source.box == state_box
                        push!(targets, (w.target.box, w.target.port))
                        push!(wires_to_remove, w)
                    end
                end
            end

            # Delete these wires coming out of the sharp state
            for w in wires_to_remove
                WD.rem_wire!(wd, w)
            end

            # Delete this sharp state box itself (doX is no longer needed in this fragment)
            WD.rem_box!(wd, state_box)

            # Find the X-type output port index of cX
            x_out_ports = Int[]
            for w in out_wires(wd, cX)
                if port_value(wd, w.source) == var_type
                    push!(x_out_ports, w.source.port)
                end
            end
            var_out_port = isempty(x_out_ports) ? 1 : first(x_out_ports)

            # Fan out the X of cX to all the previous targets (where it was originally fed by state x)
            for (tgt_box, tgt_in) in targets
                pair = (cX, var_out_port) => (tgt_box, tgt_in)
                has_wire(wd, pair) || WD.add_wire!(wd, pair)
            end

            return :case1b_rewritten
        else
            error("FAIL (case1, var=$(var_type)): multiple different sharp states for X in fragment")
        end
    end

    # else FAIL
    error("FAIL (case1, var=$(var_type)): does not satisfy 1.a or 1.b conditions")
end

"""
Case 2: cX is not in the fragment.

- If all the 'in-fragment lines' of type X within the fragment come from the same cX (in another fragment),
and the only thing in between is fan-out (multiple wires) → 2.a: do nothing;

- If all the 'in-fragment lines' of type X within the fragment come from sharp states (e.g., :doX),
  and the names of these sharp states are the same → merge the multiple states into one, then fan-out;

- Otherwise → 2.c: FAIL.
"""
function handle_case2_for_var!(
    wd,
    frag_boxes::Vector{Int},
    var_type::Symbol
)
    #  If the fragment itself has cX, that’s case1, not case2, skip it directly
    if any(b -> is_cX_box(wd, b, var_type), frag_boxes)
        return :skip_case1
    end

    # Collect all X-type wires entering the fragment
    X_input_wires = Wire[]
    for b in frag_boxes
        for w in in_wires(wd, b)
            tgt = w.target
            src = w.source
            if !(src.box in frag_boxes) && port_value(wd, tgt) == var_type
                push!(X_input_wires, w)
            end
        end
    end

    isempty(X_input_wires) && return :skip_no_X_input

    # check the source
    src_boxes = [w.source.box for w in X_input_wires]
    unique_srcs = unique(src_boxes)

    # All sources are the same cX (in another fragment)
    if length(unique_srcs) == 1
        s = unique_srcs[1]
        if is_cX_box(wd, s, var_type)
            # do nothing
            return :case2a_do_nothing
        end
    end

    # All sources are sharp state (doX/doZ...), and it is the 'same x'
    all_sharp = all(s -> is_sharp_state(wd, s), unique_srcs)
    if all_sharp
        # check if they are the same state
        names = Set{Symbol}()
        for s in unique_srcs
            n = box_name(wd, s)
            n isa Symbol || error("FAIL: sharp state without Symbol name")
            push!(names, n)
        end

        if length(names) == 1
            # rewirte
            # Merge multiple states into one (keeping the first), and feed all X_input_wires from this state.
            keep = unique_srcs[1]
            keep_name = box_name(wd, keep)

            # Delete all X_input_wires
            for w in X_input_wires
                WD.rem_wire!(wd, w)
            end

            # Fan out FIRST (while box ids are still valid)
            for w in X_input_wires
                tgt = w.target
                WD.add_wire!(wd, (keep, 1) => (tgt.box, tgt.port))
            end

            # remove extra sharp state boxes (delete in descending id order is safer)
            extras = sort([s for s in unique_srcs if s != keep], rev=true)
            for s in extras
                WD.rem_box!(wd, s)
            end


            return :case2b_rewritten
        end
    end

    # else fail
    error("FAIL (case2, var=$(var_type)): X-input wires come from mixed or invalid sources")
end


struct WireInfo
    var::Symbol
    src_box::Int
    src_port::Int
    tgt_box::Int
    tgt_port::Int
end

"""
Generate the name of a probability box based on 4 sets of variable types:
P(C_l, F_l^out ; do(D_l), do(F_l^in))

For example:
Cl_vars = [:Y]
Fout_vars = [:Z, :W]
Dl_vars = [:X]
Fin_vars = [:D]

This gives Symbol("P(Y,Z,W ; do(X),do(D))").

frag_id is optional and can be used to mark which fragment it is, e.g., P_1(...)
"""
function make_prob_box_name(
    Cl_vars::Vector{Symbol},
    Fout_vars::Vector{Symbol},
    Dl_vars::Vector{Symbol},
    Fin_vars::Vector{Symbol};
    frag_id::Union{Nothing,Int}=nothing,
)
    # sort
    Cls   = sort!(copy(Cl_vars))
    Fouts = sort!(copy(Fout_vars))
    Dls   = sort!(copy(Dl_vars))
    Fins  = sort!(copy(Fin_vars))

    # Convert a set of variables into a string in the form "X,Y"
    set_str(vs::Vector{Symbol}) =
        isempty(vs) ? "∅" : join(string.(vs), ",")

    C_str    = set_str(Cls)
    Fout_str = set_str(Fouts)
    D_str    = set_str(Dls)
    Fin_str  = set_str(Fins)

    left_parts  = String[]   # left: C_l, F_l^out
    right_parts = String[]   # right: do(D_l), do(F_l^in)

    if !isempty(Cls)
        push!(left_parts, C_str)
    end
    if !isempty(Fouts)
        push!(left_parts, Fout_str)
    end

    if !isempty(Dls)
        push!(right_parts, "do(" * D_str * ")")
    end
    if !isempty(Fins)
        push!(right_parts, "do(" * Fin_str * ")")
    end

    left  = isempty(left_parts)  ? "∅" : join(left_parts, ",")
    right = isempty(right_parts) ? "∅" : join(right_parts, ",")

    label = "P(" * left * " ; " * right * ")"

    if frag_id !== nothing
        label = "P_$frag_id" * label[2:end]
    end

    return Symbol(label)
end

"""
Step 4.2:
Replace the rewritten R-fragment `frag_boxes` with a probability box P(C_l, F_l^out ; do(D_l), do(F_l^in)).

- Ports are organized by "variable sets":
Input = all variable types of D_l + all variable types of F_l^in
Output = all variable types of C_l + all variable types of F_l^out
- Each variable type occupies only one port in each group,
  multiple wires connect to this port through fan-out.

Parameters:
wd : WiringDiagram
frag_boxes : List of box ids in the current fragment
frag_id : Number of the fragment (only affects naming)

Returns:
The id of the newly created probability box
"""
function replace_fragment_with_prob_box_grouped!(
    wd,
    frag_boxes::Vector{Int},
    frag_id::Int,
)
    # Collect the boundary wires of the fragment
    incoming = Wire[]   # outside -> fragment
    outgoing = Wire[]   # fragment -> outside

    for b in frag_boxes
        for w in in_wires(wd, b)
            if !(w.source.box in frag_boxes)
                push!(incoming, w)
            end
        end
        for w in out_wires(wd, b)
            if !(w.target.box in frag_boxes)
                push!(outgoing, w)
            end
        end
    end

    # Group into 4 Dicts by Variable Type
    Dl_map   = Dict{Symbol, Vector{WireInfo}}()  # D_l: input from sharp state
    Fin_map  = Dict{Symbol, Vector{WireInfo}}()  # F_l^in: other input
    Cl_map   = Dict{Symbol, Vector{WireInfo}}()  # C_l: output to sharp effect 
    Fout_map = Dict{Symbol, Vector{WireInfo}}()  # F_l^out: other output

    # Input edge: sharp state vs non-state
    for w in incoming
        vt = port_value(wd, w.target)
        sb, sp = w.source.box, w.source.port
        tb, tp = w.target.box, w.target.port
        info = WireInfo(vt, sb, sp, tb, tp)

        if is_sharp_state(wd, sb)
            push!(get!(Dl_map, vt, WireInfo[]), info)
        else
            push!(get!(Fin_map, vt, WireInfo[]), info)
        end
    end

    # Output edge: sharp effect vs non-effect
    for w in outgoing
        vt = port_value(wd, w.source)
        sb, sp = w.source.box, w.source.port
        tb, tp = w.target.box, w.target.port
        info = WireInfo(vt, sb, sp, tb, tp)

        if is_sharp_effect(wd, tb)
            push!(get!(Cl_map, vt, WireInfo[]), info)
        else
            push!(get!(Fout_map, vt, WireInfo[]), info)
        end
    end

    # delete all
    for b in frag_boxes
        for w in in_wires(wd, b)
            WD.rem_wire!(wd, w)
        end
        for w in out_wires(wd, b)
            WD.rem_wire!(wd, w)
        end
    end
    # for b in frag_boxes
    #     WD.rem_box!(wd, b)
    # end

    Dl_vars   = sort(collect(keys(Dl_map)))    # D_l 
    Fin_vars  = sort(collect(keys(Fin_map)))   # F_l^in
    Cl_vars   = sort(collect(keys(Cl_map)))    # C_l
    Fout_vars = sort(collect(keys(Fout_map)))  # F_l^out

    in_types  = vcat(Dl_vars,  Fin_vars)
    out_types = vcat(Cl_vars, Fout_vars)

    # Generate name: P(C_l,F_l^out ; do(D_l),do(F_l^in))
    P_name = make_prob_box_name(Cl_vars, Fout_vars, Dl_vars, Fin_vars;
                                frag_id = frag_id)

    P_box = Box(P_name, Any[in_types...], Any[out_types...])
    P_id  = WD.add_box!(wd, P_box)

    Dl_index   = Dict{Symbol, Int}()
    Fin_index  = Dict{Symbol, Int}()
    Cl_index   = Dict{Symbol, Int}()
    Fout_index = Dict{Symbol, Int}()

    for (i, vt) in enumerate(Dl_vars)
        Dl_index[vt] = i
    end
    for (k, vt) in enumerate(Fin_vars)
        Fin_index[vt] = length(Dl_vars) + k
    end

    for (i, vt) in enumerate(Cl_vars)
        Cl_index[vt] = i
    end
    for (k, vt) in enumerate(Fout_vars)
        Fout_index[vt] = length(Cl_vars) + k
    end

    # Connect the external endpoint to the P box

    # D_l
    for (vt, infos) in Dl_map
        info = first(infos)
        in_idx = Dl_index[vt]
        WD.add_wire!(wd, (info.src_box, info.src_port) => (P_id, in_idx))
    end

    # F_l^in
    for (vt, infos) in Fin_map
        info = first(infos)
        in_idx = Fin_index[vt]
        WD.add_wire!(wd, (info.src_box, info.src_port) => (P_id, in_idx))
    end

    # C_l
    for (vt, infos) in Cl_map
        out_idx = Cl_index[vt]
        for info in infos
            WD.add_wire!(wd, (P_id, out_idx) => (info.tgt_box, info.tgt_port))
        end
    end

    # F_l^out
    for (vt, infos) in Fout_map
        out_idx = Fout_index[vt]
        for info in infos
            WD.add_wire!(wd, (P_id, out_idx) => (info.tgt_box, info.tgt_port))
        end
    end

    return P_id
end

# Step 4.1 only: rewrite within each R-fragment
# Returns r_frags so Step 4.2 can reuse the fragments computed on the ORIGINAL diagram.
function id_cf_step41!(wd, display_vars; r_frags=nothing, verbose::Bool=true)
    # R-fragments computed once on the original diagram (paper requirement)
    r_frags === nothing && (r_frags = find_r_fragments(wd))

    for (l, frag_boxes) in enumerate(r_frags)
        verbose && println("---- Step 4.1 on R-fragment #$l ----")

        frag_vars = vars_in_fragment(wd, frag_boxes, display_vars)

        for v in frag_vars
            try
                status1 = handle_case1_for_var!(wd, frag_boxes, v)

                if status1 == :no_cX
                    status2 = handle_case2_for_var!(wd, frag_boxes, v)
                    verbose && println("  Fragment $l, var $v, case2 = $status2")
                else
                    verbose && println("  Fragment $l, var $v, case1 = $status1")
                end
            catch e
                verbose && println("  FAIL in fragment $l, var $v: ", e)
                rethrow(e)
            end
        end
    end

    return r_frags
end


# Step 4.2 only: collapse each fragment into a probability box
function id_cf_step42!(wd, r_frags; verbose::Bool=true)
    for (l, frag_boxes) in enumerate(r_frags)
        verbose && println("---- Step 4.2: collapsing R-fragment #$l ----")
        replace_fragment_with_prob_box_grouped!(wd, frag_boxes, l)
    end

    remove_lonely_boxes!(wd)
    return wd
end

function id_cf_step3_check!(wd)
    remaining = unabsorbed_lambda_boxes(wd)
    isempty(remaining) && return wd

    names = [string(safe_box_name(wd, b)) for b in remaining]
    error("FAIL (step3): λ_X was not absorbed into c_X for " * join(names, ", "))
end


function id_cf_step4!(wd, display_vars; verbose::Bool=true)
    id_cf_step3_check!(wd)
    r_frags = id_cf_step41!(wd, display_vars; r_frags=nothing, verbose=verbose)
    id_cf_step42!(wd, r_frags; verbose=verbose)
    return wd
end

# function id_cf_step4!(wd, display_vars)

#     # R-fragments are computed once on the original diagram
#     r_frags = find_r_fragments(wd)

#     # Step 4.1 on ALL fragments
#     for (l, frag_boxes) in enumerate(r_frags)
#         println("---- Step 4.1 on R-fragment #$l ----")

#         # variables in this fragment ∩ display_vars
#         frag_vars = vars_in_fragment(wd, frag_boxes, display_vars)

#         for v in frag_vars
#             try
#                 # Case 1: c_X is inside the fragment
#                 status1 = handle_case1_for_var!(wd, frag_boxes, v)

#                 if status1 == :no_cX
#                     # Case 2: c_X is not inside the fragment
#                     status2 = handle_case2_for_var!(wd, frag_boxes, v)
#                     println("  Fragment $l, var $v, case2 = $status2")
#                 else
#                     println("  Fragment $l, var $v, case1 = $status1")
#                 end

#             catch e
#                 # corresponds to “Else output FAIL”
#                 println("  FAIL in fragment $l, var $v: ", e)
#                 rethrow(e)
#             end
#         end
#     end


#     # Step 4.2 — collapse each fragment into a prob box
#     for (l, frag_boxes) in enumerate(r_frags)
#         println("---- Step 4.2: collapsing R-fragment #$l ----")
#         replace_fragment_with_prob_box_grouped!(wd, frag_boxes, l)
#     end
    
#     remove_lonely_boxes!(wd)
#     return wd
# end

# ============================================================
# Step 5 (+ Step 4.5 rewrite_for_data): diagram -> LaTeX fraction
# Generic (no starred-style naming). Optional data-mode rewriting.
# ============================================================

const WD = Catlab.WiringDiagrams

# factor representation
struct ProbFactor
    left::Vector{String}    # vars on left of P(...)
    dovars::Vector{String}  # vars inside do(...)
end

# small expr AST
abstract type CFExpr end

struct CFAtom <: CFExpr
    s::String
    left::Vector{String}       # var-names (e.g. "X","Y")
    dovars::Vector{String}     # var-names inside do(...)
    left_disp::Vector{String}  # printed tokens (e.g. "x'","y")
    do_disp::Vector{String}    # printed tokens (e.g. "z","w","d")
end

struct CFProd <: CFExpr
    terms::Vector{CFExpr}
end

struct CFSum <: CFExpr
    vars::Vector{String}   # printed indices, e.g. ["w","d^*"]
    body::CFExpr
end

struct CFFrac <: CFExpr
    num::CFExpr
    den::CFExpr
end

# LaTeX printing
function cf_latex(e::CFExpr)
    if e isa CFAtom
        return (e::CFAtom).s
    elseif e isa CFProd
        terms = (e::CFProd).terms
        isempty(terms) && return "1"
        return join(cf_latex.(terms), "")
    elseif e isa CFSum
        ee = e::CFSum
        return "\\sum_{" * join(ee.vars, ",") * "} " * cf_latex(ee.body)
    elseif e isa CFFrac
        ee = e::CFFrac
        return "\\frac{" * cf_latex(ee.num) * "}{" * cf_latex(ee.den) * "}"
    else
        error("unknown expr")
    end
end

# parse naming conventions
function parse_prob_box_name(sym::Symbol)
    s = String(sym)
    startswith(s, "P") || return nothing
    s2 = replace(s, r"^P_\d+" => "P")     # drop fragment id
    m = match(r"P\((.*)\)", s2); m === nothing && return nothing
    inside = m.captures[1]               # e.g. "X,Y ; do(Z),do(W)"
    parts = split(inside, " ; ")
    length(parts) == 2 || return nothing

    left  = strip(parts[1])
    right = strip(parts[2])

    left_vars = left == "∅" ? String[] : split(replace(left, " " => ""), ",")
    do_vars = String[]      # NOTE: store var-names inside do(...)
    if right != "∅"
        compact = replace(right, " " => "")
        for mm in eachmatch(r"do\(([^)]*)\)", compact)
            body = mm.captures[1]
            isempty(body) && continue
            append!(do_vars, [v for v in split(body, ",") if !isempty(v)])
        end
    end
    return ProbFactor(left_vars, do_vars)
end

parse_obs(sym::Symbol) = begin
    s = String(sym)
    m = match(r"^obs_(.*)=(.*)$", s)
    m === nothing ? nothing : (m.captures[1], m.captures[2])  # ("X","x_hat")
end

parse_do(sym::Symbol) = begin
    s = String(sym)
    m = match(r"^do_(.*)=(.*)$", s)
    m === nothing ? nothing : (m.captures[1], m.captures[2])  # ("X","x")
end

# ---------- naming policy ----------
default_var_symbol(v::String) = lowercase(v)

# build atom
function factor_to_atom(
    f::ProbFactor,
    fixed_obs::Dict{String,String},
    fixed_do::Dict{String,String},
    var2sym::Dict{String,String};
    display_mode::Symbol = :legacy,
    star_left_vars::Set{String} = Set{String}()
)
    display_mode in (:legacy, :starred_style) ||
        error("factor_to_atom: unsupported display_mode=$display_mode")

    base_disp(V::String) = get(var2sym, V, default_var_symbol(V))
    star_disp(V::String) = base_disp(V) * "^*"

    left_disp = String[]
    for V in f.left
        if display_mode == :legacy
            push!(left_disp, get(fixed_obs, V, base_disp(V)))
        else
            if isempty(f.dovars) && length(f.left) == 1 && (V in star_left_vars)
                push!(left_disp, star_disp(V))
            else
                push!(left_disp, base_disp(V))
            end
        end
    end
    do_disp = String[]
    for V in f.dovars
        if display_mode == :legacy
            push!(do_disp, get(fixed_do, V, base_disp(V)))
        else
            push!(do_disp, base_disp(V))
        end
    end

    evt_s = join(left_disp, ",")
    s = isempty(f.dovars) ? "P(" * evt_s * ")" :
                            "P(" * evt_s * "|do(" * join(do_disp, ",") * "))"

    return CFAtom(s, copy(f.left), copy(f.dovars), left_disp, do_disp)
end

# collect all var-names from factors
function collect_all_vars(factors::Vector{ProbFactor})
    vars = Set{String}()
    for f in factors
        foreach(v -> push!(vars, v), f.left)
        foreach(v -> push!(vars, v), f.dovars)
    end
    return vars
end

"""
Simplify Σ_{sumvars} ∏ atoms:

- drop summation vars that don't appear anymore
- if v appears in exactly one atom, only in its LEFT list (not do list),
  marginalize it out: Σ_v P(v,rest | do(...)) -> P(rest | do(...))
- remove trivial factor "1"

IMPORTANT: sumvars are VAR-NAMES (e.g. "Y","W"), not printed tokens.
"""
function simplify_sum(prod::CFProd, sumvars::Vector{String})
    atoms = CFAtom[t for t in prod.terms if t isa CFAtom]

    appears_in_atom(v::String, a::CFAtom) = (v in a.left) || (v in a.dovars)

    changed = true
    while changed
        changed = false

        # drop unused sum vars
        used = Set{String}()
        for a in atoms
            foreach(v -> push!(used, v), a.left)
            foreach(v -> push!(used, v), a.dovars)
        end
        new_sumvars = [v for v in sumvars if v in used]
        if length(new_sumvars) != length(sumvars)
            sumvars = new_sumvars
            changed = true
        end

        # marginalize v if it appears in exactly one atom and only in left
        for v in copy(sumvars)
            hit = Int[]
            for (i,a) in enumerate(atoms)
                appears_in_atom(v, a) && push!(hit, i)
            end
            length(hit) == 1 || continue
            i = hit[1]
            a = atoms[i]
            (v in a.left) || continue
            (v in a.dovars) && continue

            idxs = findall(==(v), a.left)
            length(idxs) == 1 || continue
            j = idxs[1]

            new_left      = [a.left[k] for k in eachindex(a.left) if k != j]
            new_left_disp = [a.left_disp[k] for k in eachindex(a.left_disp) if k != j]

            new_evt = join(new_left_disp, ",")
            new_s = isempty(new_left_disp) ? "1" :
                    (isempty(a.dovars) ? "P($new_evt)" :
                     "P($new_evt|do(" * join(a.do_disp, ",") * "))")

            atoms[i] = CFAtom(new_s, new_left, a.dovars, new_left_disp, a.do_disp)

            sumvars = [u for u in sumvars if u != v]
            changed = true
        end

        # remove trivial "1"
        before = length(atoms)
        atoms = [a for a in atoms if a.s != "1"]
        if length(atoms) != before
            changed = true
        end
    end

    return CFProd(CFExpr[a for a in atoms]), sumvars
end

# drop-do rule helpers
function build_adj(directed_edges::Vector{Tuple{String,String}})
    adj = Dict{String,Vector{String}}()
    for (u,v) in directed_edges
        push!(get!(adj, u, String[]), v)
        get!(adj, v, String[])
    end
    return adj
end

function is_descendant(target::String, sources::Vector{String}, adj::Dict{String,Vector{String}})
    seen = Set{String}()
    stack = copy(sources)
    while !isempty(stack)
        u = pop!(stack)
        u in seen && continue
        push!(seen, u)
        for v in get(adj, u, String[])
            v == target && return true
            v in seen || push!(stack, v)
        end
    end
    return false
end

"""
Drop-do rule (safe):
If ALL left vars are NOT descendants of any do-var, then P(left | do(do_vars)) -> P(left).
"""
function rule_drop_do!(atoms::Vector{CFAtom}, directed_edges::Vector{Tuple{String,String}})
    isempty(directed_edges) && return atoms
    adj = build_adj(directed_edges)

    for i in eachindex(atoms)
        a = atoms[i]
        isempty(a.dovars) && continue

        all_non_desc = all(V -> !is_descendant(V, a.dovars, adj), a.left)
        all_non_desc || continue

        new_s = "P(" * join(a.left_disp, ",") * ")"
        atoms[i] = CFAtom(new_s, a.left, String[], a.left_disp, String[])
    end
    return atoms
end

"""
Eliminate normalized sums:
If v is a sum-var and only appears in a single atom whose left is exactly [v],
and v does not appear anywhere else, then Σ_v P(v|...) = 1 and we can drop both.
"""
function rule_eliminate_sum_normalized!(prod::CFProd, sumvars::Vector{String})
    atoms = CFAtom[t for t in prod.terms if t isa CFAtom]

    function var_appears_elsewhere(v::String, atoms::Vector{CFAtom}, skip::Int)
        for (i,a) in enumerate(atoms)
            i == skip && continue
            (v in a.left || v in a.dovars) && return true
        end
        return false
    end

    changed = true
    while changed
        changed = false
        for v in copy(sumvars)
            idx = nothing
            for (i,a) in enumerate(atoms)
                if length(a.left) == 1 && a.left[1] == v && !(v in a.dovars)
                    idx = i
                    break
                end
            end
            idx === nothing && continue

            if !var_appears_elsewhere(v, atoms, idx)
                deleteat!(atoms, idx)
                sumvars = [u for u in sumvars if u != v]
                changed = true
                break
            end
        end
    end

    return CFProd(CFExpr[a for a in atoms]), sumvars
end


_to_prod(e::CFExpr) = e isa CFProd ? (e::CFProd) : CFProd(CFExpr[e])

function find_atom(e::CFExpr, pred)::Union{Nothing,CFAtom}
    if e isa CFAtom
        return pred(e::CFAtom) ? (e::CFAtom) : nothing
    elseif e isa CFProd
        for t in (e::CFProd).terms
            a = find_atom(t, pred)
            a === nothing || return a
        end
        return nothing
    elseif e isa CFSum
        return find_atom((e::CFSum).body, pred)
    elseif e isa CFFrac
        # search both sides just in case
        f = e::CFFrac
        a = find_atom(f.num, pred)
        a === nothing ? find_atom(f.den, pred) : a
    else
        return nothing
    end
end

"""
rewrite_for_data(fr; data_mode=:none | :interventions)

- :none -> no rewrite
- :interventions -> introduce a fresh dummy variable (mix_var*) with printed token mix_sym (e.g. d^*),
  and rewrite denominator into a mixture form:
      den := Σ_{mix_sym} P(obs_token) P(z|do(mix_sym)) P(mix_sym)
  and multiply numerator by Σ_{mix_sym} P(mix_sym) (merging with existing Σ_w if present).

This is a presentation/canonicalization pass.
"""
function rewrite_for_data(
    fr::CFFrac;
    data_mode::Symbol = :none,
    mix_var::String = "D",
    mix_sym::String = "d^*",
    z_var::String = "Z",
    z_sym::String = "z",
    obs_token::Union{Nothing,String} = nothing,  # NEW: anchor token, e.g. "x'"
    mix_token::Union{Nothing,String} = nothing,  # NEW: printed token for original mix var, e.g. "d"
    anchor_is_mix::Bool = false,                 # NEW: skip explicit anchor factor when anchor var == mix var
)::CFFrac
    data_mode == :none && return fr
    data_mode == :interventions || error("rewrite_for_data: unknown data_mode=$data_mode")

    mix_star = mix_var * "*"   # internal var-name for the dummy

    # Find an atom of the form P(<tok>) with no do()
    function is_P_of_token(a::CFAtom, tok::String)
        isempty(a.dovars) || return false
        length(a.left_disp) == 1 && a.left_disp[1] == tok
    end

    tok = isnothing(obs_token) ? "x'" : obs_token

    Px = nothing
    if !anchor_is_mix
        # recursive search (works even if denominator is a sum/product nest)
        Px = find_atom(fr.den, a -> is_P_of_token(a, tok))

        if Px === nothing
            # do not fail: this pass is a canonicalization for "data mode"
            println("[rewrite_for_data] WARNING: could not find P($tok) in denominator; synthesizing it.")
            println("[rewrite_for_data] den = ", cf_latex(fr.den))
            Px = CFAtom("P($tok)", ["ANCHOR"], String[], [tok], String[])
        end
    end

    # atoms: P(mix*) and P(z|do(mix*))
    Pmix = CFAtom("P($mix_sym)", [mix_star], String[], [mix_sym], String[])
    Pz_do_mix = CFAtom("P($z_sym|do($mix_sym))",
                       [z_var], [mix_star],
                       [z_sym], [mix_sym])

    # numerator: replace plain P(mix_token) by P(mix*) if present; otherwise append P(mix*).
    # This avoids leaving a stale P(d) term while introducing the mixture index d^*.
    function replace_or_append_mix_atom(terms::Vector{CFExpr}, tok::String, mix_atom::CFAtom)
        replaced = false
        new_terms = CFExpr[]
        for t in terms
            if t isa CFAtom
                a = t::CFAtom
                if isempty(a.dovars) && length(a.left_disp) == 1 && a.left_disp[1] == tok
                    push!(new_terms, mix_atom)
                    replaced = true
                    continue
                end
            end
            push!(new_terms, t)
        end
        if !replaced
            push!(new_terms, mix_atom)
        end
        return new_terms
    end

    mix_tok = isnothing(mix_token) ? lowercase(mix_var) : mix_token

    # numerator: merge sums if numerator already has Σ_{...}
    old_num = fr.num
    if old_num isa CFSum
        s = old_num::CFSum
        inner_terms = _to_prod(s.body).terms
        num_terms = replace_or_append_mix_atom(inner_terms, mix_tok, Pmix)
        num_inner = CFProd(num_terms)
        new_sumvars = mix_sym in s.vars ? s.vars : vcat(s.vars, [mix_sym])
        num2 = CFSum(new_sumvars, num_inner)
    else
        num_terms = replace_or_append_mix_atom(_to_prod(old_num).terms, mix_tok, Pmix)
        num_inner = CFProd(num_terms)
        num2 = CFSum([mix_sym], num_inner)
    end

    function is_mix_related_atom(a::CFAtom, mv::String, mt::String)
        (mv in a.left) || (mv in a.dovars) || (mt in a.left_disp) || (mt in a.do_disp)
    end

    function star_mix_atom_for_den(a::CFAtom, mv::String, mt::String, mv_star::String, msym::String)
        new_left = [v == mv ? mv_star : v for v in a.left]
        new_left_disp = [d == mt ? msym : d for d in a.left_disp]

        mix_in_do = any(v -> v == mv, a.dovars) || any(d -> d == mt, a.do_disp)
        new_do = mix_in_do ? [mv_star] : [v == mv ? mv_star : v for v in a.dovars]
        new_do_disp = mix_in_do ? [msym] : [d == mt ? msym : d for d in a.do_disp]

        evt_s = join(new_left_disp, ",")
        s = isempty(new_do_disp) ? "P(" * evt_s * ")" :
            "P(" * evt_s * "|do(" * join(new_do_disp, ",") * "))"
        return CFAtom(s, new_left, new_do, new_left_disp, new_do_disp)
    end

    den_body = fr.den isa CFSum ? (fr.den::CFSum).body : fr.den
    den_terms_old = _to_prod(den_body).terms
    mix_terms = CFExpr[]
    for t in den_terms_old
        t isa CFAtom || continue
        a = t::CFAtom
        is_mix_related_atom(a, mix_var, mix_tok) || continue
        push!(mix_terms, star_mix_atom_for_den(a, mix_var, mix_tok, mix_star, mix_sym))
    end

    den_terms = CFExpr[]
    if !isempty(mix_terms)
        if !anchor_is_mix
            push!(den_terms, Px::CFAtom)
        end
        append!(den_terms, mix_terms)
        has_pmix = any(t -> begin
            t isa CFAtom || return false
            a = t::CFAtom
            length(a.left) == 1 && a.left[1] == mix_star && isempty(a.dovars)
        end, den_terms)
        has_pmix || push!(den_terms, Pmix)
    else
        # fallback canonical form when no mix-related denominator atom is detected
        if !anchor_is_mix
            push!(den_terms, Px::CFAtom)
        end
        push!(den_terms, Pz_do_mix, Pmix)
    end

    den_inner = CFProd(den_terms)
    den2 = CFSum([mix_sym], den_inner)

    return CFFrac(num2, den2)
end

_to_prod_expr(e::CFExpr) = e isa CFProd ? (e::CFProd) : CFProd(CFExpr[e])

function map_atoms(e::CFExpr, f)::CFExpr
    if e isa CFAtom
        return f(e::CFAtom)
    elseif e isa CFProd
        ee = e::CFProd
        new_terms = CFExpr[]
        for t in ee.terms
            mapped = map_atoms(t, f)
            if mapped isa CFAtom && (mapped::CFAtom).s == "1"
                continue
            end
            push!(new_terms, mapped)
        end
        return CFProd(new_terms)
    elseif e isa CFSum
        ee = e::CFSum
        return CFSum(copy(ee.vars), map_atoms(ee.body, f))
    elseif e isa CFFrac
        ee = e::CFFrac
        return CFFrac(map_atoms(ee.num, f), map_atoms(ee.den, f))
    end
    return e
end

function rebuild_atom(
    a::CFAtom;
    left_disp::Vector{String}=copy(a.left_disp),
    do_disp::Vector{String}=copy(a.do_disp),
)::CFAtom
    if length(a.left) == 1 && isempty(a.dovars)
        s = "P(" * join(left_disp, ",") * ")"
    else
        evt_s = join(left_disp, ",")
        s = isempty(do_disp) ? "P(" * evt_s * ")" :
            "P(" * evt_s * "|do(" * join(do_disp, ",") * "))"
    end
    return CFAtom(s, copy(a.left), copy(a.dovars), left_disp, do_disp)
end

function wrap_sumvar(expr::CFExpr, var::String)::CFExpr
    if expr isa CFSum
        ee = expr::CFSum
        return var in ee.vars ? expr : CFSum(vcat(ee.vars, [var]), ee.body)
    end
    return CFSum([var], expr)
end

function rewrite_conditional_queries(
    fr::CFFrac,
    factors::Vector{ProbFactor},
    fixed_obs::Dict{String,String},
    fixed_do::Dict{String,String},
    query_obs::Dict{String,String},
    query_do::Dict{String,String},
    output_vars::Vector{String},
    queries,
    var2sym::Dict{String,String},
)::CFFrac
    queries === nothing && return fr
    isempty(output_vars) && return fr

    factor_left = Set{String}()
    factor_do = Set{String}()
    for f in factors
        union!(factor_left, f.left)
        union!(factor_do, f.dovars)
    end

    base_disp(V::String) = get(var2sym, V, default_var_symbol(V))
    star_disp(V::String) = base_disp(V) * "^*"

    lower_to_known = Dict{String,Vector{String}}()
    for V in union(factor_left, factor_do)
        push!(get!(lower_to_known, lowercase(V), String[]), V)
    end

    canonical_var(V::String) = begin
        if (V in factor_left) || (V in factor_do)
            V
        else
            matches = get(lower_to_known, lowercase(V), String[])
            length(matches) == 1 ? only(matches) : V
        end
    end

    event_left_vars = Set{String}(canonical_var.(output_vars))
    for q in queries
        !isempty(getproperty(q, :interventions)) || continue
        for (V, _) in getproperty(q, :observations)
            push!(event_left_vars, canonical_var(string(V)))
        end
    end

    dual_fixed_vars = Set{String}(intersect(keys(query_obs), keys(query_do)))
    mix_vars = Set{String}()
    for V in keys(query_obs)
        (V in dual_fixed_vars) && continue
        (V in event_left_vars) && continue
        # Lift observed roots that are only used as conditioning/intervention context.
        if (V in factor_left) && any(f -> (V in f.left) && isempty(f.dovars), factors) && (V in factor_do)
            push!(mix_vars, V)
        end
    end

    function rewrite_one(a::CFAtom, in_den::Bool)::CFAtom
        if length(a.left) == 1 && isempty(a.dovars)
            V = a.left[1]
            if V in dual_fixed_vars
                return CFAtom("1", String[], String[], String[], String[])
            elseif V in mix_vars
                return CFAtom(
                    "P(" * star_disp(V) * ")",
                    copy(a.left),
                    String[],
                    [star_disp(V)],
                    String[],
                )
            end
        end

        new_left = copy(a.left_disp)
        for i in eachindex(a.left)
            V = a.left[i]
            if V in event_left_vars
                new_left[i] = base_disp(V)
            end
        end

        new_do = copy(a.do_disp)
        for i in eachindex(a.dovars)
            V = a.dovars[i]
            if V in dual_fixed_vars
                new_do[i] = base_disp(V)
            elseif V in mix_vars
                if !in_den && any(L -> L in event_left_vars, a.left)
                    new_do[i] = base_disp(V)
                else
                    new_do[i] = star_disp(V)
                end
            end
        end

        return rebuild_atom(a; left_disp=new_left, do_disp=new_do)
    end

    num = map_atoms(fr.num, a -> rewrite_one(a, false))
    den = map_atoms(fr.den, a -> rewrite_one(a, true))

    for V in sort!(collect(mix_vars))
        num = wrap_sumvar(num, star_disp(V))
        den = wrap_sumvar(den, star_disp(V))
    end

    return CFFrac(num, den)
end

function remember_fixed_assignment!(
    assignments::Dict{String,String},
    V::String,
    token::String,
    kind::String,
)
    if haskey(assignments, V) && assignments[V] != token
        error("id_cf_step5: conflicting fixed $kind assignments for $V: $(assignments[V]) vs $token")
    end
    assignments[V] = token
end

function canonicalize_fixed_assignment_keys(
    assignments::Dict{String,String},
    known_vars,
)::Dict{String,String}
    isempty(assignments) && return assignments

    lower_to_known = Dict{String,Vector{String}}()
    for V in known_vars
        push!(get!(lower_to_known, lowercase(V), String[]), V)
    end

    canon = Dict{String,String}()
    for (V, token) in assignments
        target = V
        if !(V in known_vars)
            matches = get(lower_to_known, lowercase(V), String[])
            if length(matches) == 1
                target = only(matches)
            end
        end
        remember_fixed_assignment!(canon, target, token, "assignment")
    end
    return canon
end

Base.@kwdef struct Step5DisplayConfig
    symbols::Dict{String,String} = Dict{String,String}()
    value_rename::Dict{String,String} = Dict("x_hat" => "x'")
    mode::Symbol = :legacy  # :legacy | :starred_style
end

Base.@kwdef struct Step5LatentConfig
    mode::Symbol = :exclude_fixed  # :exclude_fixed | :include_fixed
    include::Vector{String} = String[]
end

Base.@kwdef struct Step5RuleConfig
    enabled::Bool = false
    directed_edges::Vector{Tuple{String,String}} = Tuple{String,String}[]
end

Base.@kwdef struct Step5DataConfig
    mode::Symbol = :none  # :none | :interventions | :conditional_queries
    mix_var::String = "D"
    mix_sym::String = "d^*"
    mix_target_var::String = "Z"
    anchor_var::Union{Nothing,String} = "X"
    anchor_token::Union{Nothing,String} = nothing
    mix_token::Union{Nothing,String} = nothing
end

Base.@kwdef struct Step5FixedConfig
    obs::Dict{String,String} = Dict{String,String}()
    do_assignments::Dict{String,String} = Dict{String,String}()
end

function scan_step5_boxes(
    wd::WiringDiagram,
    value_rename::Dict{String,String},
)
    fixed_obs = Dict{String,String}()
    fixed_do  = Dict{String,String}()
    factors = ProbFactor[]

    for b in WD.box_ids(wd)
        nm = box_name(wd, b)

        if (p = parse_obs(nm)) !== nothing
            V, v = p
            fixed_obs[V] = get(value_rename, v, v)
        elseif (p = parse_do(nm)) !== nothing
            V, v = p
            fixed_do[V] = get(value_rename, v, v)
        end

        f = parse_prob_box_name(nm)
        f === nothing && continue
        isempty(f.left) && continue  # P(∅ ; do(...)) = 1
        push!(factors, f)
    end

    return fixed_obs, fixed_do, factors
end

function merge_query_fixed_assignments!(
    fixed_obs::Dict{String,String},
    fixed_do::Dict{String,String},
    query_obs_assignments::Dict{String,String},
    query_do_assignments::Dict{String,String},
    queries,
    value_rename::Dict{String,String},
    manual_fixed_obs::Dict{String,String},
    manual_fixed_do::Dict{String,String},
)
    if queries !== nothing
        for q in queries
            for (V, v) in getproperty(q, :observations)
                token = get(value_rename, string(v), string(v))
                remember_fixed_assignment!(fixed_obs, string(V), token, "obs")
                remember_fixed_assignment!(query_obs_assignments, string(V), token, "obs")
            end
            for (V, v) in getproperty(q, :interventions)
                token = get(value_rename, string(v), string(v))
                remember_fixed_assignment!(fixed_do, string(V), token, "do")
                remember_fixed_assignment!(query_do_assignments, string(V), token, "do")
            end
        end
    end

    for (V, v) in manual_fixed_obs
        remember_fixed_assignment!(fixed_obs, V, v, "obs")
    end
    for (V, v) in manual_fixed_do
        remember_fixed_assignment!(fixed_do, V, v, "do")
    end
end

function canonicalize_all_assignments(
    all_vars,
    fixed_obs::Dict{String,String},
    fixed_do::Dict{String,String},
    query_obs_assignments::Dict{String,String},
    query_do_assignments::Dict{String,String},
)
    fixed_obs = canonicalize_fixed_assignment_keys(fixed_obs, all_vars)
    fixed_do = canonicalize_fixed_assignment_keys(fixed_do, all_vars)
    query_obs_assignments = canonicalize_fixed_assignment_keys(query_obs_assignments, all_vars)
    query_do_assignments = canonicalize_fixed_assignment_keys(query_do_assignments, all_vars)
    return fixed_obs, fixed_do, query_obs_assignments, query_do_assignments
end

function collect_factor_var_sets(factors::Vector{ProbFactor})
    left_var_set = Set{String}()
    do_var_set = Set{String}()
    for f in factors
        union!(left_var_set, f.left)
        union!(do_var_set, f.dovars)
    end
    return left_var_set, do_var_set
end

function prepare_display_sets!(
    factors::Vector{ProbFactor},
    fixed_obs::Dict{String,String},
    fixed_do::Dict{String,String},
    display_mode::Symbol,
)
    left_var_set, do_var_set = collect_factor_var_sets(factors)

    if display_mode == :starred_style
        # If an intervened root variable only appears in do(...) positions,
        # synthesize its marginal root factor P(v^*) for starred display mode.
        for V in sort(collect(keys(fixed_do)))
            if (V in do_var_set) && !(V in left_var_set)
                push!(factors, ProbFactor([V], String[]))
                push!(left_var_set, V)
            end
        end
    end

    star_left_vars = Set{String}()
    if display_mode == :starred_style
        for V in left_var_set
            if (V in do_var_set) || haskey(fixed_obs, V) || haskey(fixed_do, V)
                push!(star_left_vars, V)
            end
        end
    end

    return star_left_vars
end

function choose_latent_vars(
    all_vars,
    output_vars::Vector{String},
    fixed_obs::Dict{String,String},
    fixed_do::Dict{String,String},
    force_latent::Vector{String},
    latent_mode::Symbol,
)
    latent_mode in (:exclude_fixed, :include_fixed) ||
        error("id_cf_step5: unsupported latent.mode=$latent_mode")

    fixed_vars = Set{String}()
    if latent_mode == :exclude_fixed
        union!(fixed_vars, keys(fixed_obs))
        union!(fixed_vars, keys(fixed_do))
    end

    out_set = Set{String}(output_vars)
    latent = String[]
    for v in all_vars
        if !(v in out_set) && !(v in fixed_vars)
            push!(latent, v)
        end
    end
    for v in force_latent
        v in latent || push!(latent, v)
    end
    sort!(latent)
    return latent
end

function maybe_apply_algebraic_rules(
    prod::CFProd,
    sumvars::Vector{String},
    enable_rules::Bool,
    directed_edges::Vector{Tuple{String,String}},
)
    enable_rules || return prod, sumvars

    if !isempty(directed_edges)
        atoms = CFAtom[t for t in prod.terms if t isa CFAtom]
        atoms = rule_drop_do!(atoms, directed_edges)
        prod = CFProd(CFExpr[a for a in atoms])
    end
    return rule_eliminate_sum_normalized!(prod, sumvars)
end

# =========================
# Step 5 main
# =========================
"""
id_cf_step5:

- `display` controls symbol mapping/value labels/rendering mode
- `latent` controls latent-variable selection policy
- `rules` controls optional algebraic rewrite rules
- `data` controls optional data-availability rewrites
- `fixed` injects explicit fixed obs/do assignments

Core workflow:
1. Extract P(...) factors and fixed obs/do labels from Step4 output.
2. Build normalized fraction:
     num = Σ_latent ∏ factors
     den = Σ_{output_vars, latent} ∏ factors
3. Apply optional simplification/rules/rewrite layers.
"""
function id_cf_step5(
    wd::WiringDiagram;
    output_vars::Vector{String},
    simplify::Bool = true,
    queries=nothing,
    display::Step5DisplayConfig = Step5DisplayConfig(),
    latent::Step5LatentConfig = Step5LatentConfig(),
    rules::Step5RuleConfig = Step5RuleConfig(),
    data::Step5DataConfig = Step5DataConfig(),
    fixed::Step5FixedConfig = Step5FixedConfig(),
)
    display.mode in (:legacy, :starred_style) ||
        error("id_cf_step5: unsupported display.mode=$(display.mode)")
    latent.mode in (:exclude_fixed, :include_fixed) ||
        error("id_cf_step5: unsupported latent.mode=$(latent.mode)")
    data.mode in (:none, :interventions, :conditional_queries) ||
        error("id_cf_step5: unsupported data.mode=$(data.mode)")

    var2sym = display.symbols
    value_rename = display.value_rename

    query_obs_assignments = Dict{String,String}()
    query_do_assignments = Dict{String,String}()

    fixed_obs, fixed_do, factors = scan_step5_boxes(wd, value_rename)
    merge_query_fixed_assignments!(
        fixed_obs,
        fixed_do,
        query_obs_assignments,
        query_do_assignments,
        queries,
        value_rename,
        fixed.obs,
        fixed.do_assignments,
    )

    resolved_anchor_token = data.anchor_token
    if resolved_anchor_token === nothing && data.anchor_var !== nothing
        if haskey(fixed_obs, data.anchor_var)
            resolved_anchor_token = fixed_obs[data.anchor_var]
        end
    end

    resolved_mix_token = data.mix_token
    if resolved_mix_token === nothing
        if haskey(fixed_do, data.mix_var)
            resolved_mix_token = fixed_do[data.mix_var]
        else
            resolved_mix_token = get(var2sym, data.mix_var, lowercase(data.mix_var))
        end
    end

    all_vars = collect_all_vars(factors)
    fixed_obs, fixed_do, query_obs_assignments, query_do_assignments =
        canonicalize_all_assignments(all_vars, fixed_obs, fixed_do, query_obs_assignments, query_do_assignments)

    star_left_vars = prepare_display_sets!(factors, fixed_obs, fixed_do, display.mode)

    atoms = CFExpr[]
    for f in factors
        push!(atoms, factor_to_atom(
            f, fixed_obs, fixed_do, var2sym;
            display_mode=display.mode,
            star_left_vars=star_left_vars,
        ))
    end
    prod = CFProd(atoms)

    latent_vars = choose_latent_vars(
        all_vars,
        output_vars,
        fixed_obs,
        fixed_do,
        latent.include,
        latent.mode,
    )

    num_sumvars_vars = copy(latent_vars)
    num_prod = prod

    den_sumvars_vars = vcat(copy(output_vars), latent_vars)
    den_prod = prod

    if simplify
        num_prod, num_sumvars_vars = simplify_sum(num_prod, num_sumvars_vars)
        den_prod, den_sumvars_vars = simplify_sum(den_prod, den_sumvars_vars)
    end

    num_prod, num_sumvars_vars = maybe_apply_algebraic_rules(
        num_prod,
        num_sumvars_vars,
        rules.enabled,
        rules.directed_edges,
    )
    den_prod, den_sumvars_vars = maybe_apply_algebraic_rules(
        den_prod,
        den_sumvars_vars,
        rules.enabled,
        rules.directed_edges,
    )

    sum_symbol(v::String) = get(var2sym, v, default_var_symbol(v))
    num_sumvars = sum_symbol.(num_sumvars_vars)
    den_sumvars = sum_symbol.(den_sumvars_vars)

    num_expr = isempty(num_sumvars) ? num_prod : CFSum(num_sumvars, num_prod)
    den_expr = isempty(den_sumvars) ? den_prod : CFSum(den_sumvars, den_prod)

    event_mode = isempty(output_vars)
    raw_expr = event_mode ? num_expr : CFFrac(num_expr, den_expr)
    raw_tex  = cf_latex(raw_expr)

    simp_expr = raw_expr
    simp_tex  = raw_tex

    zv = data.mix_target_var
    z_sym = get(var2sym, zv, lowercase(zv))

    if event_mode || data.mode == :none
        data_expr = simp_expr
    elseif data.mode == :conditional_queries
        simp_expr isa CFFrac || error("id_cf_step5: :conditional_queries requires a fractional expression")
        data_expr = rewrite_conditional_queries(
            simp_expr::CFFrac,
            factors,
            fixed_obs,
            fixed_do,
            query_obs_assignments,
            query_do_assignments,
            output_vars,
            queries,
            var2sym,
        )
    else
        data_expr = rewrite_for_data(
            simp_expr;
            data_mode=data.mode,
            mix_var=data.mix_var,
            mix_sym=data.mix_sym,
            z_var=zv,
            z_sym=z_sym,
            obs_token=resolved_anchor_token,
            mix_token=resolved_mix_token,
            anchor_is_mix=(
                data.anchor_var !== nothing &&
                data.anchor_var == data.mix_var
            )
        )
    end
    data_tex = cf_latex(data_expr)

    return (
        raw_frac = raw_expr,
        raw_tex = raw_tex,
        simplified_tex = simp_tex,
        data_tex = data_tex
    )
end


function _infer_display_syms(
    queries::Vector{CounterfactualQuery},
    output_vars,
)
    display = Set{Symbol}()
    for q in queries
        for (v, _) in q.observations
            push!(display, v)
        end
        for (v, _) in q.interventions
            push!(display, v)
        end
        for v in q.outputs
            push!(display, v)
        end
    end
    for v in output_vars
        push!(display, Symbol(v))
    end
    return display
end

function _trace_filename(prefix::String, idx::Int, stage::String)
    base = lpad(string(idx), 2, '0') * "_" * stage * ".chyp"
    return isempty(prefix) ? base : prefix * "_" * base
end

function _maybe_write_trace!(
    trace_paths::Dict{String,String},
    trace_dir::Union{Nothing,AbstractString},
    file::String,
    wd::WiringDiagram;
    name::String,
)
    trace_dir === nothing && return nothing
    if !@isdefined(write_chyp)
        error("trace_dir was provided but write_chyp is not defined. include(\"chyp_export.jl\") first.")
    end
    path = joinpath(trace_dir, file)
    write_chyp(path, wd; name=name)
    trace_paths[file] = path
    return path
end

"""
Run the full counterfactual-identification pipeline.

Input:
- `model_or_base`: either a model accepted by `graph_b_to_scm` (for example
  `ConfoundedModel` / `ADMGModel`) or a prebuilt base `WiringDiagram`.
- `queries`: `Vector{CounterfactualQuery}`.

Output (NamedTuple):
- `identifiable`: whether Step 4 succeeded.
- `formula_available`: whether Step 5 produced a formula.
- `raw_tex` / `simplified_tex` / `data_tex`: rendered expressions (or `nothing`).
- `failure_stage`: `nothing`, `:build`, `:simplify`, `:step4`, or `:step5`.
- `error`: error string when failed.
- `step3_blockers`: unabsorbed lambda-box names before Step 4.
- `base_scm`, `wd`: intermediate diagrams.
- `trace_paths`: generated `.chyp` files when `trace_dir` is set.
- `timings_ms`: stage timings in milliseconds.
"""
function identify_counterfactual(
    model_or_base,
    queries::Vector{CounterfactualQuery};
    display_syms::Vector{Symbol}=Symbol[],
    output_vars=String[],
    run_simplify::Bool=true,
    simplify::Bool=true,
    step4_verbose::Bool=false,
    display::Step5DisplayConfig=Step5DisplayConfig(),
    latent::Step5LatentConfig=Step5LatentConfig(),
    rules::Step5RuleConfig=Step5RuleConfig(),
    data::Step5DataConfig=Step5DataConfig(),
    fixed::Step5FixedConfig=Step5FixedConfig(),
    trace_dir::Union{Nothing,AbstractString}=nothing,
    trace_prefix::String="",
)
    total_t0 = time_ns()

    base_scm = nothing
    wd = nothing
    step3_blockers = String[]
    step5_result = nothing
    raw_tex = nothing
    simplified_tex = nothing
    data_tex = nothing

    identifiable = false
    formula_available = false
    failure_stage = nothing
    error_msg = nothing

    build_ms = 0.0
    simplify_ms = 0.0
    step4_ms = 0.0
    step5_ms = 0.0

    trace_paths = Dict{String,String}()

    output_vars_str = String.(output_vars)
    inferred_display = _infer_display_syms(queries, output_vars_str)
    display_var = union(Set(display_syms), inferred_display)

    finalize = () -> begin
        total_ms = (time_ns() - total_t0) / 1e6
        return (
            identifiable = identifiable,
            formula_available = formula_available,
            raw_tex = raw_tex,
            simplified_tex = simplified_tex,
            data_tex = data_tex,
            result = step5_result,
            failure_stage = failure_stage,
            error = error_msg,
            display_vars = sort!(collect(display_var)),
            step3_blockers = step3_blockers,
            base_scm = base_scm,
            wd = wd,
            trace_paths = trace_paths,
            timings_ms = (
                build = build_ms,
                simplify = simplify_ms,
                step4 = step4_ms,
                step5 = step5_ms,
                total = total_ms,
            ),
        )
    end

    build_t0 = time_ns()
    try
        if model_or_base isa WiringDiagram
            base_scm = deepcopy(model_or_base)
        else
            if isempty(display_var)
                base_scm = graph_b_to_scm(model_or_base)
            else
                base_scm = graph_b_to_scm(model_or_base; outputs=display_var)
            end
        end
        _maybe_write_trace!(
            trace_paths,
            trace_dir,
            _trace_filename(trace_prefix, 1, "base_scm"),
            base_scm;
            name="base_scm",
        )

        wd = build_multiverse(base_scm, queries)
        _maybe_write_trace!(
            trace_paths,
            trace_dir,
            _trace_filename(trace_prefix, 2, "multiverse"),
            wd;
            name="multiverse",
        )
        build_ms = (time_ns() - build_t0) / 1e6
    catch err
        build_ms = (time_ns() - build_t0) / 1e6
        failure_stage = :build
        error_msg = sprint(showerror, err)
        return finalize()
    end

    simplify_t0 = time_ns()
    try
        if run_simplify
            drop_discard_branches!(wd)
            _maybe_write_trace!(
                trace_paths,
                trace_dir,
                _trace_filename(trace_prefix, 3, "simplify1"),
                wd;
                name="s1",
            )
            separate_sharp_effect_copy_maps!(wd)
            _maybe_write_trace!(
                trace_paths,
                trace_dir,
                _trace_filename(trace_prefix, 4, "simplify2"),
                wd;
                name="s2",
            )
            merge_identical_deterministic_boxes(wd)
            _maybe_write_trace!(
                trace_paths,
                trace_dir,
                _trace_filename(trace_prefix, 5, "simplify3"),
                wd;
                name="s3",
            )
            replace_fx_with_cx!(wd)
            _maybe_write_trace!(
                trace_paths,
                trace_dir,
                _trace_filename(trace_prefix, 6, "simplify4"),
                wd;
                name="s4",
            )
        end
        simplify_ms = (time_ns() - simplify_t0) / 1e6
    catch err
        simplify_ms = (time_ns() - simplify_t0) / 1e6
        failure_stage = :simplify
        error_msg = sprint(showerror, err)
        return finalize()
    end

    step3_blockers = sort!(unique(String(box(wd, b).value) for b in unabsorbed_lambda_boxes(wd)))

    step4_t0 = time_ns()
    try
        id_cf_step4!(wd, display_var; verbose=step4_verbose)
        identifiable = true
        _maybe_write_trace!(
            trace_paths,
            trace_dir,
            _trace_filename(trace_prefix, 7, "step4"),
            wd;
            name="step4",
        )
        step4_ms = (time_ns() - step4_t0) / 1e6
    catch err
        step4_ms = (time_ns() - step4_t0) / 1e6
        failure_stage = :step4
        error_msg = sprint(showerror, err)
        return finalize()
    end

    step5_t0 = time_ns()
    try
        step5_result = id_cf_step5(
            wd;
            output_vars=output_vars_str,
            simplify=simplify,
            queries=queries,
            display=display,
            latent=latent,
            rules=rules,
            data=data,
            fixed=fixed,
        )
        raw_tex = step5_result.raw_tex
        simplified_tex = step5_result.simplified_tex
        data_tex = step5_result.data_tex
        formula_available = true
        step5_ms = (time_ns() - step5_t0) / 1e6
    catch err
        step5_ms = (time_ns() - step5_t0) / 1e6
        failure_stage = :step5
        error_msg = sprint(showerror, err)
    end

    return finalize()
end

run_cf_pipeline(args...; kwargs...) = identify_counterfactual(args...; kwargs...)
