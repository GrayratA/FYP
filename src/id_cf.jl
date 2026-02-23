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


function id_cf_step4!(wd, display_vars; verbose::Bool=true)
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
# Generic (no cfid-style naming). Optional data-mode rewriting.
# ============================================================

const WD = Catlab.WiringDiagrams

# ---------- factor representation ----------
struct ProbFactor
    left::Vector{String}    # vars on left of P(...)
    dovars::Vector{String}  # vars inside do(...)
end

# ---------- small expr AST ----------
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

# ---------- LaTeX printing ----------
function cf_latex(e::CFExpr)
    if e isa CFAtom
        return (e::CFAtom).s
    elseif e isa CFProd
        return join(cf_latex.((e::CFProd).terms), "")
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

# ---------- parse naming conventions ----------
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
        for tok in split(replace(right, " " => ""), ",")
            mm = match(r"do\((.*)\)", tok)
            mm !== nothing && push!(do_vars, mm.captures[1])
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
    var2sym::Dict{String,String}
)
    left_disp = String[]
    for V in f.left
        push!(left_disp, get(fixed_obs, V, get(var2sym, V, default_var_symbol(V))))
    end
    do_disp = String[]
    for V in f.dovars
        push!(do_disp, get(fixed_do, V, get(var2sym, V, default_var_symbol(V))))
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

# ---------- generic simplification on sums ----------
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

# ---------- drop-do rule helpers ----------
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

# ---------- optional: data availability rewrite (generic Step 4.5) ----------
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
    obs_token::Union{Nothing,String} = nothing  # NEW: anchor token, e.g. "x'"
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

    # recursive search (works even if denominator is a sum/product nest)
    Px = find_atom(fr.den, a -> is_P_of_token(a, tok))

    if Px === nothing
        # do not fail: this pass is a canonicalization for "data mode"
        println("[rewrite_for_data] WARNING: could not find P($tok) in denominator; synthesizing it.")
        println("[rewrite_for_data] den = ", cf_latex(fr.den))
        Px = CFAtom("P($tok)", ["ANCHOR"], String[], [tok], String[])
    end

    # atoms: P(mix*) and P(z|do(mix*))
    Pmix = CFAtom("P($mix_sym)", [mix_star], String[], [mix_sym], String[])
    Pz_do_mix = CFAtom("P($z_sym|do($mix_sym))",
                       [z_var], [mix_star],
                       [z_sym], [mix_sym])

    # numerator: merge sums if numerator already has Σ_{...}
    old_num = fr.num
    if old_num isa CFSum
        s = old_num::CFSum
        innerP = _to_prod(s.body)
        num_inner = CFProd(vcat(innerP.terms, CFExpr[Pmix]))
        num2 = CFSum(vcat(s.vars, [mix_sym]), num_inner)
    else
        numP = _to_prod(old_num)
        num_inner = CFProd(vcat(numP.terms, CFExpr[Pmix]))
        num2 = CFSum([mix_sym], num_inner)
    end

    # denominator: Σ_{mix_sym} P(tok) P(z|do(mix*)) P(mix*)
    den_inner = CFProd(CFExpr[Px, Pz_do_mix, Pmix])
    den2 = CFSum([mix_sym], den_inner)

    return CFFrac(num2, den2)
end

# =========================
# Step 5 main
# =========================
"""
id_cf_step5:

- Extract P(...) probability boxes from Step4 output diagram
- Read obs_/do_ boxes as fixed values (PRINTING only)
- Choose latent vars = all_vars_in_factors - output_vars - fixed_vars
- Build RAW normalized fraction:
    num = Σ_latent ∏ factors
    den = Σ_{output_vars, latent} ∏ factors
- simplify_sum does algebraic marginalization + drop unused sums
- enable_cfid_rules applies: drop-do (with directed_edges), then eliminate Σ_v P(v|...) = 1
- data_mode (Step 4.5): optional canonicalization based on what distributions are "available"
"""
function id_cf_step5(
    wd::WiringDiagram;
    output_vars::Vector{String},
    force_latent::Vector{String} = String[],
    var2sym::Dict{String,String} = Dict{String,String}(),
    obs_value_rename::Dict{String,String} = Dict("x_hat"=>"x'"),
    do_value_rename::Dict{String,String}  = Dict{String,String}(),  # NEW
    simplify::Bool = true,
    enable_cfid_rules::Bool = false,
    directed_edges::Vector{Tuple{String,String}} = Tuple{String,String}[],
    # generic data rewrite
    data_mode::Symbol = :none,       # :none | :interventions
    mix_var::String = "D",           # dummy mixture index base var-name
    mix_sym::String = "d^*",         # printed dummy index
    mix_target_var::String = "Z",    # NEW: which variable is z in P(z|do(d^*))
    anchor_var::Union{Nothing,String} = "X",        # NEW: which var is observed anchor (e.g. X)
    anchor_token::Union{Nothing,String} = nothing   # NEW: printed token for anchor (e.g. "x'"); if nothing, derive from fixed_obs[anchor_var]
)
    # fixed obs/do (for printing)
    fixed_obs = Dict{String,String}()
    fixed_do  = Dict{String,String}()

    for b in WD.box_ids(wd)
        nm = box_name(wd, b)
        if (p = parse_obs(nm)) !== nothing
            V, v = p
            fixed_obs[V] = get(obs_value_rename, v, v)
        elseif (p = parse_do(nm)) !== nothing
            V, v = p
            fixed_do[V] = get(do_value_rename, v, v)  # APPLY rename
        end
    end

    # resolve anchor token for rewrite_for_data
    resolved_anchor_token = anchor_token
    if resolved_anchor_token === nothing && anchor_var !== nothing
        if haskey(fixed_obs, anchor_var)
            resolved_anchor_token = fixed_obs[anchor_var]
        end
    end

    # collect factors
    factors = ProbFactor[]
    for b in WD.box_ids(wd)
        nm = box_name(wd, b)
        f = parse_prob_box_name(nm)
        f === nothing && continue
        push!(factors, f)
    end

    # build product
    atoms = CFExpr[]
    for f in factors
        push!(atoms, factor_to_atom(f, fixed_obs, fixed_do, var2sym))
    end
    prod = CFProd(atoms)

    # latent var-names
    all_vars = collect_all_vars(factors)
    fixed_vars = Set{String}(keys(fixed_obs))
    union!(fixed_vars, keys(fixed_do))
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

    # numerator = Σ_latent ∏ factors
    num_sumvars_vars = copy(latent)
    num_prod = prod

    # denominator = Σ_{output_vars, latent} ∏ factors
    den_sumvars_vars = vcat(copy(output_vars), latent)
    den_prod = prod

    if simplify
        num_prod, num_sumvars_vars = simplify_sum(num_prod, num_sumvars_vars)
        den_prod, den_sumvars_vars = simplify_sum(den_prod, den_sumvars_vars)
    end

    if enable_cfid_rules
        if !isempty(directed_edges)
            num_atoms = CFAtom[t for t in num_prod.terms if t isa CFAtom]
            den_atoms = CFAtom[t for t in den_prod.terms if t isa CFAtom]
            num_atoms = rule_drop_do!(num_atoms, directed_edges)
            den_atoms = rule_drop_do!(den_atoms, directed_edges)
            num_prod = CFProd(CFExpr[a for a in num_atoms])
            den_prod = CFProd(CFExpr[a for a in den_atoms])
        end
        num_prod, num_sumvars_vars = rule_eliminate_sum_normalized!(num_prod, num_sumvars_vars)
        den_prod, den_sumvars_vars = rule_eliminate_sum_normalized!(den_prod, den_sumvars_vars)
    end

    sum_symbol(v::String) = get(var2sym, v, default_var_symbol(v))
    num_sumvars = sum_symbol.(num_sumvars_vars)
    den_sumvars = sum_symbol.(den_sumvars_vars)

    num_expr = isempty(num_sumvars) ? num_prod : CFSum(num_sumvars, num_prod)
    den_expr = isempty(den_sumvars) ? den_prod : CFSum(den_sumvars, den_prod)

    raw_frac = CFFrac(num_expr, den_expr)
    raw_tex  = cf_latex(raw_frac)

    # simplified currently same as raw
    simp_frac = raw_frac
    simp_tex  = raw_tex

    # Optional data rewrite (generic)
    zv = mix_target_var
    z_sym = get(var2sym, zv, lowercase(zv))

    data_frac = rewrite_for_data(
        simp_frac;
        data_mode=data_mode,
        mix_var=mix_var,
        mix_sym=mix_sym,
        z_var=zv,
        z_sym=z_sym,
        obs_token=resolved_anchor_token
    )
    data_tex = cf_latex(data_frac)

    return (
        raw_frac = raw_frac,
        raw_tex = raw_tex,
        simplified_tex = simp_tex,
        data_tex = data_tex
    )
end