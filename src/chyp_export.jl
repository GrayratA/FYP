# chyp_export.jl
# Robust exporter: Catlab WiringDiagram -> Chyp via add_junctions
# - Avoids "live_wires" incremental compilation pitfalls
# - Uses explicit junction generators: copy(1->2), merge(2->1), del(1->0)
#
# Requires utils.jl for:
#   - safe_box_name(wd,b)
#   - nin(wd,b), nout(wd,b)
#   - topological_order(wd)
#   - n_in_boundary(wd), n_out_boundary(wd)
#
const WD = Catlab.WiringDiagrams
using Catlab
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams: add_junctions
include("utils.jl")

# ---------- string utils ----------
sanitize_ident(s::AbstractString) = replace(s, r"[^A-Za-z0-9_]" => "_")

tensor_ids(n::Int) =
    n <= 0 ? "" :
    n == 1 ? "id" :
    join(fill("id", n), " * ")

compose_term(a::String, b::String) = a == "" ? b : "$(a) ; $(b)"

function sw_term(perm0::Vector{Int})
    # perm0 is 0-based mapping: new[i] = old[perm0[i]]
    is_id = all(perm0[i] == i-1 for i in 1:length(perm0))
    return is_id ? "" : "sw[" * join(string.(perm0), ", ") * "]"
end

# move selected positions pos (1-based) into a contiguous block starting at k (1-based)
function make_block_perm(N::Int, pos::Vector{Int}, k::Int)
    sel = sort(pos)
    others = [i for i in 1:N if !(i in sel)]
    new_order = Int[]
    append!(new_order, others[1:k-1])
    append!(new_order, sel)
    append!(new_order, others[k:end])
    return [i-1 for i in new_order]  # 0-based
end

# ---------- junction naming ----------
function chyp_box_name(wd::WiringDiagram, bid::Int)
    nm = safe_box_name(wd, bid)
    if nm === nothing
        return "box_$(bid)"
    else
        return sanitize_ident(String(nm))
    end
end

function is_junction_box(wd::WiringDiagram, bid::Int)
    # Heuristic: junctions inserted by add_junctions are typically anonymous
    safe_box_name(wd, bid) === nothing
end

function junction_gen_name(wd::WiringDiagram, bid::Int)
    m = nin(wd, bid)
    n = nout(wd, bid)
    if m == 1 && n == 2
        return "copy"
    elseif m == 2 && n == 1
        return "merge"
    elseif m == 1 && n == 0
        return "del"
    else
        # fallback: keep a unique name
        return chyp_box_name(wd, bid)
    end
end

# ---------- core exporter ----------
"""
Export a Catlab WiringDiagram to a Chyp term by:
1) Converting to a junctionized diagram: wdJ = add_junctions(wd)
2) Compiling boxes in topological order using a "live wire list" strategy,
   but now the diagram is strictly compositional due to explicit junctions.
3) Mapping anonymous junction boxes to copy/merge/del gens.
"""
function wd_to_chyp(wd::WiringDiagram; name::String="main")
    # 1) Make structure explicit
    wdJ = add_junctions(wd)

    in_id  = input_id(wdJ)
    out_id = output_id(wdJ)

    # 2) Topological order (only internal boxes > 0)
    bs = [b for b in topological_order(wdJ) if b > 0]

    # sanity: internal boxes are 1..nboxes
    all_internal = collect(1:WD.nboxes(wdJ))
    if Set(bs) != Set(all_internal)
        missing = collect(setdiff(Set(all_internal), Set(bs)))
        extra   = collect(setdiff(Set(bs), Set(all_internal)))
        error("topo failed on junctionized wd: missing=$(missing) extra=$(extra)")
    end

    # 3) Generator declarations
    gen_specs = Dict{Tuple{String,Int,Int},Bool}()

    # Always include junction structure gens
    gen_specs[("copy", 1, 2)] = true
    gen_specs[("merge", 2, 1)] = true
    gen_specs[("del", 1, 0)] = true

    for bid in bs
        g = is_junction_box(wdJ, bid) ? junction_gen_name(wdJ, bid) : chyp_box_name(wdJ, bid)
        m = nin(wdJ, bid)
        n = nout(wdJ, bid)
        gen_specs[(g, m, n)] = true
    end

    gen_lines = String[]
    for ((g,m,n), _) in sort!(collect(gen_specs); by=x->x[1])
        push!(gen_lines, "gen $(g) : $(m) -> $(n)")
    end

    # 4) Live wires represent the current ordered tensor of available outputs.
    # Boundary input wires count:
    live_wires = [(in_id, i) for i in 1:n_in_boundary(wdJ)]
    term = ""

    # Find unique source feeding (tgt_box, tgt_inport)
    function find_source_of_input(tgt_box::Int, tgt_inport::Int)
        for w in WD.wires(wdJ)
            if w.target.box == tgt_box && w.target.port == tgt_inport
                return (w.source.box, w.source.port)
            end
        end
        error("No source found for input (box=$tgt_box, inport=$tgt_inport)")
    end

    # 5) Compile boxes
    for bid in bs
        g = is_junction_box(wdJ, bid) ? junction_gen_name(wdJ, bid) : chyp_box_name(wdJ, bid)
        m = nin(wdJ, bid)
        n = nout(wdJ, bid)

        # 0-input boxes: tensor on the left
        if m == 0
            gen_term = join(filter(!isempty, [g, tensor_ids(length(live_wires))]), " * ")
            term = compose_term(term, gen_term)
            produced = [(bid, j) for j in 1:n]
            live_wires = vcat(produced, live_wires)
            continue
        end

        # required inputs
        req = [find_source_of_input(bid, j) for j in 1:m]

        # locate them
        pos = Int[]
        for r in req
            p = findfirst(==(r), live_wires)
            if p === nothing
                # Debug print to pinpoint structure
                println("---- EXPORT DEBUG ----")
                println("box bid=$bid g=$g m=$m n=$n")
                println("need wire r=$r  (producer name=", safe_box_name(wdJ, r[1]), ")")
                println("live_wires=", live_wires)
                println("topo bs=", bs)
                error("Wire $r not available when compiling box $bid ($g).")
            end
            push!(pos, p)
        end

        # make contiguous
        k = minimum(pos)
        perm0 = make_block_perm(length(live_wires), pos, k)
        @assert sort(perm0) == collect(0:length(live_wires)-1) "perm0 not a permutation!"
        sw = sw_term(perm0)
        if sw != ""
            term = compose_term(term, sw)
            live_wires = live_wires[perm0 .+ 1]
        end

        # reorder within block if necessary
        block = live_wires[k:k+m-1]
        if block != req
            local_pos = [findfirst(==(r), block) for r in req]
            any(x -> x === nothing, local_pos) && error("Internal alignment failed for box $bid")

            perm = collect(1:length(live_wires))
            for j in 1:m
                perm[k + j - 1] = k + local_pos[j] - 1
            end
            perm2_0 = [p-1 for p in perm]
            @assert sort(perm2_0) == collect(0:length(live_wires)-1) "perm2 not a permutation!"
            sw2 = sw_term(perm2_0)
            if sw2 != ""
                term = compose_term(term, sw2)
                live_wires = live_wires[perm2_0 .+ 1]
            end
        end

        # apply generator at position k
        left  = k-1
        right = length(live_wires) - left - m
        gen_term = join(filter(!isempty, [tensor_ids(left), g, tensor_ids(right)]), " * ")
        term = compose_term(term, gen_term)

        # update wires
        produced = [(bid, j) for j in 1:n]
        live_wires = vcat(live_wires[1:k-1], produced, live_wires[k+m:end])
    end

    # 6) Align boundary outputs (allow extra wires)
    out_m = n_out_boundary(wdJ)
    if out_m > 0
        desired = [(find_source_of_input(out_id, j)) for j in 1:out_m]

        used = falses(length(live_wires))
        perm_front = Int[]
        for d in desired
            p = findfirst(==(d), live_wires)
            p === nothing && error("Desired output wire $d not found in final wire list.")
            push!(perm_front, p)
            used[p] = true
        end
        perm_tail = [i for i in 1:length(live_wires) if !used[i]]
        perm = vcat(perm_front, perm_tail)

        perm0 = [p-1 for p in perm]
        @assert sort(perm0) == collect(0:length(live_wires)-1) "final perm not a permutation!"
        sw = sw_term(perm0)
        if sw != ""
            term = compose_term(term, sw)
            live_wires = live_wires[perm0 .+ 1]
        end
        # ---- drop extra wires by deleting the tail ----
        # After alignment, the first out_m wires correspond to boundary outputs.
        # Any remaining wires are "garbage" and will appear as extra horizontal lines.
        extra = length(live_wires) - out_m
        if extra > 0
            # delete from the end, one wire at a time
            for _ in 1:extra
                # apply del on the last wire: id^(k-1) * del
                k = length(live_wires)
                del_term = join(filter(!isempty, [tensor_ids(k-1), "del"]), " * ")
                term = compose_term(term, del_term)

                # update live_wires: remove the last wire
                pop!(live_wires)
            end
        end
    end

    chyp = join(gen_lines, "\n") * "\n\n" *
           "let $(sanitize_ident(name)) = $(term == "" ? "id0" : term)\n"
    return chyp
end

function write_chyp(path::AbstractString, wd::WiringDiagram; name::String="main")
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
    open(path, "w") do io
        write(io, wd_to_chyp(wd; name=name))
    end
    return path
end