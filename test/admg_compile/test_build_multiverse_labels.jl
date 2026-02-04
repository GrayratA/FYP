using Test

# Fig.9 base model
model = ConfoundedModel(
    [
        :X => :W,
        :W => :Y,
        :D => :Z,
        :Z => :Y
    ],
    Dict(:R => [:X, :Y])
)

base = graph_b_to_scm(model)

queries = [
    CounterfactualQuery(:World1, Dict(:X => :x), Dict{Symbol,Symbol}(), [:Y]),
    CounterfactualQuery(:World2, Dict(), Dict(:X => :x0, :D => :d), Symbol[]),
    CounterfactualQuery(:World3, Dict(:D => :d), Dict(:Z => :z), Symbol[]),
]

wd = build_multiverse(base, queries)


vals = String[]
for bid in box_ids(wd)
    if bid == input_id(wd) || bid == output_id(wd); continue; end
    push!(vals, string(box(wd, bid).value))
end

@test any(s -> occursin("do_X", s) && occursin("=x", s), vals)
@test any(s -> occursin("obs_X", s) && occursin("=x0", s), vals)
@test any(s -> occursin("obs_D", s) && occursin("=d", s), vals)
@test any(s -> occursin("obs_Z", s) && occursin("=z", s), vals)

@test length(output_ports(wd)) == sum(length(q.outputs) for q in queries)
