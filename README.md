# FYP Counterfactual

Counterfactual identification prototype based on Julia/Catlab, with runtime comparisons against R `cfid`.

## Repository Layout

- `src/`: core Julia implementation (`admg_compile`, `bn_import`, `simplify_cf`, `id_cf`, utilities)
  - `src/123.ipynb` and `src/test.ipynb` are development-time testing notebooks and can be ignored.
- `test/`: Julia test suites (`admg_compile`, `simplify_cf`, `id_cf`)
- `net/`: Bayesian network structure files (for example `hepar2.net`)
- `tools/`: helper scripts (for example BN-to-ADMG conversion)
- `comparisons/`: benchmark scripts and comparison artifacts
  - `comparisons/runtime/`: primary Julia vs R runtime benchmarks
  - `comparisons/runtime/structural_tests/`: structural scaling scripts/data/profiles
- `trace/`: intermediate trace artifacts

## Setup

Prerequisites:

- Julia (project uses `Project.toml` / `Manifest.toml`)
- R with package `cfid` (for R-side benchmarks)

Install Julia dependencies:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Run Tests

Run full Julia tests:

```powershell
julia --project=. test/runtest.jl
```

## Counterfactual API (`identify_counterfactual`)

Use `identify_counterfactual` as the single entry point for:
- building the SCM/multiverse from a graph + queries,
- running simplify + Step4/Step5,
- returning identifiability status and formula text.

Load core files:

```julia
using Catlab
include("src/admg_compile.jl")
include("src/simplify_cf.jl")
include("src/id_cf.jl")
include("src/bn_import.jl")     # optional: only if you read .bif/.net
include("src/chyp_export.jl")   # optional: only if you set trace_dir
```

HEPAR2 three-world conditional example:

```julia
hepar2 = read_bn_structure("net/hepar2.net")
model = hepar2.model

queries = [
    CounterfactualQuery(:Real, Dict{Symbol, Symbol}(), Dict(:age => :age_obs, :sex => :sex_obs), Symbol[]),
    CounterfactualQuery(:CF_pbc, Dict(:age => :age_do), Dict(:pbc => :pbc_obs), Symbol[]),
    CounterfactualQuery(:CF_carc, Dict(:age => :age_do), Dict(:carcinoma => :carc_obs), Symbol[]),
]

res = identify_counterfactual(
    model,
    queries;
    display_syms=[:age, :sex, :pbc, :carcinoma],
    output_vars=["carcinoma"],
    data=Step5DataConfig(mode=:conditional_queries),
    display=Step5DisplayConfig(
        symbols=Dict("age" => "age", "sex" => "sex", "pbc" => "pbc", "carcinoma" => "carcinoma"),
        value_rename=Dict("age_do" => "age"),
    ),
    trace_dir="trace/HEPAR2", # optional
)

println("identifiable = ", res.identifiable)
println("formula      = ", res.data_tex)
println("failure_stage= ", res.failure_stage)
println("error        = ", res.error)
```

Key returned fields:
- `res.identifiable`: `true` if Step4 succeeds.
- `res.formula_available`: `true` if Step5 produced formula text.
- `res.raw_tex`, `res.simplified_tex`, `res.data_tex`: formula strings.
- `res.failure_stage`: `nothing` or one of `:build`, `:simplify`, `:step4`, `:step5`.
- `res.error`: error message string when failed.
- `res.step3_blockers`: unabsorbed lambda boxes before Step4.
- `res.trace_paths`: written `.chyp` paths when `trace_dir` is set.
- `res.timings_ms`: timing breakdown by stage.

Alias:
- `run_cf_pipeline(...)` is an alias of `identify_counterfactual(...)`.

Troubleshooting:
- If `res.identifiable == false`, check `res.failure_stage` and `res.error`.
- If you pass `trace_dir`, ensure `include("src/chyp_export.jl")` was loaded.

## Runtime Benchmarks

Run main runtime comparison:

```powershell
julia --project=. comparisons/runtime/julia_benchmark.jl 30 5
Rscript comparisons/runtime/r_cfid_benchmark.R 30 5
```

For benchmark fields, alignment notes, and structural-test workflow, see:

- `comparisons/runtime/README.md`

## Structural Scaling Benchmarks

Example:

```powershell
julia comparisons/runtime/structural_tests/scripts/struct_scaling_julia.jl 4 2 8,16,24,32 `
  | Tee-Object -FilePath comparisons/runtime/structural_tests/raw/struct_julia_r4w2.txt

Rscript comparisons/runtime/structural_tests/scripts/struct_scaling_r.R 4 2 8,16,24,32 `
  | Tee-Object -FilePath comparisons/runtime/structural_tests/raw/struct_r_r4w2.txt

pwsh -File comparisons/runtime/structural_tests/scripts/make_struct_csv.ps1 -Tag r4w2
```

## Utility Script

Convert `.bif`/`.net` Bayesian network structure to ADMG Julia literal:

```powershell
julia --project=. tools/import_bn_to_admg.jl <input.bif|input.net> [output.jl] [model_name]
```
