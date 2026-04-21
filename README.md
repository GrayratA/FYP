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
