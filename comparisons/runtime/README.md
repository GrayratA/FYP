Runtime comparison scripts for the Julia prototype and R `cfid`.

These benchmarks measure in-process runtime after startup.
They report:

- `cold_total_ms`: first full call in a fresh process
- `warm_total_min_ms`, `warm_total_median_ms`, `warm_total_mean_ms`, `warm_total_max_ms`: repeated full-call timings in the same process
- Stage medians for warm runs:
  - Julia: `warm_setup_median_ms`, `warm_identify_median_ms`, `warm_build_median_ms`, `warm_simplify_median_ms`, `warm_step4_median_ms`, `warm_step5_median_ms`
  - R `cfid`: `warm_setup_median_ms`, `warm_identifiable_median_ms` and alias `warm_identify_median_ms`

Example usage:

```powershell
julia --project=. comparisons/runtime/julia_benchmark.jl 30 5
Rscript comparisons/runtime/r_cfid_benchmark.R 30 5
```

The current scripts benchmark three examples:

- `drug`
- `party`
- `hepar2_conditional`

Notes on alignment:

- Both implementations include all preprocessing done inside each script's run function in `total`.
- For `party`, Julia uses `data_mode=:none` and R `cfid` omits the `data` argument to avoid a data-mode mismatch.

## Structural Tests

Structural-scaling and profiling artifacts are under:

- `comparisons/runtime/structural_tests/scripts`: runnable scripts
- `comparisons/runtime/structural_tests/raw`: raw text outputs
- `comparisons/runtime/structural_tests/processed`: generated CSV tables
- `comparisons/runtime/structural_tests/profiles`: Rprof output files

Example usage:

```powershell
julia comparisons/runtime/structural_tests/scripts/struct_scaling_julia.jl 4 2 8,16,24,32 `
  | Tee-Object -FilePath comparisons/runtime/structural_tests/raw/struct_julia_r4w2.txt

Rscript comparisons/runtime/structural_tests/scripts/struct_scaling_r.R 4 2 8,16,24,32 `
  | Tee-Object -FilePath comparisons/runtime/structural_tests/raw/struct_r_r4w2.txt

pwsh -File comparisons/runtime/structural_tests/scripts/make_struct_csv.ps1 -Tag r4w2
```
