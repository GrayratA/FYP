Runtime comparison scripts for the Julia prototype and R `cfid`.

These benchmarks measure in-process runtime after startup.
They report:

- `cold_ms`: the first call in a fresh process
- `min_ms`, `median_ms`, `mean_ms`, `max_ms`: repeated-call timings in the same process

Example usage:

```powershell
julia --project=. comparisons/runtime/julia_benchmark.jl 30
Rscript comparisons/runtime/r_cfid_benchmark.R 30
```

The current scripts benchmark two examples:

- `drug`
- `party`
