# Performance

JuFitter optimizes the mathematically important paths first: avoid explicit
matrix inverses, reuse static covariance factorizations, keep simple
least-squares fits on the fast backend, and make expensive choices explicit.

## Running Benchmarks

The benchmark entry point is:

```bash
julia --project=. benchmarks/runbenchmarks.jl
```

For quick local checks, load the suite and run selected cases:

```julia
using BenchmarkTools
include("benchmarks/runbenchmarks.jl")

run(SUITE["fit"]["linear_100"]; seconds=0.2)
run(SUITE["fit"]["full_covariance_500_bounded"]; seconds=0.2)
run(SUITE["likelihood"]["poisson_5000"]; seconds=0.2)
```

The current suite covers:

- linear and nonlinear least-squares fits,
- no-op bounds that must preserve the fast path,
- dense covariance fits,
- bounded dense-covariance fits through `Optimization.jl`,
- Poisson likelihood fits,
- plot export,
- profile scans.

Benchmarks are not release claims until they are compared against saved
baselines in CI.

## Fast Paths

Unbounded static chi-square fits use `LsqFit`. This includes fits without
active bounds, constraints, priors, or parameter-dependent covariance. Bounds of
the form `[-Inf, Inf]` are treated as no-op bounds and keep the fast path open.

When `LsqFit` computes a weighted Jacobian, JuFitter reuses it during
`FitResult` construction. This avoids an unnecessary second automatic
differentiation pass for common least-squares workflows.

## Covariance Costs

Diagonal covariance stores precomputed inverse standard deviations and log
determinants. Dense static covariance stores a Cholesky factor and log
determinant. Residuals are whitened using linear solves.

Parameter-dependent covariance is intentionally more expensive. X
uncertainties and model-relative y uncertainties can change the effective
covariance at every parameter point, so the cost function must recompute the
relevant covariance terms.

## Known Limits

Dense covariance matrices are correct and tested, but they are not the right
representation for huge correlated datasets. They require `O(n^2)` memory and
`O(n^3)` factorization time.

If a large dataset truly has correlated uncertainties, the next necessary step
is a structured covariance API: banded matrices, Toeplitz kernels,
low-rank-plus-diagonal structure, sparse precision matrices, or custom
whitening operators. Dense-matrix micro-optimization cannot fix the asymptotic
scaling.

Makie/CairoMakie dominates first-use plotting latency. Plotting is therefore
loaded through an optional package extension: users who only fit and print
reports do not load Makie, while users who call plotting functions opt into
`using CairoMakie`. This separates fitting-engine latency from rendering
latency, but it does not remove CairoMakie's first-use compilation cost for
plot-heavy workflows.

## Rules For Future Optimizations

- Do not make the default faster by silently changing statistical semantics.
- Do not repair invalid covariance matrices without an explicit user-visible
  policy and diagnostic.
- Do not introduce explicit matrix inverses in production calculations.
- Add a torture test for every robustness fix.
- Add a benchmark before changing a hot path.
- Keep fast and full test gates separate: fast gates should exclude expensive
  plot generation; full gates should include plots and documentation assets.
