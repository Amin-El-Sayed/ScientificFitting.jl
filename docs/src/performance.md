# Performance

JuFitter optimizes the mathematically important paths first: avoid explicit
matrix inverses, reuse static covariance factorizations, keep simple
least-squares fits on the fast backend, and make expensive choices explicit.

## Startup Probe

To check the first requirement for responsive command-line use, run:

```bash
julia --project=. --startup-file=no benchmarks/startup_probe.jl \
  --save=/tmp/jufitter-startup-probe.toml
```

The probe starts a fresh Julia process, loads `JuFitter`, verifies that neither
`Makie` nor `CairoMakie` was loaded, records elapsed wall time, and can write a
small TOML summary. This is a startup smoke test, not a stable latency claim.
Use it to catch regressions such as accidentally moving plotting dependencies
back into the fitting/reporting core.

## Running BenchmarkTools Benchmarks

The benchmark entry point is:

```bash
julia --project=benchmarks benchmarks/runbenchmarks.jl --seconds=1
```

The benchmark script develops the repository package into the benchmark
environment on first run. The generated `benchmarks/Manifest.toml` is local
machine state and is intentionally ignored by git.

For quick local checks, load the suite and run selected cases:

```julia
using BenchmarkTools
include("benchmarks/runbenchmarks.jl")

run(SUITE["fit"]["linear_100"]; seconds=0.2)
run(SUITE["fit"]["full_covariance_500_bounded"]; seconds=0.2)
run(SUITE["likelihood"]["poisson_5000"]; seconds=0.2)
```

To write a machine-readable baseline:

```bash
julia --project=benchmarks benchmarks/runbenchmarks.jl \
  --seconds=1 \
  --save=benchmarks/output/local-baseline.toml
```

The saved TOML includes Julia, JuFitter, OS, CPU, Julia-thread, BLAS-thread,
git-commit, benchmark-count, timing, memory, and allocation metadata. Keep the
file with the reviewed release evidence, not as a generic repository artifact.

To compare against a saved baseline:

```bash
julia --project=benchmarks benchmarks/runbenchmarks.jl \
  --seconds=1 \
  --compare=benchmarks/output/local-baseline.toml \
  --tolerance=0.25
```

The comparison fails if a benchmark case is missing from either side. This is
intentional: adding, removing, or renaming a benchmark changes the evidence set
and requires a new reviewed baseline.

The comparison also checks release-relevant metadata: Julia version, operating
system, CPU name, machine target, Julia thread count, BLAS thread count, and
reported units must match. This prevents accidentally treating a run on another
machine or thread configuration as release evidence. For exploratory local
comparisons only, add `--allow-metadata-mismatch`; do not use that flag for
public performance claims.

Plot export benchmarks are opt-in because they require CairoMakie:

```bash
julia --project=benchmarks benchmarks/runbenchmarks.jl --plot --seconds=1
```

The current suite covers:

- linear and nonlinear least-squares fits,
- out-of-place and in-place 10k-point linear fits under identical data,
- no-op bounds that must preserve the fast path,
- dense covariance fits,
- bounded dense-covariance fits through `Optimization.jl`,
- Poisson likelihood fits,
- plot export,
- profile scans,
- profile-matrix diagnostics with profile and pairwise contour refits.

Benchmarks are not release claims until they are compared against saved
baselines on stated hardware or CI runners. Local baselines belong under
`benchmarks/output/` until a specific machine or CI environment is chosen as
the release reference.

## Performance Budget Gate

The repository also has a small steady-state gate:

```bash
julia --project=. --startup-file=no test/performance_budget_gate.jl
```

It is not a benchmark report. The gate warms compilation first and then checks
only broad budgets for representative hot paths: a 10k-point linear
least-squares fit, the same fit with no-op bounds, and a 300-point dense
covariance fit. Its job is to catch accidental slowdowns such as losing the
`LsqFit` fast path or recomputing static covariance work inside the objective.

For slow shared runners, the budgets can be scaled with
`JUFITTER_PERFORMANCE_BUDGET_SCALE`. Any release claim still needs a proper
`BenchmarkTools` run from `benchmarks/runbenchmarks.jl` and a recorded baseline.

## Fast Paths

Unbounded static chi-square fits use `LsqFit`. This includes fits without
active bounds, constraints, priors, or parameter-dependent covariance. Bounds of
the form `[-Inf, Inf]` are treated as no-op bounds and keep the fast path open.

When `LsqFit` computes a weighted Jacobian, JuFitter reuses it during
`FitResult` construction. This avoids an unnecessary second automatic
differentiation pass for common least-squares workflows.

For models whose output buffer is large or whose implementation would otherwise
allocate temporaries, use the in-place contract:

```julia
function model!(out, x, p)
    @. out = p[1] * exp(-p[2] * x) + p[3]
    return nothing
end

result = fit_model(
    model!,
    x,
    y;
    p0=[1.0, 0.5, 0.0],
    sigma_y=sigma_y,
    inplace=true,
)
```

An analytic in-place Jacobian uses `jacobian!(J, x, p)` and is passed with
`jacobian=jacobian!`. The unbounded least-squares backend forwards these
functions to LsqFit's native in-place interface. Bounds, constraints, and
parameter-dependent covariance still use the general optimizer, but preserve
the same mutating model contract with an output buffer of the correct AD type.
JuFitter evaluates both mutating functions once at `p0` and rejects incomplete
or non-finite output buffers before solver dispatch. Keep their buffer and
parameter signatures generic rather than restricting them to `Float64`; AD-based
optimizer paths evaluate them with dual-number element types.

## Covariance Costs

Diagonal covariance stores precomputed inverse standard deviations and log
determinants. Dense static covariance stores a Cholesky factor and log
determinant. Sparse static `cov_y` matrices are kept sparse on the unbounded
least-squares path and whitened with CHOLMOD's permuted factor. Residuals are
whitened using linear solves, not explicit covariance inverses.

Parameter-dependent covariance is intentionally more expensive. X
uncertainties and model-relative y uncertainties can change the effective
covariance at every parameter point, so the cost function must recompute the
relevant covariance terms.

For x uncertainties, the default path differentiates the model with respect to
each x value by automatic differentiation. That keeps the simple API correct,
but it is not the preferred route for large datasets. If the derivative is
known, pass a vectorized function:

```julia
x_derivative(x, p) = @. p[1] * exp(p[2] * x) * p[2]

result = fit_model(
    model,
    x,
    y;
    p0=[1.0, -0.2],
    sigma_y=sigma_y,
    sigma_x=sigma_x,
    x_derivative=x_derivative,
)
```

JuFitter validates that `x_derivative(x, p)` has the same length as `x` and
contains finite values. The derivative may depend on `p`; AD information is
preserved when the Gaussian NLL needs gradients or Hessians.

## Known Limits

Dense covariance matrices are correct and tested, but they are not the right
representation for huge correlated datasets. They require `O(n^2)` memory and
`O(n^3)` factorization time. Sparse static `cov_y` is the current first
structured path for least-squares problems; constrained sparse fits, sparse
Gaussian-NLL gradients, structured `cov_x`, and matrix-free whitening operators
remain future work.

If a large dataset truly has correlated uncertainties, the next necessary step
is a structured covariance API: banded matrices, Toeplitz kernels,
low-rank-plus-diagonal structure, sparse precision matrices, or custom
whitening operators. Dense-matrix micro-optimization cannot fix the asymptotic
scaling.

This is especially relevant for long time series, images, spectra, and detector
arrays. In those cases the covariance structure is usually the scientific model;
materializing a dense matrix is often both slower and less expressive than a
dedicated whitening operator.

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
