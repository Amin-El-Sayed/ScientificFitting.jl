# Performance

ScientificFitting optimizes the mathematically important paths first: avoid explicit
matrix inverses, reuse static covariance factorizations, keep simple
least-squares fits on the fast backend, and make expensive choices explicit.

## Startup Probe

To check the first requirement for responsive command-line use, run:

```bash
julia --project=. --startup-file=no benchmarks/startup_probe.jl \
  --save=/tmp/scientificfitting-startup-probe.toml
```

The probe starts a fresh Julia process, loads `ScientificFitting`, verifies that neither
`Makie` nor `CairoMakie` was loaded, records elapsed wall time, and can write a
small TOML summary. This is a startup smoke test, not a stable latency claim.
Use it to catch regressions such as accidentally moving plotting dependencies
back into the fitting/reporting core.

## Python Startup

After installing the Python package and resolving its Julia environment, run
this to measure each fit family in a **separate fresh Python process**:

```bash
python python/benchmarks/startup.py --case all \
  --output=/tmp/scientificfitting-python-startup.json
```

It separates Python imports, JuliaCall initialization, bridge loading, the first
public fit call (including result conversion), and report/diagnostic generation.
It then fits a newly created Python callable in the same process. Parameters,
covariances, and costs must match analytic references for both calls. The cases
cover every high-level fit family plus the in-place Gaussian interface, using
small linear, constant-rate, and exponential models. They isolate restart
overhead, not realistic large-data throughput. `--case histogram_density`, for
example, runs only that family; the default is `gaussian`.

Julia and Python package versions, package source, and platform accompany the
timings. The first fit includes outstanding compilation, not just optimization.
With an empty environment, JuliaCall initialization also includes provisioning;
do not compare that run with a previously installed environment.

One run per family on macOS ARM64, with Python 3.12.4, NumPy 2.5.3,
JuliaCall 0.9.35, Julia 1.12.7, and the installed 0.2.0.dev0 wrapper against
the development core:

| Fit case | Imports + bridge + first fit / s | Next fit / ms |
|---|---:|---:|
| Gaussian | 16.8 | 1.56 |
| Gaussian, in-place | 21.2 | 1.12 |
| Custom objective | 16.1 | 1.10 |
| User-defined observation likelihood | 15.9 | 1.60 |
| Poisson counts | 14.2 | 1.69 |
| Integrated bin-count model | 18.3 | 3.26 |
| Histogram density with quadrature | 22.5 | 2.60 |
| Unbinned density | 16.3 | 1.77 |
| Extended-unbinned intensity | 19.3 | 2.24 |
| Indexed observations | 15.9 | 1.67 |
| Multiple datasets | 16.3 | 2.53 |

The time in seconds sums the measured import, bridge, and first-fit stages; it
excludes data/reference preparation and reports. Initial report/diagnostic
formatting added 1.5-3.3 s. The next fit includes Python input conversion and
result snapshots, but no imports or report formatting. Density callbacks use
`vectorized=True`. These are individual local observations after installation
and precompilation, not cross-platform guarantees or cross-library rankings.
Millisecond repeated fits do not imply millisecond startup.

The core uses [PrecompileTools](https://julialang.github.io/PrecompileTools.jl/stable/)
with tiny deterministic Gaussian/Poisson fits and reports. The Python entry
uses [`invokelatest`](https://docs.julialang.org/en/v1/base/base/#Base.invokelatest)
once per fit to separate dynamic language conversion from numerical specialization.
Foreign callbacks have a fixed return-type boundary; native Julia callbacks
retain their normal specialization and automatic differentiation. These choices
cache reusable code without executing Python during package precompilation or
building platform-specific sysimages. They increase the installation cache and
do not eliminate Julia's load time; see [Python installation](python.md#Installation-And-Packaging).

## Python Model Callbacks

Run `python python/benchmarks/callbacks.py --output=/tmp/callbacks.json` with
NumPy, SciPy, and the installed Python preview. The probe compares Gaussian,
unbinned exponential, and histogram fits with the same data, starts, bounds,
solver tolerance, and local covariance calculation. Gaussian/event references
are closed form; histogram bin probabilities and curvature are analytic, with
an independent SciPy minimization. Every timed result must pass those checks.

Local medians from three warmed repetitions on macOS ARM64, Python 3.12.4,
NumPy 2.5.3, JuliaCall 0.9.35, and Julia 1.12.7, with 2,000 observations and
40 histogram bins:

| Case | Python scalar density | Python batched model | Julia finite derivatives |
|---|---:|---:|---:|
| Gaussian regression | not applicable | 1.11 ms | 0.45 ms |
| Unbinned exponential | 1,020 ms | 5.75 ms | 2.59 ms |
| Histogram density with quadrature | 216 ms | 46.8 ms | 0.78 ms |

Python times include input conversion and result snapshots; native Julia times
exclude Python conversion. All exclude first-use compilation, plotting, and
report formatting. The JSON also records Julia's default derivative mode;
these measurements are not cross-platform latency guarantees.

`vectorized=True` reduces the unbinned case from 332,000 to 166 Python density
calls per fit. Histogram quadrature still evaluates each bin adaptively and
makes 4,640 batch calls here, so it remains noticeably slower than native Julia.
Where available, supply analytic bin expectations to `fit_histogram_model`:
the same exponential-bin likelihood then takes about 1.64 ms in Python without
quadrature. See [batched density callbacks](python.md#Batched-Event-And-Histogram-Densities)
for the array contract. Batching changes evaluation granularity, not the
likelihood or the requested integration tolerance.

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

The saved TOML includes Julia, ScientificFitting, OS, CPU, Julia-thread, BLAS-thread,
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
- a 100k-point matrix-free structured-whitening fit,
- Poisson likelihood fits,
- plot export,
- profile scans,
- profile-matrix diagnostics with profile and pairwise contour refits.

Benchmarks are not release claims until they are compared against saved
baselines on stated hardware or CI runners. Local baselines belong under
`benchmarks/output/` until a specific machine or CI environment is chosen as
the release reference.

## Ecosystem Throughput Probe

From the repository checkout with Julia 1.11+, prepare the optional dependencies
once, then run the probe:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs --startup-file=no benchmarks/ecosystem.jl 10000 100000
```

This opt-in probe fits a controlled truncated-normal peak plus uniform background
with four free parameters. It compares MIGRAD and bounded L-BFGS, a distribution
factory and BuildConstructors, and the upstream DistributionsHEP objective.
All paths include ScientificFitting's local covariance and result construction.

The script reports warmed median time and allocated bytes over three repetitions.
It requires convergence, parameter differences below 0.01 standard errors,
cost agreement within 0.0002 and relative covariance-norm agreement within 0.001.
Data generation and compilation are excluded;
allocated bytes are cumulative, not peak memory. Yields are scaled by the sample
size for both solvers. Their stopping criteria differ, so agreement of the
results, not equal numeric tolerances, is the comparison requirement.

Use the reported versions and thread settings when repeating a measurement.
This probes integration overhead and conditioning; it is neither a CERN-data
analysis nor a universal solver ranking.

Heterogeneous mixtures evaluate each upstream component on the event array, then
combine the log densities with a scaled weighted sum. Homogeneous mixtures keep
their specialized native loop. Zero weights retain their parameter derivatives;
buffers are local to an objective call rather than shared between fits or stored
in an event-by-component matrix. Other distributions keep their upstream
`loglikelihood` implementation.

## Performance Budget Gate

The repository also has a small steady-state gate:

```bash
julia --project=. --startup-file=no test/performance_budget_gate.jl
```

It is not a benchmark report. The gate warms compilation first and then checks
only broad budgets for representative hot paths: a 10k-point linear
least-squares fit, the same fit with no-op bounds, and a 300-point dense
covariance fit. It also checks a 50k-point structured-whitening fit and verifies
that allocations remain linear between 10k and 50k points. Its job is to catch
accidental slowdowns such as losing the `LsqFit` fast path or recomputing static
covariance work inside the objective.

For slow shared runners, the budgets can be scaled with
`SCIENTIFICFITTING_PERFORMANCE_BUDGET_SCALE`. Any release claim still needs a proper
`BenchmarkTools` run from `benchmarks/runbenchmarks.jl` and a recorded baseline.

## Fast Paths

Unbounded static chi-square fits use `LsqFit`. This includes fits without
active bounds, constraints, priors, or parameter-dependent covariance. Bounds of
the form `[-Inf, Inf]` are treated as no-op bounds and keep the fast path open.
ScientificFitting rejects an explicit `backend=:lsqfit` request when that backend would
discard any part of the cost or parameter controls; selecting a fast backend
must never change the statistical problem.

When `LsqFit` computes a weighted Jacobian, ScientificFitting reuses it during
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
ScientificFitting evaluates both mutating functions once at `p0` and rejects incomplete
or non-finite output buffers before solver dispatch. Keep their buffer and
parameter signatures generic rather than restricting them to `Float64`; AD-based
optimizer paths evaluate them with dual-number element types.

## Covariance Costs

Diagonal covariance stores precomputed inverse standard deviations and log
determinants. Dense static covariance stores a Cholesky factor and log
determinant. Sparse static `cov_y` matrices are kept sparse on the unbounded
least-squares path and whitened with CHOLMOD's permuted factor. Residuals are
whitened using linear solves, not explicit covariance inverses.

A `WhiteningOperator` stores neither covariance nor precision matrix. It applies
the user-supplied static whitening transformation directly to residuals and to
the columns of analytic Jacobians. This keeps the asymptotic cost equal to the
operator supplied by the user.

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

ScientificFitting validates that `x_derivative(x, p)` has the same length as `x` and
contains finite values. The derivative may depend on `p`; AD information is
preserved when the Gaussian likelihood cost needs gradients or Hessians.

## Matrix-Free Structured Whitening

Suppose a time series has stationary AR(1) covariance

```math
C_{ij}=\sigma^2\rho^{|i-j|}, \qquad |\rho|<1.
```

Materializing ``C`` costs ``O(n^2)`` memory even though its whitening operation
is a one-pass recurrence. If ``r`` is the residual vector, define

```math
z_1=\frac{r_1}{\sigma}, \qquad
z_i=\frac{r_i-\rho r_{i-1}}{\sigma\sqrt{1-\rho^2}}.
```

Then ``\lVert z\rVert^2=r^T C^{-1}r``. The corresponding determinant is

```math
\log\det C=2n\log\sigma+(n-1)\log(1-\rho^2).
```

The complete ScientificFitting contract is therefore compact:

```julia
sigma = 0.20
rho = 0.65
innovation_sigma = sigma * sqrt(1 - rho^2)

function whiten_ar1!(out, residual)
    out[1] = residual[1] / sigma
    @inbounds for i in 2:length(residual)
        out[i] = (residual[i] - rho * residual[i - 1]) / innovation_sigma
    end
    return nothing
end

n = length(y)
whitening = WhiteningOperator(
    whiten_ar1!;
    logdet_covariance=2n * log(sigma) + (n - 1) * log1p(-rho^2),
    marginal_sigma=sigma,
)

result = fit_model(model, x, y; p0=p0, whitening)
```

`logdet_covariance` is required even for a chi-square fit because ScientificFitting also
reports the normalized Gaussian ``-2\log L`` value, AIC, and BIC.
`marginal_sigma` is optional and affects only pointwise error bars and prediction
bands; it does not enter the cost. Without it, use `band=:confidence` rather
than claiming a prediction band whose observation variance is unknown.

The operator represents the **complete static observation covariance**. Do not
also pass `sigma_y`, `cov_y`, x uncertainties, or active error components.
ScientificFitting rejects those combinations because adding covariance models requires
an explicit scientific derivation. The function must accept `AbstractVector`
views and generic element types. That is necessary for Jacobian whitening and
for automatic differentiation on bounded or constrained optimizer paths.

ScientificFitting can validate dimensions, finite output, and interface compatibility;
it cannot prove that a custom transformation really satisfies
``W^T W=C^{-1}`` or that its supplied determinant is correct. Reference the
operator against a small dense covariance before using it at large scale.

## Known Limits

Dense covariance matrices are correct and tested, but they are not the right
representation for huge correlated datasets. They require `O(n^2)` memory and
`O(n^3)` factorization time. Sparse static `cov_y` supports unbounded
least-squares problems. `WhiteningOperator` supports complete static observation
covariance on both least-squares and AD-compatible constrained paths. Built-in
banded, Toeplitz, low-rank-plus-diagonal, sparse-precision, structured `cov_x`,
and parameter-dependent operator types remain future work. Dense-matrix
micro-optimization cannot fix the asymptotic scaling.

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
