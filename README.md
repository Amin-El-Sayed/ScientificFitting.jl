# JuFitter

Scientific fitting for Julia: data, uncertainties, statistically explicit
results, diagnostics, and Makie plots from one coherent workflow.

JuFitter is built for laboratory, engineering, and scientific analysis where the
fit result is only useful if the uncertainty model, diagnostics, and plot are
understandable. The package is currently pre-release; the repository is being
hardened before broad public promotion or Julia package registration.

## Why JuFitter

- Weighted nonlinear least squares with diagonal, full-covariance, x/y, and
  component-based uncertainty models.
- Poisson, histogram, unbinned, extended-unbinned, indexed, custom-objective,
  and multi-dataset likelihood workflows.
- Bounds, fixed parameters, Gaussian priors, correlated parameter constraints,
  profile intervals, and two-parameter contours.
- Structured fit reports and diagnostic dashboards that tell you what to
  inspect next.
- Optional CairoMakie plotting with robust default layouts, right-side reports,
  residual/pull diagnostics, profile/contour plots, and post-fit annotation
  helpers.
- A documentation-first workflow: every serious feature should have a tested
  explanation, not only an exported function.

## Quick Start

The fitting and reporting core does not load Makie. Load CairoMakie only when
you want plots.

```julia
using JuFitter
using CairoMakie

x = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
y = [0.92, 1.80, 2.77, 3.60, 4.55, 5.49, 6.40, 7.37, 8.19, 9.18]
sigma_y = fill(0.04, length(x))

model(x, p) = @. p[1] * x + p[2]

fit = fitplot(
    model,
    x,
    y;
    p0=[1.0, 0.0],
    sigma_y=sigma_y,
    parameter_names=["m", "b"],
    xlabel="position x",
    xunit="mm",
    ylabel="voltage U",
    yunit="V",
    filename="calibration_fit.pdf",
)

result = fit.result
println(report_text(result; parameter_names=["m", "b"]))
println(diagnostic_dashboard_text(result))
```

For fitting without plotting:

```julia
using JuFitter

result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y)
report = fit_report(result; parameter_names=["m", "b"])
```

## Installation

JuFitter supports Julia 1.10 and later. The release gate tests Julia 1.10 and
the current Julia 1.12 release explicitly.

During active development, use a checked-out repository:

```julia
using Pkg
Pkg.activate("/path/to/JuFitter")
Pkg.instantiate()
```

From the repository root:

```bash
julia --project=.
```

After registration, the intended user-facing path is:

```julia
using Pkg
Pkg.add("JuFitter")
```

For plotting, add and load CairoMakie in the same environment:

```julia
using Pkg
Pkg.add(["JuFitter", "CairoMakie"])

using JuFitter
using CairoMakie
```

## Documentation

The documentation is the main entry point:

- `docs/src/quickstart.md` for the first complete fit.
- `docs/src/gallery.md` for complete scientific examples.
- `docs/src/fitting_for_practitioners.md` for practical fit judgement.
- `docs/src/statistical_foundations.md` for the statistical assumptions and
  formulas.
- `docs/src/api.md` for the generated API reference.
- `docs/src/maintenance.md` for extension points, bottlenecks, and release
  rules.

Build the local documentation from the repository root with:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Then serve `docs/build` with any static file server.

## Current Scope

JuFitter is intended to cover the common scientific fitting workflows before
v1:

- Gaussian XY fits with diagonal or dense covariance.
- Effective x-uncertainty propagation, including a vectorized derivative hook
  for large datasets.
- Likelihood fits for counts, histograms, unbinned samples, and custom
  objectives.
- Profile and contour diagnostics for cases where local covariance is not
  enough.
- Makie plots that can be used as-is or extended through returned Makie figures
  and JuFitter annotation helpers.

Known limitations are explicit:

- Dense covariance is exact but expensive: `O(n^2)` memory and `O(n^3)`
  factorization. Static sparse `cov_y` is supported for unbounded least-squares
  fits, but large correlated time series, spectra, images, and detector arrays
  still need future structured covariance or whitening operators.
- Parameter covariance is a local Hessian approximation. Nonlinear models,
  weak data, active bounds, and asymmetric likelihoods should be checked with
  profiles or contours.
- CairoMakie has a noticeable first-use compilation cost. Fitting and reporting
  remain usable without loading Makie.
- Python interoperability is planned through JuliaCall/PythonCall and has an
  opt-in gate, but it is not a release claim yet.

## Development Gates

Fast checks:

```bash
julia --project=. --startup-file=no test/docs_public_release_gate.jl
julia --project=. --startup-file=no test/docs_gallery_gate.jl
julia --project=. --startup-file=no test/docs_link_gate.jl
```

Core and package checks:

```bash
julia --project=. --startup-file=no -e 'include("test/core_runtests.jl")'
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
```

Plot checks:

```bash
julia --project=docs --startup-file=no test/plots/fitplot.jl
```

Benchmarks:

```bash
julia --project=benchmarks benchmarks/runbenchmarks.jl --seconds=1
```

Benchmark output and generated plots are ignored by git.

## Release Policy

Do not push, publish, register, deploy documentation, or make the repository
public without explicit manual approval from the maintainer. Local `codex/*`
branch commits are review checkpoints, not release actions.
