# ScientificFitting

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/scientificfitting-logo-dark.svg">
  <img alt="ScientificFitting" src="docs/src/assets/scientificfitting-logo.svg" width="560">
</picture>

[![CI](https://github.com/Amin-El-Sayed/ScientificFitting.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/Amin-El-Sayed/ScientificFitting.jl/actions/workflows/ci.yml)
[![Documentation](https://github.com/Amin-El-Sayed/ScientificFitting.jl/actions/workflows/pages.yml/badge.svg)](https://amin-el-sayed.github.io/ScientificFitting.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Scientific fitting for Julia: data, uncertainties, statistically explicit
results, diagnostics, and Makie plots from one coherent workflow.

ScientificFitting combines weighted nonlinear least squares, likelihood fits,
parameter constraints, profile diagnostics, and publication-ready Makie plots.
Version 0.1 is the first public release; feedback from real analyses is welcome.

## Install

```julia
using Pkg
Pkg.add("ScientificFitting")
```

Add `CairoMakie` separately when plots are needed; fitting and reporting do not
load Makie.

## Quickstart

```julia
using ScientificFitting
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

The same fit without plotting is simply:

```julia
using ScientificFitting

result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y)
report = fit_report(result; parameter_names=["m", "b"])
```

## Documentation

Use the [online documentation](https://amin-el-sayed.github.io/ScientificFitting.jl/)
for installation details, complete scientific examples, statistical methods,
diagnostics, plotting, performance guidance, and the API reference.

Build the local documentation from the repository root with:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Then serve `docs/build` with any static file server.

## Technical Notes

- Julia 1.10 and later are supported. CI tests Julia 1.10 and Julia 1.12.
- Gaussian fits accept diagonal, dense, sparse-static, component-based, and
  custom matrix-free static whitening uncertainty models, including x/y
  uncertainty propagation.
- Likelihood workflows cover Poisson, histogram, unbinned,
  extended-unbinned, indexed, custom-objective, and multi-dataset fits.
- Dense covariance costs `O(n^2)` memory and `O(n^3)` factorization; use a
  `WhiteningOperator` when a large structured problem has a fast whitening
  operation.
- Parameter covariance is a local Hessian approximation. Use profiles and
  contours for nonlinear, weakly constrained, bounded, or asymmetric cases.
- CairoMakie is an optional package extension. The numerical core, reports,
  and diagnostics remain usable without it.

## Contributing

Bug reports, scientific examples, documentation corrections, and focused pull
requests are welcome. Include a minimal executable example, package versions,
the uncertainty model, complete output, and the expected behavior. Numerical
changes need an analytic or independent reference; performance claims need a
reproducible benchmark on named hardware.

## License And Citation

ScientificFitting is licensed under the [MIT License](LICENSE). If it contributes to
published research, please cite the exact version using [CITATION.cff](CITATION.cff).
The documentation explains [citation, related work, and method-specific
attribution](docs/src/citation.md).
