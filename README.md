# JuFitter

JuFitter is a Julia fitting utility focused on fast scientific workflows:

- weighted nonlinear least-squares fitting
- optional full covariance for `x` and `y`
- constraint-aware fitting (hybrid backend)
- publication-ready CairoMakie plots with one-call `fitplot(...)` workflows

## Install (in this repo)

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Quick start

```julia
using JuFitter

x = collect(range(0.0, 10.0; length=200))
model(x, p) = @. p[1] * x + p[2]
sigma_y = fill(0.2, length(x))
y = model(x, [2.0, 1.0]) .+ sigma_y .* sin.(1.8 .* x)

fit = fitplot(
    model,
    x,
    y;
    p0=[1.0, 0.0],
    sigma_y=sigma_y,
    xlabel="x",
    ylabel="y",
    parameter_names=["m", "b"],
    filename="fit.pdf",
)

result = fit.result

@show result.params
@show result.param_stderr
@show result.stats.cost result.stats.cost_min result.stats.nll_min
@show result.stats.chi2 result.stats.chi2_ndf result.stats.ndf
```

## API

- `FitProblem(model, x, y; p0, sigma_y, sigma_x, cov_y, cov_x, error_components, bounds, constraints, parameter_priors, parameter_constraints, fixed_parameters, jacobian)`
- `fit(problem::FitProblem; backend=:auto, cost=:auto, maxiters=500, tol=1e-10, ci_level=0.6827, scale_covariance=:auto, initial_guesses=nothing, multistart=1)`
- `fit_model(model, x, y; kwargs...)`
- `fit_poisson_model(model, x, counts; p0, kwargs...)`
- `fit_histogram_model(expected_counts, edges, counts; p0, kwargs...)`
- `fit_histogram_density(pdf, edges, counts; p0, total_count=sum(counts), kwargs...)`
- `fit_unbinned_model(pdf, data; p0, kwargs...)`
- `fit_extended_unbinned_model(rate, data, domain; p0, kwargs...)`
- `fit_indexed_model(model, indices, y; p0, sigma_y=nothing, cov_y=nothing, kwargs...)`
- `fit_custom(objective; p0, nobs, gof=nothing, kwargs...)`
- `fit_multi_model(models, xs, ys; p0, sigma_y=nothing, kwargs...)`
- `fitplot(x, y; sigma_y=nothing, kwargs...)`
- `fitplot(model, x, y; p0, sigma_y=nothing, kwargs...)`
- `plot_fit(result; xgrid=nothing, filename=nothing, format=:pdf, theme=:clean, ...)`
- `plot_residuals(result; kind=:pull, filename=nothing, format=:pdf, ...)`
- `plot_diagnostics(result; filename=nothing, format=:pdf, ...)`
- `fit_report(result; parameter_names=nothing)`
- `report_text(result; parameter_names=nothing)`

Key `plot_fit` customization kwargs:

- labels/text:
  - `title`, `xlabel`, `ylabel`, `xunit`, `yunit`
  - `model_label` for the model formula shown above the fit summary
  - `latex_labels=true` to render string labels/titles as LaTeX
  - `parameter_names` (names in the side summary panel; can be `LaTeXString`, e.g. `L"\\lambda"`)
- plot size/layout:
  - `figure_size=(width, height)` in pixels (e.g. `(1800, 900)`)
  - `auto_limits=true` pads limits around data, error bars, fit curve, and confidence band
  - `limit_padding=0.08` controls automatic axis padding
  - `stats_panel_width` controls the right summary panel; values `<= 1` are relative, values `> 1` are pixels
  - `panel_gap` controls the gap between plot and summary panel
  - `plot_aspect` controls the axis width/height ratio (e.g. `1.0` for square plot area)
- summary display:
  - `show_stats`, `stats_sigdigits`, `stats_fontsize`, `latex_stats`
  - `stats_title` for an optional small summary title
  - `stats_position=:right` shows a compact side summary without covering the data
  - `stats_position=:inside` draws an in-axis box for compact slide-style plots
  - shows the model formula, fitted parameters ± uncertainties, `\chi^2`, `\chi^2/\mathrm{ndf}`, `P(\chi^2)`, and `\mathrm{ndf}`
  - `stats_box_color`, `stats_box_alpha`, `stats_box_strokecolor`, and `stats_box_strokewidth` control only the inside box
- styling:
  - `data_color`, `data_marker`, `data_markersize`
  - `fit_color`, `fit_linewidth`
  - `band=:confidence` or `band=:none`
  - `nsigma`, `band_color`, `band_alpha`
  - `xerr_color`, `yerr_color`, `error_whiskerwidth`
  - `show_legend`, `legend_position`
- theme/font control:
  - `theme` (`:clean`, `:minimal`, `:paper`, `:publication`, `:latex`, or custom base)
  - `theme=:minimal` for dense datasets with fine markers and precise black/white styling
  - `theme=:paper` for LaTeX-like physics publication styling
  - `theme_override=Theme(...)` for global style overrides, including fonts/font sizes
- direct Makie keyword forwarding:
  - `axis_kwargs`, `legend_kwargs`, `line_kwargs`, `scatter_kwargs`, `band_kwargs`
  - `xerrorbars_kwargs`, `yerrorbars_kwargs`
  - `stats_box_kwargs`, `stats_label_kwargs`, `stats_title_kwargs`
  - each accepts a `NamedTuple` or `Dict`, e.g. `axis_kwargs=(xscale=log10,)`

## Fit reports

Use `fit_report` when code should extract values programmatically:

```julia
report = fit_report(result; parameter_names=["m", "b"])
report.parameters[1].value
report.parameters[1].uncertainty
report.statistics.cost
report.statistics.cost_min
report.statistics.nll_min
report.statistics.chi2_ndf
report.diagnostics.warnings
```

Use `report_text` for a readable text summary:

```julia
println(report_text(result; parameter_names=["m", "b"]))
```

For nonlinear problems you can request profile-based asymmetric uncertainties:

```julia
report = fit_report(result; parameter_names=["m", "b"], errors=:profile)
interval = profile_interval(result, 1)
```

## Cost functions

JuFitter separates the statistical cost function from the numerical solver.
Use `cost=:chi2` for classical weighted least squares and `cost=:gaussian_nll`
for the full Gaussian negative log-likelihood. With `cost=:auto`, JuFitter uses
the full Gaussian NLL when the effective covariance is parameter-dependent, for
example when `sigma_x` or `cov_x` is provided.

See `docs/src/backend_design.md` and `docs/src/statistical_foundations.md` for
the mathematical details.

## Error components

Use `error_components` to define named, individually switchable uncertainty
sources. These are combined with `sigma_y`, `sigma_x`, `cov_y`, and `cov_x`.

```julia
result = fit_model(
    model,
    x,
    y;
    p0=[1.0, 0.0],
    error_components=[
        (name=:stat, target=:y, mode=:absolute, values=fill(0.1, length(y))),
        (name=:scale, target=:y, mode=:relative, values=0.02),
        (name=:model_scale, target=:y, mode=:model_relative, values=0.02),
        (name=:disabled_systematic, target=:y, mode=:absolute, values=1.0, active=false),
    ],
)
```

Supported targets are `:y` and `:x`. Supported modes are `:absolute`,
`:relative`, `:model_relative` for y only, and `:covariance`.

## Likelihood fits

Poisson, histogram, unbinned, extended-unbinned, indexed, custom objective, and simultaneous multi-dataset fits use the same
parameter-control layer as XY fits: bounds, constraints, priors, correlated
parameter constraints, fixed parameters, multistart, reports, and profiles.

```julia
count_model(x, p) = @. exp(p[1] + p[2] * x)
poisson_result = fit_poisson_model(count_model, x, counts; p0=[0.0, 0.1])

expected_counts(edges, p) = [p[1] * (edges[i + 1] - edges[i]) for i in 1:(length(edges) - 1)]
hist_result = fit_histogram_model(expected_counts, edges, counts; p0=[1.0])

normal_pdf(x, p) = exp(-0.5 * ((x - p[1]) / p[2])^2) / (p[2] * sqrt(2pi))
unbinned_result = fit_unbinned_model(normal_pdf, data; p0=[0.0, 1.0])

rate(x, p) = exp(p[1])
extended_result = fit_extended_unbinned_model(rate, data, (0.0, 1.0); p0=[0.0])

custom_result = fit_custom(p -> sum(abs2, p .- [1.0, 2.0]); p0=[0.0, 0.0], nobs=4)
```

## Fixed parameters, profiles, contours

Fixed parameters are removed from the optimized parameter vector:

```julia
result = fit_model(
    model,
    x,
    y;
    p0=[1.0, 0.2],
    sigma_y=sigma_y,
    fixed_parameters=(index=2, value=0.2, sigma_minus=0.005, sigma_plus=0.008),
)
```

The value is fixed during optimization. The optional uncertainty is reported
and included as a local fixed-parameter uncertainty in the covariance matrix.

Correlated Gaussian parameter constraints are separate from fixed parameters:

```julia
result = fit_model(
    model,
    x,
    y;
    p0=[1.0, 0.0],
    sigma_y=sigma_y,
    parameter_constraints=(
        indices=[1, 2],
        mean=[1.0, 0.0],
        covariance=[0.01 0.002; 0.002 0.04],
    ),
)
```

Profile and contour scans re-minimize the cost function with selected
parameters fixed:

```julia
prof = profile(result, 1; npoints=61, nsigma=3)
int = profile_interval(result, 1; npoints=121, nsigma=5)
cont = contour(result, 1, 2; npoints=31, nsigma=3)

plot_profile(prof; filename="profile.pdf")
plot_contour(cont; filename="contour.pdf")
```

## Robust starts

For hard nonlinear problems, provide multiple initial guesses or ask JuFitter
to generate a small deterministic candidate set:

```julia
result = fit_model(
    model,
    x,
    y;
    p0=[0.1, 5.0],
    bounds=([0.0, 0.0], [10.0, 10.0]),
    initial_guesses=[[0.1, 5.0], [2.0, 0.5]],
    multistart=2,
)
```

Covariance scaling is controlled explicitly:

```julia
scale_covariance=:auto    # default
scale_covariance=:never
scale_covariance=:always
```

## Diagnostics

Every result contains `result.diagnostics`:

```julia
result.diagnostics.warnings
result.diagnostics.covariance_condition
result.diagnostics.hessian_condition
result.diagnostics.active_bounds
```

Warnings are intentionally conservative. Active bounds, non-positive degrees of
freedom, non-convergence, and ill-conditioned covariance/Hessian estimates are
reported because local uncertainties and p-values can become unreliable in
those cases.

For XY fits, diagnostic plots are available:

```julia
plot_residuals(result; kind=:pull, filename="pulls.pdf")
plot_residuals(result; kind=:residual, filename="residuals.pdf")
plot_residuals(result; kind=:ratio, filename="ratio.pdf")
plot_diagnostics(result; filename="fit_diagnostics.pdf")
```

## Examples

Run examples from the repository root with:

```julia
julia --project=. examples/gallery/01_quickstart_linear.jl
```

Generated plots are written to `examples/output/`, which is intentionally
ignored by git.

Available examples are documented in `examples/README.md` and are organized as
a numbered gallery:

- `examples/gallery/01_quickstart_linear.jl`
- `examples/gallery/02_xy_uncertainties_photoelectric.jl`
- `examples/gallery/03_plot_customization.jl`
- `examples/gallery/04_covariance_and_effective_variance.jl`
- `examples/gallery/05_constraints_priors_profiles.jl`
- `examples/gallery/06_likelihood_workflows.jl`

## Constraint format

`constraints` can be a named tuple:

```julia
constraints = (
    ineq = p -> [p[1] - 3.0],  # <= 0
    eq = p -> [p[2] - 1.0],    # == 0
)
```

## Notes on scale

Dense covariance matrices are supported but memory-heavy for very large datasets.
They need `O(N^2)` memory and `O(N^3)` factorization time. For large `N`, prefer
diagonal uncertainties today; truly large correlated problems need future
structured covariance operators rather than materialized dense matrices.

Performance benchmarks live in `benchmarks/runbenchmarks.jl` and can be run
with:

```bash
julia --project=. benchmarks/runbenchmarks.jl
```
