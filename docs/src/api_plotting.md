# Fit Plotting

ScientificFitting's numerical core does not depend on Makie. Fitting, profiles,
diagnostics, and text reports do not require Makie; the plotting methods are
activated by CairoMakie:

```julia
using ScientificFitting
using CairoMakie
```

Without `using CairoMakie`, every plotting entry point raises an
`ArgumentError` that names the missing optional extension. This separation is
intentional: batch fits and server-side reports do not pay Makie's load or
compilation cost.

## Choose The Plotting Entry Point

| Task | Function | Return value | Runs an optimizer? |
|---|---|---|---:|
| Fit arrays and plot immediately | [`fitplot`](@ref) | `(result, figure)` | yes |
| Plot an existing x-y fit | [`plot_fit`](@ref) | `Figure` | no |
| Add content to the data axis | [`fit_axis`](@ref), `add_*!` | `Axis` or Makie plot object | no |
| Compose a custom themed figure | [`plot_theme`](@ref), [`plot_palette`](@ref), [`plot_info_panel!`](@ref) | `Theme`, style tokens, or `GridLayout` | no |

The [Plotting And Customization](plotting_design.md) guide develops complete
composition examples. This page is the exact argument and failure contract.

## Fit And Plot In One Call

```text
fitplot(model, x, y; p0, show_panel=true, print_report=false, kwargs...)
fitplot(x, y; p0=nothing, show_panel=true, print_report=false, kwargs...)
fitplot(result::FitResult; show_panel=true, print_report=false, kwargs...)
```

All methods return the named tuple `(result=result, figure=figure)`. The
two-array method fits the straight line ``y=p_1x+p_2`` and derives an initial
slope and intercept from the first and last observations unless `p0` is given.
The `FitResult` method only renders; it never repeats the fit.

### Output selection

`show_panel::Bool` and `print_report::Bool` are independent. The former controls
the right/inside numerical panel; the latter prints `report_text(result)` to the
terminal and belongs only to `fitplot`. Both styles support either panel state.

### Keyword routing

For methods that perform a fit, the following keywords are sent to
[`fit_model`](@ref):

| Concern | Fitting keywords |
|---|---|
| Observation uncertainty | `sigma_y`, `sigma_x`, `cov_y`, `cov_x`, `whitening`, `error_components` |
| Parameter information | `bounds`, `constraints`, `parameter_priors`, `parameter_constraints`, `fixed_parameters` |
| Derivatives and model evaluation | `jacobian`, `x_derivative`, `inplace` |
| Solver and covariance behavior | `backend`, `cost`, `maxiters`, `tol`, `scale_covariance`, `initial_guesses`, `multistart` |

All remaining keywords are sent to [`plot_fit`](@ref). A misspelled fitting
keyword therefore does not disappear silently: it reaches `plot_fit` and fails
as an unsupported keyword.

## Plot An Existing Result

```text
plot_fit(result::FitResult; kwargs...) -> Figure
```

`plot_fit` draws the observations, available x/y error bars, fitted model,
optional uncertainty band, and optional result panel. It does not modify
`result` or rerun the optimizer.

### Output, style, and appearance

| Keyword | Default | Contract |
|---|---:|---|
| `filename` | `nothing` | Save during construction when a path is supplied. |
| `format` | `:pdf` | Extension appended only when `filename` has no extension. |
| `theme` | `:sans` | Maintained visual style: `:sans` or `:tex`. |
| `appearance` | `:auto` | `:light`, `:dark`, or `:auto`; `:auto` currently resolves to light. |
| `theme_override` | `Theme()` | Makie theme merged after the selected ScientificFitting style. |
| `figure_size` | style/panel-dependent | Logical Makie canvas `(width, height)`; it does not set raster density. |
| `tight_layout` | `true` | Trim layout whitespace while preserving the declared figure footprint. |

The former names `:analysis`, `:presentation`, `:screen`, `:lab`, `:workbench`,
`:modern`, `:clean`, `:minimal`, and `:showcase` map to `:sans`; `:article`,
`:publication`, `:paper`, and `:latex` map to `:tex`. New code should use the
two maintained names. Unknown styles and unknown appearances raise
`ArgumentError`. Select dark output with `appearance=:dark`; visual style and
color appearance are separate arguments.

### Labels, units, model domain, and limits

| Keyword | Default | Contract |
|---|---:|---|
| `title` | `nothing` | Figure title; `nothing` produces no title. |
| `model_label` | automatic for the built-in line | Model expression shown in the right panel. |
| `xlabel`, `ylabel` | `"x"`, `"y"` | Axis quantity labels. |
| `xunit`, `yunit` | `nothing` | Appended in SI quantity-calculus form as `label / unit`; units are never inferred. |
| `latex_labels` | `false` | Convert suitable labels to Makie `LaTeXString` content. Pass explicit `L"..."` strings for mathematical notation. |
| `xgrid` | `nothing` | Explicit finite model-sampling coordinates; authoritative when supplied. |
| `fit_range` | `:axis` | `:axis` samples over padded visible x limits; `:data` samples from the first to last measured x. |
| `auto_limits` | `true` | Include data, errors, model, and displayed band in both axis limits. |
| `limit_padding` | `0.08` | Finite non-negative fractional padding around automatic content limits. |
| `plot_aspect` | `nothing` | Optional numeric `AxisAspect`; leave unset unless geometry carries meaning. |
| `axis_kwargs` | `NamedTuple()` | Makie `Axis` attributes applied after ScientificFitting's title/label defaults. |

`fit_range` must be `:axis` or `:data`. `limit_padding` must be finite and
non-negative. If `auto_limits=false`, provide limits through `axis_kwargs` or
set them on `fit_axis(figure)` after construction. Manual axis limits do not
resample the model; use a matching `xgrid` when extrapolation is intentional.

### Uncertainty band

| Keyword | Default | Contract |
|---|---:|---|
| `band` | `:confidence` | `:confidence`, `:prediction`, or `:none`. |
| `nsigma` | `1.0` | Finite positive multiplier for the displayed standard-deviation scale. |
| `band_label` | `"1-sigma band"` | Legend text; update it whenever `nsigma` or the band meaning changes. |
| `band_color`, `band_alpha` | style-dependent | Direct ScientificFitting-level overrides. |
| `band_kwargs` | `NamedTuple()` | Makie `band!` attributes applied last. |

`band=:confidence` propagates the local parameter covariance to the fitted
mean. `band=:prediction` adds pointwise observation uncertainty in y and the
effective contribution from x uncertainty. It is still a local covariance
construction, not an exact nonlinear coverage statement.

A matrix-free [`WhiteningOperator`](@ref) can define the fit without exposing
pointwise marginal errors. In that case `band=:prediction` requires
`marginal_sigma`; otherwise it raises `ArgumentError`. `band=:confidence`
remains available.

### Result panel and legend

| Keyword | Default | Contract |
|---|---:|---|
| `show_panel` | `true` | Show the structured right panel or compact in-axis panel. Independent of visual style. |
| `stats_position` | `:right` | `:right` or `:inside`. |
| `inside_stats_position` | `:lt` | `:lt`, `:lb`, `:rt`, `:rb` and their long aliases. |
| `stats_panel_width` | `:auto` | Natural Makie width, a fraction `0 < w <= 1`, or a positive pixel width. Fractions are clamped to 300--560 px. |
| `panel_gap` | style-dependent | Gap between data axis and right panel. |
| `stats_mode` | `:compact` | `:compact` or `:full`. |
| `stats_sigdigits` | `5` | Significant digits used only for displayed values. |
| `parameter_names` | `nothing` | Display names; length must equal the number of fitted parameters. |
| `stats_fontsize` | style-dependent | Explicit result-panel or in-axis text size. |
| `stats_title` | `nothing` | Optional title above the structured right panel. |
| `latex_stats` | `false` | Render structured right-panel symbols and numbers as LaTeX. |
| `stats_box_color`, `stats_box_alpha` | style-dependent, `0.95` | In-axis summary background. |
| `stats_box_strokecolor`, `stats_box_strokewidth` | style-dependent, `1.0` | In-axis summary border. |
| `show_legend` | `true` | Show data, fit, and band labels. With a right panel, the legend is placed above the report. |
| `legend_position` | `:rt` | In-axis Makie legend position when no right-side panel owns the legend. |
| `legend_kwargs` | `NamedTuple()` | Makie legend attributes applied last. |

In the right panel, `stats_mode=:full` adds cost, AIC, and BIC to the compact
parameter, chi-square, p-value, and ndf rows. In the compact in-axis box,
`:full` adds the raw chi-square. AIC and BIC are displayed values, not a license
to compare fits with different data or incompatible likelihood normalization;
see [Results And Diagnostics](api_results.md#Results-And-Diagnostics).

### Data, fit, and error-bar styling

Role defaults are used whenever a scalar keyword is `nothing`. The associated
Makie keyword container is merged last and therefore has final authority.

| Layer | ScientificFitting-level keywords | Final Makie container |
|---|---|---|
| Observations | `data_color`, `data_marker`, `data_markersize`, `data_strokecolor`, `data_strokewidth`, `data_label` | `scatter_kwargs` |
| Fit curve | `fit_color`, `fit_linewidth`, `fit_label` | `line_kwargs` |
| Band | `band_color`, `band_alpha`, `band_label` | `band_kwargs` |
| X errors | `xerr_color`, `error_whiskerwidth` | `xerrorbars_kwargs` |
| Y errors | `yerr_color`, `error_whiskerwidth` | `yerrorbars_kwargs` |

Every `*_kwargs` container accepts a `NamedTuple`, `AbstractDict`, or `nothing`.
Other container types raise `ArgumentError`. Explicit element overrides affect
only that layer; unmodified layers continue to follow the selected style. Both
styles render measurement-error lines at full contrast with the same fine
stroke. Override `linewidth` only through the corresponding Makie container.

## Extend A Finished Figure

```julia
fig = plot_fit(result; show_legend=false)
ax = fit_axis(fig)

add_vline!(ax, threshold; color=:black, linestyle=:dash)
add_vband!(ax, threshold_low, threshold_high; color=(:gray50, 0.15))
add_points!(ax, [derived_x], [derived_y]; marker=:star5)
```

| Function | Arguments and defaults | Return value | Validation |
|---|---|---|---|
| `fit_axis(figure; index=1)` | One-based axis index | Makie `Axis` | Index must exist. |
| `add_curve!(axis, f; xgrid=nothing, xspan=nothing, n=400, label=nothing, kwargs...)` | Sample on `xgrid`, `xspan`, or current x limits | Makie line plot | At least two finite x values, finite curve values, `n >= 2`. |
| `add_curve!(axis, x, y; label=nothing, kwargs...)` | Precomputed curve | Makie line plot | Equal-length finite vectors with at least two points. |
| `add_points!(axis, x, y; label=nothing, kwargs...)` | Scalar or vector coordinates | Makie scatter plot | Equal-length finite coordinates. |
| `add_vline!`, `add_hline!` | Scalar or vector coordinate, optional `label` | Makie line collection | Coordinates must be finite. |
| `add_vband!(axis, xmin, xmax; label=nothing, kwargs...)` | Ordered finite x bounds | Makie axis-relative span | Requires `xmin <= xmax`. |
| `add_hband!(axis, ymin, ymax; label=nothing, kwargs...)` | Ordered finite y bounds | Makie axis-relative span | Requires `ymin <= ymax`. |

Axis-relative bands do not inject artificial values into the orthogonal data
limits. None of these helpers changes or reruns the fit. Their `kwargs...` are
ordinary Makie plot attributes.

## Reuse The Visual Contract

```julia
theme = plot_theme(:sans; appearance=:dark)
style = plot_palette(:sans; appearance=:dark)
```

`plot_theme(style; appearance, theme_override)` returns the Makie `Theme` used
by ScientificFitting. `plot_palette(style; appearance)` returns the corresponding named
tuple of visual tokens: data/fit/band colors, multi-series colors, marker and
line sizes, typography, grids and spines, report-panel spacing, and default
figure sizes. These functions let a custom Makie layout inherit the same style
without copying private constants.

`plot_info_panel!` adds the same left-aligned information hierarchy used by
`plot_fit`:

```julia
plot_info_panel!(
    fig[1, 2];
    theme=:sans,
    appearance=:light,
    legend_plots=[data_plot, fit_plot],
    legend_labels=["data", "fit"],
    title="Fit summary",
    model_label="damped oscillator",
    parameter_lines=["A = ...", "lambda = ..."],
    statistic_lines=["chi2/ndf = ..."],
)
```

The alternative `legend_source=axis` builds the legend from labeled content on
an axis. `fontsize`, `color`, `muted_color`, and `legend_kwargs` override style
defaults. `tellwidth=true` lets the panel report its natural width;
`tellheight=false` prevents a short report from shrinking or vertically
centering the adjacent scientific axis. The function returns its `GridLayout`.

## Diagnostic Figures

Residual, pull, ratio, profile, contour, and profile-matrix figures are listed
separately in [Diagnostic Plotting](api_plotting_diagnostics.md). This page
covers fit figures and reusable Makie composition only.

## Export Semantics

Every high-level plot function returns its `Figure` even when `filename` is
provided. A filename extension takes precedence over `format`; with no
extension, `.$(format)` is appended. For explicit resolution control, save the
returned figure with Makie:

```julia
fig = plot_fit(result; theme=:tex)
save("fit.svg", fig)
save("fit.png", fig; px_per_unit=2)
```

`figure_size` controls layout size. `px_per_unit` controls raster density.
Increasing the former and scaling the image down later also scales down its
text; it is not a substitute for export resolution.

## Failure Summary

| Failure | Result |
|---|---|
| CairoMakie extension not loaded | `ArgumentError` naming the required extension |
| Invalid style, appearance, band, stats position, or fit range | `ArgumentError` |
| Non-positive/non-finite `nsigma`, negative/non-finite `limit_padding` | `ArgumentError` |
| Prediction band without matrix-free marginal errors | `ArgumentError` with the required remedy |
| Wrong number of `parameter_names` | `ArgumentError` |
| Non-finite or dimensionally inconsistent annotation data | `ArgumentError` |

## API Documentation

```@docs
ScientificFitting.fitplot
ScientificFitting.plot_fit
ScientificFitting.fit_axis
ScientificFitting.add_curve!
ScientificFitting.add_points!
ScientificFitting.add_vline!
ScientificFitting.add_hline!
ScientificFitting.add_vband!
ScientificFitting.add_hband!
ScientificFitting.plot_theme
ScientificFitting.plot_palette
ScientificFitting.plot_info_panel!
```
