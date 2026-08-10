# Plotting

JuFitter's numerical core does not depend on Makie. Fitting, profiles,
diagnostics, and text reports do not require Makie; the plotting methods are
activated by CairoMakie:

```julia
using JuFitter
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
| Inspect residuals, pulls, or ratios | [`plot_residuals`](@ref), [`plot_diagnostics`](@ref) | `Figure` | no |
| Render an existing profile or contour | [`plot_profile`](@ref), [`plot_contour`](@ref) | `Figure` | no |
| Compute and render a profile matrix | [`plot_profile_matrix`](@ref)`(result)` | `Figure` | yes, repeatedly |
| Render a stored profile matrix | [`plot_profile_matrix`](@ref)`(matrix)` | `Figure` | no |
| Compose a custom themed figure | [`plot_theme`](@ref), [`plot_palette`](@ref), [`plot_info_panel!`](@ref) | `Theme`, style tokens, or `GridLayout` | no |

The [Plotting And Customization](plotting_design.md) guide develops complete
composition examples. This page is the exact argument and failure contract.

## Fit And Plot In One Call

```julia
fitplot(model, x, y; p0, report=:plot, kwargs...)
fitplot(x, y; p0=nothing, report=:plot, kwargs...)
fitplot(result::FitResult; report=:plot, kwargs...)
```

All methods return the named tuple `(result=result, figure=figure)`. The
two-array method fits the straight line ``y=p_1x+p_2`` and derives an initial
slope and intercept from the first and last observations unless `p0` is given.
The `FitResult` method only renders; it never repeats the fit.

```julia
x = [0.0, 1.0, 2.0, 3.0]
y = [0.14, 1.08, 2.16, 3.03]
sigma_y = [0.10, 0.12, 0.11, 0.13]

out = fitplot(
    x,
    y;
    sigma_y=sigma_y,
    xlabel="position",
    xunit="mm",
    ylabel="voltage",
    yunit="V",
    report=:plot,
)

result = out.result
figure = out.figure
```

### Output selection

`report` belongs to `fitplot`, not `plot_fit`:

| `report` | Right/inside result panel | `report_text(result)` in terminal |
|---|---:|---:|
| `:plot` | yes | no |
| `:console` | no | yes |
| `:both` | yes | yes |
| `:none` | no | no |

An explicit `show_stats=true` or `false` overrides the panel default selected
by `report`. Any other `report` value raises `ArgumentError`.
For `plot_fit`, use `show_stats` directly because no terminal report is emitted.
Without an override, `:screen` includes the analysis panel and `:article`
reserves the canvas for the figure.

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

```julia
plot_fit(result::FitResult; kwargs...) -> Figure
```

`plot_fit` draws the observations, available x/y error bars, fitted model,
optional uncertainty band, and optional result panel. It does not modify
`result` or rerun the optimizer.

### Output, role, and appearance

| Keyword | Default | Contract |
|---|---:|---|
| `filename` | `nothing` | Save during construction when a path is supplied. |
| `format` | `:pdf` | Extension appended only when `filename` has no extension. |
| `theme` | `:screen` | Maintained role: `:screen` or `:article`. |
| `appearance` | `:auto` | `:light`, `:dark`, or `:auto`; `:auto` currently resolves to light. |
| `theme_override` | `Theme()` | Makie theme merged after the selected JuFitter role. |
| `figure_size` | role-dependent | Logical Makie canvas `(width, height)`; it does not set raster density. |
| `tight_layout` | `true` | Trim layout whitespace while preserving the declared figure footprint. |

The compatibility aliases are `:lab`, `:workbench`, `:modern`, `:clean`,
`:minimal`, and `:showcase => :screen`, and `:publication`, `:paper`, and
`:latex => :article`. New code should use only the two maintained names.
Unknown roles, unknown appearances, and contradictory
`theme=:dark, appearance=:light` input raise `ArgumentError`.

### Labels, units, model domain, and limits

| Keyword | Default | Contract |
|---|---:|---|
| `title` | `nothing` | Figure title; `nothing` produces no title. |
| `model_label` | automatic for the built-in line | Model expression shown in the right panel. |
| `xlabel`, `ylabel` | `"x"`, `"y"` | Axis quantity labels. |
| `xunit`, `yunit` | `nothing` | Appended as `label (unit)`; units are never inferred. |
| `latex_labels` | `false` | Convert suitable labels to Makie `LaTeXString` content. Pass explicit `L"..."` strings for mathematical notation. |
| `xgrid` | `nothing` | Explicit finite model-sampling coordinates; authoritative when supplied. |
| `fit_range` | `:axis` | `:axis` samples over padded visible x limits; `:data` samples from the first to last measured x. |
| `auto_limits` | `true` | Include data, errors, model, and displayed band in both axis limits. |
| `limit_padding` | `0.08` | Finite non-negative fractional padding around automatic content limits. |
| `plot_aspect` | `nothing` | Optional numeric `AxisAspect`; leave unset unless geometry carries meaning. |
| `axis_kwargs` | `NamedTuple()` | Makie `Axis` attributes applied after JuFitter's title/label defaults. |

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
| `band_color`, `band_alpha` | role-dependent | Direct JuFitter-level overrides. |
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
| `show_stats` | role-dependent | `true` for `:screen`, `false` for `:article`; an explicit Boolean always wins. |
| `stats_position` | `:right` | `:right` or `:inside`. |
| `inside_stats_position` | `:lt` | `:lt`, `:lb`, `:rt`, `:rb` and their long aliases. |
| `stats_panel_width` | `:auto` | Natural Makie width, a fraction `0 < w <= 1`, or a positive pixel width. Fractions are clamped to 300--560 px. |
| `panel_gap` | role-dependent | Gap between data axis and right panel. |
| `stats_mode` | `:compact` | `:compact` or `:full`. |
| `stats_sigdigits` | `5` | Significant digits used only for displayed values. |
| `parameter_names` | `nothing` | Display names; length must equal the number of fitted parameters. |
| `stats_fontsize` | role-dependent | Explicit result-panel or in-axis text size. |
| `stats_title` | `nothing` | Optional title above the structured right panel. |
| `latex_stats` | `false` | Render structured right-panel symbols and numbers as LaTeX. |
| `stats_box_color`, `stats_box_alpha` | role-dependent, `0.95` | In-axis summary background. |
| `stats_box_strokecolor`, `stats_box_strokewidth` | role-dependent, `1.0` | In-axis summary border. |
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

| Layer | JuFitter-level keywords | Final Makie container |
|---|---|---|
| Observations | `data_color`, `data_marker`, `data_markersize`, `data_strokecolor`, `data_strokewidth`, `data_label` | `scatter_kwargs` |
| Fit curve | `fit_color`, `fit_linewidth`, `fit_label` | `line_kwargs` |
| Band | `band_color`, `band_alpha`, `band_label` | `band_kwargs` |
| X errors | `xerr_color`, `error_whiskerwidth` | `xerrorbars_kwargs` |
| Y errors | `yerr_color`, `error_whiskerwidth` | `yerrorbars_kwargs` |

Every `*_kwargs` container accepts a `NamedTuple`, `AbstractDict`, or `nothing`.
Other container types raise `ArgumentError`. Explicit element overrides affect
only that layer; unmodified layers continue to follow the selected role.

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
theme = plot_theme(:screen; appearance=:dark)
style = plot_palette(:screen; appearance=:dark)
```

`plot_theme(style; appearance, theme_override)` returns the Makie `Theme` used
by JuFitter. `plot_palette(style; appearance)` returns the corresponding named
tuple of visual tokens: data/fit/band colors, multi-series colors, marker and
line sizes, typography, grids and spines, report-panel spacing, and default
figure sizes. These functions let a custom Makie layout inherit the same role
without copying private constants.

`plot_info_panel!` adds the same left-aligned information hierarchy used by
`plot_fit`:

```julia
plot_info_panel!(
    fig[1, 2];
    theme=:screen,
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
an axis. `fontsize`, `color`, `muted_color`, and `legend_kwargs` override role
defaults. `tellwidth=true` lets the panel report its natural width;
`tellheight=false` prevents a short report from shrinking or vertically
centering the adjacent scientific axis. The function returns its `GridLayout`.

## Residual, Pull, And Ratio Figures

```julia
plot_residuals(result; kind=:pull)
plot_diagnostics(result)
```

`plot_residuals` creates one axis. `kind` is:

| `kind` | Displayed value | Reference line |
|---|---|---:|
| `:residual` | ``y_i-f_i`` with available y errors | 0 |
| `:pull` | weighted or whitened residual coordinate | 0 |
| `:ratio` | ``y_i/f_i`` with propagated y-error ratio | 1 |

For dense covariance, whitened residual coordinates are not pointwise pulls in
the original measurement order. A ratio is undefined when any fitted model
value is zero and raises `ArgumentError` rather than plotting an infinity.
Non-finite coordinates, values, or errors are rejected.

Shared `plot_residuals` keywords are `filename=nothing`, `format=:pdf`,
`theme=:screen`, `appearance=:auto`, `theme_override=Theme()`,
`figure_size=(900, 520)`, `xlabel="x"`, `color=nothing`,
`reference_color=nothing`, `marker=nothing`, `markersize=nothing`,
`error_whiskerwidth=nothing`, `axis_kwargs`, `scatter_kwargs`, and
`errorbars_kwargs`.

`plot_diagnostics` stacks residual, pull, and ratio panels. It uses the same
keywords, defaults to `figure_size=(900, 900)`, and adds
`reference_line_kwargs`. All style containers are applied to every panel.

## Profile And Contour Figures

The numerical scans are defined in [Results And Diagnostics](api_results.md#Profiles-And-Contours).
Plotting a stored `ProfileResult`, `ContourResult`, or `ProfileMatrixResult`
does not perform another scan.

### One-parameter profile

```julia
plot_profile(profile_result; local_sigma=result.param_stderr[i])
```

| Keyword group | Keywords and defaults |
|---|---|
| Output and role | `filename=nothing`, `format=:pdf`, `theme=:screen`, `appearance=:auto`, `theme_override=Theme()`, `figure_size=(900, 620)` |
| Labels | `title="Profile"`, `xlabel="parameter"`, `ylabel="Delta cost"` |
| Profile | `line_color=nothing`, `line_width=nothing`, `profile_label="profile cost"`, `line_kwargs` |
| Local approximation | `local_sigma=nothing`, `local_color=nothing`, `local_linewidth=nothing`, `local_linestyle=:dash`, `local_label="local covariance parabola"`, `local_line_kwargs` |
| Threshold | `threshold_color=nothing`, `threshold_label=nothing`, `threshold_kwargs` |
| Layout | `show_legend=true`, `delta_max=nothing`, `axis_kwargs` |

`local_sigma` must be positive. `delta_max` must be positive when supplied; it
changes only the displayed y range, not the profile scan or interval.

### Two-parameter contour

```julia
plot_contour(
    contour_result;
    local_covariance=result.param_covariance,
    local_center=result.params[[i, j]],
)
```

| Keyword group | Keywords and defaults |
|---|---|
| Output and role | `filename=nothing`, `format=:pdf`, `theme=:screen`, `appearance=:auto`, `theme_override=Theme()`, `figure_size=(820, 700)` |
| Labels | `title="Contour"`, `xlabel="parameter 1"`, `ylabel="parameter 2"`, `axis_kwargs` |
| Profile surface | `show_regions=true`, `show_profile_lines=false`, `level_colors=nothing`, `region_colors=nothing`, `line_color=nothing`, `contour_kwargs` |
| Optional heatmap | `show_heatmap=false`, `colormap=:viridis`, `heatmap_kwargs` |
| Local approximation | `local_covariance=nothing`, `local_center=nothing`, `local_line_color=nothing`, `local_linewidth=nothing`, `local_linestyle=:dash`, `local_contour_kwargs` |
| Legend | `show_legend=true` |

The default uses filled profile regions and labels the common two-parameter
thresholds 2.30 and 6.18 as one- and two-sigma regions. A heatmap is opt-in.
`local_covariance` may be the relevant 2x2 matrix or the full fitted covariance;
`local_center` must contain exactly two values. Empty/non-positive contour
levels, an empty color collection, incompatible covariance dimensions, or a
surface with no finite cost value raise `ArgumentError`.

### Profile/contour matrix

```julia
plot_profile_matrix(result; parameters=[1, 2, 3])
plot_profile_matrix(matrix_result; parameter_names=["A", "lambda", "offset"])
```

The first method computes all required profile and contour refits, then renders
them. Its scan controls are `parameters=nothing`, `parameter_names=nothing`,
`npoints_profile=61`, `npoints_contour=25`, `nsigma=3`,
`profile_threshold=1.0`, `contour_levels=[2.30, 6.18]`, `adaptive=false`,
`max_refinements=2`, and `max_points=1200`.

Both methods accept `filename=nothing`, `format=:pdf`, `theme=:screen`,
`appearance=:auto`, `theme_override=Theme()`, `panel_status_mode=:issues`,
`delta_max=nothing`, and `figure_size=nothing`. The precomputed-result method
also accepts replacement `parameter_names` but no scan controls.

`panel_status_mode=:issues` labels only panels requiring attention; `:all`
labels every panel and `:none` suppresses status text. Invalid modes, mismatched
display-name counts, or inconsistent stored matrix dimensions raise
`ArgumentError`. Compute [`profile_matrix`](@ref) separately when scans should
run headlessly, be triaged before rendering, or be reused in several exports.

## Export Semantics

Every high-level plot function returns its `Figure` even when `filename` is
provided. A filename extension takes precedence over `format`; with no
extension, `.$(format)` is appended. For explicit resolution control, save the
returned figure with Makie:

```julia
fig = plot_fit(result; theme=:article)
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
| Invalid role, appearance, report mode, band, stats position, or fit range | `ArgumentError` |
| Non-positive/non-finite `nsigma`, negative/non-finite `limit_padding` | `ArgumentError` |
| Prediction band without matrix-free marginal errors | `ArgumentError` with the required remedy |
| Wrong number of `parameter_names` | `ArgumentError` |
| Non-finite or dimensionally inconsistent annotation data | `ArgumentError` |
| Ratio plot with a zero/non-finite model prediction | `ArgumentError` |
| Invalid profile/contour display geometry | `ArgumentError` rather than a misleading plot |

## API Documentation

```@docs
JuFitter.fitplot
JuFitter.plot_fit
JuFitter.fit_axis
JuFitter.add_curve!
JuFitter.add_points!
JuFitter.add_vline!
JuFitter.add_hline!
JuFitter.add_vband!
JuFitter.add_hband!
JuFitter.plot_theme
JuFitter.plot_palette
JuFitter.plot_info_panel!
JuFitter.plot_residuals
JuFitter.plot_diagnostics
JuFitter.plot_profile
JuFitter.plot_contour
JuFitter.plot_profile_matrix
```
