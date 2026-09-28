# Fit Plotting

Fitting, profiles, diagnostics, and text reports do not require Makie; the
plotting methods are activated by CairoMakie:

```julia
using ScientificFitting
using CairoMakie
```

Without `using CairoMakie`, every plotting entry point raises an
`ArgumentError` that names the missing optional extension.

## Choose The Plotting Entry Point

| Task | Function | Return value | Runs an optimizer? |
|---|---|---|---:|
| Fit arrays and plot immediately | [`fitplot`](@ref) | `(result, figure)` | yes |
| Plot an existing x-y fit | [`plot_fit`](@ref) | `Figure` | no |
| Add content to the data axis | [`fit_axis`](@ref), `add_*!` | `Axis` or Makie plot object | no |
| Compose a custom themed figure | [`plot_theme`](@ref), [`plot_palette`](@ref), [`plot_info_panel!`](@ref), [`resize_plot_to_layout!`](@ref) | `Theme`, token `NamedTuple`, `GridLayout`, or the resized `Figure` | no |

Complete composition examples live in
[Plotting And Customization](plotting_design.md); this page is the argument
and failure contract.

## Fit And Plot In One Call

```text
fitplot(model, x, y; p0, show_panel=true, print_report=false, kwargs...)
fitplot(x, y; p0=nothing, show_panel=true, print_report=false, kwargs...)
fitplot(result::FitResult; show_panel=true, print_report=false, kwargs...)
```

All methods return `(result=result, figure=figure)`. The two-array method fits
the straight line ``y=p_1x+p_2``, deriving initial slope and intercept from the
first and last observations unless `p0` is given. The `FitResult` method only
renders and never repeats the fit.

### Output selection

`show_panel::Bool` controls the right/inside numerical panel;
`print_report::Bool` prints `report_text(result)` to the terminal and exists
only on `fitplot`. The two are independent.

### Keyword routing

For methods that perform a fit, the following keywords are sent to
[`fit_model`](@ref):

| Concern | Fitting keywords |
|---|---|
| Observation uncertainty | `sigma_y`, `sigma_x`, `cov_y`, `cov_x`, `whitening`, `error_components` |
| Parameter information | `bounds`, `constraints`, `parameter_priors`, `parameter_constraints`, `fixed_parameters` |
| Derivatives and model evaluation | `jacobian`, `x_derivative`, `inplace`, `derivatives` |
| Solver and covariance behavior | `solver`, `cost`, `maxiters`, `tol`, `scale_covariance`, `initial_guesses`, `multistart` |

All remaining keywords are sent to [`plot_fit`](@ref); a misspelled fitting
keyword does not disappear silently but fails there as an unsupported keyword.

## Plot An Existing Result

```text
plot_fit(result::FitResult; kwargs...) -> Figure
```

`plot_fit` draws the observations, available x/y error bars, fitted model,
optional uncertainty band, and optional result panel, without modifying
`result` or rerunning the optimizer.

### Output, style, and appearance

| Keyword | Default | Contract |
|---|---:|---|
| `filename` | `nothing` | Save during construction when a path is supplied. |
| `format` | `:pdf` | Extension appended only when `filename` has no extension. |
| `theme` | `:sans` | Maintained visual style: `:sans` or `:tex`. |
| `appearance` | `:auto` | `:light`, `:dark`, or `:auto`; `:auto` currently resolves to light. |
| `theme_override` | `Theme()` | Makie theme merged after the selected ScientificFitting style. |
| `style` | `FitPlotStyle()` | Per-figure visual token overrides — colors, markers, line widths, `figure_size`, panel and stats-box styling. See [`FitPlotStyle`](@ref). |
| `tight_layout` | `true` | Remove empty layout rows and columns before the automatic fit-to-content pass. |

The former names `:analysis`, `:presentation`, `:screen`, `:lab`, `:workbench`,
`:modern`, `:clean`, `:minimal`, and `:showcase` map to `:sans`; `:article`,
`:publication`, `:paper`, and `:latex` map to `:tex`; new code should use the
maintained names. Unknown styles and appearances raise `ArgumentError`.

### Labels, units, model domain, and limits

| Keyword | Default | Contract |
|---|---:|---|
| `title` | `nothing` | Figure title; `nothing` produces no title. |
| `model_label` | automatic for the built-in line | Model expression shown in the right panel. |
| `xlabel`, `ylabel` | `"x"`, `"y"` | Axis quantity labels. |
| `xunit`, `yunit` | `nothing` | Appended as `label / unit` (quantity-calculus notation: a tick value 2 on an axis labeled `t / s` means t = 2 s); units are never inferred. |
| `latex_labels` | `false` | Convert suitable labels to Makie `LaTeXString` content. Pass explicit `L"..."` strings for mathematical notation. |
| `xgrid` | `nothing` | Explicit finite model-sampling coordinates; authoritative when supplied. |
| `fit_range` | `:axis` | `:axis` samples the data x range extended by `limit_padding`; `:data` samples from the smallest to the largest measured x. |
| `auto_limits` | `true` | Include data, errors, model, and displayed band in both axis limits. |
| `limit_padding` | `0.08` | Finite non-negative fractional padding around automatic content limits. |
| `plot_aspect` | `nothing` | Optional numeric `AxisAspect`; leave unset unless geometry carries meaning. |
| `axis_kwargs` | `NamedTuple()` | Makie `Axis` attributes applied after ScientificFitting's title/label defaults. |

With `auto_limits=false`, set limits through `axis_kwargs` or on
`fit_axis(figure)` after construction. Manual limits do not resample the
model; pass a matching `xgrid` for intentional extrapolation.

### Uncertainty band

| Keyword | Default | Contract |
|---|---:|---|
| `band` | `:confidence` | `:confidence`, `:prediction`, or `:none`. |
| `nsigma` | `1.0` | Finite positive multiplier for the displayed standard-deviation scale. |
| `band_label` | `nothing` | Legend text. `nothing` derives `"<nsigma>-sigma band"` from `nsigma`, e.g. `"2-sigma band"` for `nsigma=2`. The derived label does not name the band type; set the label explicitly when the meaning should be named, e.g. `band_label="2-sigma prediction band"` with `band=:prediction`. |
| `band_kwargs` | `NamedTuple()` | Makie `band!` attributes applied last. |

Band color and opacity are the `band_color` and `band_alpha` fields of
[`FitPlotStyle`](@ref).

`band=:confidence` propagates the local parameter covariance to the fitted
mean; `band=:prediction` adds pointwise observation uncertainty in y and the
effective contribution from x uncertainty. Both are
[local covariance](statistics.md#Local-Parameter-Covariance) constructions.

A matrix-free [`WhiteningOperator`](@ref) can define the fit without exposing
pointwise marginal errors; `band=:prediction` then requires `marginal_sigma`
and raises `ArgumentError` otherwise. `band=:confidence` remains available.

### Result panel and legend

| Keyword | Default | Contract |
|---|---:|---|
| `show_panel` | `true` | Show the structured right panel or compact in-axis panel. Independent of visual style. |
| `stats_position` | `:right` | `:right` or `:inside`. |
| `inside_stats_position` | `:lt` | `:lt`/`:lefttop`, `:lb`/`:leftbottom`, `:rt`/`:righttop`, `:rb`/`:rightbottom`. |
| `stats_panel_width` | `:auto` | Natural Makie width, a fraction `0 < w <= 1`, or a positive wrapping width. Fractions are clamped to 300--560 px; unbreakable TeX or legend content may expand the panel. |
| `stats_mode` | `:compact` | `:compact` or `:full`. |
| `stats_sigdigits` | `5` | Significant digits used only for displayed values. |
| `parameter_names` | `nothing` | Display names; length must equal the total number of model parameters, including fixed ones (fixed parameters are marked `(fixed)` in the panel). |
| `stats_title` | `nothing` | Optional title above the structured right panel. |
| `latex_stats` | `false` | Render structured right-panel symbols and numbers as LaTeX; requires `stats_position=:right` (the inside box renders plain text and rejects the combination with `ArgumentError`). |
| `show_legend` | `true` | Show data, fit, and band labels. With a right panel, the legend is placed above the report. |
| `legend_position` | `:rt` | In-axis Makie legend position when no right-side panel owns the legend. |
| `legend_kwargs` | `NamedTuple()` | Makie legend attributes applied last. |

Panel gap, panel text size, and the in-axis box background and border are the
`panel_gap`, `stats_fontsize`, and `stats_box_*` fields of
[`FitPlotStyle`](@ref).

In the right panel, `stats_mode=:full` adds cost, AIC, and BIC to the compact
parameter, chi-square, chi-square/ndf, p-value, and ndf rows; the in-axis box
gains the raw chi-square. AIC and BIC are displayed values; their comparison rules are in
[Results And Diagnostics](api_results.md#Results-And-Diagnostics).

### Data, fit, and error-bar styling

Per-layer visual tokens are [`FitPlotStyle`](@ref) fields, passed as one
`style` keyword; a field left at `nothing` keeps the selected style's role
default. Label keywords stay top-level, and each layer's Makie keyword
container is merged last and has final authority.

| Layer | `FitPlotStyle` fields | Label keyword | Final Makie container |
|---|---|---|---|
| Observations | `data_color`, `data_marker`, `data_markersize`, `data_strokecolor`, `data_strokewidth` | `data_label` | `scatter_kwargs` |
| Fit curve | `fit_color`, `fit_linewidth` | `fit_label` | `line_kwargs` |
| Band | `band_color`, `band_alpha` | `band_label` | `band_kwargs` |
| X errors | `xerr_color`, `error_linewidth`, `error_whiskerwidth` | — | `xerrorbars_kwargs` |
| Y errors | `yerr_color`, `error_linewidth`, `error_whiskerwidth` | — | `yerrorbars_kwargs` |

Every `*_kwargs` container accepts a `NamedTuple`, `AbstractDict`, or
`nothing`; other container types raise `ArgumentError`. Explicit overrides
affect only their layer.

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
limits; none of these helpers changes or reruns the fit, and their `kwargs...`
are ordinary Makie plot attributes.

## Reuse The Visual Contract

```julia
theme = plot_theme(:sans; appearance=:dark)
style = plot_palette(:sans; appearance=:dark)
```

`plot_theme(theme; appearance, theme_override)` returns the Makie `Theme` used
by ScientificFitting; `plot_palette(style; appearance)` returns the
corresponding named tuple of visual tokens, from data/fit/band and
multi-series colors and marker and line sizes to typography, grids and spines,
panel spacing, and default figure sizes.

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

`legend_source=axis` builds the legend from labeled axis content instead.
`fontsize`, `color`, `muted_color`, and `legend_kwargs` override style
defaults. `width=nothing` reports the panel's natural Makie width; an explicit
`width` wraps long plain-text lines. `tellwidth=true` reports that width to
the parent layout, `tellheight=true` lets a long report enlarge its row. The
function returns its `GridLayout`.

After adding every custom layout block, fit the canvas once:

```julia
resize_plot_to_layout!(
    fig;
    minimum_axis_size=(420, nothing), # keep width; retain compound row ratios
)
```

Existing explicit axis sizes remain authoritative, the current figure size is
a lower bound, and additional width goes to the first graph column. When graph
axes occupy other top-level columns, pass `flexible_columns=(...)`: each
listed `Auto` column is measured and pinned with Makie's `Auto(false, ratio)`
mode so labels and legends cannot shrink it. Explicit `Fixed` and `Relative`
tracks remain authoritative.

## Diagnostic Figures

Residual, pull, ratio, profile, contour, and profile-matrix figures are listed
in [Diagnostic Plotting](api_plotting_diagnostics.md).

## Export Semantics

Every high-level plot function returns its `Figure` even when `filename` is
provided; a filename extension takes precedence over `format`, and without one
`.$(format)` is appended. For explicit resolution control, save the returned
figure with Makie:

```julia
fig = plot_fit(result; theme=:tex)
save("fit.svg", fig)
save("fit.png", fig; px_per_unit=2)
```

`FitPlotStyle(figure_size=...)` requests the minimum logical layout size:
fixed-width content grows the canvas instead of clipping, and additional
requested width widens the flexible data axis. Raster density is controlled by
`px_per_unit`; enlarging `figure_size` and scaling the image down also scales
down its text.

## Failure Summary

| Failure | Result |
|---|---|
| CairoMakie extension not loaded | `ArgumentError` naming the required extension |
| Invalid style, appearance, band, stats position, stats mode, or fit range | `ArgumentError` |
| `latex_stats=true` with `stats_position=:inside` | `ArgumentError` |
| Non-positive/non-finite `nsigma` | `DomainError` |
| Unordered band bounds or `n < 2` curve samples | `DomainError` |
| Negative/non-finite `limit_padding` | `ArgumentError` |
| Prediction band without matrix-free marginal errors | `ArgumentError` with the required remedy |
| Wrong number of `parameter_names` | `DimensionMismatch` |
| Non-finite annotation data | `ArgumentError` |
| Dimensionally inconsistent annotation lengths | `DimensionMismatch` |
| Non-positive/non-finite layout dimensions | `DomainError` |

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
ScientificFitting.FitPlotStyle
ScientificFitting.plot_palette
ScientificFitting.plot_info_panel!
ScientificFitting.resize_plot_to_layout!
```
