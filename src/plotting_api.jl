"""
    plot_fit(result; kwargs...)

Create and return a Makie `Figure` from an existing `FitResult`. The default
layout shows data, fitted model, uncertainty band, and an optional right-side
report without requiring manual margin tuning. Use `fitplot(model, x, y; ...)`
when fitting and plotting should happen in one call.

The default `fit_range=:axis` draws the fitted model over the padded axis range.
Use `fit_range=:data` or pass `xgrid` when the curve should stop at a specific
domain boundary.

Plotting is provided by the optional CairoMakie extension. Load it with
`using CairoMakie` before calling plot functions. Fitting and text reports work
without CairoMakie.
"""
function plot_fit end

"""
    fitplot(model, x, y; p0, kwargs...)
    fitplot(x, y; kwargs...)
    fitplot(result; kwargs...)

Convenience entry point for the common notebook workflow: fit data and return a
plot in one call. Every method returns the named tuple
`(result::FitResult, figure::Figure)`. Obtain the primary axis with
`fit_axis(output.figure)` when adding custom Makie content.

Use `report=:plot`, `:console`, `:both`, or `:none` to choose where fit
statistics are shown.
"""
function fitplot end

"""
    fit_axis(figure; index=1)

Return axis number `index` from a Makie `Figure`. The default returns the first
axis, which is the primary data axis of a standard JuFitter fit figure. Invalid
indices raise `ArgumentError`. This is the stable hook for adding custom Makie
content after `plot_fit` without rebuilding the fit.
"""
function fit_axis end

"""
    add_curve!(axis, f; xgrid=nothing, xspan=nothing, n=400, label=nothing, kwargs...)
    add_curve!(axis, x, y; label=nothing, kwargs...)

Add a function-valued or precomputed curve to an existing fit axis. A function
is sampled on `xgrid`, `xspan`, or the current visible axis range. Style
defaults follow the active JuFitter plot contract unless Makie keyword
arguments such as `color` or `linewidth` are explicitly supplied.
"""
function add_curve! end

"""
    add_points!(axis, x, y; label=nothing, kwargs...)

Add extra points or derived markers to an existing fit axis. This is intended
for thresholds, calibration anchors, extrapolated intersections, and other
scientific annotations that should share the plot layout.
"""
function add_points! end

"""
    add_vline!(axis, x; label=nothing, kwargs...)

Add vertical reference line(s) to an existing fit axis.
"""
function add_vline! end

"""
    add_hline!(axis, y; label=nothing, kwargs...)

Add horizontal reference line(s) to an existing fit axis.
"""
function add_hline! end

"""
    add_vband!(axis, xmin, xmax; label=nothing, kwargs...)

Add a vertical uncertainty band, acceptance region, excluded region, or physical
threshold to an existing fit axis.
"""
function add_vband! end

"""
    add_hband!(axis, ymin, ymax; label=nothing, kwargs...)

Add a horizontal uncertainty band, acceptance region, excluded region, or
physical threshold to an existing fit axis.
"""
function add_hband! end

"""
    plot_theme(style=:analysis; appearance=:auto)

Return the Makie `Theme` used by JuFitter for a named plot style. The maintained
output roles are `:analysis`, `:presentation`, and `:article`. Analysis keeps
the numerical result panel for live work, presentation prioritizes a large
screen-readable figure, and article uses TeX typography and vector-export
conventions. The former style names remain compatibility aliases. `appearance`
is `:light`, `:dark`, or `:auto`.
"""
function plot_theme end

"""
    plot_palette(style=:analysis; appearance=:auto)

Return the visual tokens used by JuFitter's plot helpers, including color-safe
series colors, markers, line weights, typography, and layout defaults. Use this
when building compound Makie figures that should respond to a selected style.
"""
function plot_palette end

"""
    plot_info_panel!(cell; theme=:analysis, appearance=:auto,
                     legend_source=nothing, legend_plots=nothing,
                     legend_labels=nothing, model_label=nothing,
                     parameter_lines=Any[], statistic_lines=Any[], kwargs...)

Add a compact, left-aligned information panel to a Makie layout `cell` and
return its `GridLayout`. `theme` and `appearance` use the same readable panel
defaults as `plot_fit`; explicit keywords override them. Supply either
`legend_source` or matching
`legend_plots`/`legend_labels`, plus already formatted model, parameter, and
statistic lines. This low-level helper lets compound figures follow JuFitter's
information hierarchy without coupling the panel to one `FitResult`.
"""
function plot_info_panel! end

"""
    plot_residuals(result; kind=:pull, theme=:analysis, kwargs...)

Plot residuals, pulls, or data/fit ratios for a fitted model. Use this when the
main fit plot looks plausible but the noise model or model structure needs
inspection. Marker and error-bar defaults follow the selected plot style;
explicit Makie keyword containers override them.
"""
function plot_residuals end

"""
    plot_diagnostics(result; kwargs...)

Create a residual, pull, and ratio diagnostic figure. All panels inherit the
selected JuFitter style; `scatter_kwargs`, `errorbars_kwargs`, and
`reference_line_kwargs` provide explicit Makie overrides.
"""
function plot_diagnostics end

"""
    plot_profile(profile_result; kwargs...)

Plot a one-parameter profile scan, including the fitted minimum, local
parabolic approximation when available, and configured threshold levels.
Profile, approximation, and threshold line styles can be overridden
independently without rebuilding the scan.
"""
function plot_profile end

"""
    plot_contour(contour_result; kwargs...)

Plot a two-parameter profile contour with labeled confidence regions and optional
local covariance overlay. Use this to detect non-elliptic likelihood geometry or
strong parameter correlations.
"""
function plot_contour end

"""
    plot_profile_matrix(result; parameters=nothing, kwargs...)
    plot_profile_matrix(matrix_result::ProfileMatrixResult; kwargs...)

Create a kafe2/Minuit-style overview matrix: profile scans on the diagonal and
pairwise contours below the diagonal. This is the quick diagnostic view for
correlation, non-parabolicity, active bounds, and failed refits.

Pass a precomputed `ProfileMatrixResult` to render an existing Makie-free
diagnostic object without repeating its profile and contour refits. The render
method accepts `parameter_names` to replace display labels without recomputing
the matrix.
"""
function plot_profile_matrix end

function _plotting_unavailable(name::Symbol)
    throw(ArgumentError(
        "`$name` requires the optional CairoMakie plotting extension. " *
        "Run `using CairoMakie` before calling JuFitter plotting functions, " *
        "or use the fitting/reporting APIs without plotting.",
    ))
end

fitplot(args...; kwargs...) = _plotting_unavailable(:fitplot)
plot_fit(args...; kwargs...) = _plotting_unavailable(:plot_fit)
fit_axis(args...; kwargs...) = _plotting_unavailable(:fit_axis)
add_curve!(args...; kwargs...) = _plotting_unavailable(:add_curve!)
add_points!(args...; kwargs...) = _plotting_unavailable(:add_points!)
add_vline!(args...; kwargs...) = _plotting_unavailable(:add_vline!)
add_hline!(args...; kwargs...) = _plotting_unavailable(:add_hline!)
add_vband!(args...; kwargs...) = _plotting_unavailable(:add_vband!)
add_hband!(args...; kwargs...) = _plotting_unavailable(:add_hband!)
plot_theme(args...; kwargs...) = _plotting_unavailable(:plot_theme)
plot_palette(args...; kwargs...) = _plotting_unavailable(:plot_palette)
plot_info_panel!(args...; kwargs...) = _plotting_unavailable(:plot_info_panel!)
plot_residuals(args...; kwargs...) = _plotting_unavailable(:plot_residuals)
plot_diagnostics(args...; kwargs...) = _plotting_unavailable(:plot_diagnostics)
plot_profile(args...; kwargs...) = _plotting_unavailable(:plot_profile)
plot_contour(args...; kwargs...) = _plotting_unavailable(:plot_contour)
plot_profile_matrix(args...; kwargs...) = _plotting_unavailable(:plot_profile_matrix)
