"""
    plot_fit(result; kwargs...)
    plot_fit(model, x, y; p0, kwargs...)

Create a Makie fit figure from an existing `FitResult` or directly from data and
a model. The default layout shows data, fitted model, uncertainty band, and an
optional right-side report without requiring manual margin tuning.

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
plot object in one call. The returned object contains the `FitResult`, Makie
figure, axis, and generated plot handles when CairoMakie is loaded.

Use `report=:plot`, `:console`, `:both`, or `:none` to choose where fit
statistics are shown.
"""
function fitplot end

"""
    fit_axis(plot_or_axis)

Return the primary Makie `Axis` from a JuFitter plot object or pass an existing
axis through unchanged. This is the stable hook for adding custom Makie content
after `plot_fit` without rebuilding the fit.
"""
function fit_axis end

"""
    add_curve!(axis, model, x; params, label=nothing, kwargs...)
    add_curve!(axis, x, y; label=nothing, kwargs...)

Add a model curve or precomputed curve to an existing fit axis. Style defaults
follow the active JuFitter plot contract unless Makie keyword arguments such as
`color` or `linewidth` are explicitly supplied.
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
    plot_theme(style=:lab; appearance=:auto)

Return the Makie `Theme` used by JuFitter for a named plot style. Supported
styles are `:lab`, `:modern`, and `:article`; legacy aliases such as
`:workbench`, `:showcase`, and `:publication` remain accepted. `appearance` is
`:light`, `:dark`, or `:auto`.
"""
function plot_theme end

"""
    plot_palette(style=:lab; appearance=:auto)

Return the color and marker defaults used by JuFitter's plot helpers. Use this
when building compound Makie figures that should respond to the documentation
or user-selected plot style.
"""
function plot_palette end

"""
    plot_info_panel!(figure_or_grid, result; kwargs...)

Add a compact right-side information panel with model, parameters, and fit
statistics. The panel is the reusable report component used by `plot_fit`.
"""
function plot_info_panel! end

"""
    plot_residuals(result; kwargs...)

Plot residuals or pulls for a fitted model. Use this when the main fit plot
looks plausible but the noise model or model structure needs inspection.
"""
function plot_residuals end

"""
    plot_diagnostics(result; kwargs...)

Create a diagnostic figure from structured JuFitter findings. The visual
diagnostic layer is meant to point to the next inspection step, not only to draw
residuals.
"""
function plot_diagnostics end

"""
    plot_profile(profile_result; kwargs...)

Plot a one-parameter profile scan, including the fitted minimum, local
parabolic approximation when available, and configured threshold levels.
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
    plot_profile_matrix(result, parameters; kwargs...)

Create a kafe2/Minuit-style overview matrix: profile scans on the diagonal and
pairwise contours below the diagonal. This is the quick diagnostic view for
correlation, non-parabolicity, active bounds, and failed refits.
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
