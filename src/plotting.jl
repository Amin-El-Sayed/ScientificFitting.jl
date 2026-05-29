function _publication_theme()
    return Theme(
        fontsize=22,
        figure_padding=16,
        Axis=(
            xlabelsize=26,
            ylabelsize=26,
            xticklabelsize=18,
            yticklabelsize=18,
            xgridvisible=true,
            ygridvisible=true,
            xminorgridvisible=true,
            yminorgridvisible=true,
        ),
        Lines=(linewidth=3,),
        Scatter=(markersize=10,),
    )
end

function _clean_theme()
    return Theme(
        fontsize=20,
        figure_padding=(18, 22, 12, 16),
        Axis=(
            xlabelsize=24,
            ylabelsize=24,
            titlesize=25,
            xticklabelsize=17,
            yticklabelsize=17,
            xgridvisible=true,
            ygridvisible=true,
            xminorgridvisible=false,
            yminorgridvisible=false,
            xgridcolor=(:gray70, 0.35),
            ygridcolor=(:gray70, 0.35),
            topspinevisible=false,
            rightspinevisible=false,
        ),
        Lines=(linewidth=3,),
        Scatter=(markersize=10,),
    )
end

function _latex_theme()
    return Theme(
        fontsize=22,
        figure_padding=16,
        Axis=(
            xlabelsize=26,
            ylabelsize=26,
            titlesize=28,
            xticklabelsize=18,
            yticklabelsize=18,
            xgridvisible=true,
            ygridvisible=true,
            xminorgridvisible=true,
            yminorgridvisible=true,
        ),
        font="CMU Serif",
        Lines=(linewidth=3,),
        Scatter=(markersize=10,),
    )
end

function _normalize_kwargs(kwargs)
    kwargs === nothing && return NamedTuple()
    kwargs isa NamedTuple && return kwargs
    kwargs isa AbstractDict && return (; (Symbol(k) => v for (k, v) in kwargs)...)
    throw(ArgumentError("plot keyword containers must be NamedTuple, Dict, or nothing"))
end

function _merged_kwargs(defaults::NamedTuple, overrides)
    return merge(defaults, _normalize_kwargs(overrides))
end

function _prediction_band_sigma(result::FitResult, xgrid::AbstractVector)
    J = _parameter_jacobian(result.problem, result.params; x=xgrid)
    cov = result.param_covariance
    tmp = J * cov
    variances = vec(sum(tmp .* J; dims=2))
    return sqrt.(clamp.(variances, 0.0, Inf))
end

function _panel_width_px(stats_panel_width::Real, fig_width::Int)
    stats_panel_width > 0 || throw(ArgumentError("stats_panel_width must be positive"))
    if stats_panel_width <= 1
        return max(280, Int(round(fig_width * stats_panel_width)))
    end
    return Int(round(stats_panel_width))
end

function _label_with_unit(label, unit)
    unit === nothing && return label
    unit_text = string(unit)
    isempty(unit_text) && return label
    label_text = string(label)
    return isempty(label_text) ? unit_text : string(label_text, " (", unit_text, ")")
end

function _finite_extrema(values)
    finite_values = Float64[v for v in values if isfinite(v)]
    isempty(finite_values) && return nothing
    return minimum(finite_values), maximum(finite_values)
end

function _padded_limits(values; padding::Real=0.08)
    extrema = _finite_extrema(values)
    extrema === nothing && return nothing
    lo, hi = extrema
    if lo == hi
        delta = max(abs(lo), 1.0)
        return lo - 0.5 * delta, hi + 0.5 * delta
    end
    pad = Float64(padding) * (hi - lo)
    return lo - pad, hi + pad
end

function _fit_plot_limits(x, y, xerr, yerr, xgrid, ygrid, band_sigma; padding::Real)
    x_values = Float64[]
    append!(x_values, Float64.(x))
    append!(x_values, Float64.(xgrid))
    if xerr !== nothing
        append!(x_values, Float64.(x .- xerr))
        append!(x_values, Float64.(x .+ xerr))
    end

    y_values = Float64[]
    append!(y_values, Float64.(y))
    append!(y_values, Float64.(ygrid))
    if yerr !== nothing
        append!(y_values, Float64.(y .- yerr))
        append!(y_values, Float64.(y .+ yerr))
    end
    if band_sigma !== nothing
        append!(y_values, Float64.(ygrid .- band_sigma))
        append!(y_values, Float64.(ygrid .+ band_sigma))
    end

    return _padded_limits(x_values; padding=padding), _padded_limits(y_values; padding=padding)
end

function _default_grid(x::AbstractVector)
    xlo = minimum(x)
    xhi = maximum(x)
    if xlo == xhi
        xlo -= 0.5
        xhi += 0.5
    end
    return collect(range(xlo, xhi; length=800))
end

function _fmt_value(x::Real; sigdigits::Int=5)
    if isnan(x)
        return "NaN"
    elseif isinf(x)
        return signbit(x) ? "-Inf" : "Inf"
    end
    return string(round(Float64(x); sigdigits=sigdigits))
end

function _strip_math_delims(s::AbstractString)
    if startswith(s, "\$") && endswith(s, "\$") && ncodeunits(s) >= 2
        return s[2:(end - 1)]
    end
    return s
end

function _to_latex_text(s::AbstractString)
    escaped = replace(s, "_" => "\\_")
    return LaTeXString("\\text{" * escaped * "}")
end

function _latex_symbol_expr(name)
    raw = name isa LaTeXString ? String(name) : string(name)
    s = _strip_math_delims(raw)
    if occursin("\\", s) || occursin("^", s) || occursin("_", s)
        return s
    end
    return "\\mathrm{" * replace(s, "_" => "\\_") * "}"
end

function _stats_panel_lines(
    result::FitResult;
    parameter_names::Union{Nothing, AbstractVector}=nothing,
    sigdigits::Int=5,
    latex_stats::Bool=true,
)
    n = length(result.params)
    names = if parameter_names === nothing
        Any["p$i" for i in 1:n]
    else
        length(parameter_names) == n || throw(ArgumentError("parameter_names length must match parameter count"))
        collect(parameter_names)
    end

    lines = Any[latex_stats ? LaTeXString("\\textbf{Fit\\ Summary}") : "Fit Summary"]
    for i in 1:n
        v = _fmt_value(result.params[i]; sigdigits=sigdigits)
        e = _fmt_value(result.param_stderr[i]; sigdigits=sigdigits)
        if latex_stats
            pexpr = _latex_symbol_expr(names[i])
            push!(lines, LaTeXString(pexpr * " = " * v * " \\pm " * e))
        else
            push!(lines, string(names[i], " = ", v, " ± ", e))
        end
    end

    chi2_text = _fmt_value(result.stats.chi2; sigdigits=sigdigits)
    chi2_ndf_text = _fmt_value(result.stats.chi2_ndf; sigdigits=sigdigits)
    ndf_text = string(result.stats.ndf)

    if latex_stats
        push!(lines, L"\chi^2 = %$chi2_text")
        push!(lines, L"\mathrm{n_{dof}} = %$ndf_text")
        push!(lines, L"\chi^2/\mathrm{n_{dof}} = %$chi2_ndf_text")
    else
        push!(lines, string("chi2 = ", chi2_text))
        push!(lines, string("ndf = ", ndf_text))
        push!(lines, string("chi2/ndf = ", chi2_ndf_text))
    end

    return lines
end

function _as_label_text(value, latex_labels::Bool)
    if !latex_labels
        return value
    end
    if value isa LaTeXString
        return value
    end
    if value isa AbstractString
        s = _strip_math_delims(value)
        if occursin("\\", s) || occursin("^", s) || occursin("_", s)
            return LaTeXString(s)
        end
        return _to_latex_text(s)
    end
    return value
end

const _FITPLOT_FIT_KWARGS = Set([
    :sigma_y,
    :sigma_x,
    :cov_y,
    :cov_x,
    :error_components,
    :bounds,
    :constraints,
    :parameter_priors,
    :parameter_constraints,
    :fixed_parameters,
    :jacobian,
    :backend,
    :cost,
    :maxiters,
    :tol,
    :ci_level,
    :scale_covariance,
    :initial_guesses,
    :multistart,
])

function _split_fitplot_kwargs(kwargs)
    fit_kwargs = Dict{Symbol, Any}()
    plot_kwargs = Dict{Symbol, Any}()
    for (key, value) in pairs(kwargs)
        if key in _FITPLOT_FIT_KWARGS
            fit_kwargs[key] = value
        else
            plot_kwargs[key] = value
        end
    end
    return (; fit_kwargs...), (; plot_kwargs...)
end

function _linear_initial_guess(x::AbstractVector, y::AbstractVector)
    x1 = Float64(first(x))
    x2 = Float64(last(x))
    y1 = Float64(first(y))
    y2 = Float64(last(y))
    slope = x1 == x2 ? 0.0 : (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1
    return [slope, intercept]
end

_default_linear_model(x, p) = @. p[1] * x + p[2]

function _normalize_fitplot_report(report::Symbol)
    report in (:plot, :console, :both, :none) || throw(ArgumentError("report must be :plot, :console, :both, or :none"))
    return report
end

function _fitplot_result(result::FitResult; report::Symbol=:plot, kwargs...)
    report = _normalize_fitplot_report(report)
    plot_kwargs = Dict{Symbol, Any}(pairs(kwargs))
    parameter_names = get(plot_kwargs, :parameter_names, nothing)

    if !haskey(plot_kwargs, :show_stats)
        plot_kwargs[:show_stats] = report in (:plot, :both)
    end

    if report in (:console, :both)
        println(report_text(result; parameter_names=parameter_names))
    end

    fig = plot_fit(result; (; plot_kwargs...)...)
    return (result=result, figure=fig)
end

"""
    fitplot(result::FitResult; report=:plot, kwargs...)
    fitplot(model, x, y; p0, report=:plot, kwargs...)
    fitplot(x, y; p0=nothing, report=:plot, kwargs...)

Fit and plot in one call. The `model, x, y` method forwards fitting keywords
such as `sigma_y`, `sigma_x`, `bounds`, `parameter_priors`, and `backend` to
`fit_model`; plotting keywords such as `xlabel`, `ylabel`, `theme`, `nsigma`,
`report`, and `filename` are forwarded to `plot_fit`.

The `x, y` method uses a linear model by default. All methods return a named
tuple `(result, figure)` so the numerical result is not lost.
"""
function fitplot(result::FitResult; report::Symbol=:plot, kwargs...)
    return _fitplot_result(result; report=report, kwargs...)
end

function fitplot(model, x::AbstractVector, y::AbstractVector; p0::AbstractVector, report::Symbol=:plot, kwargs...)
    fit_kwargs, plot_kwargs = _split_fitplot_kwargs(kwargs)
    result = fit_model(model, x, y; p0=p0, fit_kwargs...)
    return _fitplot_result(result; report=report, plot_kwargs...)
end

function fitplot(x::AbstractVector, y::AbstractVector; p0=nothing, report::Symbol=:plot, kwargs...)
    initial = p0 === nothing ? _linear_initial_guess(x, y) : p0
    return fitplot(_default_linear_model, x, y; p0=initial, report=report, kwargs...)
end

"""
    plot_fit(
        result::FitResult;
        xgrid=nothing,
        filename=nothing,
        format=:pdf,
        theme=:publication,
        theme_override=Theme(),
        title="Fit Result",
        xlabel="x",
        ylabel="y",
        xunit=nothing,
        yunit=nothing,
        auto_limits=true,
        limit_padding=0.08,
        plot_aspect=nothing,
        figure_size=nothing,
        stats_panel_width=0.42,
        panel_gap=8,
        latex_labels=false,
        latex_stats=true,
        show_stats=true,
        stats_sigdigits=5,
        parameter_names=nothing,
        stats_fontsize=18,
        stats_box_color=:white,
        stats_box_alpha=0.95,
        stats_box_strokecolor=:black,
        stats_box_strokewidth=1.0,
        stats_linegap=2,
        stats_label_kwargs=NamedTuple(),
        stats_title_kwargs=NamedTuple(),
        stats_box_kwargs=NamedTuple(),
        show_legend=true,
        legend_position=:rt,
        axis_kwargs=NamedTuple(),
        legend_kwargs=NamedTuple(),
        data_color=:black,
        data_marker=:circle,
        data_markersize=10,
        scatter_kwargs=NamedTuple(),
        fit_color=:dodgerblue4,
        fit_linewidth=3,
        fit_label="fit",
        line_kwargs=NamedTuple(),
        band_color=:dodgerblue,
        band_alpha=0.20,
        band=:confidence,
        nsigma=1,
        band_label="1-sigma band",
        band_kwargs=NamedTuple(),
        xerr_color=:black,
        yerr_color=:black,
        error_whiskerwidth=10,
        xerrorbars_kwargs=NamedTuple(),
        yerrorbars_kwargs=NamedTuple(),
        data_label="data",
    )

Create a publication-style fit plot with data, error bars, best-fit curve, 1-sigma band,
and an optional right-side statistics panel. Makie keyword containers can be passed as
`NamedTuple`s or `Dict`s via the `*_kwargs` arguments.
"""
function plot_fit(
    result::FitResult;
    xgrid=nothing,
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:publication,
    theme_override::Theme=Theme(),
    title="Fit Result",
    xlabel="x",
    ylabel="y",
    xunit=nothing,
    yunit=nothing,
    auto_limits::Bool=true,
    limit_padding::Real=0.08,
    plot_aspect::Union{Nothing, Real}=nothing,
    figure_size::Union{Nothing, Tuple{<:Real, <:Real}}=nothing,
    stats_panel_width::Real=0.42,
    panel_gap::Real=8,
    latex_labels::Bool=false,
    latex_stats::Bool=true,
    show_stats::Bool=true,
    stats_sigdigits::Int=5,
    parameter_names::Union{Nothing, AbstractVector}=nothing,
    stats_fontsize::Real=18,
    stats_box_color=:white,
    stats_box_alpha::Real=0.95,
    stats_box_strokecolor=:black,
    stats_box_strokewidth::Real=1.0,
    stats_linegap::Real=2,
    stats_label_kwargs=NamedTuple(),
    stats_title_kwargs=NamedTuple(),
    stats_box_kwargs=NamedTuple(),
    show_legend::Bool=true,
    legend_position=:rt,
    axis_kwargs=NamedTuple(),
    legend_kwargs=NamedTuple(),
    data_color=:black,
    data_marker=:circle,
    data_markersize::Real=10,
    scatter_kwargs=NamedTuple(),
    fit_color=:dodgerblue4,
    fit_linewidth::Real=3,
    fit_label="fit",
    line_kwargs=NamedTuple(),
    band_color=:dodgerblue,
    band_alpha::Real=0.20,
    band::Symbol=:confidence,
    nsigma::Real=1.0,
    band_label="1-sigma band",
    band_kwargs=NamedTuple(),
    xerr_color=:black,
    yerr_color=:black,
    error_whiskerwidth::Real=10,
    xerrorbars_kwargs=NamedTuple(),
    yerrorbars_kwargs=NamedTuple(),
    data_label="data",
)
    band in (:confidence, :none) || throw(ArgumentError("band must be :confidence or :none"))

    thm = theme == :clean ? _clean_theme() : theme == :publication ? _publication_theme() : theme == :latex ? _latex_theme() : Theme()
    thm = merge(thm, theme_override)

    base_size = show_stats ? (1250, 720) : (980, 720)
    fig_size = figure_size === nothing ? base_size : (Int(round(figure_size[1])), Int(round(figure_size[2])))
    fig = with_theme(thm) do
        Figure(size=fig_size)
    end

    if show_stats
        colgap!(fig.layout, Int(round(panel_gap)))
    end

    axis_defaults = (
        title=_as_label_text(title, latex_labels),
        xlabel=_as_label_text(_label_with_unit(xlabel, xunit), latex_labels),
        ylabel=_as_label_text(_label_with_unit(ylabel, yunit), latex_labels),
    )
    if plot_aspect !== nothing
        axis_defaults = merge(axis_defaults, (aspect=AxisAspect(Float64(plot_aspect)),))
    end
    ax = Axis(fig[1, 1]; _merged_kwargs(axis_defaults, axis_kwargs)...)

    x = result.problem.x
    y = result.problem.y
    xerr = _xerror_for_plot(result.problem, result.params)
    yerr = _yerror_for_plot(result.problem, result.params)

    xg = xgrid === nothing ? _default_grid(x) : collect(Float64, xgrid)
    yg = _model_values(result.problem, result.params; x=xg)
    sg = Float64(nsigma) .* _prediction_band_sigma(result, xg)

    bplot = nothing
    if band == :confidence
        bplot = band!(
            ax,
            xg,
            yg .- sg,
            yg .+ sg;
            _merged_kwargs(
                (color=(band_color, band_alpha), label=_as_label_text(band_label, latex_labels)),
                band_kwargs,
            )...,
        )
    end
    fplot = lines!(
        ax,
        xg,
        yg;
        _merged_kwargs(
            (color=fit_color, linewidth=fit_linewidth, label=_as_label_text(fit_label, latex_labels)),
            line_kwargs,
        )...,
    )

    if yerr !== nothing
        errorbars!(
            ax,
            x,
            y,
            yerr;
            _merged_kwargs((color=yerr_color, whiskerwidth=error_whiskerwidth), yerrorbars_kwargs)...,
        )
    end

    if xerr !== nothing
        errorbars!(
            ax,
            x,
            y,
            xerr;
            _merged_kwargs(
                (direction=:x, color=xerr_color, whiskerwidth=error_whiskerwidth),
                xerrorbars_kwargs,
            )...,
        )
    end

    dplot = scatter!(
        ax,
        x,
        y;
        _merged_kwargs(
            (color=data_color, marker=data_marker, markersize=data_markersize, label=_as_label_text(data_label, latex_labels)),
            scatter_kwargs,
        )...,
    )

    if auto_limits
        xlims, ylims = _fit_plot_limits(x, y, xerr, yerr, xg, yg, band == :confidence ? sg : nothing; padding=limit_padding)
        xlims !== nothing && ylims !== nothing && limits!(ax, xlims..., ylims...)
    end

    if show_legend
        legend_plots = band == :confidence ? [bplot, fplot, dplot] : [fplot, dplot]
        legend_labels = band == :confidence ?
            [
                _as_label_text(band_label, latex_labels),
                _as_label_text(fit_label, latex_labels),
                _as_label_text(data_label, latex_labels),
            ] :
            [
                _as_label_text(fit_label, latex_labels),
                _as_label_text(data_label, latex_labels),
            ]
        axislegend(
            ax,
            legend_plots,
            legend_labels,
            ;
            _merged_kwargs((position=legend_position,), legend_kwargs)...,
        )
    end

    if show_stats
        panel_width_px = _panel_width_px(stats_panel_width, fig_size[1])
        stats_lines = _stats_panel_lines(
            result;
            parameter_names=parameter_names,
            sigdigits=stats_sigdigits,
            latex_stats=latex_stats,
        )
        Box(
            fig[1, 2];
            _merged_kwargs(
                (
                    color=(stats_box_color, stats_box_alpha),
                    strokecolor=stats_box_strokecolor,
                    strokewidth=stats_box_strokewidth,
                ),
                stats_box_kwargs,
            )...,
        )
        stats_grid = GridLayout(fig[1, 2], alignmode=Inside())
        rowgap!(stats_grid, Int(round(stats_linegap)))
        label_defaults = (
            halign=:left,
            valign=:top,
            justification=:left,
            lineheight=1.0,
            word_wrap=true,
            tellwidth=false,
            tellheight=true,
            padding=(12, 12, 1, 1),
            fontsize=stats_fontsize,
        )
        title_defaults = merge(label_defaults, (padding=(12, 12, 8, 12),))
        for (row, line) in enumerate(stats_lines)
            defaults = row == 1 ? title_defaults : label_defaults
            Label(
                stats_grid[row, 1],
                line;
                _merged_kwargs(defaults, row == 1 ? stats_title_kwargs : stats_label_kwargs)...,
            )
        end
        Label(
            stats_grid[length(stats_lines) + 1, 1],
            "";
            halign=:left,
            valign=:top,
            justification=:left,
            padding=(0, 0, 12, 0),
            fontsize=1,
        )
        colsize!(fig.layout, 2, Fixed(panel_width_px))
        colsize!(fig.layout, 1, Auto(1))
    end

    if filename !== nothing
        outpath = String(filename)
        if isempty(splitext(outpath)[2])
            outpath *= ".$(String(format))"
        end
        save(outpath, fig)
    end

    return fig
end

function _theme_from_symbol(theme::Symbol, theme_override::Theme)
    thm = theme == :clean ? _clean_theme() : theme == :publication ? _publication_theme() : theme == :latex ? _latex_theme() : Theme()
    return merge(thm, theme_override)
end

"""
    plot_profile(profile_result; filename=nothing, format=:pdf, theme=:publication, ...)

Plot a one-dimensional profile-likelihood scan.
"""
function plot_profile(
    profile_result::ProfileResult;
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:publication,
    theme_override::Theme=Theme(),
    title="Profile",
    xlabel="parameter",
    ylabel="Delta cost",
    line_color=:dodgerblue4,
    line_width::Real=3,
    threshold_color=:black,
    figure_size::Tuple{<:Real, <:Real}=(900, 620),
    axis_kwargs=NamedTuple(),
    line_kwargs=NamedTuple(),
)
    fig = with_theme(_theme_from_symbol(theme, theme_override)) do
        Figure(size=(Int(round(figure_size[1])), Int(round(figure_size[2]))))
    end
    ax = Axis(fig[1, 1]; _merged_kwargs((title=title, xlabel=xlabel, ylabel=ylabel), axis_kwargs)...)
    lines!(
        ax,
        profile_result.values,
        profile_result.delta_cost;
        _merged_kwargs((color=line_color, linewidth=line_width), line_kwargs)...,
    )
    hlines!(ax, [profile_result.threshold]; color=threshold_color, linestyle=:dash)

    if filename !== nothing
        outpath = String(filename)
        if isempty(splitext(outpath)[2])
            outpath *= ".$(String(format))"
        end
        save(outpath, fig)
    end

    return fig
end

"""
    plot_contour(contour_result; filename=nothing, format=:pdf, theme=:publication, ...)

Plot a two-dimensional profile-likelihood contour grid.
"""
function plot_contour(
    contour_result::ContourResult;
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:publication,
    theme_override::Theme=Theme(),
    title="Contour",
    xlabel="parameter 1",
    ylabel="parameter 2",
    colormap=:viridis,
    line_color=:black,
    figure_size::Tuple{<:Real, <:Real}=(820, 700),
    axis_kwargs=NamedTuple(),
    heatmap_kwargs=NamedTuple(),
    contour_kwargs=NamedTuple(),
)
    fig = with_theme(_theme_from_symbol(theme, theme_override)) do
        Figure(size=(Int(round(figure_size[1])), Int(round(figure_size[2]))))
    end
    ax = Axis(fig[1, 1]; _merged_kwargs((title=title, xlabel=xlabel, ylabel=ylabel), axis_kwargs)...)
    hm = heatmap!(
        ax,
        contour_result.x_values,
        contour_result.y_values,
        contour_result.delta_cost;
        _merged_kwargs((colormap=colormap,), heatmap_kwargs)...,
    )
    contour!(
        ax,
        contour_result.x_values,
        contour_result.y_values,
        contour_result.delta_cost;
        _merged_kwargs((levels=contour_result.levels, color=line_color, linewidth=2), contour_kwargs)...,
    )
    Colorbar(fig[1, 2], hm; label="Delta cost")

    if filename !== nothing
        outpath = String(filename)
        if isempty(splitext(outpath)[2])
            outpath *= ".$(String(format))"
        end
        save(outpath, fig)
    end

    return fig
end

function _diagnostic_values(result::FitResult, kind::Symbol)
    x = result.problem.x
    yhat = result.model_y
    if kind == :residual
        return x, result.residuals, _yerror_for_plot(result.problem, result.params), "Residuals", "y - fit", 0.0
    elseif kind == :pull
        return x, _weighted_data_residual(result.problem, result.params), nothing, "Pulls", "pull", 0.0
    elseif kind == :ratio
        ratio = result.problem.y ./ yhat
        yerr = _yerror_for_plot(result.problem, result.params)
        ratio_err = yerr === nothing ? nothing : yerr ./ abs.(yhat)
        return x, ratio, ratio_err, "Ratio", "data / fit", 1.0
    end
    throw(ArgumentError("diagnostic plot kind must be :residual, :pull, or :ratio"))
end

"""
    plot_residuals(result; kind=:pull, filename=nothing, format=:pdf, ...)

Plot residuals, pulls, or data/fit ratios for an XY fit.
"""
function plot_residuals(
    result::FitResult;
    kind::Symbol=:pull,
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:publication,
    theme_override::Theme=Theme(),
    figure_size::Tuple{<:Real, <:Real}=(900, 520),
    xlabel="x",
    color=:black,
    marker=:circle,
    markersize::Real=9,
    axis_kwargs=NamedTuple(),
    scatter_kwargs=NamedTuple(),
    errorbars_kwargs=NamedTuple(),
)
    x, values, errors, title, ylabel, reference = _diagnostic_values(result, kind)
    fig = with_theme(_theme_from_symbol(theme, theme_override)) do
        Figure(size=(Int(round(figure_size[1])), Int(round(figure_size[2]))))
    end
    ax = Axis(fig[1, 1]; _merged_kwargs((title=title, xlabel=xlabel, ylabel=ylabel), axis_kwargs)...)
    hlines!(ax, [reference]; color=:gray35, linestyle=:dash)
    if errors !== nothing
        errorbars!(ax, x, values, errors; _merged_kwargs((color=color, whiskerwidth=8), errorbars_kwargs)...)
    end
    scatter!(ax, x, values; _merged_kwargs((color=color, marker=marker, markersize=markersize), scatter_kwargs)...)

    if filename !== nothing
        outpath = String(filename)
        if isempty(splitext(outpath)[2])
            outpath *= ".$(String(format))"
        end
        save(outpath, fig)
    end

    return fig
end

"""
    plot_diagnostics(result; filename=nothing, format=:pdf, ...)

Create a compact residual, pull, and ratio diagnostic figure for an XY fit.
"""
function plot_diagnostics(
    result::FitResult;
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:publication,
    theme_override::Theme=Theme(),
    figure_size::Tuple{<:Real, <:Real}=(900, 900),
    xlabel="x",
    color=:black,
    marker=:circle,
    markersize::Real=8,
    axis_kwargs=NamedTuple(),
)
    fig = with_theme(_theme_from_symbol(theme, theme_override)) do
        Figure(size=(Int(round(figure_size[1])), Int(round(figure_size[2]))))
    end

    for (row, kind) in enumerate((:residual, :pull, :ratio))
        x, values, errors, title, ylabel, reference = _diagnostic_values(result, kind)
        ax = Axis(fig[row, 1]; _merged_kwargs((title=title, xlabel=row == 3 ? xlabel : "", ylabel=ylabel), axis_kwargs)...)
        hlines!(ax, [reference]; color=:gray35, linestyle=:dash)
        if errors !== nothing
            errorbars!(ax, x, values, errors; color=color, whiskerwidth=8)
        end
        scatter!(ax, x, values; color=color, marker=marker, markersize=markersize)
    end

    if filename !== nothing
        outpath = String(filename)
        if isempty(splitext(outpath)[2])
            outpath *= ".$(String(format))"
        end
        save(outpath, fig)
    end

    return fig
end
