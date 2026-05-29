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
    band_label="1-sigma band",
    band_kwargs=NamedTuple(),
    xerr_color=:black,
    yerr_color=:black,
    error_whiskerwidth::Real=10,
    xerrorbars_kwargs=NamedTuple(),
    yerrorbars_kwargs=NamedTuple(),
    data_label="data",
)
    thm = theme == :publication ? _publication_theme() : theme == :latex ? _latex_theme() : Theme()
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
        xlabel=_as_label_text(xlabel, latex_labels),
        ylabel=_as_label_text(ylabel, latex_labels),
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
    sg = _prediction_band_sigma(result, xg)

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

    if show_legend
        axislegend(
            ax,
            [bplot, fplot, dplot],
            [
                _as_label_text(band_label, latex_labels),
                _as_label_text(fit_label, latex_labels),
                _as_label_text(data_label, latex_labels),
            ],
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
    thm = theme == :publication ? _publication_theme() : theme == :latex ? _latex_theme() : Theme()
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
