const _JF_PAPER = "#ffffff"
const _JF_PAPER_SOFT = "#ffffff"
const _JF_INK = "#14151a"
const _JF_MUTED = "#14151a"
const _JF_GRID = "#edf1f7"
const _JF_TEAL = "#0077b6"
const _JF_TEAL_SOFT = "#90c9e8"
const _JF_MINIMAL_PAPER = "#ffffff"
const _JF_MINIMAL_INK = "#101010"
const _JF_MINIMAL_MUTED = "#101010"
const _JF_MINIMAL_GRID = "#eef2f7"
const _JF_MINIMAL_FIT = "#0081a7"
const _JF_MINIMAL_BAND = "#a8dadc"
const _JF_PAPER_FIT = "#000000"
const _JF_PAPER_BAND = "#000000"
const _JF_DARK_PAPER = "#111318"
const _JF_DARK_PANEL = "#171b22"
const _JF_DARK_INK = "#edf2f4"
const _JF_DARK_MUTED = "#b8c1ca"
const _JF_DARK_GRID = "#2a313a"
const _JF_DARK_FIT = "#66d9ef"
const _JF_DARK_BAND = "#66d9ef"

function _fitplot_background(theme::Symbol)
    theme == :dark && return _JF_DARK_PAPER
    return :white
end

function _style_preset(theme::Symbol)
    if theme == :dark
        return (
            data_color=(_JF_DARK_INK, 0.74),
            data_marker=:circle,
            data_markersize=6.8,
            data_strokecolor=_JF_DARK_PAPER,
            data_strokewidth=0.35,
            fit_color=_JF_DARK_FIT,
            fit_linewidth=3.0,
            band_color=_JF_DARK_BAND,
            band_alpha=0.20,
            xerr_color=(_JF_DARK_MUTED, 0.42),
            yerr_color=(_JF_DARK_MUTED, 0.42),
            error_whiskerwidth=6.0,
            stats_color=_JF_DARK_INK,
            stats_muted_color=_JF_DARK_MUTED,
            stats_accent_color=_JF_DARK_FIT,
        )
    elseif theme == :minimal
        return (
            data_color=(_JF_MINIMAL_INK, 0.62),
            data_marker=:circle,
            data_markersize=6.5,
            data_strokecolor=_JF_MINIMAL_PAPER,
            data_strokewidth=0.25,
            fit_color=_JF_MINIMAL_FIT,
            fit_linewidth=2.8,
            band_color=_JF_MINIMAL_BAND,
            band_alpha=0.16,
            xerr_color=(_JF_MINIMAL_INK, 0.32),
            yerr_color=(_JF_MINIMAL_INK, 0.32),
            error_whiskerwidth=6.0,
            stats_color=_JF_MINIMAL_INK,
            stats_muted_color=_JF_MINIMAL_MUTED,
            stats_accent_color=_JF_MINIMAL_INK,
        )
    elseif theme in (:paper, :latex, :publication)
        return (
            data_color=:black,
            data_marker=:circle,
            data_markersize=5.6,
            data_strokecolor=:white,
            data_strokewidth=0.0,
            fit_color=_JF_PAPER_FIT,
            fit_linewidth=1.8,
            band_color=_JF_PAPER_BAND,
            band_alpha=0.09,
            xerr_color=(:black, 0.70),
            yerr_color=(:black, 0.70),
            error_whiskerwidth=4.0,
            stats_color=:black,
            stats_muted_color=:black,
            stats_accent_color=:black,
        )
    end

    return (
        data_color=(_JF_INK, 0.58),
        data_marker=:circle,
        data_markersize=8.2,
        data_strokecolor=_JF_INK,
        data_strokewidth=0.9,
        fit_color=_JF_TEAL,
        fit_linewidth=3.2,
        band_color=_JF_TEAL_SOFT,
        band_alpha=0.18,
        xerr_color=(_JF_INK, 0.38),
        yerr_color=(_JF_INK, 0.38),
        error_whiskerwidth=7,
        stats_color=_JF_INK,
        stats_muted_color=_JF_MUTED,
        stats_accent_color=_JF_INK,
    )
end

function _dark_theme()
    return Theme(
        fontsize=16,
        font="TeX Gyre Heros",
        figure_padding=(14, 18, 12, 14),
        Figure=(backgroundcolor=_JF_DARK_PAPER,),
        Axis=(
            xlabelsize=22,
            ylabelsize=22,
            titlesize=22,
            xticklabelsize=16,
            yticklabelsize=16,
            backgroundcolor=_JF_DARK_PAPER,
            xlabelcolor=_JF_DARK_INK,
            ylabelcolor=_JF_DARK_INK,
            titlecolor=_JF_DARK_INK,
            xticklabelcolor=_JF_DARK_MUTED,
            yticklabelcolor=_JF_DARK_MUTED,
            xtickcolor=_JF_DARK_MUTED,
            ytickcolor=_JF_DARK_MUTED,
            xgridvisible=true,
            ygridvisible=true,
            xminorgridvisible=false,
            yminorgridvisible=false,
            xgridcolor=(_JF_DARK_GRID, 0.82),
            ygridcolor=(_JF_DARK_GRID, 0.82),
            xgridwidth=0.55,
            ygridwidth=0.55,
            topspinevisible=false,
            rightspinevisible=false,
            leftspinecolor=_JF_DARK_MUTED,
            bottomspinecolor=_JF_DARK_MUTED,
        ),
        Legend=(framevisible=false, labelsize=16, labelcolor=_JF_DARK_INK, patchsize=(24, 13), rowgap=6),
        Lines=(linewidth=1.75,),
        Scatter=(markersize=4.2,),
    )
end

function _publication_theme()
    return Theme(
        fontsize=12,
        font="CMU Serif",
        figure_padding=(10, 12, 8, 10),
        Figure=(backgroundcolor=_JF_PAPER_SOFT,),
        Axis=(
            xlabelsize=15,
            ylabelsize=15,
            titlesize=15,
            xticklabelsize=12,
            yticklabelsize=12,
            backgroundcolor=_JF_PAPER_SOFT,
            xlabelcolor=_JF_INK,
            ylabelcolor=_JF_INK,
            titlecolor=_JF_INK,
            xticklabelcolor=_JF_INK,
            yticklabelcolor=_JF_INK,
            xtickcolor=_JF_MUTED,
            ytickcolor=_JF_MUTED,
            xgridvisible=false,
            ygridvisible=false,
            topspinevisible=false,
            rightspinevisible=false,
            leftspinecolor=_JF_INK,
            bottomspinecolor=_JF_INK,
        ),
        Legend=(
            framevisible=false,
            labelsize=12,
            labelcolor=_JF_INK,
            patchsize=(22, 12),
        ),
        Lines=(linewidth=1.8,),
        Scatter=(markersize=6,),
    )
end

function _clean_theme()
    return Theme(
        fontsize=15,
        font="TeX Gyre Heros",
        figure_padding=(12, 14, 10, 12),
        Figure=(backgroundcolor=_JF_PAPER,),
        Axis=(
            xlabelsize=19,
            ylabelsize=19,
            titlesize=19,
            xticklabelsize=14,
            yticklabelsize=14,
            backgroundcolor=_JF_PAPER,
            xlabelcolor=_JF_INK,
            ylabelcolor=_JF_INK,
            titlecolor=_JF_INK,
            xticklabelcolor=_JF_INK,
            yticklabelcolor=_JF_INK,
            xtickcolor=_JF_MUTED,
            ytickcolor=_JF_MUTED,
            xgridvisible=false,
            ygridvisible=false,
            xminorgridvisible=false,
            yminorgridvisible=false,
            xgridcolor=(_JF_GRID, 0.78),
            ygridcolor=(_JF_GRID, 0.78),
            xgridwidth=0.65,
            ygridwidth=0.65,
            topspinevisible=false,
            rightspinevisible=false,
            leftspinecolor=_JF_INK,
            bottomspinecolor=_JF_INK,
        ),
        Legend=(
            framevisible=false,
            labelsize=15,
            labelcolor=_JF_INK,
            patchsize=(24, 13),
            rowgap=6,
        ),
        Lines=(linewidth=2.4,),
        Scatter=(markersize=5.8,),
    )
end

function _minimal_theme()
    return Theme(
        fontsize=16,
        font="TeX Gyre Heros",
        figure_padding=(14, 18, 12, 14),
        Figure=(backgroundcolor=_JF_MINIMAL_PAPER,),
        Axis=(
            xlabelsize=22,
            ylabelsize=22,
            titlesize=22,
            xticklabelsize=16,
            yticklabelsize=16,
            backgroundcolor=_JF_MINIMAL_PAPER,
            xlabelcolor=_JF_MINIMAL_INK,
            ylabelcolor=_JF_MINIMAL_INK,
            titlecolor=_JF_MINIMAL_INK,
            xticklabelcolor=_JF_MINIMAL_INK,
            yticklabelcolor=_JF_MINIMAL_INK,
            xtickcolor=_JF_MINIMAL_MUTED,
            ytickcolor=_JF_MINIMAL_MUTED,
            xgridvisible=true,
            ygridvisible=true,
            xminorgridvisible=false,
            yminorgridvisible=false,
            xgridcolor=(_JF_MINIMAL_GRID, 0.78),
            ygridcolor=(_JF_MINIMAL_GRID, 0.78),
            xgridwidth=0.55,
            ygridwidth=0.55,
            topspinevisible=false,
            rightspinevisible=false,
            leftspinecolor=_JF_MINIMAL_INK,
            bottomspinecolor=_JF_MINIMAL_INK,
        ),
        Legend=(framevisible=false, labelsize=16, labelcolor=_JF_MINIMAL_INK, patchsize=(24, 13), rowgap=6),
        Lines=(linewidth=1.55,),
        Scatter=(markersize=4.2,),
    )
end

function _paper_theme()
    return Theme(
        fontsize=11,
        font="CMU Serif",
        figure_padding=(8, 10, 6, 8),
        Figure=(backgroundcolor=:white,),
        Axis=(
            xlabelsize=13,
            ylabelsize=13,
            titlesize=13,
            xticklabelsize=10,
            yticklabelsize=10,
            backgroundcolor=:white,
            xlabelcolor=:black,
            ylabelcolor=:black,
            titlecolor=:black,
            xticklabelcolor=:black,
            yticklabelcolor=:black,
            xtickcolor=:black,
            ytickcolor=:black,
            xgridvisible=false,
            ygridvisible=false,
            topspinevisible=false,
            rightspinevisible=false,
            leftspinecolor=:black,
            bottomspinecolor=:black,
        ),
        Legend=(framevisible=false, labelsize=11, labelcolor=:black, patchsize=(18, 10)),
        Lines=(linewidth=1.45,),
        Scatter=(markersize=4.8,),
    )
end

function _latex_theme()
    return Theme(
        fontsize=13,
        figure_padding=(10, 12, 8, 10),
        Figure=(backgroundcolor=_JF_PAPER_SOFT,),
        Axis=(
            xlabelsize=16,
            ylabelsize=16,
            titlesize=16,
            xticklabelsize=12,
            yticklabelsize=12,
            backgroundcolor=_JF_PAPER_SOFT,
            xlabelcolor=_JF_INK,
            ylabelcolor=_JF_INK,
            titlecolor=_JF_INK,
            xticklabelcolor=_JF_INK,
            yticklabelcolor=_JF_INK,
            xtickcolor=_JF_MUTED,
            ytickcolor=_JF_MUTED,
            xgridvisible=false,
            ygridvisible=false,
            topspinevisible=false,
            rightspinevisible=false,
            leftspinecolor=_JF_INK,
            bottomspinecolor=_JF_INK,
        ),
        font="CMU Serif",
        Legend=(framevisible=false, labelsize=12, labelcolor=_JF_INK, patchsize=(22, 12)),
        Lines=(linewidth=1.8,),
        Scatter=(markersize=6,),
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

function _interpolate_sigma_to_grid(x::AbstractVector, sigma::AbstractVector, xgrid::AbstractVector)
    length(x) == length(sigma) || throw(ArgumentError("uncertainty length must match x length"))
    order = sortperm(x)
    xs = Float64.(x[order])
    ss = Float64.(sigma[order])

    length(xs) == 1 && return fill(ss[1], length(xgrid))
    out = Vector{Float64}(undef, length(xgrid))
    j = 1
    for (i, xg) in pairs(xgrid)
        if xg <= xs[1]
            out[i] = ss[1]
        elseif xg >= xs[end]
            out[i] = ss[end]
        else
            while j < length(xs) - 1 && xs[j + 1] < xg
                j += 1
            end
            if xs[j + 1] == xs[j]
                out[i] = ss[j]
                continue
            end
            weight = (xg - xs[j]) / (xs[j + 1] - xs[j])
            out[i] = (1 - weight) * ss[j] + weight * ss[j + 1]
        end
    end
    return out
end

function _observation_band_sigma(result::FitResult, xgrid::AbstractVector)
    problem = result.problem
    variance = zeros(Float64, length(xgrid))

    yerr = _yerror_for_plot(problem, result.params)
    if yerr !== nothing
        ygrid_err = _interpolate_sigma_to_grid(problem.x, yerr, xgrid)
        variance .+= ygrid_err .^ 2
    end

    xerr = _xerror_for_plot(problem, result.params)
    if xerr !== nothing
        xgrid_err = _interpolate_sigma_to_grid(problem.x, xerr, xgrid)
        dydx = _model_dydx(problem, result.params; x=xgrid)
        variance .+= (dydx .* xgrid_err) .^ 2
    end

    return sqrt.(clamp.(variance, 0.0, Inf))
end

function _fit_band_sigma(result::FitResult, xgrid::AbstractVector, band::Symbol)
    parameter_sigma = _prediction_band_sigma(result, xgrid)
    band == :confidence && return parameter_sigma
    band == :prediction && return sqrt.(parameter_sigma .^ 2 .+ _observation_band_sigma(result, xgrid) .^ 2)
    return zeros(Float64, length(xgrid))
end

function _panel_width_px(stats_panel_width, fig_width::Int, stats_lines)
    if stats_panel_width === :auto
        max_chars = maximum(length(string(line)) for line in stats_lines; init=24)
        return clamp(130 + 8 * max_chars, 300, 520)
    end

    stats_panel_width isa Real || throw(ArgumentError("stats_panel_width must be :auto or a positive number"))
    stats_panel_width > 0 || throw(ArgumentError("stats_panel_width must be positive"))
    if stats_panel_width <= 1
        return clamp(Int(round(fig_width * stats_panel_width)), 300, 560)
    end
    return Int(round(stats_panel_width))
end

function _label_with_unit(label, unit)
    label === nothing && return ""
    unit === nothing && return label
    unit_text = string(unit)
    isempty(unit_text) && return label
    label_text = string(label)
    return isempty(label_text) ? unit_text : string(label_text, " (", unit_text, ")")
end

function _plot_title(title)
    title === nothing && return ""
    return title
end

function _default_model_label(result::FitResult)
    result.problem.model === _default_linear_model && return L"y = p_1 x + p_2"
    return nothing
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
    latex_stats::Bool=false,
    stats_mode::Symbol=:compact,
)
    stats_mode in (:compact, :full) || throw(ArgumentError("stats_mode must be :compact or :full"))
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
        push!(lines, L"\chi^2/\mathrm{n_{dof}} = %$chi2_ndf_text")
        stats_mode == :full && push!(lines, L"\chi^2 = %$chi2_text")
        push!(lines, L"\mathrm{n_{dof}} = %$ndf_text")
    else
        push!(lines, string("chi2/ndf = ", chi2_ndf_text))
        stats_mode == :full && push!(lines, string("chi2 = ", chi2_text))
        push!(lines, string("ndf = ", ndf_text))
    end

    return lines
end

function _plain_stats_rows(
    result::FitResult;
    parameter_names::Union{Nothing, AbstractVector}=nothing,
    sigdigits::Int=4,
    latex_stats::Bool=false,
    stats_mode::Symbol=:compact,
)
    stats_mode in (:compact, :full) || throw(ArgumentError("stats_mode must be :compact or :full"))
    n = length(result.params)
    names = if parameter_names === nothing
        Any["p$i" for i in 1:n]
    else
        length(parameter_names) == n || throw(ArgumentError("parameter_names length must match parameter count"))
        collect(parameter_names)
    end

    rows = Tuple{Any, Any}[]
    for i in 1:n
        value = _fmt_value(result.params[i]; sigdigits=sigdigits)
        err = _fmt_value(result.param_stderr[i]; sigdigits=sigdigits)
        name = latex_stats ? LaTeXString(_latex_symbol_expr(names[i])) : string(names[i])
        uncertainty = latex_stats ? LaTeXString(value * " \\pm " * err) : string(value, " ± ", err)
        push!(rows, (name, uncertainty))
    end

    push!(rows, (latex_stats ? L"\chi^2" : "χ²", _fmt_value(result.stats.chi2; sigdigits=sigdigits)))
    push!(rows, (latex_stats ? L"\chi^2/\mathrm{ndf}" : "χ²/ndf", _fmt_value(result.stats.chi2_ndf; sigdigits=sigdigits)))
    push!(rows, (latex_stats ? L"P(\chi^2)" : "χ² prob.", _fmt_value(result.stats.pvalue; sigdigits=sigdigits)))
    push!(rows, ("ndf", string(result.stats.ndf)))
    if stats_mode == :full
        push!(rows, ("cost", _fmt_value(result.stats.cost_min; sigdigits=sigdigits)))
        push!(rows, ("AIC", _fmt_value(result.stats.aic; sigdigits=sigdigits)))
        push!(rows, ("BIC", _fmt_value(result.stats.bic; sigdigits=sigdigits)))
    end
    return rows
end

function _stats_row_texts(stats_rows)
    return [string(name, " ", value) for (name, value) in stats_rows]
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

function _stats_box_geometry(stats_lines; position::Symbol=:lt)
    max_chars = maximum(length(string(line)) for line in stats_lines; init=20)
    width = clamp(0.013 * max_chars, 0.22, 0.34)
    height = clamp(0.031 * length(stats_lines) + 0.035, 0.13, 0.36)
    margin = 0.025

    if position in (:lt, :lefttop)
        return margin, 1.0 - margin, width, height, :left, :top
    elseif position in (:lb, :leftbottom)
        return margin, margin + height, width, height, :left, :top
    elseif position in (:rt, :righttop)
        return 1.0 - margin - width, 1.0 - margin, width, height, :left, :top
    elseif position in (:rb, :rightbottom)
        return 1.0 - margin - width, margin + height, width, height, :left, :top
    end

    throw(ArgumentError("inside stats position must be :lt, :lb, :rt, or :rb"))
end

function _draw_inside_stats!(
    ax,
    stats_lines;
    position::Symbol=:lt,
    fontsize::Real=14,
    box_color=_JF_PAPER_SOFT,
    box_alpha::Real=0.90,
    box_strokecolor=_JF_GRID,
    box_strokewidth::Real=1.0,
    text_color=_JF_INK,
)
    x, y, width, height, halign, valign = _stats_box_geometry(stats_lines; position=position)
    rect = Point2f[
        (x, y),
        (x + width, y),
        (x + width, y - height),
        (x, y - height),
    ]
    poly!(
        ax,
        rect;
        space=:relative,
        color=(box_color, box_alpha),
        strokecolor=box_strokecolor,
        strokewidth=box_strokewidth,
    )
    text!(
        ax,
        x + 0.014,
        y - 0.018;
        text=join(string.(stats_lines), "\n"),
        space=:relative,
        align=(halign, valign),
        fontsize=fontsize,
        color=text_color,
        justification=:left,
        lineheight=0.95,
    )
    return nothing
end

function _draw_right_stats!(
    fig,
    stats_rows;
    panel_width_px::Int,
    fontsize::Real=15,
    title=nothing,
    model_label=nothing,
    color=_JF_INK,
    muted_color=_JF_MUTED,
    accent_color=_JF_TEAL,
    legend_plots=nothing,
    legend_labels=nothing,
    legend_kwargs=NamedTuple(),
)
    panel_cell = if legend_plots !== nothing && legend_labels !== nothing
        panel_grid = GridLayout(fig[1, 2])
        Legend(
            panel_grid[1, 1],
            legend_plots,
            legend_labels;
            _merged_kwargs(
                (
                    framevisible=false,
                    tellwidth=false,
                    tellheight=true,
                    labelsize=fontsize + 1,
                    patchsize=(26, 14),
                    rowgap=8,
                ),
                legend_kwargs,
            )...,
        )
        rowgap!(panel_grid, 14)
        panel_grid[2, 1]
    else
        fig[1, 2]
    end

    panel_ax = Axis(
        panel_cell;
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=false,
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        bottomspinevisible=false,
        leftspinevisible=false,
        backgroundcolor=:transparent,
    )
    xlims!(panel_ax, 0, 1)
    ylims!(panel_ax, 0, 1)
    y = 0.965
    if title !== nothing && !isempty(string(title))
        text!(
            panel_ax,
            0.0,
            y;
            text=title,
            space=:relative,
            align=(:left, :top),
            fontsize=fontsize + 1,
            color=color,
        )
        y -= 0.070
    end
    if model_label !== nothing && !isempty(string(model_label))
        text!(
            panel_ax,
            0.0,
            y;
            text=model_label,
            space=:relative,
            align=(:left, :top),
            fontsize=fontsize + 1,
            color=color,
        )
        y -= 0.095
    end
    y -= 0.012

    step = clamp(0.0042 * Float64(fontsize), 0.062, 0.084)
    value_x = 0.50
    for (i, (name, value)) in enumerate(stats_rows)
        row_y = y - (i - 1) * step
        text!(
            panel_ax,
            0.0,
            row_y;
            text=name,
            space=:relative,
            align=(:left, :top),
            fontsize=fontsize,
            color=color,
        )
        text!(
            panel_ax,
            value_x,
            row_y;
            text=value,
            space=:relative,
            align=(:left, :top),
            fontsize=fontsize,
            color=color,
        )
    end
    colsize!(fig.layout, 2, Fixed(panel_width_px))
    colsize!(fig.layout, 1, Auto(1))
    return nothing
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
        theme=:clean,
        theme_override=Theme(),
        title=nothing,
        model_label=nothing,
        model_label=nothing,
        xlabel="x",
        ylabel="y",
        xunit=nothing,
        yunit=nothing,
        auto_limits=true,
        limit_padding=0.08,
        plot_aspect=nothing,
        figure_size=nothing,
        stats_panel_width=:auto,
        stats_position=:right,
        inside_stats_position=:lt,
        panel_gap=22,
        latex_labels=false,
        latex_stats=false,
        show_stats=true,
        stats_mode=:compact,
        tight_layout=true,
        stats_sigdigits=5,
        parameter_names=nothing,
        stats_fontsize=15,
        stats_title=nothing,
        stats_title=nothing,
        stats_box_color=_JF_PAPER_SOFT,
        stats_box_alpha=0.95,
        stats_box_strokecolor=_JF_GRID,
        stats_box_strokewidth=1.0,
        stats_linegap=2,
        stats_label_kwargs=NamedTuple(),
        stats_title_kwargs=NamedTuple(),
        stats_box_kwargs=NamedTuple(),
        show_legend=false,
        legend_position=:rt,
        axis_kwargs=NamedTuple(),
        legend_kwargs=NamedTuple(),
        data_color=nothing,
        data_marker=nothing,
        data_markersize=nothing,
        data_strokecolor=nothing,
        data_strokewidth=nothing,
        data_strokecolor=nothing,
        data_strokewidth=nothing,
        scatter_kwargs=NamedTuple(),
        fit_color=nothing,
        fit_linewidth=nothing,
        fit_label="fit",
        line_kwargs=NamedTuple(),
        band_color=nothing,
        band_alpha=nothing,
        band=:confidence,
        nsigma=1,
        band_label="1-sigma band",
        band_kwargs=NamedTuple(),
        xerr_color=nothing,
        yerr_color=nothing,
        error_whiskerwidth=nothing,
        xerrorbars_kwargs=NamedTuple(),
        yerrorbars_kwargs=NamedTuple(),
        data_label="data",
    )

Create a publication-style fit plot with data, error bars, best-fit curve, optional
uncertainty band, and an optional right-side statistics panel. `band=:confidence`
shows the propagated parameter-covariance band. `band=:prediction` additionally
includes observation uncertainty in y and effective x uncertainty. Makie keyword
containers can be passed as `NamedTuple`s or `Dict`s via the `*_kwargs` arguments.
"""
function plot_fit(
    result::FitResult;
    xgrid=nothing,
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:clean,
    theme_override::Theme=Theme(),
    title=nothing,
    model_label=nothing,
    xlabel="x",
    ylabel="y",
    xunit=nothing,
    yunit=nothing,
    auto_limits::Bool=true,
    limit_padding::Real=0.08,
    plot_aspect::Union{Nothing, Real}=nothing,
    figure_size::Union{Nothing, Tuple{<:Real, <:Real}}=nothing,
    stats_panel_width=:auto,
    stats_position::Symbol=:right,
    inside_stats_position::Symbol=:lt,
    panel_gap::Real=22,
    latex_labels::Bool=false,
    latex_stats::Bool=false,
    show_stats::Bool=true,
    stats_mode::Symbol=:compact,
    tight_layout::Bool=true,
    stats_sigdigits::Int=5,
    parameter_names::Union{Nothing, AbstractVector}=nothing,
    stats_fontsize::Real=15,
    stats_title=nothing,
    stats_box_color=_JF_PAPER_SOFT,
    stats_box_alpha::Real=0.95,
    stats_box_strokecolor=_JF_GRID,
    stats_box_strokewidth::Real=1.0,
    stats_linegap::Real=2,
    stats_label_kwargs=NamedTuple(),
    stats_title_kwargs=NamedTuple(),
    stats_box_kwargs=NamedTuple(),
    show_legend::Bool=false,
    legend_position=:rt,
    axis_kwargs=NamedTuple(),
    legend_kwargs=NamedTuple(),
    data_color=nothing,
    data_marker=nothing,
    data_markersize=nothing,
    data_strokecolor=nothing,
    data_strokewidth=nothing,
    scatter_kwargs=NamedTuple(),
    fit_color=nothing,
    fit_linewidth=nothing,
    fit_label="fit",
    line_kwargs=NamedTuple(),
    band_color=nothing,
    band_alpha=nothing,
    band::Symbol=:confidence,
    nsigma::Real=1.0,
    band_label="1-sigma band",
    band_kwargs=NamedTuple(),
    xerr_color=nothing,
    yerr_color=nothing,
    error_whiskerwidth=nothing,
    xerrorbars_kwargs=NamedTuple(),
    yerrorbars_kwargs=NamedTuple(),
    data_label="data",
)
    band in (:confidence, :prediction, :none) || throw(ArgumentError("band must be :confidence, :prediction, or :none"))
    stats_position in (:inside, :right) || throw(ArgumentError("stats_position must be :inside or :right"))

    thm = theme == :clean ? _clean_theme() :
        theme == :minimal ? _minimal_theme() :
        theme == :dark ? _dark_theme() :
        theme == :paper ? _paper_theme() :
        theme == :publication ? _publication_theme() :
        theme == :latex ? _latex_theme() :
        Theme()
    thm = merge(thm, theme_override)
    style = _style_preset(theme)
    data_color = data_color === nothing ? style.data_color : data_color
    data_marker = data_marker === nothing ? style.data_marker : data_marker
    data_markersize = data_markersize === nothing ? style.data_markersize : data_markersize
    data_strokecolor = data_strokecolor === nothing ? style.data_strokecolor : data_strokecolor
    data_strokewidth = data_strokewidth === nothing ? style.data_strokewidth : data_strokewidth
    fit_color = fit_color === nothing ? style.fit_color : fit_color
    fit_linewidth = fit_linewidth === nothing ? style.fit_linewidth : fit_linewidth
    band_color = band_color === nothing ? style.band_color : band_color
    band_alpha = band_alpha === nothing ? style.band_alpha : band_alpha
    xerr_color = xerr_color === nothing ? style.xerr_color : xerr_color
    yerr_color = yerr_color === nothing ? style.yerr_color : yerr_color
    error_whiskerwidth = error_whiskerwidth === nothing ? style.error_whiskerwidth : error_whiskerwidth
    model_label = model_label === nothing ? _default_model_label(result) : model_label

    base_size = show_stats && stats_position == :right ? (1220, 720) : (980, 640)
    fig_size = figure_size === nothing ? base_size : (Int(round(figure_size[1])), Int(round(figure_size[2])))
    fig = with_theme(thm) do
        Figure(size=fig_size, backgroundcolor=_fitplot_background(theme))
    end

    if show_stats
        colgap!(fig.layout, Int(round(panel_gap)))
    end

    axis_defaults = (
        title=_as_label_text(_plot_title(title), latex_labels),
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
    sg = Float64(nsigma) .* _fit_band_sigma(result, xg, band)

    bplot = nothing
    if band != :none
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
            (
                color=data_color,
                marker=data_marker,
                markersize=data_markersize,
                strokecolor=data_strokecolor,
                strokewidth=data_strokewidth,
                label=_as_label_text(data_label, latex_labels),
            ),
            scatter_kwargs,
        )...,
    )

    if auto_limits
        xlims, ylims = _fit_plot_limits(x, y, xerr, yerr, xg, yg, band != :none ? sg : nothing; padding=limit_padding)
        xlims !== nothing && ylims !== nothing && limits!(ax, xlims..., ylims...)
    end

    legend_plots = band != :none ? [dplot, fplot, bplot] : [dplot, fplot]
    legend_labels = band != :none ?
        [
            _as_label_text(data_label, latex_labels),
            _as_label_text(fit_label, latex_labels),
            _as_label_text(band_label, latex_labels),
        ] :
        [
            _as_label_text(data_label, latex_labels),
            _as_label_text(fit_label, latex_labels),
        ]

    if show_legend && !(show_stats && stats_position == :right)
        axislegend(
            ax,
            legend_plots,
            legend_labels,
            ;
            _merged_kwargs((position=legend_position,), legend_kwargs)...,
        )
    end

    stats_lines = show_stats ? _stats_panel_lines(
        result;
        parameter_names=parameter_names,
        sigdigits=stats_sigdigits,
        latex_stats=stats_position == :right && latex_stats,
        stats_mode=stats_mode,
    ) : nothing
    right_stats_rows = show_stats ? _plain_stats_rows(
        result;
        parameter_names=parameter_names,
        sigdigits=stats_sigdigits,
        latex_stats=latex_stats,
        stats_mode=stats_mode,
    ) : nothing

    if show_stats && stats_position == :inside
        _draw_inside_stats!(
            ax,
            stats_lines;
            position=inside_stats_position,
            fontsize=stats_fontsize,
            box_color=stats_box_color,
            box_alpha=stats_box_alpha,
            box_strokecolor=stats_box_strokecolor,
            box_strokewidth=stats_box_strokewidth,
        )
    elseif show_stats && stats_position == :right
        panel_width_px = _panel_width_px(stats_panel_width, fig_size[1], _stats_row_texts(right_stats_rows))
        _draw_right_stats!(
            fig,
            right_stats_rows;
            panel_width_px=panel_width_px,
            fontsize=stats_fontsize,
            title=stats_title,
            model_label=model_label,
            color=style.stats_color,
            muted_color=style.stats_muted_color,
            accent_color=style.stats_accent_color,
            legend_plots=show_legend ? legend_plots : nothing,
            legend_labels=show_legend ? legend_labels : nothing,
            legend_kwargs=legend_kwargs,
        )
    end

    tight_layout && resize_to_layout!(fig)

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
    thm = theme == :clean ? _clean_theme() :
        theme == :minimal ? _minimal_theme() :
        theme == :dark ? _dark_theme() :
        theme == :paper ? _paper_theme() :
        theme == :publication ? _publication_theme() :
        theme == :latex ? _latex_theme() :
        Theme()
    return merge(thm, theme_override)
end

"""
    plot_profile(profile_result; filename=nothing, format=:pdf, theme=:publication, ...)

Plot a one-dimensional profile-likelihood scan.

Pass `local_sigma=result.param_stderr[i]` to overlay the local parabolic
covariance approximation. If the profile and parabola disagree visibly, local
symmetric errors should not be treated as the final uncertainty statement.
Use `delta_max` to focus the view on scientifically relevant interval
thresholds when a strongly non-parabolic tail would otherwise compress the
minimum.
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
    local_sigma=nothing,
    local_color=:gray35,
    local_linewidth::Real=2,
    local_linestyle=:dash,
    threshold_color=:black,
    show_legend::Bool=true,
    profile_label="profile cost",
    local_label="local covariance parabola",
    threshold_label=nothing,
    delta_max::Union{Nothing, Real}=nothing,
    figure_size::Tuple{<:Real, <:Real}=(900, 620),
    axis_kwargs=NamedTuple(),
    line_kwargs=NamedTuple(),
    local_line_kwargs=NamedTuple(),
)
    fig = with_theme(_theme_from_symbol(theme, theme_override)) do
        Figure(size=(Int(round(figure_size[1])), Int(round(figure_size[2]))), backgroundcolor=_fitplot_background(theme))
    end
    ax = Axis(fig[1, 1]; _merged_kwargs((title=title, xlabel=xlabel, ylabel=ylabel), axis_kwargs)...)
    lines!(
        ax,
        profile_result.values,
        profile_result.delta_cost;
        _merged_kwargs((color=line_color, linewidth=line_width), line_kwargs)...,
    )
    if local_sigma !== nothing
        sigma = Float64(local_sigma)
        sigma > 0 || throw(ArgumentError("local_sigma must be positive"))
        local_delta = @. abs2((profile_result.values - profile_result.best_value) / sigma)
        lines!(
            ax,
            profile_result.values,
            local_delta;
            _merged_kwargs((color=local_color, linewidth=local_linewidth, linestyle=local_linestyle), local_line_kwargs)...,
        )
    end
    hlines!(ax, [profile_result.threshold]; color=threshold_color, linestyle=:dash)

    if delta_max !== nothing
        delta_limit = Float64(delta_max)
        delta_limit > 0 || throw(ArgumentError("delta_max must be positive"))
        ylims!(ax, 0, delta_limit)
    end

    if show_legend
        handles = Any[LineElement(color=line_color, linewidth=line_width)]
        labels = Any[profile_label]
        if local_sigma !== nothing
            push!(handles, LineElement(color=local_color, linewidth=local_linewidth, linestyle=local_linestyle))
            push!(labels, local_label)
        end
        push!(handles, LineElement(color=threshold_color, linewidth=2, linestyle=:dash))
        push!(
            labels,
            threshold_label === nothing ?
            "interval threshold (Δcost = $(_fmt_value(profile_result.threshold; sigdigits=3)))" :
            threshold_label,
        )
        Legend(fig[1, 2], handles, labels; framevisible=false)
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

"""
    plot_contour(contour_result; filename=nothing, format=:pdf, theme=:publication, ...)

Plot a two-dimensional profile-likelihood contour grid.

The default emphasizes directly interpretable profile regions and labeled
two-parameter 1-sigma/2-sigma contour levels. Set `show_heatmap=true` for an
explicit delta-cost surface view; heatmaps are not the default diagnostic
because they make confidence thresholds harder to read quickly.

Pass `local_covariance=result.param_covariance` and
`local_center=result.params[[i, j]]` to overlay the local covariance ellipse in
the same parameter plane. Non-elliptic profile contours indicate that local
Gaussian covariance errors are not sufficient.
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
    show_heatmap::Bool=false,
    show_regions::Bool=true,
    show_legend::Bool=true,
    level_colors=("#4c78a8", "#72b7b2", "#f2cf5b"),
    region_colors=(("#4c78a8", 0.24), ("#72b7b2", 0.16), ("#f2cf5b", 0.12)),
    line_color=nothing,
    local_covariance=nothing,
    local_center=nothing,
    local_line_color=:gray35,
    local_linewidth::Real=2,
    local_linestyle=:dash,
    figure_size::Tuple{<:Real, <:Real}=(820, 700),
    axis_kwargs=NamedTuple(),
    heatmap_kwargs=NamedTuple(),
    contour_kwargs=NamedTuple(),
    local_contour_kwargs=NamedTuple(),
)
    fig = with_theme(_theme_from_symbol(theme, theme_override)) do
        Figure(size=(Int(round(figure_size[1])), Int(round(figure_size[2]))), backgroundcolor=_fitplot_background(theme))
    end
    ax = Axis(fig[1, 1]; _merged_kwargs((title=title, xlabel=xlabel, ylabel=ylabel), axis_kwargs)...)
    hm = if show_heatmap
        heatmap!(
            ax,
            contour_result.x_values,
            contour_result.y_values,
            contour_result.delta_cost;
            _merged_kwargs((colormap=colormap, interpolate=false, rasterize=true), heatmap_kwargs)...,
        )
    end

    plot_levels = sort(unique(contour_result.levels))
    isempty(plot_levels) && throw(ArgumentError("contour_result levels must not be empty"))
    all(isfinite, plot_levels) || throw(ArgumentError("contour_result levels must be finite"))
    all(>(0.0), plot_levels) || throw(ArgumentError("contour_result levels must be positive delta-cost thresholds"))

    colors = collect(level_colors)
    isempty(colors) && throw(ArgumentError("level_colors must contain at least one color"))
    regions = collect(region_colors)
    if show_regions && !show_heatmap
        isempty(regions) && throw(ArgumentError("region_colors must contain at least one color when show_regions=true"))
        contourf!(
            ax,
            contour_result.x_values,
            contour_result.y_values,
            contour_result.delta_cost;
            levels=vcat(0.0, plot_levels),
            colormap=regions,
        )
    end

    contour_level_label(level) =
        isapprox(level, 2.30; atol=0.015) ? "profile 1σ (2 parameters)" :
        isapprox(level, 6.18; atol=0.015) ? "profile 2σ (2 parameters)" :
        "profile Δcost = $(_fmt_value(level; sigdigits=3))"
    contour_handles = Any[]
    contour_labels = String[]
    for (index, level) in enumerate(plot_levels)
        color = line_color === nothing ? colors[mod1(index, length(colors))] : line_color
        contour!(
            ax,
            contour_result.x_values,
            contour_result.y_values,
            contour_result.delta_cost;
            _merged_kwargs((levels=[level], color=color, linewidth=3), contour_kwargs)...,
        )
        push!(contour_handles, LineElement(color=color, linewidth=3))
        push!(contour_labels, contour_level_label(level))
    end

    local_handles = Any[]
    local_labels = String[]
    if local_covariance !== nothing
        cov = Matrix{Float64}(local_covariance)
        if size(cov) != (2, 2)
            i, j = contour_result.parameter_indices
            size(cov, 1) >= max(i, j) && size(cov, 2) >= max(i, j) ||
                throw(ArgumentError("local_covariance must be 2x2 or the full parameter covariance matrix"))
            cov = cov[[i, j], [i, j]]
        end
        center = if local_center === nothing
            (mean(contour_result.x_values), mean(contour_result.y_values))
        else
            raw_center = collect(local_center)
            length(raw_center) == 2 || throw(ArgumentError("local_center must contain exactly two values"))
            (Float64(raw_center[1]), Float64(raw_center[2]))
        end
        precision = Symmetric(cov) \ Matrix{Float64}(I, 2, 2)
        local_delta = Matrix{Float64}(undef, length(contour_result.x_values), length(contour_result.y_values))
        for ix in eachindex(contour_result.x_values), iy in eachindex(contour_result.y_values)
            delta = [contour_result.x_values[ix] - center[1], contour_result.y_values[iy] - center[2]]
            local_delta[ix, iy] = dot(delta, precision * delta)
        end
        for level in plot_levels
            contour!(
                ax,
                contour_result.x_values,
                contour_result.y_values,
                local_delta;
                _merged_kwargs(
                    (levels=[level], color=local_line_color, linewidth=local_linewidth, linestyle=local_linestyle),
                    local_contour_kwargs,
                )...,
            )
            push!(local_handles, LineElement(color=local_line_color, linewidth=local_linewidth, linestyle=local_linestyle))
            push!(local_labels, "local covariance")
        end
    end
    finite_surface = replace(vec(contour_result.delta_cost), NaN => Inf)
    any(isfinite, finite_surface) ||
        throw(ArgumentError("contour_result must contain at least one finite cost value"))
    center_index = argmin(finite_surface)
    center_cartesian = CartesianIndices(contour_result.delta_cost)[center_index]
    scatter!(
        ax,
        [contour_result.x_values[center_cartesian[1]]],
        [contour_result.y_values[center_cartesian[2]]];
        marker=:cross,
        markersize=18,
        strokewidth=3,
        color=line_color === nothing ? colors[1] : line_color,
    )

    if show_heatmap
        Colorbar(fig[1, 2], hm; label="Δ cost")
    elseif show_legend
        handles = Any[MarkerElement(
            marker=:cross,
            markersize=18,
            strokewidth=3,
            color=line_color === nothing ? colors[1] : line_color,
        )]
        labels = ["profile minimum"]
        append!(handles, contour_handles)
        append!(labels, contour_labels)
        if !isempty(local_handles)
            push!(handles, first(local_handles))
            push!(labels, first(local_labels))
        end
        Legend(fig[1, 2], handles, labels; framevisible=false)
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
