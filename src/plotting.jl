const _JF_PAPER = "#ffffff"
const _JF_PAPER_SOFT = "#fbfcfd"
const _JF_INK = "#17191f"
const _JF_MUTED = "#303842"
const _JF_GRID = "#c8cdd3"
const _JF_DARK_PAPER = "#111318"
const _JF_DARK_INK = "#edf2f4"
const _JF_DARK_MUTED = "#d3dae0"
const _JF_DARK_GRID = "#73808d"

const _JF_STYLE_ALIASES = Dict(
    :analysis => :sans,
    :presentation => :sans,
    :screen => :sans,
    :lab => :sans,
    :workbench => :sans,
    :modern => :sans,
    :clean => :sans,
    :minimal => :sans,
    :showcase => :sans,
    :article => :tex,
    :publication => :tex,
    :paper => :tex,
    :latex => :tex,
)

function _resolve_plot_style(theme::Symbol, appearance::Symbol)
    appearance in (:auto, :light, :dark) ||
        throw(ArgumentError("appearance must be :auto, :light, or :dark"))

    style = get(_JF_STYLE_ALIASES, theme, theme)
    style in (:sans, :tex) ||
        throw(ArgumentError(
            "theme must be :sans, :tex, or a supported legacy alias",
        ))
    return style, appearance == :auto ? :light : appearance
end

function _style_preset(style::Symbol, appearance::Symbol)
    dark = appearance == :dark
    paper = dark ? _JF_DARK_PAPER : _JF_PAPER
    ink = dark ? _JF_DARK_INK : _JF_INK
    grid = dark ? _JF_DARK_GRID : _JF_GRID
    box = dark ? "#1b2027" : _JF_PAPER_SOFT
    box_stroke = dark ? "#65717d" : _JF_GRID

    if style == :sans
        # Screen-oriented sans typography, open axes, and light grid guides.
        # Information density is controlled independently by `show_panel`.
        fit = :dodgerblue
        return (
            name=:sans,
            diagnostic_scale=0.84,
            background_color=paper,
            axis_color=ink,
            grid_color=(grid, dark ? 0.34 : 0.46),
            data_color=ink,
            data_marker=:circle,
            data_markersize=11.0,
            data_strokecolor=paper,
            data_strokewidth=1.2,
            fit_color=fit,
            fit_linewidth=3.3,
            band_color=fit,
            band_alpha=dark ? 0.24 : 0.20,
            xerr_color=(ink, dark ? 0.84 : 0.76),
            yerr_color=(ink, dark ? 0.84 : 0.76),
            error_whiskerwidth=4.0,
            secondary_color=dark ? "#ff6b6b" : :red,
            reference_color=dark ? :slategray1 : :gray35,
            series_colors=dark ?
                (fit, "#ff6b6b", "#70cfa8", :mediumpurple1, :slategray1) :
                (fit, :red, :seagreen4, :slateblue, :gray35),
            stats_color=ink,
            stats_muted_color=ink,
            stats_fontsize=24,
            stats_box_color=box,
            stats_box_strokecolor=box_stroke,
            fontsize=23,
            xlabelsize=29,
            ylabelsize=29,
            titlesize=36,
            titlegap=18,
            subplot_titlesize=28,
            subplot_titlegap=8,
            titlealign=:left,
            ticklabelsize=23,
            legend_labelsize=23,
            legend_patchsize=(34, 20),
            legend_rowgap=5,
            panel_rowgap=2,
            panel_sectiongap=10,
            panel_gap=22,
            figure_padding=(14, 18, 14, 14),
            figure_size_with_panel=(1040, 640),
            figure_size_without_panel=(860, 560),
            xgridvisible=true,
            ygridvisible=true,
            gridwidth=1.0,
            topspinevisible=false,
            rightspinevisible=false,
            spinewidth=1.9,
            tickwidth=1.7,
            ticksize=7.0,
            tickalign=0.0,
        )
    elseif style == :tex
        fit = dark ? "#77bce6" : "#0072b2"
        return (
            name=:tex,
            diagnostic_scale=0.90,
            background_color=paper,
            axis_color=ink,
            grid_color=(grid, 0.0),
            # Hollow observations remain separable from fitted curves in
            # grayscale and in dense vector exports.
            data_color=paper,
            data_marker=:circle,
            data_markersize=11.0,
            data_strokecolor=ink,
            data_strokewidth=1.7,
            fit_color=fit,
            fit_linewidth=3.0,
            band_color=fit,
            band_alpha=dark ? 0.18 : 0.15,
            xerr_color=(ink, dark ? 0.82 : 0.76),
            yerr_color=(ink, dark ? 0.82 : 0.76),
            error_whiskerwidth=4.5,
            secondary_color=dark ? "#f06b4f" : "#d55e00",
            reference_color=dark ? "#d3dae0" : "#444a50",
            series_colors=dark ?
                (fit, "#f06b4f", "#70cfa8", "#d69acb", "#a9b4bf") :
                (fit, "#d55e00", "#009e73", "#cc79a7", "#4d4d4d"),
            stats_color=ink,
            stats_muted_color=ink,
            stats_fontsize=30,
            stats_box_color=box,
            stats_box_strokecolor=box_stroke,
            fontsize=26,
            xlabelsize=32,
            ylabelsize=32,
            titlesize=38,
            titlegap=22,
            subplot_titlesize=30,
            subplot_titlegap=10,
            # Long scientific titles must never clip against a narrow axis.
            titlealign=:left,
            ticklabelsize=26,
            legend_labelsize=26,
            legend_patchsize=(36, 21),
            legend_rowgap=5,
            panel_rowgap=2,
            panel_sectiongap=10,
            panel_gap=20,
            figure_padding=(14, 18, 13, 13),
            figure_size_with_panel=(1000, 640),
            figure_size_without_panel=(760, 520),
            xgridvisible=false,
            ygridvisible=false,
            gridwidth=0.0,
            topspinevisible=true,
            rightspinevisible=true,
            spinewidth=1.7,
            tickwidth=1.6,
            ticksize=7.0,
            tickalign=1.0,
        )
    end

    throw(ArgumentError(
        "style presets are defined only for :sans and :tex",
    ))
end

function _theme_from_style(style::Symbol, appearance::Symbol, theme_override::Theme)
    preset = _style_preset(style, appearance)
    font_theme = style == :tex ? theme_latexfonts() : Theme(font="TeX Gyre Heros")
    visual_theme = Theme(
        fontsize=preset.fontsize,
        figure_padding=preset.figure_padding,
        Figure=(backgroundcolor=preset.background_color,),
        Axis=(
            xlabelsize=preset.xlabelsize,
            ylabelsize=preset.ylabelsize,
            titlesize=preset.titlesize,
            titlegap=preset.titlegap,
            titlealign=preset.titlealign,
            xticklabelsize=preset.ticklabelsize,
            yticklabelsize=preset.ticklabelsize,
            backgroundcolor=preset.background_color,
            xlabelcolor=preset.axis_color,
            ylabelcolor=preset.axis_color,
            titlecolor=preset.axis_color,
            xticklabelcolor=preset.axis_color,
            yticklabelcolor=preset.axis_color,
            xtickcolor=preset.axis_color,
            ytickcolor=preset.axis_color,
            xtickwidth=preset.tickwidth,
            ytickwidth=preset.tickwidth,
            xticksize=preset.ticksize,
            yticksize=preset.ticksize,
            xtickalign=preset.tickalign,
            ytickalign=preset.tickalign,
            xgridvisible=preset.xgridvisible,
            ygridvisible=preset.ygridvisible,
            xminorgridvisible=false,
            yminorgridvisible=false,
            xgridcolor=preset.grid_color,
            ygridcolor=preset.grid_color,
            xgridwidth=preset.gridwidth,
            ygridwidth=preset.gridwidth,
            topspinevisible=preset.topspinevisible,
            rightspinevisible=preset.rightspinevisible,
            leftspinecolor=preset.axis_color,
            rightspinecolor=preset.axis_color,
            topspinecolor=preset.axis_color,
            bottomspinecolor=preset.axis_color,
            spinewidth=preset.spinewidth,
        ),
        Legend=(
            framevisible=false,
            labelsize=preset.legend_labelsize,
            labelcolor=preset.axis_color,
            patchsize=preset.legend_patchsize,
            rowgap=preset.legend_rowgap,
        ),
        Lines=(linewidth=preset.fit_linewidth,),
        Scatter=(
            marker=preset.data_marker,
            markersize=preset.data_markersize,
            strokecolor=preset.data_strokecolor,
            strokewidth=preset.data_strokewidth,
        ),
    )
    return merge(font_theme, visual_theme, theme_override)
end

"""
    plot_theme(theme=:sans; appearance=:auto, theme_override=Theme())

Return the Makie theme used by JuFitter plots. Use this when composing a custom
Makie figure that should remain visually consistent with `plot_fit`.
"""
function plot_theme(
    theme::Symbol=:sans;
    appearance::Symbol=:auto,
    theme_override::Theme=Theme(),
)
    style, resolved_appearance = _resolve_plot_style(theme, appearance)
    return _theme_from_style(style, resolved_appearance, theme_override)
end

"""
    plot_palette(theme=:sans; appearance=:auto)

Return the visual tokens used by a JuFitter plot style. Besides data, fit,
uncertainty-band, and error-bar defaults, the result exposes typography,
layout, and color-safe multi-series tokens for custom Makie figures.
"""
function plot_palette(theme::Symbol=:sans; appearance::Symbol=:auto)
    style, resolved_appearance = _resolve_plot_style(theme, appearance)
    return _style_preset(style, resolved_appearance)
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
    if band == :prediction && result.problem.whitening !== nothing &&
       result.problem.whitening.marginal_sigma === nothing
        throw(ArgumentError(
            "band=:prediction requires marginal_sigma in WhiteningOperator; " *
            "provide marginal standard deviations or use band=:confidence",
        ))
    end
    band == :prediction && return sqrt.(parameter_sigma .^ 2 .+ _observation_band_sigma(result, xgrid) .^ 2)
    return zeros(Float64, length(xgrid))
end

function _panel_width_px(stats_panel_width, fig_width::Int)
    if stats_panel_width === :auto
        return nothing
    end

    stats_panel_width isa Real || throw(ArgumentError("stats_panel_width must be :auto or a positive number"))
    stats_panel_width > 0 || throw(ArgumentError("stats_panel_width must be positive"))
    if stats_panel_width <= 1
        return clamp(Int(round(fig_width * stats_panel_width)), 300, 560)
    end
    return Int(round(stats_panel_width))
end

function _apply_right_panel_sizing!(fig, panel_width_px)
    if panel_width_px === nothing
        colsize!(fig.layout, 2, Auto())
    else
        colsize!(fig.layout, 2, Fixed(panel_width_px))
    end
    colsize!(fig.layout, 1, Auto(1))
    return nothing
end

function _label_with_unit(label, unit)
    label === nothing && return ""
    unit === nothing && return label
    unit_text = string(unit)
    isempty(unit_text) && return label
    label_text = string(label)
    isempty(label_text) && return unit

    # The SI quantity-calculus form makes tick values explicit: a tick at 2 on
    # an axis `t / s` means t/s = 2. Preserve math input as one LaTeX expression
    # so mixed quantity and unit strings do not produce nested `$...$` pairs.
    if label isa LaTeXString || unit isa LaTeXString
        quantity = label isa LaTeXString ? _strip_math_delims(label_text) :
            "\\text{" * replace(label_text, "_" => "\\_") * "}"
        unit_expr = unit isa LaTeXString ? _strip_math_delims(unit_text) :
            "\\mathrm{" * replace(unit_text, "_" => "\\_") * "}"
        return LaTeXString(quantity * "\\,/\\," * unit_expr)
    end
    return string(label_text, " / ", unit_text)
end

function _plot_title(title)
    title === nothing && return ""
    return title
end

function _default_model_label(result::FitResult, latex_labels::Bool)
    if result.problem.model === _default_linear_model
        return latex_labels ? L"y = p_1 x + p_2" : "y = p₁ x + p₂"
    end
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

function _fit_x_limits(x::AbstractVector, xerr; padding::Real)
    x_values = Float64[]
    append!(x_values, Float64.(x))
    if xerr !== nothing
        append!(x_values, Float64.(x .- xerr))
        append!(x_values, Float64.(x .+ xerr))
    end
    return _padded_limits(x_values; padding=padding)
end

function _default_grid(x::AbstractVector; xlimits=nothing)
    xlo, xhi = xlimits === nothing ? (minimum(x), maximum(x)) : xlimits
    if xlo == xhi
        xlo -= 0.5
        xhi += 0.5
    end
    return collect(range(xlo, xhi; length=800))
end

function _latex_number_expr(value::AbstractString)
    value == "NaN" && return "\\mathrm{NaN}"
    value == "Inf" && return "\\infty"
    value == "-Inf" && return "-\\infty"
    scientific = match(r"^(.+)[eE]([+-]?\d+)$", value)
    scientific === nothing && return value
    mantissa, exponent = scientific.captures
    return mantissa * "\\times10^{" * string(parse(Int, exponent)) * "}"
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
            push!(lines, LaTeXString(
                pexpr * " = " * _latex_number_expr(v) * " \\pm " * _latex_number_expr(e),
            ))
        else
            push!(lines, string(names[i], " = ", v, " ± ", e))
        end
    end

    chi2_text = _fmt_value(result.stats.chi2; sigdigits=sigdigits)
    chi2_ndf_text = _fmt_value(result.stats.chi2_ndf; sigdigits=sigdigits)
    ndf_text = string(result.stats.ndf)

    if latex_stats
        push!(lines, LaTeXString(
            "\\chi^2/\\mathrm{n_{dof}} = " * _latex_number_expr(chi2_ndf_text),
        ))
        stats_mode == :full && push!(lines, LaTeXString(
            "\\chi^2 = " * _latex_number_expr(chi2_text),
        ))
        push!(lines, LaTeXString("\\mathrm{n_{dof}} = " * ndf_text))
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
        uncertainty = latex_stats ?
            LaTeXString(_latex_number_expr(value) * " \\pm " * _latex_number_expr(err)) :
            string(value, " ± ", err)
        push!(rows, (name, uncertainty))
    end

    stat_value(value) = latex_stats ? LaTeXString(_latex_number_expr(string(value))) : string(value)
    push!(rows, (
        latex_stats ? L"\chi^2" : "χ²",
        stat_value(_fmt_value(result.stats.chi2; sigdigits=sigdigits)),
    ))
    push!(rows, (
        latex_stats ? L"\chi^2/\mathrm{ndf}" : "χ²/ndf",
        stat_value(_fmt_value(result.stats.chi2_ndf; sigdigits=sigdigits)),
    ))
    push!(rows, (
        latex_stats ? L"P(\chi^2)" : "χ² probability",
        stat_value(_fmt_value(result.stats.pvalue; sigdigits=sigdigits)),
    ))
    push!(rows, (latex_stats ? L"\mathrm{ndf}" : "ndf", stat_value(result.stats.ndf)))
    if stats_mode == :full
        push!(rows, (
            latex_stats ? L"\mathrm{cost}" : "cost",
            stat_value(_fmt_value(result.stats.cost_min; sigdigits=sigdigits)),
        ))
        push!(rows, (
            latex_stats ? L"\mathrm{AIC}" : "AIC",
            stat_value(_fmt_value(result.stats.aic; sigdigits=sigdigits)),
        ))
        push!(rows, (
            latex_stats ? L"\mathrm{BIC}" : "BIC",
            stat_value(_fmt_value(result.stats.bic; sigdigits=sigdigits)),
        ))
    end
    return rows
end

function _stats_row_line(name, value, latex_stats::Bool)
    latex_stats || return string(name, " = ", value)
    name_expr = _strip_math_delims(String(name))
    value_expr = _strip_math_delims(String(value))
    return LaTeXString(name_expr * " = " * value_expr)
end

function _stats_row_lines(stats_rows, latex_stats::Bool)
    return [_stats_row_line(name, value, latex_stats) for (name, value) in stats_rows]
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
        isempty(s) && return ""
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
    fontsize::Real=16,
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

"""
    plot_info_panel!(
        cell;
        theme=:sans,
        appearance=:auto,
        legend_source=nothing,
        legend_plots=nothing,
        legend_labels=nothing,
        model_label=nothing,
        parameter_lines=Any[],
        statistic_lines=Any[],
        ...
    )

Add a compact, left-aligned information panel to a Makie layout cell. The
panel is intended for custom scientific figures that should use the same
legend, model, parameter, and statistic hierarchy as `plot_fit`. `theme` and
`appearance` supply readable panel defaults from the same central style
contract; explicit panel keywords remain authoritative.
"""
function plot_info_panel!(
    cell;
    theme::Symbol=:sans,
    appearance::Symbol=:auto,
    legend_source=nothing,
    legend_plots=nothing,
    legend_labels=nothing,
    parameter_lines=Any[],
    statistic_lines=Any[],
    fontsize::Union{Nothing, Real}=nothing,
    title=nothing,
    model_label=nothing,
    color=nothing,
    muted_color=nothing,
    legend_kwargs=NamedTuple(),
    tellwidth::Bool=true,
    tellheight::Bool=false,
)
    style, resolved_appearance = _resolve_plot_style(theme, appearance)
    preset = _style_preset(style, resolved_appearance)
    fontsize = fontsize === nothing ? preset.stats_fontsize : Float64(fontsize)
    color = color === nothing ? preset.stats_color : color
    muted_color = muted_color === nothing ? preset.stats_muted_color : muted_color

    # The panel must not dictate the height of its parent row. Otherwise a
    # compact report vertically centers and shrinks the adjacent data axis.
    panel_grid = GridLayout(cell; tellwidth=tellwidth, tellheight=tellheight, valign=:top)
    row = 1
    generous_gaps = Int[]

    legend_defaults = (
        framevisible=false,
        tellwidth=true,
        tellheight=true,
        halign=:left,
        valign=:top,
        labelsize=fontsize,
        patchsize=preset.legend_patchsize,
        rowgap=preset.legend_rowgap,
    )
    if legend_source !== nothing
        Legend(panel_grid[row, 1], legend_source; _merged_kwargs(legend_defaults, legend_kwargs)...)
        push!(generous_gaps, row)
        row += 1
    elseif legend_plots !== nothing && legend_labels !== nothing
        Legend(
            panel_grid[row, 1],
            legend_plots,
            legend_labels;
            _merged_kwargs(legend_defaults, legend_kwargs)...,
        )
        push!(generous_gaps, row)
        row += 1
    end

    if title !== nothing && !isempty(string(title))
        Label(
            panel_grid[row, 1],
            title;
            halign=:left,
            tellwidth=true,
            fontsize=fontsize + 1,
            color=color,
        )
        row += 1
    end

    if model_label !== nothing && !isempty(string(model_label))
        Label(
            panel_grid[row, 1],
            model_label;
            halign=:left,
            tellwidth=true,
            fontsize=fontsize + 1,
            color=color,
        )
        push!(generous_gaps, row)
        row += 1
    end

    for line in parameter_lines
        Label(
            panel_grid[row, 1],
            line;
            halign=:left,
            tellwidth=true,
            fontsize=fontsize,
            color=color,
        )
        row += 1
    end

    !isempty(parameter_lines) && !isempty(statistic_lines) && push!(generous_gaps, row - 1)
    for line in statistic_lines
        Label(
            panel_grid[row, 1],
            line;
            halign=:left,
            tellwidth=true,
            fontsize=fontsize,
            color=muted_color,
        )
        row += 1
    end

    row > 2 && rowgap!(panel_grid, preset.panel_rowgap)
    for gap_index in generous_gaps
        gap_index < row - 1 || continue
        rowgap!(panel_grid, gap_index, preset.panel_sectiongap)
    end
    return panel_grid
end

function _draw_right_stats!(
    fig,
    stats_lines;
    panel_width_px::Union{Nothing, Int},
    parameter_count::Int,
    kwargs...,
)
    plot_info_panel!(
        fig[1, 2];
        parameter_lines=stats_lines[1:parameter_count],
        statistic_lines=stats_lines[(parameter_count + 1):end],
        kwargs...,
    )
    _apply_right_panel_sizing!(fig, panel_width_px)
    return nothing
end

const _FITPLOT_FIT_KWARGS = Set([
    :sigma_y,
    :sigma_x,
    :cov_y,
    :cov_x,
    :whitening,
    :error_components,
    :bounds,
    :constraints,
    :parameter_priors,
    :parameter_constraints,
    :fixed_parameters,
    :jacobian,
    :x_derivative,
    :inplace,
    :backend,
    :cost,
    :maxiters,
    :tol,
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

function _fitplot_result(result::FitResult; print_report::Bool=false, kwargs...)
    plot_kwargs = Dict{Symbol, Any}(pairs(kwargs))
    parameter_names = get(plot_kwargs, :parameter_names, nothing)

    if print_report
        println(report_text(result; parameter_names=parameter_names))
    end

    fig = plot_fit(result; (; plot_kwargs...)...)
    return (result=result, figure=fig)
end

"""
    fitplot(result::FitResult; show_panel=true, print_report=false, kwargs...)
    fitplot(model, x, y; p0, show_panel=true, print_report=false, kwargs...)
    fitplot(x, y; p0=nothing, show_panel=true, print_report=false, kwargs...)

Fit and plot in one call. The `model, x, y` method forwards fitting keywords
such as `sigma_y`, `sigma_x`, `whitening`, `x_derivative`, `bounds`,
`parameter_priors`, and `backend` to `fit_model`; plotting keywords such as
`xlabel`, `ylabel`, `theme`, `nsigma`, and `filename` are forwarded to
`plot_fit`. `show_panel` controls the numerical panel in the figure, while
`print_report` independently controls terminal output.

The `x, y` method uses a linear model by default. All methods return a named
tuple `(result, figure)` so the numerical result is not lost.
"""
function fitplot(
    result::FitResult;
    show_panel::Bool=true,
    print_report::Bool=false,
    kwargs...,
)
    return _fitplot_result(
        result;
        show_panel=show_panel,
        print_report=print_report,
        kwargs...,
    )
end

function fitplot(
    model,
    x::AbstractVector,
    y::AbstractVector;
    p0::AbstractVector,
    show_panel::Bool=true,
    print_report::Bool=false,
    kwargs...,
)
    fit_kwargs, plot_kwargs = _split_fitplot_kwargs(kwargs)
    result = fit_model(model, x, y; p0=p0, fit_kwargs...)
    return _fitplot_result(
        result;
        show_panel=show_panel,
        print_report=print_report,
        plot_kwargs...,
    )
end

function fitplot(
    x::AbstractVector,
    y::AbstractVector;
    p0=nothing,
    show_panel::Bool=true,
    print_report::Bool=false,
    kwargs...,
)
    initial = p0 === nothing ? _linear_initial_guess(x, y) : p0
    return fitplot(
        _default_linear_model,
        x,
        y;
        p0=initial,
        show_panel=show_panel,
        print_report=print_report,
        kwargs...,
    )
end

"""
    plot_fit(
        result::FitResult;
        xgrid=nothing,
        filename=nothing,
        format=:pdf,
        theme=:sans,
        appearance=:auto,
        theme_override=Theme(),
        title=nothing,
        model_label=nothing,
        xlabel="x",
        ylabel="y",
        xunit=nothing,
        yunit=nothing,
        auto_limits=true,
        limit_padding=0.08,
        fit_range=:axis,
        plot_aspect=nothing,
        figure_size=nothing,
        stats_panel_width=:auto,
        stats_position=:right,
        inside_stats_position=:lt,
        panel_gap=nothing,
        latex_labels=false,
        latex_stats=false,
        show_panel=true,
        stats_mode=:compact,
        tight_layout=true,
        stats_sigdigits=5,
        parameter_names=nothing,
        stats_fontsize=nothing,
        stats_title=nothing,
        stats_box_color=nothing,
        stats_box_alpha=0.95,
        stats_box_strokecolor=nothing,
        stats_box_strokewidth=1.0,
        show_legend=true,
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

Create a scientific fit plot with data, error bars, best-fit curve, optional
uncertainty band, and an optional right-side information panel. Use
`theme=:sans` for open axes, sans-serif typography, and grid guides, or
`theme=:tex` for TeX typography, a full frame, and no grid. `show_panel`
independently controls the numerical result panel and defaults to `true` for
both styles. Legacy style names remain compatibility aliases rather than
additional cosmetic presets. `appearance=:light` or `:dark` controls the color
scheme independently.
`band=:confidence`
shows the propagated parameter-covariance band. `band=:prediction` additionally
includes observation uncertainty in y and effective x uncertainty. With the
default `fit_range=:axis`, the automatically generated model grid extends to
the padded axis range; use `fit_range=:data` or pass `xgrid` to draw only over a
specific domain. Makie keyword containers can be passed as `NamedTuple`s or
`Dict`s via the `*_kwargs` arguments.
"""
function plot_fit(
    result::FitResult;
    xgrid=nothing,
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:sans,
    appearance::Symbol=:auto,
    theme_override::Theme=Theme(),
    title=nothing,
    model_label=nothing,
    xlabel="x",
    ylabel="y",
    xunit=nothing,
    yunit=nothing,
    auto_limits::Bool=true,
    limit_padding::Real=0.08,
    fit_range::Symbol=:axis,
    plot_aspect::Union{Nothing, Real}=nothing,
    figure_size::Union{Nothing, Tuple{<:Real, <:Real}}=nothing,
    stats_panel_width=:auto,
    stats_position::Symbol=:right,
    inside_stats_position::Symbol=:lt,
    panel_gap::Union{Nothing, Real}=nothing,
    latex_labels::Bool=false,
    latex_stats::Bool=false,
    show_panel::Bool=true,
    stats_mode::Symbol=:compact,
    tight_layout::Bool=true,
    stats_sigdigits::Int=5,
    parameter_names::Union{Nothing, AbstractVector}=nothing,
    stats_fontsize::Union{Nothing, Real}=nothing,
    stats_title=nothing,
    stats_box_color=nothing,
    stats_box_alpha::Real=0.95,
    stats_box_strokecolor=nothing,
    stats_box_strokewidth::Real=1.0,
    show_legend::Bool=true,
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
    fit_range in (:axis, :data) || throw(ArgumentError("fit_range must be :axis or :data"))
    isfinite(nsigma) && nsigma > 0 || throw(ArgumentError("nsigma must be finite and positive"))
    isfinite(limit_padding) && limit_padding >= 0 ||
        throw(ArgumentError("limit_padding must be finite and non-negative"))

    resolved_style, resolved_appearance = _resolve_plot_style(theme, appearance)
    thm = _theme_from_style(resolved_style, resolved_appearance, theme_override)
    style = _style_preset(resolved_style, resolved_appearance)
    panel_gap = panel_gap === nothing ? style.panel_gap : Float64(panel_gap)
    stats_fontsize = stats_fontsize === nothing ? style.stats_fontsize : Float64(stats_fontsize)
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
    stats_box_color = stats_box_color === nothing ? style.stats_box_color : stats_box_color
    stats_box_strokecolor = stats_box_strokecolor === nothing ?
        style.stats_box_strokecolor : stats_box_strokecolor
    model_label = model_label === nothing ? _default_model_label(result, latex_labels) : model_label

    base_size = show_panel && stats_position == :right ?
        style.figure_size_with_panel : style.figure_size_without_panel
    fig_size = figure_size === nothing ? base_size : (Int(round(figure_size[1])), Int(round(figure_size[2])))
    fig = with_theme(thm) do
        Figure(size=fig_size, backgroundcolor=style.background_color)
    end

    if show_panel
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

    grid_limits = xgrid === nothing && fit_range == :axis ? _fit_x_limits(x, xerr; padding=limit_padding) : nothing
    xg = xgrid === nothing ? _default_grid(x; xlimits=grid_limits) : collect(Float64, xgrid)
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

    if show_legend && !(show_panel && stats_position == :right)
        axislegend(
            ax,
            legend_plots,
            legend_labels,
            ;
            _merged_kwargs((position=legend_position,), legend_kwargs)...,
        )
    end

    stats_lines = show_panel ? _stats_panel_lines(
        result;
        parameter_names=parameter_names,
        sigdigits=stats_sigdigits,
        latex_stats=stats_position == :right && latex_stats,
        stats_mode=stats_mode,
    ) : nothing
    right_stats_rows = show_panel ? _plain_stats_rows(
        result;
        parameter_names=parameter_names,
        sigdigits=stats_sigdigits,
        latex_stats=latex_stats,
        stats_mode=stats_mode,
    ) : nothing
    right_stats_lines = show_panel ? _stats_row_lines(right_stats_rows, latex_stats) : nothing

    if show_panel && stats_position == :inside
        _draw_inside_stats!(
            ax,
            stats_lines;
            position=inside_stats_position,
            fontsize=stats_fontsize,
            box_color=stats_box_color,
            box_alpha=stats_box_alpha,
            box_strokecolor=stats_box_strokecolor,
            box_strokewidth=stats_box_strokewidth,
            text_color=style.stats_color,
        )
    elseif show_panel && stats_position == :right
        panel_width_px = _panel_width_px(stats_panel_width, fig_size[1])
        _draw_right_stats!(
            fig,
            right_stats_lines;
            panel_width_px=panel_width_px,
            parameter_count=length(result.params),
            fontsize=stats_fontsize,
            title=stats_title,
            model_label=model_label,
            color=style.stats_color,
            muted_color=style.stats_muted_color,
            legend_plots=show_legend ? legend_plots : nothing,
            legend_labels=show_legend ? legend_labels : nothing,
            legend_kwargs=legend_kwargs,
            theme=resolved_style,
            appearance=resolved_appearance,
        )
    end

    # Preserve the declared output footprint. Resizing to the layout makes a
    # compact side panel collapse the figure height and destabilizes custom
    # elements added after `plot_fit` returns.
    tight_layout && trim!(fig.layout)

    if filename !== nothing
        outpath = String(filename)
        if isempty(splitext(outpath)[2])
            outpath *= ".$(String(format))"
        end
        save(outpath, fig)
    end

    return fig
end

function _axis_limits(axis::Axis)
    rect = axis.finallimits[]
    xmin = rect.origin[1]
    ymin = rect.origin[2]
    xmax = xmin + rect.widths[1]
    ymax = ymin + rect.widths[2]
    return (xmin, xmax), (ymin, ymax)
end

"""
    fit_axis(figure; index=1)

Return the `index`-th Makie `Axis` stored in a JuFitter figure.

This is a small convenience for post-fit annotation workflows:
`fig = plot_fit(result); ax = fit_axis(fig); add_vline!(ax, x0)`. It searches
the figure layout instead of relying on manual cell indices, so the same call
works for ordinary fit plots with or without a right-side report.
"""
function fit_axis(figure::Figure; index::Integer=1)
    index >= 1 || throw(ArgumentError("index must be >= 1"))
    axes = [object for object in contents(figure.layout) if object isa Axis]
    index <= length(axes) || throw(ArgumentError("figure contains only $(length(axes)) axis object(s)"))
    return axes[Int(index)]
end

function _finite_annotation_vector(name::AbstractString, values; min_length::Int)
    vector = values isa Number ? [Float64(values)] : collect(Float64, values)
    length(vector) >= min_length || throw(ArgumentError("$name must contain at least $min_length value(s)"))
    all(isfinite, vector) || throw(ArgumentError("$name must contain only finite values"))
    return vector
end

"""
    add_curve!(axis, f; xgrid=nothing, xspan=nothing, n=400, label=nothing, kwargs...)
    add_curve!(axis, x, y; label=nothing, kwargs...)

Add a curve to an existing Makie axis and return the created plot object.

The function-valued method samples `f` either on `xgrid`, on `xspan=(xmin,
xmax)`, or on the current visible x-range of `axis`. This is useful for adding
extrapolations, reference models, or derived physical relationships after a fit
has already been computed.
"""
function add_curve!(
    axis::Axis,
    f;
    xgrid=nothing,
    xspan=nothing,
    n::Integer=400,
    label=nothing,
    kwargs...,
)
    n >= 2 || throw(ArgumentError("n must be >= 2"))
    xs = if xgrid !== nothing
        _finite_annotation_vector("xgrid", xgrid; min_length=2)
    elseif xspan !== nothing
        length(xspan) == 2 || throw(ArgumentError("xspan must contain exactly two values"))
        span = _finite_annotation_vector("xspan", xspan; min_length=2)
        collect(range(span[1], span[2]; length=Int(n)))
    else
        limits, _ = _axis_limits(axis)
        collect(range(limits[1], limits[2]; length=Int(n)))
    end
    ys = collect(Float64, [f(x) for x in xs])
    all(isfinite, ys) || throw(ArgumentError("curve values must contain only finite values"))
    return lines!(axis, xs, ys; label=label, kwargs...)
end

function add_curve!(axis::Axis, x::AbstractVector, y::AbstractVector; label=nothing, kwargs...)
    xs = _finite_annotation_vector("x", x; min_length=2)
    ys = _finite_annotation_vector("y", y; min_length=2)
    length(xs) == length(ys) || throw(ArgumentError("x and y must have equal length"))
    return lines!(axis, xs, ys; label=label, kwargs...)
end

"""
    add_points!(axis, x, y; label=nothing, kwargs...)

Add marker points to an existing fit axis. This is intended for derived
quantities, thresholds, extrapolated intersections, or highlighted data points;
it does not rerun or modify the fit.
"""
function add_points!(axis::Axis, x, y; label=nothing, kwargs...)
    xs = _finite_annotation_vector("x", x; min_length=1)
    ys = _finite_annotation_vector("y", y; min_length=1)
    length(xs) == length(ys) || throw(ArgumentError("x and y must have equal length"))
    return scatter!(axis, xs, ys; label=label, kwargs...)
end

"""
    add_vline!(axis, x; label=nothing, kwargs...)

Add vertical reference line(s) to an existing axis. This wraps Makie's
`vlines!` with JuFitter-style argument validation.
"""
function add_vline!(axis::Axis, x; label=nothing, kwargs...)
    xs = _finite_annotation_vector("x", x; min_length=1)
    return vlines!(axis, xs; label=label, kwargs...)
end

"""
    add_hline!(axis, y; label=nothing, kwargs...)

Add horizontal reference line(s) to an existing axis. This wraps Makie's
`hlines!` with JuFitter-style argument validation.
"""
function add_hline!(axis::Axis, y; label=nothing, kwargs...)
    ys = _finite_annotation_vector("y", y; min_length=1)
    return hlines!(axis, ys; label=label, kwargs...)
end

"""
    add_vband!(axis, xmin, xmax; label=nothing, kwargs...)

Add a vertical uncertainty/reference band to an existing axis. The band spans
the full axis height without contributing artificial y values to automatic
limits, so it remains an annotation layer rather than a new data model.
"""
function add_vband!(axis::Axis, xmin::Real, xmax::Real; label=nothing, kwargs...)
    isfinite(xmin) && isfinite(xmax) || throw(ArgumentError("xmin and xmax must be finite"))
    xmin <= xmax || throw(ArgumentError("xmin must be <= xmax"))
    return vspan!(axis, Float64(xmin), Float64(xmax); label=label, kwargs...)
end

"""
    add_hband!(axis, ymin, ymax; label=nothing, kwargs...)

Add a horizontal uncertainty/reference band to an existing axis. The band spans
the full axis width without contributing artificial x values to automatic
limits, so it remains an annotation layer rather than a new data model.
"""
function add_hband!(axis::Axis, ymin::Real, ymax::Real; label=nothing, kwargs...)
    isfinite(ymin) && isfinite(ymax) || throw(ArgumentError("ymin and ymax must be finite"))
    ymin <= ymax || throw(ArgumentError("ymin must be <= ymax"))
    return hspan!(axis, Float64(ymin), Float64(ymax); label=label, kwargs...)
end

function _diagnostic_colors(style::Symbol, appearance::Symbol)
    preset = _style_preset(style, appearance)
    alpha = style == :tex ? (0.20, 0.11, 0.06) : (0.28, 0.16, 0.08)
    return (
        levels=preset.series_colors[1:3],
        regions=ntuple(index -> (preset.fit_color, alpha[index]), 3),
        local_color=preset.reference_color,
    )
end

_contour_level_name(level) =
    isapprox(level, 2.30; atol=0.015) ? "1σ" :
    isapprox(level, 6.18; atol=0.015) ? "2σ" : "Δcost"

function _contour_level_label(level)
    name = _contour_level_name(level)
    level_text = _fmt_value(level; sigdigits=3)
    return name == "Δcost" ? "profile contour Δcost = $level_text" :
           "profile $name region (2 params, Δcost = $level_text)"
end

function _local_contour_label(levels)
    names = _contour_level_name.(levels)
    level_text = all(!=("Δcost"), names) ? "$(join(names, "/")) " : ""
    return "local covariance $(level_text)contours (parabolic approximation)"
end

function _panel_status_color(status::Symbol, appearance::Symbol)
    dark = appearance == :dark
    status == :stop && return dark ? "#ff8a80" : "#a33a2b"
    status == :review && return dark ? "#d3dae0" : "#4b5560"
    return dark ? _JF_DARK_MUTED : _JF_MUTED
end

function _panel_status_label(status::Symbol)
    status == :stop && return "fix"
    status == :review && return "inspect"
    status == :ok && return "ok"
    return string(status)
end

function _validate_panel_status_mode(mode::Symbol)
    mode in (:issues, :all, :none) ||
        throw(ArgumentError("panel_status_mode must be :issues, :all, or :none"))
    return nothing
end

function _draw_panel_status!(
    axis::Axis,
    status::Symbol,
    style,
    appearance::Symbol;
    mode::Symbol=:issues,
    fontsize::Union{Nothing, Real}=nothing,
)
    _validate_panel_status_mode(mode)
    mode == :none && return nothing
    mode == :issues && status == :ok && return nothing
    text!(
        axis,
        0.035,
        0.94;
        text=_panel_status_label(status),
        space=:relative,
        align=(:left, :top),
        color=_panel_status_color(status, appearance),
        fontsize=fontsize === nothing ? max(16, style.ticklabelsize - 2) : Float64(fontsize),
        font=:bold,
    )
    return nothing
end

"""
    plot_profile(profile_result; filename=nothing, format=:pdf, theme=:sans, ...)

Plot a one-dimensional profile-likelihood scan.

Pass `local_sigma=result.param_stderr[i]` to overlay the local parabolic
covariance approximation. If the profile and parabola disagree visibly, local
symmetric errors should not be treated as the final uncertainty statement.
Use `delta_max` to focus the view on scientifically relevant interval
thresholds when a strongly non-parabolic tail would otherwise compress the
minimum. Line weights and colors follow `theme`; explicit `line_kwargs`,
`local_line_kwargs`, and `threshold_kwargs` take precedence.
"""
function plot_profile(
    profile_result::ProfileResult;
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:sans,
    appearance::Symbol=:auto,
    theme_override::Theme=Theme(),
    title="Profile",
    xlabel="parameter",
    ylabel="Delta cost",
    line_color=nothing,
    line_width::Union{Nothing, Real}=nothing,
    local_sigma=nothing,
    local_color=nothing,
    local_linewidth::Union{Nothing, Real}=nothing,
    local_linestyle=:dash,
    threshold_color=nothing,
    show_legend::Bool=true,
    profile_label="profile cost",
    local_label="local covariance parabola",
    threshold_label=nothing,
    delta_max::Union{Nothing, Real}=nothing,
    figure_size::Tuple{<:Real, <:Real}=(900, 620),
    axis_kwargs=NamedTuple(),
    line_kwargs=NamedTuple(),
    local_line_kwargs=NamedTuple(),
    threshold_kwargs=NamedTuple(),
)
    resolved_style, resolved_appearance = _resolve_plot_style(theme, appearance)
    style = _style_preset(resolved_style, resolved_appearance)
    line_color = line_color === nothing ? style.fit_color : line_color
    line_width = line_width === nothing ? style.fit_linewidth : Float64(line_width)
    local_color = local_color === nothing ? style.reference_color : local_color
    local_linewidth = local_linewidth === nothing ? max(1.8, 0.7 * style.fit_linewidth) : Float64(local_linewidth)
    threshold_color = threshold_color === nothing ? style.secondary_color : threshold_color
    fig = with_theme(_theme_from_style(resolved_style, resolved_appearance, theme_override)) do
        Figure(
            size=(Int(round(figure_size[1])), Int(round(figure_size[2]))),
            backgroundcolor=style.background_color,
        )
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
    hlines!(
        ax,
        [profile_result.threshold];
        _merged_kwargs((color=threshold_color, linewidth=max(1.8, 0.65 * line_width), linestyle=:dash), threshold_kwargs)...,
    )

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
    plot_contour(contour_result; filename=nothing, format=:pdf, theme=:sans, ...)

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
    theme::Symbol=:sans,
    appearance::Symbol=:auto,
    theme_override::Theme=Theme(),
    title="Contour",
    xlabel="parameter 1",
    ylabel="parameter 2",
    colormap=:viridis,
    show_heatmap::Bool=false,
    show_regions::Bool=true,
    show_legend::Bool=true,
    show_profile_lines::Bool=false,
    level_colors=nothing,
    region_colors=nothing,
    line_color=nothing,
    local_covariance=nothing,
    local_center=nothing,
    local_line_color=nothing,
    local_linewidth::Union{Nothing, Real}=nothing,
    local_linestyle=:dash,
    figure_size::Tuple{<:Real, <:Real}=(820, 700),
    axis_kwargs=NamedTuple(),
    heatmap_kwargs=NamedTuple(),
    contour_kwargs=NamedTuple(),
    local_contour_kwargs=NamedTuple(),
)
    resolved_style, resolved_appearance = _resolve_plot_style(theme, appearance)
    style = _style_preset(resolved_style, resolved_appearance)
    diagnostic_colors = _diagnostic_colors(resolved_style, resolved_appearance)
    level_colors = level_colors === nothing ? diagnostic_colors.levels : level_colors
    region_colors = region_colors === nothing ? diagnostic_colors.regions : region_colors
    local_line_color = local_line_color === nothing ? diagnostic_colors.local_color : local_line_color
    profile_linewidth = style.fit_linewidth
    local_linewidth = local_linewidth === nothing ? max(1.8, 0.7 * profile_linewidth) : Float64(local_linewidth)
    fig = with_theme(_theme_from_style(resolved_style, resolved_appearance, theme_override)) do
        Figure(
            size=(Int(round(figure_size[1])), Int(round(figure_size[2]))),
            backgroundcolor=style.background_color,
        )
    end
    # Left alignment keeps descriptive titles inside the scientific column when
    # a natural-width legend occupies the adjacent layout cell.
    ax = Axis(
        fig[1, 1];
        _merged_kwargs((title=title, titlealign=:left, xlabel=xlabel, ylabel=ylabel), axis_kwargs)...,
    )
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
    draw_profile_lines = show_profile_lines || show_heatmap || !show_regions
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

    contour_handles = Any[]
    contour_labels = String[]
    for (index, level) in enumerate(plot_levels)
        color = line_color === nothing ? colors[mod1(index, length(colors))] : line_color
        if draw_profile_lines
            contour!(
                ax,
                contour_result.x_values,
                contour_result.y_values,
                contour_result.delta_cost;
                _merged_kwargs((levels=[level], color=color, linewidth=profile_linewidth), contour_kwargs)...,
            )
            push!(contour_handles, LineElement(color=color, linewidth=profile_linewidth))
        else
            region_color = regions[mod1(index, length(regions))]
            push!(contour_handles, PolyElement(color=region_color, strokecolor=:transparent))
        end
        push!(contour_labels, _contour_level_label(level))
    end

    local_handle = nothing
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
        end
        local_handle = LineElement(
            color=local_line_color,
            linewidth=local_linewidth,
            linestyle=local_linestyle,
        )
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
        if local_handle !== nothing
            push!(handles, local_handle)
            push!(labels, _local_contour_label(plot_levels))
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

function _local_contour_delta(covariance::AbstractMatrix, center, xs, ys)
    cov = Matrix{Float64}(covariance)
    size(cov) == (2, 2) || throw(ArgumentError("local contour covariance must be 2x2"))
    precision = Symmetric(cov) \ Matrix{Float64}(I, 2, 2)
    out = Matrix{Float64}(undef, length(xs), length(ys))
    @inbounds for ix in eachindex(xs), iy in eachindex(ys)
        delta = [xs[ix] - center[1], ys[iy] - center[2]]
        out[ix, iy] = dot(delta, precision * delta)
    end
    return out
end

"""
    plot_profile_matrix(result; parameters=nothing, parameter_names=nothing, ...)
    plot_profile_matrix(matrix_result::ProfileMatrixResult; ...)

Create a kafe2-inspired profile/contour overview for several fitted
parameters.

Diagonal panels show one-parameter profile scans against the local covariance
parabola. Lower-triangle panels show two-parameter profile contours against the
local covariance ellipse. Upper-triangle panels show the local correlation
coefficient. The plot is intended as a fast diagnostic: if profile curves are
not parabolic, or profile contours do not resemble the local ellipse, symmetric
local covariance errors should not be treated as the final uncertainty
statement.

By default, panels with warnings or critical findings are marked. Set
`panel_status_mode=:issues`, `:all`, or `:none` independently of visual style.
Passing a precomputed `ProfileMatrixResult` renders the stored scans without
repeating the profile and contour refits.
"""
function plot_profile_matrix(
    result;
    parameters=nothing,
    parameter_names=nothing,
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:sans,
    appearance::Symbol=:auto,
    theme_override::Theme=Theme(),
    npoints_profile::Int=61,
    npoints_contour::Int=25,
    nsigma::Real=3,
    profile_threshold::Real=1.0,
    contour_levels::AbstractVector=[2.30, 6.18],
    adaptive::Bool=false,
    max_refinements::Int=2,
    max_points::Int=1200,
    panel_status_mode::Symbol=:issues,
    delta_max::Union{Nothing, Real}=nothing,
    figure_size=nothing,
)
    _validate_panel_status_mode(panel_status_mode)
    matrix = profile_matrix(
        result;
        parameters=parameters,
        parameter_names=parameter_names,
        npoints_profile=npoints_profile,
        npoints_contour=npoints_contour,
        nsigma=nsigma,
        profile_threshold=profile_threshold,
        contour_levels=contour_levels,
        adaptive=adaptive,
        max_refinements=max_refinements,
        max_points=max_points,
    )
    return plot_profile_matrix(
        matrix;
        filename=filename,
        format=format,
        theme=theme,
        appearance=appearance,
        theme_override=theme_override,
        panel_status_mode=panel_status_mode,
        delta_max=delta_max,
        figure_size=figure_size,
    )
end

function plot_profile_matrix(
    matrix::ProfileMatrixResult;
    parameter_names=nothing,
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:sans,
    appearance::Symbol=:auto,
    theme_override::Theme=Theme(),
    panel_status_mode::Symbol=:issues,
    delta_max::Union{Nothing, Real}=nothing,
    figure_size=nothing,
)
    selected = matrix.parameters
    names = parameter_names === nothing ? matrix.parameter_names : collect(parameter_names)
    n = length(selected)
    length(names) == n || throw(ArgumentError(
        "parameter_names must match the number of selected parameters",
    ))
    length(matrix.best_values) == n || throw(ArgumentError(
        "profile matrix best_values do not match selected parameters",
    ))
    length(matrix.local_stderr) == n || throw(ArgumentError(
        "profile matrix local_stderr do not match selected parameters",
    ))
    size(matrix.local_covariance) == (n, n) || throw(ArgumentError(
        "profile matrix local_covariance has incompatible dimensions",
    ))
    size(matrix.local_correlation) == (n, n) || throw(ArgumentError(
        "profile matrix local_correlation has incompatible dimensions",
    ))

    matrix_levels = Float64[]
    for cont in values(matrix.contours)
        append!(matrix_levels, cont.levels)
    end
    sort!(unique!(matrix_levels))

    resolved_style, resolved_appearance = _resolve_plot_style(theme, appearance)
    style = _style_preset(resolved_style, resolved_appearance)
    _validate_panel_status_mode(panel_status_mode)
    diagnostic_colors = _diagnostic_colors(resolved_style, resolved_appearance)
    cell = n <= 3 ? 285 : 235
    fig_size = figure_size === nothing ? (cell * n + 80, cell * n + 70) :
        (Int(round(figure_size[1])), Int(round(figure_size[2])))

    fig = with_theme(_theme_from_style(resolved_style, resolved_appearance, theme_override)) do
        Figure(size=fig_size, backgroundcolor=style.background_color)
    end

    profile_color = style.fit_color
    profile_linewidth = max(2.0, 0.8 * style.fit_linewidth)
    local_linewidth = max(1.6, 0.65 * style.fit_linewidth)
    local_color = diagnostic_colors.local_color
    region_colors = collect(diagnostic_colors.regions)
    isempty(region_colors) && (region_colors = [(style.fit_color, 0.20)])
    corr_color = style.stats_muted_color
    # Dense matrices scale the selected role down, but never below a readable
    # final-size floor. This keeps one typography contract across all plots.
    density_scale = style.diagnostic_scale * (n <= 3 ? 1.0 : 0.84)
    matrix_titlesize = max(n <= 3 ? 26 : 24, round(Int, density_scale * style.titlesize))
    matrix_labelsize = max(n <= 3 ? 23 : 21, round(Int, density_scale * style.xlabelsize))
    matrix_ticklabelsize = max(n <= 3 ? 20 : 18, round(Int, density_scale * style.ticklabelsize))
    matrix_legend_size = max(n <= 3 ? 20 : 19, round(Int, density_scale * style.legend_labelsize))
    matrix_status_size = max(17, matrix_ticklabelsize - 1)
    delta_label = resolved_style == :tex ? L"\Delta\mathrm{cost}" : "Δcost"

    for row in 1:n, col in 1:n
        ax = Axis(
            fig[row, col];
            xlabel=row == n ? names[col] : "",
            ylabel=col == 1 ? (row == col ? delta_label : names[row]) : "",
            title=row == 1 ? names[col] : "",
            titlealign=:center,
            titlesize=matrix_titlesize,
            titlegap=10,
            xlabelsize=matrix_labelsize,
            ylabelsize=matrix_labelsize,
            xticklabelsize=matrix_ticklabelsize,
            yticklabelsize=matrix_ticklabelsize,
        )

        if row == col
            index = selected[row]
            prof = matrix.profiles[index]
            lines!(ax, prof.values, prof.delta_cost; color=profile_color, linewidth=profile_linewidth)
            sigma = matrix.local_stderr[row]
            if isfinite(sigma) && sigma > 0
                local_delta = @. abs2((prof.values - matrix.best_values[row]) / sigma)
                lines!(ax, prof.values, local_delta; color=local_color, linewidth=local_linewidth, linestyle=:dash)
            end
            hlines!(ax, [prof.threshold]; color=style.stats_color, linestyle=:dot)
            ylimit = delta_max === nothing ?
                max(4.0 * prof.threshold, maximum(matrix_levels; init=0.0)) :
                Float64(delta_max)
            isfinite(ylimit) && ylimit > 0 && ylims!(ax, 0, ylimit)
            _draw_panel_status!(
                ax,
                matrix.panel_status[(index, index)],
                style,
                resolved_appearance;
                mode=panel_status_mode,
                fontsize=matrix_status_size,
            )
        elseif row > col
            xindex = selected[col]
            yindex = selected[row]
            cont = matrix.contours[(xindex, yindex)]
            contourf!(
                ax,
                cont.x_values,
                cont.y_values,
                cont.delta_cost;
                levels=vcat(0.0, sort(unique(cont.levels))),
                colormap=region_colors,
            )
            local_cov = matrix.local_covariance[[col, row], [col, row]]
            if all(isfinite, local_cov)
                local_delta = _local_contour_delta(
                    local_cov,
                    matrix.best_values[[col, row]],
                    cont.x_values,
                    cont.y_values,
                )
                for level in cont.levels
                    contour!(
                        ax,
                        cont.x_values,
                        cont.y_values,
                        local_delta;
                        levels=[level],
                        color=local_color,
                        linewidth=local_linewidth,
                        linestyle=:dash,
                    )
                end
            end
            scatter!(
                ax,
                [matrix.best_values[col]],
                [matrix.best_values[row]];
                marker=:cross,
                markersize=12,
                strokewidth=2.5,
                color=style.stats_color,
            )
            _draw_panel_status!(
                ax,
                matrix.panel_status[(xindex, yindex)],
                style,
                resolved_appearance;
                mode=panel_status_mode,
                fontsize=matrix_status_size,
            )
        else
            hidedecorations!(ax)
            hidespines!(ax)
            corr = matrix.local_correlation[row, col]
            label = isfinite(corr) ? "ρ = $(_fmt_value(corr; sigdigits=3))" : "ρ unavailable"
            text!(
                ax,
                0.5,
                0.5;
                text=label,
                space=:relative,
                align=(:center, :center),
                color=corr_color,
                fontsize=max(18, matrix_ticklabelsize),
            )
        end

        row < n && hidexdecorations!(ax; grid=false, label=false)
        col > 1 && row >= col && hideydecorations!(ax; grid=false, label=false)
    end

    # Hidden upper-triangle axes have little determinable content. Equal Auto
    # weights keep every diagnostic cell aligned without guessing pixel sizes.
    for index in 1:n
        colsize!(fig.layout, index, Auto(false, 1))
        rowsize!(fig.layout, index, Auto(false, 1))
    end

    handles = Any[LineElement(color=profile_color, linewidth=profile_linewidth)]
    labels = ["profile Δcost scan"]
    for (index, level) in enumerate(matrix_levels)
        region_name =
            isapprox(level, 2.30; atol=0.015) ? "1σ profile region" :
            isapprox(level, 6.18; atol=0.015) ? "2σ profile region" :
            "profile Δcost=$(_fmt_value(level; sigdigits=3)) region"
        push!(handles, PolyElement(color=region_colors[mod1(index, length(region_colors))], strokecolor=:transparent))
        push!(labels, region_name)
    end
    push!(handles, LineElement(color=local_color, linewidth=local_linewidth, linestyle=:dash))
    push!(labels, "local parabolic covariance")
    push!(handles, MarkerElement(marker=:cross, markersize=12, strokewidth=2.5, color=style.stats_color))
    push!(labels, "fit minimum")
    Legend(
        fig[n + 1, 1:n],
        handles,
        labels;
        framevisible=false,
        orientation=:horizontal,
        tellheight=true,
        nbanks=2,
        colgap=18,
        rowgap=4,
        labelsize=matrix_legend_size,
    )

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
        values = result.residuals
        errors = _yerror_for_plot(result.problem, result.params)
        _validate_diagnostic_plot_values(kind, x, values, errors)
        return x, values, errors, "Residuals", "y - fit", 0.0
    elseif kind == :pull
        values = _weighted_data_residual(result.problem, result.params)
        _validate_diagnostic_plot_values(kind, x, values, nothing)
        return x, values, nothing, "Pulls", "pull", 0.0
    elseif kind == :ratio
        all(isfinite, yhat) || throw(ArgumentError("ratio diagnostic requires finite model predictions"))
        all(!iszero, yhat) || throw(ArgumentError("ratio diagnostic is undefined when a model prediction is zero"))
        ratio = result.problem.y ./ yhat
        yerr = _yerror_for_plot(result.problem, result.params)
        ratio_err = yerr === nothing ? nothing : yerr ./ abs.(yhat)
        _validate_diagnostic_plot_values(kind, x, ratio, ratio_err)
        return x, ratio, ratio_err, "Ratio", "data / fit", 1.0
    end
    throw(ArgumentError("diagnostic plot kind must be :residual, :pull, or :ratio"))
end

function _validate_diagnostic_plot_values(kind::Symbol, x, values, errors)
    all(isfinite, x) || throw(ArgumentError("diagnostic plot x values must be finite"))
    all(isfinite, values) || throw(ArgumentError("$(kind) diagnostic values must be finite"))
    if errors !== nothing
        all(isfinite, errors) || throw(ArgumentError("$(kind) diagnostic errors must be finite"))
        all(>=(0.0), errors) || throw(ArgumentError("$(kind) diagnostic errors must be non-negative"))
    end
    return nothing
end

"""
    plot_residuals(result; kind=:pull, filename=nothing, format=:pdf, ...)

Plot residuals, pulls, or data/fit ratios for an XY fit. Marker shape, marker
size, and error-bar caps follow `theme` unless explicitly overridden.
"""
function plot_residuals(
    result::FitResult;
    kind::Symbol=:pull,
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:sans,
    appearance::Symbol=:auto,
    theme_override::Theme=Theme(),
    figure_size::Tuple{<:Real, <:Real}=(900, 520),
    xlabel="x",
    color=nothing,
    reference_color=nothing,
    marker=nothing,
    markersize::Union{Nothing, Real}=nothing,
    error_whiskerwidth::Union{Nothing, Real}=nothing,
    axis_kwargs=NamedTuple(),
    scatter_kwargs=NamedTuple(),
    errorbars_kwargs=NamedTuple(),
)
    x, values, errors, title, ylabel, reference = _diagnostic_values(result, kind)
    resolved_style, resolved_appearance = _resolve_plot_style(theme, appearance)
    style = _style_preset(resolved_style, resolved_appearance)
    color = color === nothing ? style.data_color : color
    reference_color = reference_color === nothing ? style.fit_color : reference_color
    marker = marker === nothing ? style.data_marker : marker
    markersize = markersize === nothing ? style.data_markersize : Float64(markersize)
    error_whiskerwidth = error_whiskerwidth === nothing ? style.error_whiskerwidth : Float64(error_whiskerwidth)
    fig = with_theme(_theme_from_style(resolved_style, resolved_appearance, theme_override)) do
        Figure(
            size=(Int(round(figure_size[1])), Int(round(figure_size[2]))),
            backgroundcolor=style.background_color,
        )
    end
    ax = Axis(fig[1, 1]; _merged_kwargs((title=title, xlabel=xlabel, ylabel=ylabel), axis_kwargs)...)
    hlines!(ax, [reference]; color=reference_color, linestyle=:dash)
    if errors !== nothing
        errorbars!(
            ax,
            x,
            values,
            errors;
            _merged_kwargs((color=color, whiskerwidth=error_whiskerwidth), errorbars_kwargs)...,
        )
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
`scatter_kwargs`, `errorbars_kwargs`, and `reference_line_kwargs` are applied
to every panel after the selected style defaults.
"""
function plot_diagnostics(
    result::FitResult;
    filename::Union{Nothing, AbstractString}=nothing,
    format::Symbol=:pdf,
    theme::Symbol=:sans,
    appearance::Symbol=:auto,
    theme_override::Theme=Theme(),
    figure_size::Tuple{<:Real, <:Real}=(900, 900),
    xlabel="x",
    color=nothing,
    reference_color=nothing,
    marker=nothing,
    markersize::Union{Nothing, Real}=nothing,
    error_whiskerwidth::Union{Nothing, Real}=nothing,
    axis_kwargs=NamedTuple(),
    scatter_kwargs=NamedTuple(),
    errorbars_kwargs=NamedTuple(),
    reference_line_kwargs=NamedTuple(),
)
    resolved_style, resolved_appearance = _resolve_plot_style(theme, appearance)
    style = _style_preset(resolved_style, resolved_appearance)
    color = color === nothing ? style.data_color : color
    reference_color = reference_color === nothing ? style.fit_color : reference_color
    marker = marker === nothing ? style.data_marker : marker
    markersize = markersize === nothing ? style.data_markersize : Float64(markersize)
    error_whiskerwidth = error_whiskerwidth === nothing ? style.error_whiskerwidth : Float64(error_whiskerwidth)
    fig = with_theme(_theme_from_style(resolved_style, resolved_appearance, theme_override)) do
        Figure(
            size=(Int(round(figure_size[1])), Int(round(figure_size[2]))),
            backgroundcolor=style.background_color,
        )
    end

    for (row, kind) in enumerate((:residual, :pull, :ratio))
        x, values, errors, title, ylabel, reference = _diagnostic_values(result, kind)
        ax = Axis(fig[row, 1]; _merged_kwargs((title=title, xlabel=row == 3 ? xlabel : "", ylabel=ylabel), axis_kwargs)...)
        hlines!(
            ax,
            [reference];
            _merged_kwargs((color=reference_color, linestyle=:dash), reference_line_kwargs)...,
        )
        if errors !== nothing
            errorbars!(
                ax,
                x,
                values,
                errors;
                _merged_kwargs((color=color, whiskerwidth=error_whiskerwidth), errorbars_kwargs)...,
            )
        end
        scatter!(
            ax,
            x,
            values;
            _merged_kwargs((color=color, marker=marker, markersize=markersize), scatter_kwargs)...,
        )
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
