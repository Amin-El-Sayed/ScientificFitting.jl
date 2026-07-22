const SNAPSHOT_ONLY = get(ENV, "JUFITTER_DOC_SNAPSHOT_ONLY", "0") == "1"
const RENDER_DOC_ASSETS = !SNAPSHOT_ONLY

if RENDER_DOC_ASSETS
    using CairoMakie
end

using Distributions
using JuFitter
using LaTeXStrings
using LinearAlgebra
using Printf
using SpecialFunctions

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output")
const DOC_ASSET_DIR = joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "gallery")
const EMIT_DOC_OUTPUT_SNAPSHOTS = get(ENV, "JUFITTER_DOC_OUTPUT_SNAPSHOTS", "0") == "1"

if RENDER_DOC_ASSETS
    mkpath(OUTPUT_DIR)
    mkpath(DOC_ASSET_DIR)
end

function gallery_path(name)
    return joinpath(DOC_ASSET_DIR, name)
end

function emit_doc_output_snapshot(body::Function, id::AbstractString)
    EMIT_DOC_OUTPUT_SNAPSHOTS || return nothing

    println("=== JUFITTER_DOC_OUTPUT_BEGIN ", id, " ===")
    body()
    println("=== JUFITTER_DOC_OUTPUT_END ", id, " ===")
    return nothing
end

function style_asset_suffix(dark, style::Symbol, appearance::Symbol)
    return dark === nothing ? "$(style)_$(appearance)" : (appearance == :dark ? "dark" : "light")
end

fmt_sig(x, digits::Integer=4) = @sprintf("%.*g", digits, x)
function fmt_tex(x, digits::Integer=4)
    value = fmt_sig(x, digits)
    scientific = match(r"^(.+)[eE]([+-]?\d+)$", value)
    scientific === nothing && return value
    mantissa, exponent = scientific.captures
    return mantissa * "\\times10^{" * string(parse(Int, exponent)) * "}"
end

function style_variant_plot(
    result,
    name;
    styles=(:lab, :modern, :article),
    plain=NamedTuple(),
    latex=plain,
    kwargs...,
)
    RENDER_DOC_ASSETS || return nothing

    for style in styles, appearance in (:light, :dark)
        typography = style == :article ? latex : plain
        plot_fit(
            result;
            kwargs...,
            typography...,
            theme=style,
            appearance=appearance,
            filename=gallery_path("$(name)_$(style)_$(appearance).png"),
            format=:png,
        )
    end
end

function poisson_deviance_residuals(counts, expected)
    return [
        sign(n - mu) * sqrt(max(2 * (n == 0 ? mu : mu - n + n * log(n / mu)), 0.0))
        for (n, mu) in zip(counts, expected)
    ]
end

function save_poisson_counts(
    result,
    x,
    counts,
    model,
    name;
    dark::Union{Nothing, Bool}=nothing,
    style::Symbol=:modern,
    appearance::Symbol=dark === nothing ? :light : (dark ? :dark : :light),
)
    RENDER_DOC_ASSETS || return nothing

    dark_mode = appearance == :dark
    palette = plot_palette(style; appearance=appearance)
    foreground = palette.stats_color
    muted = palette.stats_muted_color
    data_color = palette.data_color
    fit_color = palette.fit_color
    band_color = (palette.band_color, max(palette.band_alpha, 0.16))
    residual_positive = fit_color
    residual_negative = style == :article ? (dark_mode ? "#d7d7d7" : "#5f6873") :
                        (dark_mode ? "#b7c8dc" : "#52606f")
    fig = with_theme(plot_theme(style; appearance=appearance)) do
        Figure(size=(1460, 850), backgroundcolor=dark_mode ? "#111318" : "#ffffff")
    end
    ax = Axis(fig[1, 1]; title="Radioactive decay with detector background", ylabel="counts per 10 s")
    xg = collect(range(minimum(x), maximum(x); length=300))
    yg = model(xg, result.params)
    lower = [quantile(Poisson(mu), 0.16) for mu in yg]
    upper = [quantile(Poisson(mu), 0.84) for mu in yg]
    band!(ax, xg, lower, upper; color=band_color, label="central 68% count interval")
    lines!(ax, xg, yg; color=fit_color, linewidth=3, label="expected counts")
    scatter!(ax, x, counts; color=data_color, markersize=palette.data_markersize + 2.0, label="observed counts")
    hidexdecorations!(ax; grid=false)

    expected = model(x, result.params)
    residuals = poisson_deviance_residuals(counts, expected)
    residual_colors = ifelse.(residuals .>= 0, residual_positive, residual_negative)
    residual_ax = Axis(fig[2, 1]; xlabel="elapsed time (min)", ylabel="deviance residual")
    barplot!(residual_ax, x, residuals; width=0.66, color=residual_colors)
    hlines!(residual_ax, [0.0]; color=(foreground, 0.70), linewidth=1.5)
    hlines!(residual_ax, [-2.0, 2.0]; color=(foreground, 0.32), linestyle=:dash, linewidth=1.5)
    linkxaxes!(ax, residual_ax)

    half_life = log(2) / result.params[2]
    sigma_half_life = log(2) * result.param_stderr[2] / result.params[2]^2
    deviance_pvalue = ccdf(Chisq(result.stats.ndf), result.stats.chi2)
    article = style == :article
    parameter_lines = article ? Any[
        LaTeXString("S_0 = $(fmt_tex(result.params[1], 5)) \\pm $(fmt_tex(result.param_stderr[1], 2))\\;\\mathrm{counts}"),
        LaTeXString("\\lambda = $(fmt_tex(result.params[2], 5)) \\pm $(fmt_tex(result.param_stderr[2], 2))\\;\\mathrm{min^{-1}}"),
        LaTeXString("B = $(fmt_tex(result.params[3], 4)) \\pm $(fmt_tex(result.param_stderr[3], 2))\\;\\mathrm{counts}"),
        LaTeXString("t_{1/2} = $(fmt_tex(half_life, 4)) \\pm $(fmt_tex(sigma_half_life, 2))\\;\\mathrm{min}"),
    ] : Any[
        "S₀ = $(fmt_sig(result.params[1], 5)) ± $(fmt_sig(result.param_stderr[1], 2)) counts",
        "λ = $(fmt_sig(result.params[2], 5)) ± $(fmt_sig(result.param_stderr[2], 2)) min⁻¹",
        "B = $(fmt_sig(result.params[3], 4)) ± $(fmt_sig(result.param_stderr[3], 2)) counts",
        "t₁/₂ = $(fmt_sig(half_life, 4)) ± $(fmt_sig(sigma_half_life, 2)) min",
    ]
    statistic_lines = article ? Any[
        LaTeXString("D/\\mathrm{ndf} = $(fmt_tex(result.stats.chi2, 4))/$(result.stats.ndf) = $(fmt_tex(result.stats.chi2_ndf, 4))"),
        LaTeXString("P(D) = $(fmt_tex(deviance_pvalue, 4))"),
    ] : Any[
        "D/ndf = $(fmt_sig(result.stats.chi2, 4))/$(result.stats.ndf) = $(fmt_sig(result.stats.chi2_ndf, 4))",
        "P(D) = $(fmt_sig(deviance_pvalue, 4))",
    ]
    plot_info_panel!(
        fig[1:2, 2];
        legend_source=ax,
        title="Poisson likelihood fit",
        model_label=article ? L"\mu(t)=S_0 e^{-\lambda t}+B" : "μ(t) = S₀ exp(−λt) + B",
        parameter_lines=parameter_lines,
        statistic_lines=statistic_lines,
        color=foreground,
        muted_color=muted,
        fontsize=palette.stats_fontsize + 2,
    )
    rowsize!(fig.layout, 1, Relative(0.72))
    rowsize!(fig.layout, 2, Relative(0.28))
    colgap!(fig.layout, 24)
    save(gallery_path("$(name)_$(style_asset_suffix(dark, style, appearance)).png"), fig)
end

function save_histogram_fit(
    result,
    edges,
    counts,
    expected_counts,
    name;
    dark::Union{Nothing, Bool}=nothing,
    style::Symbol=:modern,
    appearance::Symbol=dark === nothing ? :light : (dark ? :dark : :light),
)
    RENDER_DOC_ASSETS || return nothing

    dark_mode = appearance == :dark
    palette = plot_palette(style; appearance=appearance)
    foreground = palette.stats_color
    muted = palette.stats_muted_color
    data_color = palette.data_color
    fit_color = palette.fit_color
    observed_color = style == :article ? (foreground, dark_mode ? 0.20 : 0.14) :
                     (palette.fit_color, dark_mode ? 0.20 : 0.22)
    background_base = style == :article ? foreground : (dark_mode ? "#b7c8dc" : "#52606f")
    background_color = (background_base, style == :article ? 0.10 : 0.20)
    residual_positive = fit_color
    residual_negative = style == :article ? (dark_mode ? "#d7d7d7" : "#5f6873") :
                        (dark_mode ? "#b7c8dc" : "#52606f")
    fig = with_theme(plot_theme(style; appearance=appearance)) do
        Figure(size=(1460, 850), backgroundcolor=dark_mode ? "#111318" : "#ffffff")
    end
    ax = Axis(fig[1, 1]; title="Binned detector spectrum", ylabel="events per bin")
    centers = [(edges[i] + edges[i + 1]) / 2 for i in 1:(length(edges) - 1)]
    widths = diff(edges)
    expected = expected_counts(edges, result.params)
    background = result.params[4] .* widths
    barplot!(ax, centers, background; width=widths, color=background_color, label="fitted background")
    barplot!(
        ax,
        centers,
        expected;
        width=widths,
        color=:transparent,
        strokecolor=fit_color,
        strokewidth=3,
        label="expected bin counts",
    )
    barplot!(
        ax,
        centers,
        counts;
        width=0.82 .* widths,
        color=observed_color,
        strokecolor=(foreground, 0.65),
        strokewidth=1.3,
        label="observed bin counts",
    )
    scatter!(ax, centers, counts; color=data_color, markersize=palette.data_markersize + 1.4)
    hidexdecorations!(ax; grid=false)

    residuals = poisson_deviance_residuals(counts, expected)
    residual_colors = ifelse.(residuals .>= 0, residual_positive, residual_negative)
    residual_ax = Axis(fig[2, 1]; xlabel="pulse amplitude (V)", ylabel="deviance residual")
    barplot!(residual_ax, centers, residuals; width=0.82 .* widths, color=residual_colors)
    hlines!(residual_ax, [0.0]; color=(foreground, 0.70), linewidth=1.5)
    hlines!(residual_ax, [-2.0, 2.0]; color=(foreground, 0.32), linestyle=:dash, linewidth=1.5)
    linkxaxes!(ax, residual_ax)

    deviance_pvalue = ccdf(Chisq(result.stats.ndf), result.stats.chi2)
    article = style == :article
    parameter_lines = article ? Any[
        LaTeXString("N_{\\mathrm{peak}} = $(fmt_tex(result.params[1], 5)) \\pm $(fmt_tex(result.param_stderr[1], 2))"),
        LaTeXString("\\mu = $(fmt_tex(result.params[2], 5)) \\pm $(fmt_tex(result.param_stderr[2], 2))\\;\\mathrm{V}"),
        LaTeXString("\\sigma = $(fmt_tex(result.params[3], 5)) \\pm $(fmt_tex(result.param_stderr[3], 2))\\;\\mathrm{V}"),
        LaTeXString("b = $(fmt_tex(result.params[4], 4)) \\pm $(fmt_tex(result.param_stderr[4], 2))\\;\\mathrm{events/V}"),
    ] : Any[
        "Nₚₑₐₖ = $(fmt_sig(result.params[1], 5)) ± $(fmt_sig(result.param_stderr[1], 2)) events",
        "μ = $(fmt_sig(result.params[2], 5)) ± $(fmt_sig(result.param_stderr[2], 2)) V",
        "σ = $(fmt_sig(result.params[3], 5)) ± $(fmt_sig(result.param_stderr[3], 2)) V",
        "b = $(fmt_sig(result.params[4], 4)) ± $(fmt_sig(result.param_stderr[4], 2)) events/V",
    ]
    statistic_lines = article ? Any[
        LaTeXString("D/\\mathrm{ndf} = $(fmt_tex(result.stats.chi2, 4))/$(result.stats.ndf) = $(fmt_tex(result.stats.chi2_ndf, 4))"),
        LaTeXString("P(D) = $(fmt_tex(deviance_pvalue, 4))"),
    ] : Any[
        "D/ndf = $(fmt_sig(result.stats.chi2, 4))/$(result.stats.ndf) = $(fmt_sig(result.stats.chi2_ndf, 4))",
        "P(D) = $(fmt_sig(deviance_pvalue, 4))",
    ]
    plot_info_panel!(
        fig[1:2, 2];
        legend_source=ax,
        title="Binned Poisson likelihood",
        model_label=article ? L"n_i \sim \mathrm{Poisson}(\mu_i)" : "nᵢ ~ Poisson(μᵢ)",
        parameter_lines=parameter_lines,
        statistic_lines=statistic_lines,
        color=foreground,
        muted_color=muted,
        fontsize=palette.stats_fontsize + 2,
    )
    rowsize!(fig.layout, 1, Relative(0.72))
    rowsize!(fig.layout, 2, Relative(0.28))
    colgap!(fig.layout, 24)
    save(gallery_path("$(name)_$(style_asset_suffix(dark, style, appearance)).png"), fig)
end

function line_intersection(emission_result, baseline_result)
    emission_slope, emission_intercept = emission_result.params
    baseline_slope, baseline_intercept = baseline_result.params
    denominator = emission_slope - baseline_slope
    threshold = (baseline_intercept - emission_intercept) / denominator

    threshold_gradient_emission = [-threshold / denominator, -1 / denominator]
    threshold_gradient_baseline = [threshold / denominator, 1 / denominator]
    threshold_variance =
        dot(threshold_gradient_emission, emission_result.param_covariance * threshold_gradient_emission) +
        dot(threshold_gradient_baseline, baseline_result.param_covariance * threshold_gradient_baseline)

    work_function_eV = emission_slope * threshold
    work_gradient_emission = [
        threshold + emission_slope * threshold_gradient_emission[1],
        emission_slope * threshold_gradient_emission[2],
    ]
    work_gradient_baseline = emission_slope .* threshold_gradient_baseline
    work_variance =
        dot(work_gradient_emission, emission_result.param_covariance * work_gradient_emission) +
        dot(work_gradient_baseline, baseline_result.param_covariance * work_gradient_baseline)

    return (
        threshold=threshold,
        sigma_threshold=sqrt(max(threshold_variance, 0.0)),
        work_function_eV=work_function_eV,
        sigma_work_function_eV=sqrt(max(work_variance, 0.0)),
    )
end

function save_photoelectric_work_function(
    emission_result,
    baseline_result,
    frequency,
    voltage,
    sigma_frequency,
    sigma_voltage,
    emission_mask,
    name;
    dark::Union{Nothing, Bool}=nothing,
    style::Symbol=:modern,
    appearance::Symbol=dark === nothing ? :light : (dark ? :dark : :light),
)
    RENDER_DOC_ASSETS || return nothing

    dark_mode = appearance == :dark
    palette = plot_palette(style; appearance=appearance)
    color = palette.data_color
    muted = palette.stats_muted_color
    emission_color = palette.fit_color
    baseline_color = style == :article ? (dark_mode ? "#d7d7d7" : "#555555") :
        (dark_mode ? "#b7c8dc" : "#52606f")
    threshold_color = style == :article ? (dark_mode ? "#edf2f4" : "#111111") :
        (dark_mode ? "#cbd5e1" : "#374151")
    emission_band = (palette.band_color, max(palette.band_alpha, 0.16))
    baseline_band = (baseline_color, style == :article ? 0.09 : 0.18)
    error_whiskerwidth = palette.error_whiskerwidth
    fit_linewidth = palette.fit_linewidth
    article = style == :article

    derived = line_intersection(emission_result, baseline_result)
    threshold = derived.threshold
    sigma_threshold = derived.sigma_threshold
    work_function_eV = derived.work_function_eV
    sigma_work_function_eV = derived.sigma_work_function_eV
    emission_slope, emission_intercept = emission_result.params
    baseline_slope, baseline_intercept = baseline_result.params

    fig = with_theme(plot_theme(style; appearance=appearance)) do
        Figure(size=(1360, 800), backgroundcolor=dark_mode ? "#111318" : "#ffffff")
    end
    ax = Axis(
        fig[1, 1];
        title="Photoelectric work-function extraction",
        xlabel=article ? L"\nu\;(\mathrm{THz})" : "frequency ν (THz)",
        ylabel=article ? L"U_0\;(\mathrm{V})" : "stopping voltage U₀ (V)",
    )

    errorbars!(ax, frequency, voltage, sigma_voltage; color=(muted, 0.44), whiskerwidth=error_whiskerwidth)
    errorbars!(ax, frequency, voltage, sigma_frequency; direction=:x, color=(muted, 0.30), whiskerwidth=error_whiskerwidth)
    scatter!(
        ax,
        frequency[emission_mask],
        voltage[emission_mask];
        color=color,
        markersize=palette.data_markersize,
        label=article ? L"\mathrm{emission\ regime}" : "emission regime",
    )
    scatter!(
        ax,
        frequency[.!emission_mask],
        voltage[.!emission_mask];
        color=baseline_color,
        marker=:diamond,
        markersize=palette.data_markersize,
        label=article ? L"\mathrm{baseline\ regime}" : "baseline regime",
    )

    xmin, xmax = extrema(frequency)
    xg = collect(range(xmin - 15, xmax + 15; length=500))
    emission_y = @. emission_slope * xg + emission_intercept
    baseline_y = @. baseline_slope * xg + baseline_intercept
    J = hcat(xg, ones(length(xg)))
    emission_sigma = sqrt.(clamp.(vec(sum((J * emission_result.param_covariance) .* J; dims=2)), 0.0, Inf))
    baseline_sigma = sqrt.(clamp.(vec(sum((J * baseline_result.param_covariance) .* J; dims=2)), 0.0, Inf))
    band!(
        ax,
        xg,
        emission_y .- emission_sigma,
        emission_y .+ emission_sigma;
        color=emission_band,
        label=article ? L"\mathrm{emission\ }1\sigma\mathrm{\ fit\ band}" : "emission 1σ fit band",
    )
    band!(
        ax,
        xg,
        baseline_y .- baseline_sigma,
        baseline_y .+ baseline_sigma;
        color=baseline_band,
        label=article ? L"\mathrm{baseline\ }1\sigma\mathrm{\ fit\ band}" : "baseline 1σ fit band",
    )
    lines!(ax, xg, emission_y; color=emission_color, linewidth=fit_linewidth, label=article ? L"\mathrm{emission\ fit}" : "emission fit")
    lines!(ax, xg, baseline_y; color=baseline_color, linewidth=fit_linewidth, label=article ? L"\mathrm{baseline\ fit}" : "baseline fit")

    vspan!(
        ax,
        threshold - sigma_threshold,
        threshold + sigma_threshold;
        color=(threshold_color, 0.14),
        label=article ? L"\mathrm{intersection\ }1\sigma\mathrm{\ interval}" : "intersection 1σ interval",
    )
    vlines!(ax, [threshold]; color=(threshold_color, 0.85), linestyle=:dash, linewidth=max(1.8, fit_linewidth - 0.2))
    threshold_y = emission_slope * threshold + emission_intercept
    scatter!(
        ax,
        [threshold],
        [threshold_y];
        color=threshold_color,
        marker=:star5,
        markersize=style == :article ? 13 : 16,
        label=article ? L"\mathrm{line\ intersection}" : "line intersection",
    )
    limits!(ax, xmin - 20, xmax + 20, minimum(voltage .- sigma_voltage) - 0.16, maximum(voltage .+ sigma_voltage) + 0.25)

    h_fit = emission_slope * 1.602176634e-19 / 1e12
    sigma_h = sqrt(max(emission_result.param_covariance[1, 1], 0.0)) * 1.602176634e-19 / 1e12
    parameter_lines = if article
        [
            LaTeXString("m_{\\mathrm{emit}} = $(fmt_tex(emission_slope, 5))\\;\\mathrm{V/THz}"),
            LaTeXString("h = $(fmt_tex(h_fit, 4)) \\pm $(fmt_tex(sigma_h, 2))\\;\\mathrm{J\\,s}"),
            LaTeXString("m_{\\mathrm{base}} = $(fmt_tex(baseline_slope, 4))\\;\\mathrm{V/THz}"),
            LaTeXString("\\nu_0 = $(fmt_tex(threshold, 5)) \\pm $(fmt_tex(sigma_threshold, 2))\\;\\mathrm{THz}"),
            LaTeXString("\\Phi = $(fmt_tex(work_function_eV, 5)) \\pm $(fmt_tex(sigma_work_function_eV, 2))\\;\\mathrm{eV}"),
        ]
    else
        [
            "mₑ = $(fmt_sig(emission_slope, 5)) V/THz",
            "h = $(fmt_sig(h_fit, 4)) ± $(fmt_sig(sigma_h, 2)) J s",
            "mᵦ = $(fmt_sig(baseline_slope, 4)) V/THz",
            "ν₀ = $(fmt_sig(threshold, 5)) ± $(fmt_sig(sigma_threshold, 2)) THz",
            "Φ = $(fmt_sig(work_function_eV, 5)) ± $(fmt_sig(sigma_work_function_eV, 2)) eV",
        ]
    end
    statistic_lines = if article
        [
            LaTeXString("\\chi^2_{\\mathrm{emit}}/\\mathrm{ndf} = $(fmt_tex(emission_result.stats.chi2_ndf, 4))"),
            LaTeXString("\\chi^2_{\\mathrm{base}}/\\mathrm{ndf} = $(fmt_tex(baseline_result.stats.chi2_ndf, 4))"),
        ]
    else
        [
            "χ²ₑ/ndf = $(fmt_sig(emission_result.stats.chi2_ndf, 4))",
            "χ²ᵦ/ndf = $(fmt_sig(baseline_result.stats.chi2_ndf, 4))",
        ]
    end
    plot_info_panel!(
        fig[1, 2];
        legend_source=ax,
        model_label=article ? L"U_0(\nu)=m_r\nu+b_r" : "U₀(ν) = mᵣν + bᵣ",
        parameter_lines=parameter_lines,
        statistic_lines=statistic_lines,
        color=palette.stats_color,
        fontsize=palette.stats_fontsize + 5,
        muted_color=muted,
    )
    colgap!(fig.layout, 18)
    suffix = dark === nothing ? "$(style)_$(appearance)" : (appearance == :dark ? "dark" : "light")
    save(gallery_path("$(name)_$(suffix).png"), fig)
end

# 0. Quickstart plot matching docs/src/quickstart.md.
quick_x = [0.0, 0.4348, 0.8696, 1.3043, 1.7391, 2.1739, 2.6087, 3.0435,
           3.4783, 3.9130, 4.3478, 4.7826, 5.2174, 5.6522, 6.0870, 6.5217,
           6.9565, 7.3913, 7.8261, 8.2609, 8.6957, 9.1304, 9.5652, 10.0]
quick_y = [0.7000, 1.6125, 2.4832, 3.2749, 3.9858, 4.6545, 5.3439, 6.1123,
           6.9838, 7.9338, 8.8975, 9.7983, 10.5849, 11.2581, 11.8738,
           12.5193, 13.2732, 14.1663, 15.1641, 16.1797, 17.1125, 17.8965,
           18.5336, 19.0964]
quick_sigma_y = [0.1600, 0.1687, 0.1774, 0.1861, 0.1948, 0.2035, 0.2122,
                 0.2209, 0.2296, 0.2383, 0.2470, 0.2557, 0.2643, 0.2730,
                 0.2817, 0.2904, 0.2991, 0.3078, 0.3165, 0.3252, 0.3339,
                 0.3426, 0.3513, 0.3600]
quick_model(x, p) = @. p[1] * x + p[2]
quick_result = fit_model(quick_model, quick_x, quick_y; p0=[1.0, 0.0], sigma_y=quick_sigma_y)
emit_doc_output_snapshot("quickstart") do
    println(report_text(quick_result; parameter_names=["m", "b"]))
    println()
    println(diagnostic_dashboard_text(quick_result))
end
style_variant_plot(
    quick_result,
    "quickstart_linear";
    plain=(
        title="Quickstart calibration",
        model_label="U(x) = m x + b",
        xlabel="x",
        xunit="mm",
        ylabel="U",
        yunit="V",
        parameter_names=["m", "b"],
        latex_labels=false,
        latex_stats=false,
        band_label="1σ prediction band",
    ),
    latex=(
        title=L"\mathrm{Quickstart\ calibration}",
        model_label=L"U(x)=m x+b",
        xlabel=L"x",
        xunit=L"\mathrm{mm}",
        ylabel=L"U",
        yunit=L"\mathrm{V}",
        parameter_names=[L"m", L"b"],
        latex_labels=true,
        latex_stats=true,
        band_label=L"1\sigma\ \mathrm{prediction\ band}",
    ),
    band=:prediction,
    nsigma=1,
    show_legend=true,
    stats_position=:right,
    stats_mode=:full,
    stats_fontsize=20,
    figure_size=(1200, 760),
)

# The same scientific content rendered through every public style preset.
if RENDER_DOC_ASSETS
    for style in (:lab, :modern, :article)
        typography = if style == :article
            (
                title=L"\mathrm{Quickstart\ calibration}",
                model_label=L"U(x)=m x+b",
                xlabel=L"x",
                xunit=L"\mathrm{mm}",
                ylabel=L"U",
                yunit=L"\mathrm{V}",
                parameter_names=[L"m", L"b"],
                latex_labels=true,
                latex_stats=true,
            )
        else
            (
                title="Quickstart calibration",
                model_label="U(x) = m x + b",
                xlabel="x",
                xunit="mm",
                ylabel="U",
                yunit="V",
                parameter_names=["m", "b"],
                latex_labels=false,
                latex_stats=false,
            )
        end
        plot_fit(
            quick_result;
            typography...,
            theme=style,
            filename=gallery_path("plot_style_$(style).png"),
            format=:png,
            band=:prediction,
            nsigma=1,
            band_label=style == :article ? L"1\sigma\ \mathrm{prediction\ band}" : "1σ prediction band",
            show_legend=true,
            stats_position=:right,
            stats_mode=:full,
            figure_size=(1200, 760),
        )
    end
end

# 1. Linear calibration with visible heteroscedastic uncertainties.
x = collect(range(0.0, 10.0; length=28))
scatter_scale = @. 0.22 + 0.025 * x
sigma_y = @. 0.10 + 0.012 * x
y = @. 0.82 + 1.72 * x + scatter_scale * (0.55 * sin(1.25 * x) + 0.18 * cos(3.7 * x))
calibration_model(x, p) = @. p[1] * x + p[2]
linear_result = fit_model(calibration_model, x, y; p0=[1.5, 0.5], sigma_y=sigma_y)
emit_doc_output_snapshot("linear_calibration") do
    println(report_text(linear_result; parameter_names=["m", "b"]))
    println(diagnostic_dashboard_text(linear_result))
end
style_variant_plot(
    linear_result,
    "linear_calibration";
    plain=(
        title="Sensor calibration",
        model_label="U(x) = m x + b",
        xlabel="x",
        xunit="mm",
        ylabel="U",
        yunit="V",
        parameter_names=["m", "b"],
        latex_labels=false,
        latex_stats=false,
        band_label="1σ prediction band",
    ),
    latex=(
        title=L"\mathrm{Sensor\ calibration}",
        model_label=L"U(x)=m x + b",
        xlabel=L"x",
        xunit=L"\mathrm{mm}",
        ylabel=L"U",
        yunit=L"\mathrm{V}",
        parameter_names=[L"m", L"b"],
        latex_labels=true,
        latex_stats=true,
        band_label=L"1\sigma\ \mathrm{prediction\ band}",
    ),
    band=:prediction,
    nsigma=1,
    show_legend=true,
    legend_position=:lt,
    stats_position=:right,
    stats_mode=:full,
    stats_fontsize=20,
    figure_size=(1200, 760),
)

# 2. Photoelectric work-function extraction from the intersection of two regimes.
frequency_THz = [350.0, 380.0, 410.0, 440.0, 470.0, 495.0, 515.0, 532.0,
                 565.0, 590.0, 620.0, 655.0, 690.0, 730.0, 775.0, 825.0, 880.0, 940.0]
sigma_frequency_THz = [4.5, 4.2, 4.0, 3.8, 3.6, 3.4, 3.2, 3.0,
                       2.9, 2.8, 2.7, 2.6, 2.5, 2.5, 2.4, 2.4, 2.3, 2.3]
sigma_voltage = [0.038, 0.040, 0.041, 0.043, 0.045, 0.047, 0.050, 0.052,
                 0.048, 0.050, 0.052, 0.054, 0.057, 0.060, 0.064, 0.068, 0.073, 0.080]
baseline_mask = frequency_THz .<= 532.0
emission_mask = .!baseline_mask
photo_model(f, p) = @. p[1] * f + p[2]
baseline_truth = [0.00012, -0.045]
emission_truth = [0.00413, -2.247]
baseline_offsets = [0.9, -1.1, 0.35, 1.2, -0.85, 0.65, -0.45, 0.75]
emission_offsets = [0.8, -1.2, 0.45, 1.0, -0.75, 1.15, -0.55, 0.7, -1.0, 0.35]
voltage = similar(frequency_THz)
voltage[baseline_mask] = photo_model(frequency_THz[baseline_mask], baseline_truth) .+
                         sigma_voltage[baseline_mask] .* baseline_offsets
voltage[emission_mask] = photo_model(frequency_THz[emission_mask], emission_truth) .+
                         sigma_voltage[emission_mask] .* emission_offsets
baseline_result = fit_model(
    photo_model,
    frequency_THz[baseline_mask],
    voltage[baseline_mask];
    p0=[0.0, 0.0],
    sigma_y=sigma_voltage[baseline_mask],
    sigma_x=sigma_frequency_THz[baseline_mask],
)
emission_result = fit_model(
    photo_model,
    frequency_THz[emission_mask],
    voltage[emission_mask];
    p0=[0.004, -2.2],
    sigma_y=sigma_voltage[emission_mask],
    sigma_x=sigma_frequency_THz[emission_mask],
    bounds=([0.0, -20.0], [0.02, 5.0]),
    initial_guesses=[[0.004, -2.2], [0.0042, -2.6], [0.0038, -1.8]],
)
emit_doc_output_snapshot("photoelectric_threshold") do
    derived = line_intersection(emission_result, baseline_result)
    elementary_charge = 1.602176634e-19
    emission_slope = emission_result.params[1]
    h_fit = emission_slope * elementary_charge / 1e12
    sigma_h = sqrt(emission_result.param_covariance[1, 1]) * elementary_charge / 1e12

    println("h = ", h_fit, " +/- ", sigma_h, " J s")
    println("Phi = ", derived.work_function_eV, " +/- ", derived.sigma_work_function_eV, " eV")
    println("nu0 = ", derived.threshold, " +/- ", derived.sigma_threshold, " THz")
    println()
    println("baseline")
    println(diagnostic_dashboard_text(baseline_result))
    println("emission")
    println(diagnostic_dashboard_text(emission_result))
end
for style in (:lab, :modern, :article), appearance in (:light, :dark)
    save_photoelectric_work_function(
        emission_result,
        baseline_result,
        frequency_THz,
        voltage,
        sigma_frequency_THz,
        sigma_voltage,
        emission_mask,
        "photoelectric_threshold";
        style=style,
        appearance=appearance,
    )
end

# 3. Exponential decay with full covariance.
x_cov = [0.0, 0.1190, 0.2381, 0.3571, 0.4762, 0.5952, 0.7143, 0.8333,
         0.9524, 1.0714, 1.1905, 1.3095, 1.4286, 1.5476, 1.6667, 1.7857,
         1.9048, 2.0238, 2.1429, 2.2619, 2.3810, 2.5]
y_cov = [2.25155, 2.00932, 1.79547, 1.60660, 1.44018, 1.29422, 1.16695,
         1.05651, 0.96079, 0.87740, 0.80381, 0.73758, 0.67658, 0.61929,
         0.56494, 0.51354, 0.46572, 0.42253, 0.38510, 0.35431, 0.33052,
         0.31346]
decay_model(x, p) = @. p[1] * exp(p[2] * x) + p[3]
n = length(x_cov)
base_sigma = 0.055
corr_len = 2.3
cov_y = [base_sigma^2 * exp(-abs(i - j) / corr_len) for i in 1:n, j in 1:n]
cov_result = fit_model(decay_model, x_cov, y_cov; p0=[1.5, -0.7, 0.0], cov_y=cov_y)
emit_doc_output_snapshot("full_covariance") do
    amplitude, decay_rate, offset = cov_result.params
    sigma_amplitude, sigma_decay_rate, sigma_offset = cov_result.param_stderr

    println("A = ", amplitude, " +/- ", sigma_amplitude)
    println("lambda = ", -decay_rate, " +/- ", sigma_decay_rate)
    println("C = ", offset, " +/- ", sigma_offset)
    println(diagnostic_dashboard_text(cov_result))
end
style_variant_plot(
    cov_result,
    "full_covariance_decay";
    plain=(
        title="Correlated decay data",
        model_label="y(t) = A exp(-λ t) + C",
        xlabel="t",
        xunit="s",
        ylabel="signal",
        parameter_names=["A", "λ", "C"],
        latex_labels=false,
        latex_stats=false,
        band_label="1σ prediction band",
    ),
    latex=(
        title=L"\mathrm{Correlated\ decay\ data}",
        model_label=L"y(t)=A e^{-\lambda t}+C",
        xlabel=L"t",
        xunit=L"\mathrm{s}",
        ylabel=L"signal",
        parameter_names=[L"A", L"\lambda", L"C"],
        latex_labels=true,
        latex_stats=true,
        band_label=L"1\sigma\ \mathrm{prediction\ band}",
    ),
    band=:prediction,
    nsigma=1,
    show_legend=true,
    legend_position=:lt,
    stats_position=:right,
    stats_mode=:full,
    stats_fontsize=20,
    figure_size=(1200, 760),
)

# 4. Effective-variance fit with x and y uncertainties.
line_model(x, p) = @. p[1] * x + p[2]
x_obs = [0.2240, 0.3698, 0.6835, 0.8413, 0.9811, 1.2948,
         1.4586, 1.7364, 1.8641, 2.0579, 2.3356, 2.5954,
         2.7172, 2.9949, 3.1047, 3.3885, 3.6362, 3.7640]
y_obs = [1.450, 1.694, 1.838, 1.978, 2.218, 2.366,
         2.502, 2.746, 2.946, 3.066, 3.274, 3.490,
         3.618, 3.774, 4.014, 4.218, 4.322, 4.530]
sigma_x = fill(0.050, length(x_obs))
sigma_y_xy = fill(0.033, length(x_obs))
xy_result = fit_model(
    line_model,
    x_obs,
    y_obs;
    p0=[0.8, 1.2],
    sigma_y=sigma_y_xy,
    sigma_x=sigma_x,
)
emit_doc_output_snapshot("xy_uncertainties") do
    println(report_text(xy_result; parameter_names=["m", "b"]))
    println(diagnostic_dashboard_text(xy_result))
end
style_variant_plot(
    xy_result,
    "xy_uncertainties";
    plain=(
        title="Calibration with x and y uncertainty",
        model_label="U(x) = m x + b",
        xlabel="measured position",
        xunit="mm",
        ylabel="measured voltage",
        yunit="V",
        parameter_names=["m", "b"],
        latex_labels=false,
        latex_stats=false,
        band_label="1σ prediction band",
    ),
    latex=(
        title=L"\mathrm{Calibration\ with\ x\ and\ y\ uncertainty}",
        model_label=L"U(x)=m x+b",
        xlabel=L"x_\mathrm{meas}",
        xunit=L"\mathrm{mm}",
        ylabel=L"U_\mathrm{meas}",
        yunit=L"\mathrm{V}",
        parameter_names=[L"m", L"b"],
        latex_labels=true,
        latex_stats=true,
        band_label=L"1\sigma\ \mathrm{prediction\ band}",
    ),
    band=:prediction,
    nsigma=1,
    show_legend=true,
    legend_position=:lt,
    stats_position=:right,
    stats_mode=:full,
    stats_fontsize=20,
    figure_size=(1200, 760),
)

# 5. Bounds, prior, profile, and a genuinely non-elliptic contour.
time_saturation = collect(range(0.15, 2.2; length=18))
saturation_model(t, p) = @. p[1] * (1 - exp(-t / p[2])) + p[3]
sigma_time = @. 0.010 + 0.004 * time_saturation
sigma_response = @. 0.045 + 0.008 * time_saturation
saturation_residual_pattern = [
    0.50, -0.90, 0.30, 1.10, -0.70, 0.80, -1.00, 0.40, 0.90,
    -0.60, 0.70, -0.80, 1.00, -0.40, 0.55, -0.75, 0.65, -0.35,
]
response_saturation =
    saturation_model(time_saturation, [4.8, 3.4, 0.12]) .+
    sigma_response .* saturation_residual_pattern
saturation_result = fit_model(
    saturation_model,
    time_saturation,
    response_saturation;
    p0=[3.0, 2.0, 0.0],
    sigma_y=sigma_response,
    sigma_x=sigma_time,
    bounds=([0.1, 0.1, -0.5], [20.0, 20.0, 1.0]),
    parameter_priors=(index=3, mean=0.10, sigma=0.08),
    initial_guesses=[[3.0, 2.0, 0.0], [8.0, 7.0, 0.1], [2.0, 1.0, 0.2]],
)
style_variant_plot(
    saturation_result,
    "constraints_priors";
    plain=(
        title="Early saturation measurement",
        model_label="y(t) = A (1 - exp(-t / tau)) + c",
        xlabel="t",
        xunit="s",
        ylabel="y",
        yunit="V",
        parameter_names=["A", "tau", "c"],
        latex_labels=false,
        latex_stats=false,
        band_label="1σ prediction band",
    ),
    latex=(
        title=L"\mathrm{Early\ saturation\ measurement}",
        model_label=L"y(t)=A(1-e^{-t/\tau})+c",
        xlabel=L"t",
        xunit=L"\mathrm{s}",
        ylabel=L"y",
        yunit=L"\mathrm{V}",
        parameter_names=[L"A", L"\tau", L"c"],
        latex_labels=true,
        latex_stats=true,
        band_label=L"1\sigma\ \mathrm{prediction\ band}",
    ),
    band=:prediction,
    nsigma=1,
    show_legend=true,
    legend_position=:lt,
    stats_position=:right,
    stats_mode=:full,
    stats_fontsize=20,
    figure_size=(1200, 760),
)
amplitude_interval = profile_interval(saturation_result, 1; npoints=81, nsigma=4)
emit_doc_output_snapshot("constraints_profiles") do
    println("amplitude center = ", saturation_result.params[1])
    println("amplitude 1sigma lower = ", amplitude_interval.lower)
    println("amplitude 1sigma upper = ", amplitude_interval.upper)
    println("amplitude -sigma = ", amplitude_interval.uncertainty_minus)
    println("amplitude +sigma = ", amplitude_interval.uncertainty_plus)
    println(diagnostic_dashboard_text(saturation_result))
end
if RENDER_DOC_ASSETS
    prof = JuFitter.profile(saturation_result, 1; npoints=61, nsigma=4)
    cont = JuFitter.contour(saturation_result, 1, 2; npoints=121, nsigma=4)
    profile_overview_parameters = [1, 2, 3]
    profile_overview_names = ["A", "tau", "c"]

    for style in (:lab, :modern, :article), appearance in (:light, :dark)
        plot_profile(
            prof;
            theme=style,
            appearance=appearance,
            title=style == :article ? L"\mathrm{Profile\ cost\ versus\ local\ parabola}" :
                  "Profile cost versus local parabola",
            xlabel=style == :article ? L"\mathrm{amplitude}\ A" : "amplitude A",
            local_sigma=saturation_result.param_stderr[1],
            delta_max=8,
            filename=gallery_path("saturation_profile_$(style)_$(appearance).png"),
            format=:png,
        )
        plot_contour(
            cont;
            theme=style,
            appearance=appearance,
            title=style == :article ? L"\mathrm{Profile\ contours\ versus\ local\ covariance}" :
                  "Profile contours versus local covariance",
            xlabel=style == :article ? L"\mathrm{amplitude}\ A" : "amplitude A",
            ylabel=style == :article ? L"\mathrm{time\ constant}\ \tau" : "time constant tau",
            local_covariance=saturation_result.param_covariance,
            local_center=saturation_result.params[[1, 2]],
            figure_size=(980, 720),
            filename=gallery_path("amplitude_timescale_contour_$(style)_$(appearance).png"),
            format=:png,
        )
    end
    for style in (:lab, :modern, :article), appearance in (:light, :dark)
        plot_profile_matrix(
            saturation_result;
            parameters=profile_overview_parameters,
            parameter_names=profile_overview_names,
            npoints_profile=41,
            npoints_contour=21,
            nsigma=3,
            adaptive=true,
            max_refinements=1,
            panel_status_mode=:issues,
            theme=style,
            appearance=appearance,
            figure_size=(1020, 980),
            filename=gallery_path("saturation_profile_matrix_$(style)_$(appearance).png"),
            format=:png,
        )
    end
end

# 6. Poisson decay counts and a binned detector spectrum.
x_counts = collect(0.0:1.0:18.0)
counts = [48, 37, 35, 27, 27, 17, 22, 13, 16, 8, 13, 5, 11, 4, 7, 2, 6, 1, 5]
poisson_model(t, p) = @. p[1] * exp(-p[2] * t) + p[3]
poisson_result = fit_poisson_model(
    poisson_model,
    x_counts,
    counts;
    p0=[40.0, 0.15, 3.0],
    bounds=([1e-6, 1e-6, 1e-6], [200.0, 2.0, 50.0]),
    parameter_names=["initial signal", "decay constant", "background"],
    initial_guesses=[[40.0, 0.15, 3.0], [70.0, 0.30, 2.0], [25.0, 0.08, 5.0]],
)
emit_doc_output_snapshot("poisson_decay") do
    lambda = poisson_result.params[2]
    sigma_lambda = poisson_result.param_stderr[2]
    half_life = log(2) / lambda
    sigma_half_life = log(2) * sigma_lambda / lambda^2

    println("half-life = ", half_life, " +/- ", sigma_half_life, " min")
    println("deviance/ndf = ", poisson_result.stats.chi2_ndf)
    println("P(D) = ", poisson_result.stats.pvalue)
    println(diagnostic_dashboard_text(poisson_result))
end
for style in (:lab, :modern, :article), appearance in (:light, :dark)
    save_poisson_counts(
        poisson_result,
        x_counts,
        counts,
        poisson_model,
        "poisson_counts";
        style=style,
        appearance=appearance,
    )
end

edges = [0.0, 0.4, 0.9, 1.5, 2.2, 3.0, 4.0, 5.2, 6.6, 8.2, 10.0]
hist_counts = [0, 3, 9, 24, 47, 69, 51, 24, 8, 4]
function expected_counts(edges, p)
    peak_yield, centroid, width, background_density = p
    return [
        peak_yield * 0.5 * (
            erf((edges[i + 1] - centroid) / (sqrt(2) * width)) -
            erf((edges[i] - centroid) / (sqrt(2) * width))
        ) + background_density * (edges[i + 1] - edges[i])
        for i in 1:(length(edges) - 1)
    ]
end
hist_result = fit_histogram_model(
    expected_counts,
    edges,
    hist_counts;
    p0=[210.0, 3.8, 1.0, 1.0],
    bounds=([1e-6, 0.0, 0.05, 1e-6], [1000.0, 10.0, 5.0, 100.0]),
    parameter_names=["peak yield", "centroid", "width", "background density"],
    initial_guesses=[[210.0, 3.8, 1.0, 1.0], [300.0, 4.2, 1.5, 0.5], [150.0, 3.2, 0.7, 2.0]],
)
emit_doc_output_snapshot("histogram_likelihood") do
    println(
        "centroid = ",
        hist_result.params[2],
        " +/- ",
        hist_result.param_stderr[2],
        " V",
    )
    println("deviance/ndf = ", hist_result.stats.chi2_ndf)
    println("P(D) = ", hist_result.stats.pvalue)
    println(diagnostic_dashboard_text(hist_result))
end
for style in (:lab, :modern, :article), appearance in (:light, :dark)
    save_histogram_fit(
        hist_result,
        edges,
        hist_counts,
        expected_counts,
        "histogram_likelihood";
        style=style,
        appearance=appearance,
    )
end

# 7. Multi-dataset model criticism and partial parameter sharing.
include(joinpath(@__DIR__, "10_multi_dataset_calibration.jl"))

println("Generated documentation gallery assets in ", DOC_ASSET_DIR)
