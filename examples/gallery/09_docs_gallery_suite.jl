using CairoMakie
using Distributions
using JuFitter
using LaTeXStrings
using LinearAlgebra
using Printf
using SpecialFunctions

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output")
const DOC_ASSET_DIR = joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "gallery")

mkpath(OUTPUT_DIR)
mkpath(DOC_ASSET_DIR)

function gallery_path(name)
    return joinpath(DOC_ASSET_DIR, name)
end

fmt_sig(x, digits::Integer=4) = @sprintf("%.*g", digits, x)

function lightdark_plot(result, name; kwargs...)
    plot_fit(
        result;
        kwargs...,
        theme=:minimal,
        filename=gallery_path("$(name)_light.png"),
        format=:png,
    )
    plot_fit(
        result;
        kwargs...,
        theme=:dark,
        filename=gallery_path("$(name)_dark.png"),
        format=:png,
    )
end

function gallery_theme(dark::Bool)
    if dark
        return Theme(
            fontsize=20,
            font="TeX Gyre Heros",
            Figure=(backgroundcolor="#111318",),
            Axis=(
                backgroundcolor="#111318",
                xlabelsize=27,
                ylabelsize=27,
                titlesize=26,
                xticklabelsize=20,
                yticklabelsize=20,
                xlabelcolor="#edf2f4",
                ylabelcolor="#edf2f4",
                titlecolor="#edf2f4",
                xticklabelcolor="#b8c1ca",
                yticklabelcolor="#b8c1ca",
                xtickcolor="#b8c1ca",
                ytickcolor="#b8c1ca",
                xgridcolor=("#2a313a", 0.85),
                ygridcolor=("#2a313a", 0.85),
                topspinevisible=false,
                rightspinevisible=false,
                leftspinecolor="#b8c1ca",
                bottomspinecolor="#b8c1ca",
            ),
            Legend=(framevisible=false, labelcolor="#edf2f4", labelsize=20, patchsize=(30, 16)),
        )
    end
    return Theme(
        fontsize=20,
        font="TeX Gyre Heros",
        Figure=(backgroundcolor="#ffffff",),
        Axis=(
            backgroundcolor="#ffffff",
            xlabelsize=27,
            ylabelsize=27,
            titlesize=26,
            xticklabelsize=20,
            yticklabelsize=20,
            xgridcolor=("#eef2f7", 0.9),
            ygridcolor=("#eef2f7", 0.9),
            topspinevisible=false,
            rightspinevisible=false,
        ),
        Legend=(framevisible=false, labelsize=20, patchsize=(30, 16)),
    )
end

function poisson_deviance_residuals(counts, expected)
    return [
        sign(n - mu) * sqrt(max(2 * (n == 0 ? mu : mu - n + n * log(n / mu)), 0.0))
        for (n, mu) in zip(counts, expected)
    ]
end

function save_poisson_counts(result, x, counts, model, name; dark::Bool=false)
    color = dark ? "#edf2f4" : "#14151a"
    fit_color = dark ? "#66d9ef" : "#0081a7"
    band_color = dark ? ("#66d9ef", 0.16) : ("#a8dadc", 0.34)
    residual_positive = dark ? "#66d9ef" : "#0081a7"
    residual_negative = dark ? "#f4b860" : "#b45309"
    fig = with_theme(gallery_theme(dark)) do
        Figure(size=(1460, 850), backgroundcolor=dark ? "#111318" : "#ffffff")
    end
    ax = Axis(fig[1, 1]; title="Radioactive decay with detector background", ylabel="counts per 10 s")
    xg = collect(range(minimum(x), maximum(x); length=300))
    yg = model(xg, result.params)
    lower = [quantile(Poisson(mu), 0.16) for mu in yg]
    upper = [quantile(Poisson(mu), 0.84) for mu in yg]
    band!(ax, xg, lower, upper; color=band_color, label="central 68% count interval")
    lines!(ax, xg, yg; color=fit_color, linewidth=3, label="expected counts")
    scatter!(ax, x, counts; color=color, markersize=11, label="observed counts")
    hidexdecorations!(ax; grid=false)

    expected = model(x, result.params)
    residuals = poisson_deviance_residuals(counts, expected)
    residual_colors = ifelse.(residuals .>= 0, residual_positive, residual_negative)
    residual_ax = Axis(fig[2, 1]; xlabel="elapsed time (min)", ylabel="deviance residual")
    barplot!(residual_ax, x, residuals; width=0.66, color=residual_colors)
    hlines!(residual_ax, [0.0]; color=(color, 0.70), linewidth=1.5)
    hlines!(residual_ax, [-2.0, 2.0]; color=(color, 0.32), linestyle=:dash, linewidth=1.5)
    linkxaxes!(ax, residual_ax)

    side = GridLayout()
    fig[1:2, 2] = side
    Legend(side[1, 1], ax; framevisible=false, tellheight=true)
    panel = Axis(side[2, 1]; backgroundcolor=:transparent)
    hidedecorations!(panel)
    hidespines!(panel)
    half_life = log(2) / result.params[2]
    sigma_half_life = log(2) * result.param_stderr[2] / result.params[2]^2
    text!(
        panel,
        0,
        1;
        text="signal at t = 0\n" *
             "  $(fmt_sig(result.params[1], 5)) ± $(fmt_sig(result.param_stderr[1], 2)) counts\n\n" *
             "decay constant λ\n" *
             "  $(fmt_sig(result.params[2], 5)) ± $(fmt_sig(result.param_stderr[2], 2)) min⁻¹\n\n" *
             "background (local covariance)\n" *
             "  $(fmt_sig(result.params[3], 4)) ± $(fmt_sig(result.param_stderr[3], 2)) counts\n\n" *
             "half-life (local propagation)\n" *
             "  $(fmt_sig(half_life, 4)) ± $(fmt_sig(sigma_half_life, 2)) min\n\n" *
             "Poisson deviance / ndf\n" *
             "  $(fmt_sig(result.stats.chi2, 4)) / $(result.stats.ndf) = $(fmt_sig(result.stats.chi2_ndf, 4))\n" *
             "asymptotic P(D) = $(fmt_sig(result.stats.pvalue, 4))",
        space=:relative,
        align=(:left, :top),
        color=color,
        fontsize=20,
        lineheight=1.1,
    )
    rowsize!(side, 1, Auto())
    rowsize!(side, 2, Relative(1))
    rowsize!(fig.layout, 1, Relative(0.72))
    rowsize!(fig.layout, 2, Relative(0.28))
    colsize!(fig.layout, 2, Fixed(430))
    save(gallery_path("$(name)_$(dark ? "dark" : "light").png"), fig)
end

function save_histogram_fit(result, edges, counts, expected_counts, name; dark::Bool=false)
    color = dark ? "#edf2f4" : "#14151a"
    fit_color = dark ? "#66d9ef" : "#0081a7"
    observed_color = dark ? ("#edf2f4", 0.28) : ("#4c78a8", 0.28)
    background_color = dark ? ("#f4b860", 0.20) : ("#f4b183", 0.34)
    residual_positive = dark ? "#66d9ef" : "#0081a7"
    residual_negative = dark ? "#f4b860" : "#b45309"
    fig = with_theme(gallery_theme(dark)) do
        Figure(size=(1460, 850), backgroundcolor=dark ? "#111318" : "#ffffff")
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
        strokecolor=(color, 0.65),
        strokewidth=1.3,
        label="observed bin counts",
    )
    scatter!(ax, centers, counts; color=color, markersize=8)
    hidexdecorations!(ax; grid=false)

    residuals = poisson_deviance_residuals(counts, expected)
    residual_colors = ifelse.(residuals .>= 0, residual_positive, residual_negative)
    residual_ax = Axis(fig[2, 1]; xlabel="pulse amplitude (V)", ylabel="deviance residual")
    barplot!(residual_ax, centers, residuals; width=0.82 .* widths, color=residual_colors)
    hlines!(residual_ax, [0.0]; color=(color, 0.70), linewidth=1.5)
    hlines!(residual_ax, [-2.0, 2.0]; color=(color, 0.32), linestyle=:dash, linewidth=1.5)
    linkxaxes!(ax, residual_ax)

    side = GridLayout()
    fig[1:2, 2] = side
    Legend(side[1, 1], ax; framevisible=false, tellheight=true)
    panel = Axis(side[2, 1]; backgroundcolor=:transparent)
    hidedecorations!(panel)
    hidespines!(panel)
    text!(
        panel,
        0,
        1;
        text="peak yield\n" *
             "  $(fmt_sig(result.params[1], 5)) ± $(fmt_sig(result.param_stderr[1], 2)) events\n\n" *
             "centroid\n" *
             "  $(fmt_sig(result.params[2], 5)) ± $(fmt_sig(result.param_stderr[2], 2)) V\n\n" *
             "Gaussian width\n" *
             "  $(fmt_sig(result.params[3], 5)) ± $(fmt_sig(result.param_stderr[3], 2)) V\n\n" *
             "background density\n" *
             "  $(fmt_sig(result.params[4], 4)) ± $(fmt_sig(result.param_stderr[4], 2)) events/V\n\n" *
             "Poisson deviance / ndf\n" *
             "  $(fmt_sig(result.stats.chi2, 4)) / $(result.stats.ndf) = $(fmt_sig(result.stats.chi2_ndf, 4))\n" *
             "asymptotic P(D) = $(fmt_sig(result.stats.pvalue, 4))",
        space=:relative,
        align=(:left, :top),
        color=color,
        fontsize=20,
        lineheight=1.1,
    )
    rowsize!(side, 1, Auto())
    rowsize!(side, 2, Relative(1))
    rowsize!(fig.layout, 1, Relative(0.72))
    rowsize!(fig.layout, 2, Relative(0.28))
    colsize!(fig.layout, 2, Fixed(430))
    save(gallery_path("$(name)_$(dark ? "dark" : "light").png"), fig)
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
    dark::Bool=false,
)
    color = dark ? "#edf2f4" : "#14151a"
    muted = dark ? "#b8c1ca" : "#5b6270"
    emission_color = dark ? "#66d9ef" : "#007f9e"
    baseline_color = dark ? "#f4b860" : "#b85c38"
    threshold_color = dark ? "#f7e06e" : "#7a5c00"
    emission_band = dark ? ("#66d9ef", 0.20) : ("#89d5e0", 0.30)
    baseline_band = dark ? ("#f4b860", 0.18) : ("#f4b183", 0.30)

    derived = line_intersection(emission_result, baseline_result)
    threshold = derived.threshold
    sigma_threshold = derived.sigma_threshold
    work_function_eV = derived.work_function_eV
    sigma_work_function_eV = derived.sigma_work_function_eV
    emission_slope, emission_intercept = emission_result.params
    baseline_slope, baseline_intercept = baseline_result.params

    fig = with_theme(gallery_theme(dark)) do
        Figure(size=(1460, 800), backgroundcolor=dark ? "#111318" : "#ffffff")
    end
    ax = Axis(
        fig[1, 1];
        title="Photoelectric work-function extraction",
        xlabel="frequency ν (THz)",
        ylabel="stopping voltage U₀ (V)",
    )

    errorbars!(ax, frequency, voltage, sigma_voltage; color=(muted, 0.46), whiskerwidth=5)
    errorbars!(ax, frequency, voltage, sigma_frequency; direction=:x, color=(muted, 0.30), whiskerwidth=5)
    scatter!(ax, frequency[emission_mask], voltage[emission_mask]; color=color, markersize=10, label="emission regime")
    scatter!(ax, frequency[.!emission_mask], voltage[.!emission_mask]; color=baseline_color, marker=:diamond, markersize=10, label="baseline regime")

    xmin, xmax = extrema(frequency)
    xg = collect(range(xmin - 15, xmax + 15; length=500))
    emission_y = @. emission_slope * xg + emission_intercept
    baseline_y = @. baseline_slope * xg + baseline_intercept
    J = hcat(xg, ones(length(xg)))
    emission_sigma = sqrt.(clamp.(vec(sum((J * emission_result.param_covariance) .* J; dims=2)), 0.0, Inf))
    baseline_sigma = sqrt.(clamp.(vec(sum((J * baseline_result.param_covariance) .* J; dims=2)), 0.0, Inf))
    band!(ax, xg, emission_y .- emission_sigma, emission_y .+ emission_sigma; color=emission_band, label="emission 1σ fit band")
    band!(ax, xg, baseline_y .- baseline_sigma, baseline_y .+ baseline_sigma; color=baseline_band, label="baseline 1σ fit band")
    lines!(ax, xg, emission_y; color=emission_color, linewidth=3.2, label="emission fit")
    lines!(ax, xg, baseline_y; color=baseline_color, linewidth=3.2, label="baseline fit")

    vspan!(ax, threshold - sigma_threshold, threshold + sigma_threshold; color=(threshold_color, 0.14), label="intersection 1σ interval")
    vlines!(ax, [threshold]; color=(threshold_color, 0.85), linestyle=:dash, linewidth=2.5)
    threshold_y = emission_slope * threshold + emission_intercept
    scatter!(ax, [threshold], [threshold_y]; color=threshold_color, marker=:star5, markersize=18, label="line intersection")
    limits!(ax, xmin - 20, xmax + 20, minimum(voltage .- sigma_voltage) - 0.16, maximum(voltage .+ sigma_voltage) + 0.25)

    h_fit = emission_slope * 1.602176634e-19 / 1e12
    sigma_h = sqrt(max(emission_result.param_covariance[1, 1], 0.0)) * 1.602176634e-19 / 1e12
    side = GridLayout()
    fig[1, 2] = side
    Legend(side[1, 1], ax; framevisible=false, tellheight=true)
    panel = Axis(side[2, 1]; backgroundcolor=:transparent)
    hidedecorations!(panel)
    hidespines!(panel)
    text!(
        panel,
        0,
        1;
        text="emission slope h/e = $(fmt_sig(emission_slope, 5)) V/THz\n" *
             "h = $(fmt_sig(h_fit, 4)) ± $(fmt_sig(sigma_h, 2)) J s\n" *
             "baseline slope = $(fmt_sig(baseline_slope, 4)) V/THz\n" *
             "ν₀ = $(fmt_sig(threshold, 5)) ± $(fmt_sig(sigma_threshold, 2)) THz\n" *
             "Φ = $(fmt_sig(work_function_eV, 5)) ± $(fmt_sig(sigma_work_function_eV, 2)) eV\n\n" *
             "emission χ²/ndf = $(fmt_sig(emission_result.stats.chi2_ndf, 4))\n" *
             "baseline χ²/ndf = $(fmt_sig(baseline_result.stats.chi2_ndf, 4))",
        space=:relative,
        align=(:left, :top),
        color=color,
        fontsize=20,
        lineheight=1.12,
    )
    rowsize!(side, 1, Auto())
    rowsize!(side, 2, Relative(1))
    colsize!(fig.layout, 2, Fixed(480))
    save(gallery_path("$(name)_$(dark ? "dark" : "light").png"), fig)
end

function save_multi_dataset(result, x1, y1, x2, y2, sigma1, sigma2, name; dark::Bool=false)
    color = dark ? "#edf2f4" : "#14151a"
    fit_a = dark ? "#66d9ef" : "#0081a7"
    fit_b = dark ? "#f59e0b" : "#b45309"
    band_a = dark ? ("#66d9ef", 0.16) : ("#a8dadc", 0.28)
    band_b = dark ? ("#f59e0b", 0.14) : ("#fed7aa", 0.42)
    p = result.params
    cov = result.param_covariance
    fig = with_theme(gallery_theme(dark)) do
        Figure(size=(1260, 760), backgroundcolor=dark ? "#111318" : "#ffffff")
    end
    ax = Axis(fig[1, 1]; title="Shared-slope calibration transfer", xlabel="input x", ylabel="response y")
    errorbars!(ax, x1, y1, sigma1; color=(fit_a, 0.45), whiskerwidth=6)
    errorbars!(ax, x2, y2, sigma2; color=(fit_b, 0.45), whiskerwidth=6)
    scatter!(ax, x1, y1; color=fit_a, markersize=8, label="sensor A")
    scatter!(ax, x2, y2; color=fit_b, markersize=8, marker=:rect, label="sensor B")
    xg = collect(range(0, 5; length=250))
    y_a = @. p[1] * xg + p[2]
    y_b = @. p[1] * xg + p[3]
    Ja = hcat(xg, ones(length(xg)), zeros(length(xg)))
    Jb = hcat(xg, zeros(length(xg)), ones(length(xg)))
    sigma_a = sqrt.(clamp.(vec(sum((Ja * cov) .* Ja; dims=2)), 0.0, Inf))
    sigma_b = sqrt.(clamp.(vec(sum((Jb * cov) .* Jb; dims=2)), 0.0, Inf))
    band!(ax, xg, y_a .- sigma_a, y_a .+ sigma_a; color=band_a, label="1σ fit band A")
    band!(ax, xg, y_b .- sigma_b, y_b .+ sigma_b; color=band_b, label="1σ fit band B")
    lines!(ax, xg, y_a; color=fit_a, linewidth=3, label="fit A")
    lines!(ax, xg, y_b; color=fit_b, linewidth=3, label="fit B")
    axislegend(ax; position=:lt, nbanks=2)
    offset_gap = p[2] - p[3]
    offset_grad = [0.0, 1.0, -1.0]
    sigma_gap = sqrt(max(dot(offset_grad, cov * offset_grad), 0.0))
    panel = Axis(fig[1, 2]; backgroundcolor=:transparent)
    hidedecorations!(panel)
    hidespines!(panel)
    text!(
        panel,
        0,
        1;
        text="shared slope = $(round(p[1]; sigdigits=5))\n" *
             "offset A = $(round(p[2]; sigdigits=5))\n" *
             "offset B = $(round(p[3]; sigdigits=5))\n" *
             "offset gap = $(round(offset_gap; sigdigits=5)) ± $(round(sigma_gap; sigdigits=2))\n" *
             "χ²/ndf = $(round(result.stats.chi2_ndf; sigdigits=4))\n" *
             "P(χ²) = $(round(result.stats.pvalue; sigdigits=4))",
        space=:relative,
        align=(:left, :top),
        color=color,
        fontsize=20,
        lineheight=1.15,
    )
    colsize!(fig.layout, 2, Fixed(380))
    rowsize!(fig.layout, 1, Auto(1))
    save(gallery_path("$(name)_$(dark ? "dark" : "light").png"), fig)
end

# 0. Quickstart plot matching docs/src/quickstart.md.
quick_x = collect(range(0.0, 10.0; length=24))
quick_sigma_y = @. 0.16 + 0.02 * quick_x
quick_y = @. 1.85 * quick_x + 0.7 + quick_sigma_y * sin(1.6 * quick_x)
quick_model(x, p) = @. p[1] * x + p[2]
quick_result = fit_model(quick_model, quick_x, quick_y; p0=[1.0, 0.0], sigma_y=quick_sigma_y)
lightdark_plot(
    quick_result,
    "quickstart_linear";
    title=L"\mathrm{Quickstart\ calibration}",
    model_label=L"U(x)=m x+b",
    xlabel=L"t",
    xunit=L"\mathrm{s}",
    ylabel=L"U",
    yunit=L"\mathrm{V}",
    parameter_names=[L"m", L"b"],
    latex_labels=true,
    latex_stats=true,
    band=:prediction,
    nsigma=1,
    band_label=L"1\sigma\ \mathrm{prediction\ band}",
    show_legend=true,
    stats_position=:right,
    stats_mode=:full,
    stats_fontsize=20,
    figure_size=(1200, 760),
)

# The same scientific content rendered through every public style preset.
for style in (:clean, :minimal, :paper, :publication, :latex, :dark)
    plot_fit(
        quick_result;
        theme=style,
        filename=gallery_path("plot_style_$(style).png"),
        format=:png,
        title=L"\mathrm{Quickstart\ calibration}",
        model_label=L"U(t)=m t+b",
        xlabel=L"t",
        xunit=L"\mathrm{s}",
        ylabel=L"U",
        yunit=L"\mathrm{V}",
        parameter_names=[L"m", L"b"],
        latex_labels=true,
        latex_stats=true,
        band=:prediction,
        nsigma=1,
        band_label=L"1\sigma\ \mathrm{prediction\ band}",
        show_legend=true,
        stats_position=:right,
        stats_mode=:full,
        stats_fontsize=20,
        figure_size=(1200, 760),
    )
end

# 1. Linear calibration with visible heteroscedastic uncertainties.
x = collect(range(0.0, 10.0; length=28))
scatter_scale = @. 0.22 + 0.025 * x
sigma_y = @. 0.10 + 0.012 * x
y = @. 0.82 + 1.72 * x + scatter_scale * (0.55 * sin(1.25 * x) + 0.18 * cos(3.7 * x))
calibration_model(x, p) = @. p[1] * x + p[2]
linear_result = fit_model(calibration_model, x, y; p0=[1.5, 0.5], sigma_y=sigma_y)
lightdark_plot(
    linear_result,
    "linear_calibration";
    title=L"\mathrm{Sensor\ calibration}",
    model_label=L"U(x)=m x + b",
    xlabel=L"x",
    xunit=L"\mathrm{mm}",
    ylabel=L"U",
    yunit=L"\mathrm{V}",
    parameter_names=[L"m", L"b"],
    latex_labels=true,
    latex_stats=true,
    band=:prediction,
    nsigma=1,
    band_label=L"1\sigma\ \mathrm{prediction\ band}",
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
save_photoelectric_work_function(emission_result, baseline_result, frequency_THz, voltage, sigma_frequency_THz, sigma_voltage, emission_mask, "photoelectric_threshold"; dark=false)
save_photoelectric_work_function(emission_result, baseline_result, frequency_THz, voltage, sigma_frequency_THz, sigma_voltage, emission_mask, "photoelectric_threshold"; dark=true)

# 3. Exponential decay with full covariance.
x_cov = collect(range(0.0, 2.5; length=22))
decay_model(x, p) = @. p[1] * exp(p[2] * x) + p[3]
n = length(x_cov)
base_sigma = 0.055
corr_len = 2.3
cov_y = [base_sigma^2 * exp(-abs(i - j) / corr_len) for i in 1:n, j in 1:n]
y_cov = decay_model(x_cov, [2.0, -1.12, 0.24]) .+ 0.75 .* base_sigma .* (sin.(1.8 .* x_cov) .+ 0.28 .* cos.(4.1 .* x_cov))
cov_result = fit_model(decay_model, x_cov, y_cov; p0=[1.5, -0.7, 0.0], cov_y=cov_y)
lightdark_plot(
    cov_result,
    "full_covariance_decay";
    title=L"\mathrm{Correlated\ decay\ data}",
    model_label=L"y(t)=A e^{-\lambda t}+C",
    xlabel=L"t",
    xunit=L"\mathrm{s}",
    ylabel=L"signal",
    parameter_names=[L"A", L"\lambda", L"C"],
    latex_labels=true,
    latex_stats=true,
    band=:prediction,
    nsigma=1,
    band_label=L"1\sigma\ \mathrm{prediction\ band}",
    show_legend=true,
    legend_position=:lt,
    stats_position=:right,
    stats_mode=:full,
    stats_fontsize=20,
    figure_size=(1200, 760),
)

# 4. Effective-variance fit with x and y uncertainties.
x_true = collect(range(0.0, 4.0; length=18))
line_model(x, p) = @. p[1] * x + p[2]
sigma_x = fill(0.16, length(x_true))
sigma_y_xy = fill(0.10, length(x_true))
x_obs = x_true .+ sigma_x .* cos.(2.2 .* x_true)
y_obs = line_model(x_obs, [0.9, 1.2]) .+ sigma_y_xy .* sin.(3.1 .* x_obs)
xy_result = fit_model(line_model, x_obs, y_obs; p0=[0.5, 0.5], sigma_y=sigma_y_xy, sigma_x=sigma_x)
lightdark_plot(
    xy_result,
    "xy_uncertainties";
    title=L"\mathrm{XY\ uncertainty\ propagation}",
    model_label=L"y=m x+b",
    xlabel=L"x_\mathrm{meas}",
    ylabel=L"y_\mathrm{meas}",
    parameter_names=[L"m", L"b"],
    latex_labels=true,
    latex_stats=true,
    band=:prediction,
    nsigma=1,
    band_label=L"1\sigma\ \mathrm{prediction\ band}",
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
lightdark_plot(
    saturation_result,
    "constraints_priors";
    title=L"\mathrm{Early\ saturation\ measurement}",
    model_label=L"y(t)=A(1-e^{-t/\tau})+c",
    xlabel=L"t\;(\mathrm{s})",
    ylabel=L"y\;(\mathrm{V})",
    parameter_names=[L"A", L"\tau", L"c"],
    latex_labels=true,
    latex_stats=true,
    band=:prediction,
    nsigma=1,
    band_label=L"1\sigma\ \mathrm{prediction\ band}",
    show_legend=true,
    legend_position=:lt,
    stats_position=:right,
    stats_mode=:full,
    stats_fontsize=20,
    figure_size=(1200, 760),
)
prof = JuFitter.profile(saturation_result, 1; npoints=61, nsigma=4)
cont = JuFitter.contour(saturation_result, 1, 2; npoints=121, nsigma=4)
plot_profile(
    prof;
    theme=:minimal,
    title="Profile cost versus local parabola",
    xlabel="amplitude A",
    local_sigma=saturation_result.param_stderr[1],
    delta_max=8,
    filename=gallery_path("saturation_profile_light.png"),
    format=:png,
)
plot_profile(
    prof;
    theme=:dark,
    title="Profile cost versus local parabola",
    xlabel="amplitude A",
    line_color="#66d9ef",
    local_sigma=saturation_result.param_stderr[1],
    local_color="#b8c1ca",
    threshold_color="#edf2f4",
    delta_max=8,
    filename=gallery_path("saturation_profile_dark.png"),
    format=:png,
)
plot_contour(
    cont;
    theme=:minimal,
    title="Profile contours versus local covariance",
    xlabel="amplitude A",
    ylabel="time constant τ",
    local_covariance=saturation_result.param_covariance,
    local_center=saturation_result.params[[1, 2]],
    figure_size=(980, 720),
    filename=gallery_path("amplitude_timescale_contour_light.png"),
    format=:png,
)
plot_contour(
    cont;
    theme=:dark,
    title="Profile contours versus local covariance",
    xlabel="amplitude A",
    ylabel="time constant τ",
    local_covariance=saturation_result.param_covariance,
    local_center=saturation_result.params[[1, 2]],
    local_line_color="#b8c1ca",
    figure_size=(980, 720),
    filename=gallery_path("amplitude_timescale_contour_dark.png"),
    format=:png,
)

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
save_poisson_counts(poisson_result, x_counts, counts, poisson_model, "poisson_counts"; dark=false)
save_poisson_counts(poisson_result, x_counts, counts, poisson_model, "poisson_counts"; dark=true)

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
save_histogram_fit(hist_result, edges, hist_counts, expected_counts, "histogram_likelihood"; dark=false)
save_histogram_fit(hist_result, edges, hist_counts, expected_counts, "histogram_likelihood"; dark=true)

# 7. Multi-dataset fit with shared and local parameters.
x1 = collect(0.0:1.0:5.0)
x2 = collect(0.0:1.0:5.0)
local_linear(x, p) = @. p[1] * x + p[2]
sigma1 = fill(0.11, length(x1))
sigma2 = fill(0.13, length(x2))
y1 = local_linear(x1, [2.0, 1.0]) .+ sigma1 .* sin.(1.4 .* x1)
y2 = local_linear(x2, [2.0, -1.0]) .+ sigma2 .* cos.(1.7 .* x2)
multi_result = fit_multi_model(
    [local_linear, local_linear],
    [x1, x2],
    [y1, y2];
    p0=[1.7, 0.8, -0.8],
    parameter_map=[[1, 2], [1, 3]],
    sigma_y=[sigma1, sigma2],
    parameter_names=["shared slope", "offset A", "offset B"],
)
save_multi_dataset(multi_result, x1, y1, x2, y2, sigma1, sigma2, "multi_dataset_shared_slope"; dark=false)
save_multi_dataset(multi_result, x1, y1, x2, y2, sigma1, sigma2, "multi_dataset_shared_slope"; dark=true)

println("Generated documentation gallery assets in ", DOC_ASSET_DIR)
