using CairoMakie
using JuFitter
using LaTeXStrings
using LinearAlgebra
using Printf

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
            fontsize=14,
            font="TeX Gyre Heros",
            Figure=(backgroundcolor="#111318",),
            Axis=(
                backgroundcolor="#111318",
                xlabelsize=18,
                ylabelsize=18,
                titlesize=18,
                xticklabelsize=13,
                yticklabelsize=13,
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
            Legend=(framevisible=false, labelcolor="#edf2f4", labelsize=13, patchsize=(18, 10)),
        )
    end
    return Theme(
        fontsize=14,
        font="TeX Gyre Heros",
        Figure=(backgroundcolor="#ffffff",),
        Axis=(
            backgroundcolor="#ffffff",
            xlabelsize=18,
            ylabelsize=18,
            titlesize=18,
            xticklabelsize=13,
            yticklabelsize=13,
            xgridcolor=("#eef2f7", 0.9),
            ygridcolor=("#eef2f7", 0.9),
            topspinevisible=false,
            rightspinevisible=false,
        ),
        Legend=(framevisible=false, labelsize=13, patchsize=(18, 10)),
    )
end

function save_poisson_counts(result, x, counts, model, name; dark::Bool=false)
    color = dark ? "#edf2f4" : "#14151a"
    fit_color = dark ? "#66d9ef" : "#0081a7"
    band_color = dark ? ("#66d9ef", 0.18) : ("#a8dadc", 0.30)
    fig = with_theme(gallery_theme(dark)) do
        Figure(size=(1040, 610), backgroundcolor=dark ? "#111318" : "#ffffff")
    end
    ax = Axis(fig[1, 1]; title="Poisson count model", xlabel="channel", ylabel="counts")
    yerr = sqrt.(counts)
    errorbars!(ax, x, counts, yerr; color=(color, 0.45), whiskerwidth=6)
    scatter!(ax, x, counts; color=(color, 0.78), markersize=9, label="data")
    xg = collect(range(minimum(x), maximum(x); length=300))
    yg = model(xg, result.params)
    band!(ax, xg, yg .- sqrt.(yg), yg .+ sqrt.(yg); color=band_color)
    lines!(ax, xg, yg; color=fit_color, linewidth=3, label="model")
    axislegend(ax; position=:lt)
    panel = Axis(fig[1, 2]; backgroundcolor=:transparent)
    hidedecorations!(panel)
    hidespines!(panel)
    text!(
        panel,
        0,
        1;
        text="log scale = $(round(result.params[1]; sigdigits=5))\n" *
             "slope = $(round(result.params[2]; sigdigits=5))\n" *
             "χ²/ndf = $(round(result.stats.chi2_ndf; sigdigits=4))\n" *
             "P(χ²) = $(round(result.stats.pvalue; sigdigits=4))",
        space=:relative,
        align=(:left, :top),
        color=color,
        fontsize=15,
        lineheight=1.1,
    )
    colsize!(fig.layout, 2, Fixed(220))
    rowsize!(fig.layout, 1, Auto(1))
    save(gallery_path("$(name)_$(dark ? "dark" : "light").png"), fig)
end

function save_histogram_fit(result, edges, counts, expected_counts, name; dark::Bool=false)
    color = dark ? "#edf2f4" : "#14151a"
    fit_color = dark ? "#66d9ef" : "#0081a7"
    fig = with_theme(gallery_theme(dark)) do
        Figure(size=(1040, 610), backgroundcolor=dark ? "#111318" : "#ffffff")
    end
    ax = Axis(fig[1, 1]; title="Histogram likelihood fit", xlabel="observable", ylabel="bin counts")
    centers = [(edges[i] + edges[i + 1]) / 2 for i in 1:(length(edges) - 1)]
    widths = diff(edges)
    barplot!(ax, centers, counts; width=0.92 .* widths, color=(fit_color, 0.22), strokecolor=(color, 0.30), strokewidth=1)
    scatter!(ax, centers, counts; color=color, markersize=10)
    expected = expected_counts(edges, result.params)
    lines!(ax, centers, expected; color=fit_color, linewidth=3)
    panel = Axis(fig[1, 2]; backgroundcolor=:transparent)
    hidedecorations!(panel)
    hidespines!(panel)
    text!(
        panel,
        0,
        1;
        text="scale = $(round(result.params[1]; sigdigits=5))\n" *
             "slope = $(round(result.params[2]; sigdigits=5))\n" *
             "χ²/ndf = $(round(result.stats.chi2_ndf; sigdigits=4))\n" *
             "P(χ²) = $(round(result.stats.pvalue; sigdigits=4))",
        space=:relative,
        align=(:left, :top),
        color=color,
        fontsize=15,
        lineheight=1.1,
    )
    colsize!(fig.layout, 2, Fixed(220))
    rowsize!(fig.layout, 1, Auto(1))
    save(gallery_path("$(name)_$(dark ? "dark" : "light").png"), fig)
end

function save_photoelectric_work_function(result, frequency, voltage, sigma_frequency, sigma_voltage, fit_mask, name; dark::Bool=false)
    color = dark ? "#edf2f4" : "#14151a"
    muted = dark ? "#b8c1ca" : "#5b6270"
    fit_color = dark ? "#66d9ef" : "#0081a7"
    work_color = dark ? "#f59e0b" : "#b45309"
    band_color = dark ? ("#66d9ef", 0.20) : ("#a8dadc", 0.28)

    slope, intercept = result.params
    cov = result.param_covariance
    sigma_intercept = sqrt(max(cov[2, 2], 0.0))
    threshold = -intercept / slope
    threshold_grad = [intercept / slope^2, -1 / slope]
    sigma_threshold = sqrt(max(dot(threshold_grad, cov * threshold_grad), 0.0))
    work_function_eV = -intercept

    fig = with_theme(gallery_theme(dark)) do
        Figure(size=(1320, 760), backgroundcolor=dark ? "#111318" : "#ffffff")
    end
    ax = Axis(
        fig[1, 1];
        title="Photoelectric work-function extraction",
        xlabel="frequency ν (THz)",
        ylabel="stopping voltage U₀ (V)",
    )

    errorbars!(ax, frequency, voltage, sigma_voltage; color=(muted, 0.46), whiskerwidth=5)
    errorbars!(ax, frequency, voltage, sigma_frequency; direction=:x, color=(muted, 0.30), whiskerwidth=5)
    scatter!(ax, frequency[fit_mask], voltage[fit_mask]; color=color, markersize=9, label="fit data")
    scatter!(ax, frequency[.!fit_mask], voltage[.!fit_mask]; color=(muted, 0.55), marker=:diamond, markersize=9, label="below threshold")

    xg = collect(range(0.0, maximum(frequency) * 1.04; length=500))
    yg = @. slope * xg + intercept
    J = hcat(xg, ones(length(xg)))
    sigma_band = sqrt.(clamp.(vec(sum((J * cov) .* J; dims=2)), 0.0, Inf))
    band!(ax, xg, yg .- sigma_band, yg .+ sigma_band; color=band_color, label="1σ fit uncertainty")
    lines!(ax, xg, yg; color=fit_color, linewidth=3.0, label="linear extrapolation")

    hlines!(ax, [0.0]; color=(muted, 0.55), linestyle=:dash)
    vlines!(ax, [threshold]; color=(work_color, 0.75), linestyle=:dash)
    scatter!(ax, [0.0], [intercept]; color=work_color, markersize=12, label="y-intercept")
    scatter!(ax, [threshold], [0.0]; color=fit_color, marker=:utriangle, markersize=12, label="threshold frequency")
    text!(
        ax,
        35,
        intercept - 0.16;
        text="Φ = $(round(work_function_eV; sigdigits=4)) ± $(round(sigma_intercept; sigdigits=2)) eV",
        color=work_color,
        fontsize=15,
        align=(:left, :top),
    )
    text!(
        ax,
        threshold + 30,
        0.22;
        text="ν₀ = $(round(threshold; sigdigits=4)) ± $(round(sigma_threshold; sigdigits=2)) THz",
        color=fit_color,
        fontsize=15,
        align=(:left, :bottom),
    )
    axislegend(ax; position=:lt, nbanks=2)
    limits!(ax, 0, maximum(frequency) * 1.05, intercept - 0.55, maximum(voltage) + 0.65)

    h_fit = slope * 1.602176634e-19 / 1e12
    sigma_h = sqrt(max(cov[1, 1], 0.0)) * 1.602176634e-19 / 1e12
    panel = Axis(fig[1, 2]; backgroundcolor=:transparent)
    hidedecorations!(panel)
    hidespines!(panel)
    text!(
        panel,
        0,
        1;
        text="slope h/e = $(fmt_sig(slope, 5)) V/THz\n" *
             "h = $(fmt_sig(h_fit, 5)) ± $(fmt_sig(sigma_h, 2)) J s\n" *
             "intercept = $(fmt_sig(intercept, 5)) ± $(fmt_sig(sigma_intercept, 2)) V\n" *
             "Φ = $(fmt_sig(work_function_eV, 5)) ± $(fmt_sig(sigma_intercept, 2)) eV\n" *
             "ν₀ = $(fmt_sig(threshold, 5)) ± $(fmt_sig(sigma_threshold, 2)) THz\n" *
             "χ²/ndf = $(fmt_sig(result.stats.chi2_ndf, 4))",
        space=:relative,
        align=(:left, :top),
        color=color,
        fontsize=15,
        lineheight=1.12,
    )
    colsize!(fig.layout, 2, Fixed(330))
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
        Figure(size=(1320, 760), backgroundcolor=dark ? "#111318" : "#ffffff")
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
        fontsize=16,
        lineheight=1.15,
    )
    colsize!(fig.layout, 2, Fixed(330))
    rowsize!(fig.layout, 1, Auto(1))
    save(gallery_path("$(name)_$(dark ? "dark" : "light").png"), fig)
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
    stats_fontsize=13,
    figure_size=(1320, 760),
)

# 2. Photoelectric work-function extraction with x/y uncertainties.
const c = 299_792_458.0
wavelength_nm = [150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0]
voltage = [5.99, 3.87, 2.69, 1.78, 1.28, 0.77, 0.50, 0.15, 0.0, 0.0]
sigma_wavelength_nm = fill(0.01, length(wavelength_nm))
sigma_voltage = fill(0.04, length(voltage))
frequency_THz = @. c / (wavelength_nm * 1e-9) / 1e12
sigma_frequency_THz = @. c * (sigma_wavelength_nm * 1e-9) / (wavelength_nm * 1e-9)^2 / 1e12
fit_mask = voltage .> 0.0
photo_model(f, p) = @. p[1] * f + p[2]
photo_result = fit_model(
    photo_model,
    frequency_THz[fit_mask],
    voltage[fit_mask];
    p0=[0.004, -2.2],
    sigma_y=sigma_voltage[fit_mask],
    sigma_x=sigma_frequency_THz[fit_mask],
    bounds=([0.0, -20.0], [0.02, 5.0]),
    initial_guesses=[[0.004, -2.2], [0.0042, -2.6], [0.0038, -1.8]],
)
save_photoelectric_work_function(photo_result, frequency_THz, voltage, sigma_frequency_THz, sigma_voltage, fit_mask, "photoelectric_threshold"; dark=false)
save_photoelectric_work_function(photo_result, frequency_THz, voltage, sigma_frequency_THz, sigma_voltage, fit_mask, "photoelectric_threshold"; dark=true)

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
    stats_fontsize=13,
    figure_size=(1320, 760),
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
    stats_fontsize=13,
    figure_size=(1320, 760),
)

# 5. Constraints, prior, profile, and contour.
x_quad = collect(range(-2.0, 2.3; length=28))
quad_model(x, p) = @. p[1] * x^2 + p[2] * x + p[3]
sigma_quad = fill(0.08, length(x_quad))
y_quad = quad_model(x_quad, [0.65, -0.75, 0.35]) .+ sigma_quad .* cos.(2.0 .* x_quad)
quad_result = fit_model(
    quad_model,
    x_quad,
    y_quad;
    p0=[0.25, -0.2, 0.0],
    sigma_y=sigma_quad,
    bounds=([0.0, -2.0, -1.0], [2.0, 2.0, 2.0]),
    constraints=(ineq=p -> [-p[1]],),
    parameter_priors=(index=3, mean=0.3, sigma=0.2),
)
lightdark_plot(
    quad_result,
    "constraints_priors";
    title=L"\mathrm{Constrained\ quadratic\ model}",
    model_label=L"y=a x^2+b x+c",
    xlabel=L"x",
    ylabel=L"response",
    parameter_names=[L"a", L"b", L"c"],
    latex_labels=true,
    latex_stats=true,
    band=:prediction,
    nsigma=1,
    band_label=L"1\sigma\ \mathrm{prediction\ band}",
    show_legend=true,
    legend_position=:lt,
    stats_position=:right,
    stats_mode=:full,
    stats_fontsize=13,
    figure_size=(1320, 760),
)
prof = JuFitter.profile(quad_result, 1; npoints=45, nsigma=3)
cont = JuFitter.contour(quad_result, 1, 2; npoints=31, nsigma=2)
plot_profile(prof; theme=:minimal, title="Profile for curvature", xlabel="curvature", filename=gallery_path("curvature_profile_light.png"), format=:png)
plot_profile(prof; theme=:dark, title="Profile for curvature", xlabel="curvature", line_color="#66d9ef", threshold_color="#edf2f4", filename=gallery_path("curvature_profile_dark.png"), format=:png)
plot_contour(cont; theme=:minimal, title="Curvature-slope contour", xlabel="curvature", ylabel="slope", filename=gallery_path("curvature_slope_contour_light.png"), format=:png)
plot_contour(cont; theme=:dark, title="Curvature-slope contour", xlabel="curvature", ylabel="slope", line_color="#edf2f4", filename=gallery_path("curvature_slope_contour_dark.png"), format=:png)

# 6. Poisson counts and histogram likelihoods.
x_counts = collect(1.0:10.0)
counts = [2, 4, 5, 6, 9, 11, 12, 14, 15, 19]
poisson_model(x, p) = @. exp(p[1] + p[2] * x)
poisson_result = fit_poisson_model(
    poisson_model,
    x_counts,
    counts;
    p0=[0.0, 0.1],
    bounds=([-10.0, -10.0], [10.0, 10.0]),
    parameter_names=["log scale", "slope"],
)
save_poisson_counts(poisson_result, x_counts, counts, poisson_model, "poisson_counts"; dark=false)
save_poisson_counts(poisson_result, x_counts, counts, poisson_model, "poisson_counts"; dark=true)

edges = collect(0.0:1.0:5.0)
hist_counts = [4, 9, 13, 20, 31]
expected_counts(edges, p) = [p[1] * (edges[i + 1] - edges[i]) * exp(p[2] * (edges[i] + edges[i + 1]) / 2) for i in 1:(length(edges) - 1)]
hist_result = fit_histogram_model(
    expected_counts,
    edges,
    hist_counts;
    p0=[3.0, 0.3],
    bounds=([1e-6, -5.0], [100.0, 5.0]),
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
