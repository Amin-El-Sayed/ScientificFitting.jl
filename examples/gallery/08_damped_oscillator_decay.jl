const SNAPSHOT_ONLY = get(ENV, "JUFITTER_DOC_SNAPSHOT_ONLY", "0") == "1"
const RENDER_DOC_ASSETS = !SNAPSHOT_ONLY

if RENDER_DOC_ASSETS
    using CairoMakie
end

using JuFitter
using LaTeXStrings
using LinearAlgebra
using Printf

const DATA_FILE = joinpath(@__DIR__, "..", "data", "damped_oscillator", "pohl_wheel_free_decay.csv")
const OUTPUT_DIR = joinpath(@__DIR__, "..", "output")
const DOC_ASSET_DIR = joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "gallery")
const EMIT_DOC_OUTPUT_SNAPSHOTS = get(ENV, "JUFITTER_DOC_OUTPUT_SNAPSHOTS", "0") == "1"

function load_damped_oscillator(path)
    rows = readlines(path)[2:end]
    parsed = [parse.(Float64, split(row, ",")) for row in rows if !isempty(strip(row))]
    return (
        time=[row[1] for row in parsed],
        angle=[row[2] for row in parsed],
        sigma_angle=[row[3] for row in parsed],
    )
end

data = load_damped_oscillator(DATA_FILE)
time = data.time
angle = data.angle
sigma_angle = 0.5 .* data.sigma_angle
sigma_time = fill(0.0005, length(time))
time_reference = (minimum(time) + maximum(time)) / 2

# Centering time makes phase, frequency, and frequency drift easier to separate.
constant_frequency_model(t, p) = @. p[1] * exp(-p[4] * (t - time_reference)) *
                                     cos(p[2] * (t - time_reference) + p[3])

frequency_drift_model(t, p) = @. p[1] * exp(-p[4] * (t - time_reference)) *
                                  cos(
    p[2] * (t - time_reference) +
    0.5 * p[5] * (t - time_reference)^2 +
    p[3],
)

constant_result = fit_model(
    constant_frequency_model,
    time,
    angle;
    p0=[1.6, 3.26, 0.0, 0.0035],
    sigma_y=sigma_angle,
    sigma_x=sigma_time,
    bounds=([0.0, 2.0, -20.0, 0.0], [5.0, 5.0, 20.0, 0.05]),
    initial_guesses=[
        [1.6, 3.26, 0.0, 0.0035],
        [1.8, 3.20, 2.0, 0.0020],
        [1.5, 3.35, -2.0, 0.0060],
    ],
    maxiters=3000,
)

drift_result = fit_model(
    frequency_drift_model,
    time,
    angle;
    p0=[1.6, 3.26, 0.0, 0.0035, 0.0],
    sigma_y=sigma_angle,
    sigma_x=sigma_time,
    bounds=([0.0, 2.0, -20.0, 0.0, -0.01], [5.0, 5.0, 20.0, 0.05, 0.01]),
    initial_guesses=[
        [1.6, 3.26, 0.0, 0.0035, 0.0],
        [1.8, 3.20, 2.0, 0.0020, 0.0001],
        [1.5, 3.35, -2.0, 0.0060, -0.0001],
    ],
    maxiters=4000,
)

function drift_prediction_sigma(result, t, sigma_y, sigma_t)
    amplitude, omega, phase, damping, beta = result.params
    covariance = result.param_covariance
    tau = t - time_reference
    envelope = amplitude * exp(-damping * tau)
    theta = omega * tau + 0.5 * beta * tau^2 + phase

    jacobian = [
        exp(-damping * tau) * cos(theta),
        -envelope * tau * sin(theta),
        -envelope * sin(theta),
        -envelope * tau * cos(theta),
        -0.5 * envelope * tau^2 * sin(theta),
    ]
    confidence_variance = max(dot(jacobian, covariance * jacobian), 0.0)
    time_derivative = envelope * (-damping * cos(theta) - (omega + beta * tau) * sin(theta))
    return sqrt(confidence_variance + sigma_y^2 + (time_derivative * sigma_t)^2)
end

fmt(x, digits=4) = @sprintf("%.*g", digits, x)
function fmt_tex(x, digits=4)
    value = fmt(x, digits)
    scientific = match(r"^(.+)[eE]([+-]?\d+)$", value)
    scientific === nothing && return value
    mantissa, exponent = scientific.captures
    return mantissa * "\\times10^{" * string(parse(Int, exponent)) * "}"
end

function style_asset_suffix(dark, style::Symbol, appearance::Symbol)
    return dark === nothing ? "$(style)_$(appearance)" : (appearance == :dark ? "dark" : "light")
end

function emit_doc_output_snapshot(body::Function, id::AbstractString)
    EMIT_DOC_OUTPUT_SNAPSHOTS || return nothing

    println("=== JUFITTER_DOC_OUTPUT_BEGIN ", id, " ===")
    body()
    println("=== JUFITTER_DOC_OUTPUT_END ", id, " ===")
    return nothing
end

function add_pull_reference!(axis, xmin, xmax, color_1sigma, color_2sigma, zero_color)
    band!(axis, [xmin, xmax], [-2.0, -2.0], [2.0, 2.0]; color=color_2sigma)
    band!(axis, [xmin, xmax], [-1.0, -1.0], [1.0, 1.0]; color=color_1sigma)
    hlines!(axis, [0.0]; color=zero_color, linewidth=1.4)
end

function save_model_comparison(
    filename;
    dark::Union{Nothing, Bool}=nothing,
    style::Symbol=:modern,
    appearance::Symbol=dark === nothing ? :light : (dark ? :dark : :light),
)
    RENDER_DOC_ASSETS || return nothing

    dark_mode = appearance == :dark
    palette = plot_palette(style; appearance=appearance)
    foreground = palette.stats_color
    muted = palette.stats_muted_color
    constant_color = style == :article ? (dark_mode ? "#c8cdd4" : "#505761") :
                     (style == :modern ? (dark_mode ? "#b7c8dc" : "#52606f") : muted)
    drift_color = palette.fit_color
    measurement_color = palette.data_color
    prediction_color = (palette.band_color, palette.band_alpha)
    pull_1sigma = (palette.band_color, dark_mode ? 0.15 : 0.24)
    pull_2sigma = (palette.band_color, dark_mode ? 0.07 : 0.11)
    background = dark_mode ? "#111318" : (style == :modern ? "#fbfcfd" : "#ffffff")
    article = style == :article
    fit_title = "Free decay: constant frequency versus frequency drift"
    angle_label = article ? L"\varphi\;(\mathrm{rad})" : "angle φ (rad)"
    pull_label = article ? L"r_i" : "pull rᵢ"
    time_label = article ? L"t\;(\mathrm{s})" : "elapsed time t (s)"

    figure = with_theme(plot_theme(style; appearance=appearance)) do
        Figure(size=(1440, 960), backgroundcolor=background)
    end
    fit_axis = Axis(
        figure[1, 1];
        title=fit_title,
        ylabel=angle_label,
    )
    constant_pull_axis = Axis(
        figure[2, 1];
        title="Pulls: constant-frequency model",
        titlealign=:left,
        ylabel=pull_label,
    )
    drift_pull_axis = Axis(
        figure[3, 1];
        title="Pulls: frequency-drift model",
        titlealign=:left,
        xlabel=time_label,
        ylabel=pull_label,
    )

    time_grid = collect(range(minimum(time), maximum(time); length=1800))
    constant_grid = constant_frequency_model(time_grid, constant_result.params)
    drift_grid = frequency_drift_model(time_grid, drift_result.params)
    prediction_sigma = [
        drift_prediction_sigma(drift_result, t, sigma_angle[1], sigma_time[1])
        for t in time_grid
    ]
    band!(
        fit_axis,
        time_grid,
        drift_grid .- prediction_sigma,
        drift_grid .+ prediction_sigma;
        color=prediction_color,
        label="drift model: local 1σ prediction band",
    )
    lines!(
        fit_axis,
        time_grid,
        constant_grid;
        color=constant_color,
        linestyle=:dash,
        linewidth=max(1.8, palette.fit_linewidth - 0.3),
        label="constant-frequency model",
    )
    lines!(
        fit_axis,
        time_grid,
        drift_grid;
        color=drift_color,
        linewidth=palette.fit_linewidth,
        label="frequency-drift model",
    )
    errorbars!(fit_axis, time, angle, sigma_angle; color=palette.yerr_color, whiskerwidth=palette.error_whiskerwidth)
    errorbars!(
        fit_axis,
        time,
        angle,
        sigma_time;
        direction=:x,
        color=palette.xerr_color,
        whiskerwidth=palette.error_whiskerwidth,
    )
    scatter!(
        fit_axis,
        time,
        angle;
        color=measurement_color,
        markersize=max(3.8, palette.data_markersize - 1.4),
        label="measured angle",
    )
    hidexdecorations!(fit_axis; grid=false)

    xmin, xmax = extrema(time)
    add_pull_reference!(constant_pull_axis, xmin, xmax, pull_1sigma, pull_2sigma, (foreground, 0.55))
    add_pull_reference!(drift_pull_axis, xmin, xmax, pull_1sigma, pull_2sigma, (foreground, 0.55))
    lines!(constant_pull_axis, time, constant_result.weighted_residuals; color=(constant_color, 0.50), linewidth=1.2)
    scatter!(constant_pull_axis, time, constant_result.weighted_residuals; color=constant_color, markersize=4.8)
    lines!(drift_pull_axis, time, drift_result.weighted_residuals; color=(drift_color, 0.50), linewidth=1.2)
    scatter!(drift_pull_axis, time, drift_result.weighted_residuals; color=drift_color, markersize=4.8)
    hidexdecorations!(constant_pull_axis; grid=false)
    linkxaxes!(fit_axis, constant_pull_axis, drift_pull_axis)
    ylims!(constant_pull_axis, -3.2, 3.2)
    ylims!(drift_pull_axis, -3.2, 3.2)

    beta = drift_result.params[5]
    sigma_beta = drift_result.param_stderr[5]
    time_span = maximum(time) - minimum(time)
    frequency_change = beta * time_span
    sigma_frequency_change = sigma_beta * time_span
    damping_time = 1 / drift_result.params[4]
    sigma_damping_time = drift_result.param_stderr[4] / drift_result.params[4]^2
    delta_aic = constant_result.stats.aic - drift_result.stats.aic

    parameter_lines = article ? Any[
        LaTeXString("\\omega_{\\mathrm{ref}} = $(fmt_tex(drift_result.params[2], 6)) \\pm $(fmt_tex(drift_result.param_stderr[2], 2))\\;\\mathrm{rad\\,s^{-1}}"),
        LaTeXString("\\beta = $(fmt_tex(1e3 * beta, 5)) \\pm $(fmt_tex(1e3 * sigma_beta, 2))\\;\\mathrm{mrad\\,s^{-2}}"),
        LaTeXString("\\Delta\\omega_{60s} = $(fmt_tex(frequency_change, 4)) \\pm $(fmt_tex(sigma_frequency_change, 2))\\;\\mathrm{rad\\,s^{-1}}"),
        LaTeXString("\\lambda = $(fmt_tex(drift_result.params[4], 5)) \\pm $(fmt_tex(drift_result.param_stderr[4], 2))\\;\\mathrm{s^{-1}}"),
        LaTeXString("\\tau_d = $(fmt_tex(damping_time, 5)) \\pm $(fmt_tex(sigma_damping_time, 2))\\;\\mathrm{s}"),
    ] : Any[
        "ωref = $(fmt(drift_result.params[2], 6)) ± $(fmt(drift_result.param_stderr[2], 2)) rad s⁻¹",
        "β = $(fmt(1e3 * beta, 5)) ± $(fmt(1e3 * sigma_beta, 2)) mrad s⁻²",
        "Δω over 60 s = $(fmt(frequency_change, 4)) ± $(fmt(sigma_frequency_change, 2)) rad s⁻¹",
        "λ = $(fmt(drift_result.params[4], 5)) ± $(fmt(drift_result.param_stderr[4], 2)) s⁻¹",
        "τd = $(fmt(damping_time, 5)) ± $(fmt(sigma_damping_time, 2)) s",
    ]
    statistic_lines = article ? Any[
        LaTeXString("\\chi^2_{\\mathrm{const}}/\\mathrm{ndf} = $(fmt_tex(constant_result.stats.chi2_ndf, 4))"),
        LaTeXString("P_{\\mathrm{const}}(\\chi^2) = $(fmt_tex(constant_result.stats.pvalue, 3))"),
        "constant diagnostic = $(diagnostic_dashboard(constant_result).status)",
        LaTeXString("\\chi^2_{\\mathrm{drift}}/\\mathrm{ndf} = $(fmt_tex(drift_result.stats.chi2_ndf, 4))"),
        LaTeXString("P_{\\mathrm{drift}}(\\chi^2) = $(fmt_tex(drift_result.stats.pvalue, 3))"),
        "drift diagnostic = $(diagnostic_dashboard(drift_result).status)",
        LaTeXString("\\Delta\\mathrm{AIC} = $(fmt_tex(delta_aic, 5))\\;\\mathrm{in\\ favor\\ of\\ drift}"),
        LaTeXString("\\chi^2/\\mathrm{ndf} \\ll 1\\;\\mathrm{still\\ requires\\ review}"),
    ] : Any[
        "constant model: χ²/ndf = $(fmt(constant_result.stats.chi2_ndf, 4))",
        "constant model: P(χ²) = $(fmt(constant_result.stats.pvalue, 3))",
        "constant diagnostic = $(diagnostic_dashboard(constant_result).status)",
        "drift model: χ²/ndf = $(fmt(drift_result.stats.chi2_ndf, 4))",
        "drift model: P(χ²) = $(fmt(drift_result.stats.pvalue, 3))",
        "drift diagnostic = $(diagnostic_dashboard(drift_result).status)",
        "ΔAIC = $(fmt(delta_aic, 5)) in favor of drift",
        "χ²/ndf ≪ 1 still requires review",
    ]
    plot_info_panel!(
        figure[1:3, 2];
        legend_source=fit_axis,
        title="Frequency-drift result",
        model_label=article ? L"A e^{-\lambda t}\cos(\omega_{\mathrm{ref}}t+\beta t^2/2+\varphi)" :
            "A exp(-λt) cos(ωref t + βt²/2 + φ)",
        parameter_lines=parameter_lines,
        statistic_lines=statistic_lines,
        color=foreground,
        muted_color=muted,
        fontsize=palette.stats_fontsize + 1,
    )

    rowsize!(figure.layout, 1, Relative(0.62))
    rowsize!(figure.layout, 2, Relative(0.19))
    rowsize!(figure.layout, 3, Relative(0.19))
    colsize!(figure.layout, 2, Fixed(400))
    colgap!(figure.layout, 18)
    save(filename, figure)
end

if RENDER_DOC_ASSETS
    mkpath(OUTPUT_DIR)
    mkpath(DOC_ASSET_DIR)

    for (dark, suffix) in ((false, "light"), (true, "dark"))
        save_model_comparison(joinpath(OUTPUT_DIR, "08_damped_oscillator_decay_$(suffix).png"); dark=dark)
        save_model_comparison(joinpath(DOC_ASSET_DIR, "damped_oscillator_decay_$(suffix).png"); dark=dark)
    end

    for style in (:lab, :modern, :article), appearance in (:light, :dark)
        save_model_comparison(
            joinpath(DOC_ASSET_DIR, "damped_oscillator_decay_$(style_asset_suffix(nothing, style, appearance)).png");
            style=style,
            appearance=appearance,
        )
    end
end

println("Constant-frequency model")
println(report_text(constant_result; parameter_names=["A_ref", "omega_ref", "phi_ref", "lambda"]))
println(diagnostic_dashboard_text(constant_result))
emit_doc_output_snapshot("resonance_constant") do
    println(report_text(constant_result; parameter_names=["A_ref", "omega_ref", "phi_ref", "lambda"]))
    println(diagnostic_dashboard_text(constant_result))
end
println()
println("Frequency-drift model")
println(report_text(drift_result; parameter_names=["A_ref", "omega_ref", "phi_ref", "lambda", "beta"]))
println(diagnostic_dashboard_text(drift_result))
emit_doc_output_snapshot("resonance_drift") do
    println(report_text(drift_result; parameter_names=["A_ref", "omega_ref", "phi_ref", "lambda", "beta"]))
    println(diagnostic_dashboard_text(drift_result))
end
