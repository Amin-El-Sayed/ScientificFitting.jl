using CairoMakie
using JuFitter
using LinearAlgebra
using Printf

const DATA_FILE = joinpath(@__DIR__, "..", "data", "damped_oscillator", "pohl_wheel_free_decay.csv")
const OUTPUT_DIR = joinpath(@__DIR__, "..", "output")
const DOC_ASSET_DIR = joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "gallery")

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
sigma_angle = data.sigma_angle
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

function oscillator_theme(dark::Bool)
    foreground = dark ? "#edf2f4" : "#14151a"
    muted = dark ? "#b8c1ca" : "#5b6270"
    grid = dark ? ("#2a313a", 0.85) : ("#e9eef4", 0.95)
    background = dark ? "#111318" : "#ffffff"
    return Theme(
        fontsize=20,
        font="TeX Gyre Heros",
        Figure=(backgroundcolor=background,),
        Axis=(
            backgroundcolor=background,
            xlabelsize=27,
            ylabelsize=27,
            titlesize=26,
            xticklabelsize=20,
            yticklabelsize=20,
            xlabelcolor=foreground,
            ylabelcolor=foreground,
            titlecolor=foreground,
            xticklabelcolor=muted,
            yticklabelcolor=muted,
            xtickcolor=muted,
            ytickcolor=muted,
            xgridcolor=grid,
            ygridcolor=grid,
            topspinevisible=false,
            rightspinevisible=false,
            leftspinecolor=muted,
            bottomspinecolor=muted,
        ),
        Legend=(framevisible=false, labelcolor=foreground, labelsize=19, patchsize=(30, 16)),
    )
end

function add_pull_reference!(axis, xmin, xmax, color_1sigma, color_2sigma, zero_color)
    band!(axis, [xmin, xmax], [-2.0, -2.0], [2.0, 2.0]; color=color_2sigma)
    band!(axis, [xmin, xmax], [-1.0, -1.0], [1.0, 1.0]; color=color_1sigma)
    hlines!(axis, [0.0]; color=zero_color, linewidth=1.4)
end

function save_model_comparison(filename; dark::Bool=false)
    foreground = dark ? "#edf2f4" : "#14151a"
    muted = dark ? "#b8c1ca" : "#5b6270"
    constant_color = dark ? "#f4b860" : "#b45309"
    drift_color = dark ? "#66d9ef" : "#007f9e"
    measurement_color = dark ? "#edf2f4" : "#14151a"
    prediction_color = dark ? ("#66d9ef", 0.20) : ("#89d5e0", 0.34)
    pull_1sigma = dark ? ("#66d9ef", 0.13) : ("#a8dadc", 0.30)
    pull_2sigma = dark ? ("#66d9ef", 0.06) : ("#a8dadc", 0.14)
    background = dark ? "#111318" : "#ffffff"

    figure = with_theme(oscillator_theme(dark)) do
        Figure(size=(1720, 1040), backgroundcolor=background)
    end
    fit_axis = Axis(
        figure[1, 1];
        title="Free decay: constant frequency versus frequency drift",
        ylabel="angle φ (rad)",
    )
    constant_pull_axis = Axis(
        figure[2, 1];
        title="Pulls: constant-frequency model",
        titlealign=:left,
        titlesize=19,
        ylabel="pull rᵢ",
        ylabelsize=22,
    )
    drift_pull_axis = Axis(
        figure[3, 1];
        title="Pulls: frequency-drift model",
        titlealign=:left,
        titlesize=19,
        xlabel="elapsed time t (s)",
        ylabel="pull rᵢ",
        ylabelsize=22,
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
        linewidth=2.2,
        label="constant-frequency model",
    )
    lines!(
        fit_axis,
        time_grid,
        drift_grid;
        color=drift_color,
        linewidth=2.6,
        label="frequency-drift model",
    )
    errorbars!(fit_axis, time, angle, sigma_angle; color=(muted, 0.28), whiskerwidth=3)
    errorbars!(
        fit_axis,
        time,
        angle,
        sigma_time;
        direction=:x,
        color=(muted, 0.18),
        whiskerwidth=3,
    )
    scatter!(
        fit_axis,
        time,
        angle;
        color=(measurement_color, 0.72),
        markersize=4.4,
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

    side = GridLayout()
    figure[1:3, 2] = side
    Legend(side[1, 1], fit_axis; framevisible=false, tellheight=true)
    report_axis = Axis(side[2, 1]; backgroundcolor=:transparent)
    hidedecorations!(report_axis)
    hidespines!(report_axis)

    beta = drift_result.params[5]
    sigma_beta = drift_result.param_stderr[5]
    time_span = maximum(time) - minimum(time)
    frequency_change = beta * time_span
    sigma_frequency_change = sigma_beta * time_span
    damping_time = 1 / drift_result.params[4]
    sigma_damping_time = drift_result.param_stderr[4] / drift_result.params[4]^2
    delta_aic = constant_result.stats.aic - drift_result.stats.aic

    text!(
        report_axis,
        0,
        1;
        text="frequency-drift model\n" *
             "  ωref = $(fmt(drift_result.params[2], 6)) ± $(fmt(drift_result.param_stderr[2], 2)) rad s⁻¹\n" *
             "  β = $(fmt(1e3 * beta, 5)) ± $(fmt(1e3 * sigma_beta, 2)) mrad s⁻²\n" *
             "  Δω over 60 s\n" *
             "    $(fmt(frequency_change, 4)) ± $(fmt(sigma_frequency_change, 2)) rad s⁻¹\n" *
             "  λ = $(fmt(drift_result.params[4], 5)) ± $(fmt(drift_result.param_stderr[4], 2)) s⁻¹\n" *
             "  τd = $(fmt(damping_time, 5)) ± $(fmt(sigma_damping_time, 2)) s\n\n" *
             "constant-frequency model\n" *
             "  χ²/ndf = $(fmt(constant_result.stats.chi2_ndf, 4))\n" *
             "  P(χ²) = $(fmt(constant_result.stats.pvalue, 3))\n" *
             "  diagnostic status: $(diagnostic_dashboard(constant_result).status)\n\n" *
             "frequency-drift model\n" *
             "  χ²/ndf = $(fmt(drift_result.stats.chi2_ndf, 4))\n" *
             "  P(χ²) = $(fmt(drift_result.stats.pvalue, 3))\n" *
             "  diagnostic status: $(diagnostic_dashboard(drift_result).status)\n\n" *
             "ΔAIC = $(fmt(delta_aic, 5)) in favor of drift\n" *
             "but χ²/ndf ≪ 1 requires review",
        space=:relative,
        align=(:left, :top),
        color=foreground,
        fontsize=18,
        lineheight=1.08,
    )

    rowsize!(side, 1, Auto())
    rowsize!(side, 2, Relative(1))
    rowsize!(figure.layout, 1, Relative(0.62))
    rowsize!(figure.layout, 2, Relative(0.19))
    rowsize!(figure.layout, 3, Relative(0.19))
    colsize!(figure.layout, 2, Fixed(540))
    save(filename, figure)
end

mkpath(OUTPUT_DIR)
mkpath(DOC_ASSET_DIR)

for (dark, suffix) in ((false, "light"), (true, "dark"))
    save_model_comparison(joinpath(OUTPUT_DIR, "08_damped_oscillator_decay_$(suffix).png"); dark=dark)
    save_model_comparison(joinpath(DOC_ASSET_DIR, "damped_oscillator_decay_$(suffix).png"); dark=dark)
end

println("Constant-frequency model")
println(report_text(constant_result; parameter_names=["A_ref", "omega_ref", "phi_ref", "lambda"]))
println(diagnostic_dashboard_text(constant_result))
println()
println("Frequency-drift model")
println(report_text(drift_result; parameter_names=["A_ref", "omega_ref", "phi_ref", "lambda", "beta"]))
println(diagnostic_dashboard_text(drift_result))
