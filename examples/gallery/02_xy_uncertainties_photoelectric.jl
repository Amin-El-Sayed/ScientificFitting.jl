using CairoMakie
using JuFitter
using LinearAlgebra
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

# Controlled teaching data: two experimentally distinguishable linear regimes,
# with individual x/y standard uncertainties at every point.
const elementary_charge = 1.602176634e-19

frequency_THz = [350.0, 380.0, 410.0, 440.0, 470.0, 495.0, 515.0, 532.0,
                 565.0, 590.0, 620.0, 655.0, 690.0, 730.0, 775.0, 825.0,
                 880.0, 940.0]
voltage_V = [0.0312, -0.0434, 0.01855, 0.0594, -0.02685, 0.04495, -0.0057,
             0.05784, 0.12324, 0.13123, 0.34230, 0.52185, 0.57404, 0.85602,
             0.94333, 1.23891, 1.35237, 1.70871]
sigma_frequency_THz = [4.5, 4.2, 4.0, 3.8, 3.6, 3.4, 3.2, 3.0,
                       2.9, 2.8, 2.7, 2.6, 2.5, 2.5, 2.4, 2.4, 2.3, 2.3]
sigma_voltage_V = [0.038, 0.040, 0.041, 0.043, 0.045, 0.047, 0.050, 0.052,
                   0.048, 0.050, 0.052, 0.054, 0.057, 0.060, 0.064, 0.068,
                   0.073, 0.080]

baseline_mask = frequency_THz .<= 532.0
emission_mask = .!baseline_mask
reference_frequency_THz = 550.0
line_model(x, p) = @. p[1] * (x - reference_frequency_THz) + p[2]

baseline = fit_model(
    line_model,
    frequency_THz[baseline_mask],
    voltage_V[baseline_mask];
    p0=[0.0, 0.02],
    sigma_x=sigma_frequency_THz[baseline_mask],
    sigma_y=sigma_voltage_V[baseline_mask],
)
emission = fit_model(
    line_model,
    frequency_THz[emission_mask],
    voltage_V[emission_mask];
    p0=[0.0042, 0.02],
    sigma_x=sigma_frequency_THz[emission_mask],
    sigma_y=sigma_voltage_V[emission_mask],
    bounds=([0.0, -5.0], [0.02, 5.0]),
)

me, ce = emission.params
mb, cb = baseline.params
photoelectric_slope = me - mb
threshold_offset_THz = (cb - ce) / photoelectric_slope
threshold_THz = reference_frequency_THz + threshold_offset_THz
gradient_emission = [
    -threshold_offset_THz / photoelectric_slope,
    -1 / photoelectric_slope,
]
gradient_baseline = [
    threshold_offset_THz / photoelectric_slope,
    1 / photoelectric_slope,
]
threshold_variance =
    dot(gradient_emission, emission.param_covariance * gradient_emission) +
    dot(gradient_baseline, baseline.param_covariance * gradient_baseline)
sigma_threshold_THz = sqrt(threshold_variance)

h_fit = photoelectric_slope * elementary_charge / 1e12
sigma_photoelectric_slope = sqrt(
    emission.param_covariance[1, 1] + baseline.param_covariance[1, 1],
)
sigma_h = sigma_photoelectric_slope * elementary_charge / 1e12

work_function_eV = photoelectric_slope * threshold_THz
work_gradient_emission = [reference_frequency_THz, -1.0]
work_gradient_baseline = [-reference_frequency_THz, 1.0]
work_variance =
    dot(work_gradient_emission, emission.param_covariance * work_gradient_emission) +
    dot(work_gradient_baseline, baseline.param_covariance * work_gradient_baseline)
sigma_work_function_eV = sqrt(work_variance)

# JuFitter returns Makie objects, so derived quantities and annotations remain
# ordinary Makie operations rather than a special plotting mini-language.
style = :modern
appearance = :light
palette = plot_palette(style; appearance=appearance)
baseline_color = palette.stats_muted_color
threshold_color = palette.stats_color

fig = with_theme(plot_theme(style; appearance=appearance)) do
    Figure(size=(1120, 700))
end
ax = Axis(
    fig[1, 1];
    title="Photoelectric threshold from two fitted regimes",
    xlabel="frequency ν (THz)",
    ylabel="stopping voltage U₀ (V)",
)
errorbars!(
    ax,
    frequency_THz,
    voltage_V,
    sigma_voltage_V;
    color=palette.yerr_color,
    whiskerwidth=palette.error_whiskerwidth,
)
errorbars!(
    ax,
    frequency_THz,
    voltage_V,
    sigma_frequency_THz;
    direction=:x,
    color=palette.xerr_color,
    whiskerwidth=palette.error_whiskerwidth,
)
scatter!(ax, frequency_THz[baseline_mask], voltage_V[baseline_mask]; color=baseline_color, marker=:diamond, label="baseline")
scatter!(ax, frequency_THz[emission_mask], voltage_V[emission_mask]; color=palette.data_color, label="emission")

xgrid = collect(range(minimum(frequency_THz) - 15, maximum(frequency_THz) + 15; length=500))
J = hcat(xgrid .- reference_frequency_THz, ones(length(xgrid)))
for (result, color, label) in ((baseline, baseline_color, "baseline fit"), (emission, palette.fit_color, "emission fit"))
    ygrid = line_model(xgrid, result.params)
    sigma_fit = sqrt.(clamp.(vec(sum((J * result.param_covariance) .* J; dims=2)), 0.0, Inf))
    band!(ax, xgrid, ygrid .- sigma_fit, ygrid .+ sigma_fit; color=(color, 0.22), label="$label 1σ")
    lines!(ax, xgrid, ygrid; color=color, linewidth=palette.fit_linewidth, label=label)
end
add_vband!(
    ax,
    threshold_THz - sigma_threshold_THz,
    threshold_THz + sigma_threshold_THz;
    color=(threshold_color, 0.14),
    label="threshold 1σ",
)
add_vline!(ax, threshold_THz; color=threshold_color, linestyle=:dash, linewidth=2)
add_points!(
    ax,
    [threshold_THz],
    [line_model([threshold_THz], emission.params)[1]];
    marker=:star5,
    markersize=18,
    color=threshold_color,
    label="intersection",
)
plot_info_panel!(
    fig[1, 2];
    legend_source=ax,
    model_label="ΔU(ν) = mγ (ν - ν₀)",
    parameter_lines=[
        "photoelectric slope = $(round(photoelectric_slope; sigdigits=5)) V/THz",
        "h = $(round(h_fit; sigdigits=4)) ± $(round(sigma_h; sigdigits=2)) J s",
        "baseline slope = $(round(mb; sigdigits=4)) V/THz",
        "ν0 = $(round(threshold_THz; sigdigits=5)) ± $(round(sigma_threshold_THz; sigdigits=2)) THz",
        "Φ = $(round(work_function_eV; sigdigits=5)) ± $(round(sigma_work_function_eV; sigdigits=2)) eV",
    ],
    statistic_lines=[
        "emission χ²/ndf = $(round(emission.stats.chi2_ndf; sigdigits=4))",
        "baseline χ²/ndf = $(round(baseline.stats.chi2_ndf; sigdigits=4))",
    ],
    color=palette.stats_color,
    muted_color=palette.stats_muted_color,
)
save(example_output("02_photoelectric_fit.pdf"), fig)

println("h = ", h_fit, " +/- ", sigma_h, " J s")
println("threshold = ", threshold_THz, " +/- ", sigma_threshold_THz, " THz")
println("work function = ", work_function_eV, " +/- ", sigma_work_function_eV, " eV")
println("baseline")
println(diagnostic_dashboard_text(baseline))
println("emission")
println(diagnostic_dashboard_text(emission))
println("Saved plot to ", example_output("02_photoelectric_fit.pdf"))
