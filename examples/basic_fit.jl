using JuFitter
using LaTeXStrings
include(joinpath(@__DIR__, "_example_utils.jl"))

# Photoelectric-effect example with threshold kink:
#   U_0(nu) = max((h/e) * (nu - nu_0), 0)
# Above the threshold frequency nu_0 this is the Einstein equation
#   U_0 = (h/e) * nu - Phi/e,
# while below threshold no stopping voltage is observed.

const c = 299_792_458.0                   # m/s, exact
const elementary_charge = 1.602176634e-19 # C, exact
const h_reference = 6.62607015e-34         # J s, exact SI value

wavelength_nm = [150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0]
stopping_voltage_V = [5.97, 3.90, 2.66, 1.83, 1.24, 0.80, 0.46, 0.18, 0.0, 0.0]

# The table does not provide uncertainties. These are conservative readout
# assumptions for the example; replace them with the actual lab uncertainties if known.
sigma_wavelength_nm = fill(0.01, length(wavelength_nm))
sigma_voltage_V = fill(0.01, length(stopping_voltage_V))

frequency_THz = @. c / (wavelength_nm * 1e-9) / 1e12
sigma_frequency_THz = @. c * (sigma_wavelength_nm * 1e-9) / (wavelength_nm * 1e-9)^2 / 1e12

# p[1] = h/e in V/THz, p[2] = threshold frequency in THz.
photo_threshold_model(f_THz, p) = @. max(p[1] * (f_THz - p[2]), 0.0)

result = fit_model(
    photo_threshold_model,
    frequency_THz,
    stopping_voltage_V;
    p0=[0.004, 520.0],
    sigma_y=sigma_voltage_V,
    sigma_x=sigma_frequency_THz,
    bounds=([0.0, 0.0], [0.02, 2000.0]),
    initial_guesses=[[0.004, 500.0], [0.0042, 540.0], [0.0038, 450.0]],
    multistart=3,
)

slope = result.params[1]
threshold_frequency_THz = result.params[2]
sigma_slope = result.param_stderr[1]
sigma_threshold_frequency_THz = result.param_stderr[2]

h_fit = slope * elementary_charge / 1e12
sigma_h = sigma_slope * elementary_charge / 1e12
work_function_eV = slope * threshold_frequency_THz
sigma_work_function_eV = sqrt(
    (threshold_frequency_THz * sigma_slope)^2 +
    (slope * sigma_threshold_frequency_THz)^2 +
    2 * threshold_frequency_THz * slope * result.param_covariance[1, 2]
)
threshold_wavelength_nm = c / (threshold_frequency_THz * 1e12) * 1e9

println("Photoelectric-effect threshold fit")
println("slope h/e = ", slope, " ± ", sigma_slope, " V/THz")
println("threshold frequency = ", threshold_frequency_THz, " ± ", sigma_threshold_frequency_THz, " THz")
println("threshold wavelength = ", threshold_wavelength_nm, " nm")
println("h = ", h_fit, " ± ", sigma_h, " J s")
println("h reference = ", h_reference, " J s")
println("relative h deviation = ", (h_fit - h_reference) / h_reference)
println("work function Phi = ", work_function_eV, " ± ", sigma_work_function_eV, " eV")
println("chi2/ndf = ", result.stats.chi2_ndf)
println(report_text(result; parameter_names=["h/e", "nu_0"]))

plot_fit(
    result;
    filename=example_output("photoeffect_threshold_fit.pdf"),
    format=:pdf,
    theme=:latex,
    latex_labels=true,
    latex_stats=true,
    title=L"\text{Photoeffekt mit Schwellenknick: }U_0=\max\left(\frac{h}{e}(\nu-\nu_0),0\right)",
    xlabel=L"\nu\ \text{in THz}",
    ylabel=L"U_0\ \text{in V}",
    parameter_names=[L"h/e", L"\nu_0"],
    fit_label=L"\text{Schwellen-Fit}",
    data_label=L"\text{Messdaten}",
    band_label=L"1\sigma\ \text{Band}",
    fit_color=:darkblue,
    band_color=:dodgerblue,
    data_marker=:circle,
    figure_size=(1200, 720),
    stats_panel_width=340,
)

plot_diagnostics(
    result;
    filename=example_output("photoeffect_threshold_diagnostics.pdf"),
    format=:pdf,
    theme=:latex,
    xlabel=L"\nu\ \text{in THz}",
)

println("Saved plot to ", example_output("photoeffect_threshold_fit.pdf"))
println("Saved diagnostics to ", example_output("photoeffect_threshold_diagnostics.pdf"))
