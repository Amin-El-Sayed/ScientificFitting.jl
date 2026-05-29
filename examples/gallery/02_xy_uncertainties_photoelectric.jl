using JuFitter
using LaTeXStrings
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

# Photoelectric-effect threshold model:
# U_0(nu) = max((h/e) * (nu - nu_0), 0).
const c = 299_792_458.0
const elementary_charge = 1.602176634e-19
const h_reference = 6.62607015e-34

wavelength_nm = [150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0]
stopping_voltage_V = [5.97, 3.90, 2.66, 1.83, 1.24, 0.80, 0.46, 0.18, 0.0, 0.0]
sigma_wavelength_nm = fill(0.01, length(wavelength_nm))
sigma_voltage_V = fill(0.01, length(stopping_voltage_V))

frequency_THz = @. c / (wavelength_nm * 1e-9) / 1e12
sigma_frequency_THz = @. c * (sigma_wavelength_nm * 1e-9) / (wavelength_nm * 1e-9)^2 / 1e12

photo_threshold_model(f_THz, p) = @. max(p[1] * (f_THz - p[2]), 0.0)

fit = fitplot(
    photo_threshold_model,
    frequency_THz,
    stopping_voltage_V;
    p0=[0.004, 520.0],
    sigma_y=sigma_voltage_V,
    sigma_x=sigma_frequency_THz,
    bounds=([0.0, 0.0], [0.02, 2000.0]),
    initial_guesses=[[0.004, 500.0], [0.0042, 540.0], [0.0038, 450.0]],
    multistart=3,
    filename=example_output("02_photoelectric_fit.pdf"),
    theme=:latex,
    latex_labels=true,
    latex_stats=true,
    title=L"\text{Photoelectric threshold fit}",
    xlabel=L"\nu",
    xunit=L"\mathrm{THz}",
    ylabel=L"U_0",
    yunit=L"\mathrm{V}",
    parameter_names=[L"h/e", L"\nu_0"],
    fit_label=L"\text{threshold model}",
    data_label=L"\text{data}",
    band_label=L"1\sigma",
    report=:console,
)

result = fit.result
slope = result.params[1]
sigma_slope = result.param_stderr[1]
h_fit = slope * elementary_charge / 1e12
sigma_h = sigma_slope * elementary_charge / 1e12

plot_diagnostics(
    result;
    filename=example_output("02_photoelectric_diagnostics.pdf"),
    theme=:latex,
    xlabel=L"\nu\ \mathrm{in\ THz}",
)

println("h = ", h_fit, " +/- ", sigma_h, " J s")
println("relative h deviation = ", (h_fit - h_reference) / h_reference)
println("Saved plots to ", example_output("02_photoelectric_fit.pdf"), " and ", example_output("02_photoelectric_diagnostics.pdf"))
