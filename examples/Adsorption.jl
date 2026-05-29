using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))
using CairoMakie
using LaTeXStrings

#Aufbau A wird verwendet

T = 19.5
sigma_T = 0.5
pr = [0.0, 4.0, 7.0, 14.0, 29.0, 72.0, 123.0, 168.0, 302.0, 398.0, 507.0, 586.0, 636.0, 684.0]
V = [0.015, 0.21, 0.4, 0.61, 0.81, 1.02, 1.12, 1.21, 1.32, 1.43, 1.57, 1.67, 1.77, 1.9] .- 0.015
conversion_factor = (0.641 - 0.000992 * 25.0) / 72.15
n_ads = V .* conversion_factor

sigma_pr = fill(1.0, length(pr))
sigma_n_ads = fill(0.02 * conversion_factor, length(pr))
langmuir_model(pr, p) = @. p[1] * (p[2] * pr) / (1 + (p[2] * pr))
BET_model(pr, p) = @. p[1] * (exp((-p[2] - 26800) / (8.314 * (p[4] + 273.15))) * (pr / p[3])) / ((1 - (pr / p[3])) * (1 + (exp((-p[2] - 26800) / (8.314 * (p[4] + 273.15))) - 1) * (pr / p[3])))


# Initial guesses
p0_guess_langmuir = [6.0, 0.02]
bounds_langmuir = ([0.0, 0.0], [100.0, 1.0])

p0_guess_BET = [0.01, -38000.0, 1000.0, T]
# Keep p0 strictly above measured pressures to avoid the BET singularity at p/p0 = 1.
bounds_BET = ([0.0, -1.0e5, maximum(pr) + 1.0, 0.0], [1.0, -1000.0, 5000.0, 100.0])

T_prior = (index=4, mean=T, sigma=sigma_T)

# Includes BOTH y and x uncertainties
result_langmuir = fit_model(langmuir_model, pr, n_ads; p0=p0_guess_langmuir, sigma_y=sigma_n_ads, sigma_x=sigma_pr, bounds=bounds_langmuir)

result_BET = fit_model(BET_model, pr, n_ads; p0=p0_guess_BET, bounds=bounds_BET, sigma_y=sigma_n_ads, sigma_x=sigma_pr, parameter_priors=[T_prior])



latex_theme = Theme(
    font="CMU Serif",
    fontsize=22,
    Legend=(
        labelsize=20,
        framevisible=true,
        backgroundcolor=:white,
        patchsize=(28, 20),
    ),
)

plot_fit(
    result_langmuir;
    filename=example_output("adsorption_plot_langmuir"),
    format=:svg,
    theme=:latex,
    theme_override=latex_theme,
    plot_aspect=2.0,
    figure_size=(1300, 500),
    stats_panel_width=0.08,
    panel_gap=6,
    latex_labels=true,
    latex_stats=true,
    title=L"\text{Langmuir-Modell: }n_\mathrm{mono}(p)=n_\mathrm{max}\frac{Kp}{1+Kp}",
    xlabel=L"\text{Druck }p \text{ in mbar}",
    ylabel=L"\text{Adsorbierte Stoffmenge }n_\mathrm{ads} \text{ in mol}",
    parameter_names=[L"n_\mathrm{mono}", L"K"],
    fit_label=L"\text{Fit}",
    band_label=L"1\sigma\ \text{Band}",
    data_label=L"\text{Daten}",
    fit_color=:darkblue,
    fit_linewidth=1.0,
    band_color=:dodgerblue,
    band_alpha=0.22,
    data_marker=:circle,
    data_markersize=8,
    show_legend=true,
    legend_position=:lt,
)

plot_fit(
    result_BET;
    filename=example_output("adsorption_plot_BET"),
    format=:svg,
    theme=:latex,
    theme_override=latex_theme,
    plot_aspect=2.0,
    figure_size=(1300, 650),
    stats_panel_width=0.08,
    panel_gap=6,
    latex_labels=true,
    latex_stats=true,
    title=L"\text{BET-Modell: }n_\mathrm{ads}(S)=n_\mathrm{mono}\frac{exp\left(\frac{-\Delta_{\text{ads}}\text{H}-\Delta_{\text{vap}}\text{H}}{\text{RT}}\right)\cdot\frac{p}{p_0}}{(1-\frac{p}{p_0})(1+(exp\left(\frac{-\Delta_{\text{ads}}\text{H}-\Delta_{\text{vap}}\text{H}}{\text{RT}}\right)-1)\frac{p}{p_0})}",
    xlabel=L"\text{Druck }p \text{ in mbar}",
    ylabel=L"\text{Adsorbierte Stoffmenge }n_\mathrm{ads} \text{ in mol}",
    parameter_names=[L"n_\mathrm{mono}", L"\Delta_{\text{ads}}\text{H}", L"p_0", L"\text{T}"],
    fit_label=L"\text{Fit}",
    band_label=L"1\sigma\ \text{Band}",
    data_label=L"\text{Daten}",
    fit_color=:darkred,
    fit_linewidth=1.0,
    band_color=:orangered,
    band_alpha=0.22,
    data_marker=:circle,
    data_markersize=8,
    show_legend=true,
    legend_position=:lt,
)

println("Fitted params = ", result_BET.params)
println("Parameter uncertainties = ", result_BET.param_stderr)
println("backend = ", result_BET.backend)
println("chi2 = ", result_BET.stats.chi2)
println("chi2/ndf = ", result_BET.stats.chi2_ndf)
println("Saved plot to ", example_output("adsorption_plot_BET.svg"))
