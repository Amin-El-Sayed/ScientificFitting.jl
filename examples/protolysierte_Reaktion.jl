using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))
using CairoMakie
using LaTeXStrings

x = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0]
y_1 = [23.00, 22.75, 22.45, 22.40, 22.20, 22.00, 21.75, 21.55, 21.25, 21.10, 20.90, 20.50, 19.95, 19.95, 19.05, 18.80, 18.25, 17.85, 17.30]
y_2 = [22.95, 22.6, 22.40, 22.30, 22.15, 21.9, 21.65, 21.35, 21.15, 21.0, 20.85, 20.35, 19.75, 19.25, 18.85, 18.6, 18.1, 17.6, 17.15]
y = (y_1 .+ y_2) ./ 2

sigma_x = fill(1 / 6, length(x))
sigma_y_1 = fill(0.1, length(y_1))
sigma_y_2 = fill(0.1, length(y_2))
sigma_y = sqrt.(sigma_y_1.^2 .+ sigma_y_2.^2)
model(x, p) = @. (p[1]-p[3]) * exp(-p[2] * x) + p[3]# - p[4]

# Initial guess (p0) for [A, k, C]
p0_guess = [6.0, 0.02, 16.0, 11.95]
bounds = ([0.0, 0.0, 0.0, 0.0], [100.0, 1.0, 50.0, 50.0])
infinity_prior = (index=3, mean=10.575, sigma=0.1414)
offset_prior = (index=4, mean=11.95, sigma=0.1414)

# Includes BOTH y and x uncertainties
result = fit_model(model, x, y; p0=p0_guess, sigma_y=sigma_y, sigma_x=sigma_x, bounds=bounds, parameter_priors=[infinity_prior, offset_prior])

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
    result;
    filename=example_output("protolysierte_Reaktion_plot"),
    format=:svg,
    theme=:latex,
    theme_override=latex_theme,
    plot_aspect=1.5,
    figure_size=(1300, 500),
    stats_panel_width=0.06,
    panel_gap=6,
    latex_labels=true,
    latex_stats=true,
    title=L"\text{Modellfunktion: }\alpha(t)=(\alpha_0-\alpha_{\infty}) e^{-k t}+\alpha_{\infty}", #-\alpha_{\text{offset}}
    xlabel=L"\text{Zeit }t\ \text{in min}",
    ylabel=L"\text{Drehwinkel }\alpha\ \text{in }^\circ",
    parameter_names=[L"\alpha_0", L"k", L"\alpha_{\infty}", L"\alpha_{\text{offset}}"],
    fit_label=L"\text{Fit}",
    band_label=L"1\sigma\ \text{Band}",
    data_label=L"\text{Daten}",
    fit_color=:darkblue,
    fit_linewidth=1.5,
    band_color=:dodgerblue,
    band_alpha=0.22,
    data_marker=:circle,
    data_markersize=8,
    show_legend=true,
    legend_position=:rt,
)

println("Fitted params = ", result.params)
println("Parameter uncertainties = ", result.param_stderr)
println("backend = ", result.backend)
println("chi2 = ", result.stats.chi2)
println("chi2/ndf = ", result.stats.chi2_ndf)
println("Saved plot to ", example_output("protolysierte_Reaktion_plot.svg"))
