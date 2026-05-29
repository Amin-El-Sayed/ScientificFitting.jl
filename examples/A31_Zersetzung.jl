using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))
using CairoMakie
using LaTeXStrings

T = 292.15
sigma_T = 1.0

h_0 = 8.0
sigma_h_0 = 0.1

#t = [4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0]
#h = [9.0, 10.1, 11.2, 12.1, 13.8, 14.7, 15.4, 16.1, 16.8, 17.45]
#h_prime = [18.1, 18.65, 19.2, 19.7, 20.2, 20.65, 21.0, 21.7, 22.2, 22.75]
t = [4.0, 8.0, 12.0, 16.0, 24.0, 28.0, 32.0, 36.0, 40.0]
h = [9.0, 10.1, 11.2, 12.1, 13.8, 14.7, 15.4, 16.1, 16.8]
h_prime = [17.45, 18.1, 18.65, 19.2, 19.7, 20.2, 20.65, 21.10, 21.4]
y = log.(h_prime .- h)

sigma_x = fill(1 / 12, length(t))
sigma_h = fill(0.1, length(t))
sigma_h_prime = fill(0.1, length(t))
sigma_y = @. sqrt(2) * sigma_h / (h_prime - h) # Propagation of uncertainty for log(h'-h)
model(t, p) = @. (-p[1] * t) + log(abs((p[2]-p[3])*(1-exp(-p[1]*p[4]))))

p0_guess = [6.0, 12.0, 3.4, 40.0] # [k, h_infty, h_0, tau]
bounds = ([0.0, 0.0, 0.0, 0.0], [1000.0, 1000.0, 1000.0, 50.0])
tau_prior = (index=4, mean=40, sigma=1/6)
h_0_prior = (index=3, mean=h_0, sigma=sigma_h_0)

# Includes BOTH y and x uncertainties
result = fit_model(model, t, y; p0=p0_guess, sigma_y=sigma_y, sigma_x=sigma_x, bounds=bounds, parameter_priors=[tau_prior, h_0_prior])

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
    filename=example_output("Zersetzung_langsam_plot"),
    format=:svg,
    theme=:latex,
    theme_override=latex_theme,
    plot_aspect=1.5,
    figure_size=(1200, 500),
    stats_panel_width=0.06,
    panel_gap=6,
    latex_labels=true,
    latex_stats=true,
    title=L"\text{model function: }ln(h^' -h)=-kt+ln((h_{\infty}-h_0)-(1-exp(-k\tau)))", #-\alpha_{\text{offset}}
    xlabel=L"\text{Zeit }t\ \text{in min}",
    ylabel=L"ln(h^' -h)",
    parameter_names=[L"k", L"h_{\infty}", L"h_0", L"\tau"],
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
println("Saved plot to ", example_output("Zersetzung_langsam_plot.svg"))
