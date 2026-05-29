using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))
using CairoMakie
using LaTeXStrings

x = collect(range(0.0, 3.0; length=18))
model(x, p) = @. p[1] * exp(-p[2] * x) + p[3]
true_p = [2.4, 1.1, 0.2]

sigma_y = fill(0.06, length(x))
y = model(x, true_p) .+ sigma_y .* sin.(2.8 .* x)

result = fit_model(model, x, y; p0=[2.0, 0.8, 0.0], sigma_y=sigma_y)

latex_theme = Theme(
    font="CMU Serif",
    fontsize=22,
)

plot_fit(
    result;
    filename=example_output("latex_styled_plot"),
    format=:pdf,
    theme=:latex,
    theme_override=latex_theme,
    latex_labels=true,
    latex_stats=true,
    title=L"\text{Exponential Fit: } y(x)=A e^{-\lambda x}+C",
    xlabel=L"x\,\,(\mathrm{s})",
    ylabel=L"y\,\,(\mathrm{a.u.})",
    parameter_names=[L"A", L"\lambda", L"C"],
    fit_color=:darkred,
    band_color=:red,
    band_alpha=0.22,
    data_marker=:circle,
    data_markersize=11,
    show_legend=true,
    legend_position=:rt,
)

println("Saved plot to ", example_output("latex_styled_plot.pdf"))
