using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))
using CairoMakie

x = collect(range(-1.5, 2.0; length=13))
model(x, p) = @. p[1] * x^2 + p[2] * x + p[3]
true_p = [0.7, -0.4, 0.2]

sigma_y = 0.07 .+ 0.03 .* abs.(x)
y = model(x, true_p) .+ sigma_y .* sin.(3.4 .* x)

result = fit_model(model, x, y; p0=[0.3, 0.0, 0.0], sigma_y=sigma_y)

# `theme_override` can be used to tweak fonts and global style.
custom_theme = Theme(
    fontsize=21,
    Lines=(linewidth=4,),
    Scatter=(markersize=14,),
)

plot_fit(
    result;
    filename=example_output("styled_plot_with_summary"),
    format=:pdf,
    theme=:publication,
    theme_override=custom_theme,
    title="Quadratic Calibration Fit",
    xlabel="control variable u",
    ylabel="response r",
    parameter_names=["a", "b", "c"],
    fit_color=:navy,
    band_color=:skyblue,
    band_alpha=0.30,
    data_color=:black,
    data_marker=:star5,
    data_markersize=14,
    yerr_color=:gray30,
    xerr_color=:gray30,
    show_legend=true,
    legend_position=:lt,
    stats_fontsize=18,
    stats_sigdigits=6,
)

println("Saved plot to ", example_output("styled_plot_with_summary.pdf"))
