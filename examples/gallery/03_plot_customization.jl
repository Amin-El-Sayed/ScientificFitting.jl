using CairoMakie
using JuFitter
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

x = collect(range(-1.5, 2.0; length=16))
model(x, p) = @. p[1] * x^2 + p[2] * x + p[3]
sigma_y = 0.07 .+ 0.025 .* abs.(x)
y = model(x, [0.72, -0.38, 0.2]) .+ sigma_y .* sin.(3.4 .* x)

custom_theme = Theme(
    fontsize=21,
    Axis=(xgridvisible=true, ygridvisible=true),
    Lines=(linewidth=4,),
    Scatter=(markersize=13,),
)

fit = fitplot(
    model,
    x,
    y;
    p0=[0.3, 0.0, 0.0],
    sigma_y=sigma_y,
    filename=example_output("03_plot_customization.svg"),
    theme=:clean,
    theme_override=custom_theme,
    title="Quadratic calibration",
    xlabel="control variable",
    xunit="a.u.",
    ylabel="response",
    yunit="mV",
    parameter_names=["a", "b", "c"],
    nsigma=2,
    band_label="2-sigma band",
    fit_color=:navy,
    band_color=:skyblue,
    band_alpha=0.28,
    data_marker=:star5,
    data_markersize=14,
    yerr_color=:gray30,
    legend_position=:lt,
    stats_sigdigits=6,
    report=:both,
    axis_kwargs=(xtickalign=1, ytickalign=1),
)

println("Saved plot to ", example_output("03_plot_customization.svg"))
println("Figure object type: ", typeof(fit.figure))
