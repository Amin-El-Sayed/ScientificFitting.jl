using CairoMakie
using ScientificFitting
using LaTeXStrings
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

x = collect(range(-1.5, 2.0; length=16))
model(x, p) = @. p[1] * x^2 + p[2] * x + p[3]
sigma_y = 0.07 .+ 0.025 .* abs.(x)
y = model(x, [0.72, -0.38, 0.2]) .+ sigma_y .* sin.(3.4 .* x)

# Theme overrides should change only the intended Makie attributes. Typography,
# layout, and the remaining marks continue to follow ScientificFitting's style contract.
custom_theme = Theme(Axis=(xgridvisible=false, ygridvisible=true))

fit = fitplot(
    model,
    x,
    y;
    p0=[0.3, 0.0, 0.0],
    sigma_y=sigma_y,
    filename=example_output("03_plot_customization.svg"),
    theme=:sans,
    theme_override=custom_theme,
    title="Quadratic calibration",
    model_label=L"y = a x^2 + b x + c",
    xlabel="control variable",
    xunit="a.u.",
    ylabel="response",
    yunit="mV",
    parameter_names=["a", "b", "c"],
    nsigma=2,
    band_label="2-sigma band",
    fit_color="#0072b2",
    band_color="#0072b2",
    band_alpha=0.18,
    data_marker=:diamond,
    data_markersize=9,
    stats_sigdigits=6,
    show_panel=true,
    print_report=true,
    axis_kwargs=(xtickalign=1, ytickalign=1),
)

println("Saved plot to ", example_output("03_plot_customization.svg"))
println("Figure object type: ", typeof(fit.figure))
