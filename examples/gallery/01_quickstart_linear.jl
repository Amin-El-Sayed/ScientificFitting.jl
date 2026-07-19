using JuFitter
using CairoMakie
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

# Minimal workflow: x/y arrays plus y-uncertainties. Without an explicit model,
# fitplot uses a straight line with a robust initial guess from the endpoints.
x = [0.0, 0.4348, 0.8696, 1.3043, 1.7391, 2.1739, 2.6087, 3.0435,
     3.4783, 3.9130, 4.3478, 4.7826, 5.2174, 5.6522, 6.0870, 6.5217,
     6.9565, 7.3913, 7.8261, 8.2609, 8.6957, 9.1304, 9.5652, 10.0]
y = [0.7000, 1.6125, 2.4832, 3.2749, 3.9858, 4.6545, 5.3439, 6.1123,
     6.9838, 7.9338, 8.8975, 9.7983, 10.5849, 11.2581, 11.8738,
     12.5193, 13.2732, 14.1663, 15.1641, 16.1797, 17.1125, 17.8965,
     18.5336, 19.0964]
sigma_y = [0.1600, 0.1687, 0.1774, 0.1861, 0.1948, 0.2035, 0.2122,
           0.2209, 0.2296, 0.2383, 0.2470, 0.2557, 0.2643, 0.2730,
           0.2817, 0.2904, 0.2991, 0.3078, 0.3165, 0.3252, 0.3339,
           0.3426, 0.3513, 0.3600]

fit = fitplot(
    x,
    y;
    sigma_y=sigma_y,
    title="Quickstart calibration",
    model_label="U(x) = m x + b",
    xlabel="x",
    xunit="mm",
    ylabel="U",
    yunit="V",
    parameter_names=["m", "b"],
    band=:prediction,
    nsigma=1,
    band_label="1σ prediction band",
    show_legend=true,
    report=:both,
    filename=example_output("01_quickstart_linear.pdf"),
)

println()
println(diagnostic_dashboard_text(fit.result))
