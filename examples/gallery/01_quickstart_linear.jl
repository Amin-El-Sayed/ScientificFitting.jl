using JuFitter
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

# Minimal workflow: x/y arrays plus y-uncertainties. Without an explicit model,
# fitplot uses a straight line with a robust initial guess from the endpoints.
x = collect(range(0.0, 10.0; length=24))
sigma_y = 0.16 .+ 0.02 .* x
y = @. 1.85 * x + 0.7 + sigma_y * sin(1.6 * x)

fit = fitplot(
    x,
    y;
    sigma_y=sigma_y,
    filename=example_output("01_quickstart_linear.pdf"),
    xlabel="time",
    xunit="s",
    ylabel="signal",
    yunit="V",
    parameter_names=["slope", "offset"],
    report=:both,
)

println("Saved plot to ", example_output("01_quickstart_linear.pdf"))
println("Fit parameters: ", fit.result.params)
