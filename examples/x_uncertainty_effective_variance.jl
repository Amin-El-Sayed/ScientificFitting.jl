using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))

# Effective-variance example: compare fit with and without x-errors.
x_true = collect(range(0.0, 4.0; length=14))
model(x, p) = @. p[1] * x + p[2]
true_p = [0.9, 1.2]

sigma_x = fill(0.16, length(x_true))
sigma_y = fill(0.10, length(x_true))

x_obs = x_true .+ sigma_x .* cos.(2.2 .* x_true)
y_obs = model(x_obs, true_p) .+ sigma_y .* sin.(3.1 .* x_obs)

result_yonly = fit_model(model, x_obs, y_obs; p0=[0.5, 0.5], sigma_y=sigma_y)
result_xy = fit_model(model, x_obs, y_obs; p0=[0.5, 0.5], sigma_y=sigma_y, sigma_x=sigma_x)

println("backend (y-only): ", result_yonly.backend)
println("backend (x+y): ", result_xy.backend)
println("params y-only: ", result_yonly.params)
println("params x+y:    ", result_xy.params)

plot_fit(result_xy; filename=example_output("x_uncertainty.png"))
println("Saved plot to ", example_output("x_uncertainty.png"))
