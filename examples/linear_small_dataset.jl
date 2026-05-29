using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))

# Small dataset (12 points), heteroscedastic y-errors.
x = collect(range(-2.0, 3.5; length=12))
model(x, p) = @. p[1] * x + p[2]
true_p = [1.25, -0.35]

sigma_y = 0.08 .+ 0.02 .* abs.(x)
y = model(x, true_p) .+ sigma_y .* sin.(2.6 .* x)

result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y)

println("backend: ", result.backend)
println("params: ", result.params)
println("stderr: ", result.param_stderr)
println("chi2/ndf: ", result.stats.chi2_ndf)

xgrid = collect(range(minimum(x), maximum(x); length=300))
plot_fit(result; xgrid=xgrid, filename=example_output("linear_small"), format=:svg)
println("Saved plot to ", example_output("linear_small.svg"))
