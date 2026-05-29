using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))

# Quadratic model with bounds + inequality constraint.
x = collect(range(-1.8, 2.2; length=15))
model(x, p) = @. p[1] * x^2 + p[2] * x + p[3]
true_p = [0.6, -0.8, 0.4]

sigma_y = fill(0.09, length(x))
y = model(x, true_p) .+ sigma_y .* cos.(2.0 .* x)

constraints = (
    ineq = p -> [-p[1]],  # enforce p1 >= 0 (convex parabola)
)

bounds = ([0.0, -2.0, -1.0], [2.0, 2.0, 2.0])

result = fit_model(
    model,
    x,
    y;
    p0=[0.2, -0.2, 0.0],
    sigma_y=sigma_y,
    bounds=bounds,
    constraints=constraints,
)

println("backend: ", result.backend)
println("params: ", result.params)
println("p1 >= 0 check: ", result.params[1] >= 0.0)

xgrid = collect(range(-2.0, 2.4; length=500))
plot_fit(result; xgrid=xgrid, filename=example_output("constrained_quadratic"), format=:png)
println("Saved plot to ", example_output("constrained_quadratic.png"))
