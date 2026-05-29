using JuFitter
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

x = collect(range(-2.0, 2.3; length=28))
model(x, p) = @. p[1] * x^2 + p[2] * x + p[3]
sigma_y = fill(0.08, length(x))
y = model(x, [0.65, -0.75, 0.35]) .+ sigma_y .* cos.(2.0 .* x)

constraints = (
    ineq = p -> [-p[1]],  # p[1] >= 0: convex parabola
)

fit = fitplot(
    model,
    x,
    y;
    p0=[0.25, -0.2, 0.0],
    sigma_y=sigma_y,
    bounds=([0.0, -2.0, -1.0], [2.0, 2.0, 2.0]),
    constraints=constraints,
    parameter_priors=(index=3, mean=0.3, sigma=0.2),
    filename=example_output("05_constraints_priors_fit.pdf"),
    title="Constrained quadratic fit",
    xlabel="position",
    ylabel="response",
    parameter_names=["curvature", "slope", "offset"],
    report=:both,
)

result = fit.result
prof = profile(result, 1; npoints=41, nsigma=3)
interval = profile_interval(result, 1; npoints=81, nsigma=4)
cont = contour(result, 1, 2; npoints=25, nsigma=2)

plot_profile(
    prof;
    filename=example_output("05_curvature_profile.pdf"),
    xlabel="curvature",
)
plot_contour(
    cont;
    filename=example_output("05_curvature_slope_contour.pdf"),
    xlabel="curvature",
    ylabel="slope",
)

println("Profile interval for curvature: -", interval.uncertainty_minus, " +", interval.uncertainty_plus)
println("Saved profile and contour plots.")
