using JuFitter
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

x = collect(range(0.15, 2.2; length=18))
model(t, p) = @. p[1] * (1 - exp(-t / p[2])) + p[3]
sigma_x = @. 0.010 + 0.004 * x
sigma_y = @. 0.045 + 0.008 * x
residual_pattern = [
    0.50, -0.90, 0.30, 1.10, -0.70, 0.80, -1.00, 0.40, 0.90,
    -0.60, 0.70, -0.80, 1.00, -0.40, 0.55, -0.75, 0.65, -0.35,
]
y = model(x, [4.8, 3.4, 0.12]) .+ sigma_y .* residual_pattern

fit = fitplot(
    model,
    x,
    y;
    p0=[3.0, 2.0, 0.0],
    sigma_y=sigma_y,
    sigma_x=sigma_x,
    bounds=([0.1, 0.1, -0.5], [20.0, 20.0, 1.0]),
    parameter_priors=(index=3, mean=0.10, sigma=0.08),
    initial_guesses=[[3.0, 2.0, 0.0], [8.0, 7.0, 0.1], [2.0, 1.0, 0.2]],
    filename=example_output("05_constraints_priors_fit.pdf"),
    title="Early saturation measurement",
    xlabel="time (s)",
    ylabel="response (V)",
    parameter_names=["amplitude", "time constant", "offset"],
    report=:both,
)

result = fit.result
prof = profile(result, 1; npoints=61, nsigma=4)
interval = profile_interval(result, 1; npoints=81, nsigma=4)
cont = contour(result, 1, 2; npoints=121, nsigma=4)

plot_profile(
    prof;
    filename=example_output("05_saturation_profile.pdf"),
    title="Profile cost versus local parabola",
    xlabel="amplitude A",
    local_sigma=result.param_stderr[1],
    delta_max=8,
)
plot_contour(
    cont;
    filename=example_output("05_amplitude_timescale_contour.pdf"),
    title="Profile contours versus local covariance",
    xlabel="amplitude A",
    ylabel="time constant tau",
    local_covariance=result.param_covariance,
    local_center=result.params[[1, 2]],
)

println("Profile interval for amplitude: -", interval.uncertainty_minus, " +", interval.uncertainty_plus)
println("Saved profile and contour plots.")
