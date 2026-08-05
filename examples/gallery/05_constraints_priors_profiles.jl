using CairoMakie
using JuFitter
using Printf
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

x = [
    0.150000, 0.270588, 0.391176, 0.511765, 0.632353, 0.752941,
    0.873529, 0.994118, 1.114706, 1.235294, 1.355882, 1.476471,
    1.597059, 1.717647, 1.838235, 1.958824, 2.079412, 2.200000,
]
y = [
    0.345641, 0.470693, 0.593535, 0.839840, 0.969673, 1.129630,
    1.165946, 1.315719, 1.477930, 1.631653, 1.729280, 1.759687,
    1.878716, 1.988480, 2.172411, 2.176549, 2.346881, 2.447489,
]
sigma_x = [
    0.010600, 0.011082, 0.011565, 0.012047, 0.012529, 0.013012,
    0.013494, 0.013976, 0.014459, 0.014941, 0.015424, 0.015906,
    0.016388, 0.016871, 0.017353, 0.017835, 0.018318, 0.018800,
]
sigma_y = [
    0.046200, 0.047165, 0.048129, 0.049094, 0.050059, 0.051024,
    0.051988, 0.052953, 0.053918, 0.054882, 0.055847, 0.056812,
    0.057776, 0.058741, 0.059706, 0.060671, 0.061635, 0.062600,
]
model(t, p) = @. p[1] * (1 - exp(-t / p[2])) + p[3]

fit = fitplot(
    model,
    x,
    y;
    p0=[4.5, 3.0, 0.1],
    sigma_y=sigma_y,
    sigma_x=sigma_x,
    bounds=([0.1, 0.1, -0.5], [20.0, 20.0, 1.0]),
    parameter_priors=(index=3, mean=0.10, sigma=0.08),
    initial_guesses=[[6.0, 5.0, 0.1], [3.0, 2.0, 0.1]],
    maxiters=2000,
    filename=example_output("05_constraints_priors_fit.pdf"),
    title="Early saturation measurement",
    xlabel="time (s)",
    ylabel="response (V)",
    parameter_names=["amplitude", "time constant", "offset"],
    report=:both,
)

result = fit.result
interval = profile_interval(result, 1; npoints=81, nsigma=4)
prof = interval.profile_result
cont = JuFitter.contour(result, 1, 2; npoints=121, nsigma=4)

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
    title="Profile versus local covariance",
    xlabel="amplitude A",
    ylabel="time constant tau",
    local_covariance=result.param_covariance,
    local_center=result.params[[1, 2]],
)

@printf(
    "Profile interval for amplitude: %.3f -%.3f +%.3f V\n",
    result.params[1],
    interval.uncertainty_minus,
    interval.uncertainty_plus,
)
println("Saved profile and contour plots.")
