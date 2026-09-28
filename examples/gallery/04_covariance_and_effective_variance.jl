using ScientificFitting
using LinearAlgebra
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

# Full y-covariance: correlated readout noise with correlation time 0.3 s,
# built from the time axis so off-diagonal terms follow |x_i - x_j|.
x = collect(range(0.0, 2.5; length=18))
model(x, p) = @. p[1] * exp(p[2] * x) + p[3]
n = length(x)
base_sigma = 0.05
correlation_time = 0.3  # s
cov_y = [base_sigma^2 * exp(-abs(x[i] - x[j]) / correlation_time) for i in 1:n, j in 1:n]
# One fixed noise realization drawn from N(0, cov_y): the Cholesky factor of
# cov_y applied to a standard-normal vector, pasted literally so the script
# is reproducible without a random number generator and the printed chi2/ndf
# is consistent with the stated covariance.
y = [2.24655, 1.93747, 1.6677, 1.49308, 1.28078, 1.09324, 0.906305, 0.852075,
     0.750124, 0.673371, 0.581348, 0.57533, 0.508985, 0.497965, 0.403454,
     0.366872, 0.310196, 0.265903]

cov_fit = fitplot(
    model,
    x,
    y;
    p0=[1.5, -0.7, 0.0],
    cov_y=cov_y,
    filename=example_output("04_full_covariance.pdf"),
    title="Exponential fit with full y-covariance",
    xlabel="time",
    xunit="s",
    ylabel="signal",
    parameter_names=["A", "lambda", "C"],
    show_panel=false,
    print_report=true,
)

# Effective variance: x-errors enter through the model derivative,
# sigma_eff^2 = sigma_y^2 + (m * sigma_x)^2. For a straight line with constant
# sigma_x this weight is the same at every point; the x-errors still act
# through the m-dependence of sigma_eff in the likelihood, including its
# log-determinant term.
x_true = collect(range(0.0, 4.0; length=16))
linear_model(x, p) = @. p[1] * x + p[2]
sigma_x = fill(0.16, length(x_true))
sigma_y = fill(0.10, length(x_true))
# The response is generated at the TRUE abscissa; only the recorded x is
# perturbed. The declared x uncertainty therefore really is in the data.
x_obs = x_true .+ sigma_x .* cos.(2.2 .* x_true)
y_obs = linear_model(x_true, [0.9, 1.2]) .+ sigma_y .* sin.(3.1 .* x_true)

xy_fit = fitplot(
    linear_model,
    x_obs,
    y_obs;
    p0=[0.5, 0.5],
    sigma_y=sigma_y,
    sigma_x=sigma_x,
    filename=example_output("04_effective_variance.pdf"),
    title="Linear fit with x and y uncertainties",
    xlabel="measured x",
    ylabel="measured y",
    parameter_names=["m", "b"],
    show_panel=false,
    print_report=true,
)

println("Full covariance backend: ", cov_fit.result.backend)
println("Effective-variance backend: ", xy_fit.result.backend)
