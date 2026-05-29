using JuFitter
using LinearAlgebra
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

# Full y-covariance: correlated readout noise.
x = collect(range(0.0, 2.5; length=18))
model(x, p) = @. p[1] * exp(p[2] * x) + p[3]
n = length(x)
base_sigma = 0.05
corr_len = 2.0
cov_y = [base_sigma^2 * exp(-abs(i - j) / corr_len) for i in 1:n, j in 1:n]
y = model(x, [2.0, -1.1, 0.25]) .+ 0.6 .* base_sigma .* (sin.(1.7 .* x) .+ 0.3 .* cos.(2.9 .* x))

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
    report=:console,
)

# Effective variance: x-errors contribute through the local model derivative.
x_true = collect(range(0.0, 4.0; length=16))
linear_model(x, p) = @. p[1] * x + p[2]
sigma_x = fill(0.16, length(x_true))
sigma_y = fill(0.10, length(x_true))
x_obs = x_true .+ sigma_x .* cos.(2.2 .* x_true)
y_obs = linear_model(x_obs, [0.9, 1.2]) .+ sigma_y .* sin.(3.1 .* x_obs)

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
    report=:console,
)

println("Full covariance backend: ", cov_fit.result.backend)
println("Effective-variance backend: ", xy_fit.result.backend)
