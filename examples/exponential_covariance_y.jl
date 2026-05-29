using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))
using LinearAlgebra

# Exponential model with full y-covariance (correlated errors).
x = collect(range(0.0, 2.5; length=16))
model(x, p) = @. p[1] * exp(p[2] * x) + p[3]
true_p = [2.0, -1.1, 0.25]

n = length(x)
base_sigma = 0.05
corr_len = 2.0
cov_y = Matrix{Float64}(undef, n, n)
for i in 1:n, j in 1:n
    cov_y[i, j] = base_sigma^2 * exp(-abs(i - j) / corr_len)
end

noise_shape = sin.(1.7 .* x) .+ 0.3 .* cos.(2.9 .* x)
# Scale deterministic noise to be modest relative to covariance diagonal.
y = model(x, true_p) .+ 0.6 .* base_sigma .* noise_shape

result = fit_model(model, x, y; p0=[1.5, -0.7, 0.0], cov_y=cov_y)

println("backend: ", result.backend)
println("params: ", result.params)
println("stderr: ", result.param_stderr)
println("AIC/BIC: ", (result.stats.aic, result.stats.bic))

plot_fit(result; filename=example_output("exp_covariance.pdf"))
println("Saved plot to ", example_output("exp_covariance.pdf"))
