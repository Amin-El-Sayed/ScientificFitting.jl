using JuFitter
include(joinpath(@__DIR__, "_example_utils.jl"))

# Nonlinear model with user-supplied analytic Jacobian.
x = collect(range(0.0, 1.8; length=18))
model(x, p) = @. p[1] * exp(p[2] * x) + p[3]
true_p = [1.7, -1.3, 0.15]

function jacobian(x, p)
    a, b, _ = p
    ex = exp.(b .* x)
    J = Matrix{Float64}(undef, length(x), 3)
    J[:, 1] .= ex
    J[:, 2] .= a .* x .* ex
    J[:, 3] .= 1.0
    return J
end

sigma_y = fill(0.05, length(x))
y = model(x, true_p) .+ sigma_y .* sin.(4.0 .* x)

result = fit_model(
    model,
    x,
    y;
    p0=[1.2, -0.7, 0.0],
    sigma_y=sigma_y,
    jacobian=jacobian,
    backend=:lsqfit,
)

println("backend: ", result.backend)
println("params: ", result.params)
println("correlation matrix:\n", result.param_correlation)

# Demonstrate explicit x-grid, labels/title, and style options.
xgrid = collect(range(0.0, 2.2; length=600))
plot_fit(
    result;
    xgrid=xgrid,
    filename=example_output("analytic_jacobian_plot"),
    format=:svg,
    theme=:publication,
    title="Exponential Decay Fit",
    xlabel="time t (s)",
    ylabel="signal S(t)",
    parameter_names=["A", "lambda", "C"],
    fit_color=:tomato4,
    fit_linewidth=4,
    band_color=:tomato2,
    band_alpha=0.25,
    data_color=:black,
    data_marker=:utriangle,
    data_markersize=12,
    xerr_color=:gray35,
    yerr_color=:gray35,
    error_whiskerwidth=8,
    stats_sigdigits=6,
    stats_fontsize=18,
)
println("Saved plot to ", example_output("analytic_jacobian_plot.svg"))
