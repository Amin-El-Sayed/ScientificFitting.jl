# Quickstart

This is the smallest useful JuFitter workflow: data, uncertainties, model,
fit, plot, report.

```julia
using JuFitter

x = collect(range(0.0, 10.0; length=40))
model(x, p) = @. p[1] * x + p[2]

sigma_y = fill(0.25, length(x))
y = model(x, [1.7, 0.8]) .+ sigma_y .* sin.(1.2 .* x)

plot_result = fitplot(
    model,
    x,
    y;
    p0=[1.0, 0.0],
    sigma_y=sigma_y,
    parameter_names=["slope", "offset"],
    xlabel="time / s",
    ylabel="signal / V",
    filename="quickstart_fit.pdf",
)

result = plot_result.result
println(report_text(result; parameter_names=["slope", "offset"]))
println(diagnose_text(result))
```

## What Happens

`fitplot` calls `fit_model`, creates a `FitResult`, and renders the default
fit plot with uncertainty band and optional result panel. If you only need the
numbers, call `fit_model` directly:

```julia
result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y)
```

## Choosing The Cost

For ordinary Gaussian y uncertainties, `cost=:auto` uses chi-square:

```julia
result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y)
```

For parameter-dependent covariance, for example x uncertainties, `cost=:auto`
uses the full Gaussian negative log-likelihood:

```julia
sigma_x = fill(0.03, length(x))
result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_x=sigma_x, sigma_y=sigma_y)
```

For count data, use the Poisson helper:

```julia
rate(x, p) = @. exp(p[1] + p[2] * x)
counts = round.(rate(x, [1.0, 0.12]))

result = fit_poisson_model(rate, x, counts; p0=[0.8, 0.1])
```

## Reading The Result

Important fields are:

- `result.params`: fitted parameters.
- `result.param_stderr`: local one-sigma standard errors.
- `result.param_covariance`: local parameter covariance matrix.
- `result.stats`: cost, chi-square or deviance, ndf, p-value, AIC, and BIC.
- `result.diagnostics`: warnings for non-convergence, active bounds,
  ill-conditioned covariance, and unavailable goodness-of-fit statistics.
- `diagnose(result)`: actionable troubleshooting findings with evidence and
  recommended next steps.

Local errors are only local approximations. For nonlinear fits, active bounds,
or asymmetric likelihoods, use `profile`, `profile_interval`, and `contour`.
