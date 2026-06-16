# Reference Map

This page is the short map of JuFitter's public surface. Use it when you know
what kind of task you have and need the right entry point. For signatures and
field names, use the [API Reference](api.md).

## Start With The Workflow

Most analyses follow one of these paths:

| task | entry point | result |
| --- | --- | --- |
| quick x-y fit with plot | `fitplot(...)` | fit result plus Makie figure |
| fit first, plot later | `fit_model(...)` | `FitResult` |
| explicit reusable problem | `FitProblem(...)` then `fit(...)` | `FitResult` |
| counts or distributions | likelihood helpers | `LikelihoodFitResult` |
| custom objective | `fit_custom(...)` | `LikelihoodFitResult` |

The object to keep is the result object. It stores the parameters, local
covariance, residual information where available, goodness-of-fit statistics,
diagnostics, and enough metadata to regenerate reports or plots without
rewriting the fit.

## Describe The Uncertainty Model

Choose the simplest uncertainty representation that is scientifically honest:

- `sigma_y` for independent vertical standard uncertainties,
- `sigma_x` for relevant x uncertainty propagated through the model slope,
- `cov_y` or `cov_x` when measured points are correlated,
- `ErrorComponent` when named uncertainty contributions should remain visible,
- `ParameterPrior`, `ParameterConstraint`, and `FixedParameter` for external
  information about parameters.

Dense covariance is exact but expensive. It is appropriate for small and
medium correlated datasets, not for every long time series, image, spectrum, or
detector array. Large structured correlations are tracked as future whitening
operator work in the release audit.

## Read The Result Before Plotting

For terminal or notebook output, use:

```julia
println(report_text(result))
println(diagnostic_dashboard_text(result))
```

The main fields are:

- `result.params`: best-fit parameter values,
- `result.param_covariance`: local parameter covariance matrix,
- `result.param_stderr`: local symmetric parameter errors,
- `result.stats`: chi-square, likelihood, p-value, AIC, and BIC summary,
- `result.diagnostics`: numerical and statistical warnings.

Local covariance is a quadratic approximation around the minimum. It is fast
and useful, but it can be misleading for nonlinear models, weak data, active
bounds, strong correlations, or asymmetric likelihoods.

## Diagnose Before Trusting

Use `diagnose(result)` or `diagnostic_dashboard(result)` when a fit looks wrong
or when it matters scientifically. The diagnostics check convergence, degrees
of freedom, p-values, covariance/Hessian conditioning, active bounds, strong
correlations, large pulls, and residual structure.

If local covariance looks suspicious, compute profile or contour diagnostics:

```julia
prof = profile(result, 1; adaptive=true)
interval = profile_interval(result, 1; adaptive=true)

cont = contour(result, 1, 2; adaptive=true)
```

Then inspect them visually:

```julia
plot_profile(prof; local_sigma=result.param_stderr[1])
plot_contour(
    cont;
    local_covariance=result.param_covariance,
    local_center=result.params[[1, 2]],
)
plot_profile_matrix(result; parameters=[1, 2, 3])
```

Profiles and contours answer a different question than the main fit plot: they
show whether the uncertainty geometry near the minimum is close enough to a
local Gaussian approximation.

## Plot Without Re-Fitting

`plot_fit(result; ...)` creates a Makie figure from an existing result. It does
not change the numerical fit. This is the intended pattern for article
publication-quality figures: fit once, then add annotations on top.

```julia
fig = plot_fit(result; report=:plot, show_legend=true)
ax = fit_axis(fig)
add_vline!(ax, threshold; label="threshold")
add_curve!(ax, reference_model; xspan=(0, 10), label="reference")
```

Use `plot_theme` and `plot_palette` when building custom Makie layouts that
should follow JuFitter's `:lab`, `:modern`, or `:article` styles.

## When To Leave This Page

- First fit: [Quickstart](quickstart.md).
- Full examples: [Gallery](gallery.md).
- Statistical meaning: [Statistical Foundations](statistical_foundations.md).
- Exact signatures: [API Reference](api.md).
- Implementation structure: [Backend Design](backend_design.md).
