# API Overview

JuFitter's public API is organized around one workflow: define a scientific
problem, fit it, inspect diagnostics, then decide how much output you want.
Most users should start with `fitplot` or `fit_model`; lower-level objects are
available when a notebook grows into a reusable analysis script.

## Fitting Entry Points

- `fitplot(x, y; ...)` fits common array data and returns a Makie figure.
- `fit_model(model, x, y; p0, ...)` returns a `FitResult` without forcing a
  plot-first workflow.
- `FitProblem(...)` stores data, model, uncertainties, constraints, and backend
  choices explicitly for advanced workflows.
- Likelihood helpers cover Poisson counts, histograms, unbinned samples,
  indexed likelihoods, custom objectives, and multi-dataset costs.

## Result Objects

`FitResult` is the object to keep. It stores the best-fit parameters,
covariance estimate, residuals, degrees of freedom, goodness-of-fit statistics,
optimizer diagnostics, and enough metadata for later reports and plots.

Use `report_text(result)` for console output, `fit_report(result)` for a
structured report object, and `diagnostic_dashboard_text(result)` when a fit
needs a quick "what should I inspect next?" summary.

## Plotting Entry Points

`plot_fit(result; ...)` turns an existing result into a Makie figure. It does
not refit. This is deliberate: a user can add markers, thresholds, reference
bands, or additional curves without changing the numerical result.

Useful plot controls include:

- `report=:none`, `:console`, `:plot`, or `:both`,
- `show_legend`, `show_stats`, and `stats_position`,
- `style=:workbench`, `:showcase`, or `:publication`,
- `appearance=:light` or `:dark`,
- `sigma_band` and band labels for uncertainty semantics.

Post-fit helpers such as `fit_axis`, `add_curve!`, `add_points!`, `add_vline!`,
`add_hline!`, `add_vband!`, and `add_hband!` are thin Makie wrappers for common
scientific annotations.

## Diagnostics and Uncertainty Checks

Profiles and contours are not decorative plots. They test whether a local
parabolic covariance approximation is credible. Use `profile_interval`,
`profile_curve`, `contour_grid`, and `plot_profile_matrix` when parameters are
bounded, strongly correlated, weakly constrained, or visibly nonlinear.

For the conceptual pipeline, see [How JuFitter Works](@ref). For the statistical
assumptions behind the methods, see [Statistical Foundations](@ref).
