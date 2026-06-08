# API Reference

This page is the compact public reference for JuFitter. It is organized by how
the package is normally used: define a fit, describe uncertainties and
constraints, inspect the result, diagnose failure modes, and make plots. For a
worked path through the same ideas, start with [Quickstart](quickstart.md) or
the [Gallery](gallery.md).

## Fitting Workflows

Use these entry points for ordinary model fits and likelihood fits.

```@docs
JuFitter.fit
JuFitter.fit_model
JuFitter.fit_custom
JuFitter.fit_poisson_model
JuFitter.fit_histogram_model
JuFitter.fit_histogram_density
JuFitter.fit_unbinned_model
JuFitter.fit_extended_unbinned_model
JuFitter.fit_indexed_model
JuFitter.fit_multi_model
```

## Fit Problems And Options

These types describe the normalized problem JuFitter solves. Most users create
them indirectly through the fitting functions, but they are useful when a
workflow needs explicit construction, testing, or serialization.

```@docs
JuFitter.FitProblem
JuFitter.LikelihoodFitProblem
JuFitter.FitOptions
```

## Uncertainties, Bounds, And Constraints

Use these structures when pointwise `sigma_y`/`sigma_x` is not enough. Dense
covariance is exact for small and medium datasets, but large correlated
measurements should eventually use structured covariance or custom whitening
operators; see the roadmap and release audit for that limitation.

```@docs
JuFitter.ErrorComponent
JuFitter.ParameterPrior
JuFitter.FixedParameter
JuFitter.ParameterConstraint
JuFitter.ConstraintSpec
```

## Results And Reports

These objects hold fit outputs and make the numerical assumptions explicit.
Parameter covariance is local: for nonlinear models, active bounds, weak data,
or visibly asymmetric likelihoods, inspect profiles or contours before treating
symmetric errors as final.

```@docs
JuFitter.FitResult
JuFitter.LikelihoodFitResult
JuFitter.FitStatistics
JuFitter.FitDiagnostics
JuFitter.FitReport
JuFitter.ParameterEstimate
JuFitter.fit_report
JuFitter.report_text
```

## Diagnostics

Diagnostics are structured so they can be printed in a terminal, shown in a
notebook, or reused by plots. The dashboard is intentionally action-oriented:
it should tell a lab user what to inspect next.

```@docs
JuFitter.DiagnosticFinding
JuFitter.DiagnosticReport
JuFitter.DiagnosticDashboard
JuFitter.diagnose
JuFitter.diagnose_text
JuFitter.diagnostic_dashboard
JuFitter.diagnostic_dashboard_text
```

## Profiles And Contours

Use profile and contour scans when the local covariance approximation is not
enough. The one-dimensional profile fixes one parameter and refits the others;
the two-dimensional contour repeats that idea over a parameter pair.

```@docs
JuFitter.ProfileResult
JuFitter.ProfileInterval
JuFitter.ContourResult
JuFitter.ProfileMatrixResult
JuFitter.ProfileMatrixPanelTriage
JuFitter.profile
JuFitter.profile_interval
JuFitter.contour
JuFitter.profile_matrix
JuFitter.profile_matrix_triage
```

## Plotting

Plotting is loaded through the optional CairoMakie extension. Fitting and text
reports work without Makie; plot calls require `using CairoMakie`.

```@docs
JuFitter.fitplot
JuFitter.plot_fit
JuFitter.fit_axis
JuFitter.add_curve!
JuFitter.add_points!
JuFitter.add_vline!
JuFitter.add_hline!
JuFitter.add_vband!
JuFitter.add_hband!
JuFitter.plot_theme
JuFitter.plot_palette
JuFitter.plot_info_panel!
JuFitter.plot_residuals
JuFitter.plot_diagnostics
JuFitter.plot_profile
JuFitter.plot_contour
JuFitter.plot_profile_matrix
```
