# Results And Diagnostics

This page defines fitted results, local and profile uncertainty, and diagnostic
output. Fit construction is covered by [Fitting](api_fitting.md).

## Results

### Fit Result Fields

`FitResult` and `LikelihoodFitResult` use the same parameter and status field
names. `FitResult` additionally contains x-y-specific model and residual data.

| Field | Meaning |
|---|---|
| `problem` | Validated problem used for the selected candidate. |
| `options` | Normalized solver options. |
| `backend` | `:lsqfit`, `:optimization`, or `:fixed`. |
| `converged` | Whether the selected solver reported convergence. |
| `iterations` | Iteration count, or `missing` when unavailable. |
| `message` | Native solver termination message. |
| `params` | Best-fit parameter vector. |
| `param_stderr` | Local one-standard-deviation estimates from the covariance diagonal. |
| `param_covariance` | Local parameter covariance matrix. |
| `param_correlation` | Correlation matrix derived from that covariance. |
| `stats` | [`FitStatistics`](@ref). |
| `diagnostics` | Numerical checks computed during result construction. |

Only `FitResult` has `model_y`, `residuals`, `weighted_residuals`, and
`jacobian`. `weighted_residuals` are pulls only when the uncertainty model
supports that pointwise interpretation; with full whitening they are whitened
coordinates.

### Fit Statistics

| Field | Meaning |
|---|---|
| `cost` | Symbol identifying the minimized cost. |
| `cost_min` | Minimized cost including parameter terms. |
| `minus2loglik_min` | Value used for likelihood-derived summaries. It is a normalized ``-2\log L`` only when the objective follows that convention. |
| `chi2` | Chi-square or deviance goodness-of-fit statistic, otherwise `NaN`. |
| `chi2_ndf` | `chi2 / ndf` when defined. |
| `ndf` | Independent observations and Gaussian constraint dimensions minus free parameters. |
| `pvalue` | Upper-tail chi-square probability when a reference distribution exists. |
| `aic`, `bic` | Information criteria; meaningful only for compatible likelihood normalizations. |

An arbitrary custom loss has no likelihood interpretation. For indexed and
multi-dataset wrappers, `minus2loglik_min` equals the chi-square objective only
when no normalized Gaussian parameter terms are present. See
[Likelihoods and Model Comparison](likelihood_models.md) before comparing AIC
or BIC across different data or uncertainty models.

```@docs
JuFitter.FitResult
JuFitter.LikelihoodFitResult
JuFitter.FitStatistics
JuFitter.FitDiagnostics
```

## [Parameter Covariance](@id parameter-covariance-reference)

`param_covariance` is a local quadratic approximation. It can be misleading for
nonlinear models, weak data, active bounds, asymmetric likelihoods, or multiple
minima. `scale_covariance` changes only residual-scale treatment; it does not
make a non-quadratic likelihood Gaussian.

Use [`profile_interval`](@ref) for asymmetric one-parameter intervals and
[`profile_matrix`](@ref) when several parameters may be correlated or
non-parabolic.

## Profiles And Contours

Profiles fix the displayed parameter or parameter pair and re-optimize every
remaining free parameter. A scan point is therefore a fit, not merely an
evaluation of the original model. The same functions accept `FitResult` and
`LikelihoodFitResult`.

For costs on the ``-2\log L`` or chi-square scale, common asymptotic thresholds
are:

| Coverage | One profiled parameter | Two profiled parameters |
|---:|---:|---:|
| 68.27% | `threshold = 1.00` | `levels = [2.30]` |
| 95.45% | `threshold = 4.00` | `levels = [6.18]` |

Defaults are `1.00` for profiles and `[2.30, 6.18]` for contours. These are
Wilks-theorem approximations, not universal finite-sample guarantees.

`adaptive=true` refines threshold crossings or contour cells instead of making
the entire rectangular grid dense. Failed refits become `Inf` by default and
are surfaced by diagnostics; use `on_failure=:throw` to stop at the first failed
point.

| Scan control | Meaning |
|---|---|
| `values`, `xvalues`, `yvalues` | Explicit finite scan coordinates; replace the automatic range. |
| `npoints` | Resolution of an automatically generated axis. |
| `nsigma` | Half-width of the automatic range in local standard errors. |
| `threshold`, `levels` | Positive delta-cost thresholds for intervals or regions. |
| `adaptive` | Refine only threshold-crossing intervals or cells. |
| `max_refinements`, `max_points` | Bound adaptive work and total scan size. |
| `on_failure` | `:inf` records a failed refit; `:throw` stops immediately. |

`profile_interval` linearly interpolates threshold crossings.
A side that is not bracketed is returned as `NaN`, not silently extrapolated.
`profile_matrix` accepts `parameters` and `parameter_names`; `profile_tolerance`
and `contour_tolerance` compare scans with local quadratic geometry.

```@docs
JuFitter.profile
JuFitter.profile_interval
JuFitter.contour
JuFitter.profile_matrix
JuFitter.profile_matrix_triage
JuFitter.ProfileResult
JuFitter.ProfileInterval
JuFitter.ContourResult
JuFitter.ProfileMatrixResult
JuFitter.ProfileMatrixPanelTriage
```

## Diagnostics And Reports

| Need | Function | Return value |
|---|---|---|
| Programmatic findings | [`diagnose`](@ref) | [`DiagnosticReport`](@ref) |
| Full diagnostic text | [`diagnose_text`](@ref) | `String` |
| Short action list | [`diagnostic_dashboard`](@ref) | [`DiagnosticDashboard`](@ref) |
| Dashboard text | [`diagnostic_dashboard_text`](@ref) | `String` |
| Structured fit report | [`fit_report`](@ref) | [`FitReport`](@ref) |
| Console or notebook report | [`report_text`](@ref) | `String` |

Dashboard status is `:ok`, `:review`, or `:stop`. Text output renders `:stop` as
`critical - fix before use`. The dashboard summarizes implemented checks;
`:ok` is not proof that the physical model is true.

`diagnostic_dashboard(...; max_actions=5)` limits the deduplicated action list.
`fit_report(...; errors=:profile)` replaces local symmetric display errors with
profile intervals and therefore performs additional fits. Its
`profile_threshold`, `profile_npoints`, and `profile_nsigma` keywords control
those scans. `report_text(...; sigdigits=6)` controls numerical formatting only.

```@docs
JuFitter.DiagnosticFinding
JuFitter.DiagnosticReport
JuFitter.DiagnosticDashboard
JuFitter.diagnose
JuFitter.diagnose_text
JuFitter.diagnostic_dashboard
JuFitter.diagnostic_dashboard_text
JuFitter.ParameterEstimate
JuFitter.FitReport
JuFitter.fit_report
JuFitter.report_text
```
