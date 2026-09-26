# Results And Diagnostics

Fit construction is covered by [Fitting](api_fitting.md).

## Results

### Fit Result Fields

`FitResult` and `LikelihoodFitResult` use the same parameter and status field
names.

| Field | Meaning |
|---|---|
| `problem` | Validated problem used for the selected candidate. |
| `options` | Normalized solver options. |
| `backend` | Selected backend, for example `:lsqfit`, `:optimization`, `:native_minuit`, or `:fixed`. |
| `converged` | Whether the selected solver reported convergence. |
| `iterations` | Iteration count, or `missing` when unavailable. |
| `message` | Native solver termination message. |
| `params` | Best-fit parameter vector. |
| `param_stderr` | Local one-standard-deviation estimates from the covariance diagonal. |
| `param_covariance` | Local parameter covariance matrix. |
| `param_correlation` | Correlation matrix derived from that covariance. |
| `stats` | [`ScientificFitting.FitStatistics`](@ref). |
| `diagnostics` | Numerical checks computed during result construction. |
| `solver_result` | [`FitSolverResult`](@ref) with native details and free-parameter mapping for scalar adapters; otherwise `nothing`. |

Only `FitResult` has `model_y`, `residuals`, `weighted_residuals`, and
`jacobian`. With a non-diagonal covariance, `weighted_residuals` are whitened
coordinates, not pointwise pulls
([Residuals And Pulls](statistics.md#Residuals-And-Pulls)).

### Predictions Without A Plotting Backend

`predict(result, x)` evaluates a Gaussian fit on new coordinates; with
`uncertainty=true` it also returns the local standard uncertainty of the
fitted mean.

```@docs
ScientificFitting.predict
```

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

An arbitrary custom loss has no likelihood interpretation, and for indexed and
multi-dataset wrappers `minus2loglik_min` equals the chi-square objective only
when no normalized Gaussian parameter terms are present; see
[Likelihoods and Model Comparison](statistics.md#Observation-Likelihoods) before comparing AIC
or BIC.

```@docs
ScientificFitting.FitResult
ScientificFitting.LikelihoodFitResult
ScientificFitting.FitStatistics
ScientificFitting.FitDiagnostics
```

## [Parameter Covariance](@id parameter-covariance-reference)

`param_covariance` is a local quadratic approximation; its failure modes are
derived in [Local Parameter Covariance](statistics.md#Local-Parameter-Covariance)
and the `scale_covariance` policy in
[Covariance Scaling](statistics.md#Covariance-Scaling).

Use [`profile_interval`](@ref) for asymmetric one-parameter intervals and
[`profile_matrix`](@ref) when several parameters may be correlated or
non-parabolic.

## Profiles And Contours

Profiles fix the displayed parameter or parameter pair and re-optimize every
remaining free parameter
([Profiles And Contours](statistics.md#Profiles-And-Contours)). The same
functions accept `FitResult` and `LikelihoodFitResult`.

Common asymptotic thresholds on the ``-2\log L`` or chi-square scale:

| Coverage | One profiled parameter | Two profiled parameters |
|---:|---:|---:|
| 68.27% | `threshold = 1.00` | `levels = [2.30]` |
| 95.45% | `threshold = 4.00` | `levels = [6.18]` |

Defaults are `1.00` for profiles and `[2.30, 6.18]` for contours.

Failed refits become `Inf` by default and are surfaced by diagnostics; a
finite objective from a non-converged nuisance fit also counts as a failure.
Refits inherit the original solver limits and tolerances.

| Scan control | Meaning |
|---|---|
| `values`, `xvalues`, `yvalues` | Explicit finite scan coordinates; replace the automatic range. |
| `npoints` | Resolution of an automatically generated axis. |
| `nsigma` | Half-width of the automatic range in local standard errors. |
| `threshold`, `levels` | Positive delta-cost thresholds for intervals or regions. |
| `adaptive` | Refine only threshold-crossing intervals or cells. |
| `max_refinements`, `max_points` | Bound adaptive work and total scan size. |
| `on_failure` | `:inf` records a failed refit; `:throw` stops immediately. |

`profile_interval` linearly interpolates threshold crossings; a side that is
not bracketed is returned as `NaN`. The search stops at failed grid points
rather than interpolating across gaps. Pass a completed `ProfileResult` to
extract its interval without additional refits. `profile_matrix` accepts
`parameters` and `parameter_names`; `profile_tolerance` and
`contour_tolerance` compare scans with local quadratic geometry.

```@docs
ScientificFitting.profile
ScientificFitting.profile_interval
ScientificFitting.contour
ScientificFitting.profile_matrix
ScientificFitting.ProfileResult
ScientificFitting.ProfileInterval
ScientificFitting.ContourResult
ScientificFitting.ProfileMatrixResult
```

## Diagnostics And Reports

| Need | Function | Return value |
|---|---|---|
| Programmatic findings | [`diagnose`](@ref) | [`ScientificFitting.DiagnosticReport`](@ref) |
| Full diagnostic text | [`diagnose_text`](@ref) | `String` |
| Short action list | [`diagnostic_dashboard`](@ref) | [`ScientificFitting.DiagnosticDashboard`](@ref) |
| Dashboard text | [`diagnostic_dashboard_text`](@ref) | `String` |
| Structured fit report | [`fit_report`](@ref) | [`ScientificFitting.FitReport`](@ref) |
| Console or notebook report | [`report_text`](@ref) | `String` |

Text output renders `:stop` as `critical - fix before use`; `:ok` is not proof
that the physical model is true.

`max_actions=0` suppresses the deduplicated action list without removing any
findings.
`fit_report(...; errors=:profile)` replaces local symmetric display errors with
profile intervals, performs the additional fits, and leaves unbracketed sides
as `NaN`. Its `profile_threshold`, `profile_npoints`, and `profile_nsigma`
keywords control those scans. `report_text(...; sigdigits=6)` controls
numerical formatting only.

```@docs
ScientificFitting.DiagnosticFinding
ScientificFitting.DiagnosticReport
ScientificFitting.DiagnosticDashboard
ScientificFitting.diagnose
ScientificFitting.diagnose_text
ScientificFitting.diagnostic_dashboard
ScientificFitting.diagnostic_dashboard_text
ScientificFitting.ParameterEstimate
ScientificFitting.FitReport
ScientificFitting.fit_report
ScientificFitting.report_text
```
