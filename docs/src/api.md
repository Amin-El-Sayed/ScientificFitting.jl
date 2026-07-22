# API Reference

This is the technical reference for JuFitter's public API. If you are fitting a
first dataset, begin with [Quickstart](quickstart.md). If you need the reasoning
behind a statistical choice, use [Statistical Foundations](statistical_foundations.md).

## Choose An Entry Point

The public fitting functions differ by observation model, not by plotting
style or optimizer.

| Data and sampling model | Entry point | Model contract | Result |
|---|---|---|---|
| Numeric ``x`` and ``y`` with Gaussian uncertainties | [`fit_model`](@ref) | `model(x, p) -> ŷ` | [`FitResult`](@ref) |
| Independent counts | [`fit_poisson_model`](@ref) | `model(x, p) -> expected_counts` | [`LikelihoodFitResult`](@ref) |
| Histogram with precomputed expected bin counts | [`fit_histogram_model`](@ref) | `expected_counts(edges, p) -> μ` | [`LikelihoodFitResult`](@ref) |
| Histogram from a normalized density | [`fit_histogram_density`](@ref) | `pdf(x, p) -> density` | [`LikelihoodFitResult`](@ref) |
| Independent unbinned observations | [`fit_unbinned_model`](@ref) | `pdf(x, p) -> density` | [`LikelihoodFitResult`](@ref) |
| Unbinned events with a parameter-dependent rate | [`fit_extended_unbinned_model`](@ref) | `rate(x, p) -> event_rate` | [`LikelihoodFitResult`](@ref) |
| Observations addressed by non-numeric indices | [`fit_indexed_model`](@ref) | `model(indices, p) -> ŷ` | [`LikelihoodFitResult`](@ref) |
| Several datasets sharing parameters | [`fit_multi_model`](@ref) | one `model_i(x_i, p_i)` per dataset | [`LikelihoodFitResult`](@ref) |
| A custom scalar objective | [`fit_custom`](@ref) | `objective(p) -> scalar` | [`LikelihoodFitResult`](@ref) |

For reusable low-level workflows, construct [`FitProblem`](@ref) or
[`LikelihoodFitProblem`](@ref) and call [`fit`](@ref). `JuFitter.fit` extends
the standard `StatsAPI.fit` generic, so it coexists with packages such as
StatsBase and Distributions.

## Common Conventions

### Parameters And Model Functions

`p0` fixes the parameter order. Parameter indices in bounds, fixed values,
constraints, profiles, and contours are Julia's one-based indices into that
vector.

The ordinary allocating model contract is

```julia
model(x, p) -> vector with length(y)
```

For the allocation-sensitive path, use

```julia
model!(out, x, p)
jacobian!(J, x, p)  # optional
```

and pass `inplace=true`. JuFitter validates at construction time that these
functions fill every output. On automatic-differentiation paths, `p`, `out`,
and `J` may contain non-`Float64` scalar types; mutating functions must not
hard-code `Float64` buffers internally.

Input observations and starting values are copied to `Float64` storage. A
model must return finite values at every parameter point used by the solver.

### Parameter Control

The controls below are accepted by `fit_model` and all likelihood wrappers.

| Keyword | Meaning |
|---|---|
| `p0` | Required complete starting vector. |
| `bounds=(lower, upper)` | Componentwise closed bounds; use `±Inf` for an open side. |
| `fixed_parameters` | Remove parameters from the optimizer with `FixedParameter`, `i => value`, or equivalent named tuples. |
| `parameter_priors` | Independent normalized Gaussian or split-normal terms. |
| `parameter_constraints` | Correlated Gaussian terms on selected parameters. |
| `constraints` | General nonlinear constraints; `ineq(p) <= 0` and `eq(p) == 0`. |

Likelihood wrappers additionally accept `parameter_names`, which are stored in
the low-level problem for reports and diagnostics. For a `FitResult`, pass
reader-facing names to `fit_report`, `report_text`, or the plotting function.

`FixedParameter(..., sigma)` records external uncertainty but does not make the
parameter free. To propagate an uncertain external measurement through the
fit, use a `ParameterPrior` or `ParameterConstraint` instead.

### Solver Control

| Keyword | Default | Contract |
|---|---:|---|
| `maxiters` | `500` for `fit_model`, `1000` for likelihood wrappers | Positive iteration limit for each candidate. |
| `tol` | `1e-10` | Positive absolute and relative solver tolerance. |
| `initial_guesses` | `nothing` | Additional complete starting vectors. |
| `multistart` | `1` | Number of deterministic candidates generated from finite bounds or scaled versions of `p0`. |

`fit_model` additionally accepts:

| Keyword | Default | Contract |
|---|---:|---|
| `backend` | `:auto` | `:auto`, `:lsqfit`, or `:optimization`. |
| `cost` | `:auto` | `:chi2` or full `:gaussian_likelihood` on the ``-2\log L`` scale; `:auto` uses the latter for parameter-dependent covariance. |
| `scale_covariance` | `:auto` | `:auto`, `:always`, or `:never`; see [Parameter Covariance](@ref parameter-covariance-reference). |
| `jacobian` | `nothing` | Analytic model Jacobian, allocating or in-place according to `inplace`. |
| `x_derivative` | `nothing` | Vector ``\partial f/\partial x`` for efficient x-uncertainty propagation. |

`backend=:auto` uses the fast LsqFit path only when static chi-square least
squares represents the complete problem. An explicit incompatible
`backend=:lsqfit` request raises an error rather than silently dropping bounds,
constraints, priors, or parameter-dependent covariance.

## Gaussian Fits

### Observation Uncertainty

Choose one representation for each physical uncertainty source. JuFitter
rejects contradictory combinations instead of guessing how they should be
combined.

| Keyword | Accepted value | Statistical role |
|---|---|---|
| `sigma_y` | positive vector | Independent y standard deviations. |
| `cov_y` | dense or sparse SPD matrix | Complete y covariance. Mutually exclusive with `sigma_y`. |
| `sigma_x` | positive vector | Independent x standard deviations propagated through ``\partial f/\partial x``. |
| `cov_x` | dense or sparse SPD matrix | Complete x covariance. Mutually exclusive with `sigma_x`. |
| `error_components` | named [`ErrorComponent`](@ref)s | Additive absolute, relative, model-relative, or covariance contributions. |
| `whitening` | [`WhiteningOperator`](@ref) | Complete static observation covariance represented by ``W^\mathsf{T}W=C^{-1}``. |

`whitening` is intentionally exclusive with all other observation-uncertainty
keywords. It describes the complete covariance; combining it with another
source without an explicit derivation would double-count uncertainty.

With no supplied observation uncertainty, `fit_model` performs unweighted
least squares and `scale_covariance=:auto` estimates the residual scale from
``\chi^2/\mathrm{ndf}``. With physical uncertainties, `:auto` leaves their
scale unchanged.

### Error Components

An error component has a stable name and can be activated or deactivated
without rewriting the fit:

```julia
ErrorComponent(:readout, :y, :absolute, sigma_readout)
ErrorComponent(:gain, :y, :relative, 0.015)
ErrorComponent(:calibration, :y, :model_relative, 0.008)
ErrorComponent(:shared, :y, :covariance, covariance_matrix)
```

`target` is `:x` or `:y`. `mode` is `:absolute`, `:relative`,
`:model_relative`, or `:covariance`; x components do not support
`:model_relative`.

### Fit Completion And Failure

Scientific input errors such as non-finite observations, non-positive standard
deviations, invalid bounds, or non-positive-definite covariance matrices raise
`ArgumentError` before optimization.

For multistart fits, JuFitter returns the converged candidate with the lowest
finite cost. If no candidate converges but one returns a finite result, that
result is returned with `converged == false`; callers must inspect the status or
use [`diagnostic_dashboard`](@ref). If every candidate fails, the underlying
error is raised.

```@docs
JuFitter.fit(::JuFitter.FitProblem)
JuFitter.fit_model
JuFitter.FitProblem
JuFitter.FitOptions
```

## Likelihood And Count Fits

Poisson, histogram, unbinned, and extended-unbinned entry points minimize costs
on the ``-2\log L`` scale. Poisson and histogram fits also compute Poisson
deviance, so `chi2`, `chi2_ndf`, and `pvalue` are available as goodness-of-fit
summaries. Ordinary and extended unbinned fits do not invent a chi-square
statistic; those fields are `NaN`.

`fit_indexed_model` and `fit_multi_model` minimize chi-square but omit additive
Gaussian normalization constants. Their AIC/BIC values may compare models fit
to the same observations with the same uncertainty model; they must not be used
to compare different uncertainty scales or datasets.

| Entry point | Additional contract |
|---|---|
| `fit_poisson_model` | Every expected count must be finite and strictly positive; observed counts must be non-negative integers. |
| `fit_histogram_model` | `length(edges) == length(counts) + 1`; edges increase strictly; the model returns one positive expectation per bin. |
| `fit_histogram_density` | Integrates `pdf(x, p)` over every bin with Gauss-Kronrod quadrature; `total_count > 0`, `rtol > 0`. |
| `fit_unbinned_model` | The supplied density must already be normalized and positive at every observation. |
| `fit_extended_unbinned_model` | `rate` is an intensity, not a density; its integral over `domain` is the expected event count. |
| `fit_indexed_model` | Supports `sigma_y` or `cov_y`; indices may be any container accepted by the model. |
| `fit_multi_model` | Supports per-dataset `sigma_y`; `parameter_map[i]` selects global parameters passed to model `i`. |

For `fit_custom`, `objective` should be a normalized ``-2\log L`` cost if local
covariance, AIC, and BIC are to retain their standard interpretation. With an
arbitrarily scaled loss, optimization still works, but these inferential fields
are only arithmetic summaries. `nobs` must count statistically independent
observations.

```@docs
JuFitter.fit(::JuFitter.LikelihoodFitProblem)
JuFitter.fit_custom
JuFitter.fit_poisson_model
JuFitter.fit_histogram_model
JuFitter.fit_histogram_density
JuFitter.fit_unbinned_model
JuFitter.fit_extended_unbinned_model
JuFitter.fit_indexed_model
JuFitter.fit_multi_model
JuFitter.LikelihoodFitProblem
```

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
| `minus2loglik_min` | Gaussian or likelihood cost on the ``-2\log L`` scale; equal to chi-square for indexed/multi wrappers. |
| `chi2` | Chi-square or deviance goodness-of-fit statistic, otherwise `NaN`. |
| `chi2_ndf` | `chi2 / ndf` when defined. |
| `ndf` | Independent observations and Gaussian constraint dimensions minus free parameters. |
| `pvalue` | Upper-tail chi-square probability when a reference distribution exists. |
| `aic`, `bic` | Information criteria; meaningful only for compatible likelihood normalizations. |

See [Statistical Foundations](statistical_foundations.md) before comparing AIC
or BIC across different data, uncertainty models, or arbitrarily normalized
objectives.

```@docs
JuFitter.FitResult
JuFitter.LikelihoodFitResult
JuFitter.FitStatistics
JuFitter.FitDiagnostics
```

## [Parameter Covariance](@id parameter-covariance-reference)

`param_covariance` is a local quadratic approximation. It can be misleading
for nonlinear models, weak data, active bounds, asymmetric likelihoods, or
multiple minima. `scale_covariance` changes only its residual-scale treatment;
it does not make a non-quadratic likelihood Gaussian.

Use [`profile_interval`](@ref) for asymmetric one-parameter intervals and
[`profile_matrix`](@ref) when several parameters may be correlated or
non-parabolic.

## Profiles And Contours

Profiles fix the displayed parameter or parameter pair and re-optimize every
remaining free parameter. Therefore a scan point is a fit, not merely an
evaluation of the original model.

For costs on the ``-2\log L`` or chi-square scale, common asymptotic thresholds
are:

| Coverage | One profiled parameter | Two profiled parameters |
|---:|---:|---:|
| 68.27% | `threshold = 1.00` | `levels = [2.30]` |
| 95.45% | `threshold = 4.00` | `levels = [6.18]` |

The defaults use `1.00` for profiles and `[2.30, 6.18]` for contours. These are
Wilks-theorem approximations, not universal finite-sample guarantees.

`adaptive=true` refines threshold-crossing intervals or contour cells instead
of making the complete rectangular grid dense. Failed refits become `Inf` by
default and are surfaced by profile diagnostics; use `on_failure=:throw` when
the first failed point should stop the scan.

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
| Console/notebook report | [`report_text`](@ref) | `String` |

Dashboard status is `:ok`, `:review`, or `:critical`. It summarizes implemented
checks; `:ok` is not proof that the physical model is true.

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

## Constraints And Uncertainty Objects

```@docs
JuFitter.ConstraintSpec
JuFitter.ParameterPrior
JuFitter.FixedParameter
JuFitter.ParameterConstraint
JuFitter.ErrorComponent
JuFitter.WhiteningOperator
```

## Plotting

Fitting, profiles, diagnostics, and text reports do not require Makie. Plotting
is activated only after loading CairoMakie:

```julia
using JuFitter
using CairoMakie
```

[`plot_fit`](@ref) returns a `Figure`. [`fitplot`](@ref) fits and returns the
named tuple `(result, figure)`. Use `fit_axis(figure)` as the stable extension
point for native Makie calls or JuFitter's `add_*!` helpers.

`report=:plot`, `:console`, `:both`, or `:none` controls the information panel
and terminal output. `fit_range=:axis` extrapolates the fitted curve over the
padded visible x range; use `fit_range=:data` or an explicit `xgrid` at a
physical domain boundary.

The three maintained styles are `:lab`, `:modern`, and `:article`. Explicit
Makie keywords override style defaults only for the element receiving them.
See [Plotting And Customization](plotting_design.md) for the complete layout and
extension contract.

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
