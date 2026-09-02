# Fitting

This page defines the fitting inputs and observation-model contracts. Shared
parameter ordering and solver conventions are listed in the
[API overview](api.md).

## Gaussian Fits

### Observation Uncertainty

Choose one representation for each physical uncertainty source. ScientificFitting
rejects contradictory combinations instead of guessing how they combine.

| Keyword | Accepted value | Statistical role |
|---|---|---|
| `sigma_y` | positive vector | Independent y standard deviations. |
| `cov_y` | dense or sparse SPD matrix | Complete y covariance. Mutually exclusive with `sigma_y`. |
| `sigma_x` | positive vector | Independent x standard deviations propagated through ``\partial f/\partial x``. |
| `cov_x` | dense or sparse SPD matrix | Complete x covariance. Mutually exclusive with `sigma_x`. |
| `error_components` | named [`ErrorComponent`](@ref)s | Additive absolute, relative, model-relative, or covariance contributions. |
| `whitening` | [`WhiteningOperator`](@ref) | Complete static covariance represented by ``W^\mathsf{T}W=C^{-1}``. |

`whitening` is intentionally exclusive with every other observation-uncertainty
keyword. It describes the complete covariance; adding another source without an
explicit derivation would double-count uncertainty.

With no supplied observation uncertainty, `fit_model` performs unweighted least
squares and `scale_covariance=:auto` estimates residual scale from
``\chi^2/\mathrm{ndf}``. With physical uncertainties, `:auto` leaves their scale
unchanged.

### Error Components

An error component has a stable name and can be activated or deactivated without
rewriting the fit:

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

Non-finite observations, non-positive standard deviations, invalid bounds, and
non-positive-definite covariance matrices raise `ArgumentError` before
optimization.

For multistart fits, ScientificFitting returns the converged candidate with the lowest
finite cost. If no candidate converges but one returns a finite result, it is
returned with `converged == false`; inspect the status or use
[`diagnostic_dashboard`](@ref). If every candidate fails, the underlying error
is raised.

```@docs
ScientificFitting.fit(::ScientificFitting.FitProblem)
ScientificFitting.fit_model
ScientificFitting.FitProblem
ScientificFitting.FitOptions
```

## Likelihood And Count Fits

Poisson, histogram, unbinned, and extended-unbinned entry points minimize costs
on the ``-2\log L`` scale. Poisson and histogram fits also compute Poisson
deviance, so `chi2`, `chi2_ndf`, and `pvalue` are available as goodness-of-fit
summaries. Ordinary and extended unbinned fits do not invent a chi-square
statistic; those fields are `NaN`.

`fit_indexed_model` and `fit_multi_model` minimize chi-square but omit additive
Gaussian normalization constants. Their AIC/BIC values may compare models fit
to the same observations with the same uncertainty model; they must not compare
different uncertainty scales or datasets.

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
observations. If supplied, `gof(p)` is the data goodness-of-fit statistic;
ScientificFitting adds quadratic contributions and dimensions from Gaussian parameter
priors and constraints.

```@docs
ScientificFitting.fit(::ScientificFitting.LikelihoodFitProblem)
ScientificFitting.fit_custom
ScientificFitting.fit_poisson_model
ScientificFitting.fit_histogram_model
ScientificFitting.fit_histogram_density
ScientificFitting.fit_unbinned_model
ScientificFitting.fit_extended_unbinned_model
ScientificFitting.fit_indexed_model
ScientificFitting.fit_multi_model
ScientificFitting.LikelihoodFitProblem
```

## Constraints And Uncertainty Objects

```@docs
ScientificFitting.ConstraintSpec
ScientificFitting.ParameterPrior
ScientificFitting.FixedParameter
ScientificFitting.ParameterConstraint
ScientificFitting.ErrorComponent
ScientificFitting.WhiteningOperator
```
