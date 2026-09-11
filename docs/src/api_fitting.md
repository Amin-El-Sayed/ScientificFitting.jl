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
| `fit_likelihood_model` | Supply `logprob(y, prediction, p)`, or fixed additive `error` distributions. A multivariate error object models the residual vector jointly. |
| `fit_distribution` | `make_distribution(p)` builds one upstream distribution per objective evaluation; scalar events use a vector, multivariate events a matrix with explicit `obsdim`. |
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

### Minimization And Local Errors

All likelihood helpers accept these independent controls; Gaussian `fit_model`
uses `scale_covariance`. Both families accept `solver` objects as an alternative
to their legacy `backend` or `optimizer` shortcuts.

| Keyword | Choices and behavior |
|---|---|
| `optimizer` | `:auto` selects LBFGS, or IPNewton for nonlinear constraints. Explicit `:lbfgs`, `:ipnewton`, and `:nelder_mead` are available. |
| `parameter_covariance` | `:auto` selects `:none` with derivative-free solvers, `:hessian` otherwise. `:hessian` requires a locally smooth cost; `:none` leaves free-parameter errors as `NaN` and preserves explicitly supplied fixed-parameter errors. |

Nelder-Mead uses NLopt's native box bounds without numerical derivatives or a
custom penalty. Fixed values, Gaussian priors and correlated parameter terms
remain active. Nonlinear constraints require a capable solver, such as IPNewton;
incompatible methods reject them, never ignore them. All methods optimize continuous parameters and
are local searches. Begin at finite cost inside the likelihood's support.

For Nelder-Mead, `maxiters` is an **objective-evaluation budget**, not an
iteration count; `result.iterations` is `missing`. `tol` sets absolute/relative
parameter stopping tolerances, so choose parameter units/scales accordingly.
Function-value stopping is disabled: distant simplex vertices can have equal
costs without locating the minimum. Reaching the budget is
not convergence. Profiles preserve both controls; use explicit `values` grids
when no local errors exist. [Non-regular likelihoods](likelihood_models.md#A-Moving-Support-Boundary)
need more than a successful minimization to justify confidence intervals.

```@docs
ScientificFitting.fit(::ScientificFitting.LikelihoodFitProblem)
ScientificFitting.fit_custom
ScientificFitting.fit_likelihood_model
ScientificFitting.fit_distribution
ScientificFitting.fit_poisson_model
ScientificFitting.fit_histogram_model
ScientificFitting.fit_histogram_density
ScientificFitting.fit_unbinned_model
ScientificFitting.fit_extended_unbinned_model
ScientificFitting.fit_indexed_model
ScientificFitting.fit_multi_model
ScientificFitting.LikelihoodFitProblem
```

### Distribution Objects

Choose the interface according to what is random. These are different models,
not different numerical solvers:

```@example distribution_objects
using ScientificFitting, Distributions

# Fit an event distribution: both its location and scale are unknown.
events = [-0.5, 0.2, 0.8, 1.1]
make_distribution(p) = Normal(p[1], exp(p[2]))
result = fit_distribution(make_distribution, events;
    p0=[0., 0.], parameter_names=["mean", "log_scale"])
println(report_text(result))
```

```@example distribution_errors
using ScientificFitting, Distributions

# Fit a response curve: the additive measurement-error distributions are known.
x, y = [0., 1., 2., 3.], [0.1, 1.2, 1.8, 3.4]
sigma = [0.2, 0.3, 0.2, 0.4]
line(x, p) = @. p[1]*x + p[2]
result = fit_likelihood_model(line, x, y;
    error=Normal.(0., sigma), p0=[1., 0.])
println(report_text(result))
```

A single `error=TDist(4)` applies the same independent Student-t error to each
residual. `error=MvNormal(zeros(length(y)), covariance)` instead evaluates one
joint log density; it does not multiply marginal errors. For a Poisson or other
prediction-dependent observation model, retain the `logprob` callback.

Factories may return DistributionsHEP shapes or NumericalDistributions objects.
The latter's normalization runs once when the factory constructs the distribution,
then is reused for all observations. ScientificFitting differentiates that
normalization too; a constructor that cannot accept dual numbers needs explicit
`derivatives=:finite`. A `pdf`-only object uses `log(pdf)`, with the tail precision
of its implementation. No interpolated or moving-support model is assumed smooth.

With `using BuildConstructors`, `fit_distribution(constructor, events)` takes
names, starts, bounds and fixed/shared state from the constructor. Optional
`p0=(mu=0.2,)` overrides free starting values. Change bounds/fixed state on the
constructor, not through duplicate fit keywords. Metadata `uncertainty` values
do not create priors or statistical errors. The model is snapshotted for fitting
and profiling; `fitted_model(result)` reconstructs the fitted distribution or
intensity. `parameter_values(result)` provides named fitted values for an explicit
`update!` of the original constructor. The integration does not require Minuit.

```@docs
ScientificFitting.fitted_model
```

## Solver Adapters

`solver=OptimizationSolver(algorithm; native_options...)` accepts algorithms
from the corresponding Optimization.jl solver packages. For example:

```@example solver_choice
using ScientificFitting, OptimizationOptimJL

# Same data likelihood; only the numerical minimizer is selected here.
cost(p) = (p[1] - 2)^2 + (p[2] + 1)^2 / 4
result = fit_custom(cost; p0=[0., 0.], nobs=10,
                    solver=OptimizationSolver(BFGS()))
println(report_text(result))
```

With the optional NativeMinuit package installed, use
`import NativeMinuit` and `solver=NativeMinuitSolver(steps=[0.2, 0.3])` instead.
`import` activates the extension without importing NativeMinuit's own `profile`
name; use `ScientificFitting.profile` for our result-based refits and
`NativeMinuit.profile` for its native API.
NativeMinuit requires Julia 1.11 or later; the core still supports Julia 1.10.
It is an LGPL-licensed dependency, not bundled or copied into ScientificFitting.
Its MIGRAD adapter supports box bounds, fixed parameters and all statistical
parameter terms, but rejects nonlinear equality/inequality constraints.

`steps` contains numerical initial step sizes in **full parameter order**, not
measurement uncertainties. `tol` is Minuit's EDM tolerance; `maxiters` is its
function-call budget. A single MIGRAD pass runs per multistart candidate; the
adapter does not hide additional retries. Other native constructor options,
such as `strategy=2`, are retained by profile refits. No solver setting changes
the objective's ``\chi^2``/``-2\log L`` scale (`errordef=1`).

`result.solver_result.raw` exposes the native solver result. For NativeMinuit,
this is the `Minuit` object for native HESSE/MINOS/contour operations. Its vector
order is `result.solver_result.parameter_indices`; parameters already fixed by
ScientificFitting are absent. Mutating that object does not update the stored
ScientificFitting result. Local errors in `param_covariance` still follow
ScientificFitting's covariance policy, not an implicit replacement by MINOS.
For native MINOS results, inspect validity and parameter-limit flags: reaching
a bound is not finding a likelihood-threshold crossing.

Named likelihood problems carry their unique, nonempty `parameter_names` into
the native object, including during reduced nuisance-parameter refits.

Third-party adapters implement two methods: capabilities and `solve_fit`.
They receive a standard `OptimizationProblem` with the complete objective and
free-coordinate constraints. The statistical result is constructed by the
shared core. See [Backend Design](backend_design.md#The-Solver-Extension-Boundary).

```@docs
ScientificFitting.AbstractFitSolver
ScientificFitting.OptimizationSolver
ScientificFitting.NativeMinuitSolver
ScientificFitting.FitSolverResult
ScientificFitting.solver_capabilities
ScientificFitting.solve_fit
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
