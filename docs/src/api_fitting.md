# Fitting

This page defines the fitting inputs and observation-model contracts. Shared
parameter ordering and solver conventions are listed in the
[API overview](api.md).

## Gaussian Fits

### Observation Uncertainty

Choose one representation for each physical uncertainty source; contradictory
combinations are rejected.

| Keyword | Accepted value | Statistical role |
|---|---|---|
| `sigma_y` | positive vector | Independent y standard deviations. |
| `cov_y` | dense or sparse SPD matrix | Complete y covariance. Mutually exclusive with `sigma_y`. |
| `sigma_x` | positive vector | Independent x standard deviations propagated through ``\partial f/\partial x``. |
| `cov_x` | dense or sparse SPD matrix | Complete x covariance. Mutually exclusive with `sigma_x`. |
| `error_components` | named [`ErrorComponent`](@ref)s | Additive absolute, relative, model-relative, or covariance contributions. |
| `whitening` | [`WhiteningOperator`](@ref) | Complete static covariance represented by ``W^\mathsf{T}W=C^{-1}``. |

`whitening` is exclusive with every other observation-uncertainty keyword: it
already represents the complete covariance
([Structured Whitening](statistics.md#Structured-Whitening)).

With no supplied observation uncertainty, `fit_model` performs unweighted
least squares; the `scale_covariance` policy is specified in
[Covariance Scaling](statistics.md#Covariance-Scaling).

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

Non-finite observations, non-positive standard deviations, invalid bounds, and
non-positive-definite covariance matrices raise `ArgumentError` or
`DomainError` before
optimization.

Multistart fits return the candidate with the lowest finite cost — convergence
status only breaks exact ties — reported with its own `converged` flag;
inspect the status or use [`diagnostic_dashboard`](@ref). If every candidate
fails, the underlying error is raised.

```@docs
ScientificFitting.fit(::ScientificFitting.FitProblem)
ScientificFitting.fit_model
ScientificFitting.FitProblem
ScientificFitting.FitOptions
```

## Likelihood And Count Fits

Poisson, histogram, unbinned, and extended-unbinned entry points minimize
costs on the ``-2\log L`` scale. Poisson and histogram fits fill `chi2`,
`chi2_ndf`, and `pvalue` from the
[Poisson deviance](statistics.md#Poisson-Counts-And-Histograms);
[ordinary and extended unbinned fits](statistics.md#Unbinned-And-Extended-Likelihoods),
and any fit with an exactly zero Poisson expectation, leave those fields `NaN`.

`fit_indexed_model` and `fit_multi_model` minimize chi-square but omit
additive Gaussian normalization constants, so their AIC/BIC compare only
models fit to the same observations with the same uncertainty model.

| Entry point | Additional contract |
|---|---|
| `fit_likelihood_model` | Supply `logprob(y, prediction, p)`, or fixed additive `error` distributions. A multivariate error object models the residual vector jointly. |
| `fit_distribution` | Builds one upstream distribution per objective evaluation. Accepts events (multivariate matrices need `obsdim`) or `edges, counts` with an explicit expected total or extended component yields. |
| `fit_poisson_model` | Every expected count must be finite and nonnegative; observed counts must be non-negative integers. |
| `fit_histogram_model` | `length(edges) == length(counts) + 1`; edges increase strictly; the model returns one nonnegative expectation per bin. |
| `fit_histogram_density` | Integrates `pdf(x, p)` over every bin with Gauss-Kronrod quadrature; `total_count > 0`, `rtol > 0`. |
| `fit_unbinned_model` | The supplied density must already be normalized and positive at every observation. |
| `fit_extended_unbinned_model` | `rate` is an intensity, not a density; its integral over `domain` is the expected event count. |
| `fit_indexed_model` | Supports `sigma_y` or `cov_y`; indices may be any container accepted by the model. |
| `fit_multi_model` | Supports per-dataset `sigma_y`; `parameter_map[i]` selects global parameters passed to model `i`. |

For `fit_custom`, `objective` should be a normalized ``-2\log L`` cost
([The Cost Convention](statistics.md#The-Cost-Convention)); with an
arbitrarily scaled loss, local covariance, AIC, and BIC are only arithmetic
summaries. `nobs` must count statistically independent observations. An
optional `gof(p)` supplies the data goodness-of-fit statistic; Gaussian
parameter priors and constraints contribute as in
[Degrees Of Freedom](statistics.md#Degrees-Of-Freedom).

### Minimization And Local Errors

All likelihood helpers accept these independent controls; Gaussian `fit_model`
uses `scale_covariance` instead.

| Keyword | Choices and behavior |
|---|---|
| `optimizer` | `:auto` selects LBFGS, or IPNewton for nonlinear constraints. Explicit `:lbfgs`, `:ipnewton`, and `:nelder_mead` are available. |
| `parameter_covariance` | `:auto` selects `:none` with derivative-free solvers, `:hessian` otherwise. `:hessian` requires a locally smooth cost; `:none` leaves free-parameter errors as `NaN` and preserves explicitly supplied fixed-parameter errors. |

Nelder-Mead uses NLopt's native box bounds without numerical derivatives or a
custom penalty; fixed values, Gaussian priors, and correlated parameter terms
remain active. Incompatible solvers reject nonlinear constraints, never ignore
them. All methods are local searches over continuous parameters; begin at
finite cost inside the likelihood's support.

For Nelder-Mead, `maxiters` is an **objective-evaluation budget**,
`result.iterations` is `missing`, and reaching the budget is not convergence.
`tol` sets absolute/relative parameter stopping tolerances, so choose
parameter units accordingly; function-value stopping is disabled. Profiles
preserve both controls; use explicit `values` grids when no local errors
exist. [Non-regular likelihoods](statistics.md#Observation-Likelihoods) need
more than a successful minimization to justify confidence intervals.

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

Choose the interface according to what is random:

```@example distribution_objects
using ScientificFitting, Distributions

# Fit an event distribution: both its location and scale are unknown.
events = [-0.5, 0.2, 0.8, 1.1]
make_distribution(p) = Normal(p[1], exp(p[2]))
result = fit_distribution(make_distribution, events;
    p0=[0., 0.], parameter_names=["mean", "log_scale"])
println(report_text(result))
```

For binned events, pass edges and counts instead. The following example treats
50 events as an independently known expected full-support yield, **not** a
number estimated from the observed histogram, whose window may omit events; an
estimated total needs a fitted yield parameter, not `total_count=sum(counts)`
treated as known.

```@example distribution_bins
using ScientificFitting, Distributions

edges = [-2., -0.8, 0.2, 1.5, 3.]
counts = [5, 15, 18, 7]
# Fit the peak location; the width and expected full-support count are known.
result = fit_distribution(p -> Normal(p[1], 1.), edges, counts;
    p0=[0.], total_count=50., parameter_names=["mean"])
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

Model-construction functions may return DistributionsHEP or
NumericalDistributions objects; the latter's normalization runs once per
constructed distribution, is reused for all observations, and is
differentiated with the model.

With `using BuildConstructors`, `fit_distribution(constructor, events)` takes
names, starts, bounds, and fixed/shared state from the constructor
([Named Model Construction](interfaces.md#Named-Model-Construction)).

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

With the optional NativeMinuit package installed
([Minimizers](interfaces.md#Minimizers)), `import NativeMinuit` — keeping its
own `profile` name out of scope — and pass
`solver=NativeMinuitSolver(steps=[0.2, 0.3])`. The MIGRAD adapter supports box
bounds, fixed parameters, and all statistical parameter terms, but rejects
nonlinear equality/inequality constraints.

`steps` contains numerical initial step sizes in **full parameter order**, not
measurement uncertainties. `tol` is Minuit's EDM tolerance, defaulting to
`0.1`; `maxiters` is the requested MIGRAD function-call budget —
gradient/covariance checks are additional. A failed or rejected attempt may
restart once from the returned point within the remaining budget, with
unchanged tolerance and strategy; diagnostics record the restart. Native
constructor options such as `strategy=2` are retained by profile refits. No
solver setting changes the objective's ``\chi^2``/``-2\log L`` scale
(`errordef=1`).

For smooth interior fits with positive local covariance, SF also checks
``g^{\mathsf T}\operatorname{Cov}(\hat p)g/4``: the estimated remaining cost
decrease using a fresh gradient and SF's curvature. It follows Minuit's
acceptance limit, ten times the nominal EDM goal ``0.002\,\mathrm{tol}``; a
failed check sets `converged=false` and reports `not_stationary`. The check
does not establish a global minimum and is not applied at active bounds, with
nonlinear constraints, or without usable local curvature.

`result.solver_result.raw` exposes the native solver result — for
NativeMinuit, the `Minuit` object for native HESSE/MINOS/contour operations.
Its vector order is `result.solver_result.parameter_indices`; parameters
already fixed by ScientificFitting are absent, and mutating the object does
not update the stored result. `param_covariance` keeps ScientificFitting's
covariance policy. For native MINOS results, inspect validity and
parameter-limit flags: reaching a bound is not finding a likelihood-threshold
crossing.

Named likelihood problems carry their unique, nonempty `parameter_names` into
the native object, including during reduced nuisance-parameter refits.

Third-party adapters implement `solver_capabilities` and `solve_fit`, and may
specialize `default_fit_tolerance(solver, derivatives)`; the extension
contract is specified in
[Backend Design](backend_design.md#The-Solver-Extension-Boundary).

```@docs
ScientificFitting.AbstractFitSolver
ScientificFitting.OptimizationSolver
ScientificFitting.NativeMinuitSolver
ScientificFitting.FitSolverResult
ScientificFitting.solver_capabilities
ScientificFitting.default_fit_tolerance
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
