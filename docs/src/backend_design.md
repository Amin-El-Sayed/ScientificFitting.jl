# Backend Design

This page is the contributor map for ScientificFitting's numerical core. It explains
where statistical meaning is established, where numerical work happens, and
which boundaries must not be blurred. User-facing signatures and defaults are
documented in the [API Reference](api.md); the statistical derivations are in
[Mathematics and Statistics](statistical_foundations.md).

The central dependency direction is deliberately one-way:

```math
\begin{aligned}
\text{public constructor}
&\longrightarrow \text{validated problem} \\
&\longrightarrow \text{objective and compatible solver} \\
&\longrightarrow \text{result} \\
&\longrightarrow \text{reports, diagnostics, profiles, and plots}.
\end{aligned}
```

Later stages may inspect an earlier result. They must not reconstruct or mutate
the statistical problem behind it.

## Two Problem Families

ScientificFitting has two normalized problem types because a residual vector and a
general likelihood do not expose the same information.

| Normalized problem | Scientific payload and public path |
|---|---|
| `FitProblem` -> `FitResult` | Gaussian x-y data built by `fit_model` or `FitProblem`; stores the model, observations, uncertainty, and parameter controls; minimizes static ``\chi^2`` or normalized Gaussian ``-2\log L``. |
| `LikelihoodFitProblem` -> `LikelihoodFitResult` | Counts, samples, indexed data, or custom objectives built by the likelihood helpers or `fit_custom`; stores the objective, optional goodness statistic, observation count, and parameter controls. |

Both families share parameter bounds, fixed parameters, Gaussian parameter
terms, nonlinear constraints, multistart selection, local covariance,
diagnostics, and profile refits. They stay separate where their data contracts
differ: a generic likelihood need not have x-y residuals, model predictions, or
a natural fit curve.

## One Fit, Step By Step

### 1. Normalize and validate the scientific input

Convenience functions construct a `FitProblem` or `LikelihoodFitProblem` before
optimization begins. Problem construction and the public fit entry points copy
numeric inputs into stable storage and reject mismatched dimensions, non-finite
observations or starting values, non-positive standard deviations, invalid
covariance matrices, inconsistent parameter indices, and fixed values outside
declared bounds before solver dispatch.

This boundary is intentional. A solver should never be asked to discover that a
covariance matrix is not positive definite or that the model returned the wrong
number of predictions.

### 2. Map full parameters to optimizer coordinates

The scientific model always sees the complete parameter vector. Fixed
parameters are removed only from the optimizer-visible vector:

```math
q_{\mathrm{free}}
\xrightarrow{\text{expand}}
p_{\mathrm{full}}.
```

Bounds are reduced to the same free coordinates. Nonlinear constraint callbacks
are wrapped so that user code still receives `p_full`. After fitting, free
covariance and Jacobian blocks are embedded back into full parameter order.
This single mapping is reused by ordinary fits, multistart candidates, profiles,
and contours.

### 3. Prepare reusable evaluation state

Static Gaussian uncertainty is prepared outside repeated scalar objective
evaluations:

- diagonal errors become inverse standard deviations and a log determinant,
- dense covariance becomes a Cholesky factor and a log determinant,
- sparse covariance stays sparse on the compatible least-squares path,
- `WhiteningOperator` keeps the supplied matrix-free operation and determinant,
- correlated parameter constraints are factorized once.

The general Gaussian path stores this state in `FitEvaluationCache`; likelihood
fits use `LikelihoodEvaluationCache` for reusable parameter-constraint state.
The LsqFit path prepares equivalent static weights directly for its native
residual interface.

Parameter-dependent covariance is not cached as if it were static. Effective
x-error covariance and model-relative uncertainty must be recomputed at each
parameter point because that dependence is part of the probability model.

### 4. Construct exactly one objective

For a Gaussian problem, `cost=:auto` selects static chi-square when the
covariance is parameter independent:

```math
\chi^2(p)=r(p)^\mathsf{T}V^{-1}r(p).
```

When the effective covariance depends on the fitted parameters, it selects the
normalized Gaussian objective:

```math
-2\log L(p)
=n\log(2\pi)+\log\det V(p)+r(p)^\mathsf{T}V(p)^{-1}r(p).
```

The log determinant cannot be dropped in the second case. Doing so changes the
optimum, not merely the reported normalization.

Likelihood problems provide their data objective directly on the ``-2\log L``
scale. Gaussian parameter priors and correlated parameter constraints are then
added by the shared parameter layer. Bounds and fixed parameters restrict the
parameter space; they are not hidden penalty terms.

### 5. Dispatch only to a compatible solver

Solver selection follows the represented problem rather than a speed preference:

| Condition | Backend |
|---|---|
| All parameters fixed | no optimizer; evaluate the complete result once |
| Unbounded, static Gaussian chi-square without extra parameter terms | LsqFit least-squares path |
| Bounds, priors, parameter constraints, parameter-dependent covariance, or likelihood objective | Optimization.jl with LBFGS |
| Nonlinear equality or inequality constraints | Optimization.jl with IPNewton |

An explicit `backend=:lsqfit` request is rejected if it would discard any part
of the statistical problem. Backend selection may change how the same objective
is minimized; it must never change which objective is being minimized.

### 6. Build the result once

`FitResult` and `LikelihoodFitResult` are the numerical source of truth. Result
construction records the selected minimum, solver status, parameter estimates,
local parameter covariance, correlations, statistics, and diagnostics. Gaussian
x-y results additionally retain model predictions, raw residuals, weighted
residuals, and the weighted Jacobian.

For static least squares, local covariance comes from the weighted Jacobian. For
Gaussian likelihood and general likelihood fits, it comes from the objective
Hessian on the ``-2\log L`` scale. This covariance is a local approximation;
profiles and contours remain separate refit operations when the cost is not
locally quadratic.

### 7. Read the result without changing it

`report_text`, `diagnose`, and `diagnostic_dashboard` consume result fields.
`profile` and `contour` reuse the stored normalized problem, fix one or two
parameters, and refit the remaining nuisance parameters. Plotting is an optional
CairoMakie package extension and consumes the same result objects. None of these
paths reruns or alters the original fit unless the API explicitly describes a
profile or contour refit.

## Numerical Invariants

These are architectural rules, not implementation preferences.

| Invariant | Consequence |
|---|---|
| Statistical semantics precede solver choice. | An optimizer cannot silently drop bounds, priors, constraints, or covariance terms. |
| Covariance is applied by factorization, solves, or a validated whitening operation. | Production cost evaluation does not form an explicit covariance inverse. |
| Numerical repairs are visible. | Invalid inputs fail; ScientificFitting does not add hidden diagonal jitter to make a covariance appear usable. |
| Static work stays outside the hot objective. | Repeated evaluations reuse factors, determinants, and prepared constraint state. |
| Full and free parameter order have one mapping. | Fixed parameters, callbacks, covariance dimensions, ndf, profiles, and reports remain consistent. |
| ``-2\log L`` is the likelihood scale. | Hessian covariance, likelihood-ratio thresholds, AIC, and BIC use one convention. |
| A local covariance is not a coverage guarantee. | Diagnostics expose suspect curvature; profile and contour results remain first-class outputs. |
| Plotting is optional. | `using ScientificFitting` provides fitting, reports, diagnostics, and profiles without loading Makie. |

These invariants are the review contract for changes to the numerical core.

## Source Map

The core is split by responsibility rather than by feature-specific vertical
stacks.

| Source | Owns |
|---|---|
| `types.jl` | validated problem/result types, uncertainty inputs, and in-place wrappers |
| `parameters.jl` | full/free parameter mapping, fixed parameters, bounds, and multistart candidates |
| `weights.jl` | covariance preparation, whitening, weighted residuals/Jacobians, local covariance helpers, backend compatibility |
| `costs.jl` | chi-square, normalized Gaussian likelihood, priors, and correlated parameter terms |
| `fit.jl` | Gaussian solver dispatch and `FitResult` construction |
| `likelihood_fits.jl` | likelihood problem construction, wrappers, solver path, and `LikelihoodFitResult` |
| `profile.jl` | fixed-parameter refits, profile intervals, contours, matrix summaries, and their diagnostics |
| `diagnostics.jl` | structured findings, severity, evidence, and next actions |
| `report.jl` | Makie-free report objects and text formatting |
| `plotting_api.jl` | public plotting boundary and informative fallback methods |
| `ext/ScientificFittingCairoMakieExt.jl` plus `plotting.jl` | CairoMakie rendering only |

This map is also a review rule. For example, a plotting feature should not add a
second statistical calculation, and a new optimizer should not own covariance
semantics.

## Where A New Feature Belongs

Before adding a type or abstraction, first ask whether an existing problem can
already express the required statistics.

| Change | Preferred integration |
|---|---|
| New convenience fitting function | validate its domain-specific inputs, then construct an existing problem type |
| New static uncertainty representation | add validation, preparation/whitening, determinant semantics, and an analytic covariance reference |
| New likelihood family | provide a ``-2\log L`` objective, a justified goodness statistic when one exists, and an explicit observation count |
| New numerical backend | add a compatibility predicate and prove that the represented objective is unchanged |
| New diagnostic | consume a result or profile object and return structured evidence plus an action |
| New report or plot | consume existing result fields; keep rendering inside the optional extension |

Do not add a parallel result type, cache, or solver path merely to support a new
presentation. Small APIs that compose existing contracts are easier to audit
than duplicated feature stacks.

## Verification Map

Architecture changes need evidence at the layer they affect:

| Claim | Primary evidence |
|---|---|
| Gaussian values, covariance, normalization, and constraints | `test/statistics/linear_gaussian_reference.jl` and `covariance_semantics_reference.jl` |
| Poisson, histogram, unbinned, extended, indexed, and multi-fit semantics | `test/statistics/likelihood_reference.jl` |
| Profiles, contours, local approximations, and failed refits | `test/statistics/profile_contour_reference.jl` |
| Structured matrix-free covariance | `test/statistics/structured_whitening_reference.jl` |
| In-place models and Jacobians | `test/numerics/inplace_model_reference.jl` |
| Invalid scientific and numerical inputs | `test/numerics/torture_inputs.jl` |
| Public compatibility and optional plotting boundary | `test/regression/current_api.jl` |
| Steady-state hot-path budgets | `test/performance_budget_gate.jl` |
| Plot composition and extension behavior | `test/plots/fitplot.jl` |

The core gate is `julia --project=. test/core_runtests.jl`; the complete package
gate is `julia --project=. test/runtests.jl`. Performance methodology and the
benchmark runner are documented on the [Performance](performance.md) page.
