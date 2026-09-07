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
- sparse covariance keeps its sparse Cholesky factor,
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

Derivative selection is stored in the problem, not in the renderer.
`derivatives=:auto` keeps the Julia defaults; `derivatives=:finite` uses central
finite differences for foreign callbacks that accept ordinary floating-point
values but not ForwardDiff dual numbers. It applies to optimization, constraint
derivatives, post-fit covariance, predictions, and profile/contour refits.
Explicit model Jacobians and x-derivatives take precedence where they apply.
The solver stopping tolerance defaults to `1e-6` in finite mode and `1e-10`
otherwise. This avoids demanding convergence below the noise floor of
numerically differenced gradients; it is not an error bound on fitted
parameters. Explicit `tol` values are preserved, including when they lead to
a reported convergence failure.
For LsqFit, `maxiters` maps to `maxIter` and `tol` to both its step (`x_tol`)
and gradient (`g_tol`) criteria. The result keeps the solver's actual
convergence flag; reaching an iteration limit does not imply success.

The finite mode differentiates the **whole** objective, including any
parameter-dependent covariance and its log determinant. It is an approximation,
not a claim of exact derivatives: noisy models, badly scaled parameters, and
non-smooth/domain-limited callbacks still need care. Models must be evaluable
in a neighborhood of each evaluation point, including near declared bounds.
For pointwise x-error propagation, two vectorized model calls estimate all
`df/dx` values instead of crossing a Python boundary once per observation.

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

CHOLMOD's sparse solves do not accept ForwardDiff dual numbers. Static sparse
covariance therefore requires `derivatives=:finite` with the general optimizer,
or an AD-compatible `WhiteningOperator` instead. The least-squares path remains
available without that override. Python uses finite derivatives consistently,
including sparse fits with bounds and profile refits. Sparse covariance
components are validated through their stored entries, without a dense copy.

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

Automatic scan ranges intersect the local-covariance range with declared bounds;
they do not generate known-infeasible trial values. Explicit grids are preserved,
and missing threshold crossings remain missing rather than being replaced by a
bound. This policy is shared by Gaussian and general-likelihood profiles,
contours, and matrices.

The Python renderer uses native Matplotlib objects, not Makie. Both renderers
share core residual/ratio preparation; observation pulls use stored whitened
residuals without auxiliary parameter terms. Completed profile/contour snapshots
are rendered without further scans. Matplotlib's constrained layout reserves
outside legends and reports; the wrapper does not install a separate layout
engine or resize callbacks.

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
| `derivatives.jl` and `prediction.jl` | shared derivative policy and model-mean uncertainty propagation |
| `likelihood_fits.jl` | likelihood problem construction, wrappers, solver path, and `LikelihoodFitResult` |
| `profile.jl` | fixed-parameter refits, profile intervals, contours, matrix summaries, and their diagnostics |
| `diagnostics.jl` | structured findings, severity, evidence, next actions, and renderer-independent residual values |
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
| Solver limits, convergence status, and stopped profile refits | `test/numerics/solver_status_reference.jl` |
| Invalid scientific and numerical inputs | `test/numerics/torture_inputs.jl` |
| Public compatibility and optional plotting boundary | `test/regression/current_api.jl` |
| Steady-state hot-path budgets | `test/performance_budget_gate.jl` |
| Plot composition and extension behavior | `test/plots/fitplot.jl` |
| NumPy callback parity, native Matplotlib panels, profiles, and ownership | `python/tests` |

The core gate is `julia --project=. test/core_runtests.jl`; the complete package
gate is `julia --project=. test/runtests.jl`. Performance methodology and the
benchmark runner are documented on the [Performance](performance.md) page.

## Planned Work

- [x] **v0.2: concrete statistical scope at the entry point.** README and
  documentation entry now distinguish observation models, likelihood
  optimization, local parameter errors, profiles, and model bands from
  posterior sampling. Source checks guard that distinction.
- [ ] **User-defined measurement-error distributions.** `fit_likelihood_model`
  now accepts batched log densities/masses in Julia and Python, reusing the
  common likelihood engine. Gaussian, fitted-scale, Binomial, and Student-t
  references cover normalization and curvature. Remaining: suitable explicit
  solver/inference controls for non-smooth or support-limited objectives and
  worked guidance beyond the current smooth, independent-observation contract.
- [ ] **v0.2: complete the native Python interface.** The preview now wraps
  every high-level fitting family, including correlated parameter constraints
  and named multi-dataset sharing, sparse/structured covariance, error
  components, and in-place callbacks. Named result/report/diagnostic snapshots,
  asymmetric profile intervals, and profile matrices now preserve core findings
  and missing values. Native Matplotlib reports, Gaussian residual diagnostics,
  and profile/contour matrices are implemented with panel visibility independent
  of typography. Executable Student-t guidance and Poisson, unequal-bin, and
  multi-dataset examples now use ordinary Matplotlib with numerical references
  and export checks. Named additional starts follow `p0` order, and multi-dataset
  errors use the same real-input validation as single-dataset fits.
  Remaining: the final numerical API-parity review, including explicit multistart
  budgets in the Julia gallery (supplying `initial_guesses` alone does not run them).
  Do not add parallel plotting abstractions where native composition suffices. Deliver
  documented Python APIs, pip/conda installation without manual Julia setup,
  and clean-install tests on supported platforms. Measure numerical parity,
  callback overhead (including scalar quadrature), startup time, and installed
  size; automatic Julia provisioning does not eliminate its runtime footprint.
