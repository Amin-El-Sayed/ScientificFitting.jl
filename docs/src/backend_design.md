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
non-smooth/domain-limited callbacks still need care. Differentiated models must be evaluable
in a neighborhood of each evaluation point, including near declared bounds.
For pointwise x-error propagation, four vectorized model calls give a central
fourth-order estimate of every `df/dx`, rather than one Python call per
observation. The step ``h_i = \epsilon^{1/5}\max(|x_i|,1)`` balances truncation
and roundoff; the larger step also reduces noise when the likelihood gradient
differentiates this estimate again. Analytic `x_derivative` callbacks bypass it.

The Python bridge converts inputs once and enters the selected fit through
`invokelatest`, an inference boundary outside the numerical loop. This prevents
compiling all solver branches through Python's dynamic call dispatcher. Foreign
callbacks use a private `_TypedCallback{R}`: the result type is concrete, while
one dynamic dispatch calls the runtime-created closure. Different Python models
therefore share Julia solver specializations. Native Julia models bypass it.
The Python adapter declares the output rank: models and x-derivatives return
vectors, allocating parameter Jacobians return matrices. Input validation must
preserve that distinction before Julia conversion.
A small `PrecompileTools` workload caches Gaussian/Poisson numerical kernels,
prediction, and reporting with ordinary Julia callbacks. It neither starts
Python nor loads plotting during installation or package import.

Event and histogram density helpers accept `vectorized=true` without changing
their scalar default. Event log likelihoods use one density call per objective;
integrals use QuadGK's `BatchIntegrand`, retaining adaptive error control for
each bin. Batch buffers retain the parameter element type, including dual
numbers, and are reused across bins within one objective evaluation. Profile
refits retain the same density closure. An analytic bin integral can instead
be supplied through `fit_histogram_model`, avoiding quadrature entirely.

Solver selection follows the represented problem rather than a speed preference:

| Condition | Backend |
|---|---|
| All parameters fixed | no optimizer; evaluate the complete result once |
| Unbounded, static Gaussian chi-square without extra parameter terms | LsqFit least-squares path |
| Bounds, priors, parameter constraints, parameter-dependent covariance, or likelihood objective | Optimization.jl with LBFGS |
| Nonlinear equality or inequality constraints | Optimization.jl with IPNewton |
| Likelihood with explicit `optimizer=:nelder_mead` and no nonlinear constraints | OptimizationNLopt with native bounded Nelder-Mead, without derivatives |

An explicit `backend=:lsqfit` request is rejected if it would discard any part
of the statistical problem. Backend selection may change how the same objective
is minimized; it must never change which objective is being minimized.

Likelihood `optimizer` and `parameter_covariance` are independent options stored
in `FitOptions` and retained by profile refits. Nelder-Mead defaults to no
Hessian calculation, not a numerical Hessian across a kink or support boundary.
It uses the same objective/cache and parameter controls, without a separate
statistics implementation or bound penalty. Its `maxiters` is NLopt's objective
evaluation budget; the unavailable iteration count remains `missing`.

Both fit families use one multistart ranking rule: a converged finite result
outranks an unconverged one; within the same status, the lower cost wins. If
every run stops early, return the best finite result with `converged=false`,
not merely the first candidate.

Both solver objectives retain the model/objective type even when their evaluation
cache is passed as solver context. ForwardDiff tags must not be shared between
precompiled fits and unseen models with nested x derivatives
([ForwardDiff #714](https://github.com/JuliaDiff/ForwardDiff.jl/issues/714)).
The startup gate compares automatic and analytic x derivatives in a fresh
process, for diagonal and dense x covariance, and checks a custom objective
with a nested derivative against its exact minimum and covariance. Python's
typed finite-difference callbacks still share their precompiled solver path.

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
Hessian on the ``-2\log L`` scale. `LikelihoodFitResult` construction evaluates
it once and reuses it for diagnostics.
For likelihoods with `parameter_covariance=:none`, free errors/covariances are
`NaN`, fixed errors remain zero, and diagnostics explain the omission rather
than claiming a curvature failure. Covariance conditioning uses free coordinates;
a fixed parameter's zero variance is not a singularity. This covariance is a local approximation;
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
| Batched densities, adaptive integration, derivatives, and scalar parity | `test/statistics/vectorized_density_reference.jl` |
| Profiles, contours, local approximations, and failed refits | `test/statistics/profile_contour_reference.jl` |
| Structured matrix-free covariance | `test/statistics/structured_whitening_reference.jl` |
| In-place models and Jacobians | `test/numerics/inplace_model_reference.jl` |
| Solver limits, convergence status, and stopped profile refits | `test/numerics/solver_status_reference.jl` |
| Laplace median, moving support, derivative-free nuisance refits, and optional Hessian errors | `test/numerics/nonsmooth_likelihood_reference.jl` |
| Invalid scientific and numerical inputs | `test/numerics/torture_inputs.jl` |
| Public compatibility and optional plotting boundary | `test/regression/current_api.jl` |
| Steady-state hot-path budgets | `test/performance_budget_gate.jl` |
| Fresh-process loading and nested automatic derivatives | `test/startup_probe_gate.jl` |
| Plot composition and extension behavior | `test/plots/fitplot.jl` |
| NumPy callback parity, native Matplotlib panels, profiles, and ownership | `python/tests` |

The core gate is `julia --project=. test/core_runtests.jl`; the complete package
gate is `julia --project=. test/runtests.jl`. Performance methodology and the
benchmark runner are documented on the [Performance](performance.md) page.

For a page-level output check, run
`julia --project=docs test/docs_output_snapshots.jl gallery/resonance_decay.md`
(additional page paths are accepted). Without page arguments the gate executes
every documented workflow. It compares the displayed output with both the
page's code cells and the example generator; only solver iteration counts are
normalized across solver/platform versions.

## Planned Work

- [x] **v0.2: concrete statistical scope at the entry point.** README and
  documentation entry now distinguish observation models, likelihood
  optimization, local parameter errors, profiles, and model bands from
  posterior sampling. Source checks guard that distinction.
- [x] **User-defined measurement-error distributions.** `fit_likelihood_model`
  now accepts batched log densities/masses in Julia and Python, reusing the
  common likelihood engine. Gaussian, fitted-scale, Binomial, and Student-t
  references cover normalization and curvature. Bounded derivative-free
  Nelder-Mead and independent covariance controls now cover non-smooth and
  moving-support examples without fabricated Hessian errors. Laplace and
  exponential references verify minima and profile costs; worked guidance
  distinguishes successful minimization from valid interval coverage.
- [ ] **v0.2: complete the native Python interface.** The preview now wraps
  every high-level fit family, named parameter controls, sparse/structured
  covariance, in-place callbacks, and core reports/profiles/diagnostics.
  [Python examples](python.md) use NumPy models and native Matplotlib; the
  wrapper does not duplicate statistical algorithms or require Makie.
  Numerical API review and eleven fresh-process cases verify parameters,
  covariance, and cost against analytic references. [Startup and callback
  measurements](performance.md#Python-Startup) distinguish first use from warm
  fitting and record the remaining runtime and quadrature overhead.
  Local installed-wheel provisioning and Conda Python 3.14 reference fits pass
  on macOS ARM64. Wheel/sdist metadata, MIT license, and the 0.2 core pin are
  checked; installed Conda and packaged wheel runtime files match the current
  sources. Registry installation checks reject persisted development/repository
  overrides and record the actual resolved core source and tree hash.
  **Remaining:** run the configured installed-package CI on Linux, macOS, and
  Windows, then repeat clean wheel/Conda installation against the registered
  0.2 core without development overrides before publishing Python packages.
