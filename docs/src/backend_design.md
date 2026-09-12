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
LsqFit and Optimization stopping tolerances default to `1e-6` in finite mode
and `1e-10` otherwise. NativeMinuit instead keeps its native EDM tolerance
(`0.1`); equal tolerance numbers do not mean equal accuracy across solvers.
The finite-difference defaults avoid demanding convergence below the noise floor of
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

### The Solver Extension Boundary

Both problem families enter `src/solvers.jl` for scalar minimization; the
specialized LsqFit residual path remains in `src/fit.jl`. An explicit
`solver::AbstractFitSolver` overrides automatic numerical selection, not the
statistical model. An adapter implements:

- `solver_capabilities`: support for bounds/constraints and required derivatives.
- `solve_fit`: minimize the prepared `OptimizationProblem` and return
  `FitSolverResult` in its free coordinates, with native result and actual status.
- Optionally `default_fit_tolerance`: supply the solver's stopping tolerance
  when the user omits `tol`; the core stores it before multistart or profiling.

The core owns objective construction, parameter mapping, capability checks,
derivative policy, multistart ranking, statistical summaries and profile refits.
Adapters do not add priors, reinterpret error scales or silently drop constraints.
Unknown iteration counts stay `missing`. Solver settings persist when a profile
changes the number of free parameters, so NLopt adapters take algorithm enums
rather than dimension-bound native `Opt` instances.

`OptimizationSolver` uses SciML's algorithm traits and termination codes.
`ScientificFittingNativeMinuitExt` loads only with NativeMinuit; it supplies
the same complete cost and derivative policy to MIGRAD with `errordef=1`.
Its native result retains covariance and failure evidence without replacing
ScientificFitting's covariance policy. Usage and budget conventions are in
[Solver Adapters](api_fitting.md#Solver-Adapters).

Native state is retained behind a non-specializing result field. Reports,
diagnostics and plots share the same result type across user models and solvers;
they do not need a separate compilation for every native solution type. The
objective and solver kernels still specialize on their computational inputs.

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
| `solvers.jl` and `ext/ScientificFittingNativeMinuitExt.jl` | scalar solver contract, capabilities and optional MIGRAD adapter |
| `derivatives.jl` and `prediction.jl` | shared derivative policy and model-mean uncertainty propagation |
| `likelihood_fits.jl` | likelihood problem construction, wrappers, solver path, and `LikelihoodFitResult` |
| `distribution_fits.jl` and `ext/ScientificFittingDistributionsHEPExt.jl` | upstream probability objects, normalization, histogram integrals and extended yields |
| `ext/ScientificFittingBuildConstructorsExt.jl` | named constructor metadata, parameter mapping and fitted-model reconstruction |
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

## [Ecosystem Integrations](@id v03-ecosystem)

[Packages and Interfaces](interfaces.md) contains executable examples of
distribution fitting, named model construction, interchangeable solvers and
posterior inference. The [Python interface](python.md) keeps NumPy callbacks
and Matplotlib plots; it does not wrap Julia constructor or solver objects.

| Interface | Contract | Focused reference tests |
|:---|:---|:---|
| Distribution objects | Reuse upstream `logpdf`/`pdf`/`cdf`; continuous, discrete and joint observations, binned and extended likelihoods | `test/extensions/distribution_ecosystem.jl`, `composite_distributions.jl` |
| BuildConstructors | Preserve names, starts, bounds, fixed/shared parameters and reconstruction; independent of Minuit | `test/extensions/buildconstructors.jl` |
| NativeMinuit / Optimization | Same statistical objective, parameter controls and nuisance refits; explicit solver capabilities | `test/extensions/native_minuit.jl`, `model_composition.jl` |

Analytic references cover normalization, gradients, Hessians, mixture
boundaries and nuisance refits. A named extended model agrees between Optim
and NativeMinuit; the [throughput probe](performance.md#Ecosystem-Throughput-Probe)
compares both paths at matched accuracy. Construction and normalization happen
once per objective evaluation; heterogeneous event mixtures use bounded batches.

### Real-Data Reference

The [LHCb three-hadron B-decay data](https://opendata.cern.ch/record/4900)
supply the [mass-spectrum example](gallery/lhcb_mass_spectrum.md): 3,420,295
MagnetUp candidates, 9,717 after the notebook's PID cuts, 7,368 in the fit window.
`examples/data/lhcb_mass/prepare.jl` verifies the ROOT checksum and rebuilds every
bin. The CC0 collision data, DOI, cuts and model limits are documented on the page;
the full ROOT download stays outside Git and the ordinary docs build.

`benchmarks/lhcb_reference.py` checks the two peak shapes independently using
SciPy CDFs and C++ Minuit2. Executed documentation checks the minima, local
errors and signal-yield profile against those references. This is a conditional
mass-spectrum fit, not a reproduction of an LHCb paper or a CP-asymmetry result.
Keep data preparation, warm fitting and uncertainty analysis separate in timing;
3.4 million scanned candidates must not be advertised as 3.4 million fitted events.

The binned distribution adapter processes each mixture component across all
bins, rather than dynamically dispatching and allocating component vectors
inside every bin. It retains log-space probability sums and derivatives through
truncation and zero weights. Regression checks compare values, gradients and
Hessians with the scalar formulation and bound allocations for 10,000 bins.
This removes the dominant adapter overhead in the LHCb example; upstream
CDF evaluation still costs more than its model-specific analytic control.

Unbinned heterogeneous mixtures reduce native component batches in blocks of
at most 4,096 events. This bounds simultaneous event-sized gradient/Hessian
scratch arrays without binning, sampling, or rebuilding the model per block.
The observations themselves remain resident; total work and cumulative
allocations still scale with the event count. Vector and joint-event regressions
check the complete reduction, including the final partial block and zero weights.

Likelihood quadrature uses the maximum norm of the value and all nested AD
coefficients for its estimated error, so a locally constant density cannot hide
an oscillatory derivative. Moving finite integration bounds are mapped to
`[0, 1]` without discarding their derivatives. CDFs that violate range/order by
roundoff are reintegrated from the PDF in `integration=:auto`; strict `:cdf`
and larger violations still fail. No probabilities are clipped. Analytic
regressions cover values, gradients and Hessians; a DistributionsHEP tail case
reproduces the original CDF-roundoff failure. Upstream normalization routines
remain responsible for the derivatives of their own probability objects.

The example uses upper-only truncation of an exponential's natural support.
Distributions 0.25.131's redundant lower bound at zero produces undefined AD
derivatives in its log normalizer; finite differences give the same optimum,
but the equivalent upper-only constructor retains valid automatic derivatives.
The shared scalar-solver boundary rejects non-finite initial gradients and,
when required, Hessians before entering the optimizer. The error identifies
the derivative failure and the explicit finite-difference option; it never
silently changes differentiation policy. This preflight does not prove derivative
correctness throughout the parameter domain. No upstream types are patched.

### Shared Contract

1. **Statistical meaning stays in ScientificFitting.** Solvers receive the same
   validated objective/residuals, parameter mapping, and constraints. Likelihood
   costs use ``-2\log L``; the Minuit adapter must use the matching error scale
   (`errordef=1`). Probability densities, discrete masses, and event intensities
   remain distinct. Dependent observations require a joint likelihood, not a
   product of marginal probabilities; discrete data do not imply discrete fit
   parameters are supported.
2. **Capabilities are explicit.** Check bounds, nonlinear constraints, required
   derivatives, and residual versus scalar objectives before solving. Reject
   unsupported requests instead of dropping constraints or adding hidden
   penalties. Retain explicit solver/derivative options in multistart and all
   nuisance-parameter refits. Document backend-specific tolerance and budget
   meanings rather than pretending iterations and function calls are identical.
3. **Results preserve evidence.** Normalize parameters, objective values, and
   convergence status while retaining native failure details. Keep minimizer
   selection independent of local covariance, MINOS intervals, and contour
   methods. Record missing crossings, parameter boundaries, and failed scans;
   do not substitute symmetric errors for an invalid asymmetric interval.
4. **Metadata is not statistics.** Constructor metadata used for optimizer step
   sizes is not a measurement uncertainty or Gaussian prior. Parameter terms
   remain explicit. Reject conflicting names/bounds and avoid silently mutating
   the user's constructor during optimization or profiling.
5. **Reuse expensive work at the correct scope.** Prepare a parameter-dependent
   distribution and its normalization once per distinct model at each parameter
   point, not once per observation. Preserve dual-number types, batched
   evaluation, stable log densities, and analytic bin integrals/CDF differences
   where available. Fall back to controlled quadrature or an explicit finite
   derivative mode where appropriate; never cache a normalization across changed
   parameters or claim smoothness for an arbitrary interpolated density.

### Solver And Posterior Checks

The solver contract is exercised with LsqFit, Optim through Optimization.jl,
bounded NLopt Nelder-Mead, and NativeMinuit. The tests cover residual,
scalar-gradient, nonlinear-constrained, derivative-free and profile-refit
paths. C++ Minuit2 provides an independent likelihood-fit reference.

Turing is a complementary posterior-inference workflow, not a MIGRAD substitute.
Its [external-likelihood interface](https://turinglang.org/docs/usage/external-likelihoods/index.html)
reuses `-problem.objective(p)/2` for an explicitly normalized data likelihood.
The [executed example](interfaces.md#Posterior-Inference) specifies priors and
support in Turing, without copying SF parameter terms or constraints. A
Gamma-Poisson reference checks data/prior separation, values, gradients and
Hessians, including an SF fit with an auxiliary parameter term. Four NUTS chains
are checked against the analytic posterior and Turing's chain diagnostics.
There is no Turing dependency or sampler wrapper in the core; posterior credible
intervals and profile confidence intervals retain their different meanings.

Each required adapter needs analytic or independent parameter/objective/error
references, a constraint/failure case, profile-refit parity, and executable
documentation. A nested signal-plus-background model must work through
BuildConstructors with two solvers without rewriting its likelihood. Measure
adapter overhead against direct upstream calls with the same objective and
accuracy, separating startup, warm fitting, and uncertainty analysis. Extend
existing targeted tests rather than duplicating a benchmark framework.

The core supports Julia 1.10. The optional
[NativeMinuit 0.7.2](https://github.com/fkguo/NativeMinuit.jl/blob/main/Project.toml)
requires Julia 1.11 and declares LGPL-2.1-or-later. The Turing 0.48 reference
runs on Julia 1.12. Plot tests cover CairoMakie 0.13 and 0.15.14+; the newer
line permits coexistence with Turing's current chain-plotting dependencies.
