# Backend Design

This page is the contributor map for ScientificFitting's numerical core.
User-facing signatures and defaults are in the [API Reference](api.md); the
statistical derivations are in [Mathematics and Statistics](statistics.md).

The dependency direction is one-way:

```math
\begin{aligned}
\text{public constructor}
&\longrightarrow \text{validated problem} \\
&\longrightarrow \text{objective and compatible solver} \\
&\longrightarrow \text{result} \\
&\longrightarrow \text{reports, diagnostics, profiles, and plots}.
\end{aligned}
```

Later stages may inspect an earlier result, never reconstruct or mutate the
statistical problem behind it.

## Two Problem Families

A residual vector and a general likelihood do not expose the same information,
so there are two normalized problem types.

| Normalized problem | Scientific payload and public path |
|---|---|
| `FitProblem` -> `FitResult` | Gaussian x-y data built by `fit_model` or `FitProblem`; stores the model, observations, uncertainty, and parameter controls; minimizes static ``\chi^2`` or normalized Gaussian ``-2\log L``. |
| `LikelihoodFitProblem` -> `LikelihoodFitResult` | Counts, samples, indexed data, or custom objectives built by the likelihood helpers or `fit_custom`; stores the objective, optional goodness statistic, observation count, and parameter controls. |

Both families share parameter bounds, fixed parameters, Gaussian parameter
terms, nonlinear constraints, multistart selection, local covariance,
diagnostics, and profile refits. They stay separate where their data contracts
differ: a generic likelihood need not have x-y residuals, predictions, or a
fit curve.

## One Fit, Step By Step

### 1. Normalize and validate the scientific input

Problem construction and the public fit entry points copy numeric inputs into
stable storage and reject mismatched dimensions, non-finite observations or
starting values, non-positive standard deviations, invalid covariance
matrices, inconsistent parameter indices, and fixed values outside declared
bounds before solver dispatch.

### 2. Map full parameters to optimizer coordinates

The scientific model always sees the complete parameter vector; fixed
parameters are removed only from the optimizer-visible vector:

```math
q_{\mathrm{free}}
\xrightarrow{\text{expand}}
p_{\mathrm{full}}.
```

Bounds are reduced to the same free coordinates, and nonlinear constraint
callbacks still receive `p_full`. After fitting, free covariance and Jacobian
blocks are embedded back into full parameter order. Ordinary fits, multistart
candidates, profiles, and contours reuse this one mapping.

### 3. Prepare reusable evaluation state

Static Gaussian uncertainty is prepared outside repeated scalar objective
evaluations:

- diagonal errors become inverse standard deviations and a log determinant,
- dense covariance becomes a Cholesky factor and a log determinant,
- sparse covariance keeps its sparse Cholesky factor,
- `WhiteningOperator` keeps the supplied matrix-free operation and determinant,
- correlated parameter constraints are factorized once.

This state lives in `FitEvaluationCache` (Gaussian) and
`LikelihoodEvaluationCache` (likelihood); the LsqFit path prepares equivalent
static weights for its native residual interface. Parameter-dependent
covariance — [effective x errors](statistics.md#Uncertainty-In-X) and
model-relative components — is recomputed at each parameter point, not cached.

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

The determinant term then changes the optimum and cannot be dropped; see
[Full Gaussian Likelihood](statistics.md#Full-Gaussian-Likelihood).
Likelihood problems provide their data objective directly on the
[``-2\log L`` scale](statistics.md#The-Cost-Convention). The shared parameter
layer adds Gaussian priors and correlated parameter constraints; bounds and
fixed parameters restrict the parameter space, never as hidden penalty terms.

### 5. Dispatch only to a compatible solver

Derivative selection is stored in the problem. `derivatives=:auto` keeps the
Julia defaults; `derivatives=:finite` uses central finite differences for
foreign callbacks that accept ordinary floats but not ForwardDiff duals, and
applies to optimization, constraint derivatives, post-fit covariance,
predictions, and profile/contour refits. Explicit model Jacobians and
x-derivatives take precedence. Tolerance and budget defaults — the finite-mode
`1e-6` versus `1e-10` split, NativeMinuit's EDM `0.1`, and the LsqFit
`maxiters`/`maxIter` and `tol`/`x_tol`+`g_tol` mapping — are tabulated under
[Solver Control](api.md#Solver-Control); explicit `tol` values are preserved,
and the result keeps the solver's actual convergence flag, so an iteration
limit never implies success.

Finite mode differentiates the **whole** objective, including any
parameter-dependent covariance and its log determinant; differentiated models
must be evaluable in a neighborhood of each evaluation point, including near
declared bounds. Pointwise x-error propagation uses four vectorized model
calls for a central fourth-order estimate of every `df/dx` with step
``h_i = \epsilon^{1/5}\max(|x_i|,1)``; analytic `x_derivative` callbacks
bypass it.

The Python bridge converts inputs once and enters the fit through
`invokelatest`, an inference boundary outside the numerical loop; foreign
callbacks use a private `_TypedCallback{R}` (concrete result type, one dynamic
dispatch), so different Python models share Julia solver specializations,
while native Julia models bypass it. Models and x-derivatives return vectors,
allocating parameter Jacobians return matrices. A `PrecompileTools` workload
caches Gaussian/Poisson kernels, prediction, and reporting without starting
Python or loading plotting.

Density helpers accept `vectorized=true` without changing their scalar
default: one density call per objective, QuadGK's `BatchIntegrand` with
per-bin adaptive error control, dual-number-preserving batch buffers reused
across bins, and the same density closure in profile refits. An analytic bin
integral supplied through `fit_histogram_model` avoids quadrature entirely.

Solver selection follows the represented problem rather than a speed
preference:

| Condition | Backend |
|---|---|
| All parameters fixed | no optimizer; evaluate the complete result once |
| Unbounded, static Gaussian chi-square without extra parameter terms | LsqFit least-squares path |
| Bounds, priors, parameter constraints, parameter-dependent covariance, or likelihood objective | Optimization.jl with LBFGS |
| Nonlinear equality or inequality constraints | Optimization.jl with IPNewton |
| Likelihood with explicit `solver=:nelder_mead` and no nonlinear constraints | OptimizationNLopt with native bounded Nelder-Mead, without derivatives |

There is no backend forcing: features outside static least squares route to
the scalar solver automatically. Backend selection may change how the
objective is minimized, never which objective is minimized.

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
derivative policy, multistart ranking, statistical summaries, and profile
refits. Adapters do not add priors, reinterpret error scales, or silently drop
constraints; unknown iteration counts stay `missing`. Solver settings persist
when a profile changes the number of free parameters, so NLopt adapters take
algorithm enums rather than dimension-bound native `Opt` instances.

`OptimizationSolver` uses SciML's algorithm traits and termination codes.
`ScientificFittingNativeMinuitExt` loads only with NativeMinuit and supplies
the same cost and derivative policy to MIGRAD with `errordef=1`; its native
result keeps covariance and failure evidence behind a non-specializing result
field, so one result type serves all models and solvers while the numerical
kernels still specialize. Conventions:
[Solver Adapters](api_fitting.md#Solver-Adapters).

Likelihood `optimizer` and `parameter_covariance` are independent options
stored in `FitOptions` and retained by profile refits; their contracts,
including Nelder-Mead's evaluation budget and Hessian-free default, are in
[Minimization And Local Errors](api_fitting.md#Minimization-And-Local-Errors).

One multistart ranking rule serves both families: a converged finite result
outranks an unconverged one; within the same status, the lower cost wins. If
every run stops early, the best finite result is returned with
`converged=false`.

ForwardDiff tags must not be shared between precompiled fits and unseen models
with nested x derivatives
([ForwardDiff #714](https://github.com/JuliaDiff/ForwardDiff.jl/issues/714));
the startup gate checks exactly this in a fresh process. CHOLMOD's sparse
solves reject ForwardDiff duals, so static sparse covariance requires
`derivatives=:finite` with the general optimizer or an AD-compatible
`WhiteningOperator`; the least-squares path needs no override, Python uses
finite derivatives consistently (including sparse fits with bounds and profile
refits), and sparse covariance components are validated through their stored
entries, without a dense copy.

### 6. Build the result once

`FitResult` and `LikelihoodFitResult` are the numerical source of truth:
selected minimum, solver status, parameter estimates, local covariance,
correlations, statistics, and diagnostics. Gaussian x-y results additionally
retain predictions, raw and weighted residuals, and the weighted Jacobian.

For static least squares, local covariance comes from the weighted Jacobian;
otherwise from the objective Hessian on the ``-2\log L`` scale, evaluated once
and reused for diagnostics. With `parameter_covariance=:none`, free errors and
covariances are `NaN`, fixed errors remain zero, and diagnostics explain the
omission rather than claiming a curvature failure. Conditioning uses free
coordinates; a fixed parameter's zero variance is not a singularity. The
covariance is a
[local approximation](statistics.md#Local-Parameter-Covariance); profiles and
contours remain separate refit operations.

### 7. Read the result without changing it

`report_text`, `diagnose`, and `diagnostic_dashboard` consume result fields.
`profile` and `contour` reuse the stored normalized problem, fix one or two
parameters, and refit the remaining nuisance parameters. Plotting is an
optional CairoMakie package extension consuming the same result objects. None
of these paths reruns or alters the original fit unless the API explicitly
describes a refit.

Automatic scan ranges intersect the local-covariance range with declared
bounds; explicit grids are preserved, and missing threshold crossings remain
missing rather than being replaced by a bound, in Gaussian and
general-likelihood profiles, contours, and matrices alike.

The [Python renderer](python.md) uses native Matplotlib objects, not Makie.
Both renderers share residual/ratio preparation; observation pulls use stored
whitened residuals without auxiliary parameter terms, completed
profile/contour snapshots render without further scans, and Matplotlib's
constrained layout reserves outside legends and reports without a separate
layout engine or resize callbacks.

## Numerical Invariants

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

The core is split by responsibility, not by feature-specific vertical stacks.

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

This map is also a review rule: a plotting feature should not add a second
statistical calculation, and a new optimizer should not own covariance
semantics.

## Where A New Feature Belongs

Before adding a type or abstraction, ask whether an existing problem can
already express the required statistics.

| Change | Preferred integration |
|---|---|
| New convenience fitting function | validate its domain-specific inputs, then construct an existing problem type |
| New static uncertainty representation | add validation, preparation/whitening, determinant semantics, and an analytic covariance reference |
| New likelihood family | provide a ``-2\log L`` objective, a justified goodness statistic when one exists, and an explicit observation count |
| New numerical backend | add a compatibility predicate and prove that the represented objective is unchanged |
| New diagnostic | consume a result or profile object and return structured evidence plus an action |
| New report or plot | consume existing result fields; keep rendering inside the optional extension |

Do not add a parallel result type, cache, or solver path merely to support a
new presentation.

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

The core gate is `julia --project=. test/core_runtests.jl`; the complete
package gate is `julia --project=. test/runtests.jl`. Benchmark methodology is
in [Performance Checks](backend_design.md#Performance-Checks).

For a page-level output check, run
`julia --project=docs test/docs_output_snapshots.jl gallery/resonance_decay.md`
(more page paths are accepted; without arguments it executes every documented
workflow). It compares the displayed output with both the page's code cells
and the example generator; only solver iteration counts are normalized across
solver/platform versions.

## [Ecosystem Integrations](@id v03-ecosystem)

[Packages and Interfaces](interfaces.md) holds executable examples of
distribution fitting, named model construction, interchangeable solvers and
posterior inference. The [Python interface](python.md) keeps NumPy callbacks
and Matplotlib plots and does not wrap Julia constructor or solver objects.

| Interface | Contract | Focused reference tests |
|:---|:---|:---|
| Distribution objects | Reuse upstream `logpdf`/`pdf`/`cdf`; continuous, discrete and joint observations, binned and extended likelihoods | `test/extensions/distribution_ecosystem.jl`, `composite_distributions.jl` |
| BuildConstructors | Preserve names, starts, bounds, fixed/shared parameters and reconstruction; independent of Minuit | `test/extensions/buildconstructors.jl` |
| NativeMinuit / Optimization | Same statistical objective, parameter controls and nuisance refits; explicit solver capabilities | `test/extensions/native_minuit.jl`, `model_composition.jl` |

Analytic references cover normalization, gradients, Hessians, mixture
boundaries and nuisance refits. A named extended model agrees between Optim
and NativeMinuit; the opt-in probe `benchmarks/ecosystem.jl` compares both
paths at matched accuracy.

### Real-Data Reference

The [LHCb three-hadron B-decay data](https://opendata.cern.ch/record/4900)
supply the [mass-spectrum example](gallery/lhcb_mass_spectrum.md): 3,420,295
MagnetUp candidates, 9,717 after the notebook's PID cuts, 7,368 in the fit
window. `examples/data/lhcb_mass/prepare.jl` verifies the ROOT checksum and
rebuilds every bin; the CC0 collision data, DOI, cuts and model limits are
documented on the page, and the full ROOT download stays outside Git and the
ordinary docs build.

`benchmarks/lhcb_reference.py` checks the two peak shapes independently using
SciPy CDFs and C++ Minuit2; executed documentation checks the minima, local
errors and signal-yield profile against those references. This is a
conditional mass-spectrum fit, not a reproduction of an LHCb paper; keep data
preparation, warm fitting and uncertainty analysis separate in timing, and do
not advertise 3.4 million scanned candidates as 3.4 million fitted events.

The binned distribution adapter processes each mixture component across all
bins rather than dispatching per bin, retaining log-space probability sums and
derivatives through truncation and zero weights; regressions compare values,
gradients and Hessians with the scalar formulation and bound allocations for
10,000 bins. Unbinned heterogeneous mixtures reduce native component batches
in blocks of at most 4,096 events, bounding event-sized gradient/Hessian
scratch without binning, sampling, or rebuilding the model per block; total
work still scales with the event count, and regressions cover the final
partial block and zero weights.

Likelihood quadrature estimates its error from the maximum norm of the value
and all nested AD coefficients, so a locally constant density cannot hide an
oscillatory derivative. Moving finite integration bounds are mapped to
`[0, 1]` without discarding derivatives. CDFs violating range/order by
roundoff are reintegrated from the PDF in `integration=:auto`; strict `:cdf`
and larger violations still fail, and no probabilities are clipped. A
DistributionsHEP tail case reproduces the original CDF-roundoff failure.

The example uses upper-only truncation of an exponential's natural support:
Distributions 0.25.131's redundant lower bound at zero produces undefined AD
derivatives in its log normalizer, while the equivalent upper-only constructor
keeps valid automatic derivatives. The shared scalar-solver boundary rejects
non-finite initial gradients and, when required, Hessians; the error names the
derivative failure and the finite-difference option and never silently changes
differentiation policy. No upstream types are patched.

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
bounded NLopt Nelder-Mead, and NativeMinuit, across residual, scalar-gradient,
nonlinear-constrained, derivative-free and profile-refit paths; C++ Minuit2
provides an independent likelihood-fit reference.

Turing is a complementary posterior-inference workflow, not a MIGRAD
substitute. Its
[external-likelihood interface](https://turinglang.org/docs/usage/external-likelihoods/index.html)
reuses `-problem.objective(p)/2` for an explicitly normalized data likelihood;
the [executed example](interfaces.md#Posterior-Inference) specifies priors and
support in Turing, without copying SF parameter terms or constraints. A
Gamma-Poisson reference checks data/prior separation, values, gradients,
Hessians, and four NUTS chains against the analytic posterior and Turing's
chain diagnostics. The core has no Turing dependency or sampler wrapper;
posterior credible intervals and profile confidence intervals keep their
different meanings.

Each required adapter needs analytic or independent parameter/objective/error
references, a constraint/failure case, profile-refit parity, and executable
documentation. A nested signal-plus-background model must work through
BuildConstructors with two solvers without rewriting its likelihood. Measure
adapter overhead against direct upstream calls at the same objective and
accuracy, separating startup, warm fitting, and uncertainty analysis, by
extending existing targeted tests rather than duplicating a benchmark
framework.

The core supports Julia 1.10. The optional
[NativeMinuit 0.7.2](https://github.com/fkguo/NativeMinuit.jl/blob/main/Project.toml)
requires Julia 1.11 and declares LGPL-2.1-or-later; the Turing 0.48 reference
runs on Julia 1.12. Plot tests cover CairoMakie 0.13 and 0.15.14+, the newer
line permitting coexistence with Turing's chain-plotting dependencies.

## Performance Checks

The canonical benchmark run is

```bash
julia --project=benchmarks benchmarks/runbenchmarks.jl --seconds=1
```

with `--save` and `--compare` for TOML baselines; a comparison fails on missing
benchmark cases or mismatched machine metadata. Two CI gates back the claims on
this page: `test/startup_probe_gate.jl` verifies fresh-process loading —
`using ScientificFitting` without Makie — together with nested automatic x
derivatives, and `test/performance_budget_gate.jl` checks broad steady-state
budgets for representative warmed hot paths, catching regressions such as
losing the fast path or recomputing static covariance work inside the
objective. The LsqFit least-squares fast path is selected automatically for
unbounded static chi-square fits (`[-Inf, Inf]` counts as no-op bounds), and
its weighted Jacobian is reused during `FitResult` construction.
