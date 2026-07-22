# Backend Design

This page is for contributors and advanced users who need to audit how a fit is
constructed. It records the numerical structure that must stay clear as
JuFitter grows; public signatures and defaults remain in the
[API Reference](api.md).

```@raw html
<div class="jufitter-flow">
  <div class="jufitter-flow-step"><strong>fit_model</strong><span>Normalizes user data and options.</span></div>
  <div class="jufitter-flow-step"><strong>FitProblem</strong><span>Stores model, data, bounds, fixed parameters, priors, constraints.</span></div>
  <div class="jufitter-flow-step"><strong>Evaluation Cache</strong><span>Prepares static covariance factors and reusable weights.</span></div>
  <div class="jufitter-flow-step"><strong>Backend Choice</strong><span>Fast LsqFit path or general Optimization.jl path.</span></div>
  <div class="jufitter-flow-step"><strong>Result</strong><span>Parameters, statistics, covariance, diagnostics.</span></div>
  <div class="jufitter-flow-step"><strong>Output</strong><span>Reports and plots never change the numerical result.</span></div>
</div>
```

## Separation Of Responsibilities

The optimizer is not allowed to define the statistics. JuFitter separates the
fit into:

```math
\text{FitProblem} + \text{CostFunction} + \text{Optimizer}
\longrightarrow
\text{FitResult}.
```

The `FitProblem` stores data, model, uncertainty inputs, fixed parameters,
bounds, priors, and constraints. The cost layer decides what is minimized. The
optimizer only solves the numerical minimization problem.

The low-level `fit(problem)` methods extend `StatsAPI.fit`; JuFitter does not
create a competing generic with the same ecosystem-wide name.

Plotting is a package extension, not a numerical backend dependency.
`plotting_api.jl` defines the public plotting functions and informative
fallbacks, while `ext/JuFitterCairoMakieExt.jl` loads the CairoMakie
implementation from `plotting.jl`. Consequently `using JuFitter` keeps fitting,
reports, diagnostics, profiles, and profile-matrix computation Makie-free;
`using CairoMakie` activates rendering methods without changing a fit result.

## Cost Functions

For Gaussian XY fits, JuFitter currently distinguishes:

- `:chi2`: quadratic residual cost for static Gaussian uncertainties.
- `:gaussian_likelihood`: full Gaussian likelihood cost on the ``-2\log L``
  scale, including normalization and log-determinant terms.

`cost=:auto` chooses `:gaussian_likelihood` when the effective covariance
depends on the parameters, especially for x uncertainties or model-relative y
uncertainties. Otherwise it uses `:chi2`.

For residuals

```math
r(\theta)=y-m(x,\theta)
```

and covariance ``V``, the chi-square cost is

```math
\chi^2(\theta)=r(\theta)^T V^{-1}r(\theta).
```

For the full Gaussian likelihood cost:

```math
-2\log L(\theta)
=
n\log(2\pi)
+
\log\det V(\theta)
+
r(\theta)^T V(\theta)^{-1}r(\theta).
```

The log determinant is essential when ``V`` depends on ``\theta``.

## Backend Selection

The fast path uses `LsqFit` for unbounded static chi-square fits. Analytic
Jacobians are forwarded to the backend when provided, and the weighted Jacobian
computed by `LsqFit` is reused when constructing the `FitResult`.

The general path uses `Optimization.jl` for scalar objectives, bounds,
constraints, priors, parameter-dependent covariance, and likelihood workflows.

An explicit backend request cannot weaken the statistical problem.
`backend=:lsqfit` is rejected when the fit contains a non-chi-square cost,
bounds, priors, parameter constraints, nonlinear constraints, active error
components, or parameter-dependent covariance. Use `backend=:auto` unless a
specific compatible backend is needed for a controlled comparison.

No-op bounds such as `[-Inf, Inf]` are normalized and do not block the fast
least-squares path. This matters because generic APIs often pass bounds even
when they do not mathematically constrain the problem.

## Covariance and Whitening

Static uncertainty information is prepared in `FitEvaluationCache`.

Diagonal covariance stores inverse standard deviations and log determinants.
Dense static covariance stores a Cholesky factor and log determinant. Residuals
are whitened by linear solves, not by multiplying with an explicit inverse
matrix.

`WhiteningOperator` is the matrix-free static path. Its user function applies a
complete operator ``W`` with ``W^T W=V^{-1}``; JuFitter applies the same
operation to residuals and Jacobian columns. The supplied log determinant keeps
the normalized Gaussian ``-2\log L`` value and information criteria consistent.
The operator is exclusive with other observation-uncertainty inputs, so the
cache never silently combines covariance models with unknown semantics.

Parameter-dependent covariance remains dynamic by design. Recomputing the
effective covariance for x uncertainties or model-relative uncertainty is part
of the statistical model, not accidental overhead.

## Parameter Mapping

Fixed parameters are removed from the optimizer-visible vector. JuFitter maps
between:

- the full parameter vector used by the model,
- the free parameter vector used by the optimizer,
- reported parameter estimates and uncertainties.

This keeps fixed parameters, priors, bounds, covariance dimensions, and degrees
of freedom consistent. A fixed parameter is still part of the scientific model,
so it must satisfy any declared bound for that parameter. Profile and contour
refits use the same rule; a scan point outside a bound is a failed refit, not a
valid uncertainty point.

## Diagnostics

`FitDiagnostics` records warnings, condition numbers, active bounds, and
structured diagnostic findings. The diagnostics layer must remain separate from
the numerical result: it interprets the result, but it does not silently repair
or modify it.

Important cases that must stay visible:

- optimizer non-convergence,
- non-positive degrees of freedom,
- unavailable goodness-of-fit statistics,
- ill-conditioned covariance or Hessian matrices,
- active bounds,
- strong parameter correlations,
- suspicious residual structure,
- failed profile or contour refits,
- non-finite or non-positive-semidefinite local parameter covariance.

## Scaling Limits

Diagonal and uncorrelated problems are the target class for very large datasets.
Dense covariance matrices are supported and tested, but they scale as `O(n^2)`
memory and `O(n^3)` factorization time.

Large correlated datasets can provide a `WhiteningOperator` instead of a dense
matrix. Its cost and storage are defined by the application-specific operation;
an AR(1) recurrence is ``O(n)`` in both time and working memory. Built-in
banded, Toeplitz, low-rank-plus-diagonal, and sparse-precision wrappers remain
future extensions. This is a data-structure problem, not a micro-optimization
of dense matrices.

Large ordinary least-squares problems can use `fit_model(...; inplace=true)`
with `model!(out, x, p)` and, optionally, `jacobian!(J, x, p)`. The LsqFit path
then evaluates model values and residuals in reusable solver buffers. General
constrained and parameter-dependent-covariance objectives retain the same user
contract, but automatic differentiation still requires temporary dual-number
buffers.

## Extension Points

The main backend extension points are:

- built-in structured covariance wrappers on top of the whitening contract,
- parameter-dependent and structured x-covariance operators,
- in-place likelihood and custom-objective evaluation,
- analytic Jacobian hooks for likelihood workflows,
- optimizer fallback and parameter scaling policies,
- stronger profile/contour refinement,
- ODE/PDE model adapters once the base fitting API is stable.

Every extension should add reference tests before changing statistical
semantics and benchmark coverage before changing a hot path.
