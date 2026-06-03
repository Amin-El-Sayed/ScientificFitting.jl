# Maintenance Notes

This page exists so future contributors can understand JuFitter's structure,
known bottlenecks, and extension points without reverse-engineering the whole
package.

## Core Structure

- `types.jl`: public data structures and argument normalization.
- `parameters.jl`: mapping between full parameters and optimizer-visible free
  parameters.
- `weights.jl`: residuals, covariance preparation, whitening, Jacobians,
  covariance estimates, and backend selection.
- `costs.jl`: chi-square and Gaussian NLL semantics.
- `fit.jl`: XY fit orchestration, backend calls, and `FitResult` construction.
- `likelihood_fits.jl`: custom, Poisson, histogram, unbinned, indexed, and
  multi-fit likelihood workflows.
- `profile.jl`: profile likelihoods, intervals, and two-parameter contours.
- `plotting.jl`: Makie themes, fit plots, residuals, diagnostics, profiles,
  and contours.
- `report.jl`: structured fit reports and plain-text reporting.

## Hot Paths

- Static diagonal and dense covariance information is prepared in
  `FitEvaluationCache`.
- Diagonal covariance uses precomputed inverse standard deviations.
- Dense covariance uses Cholesky solves and cached log determinants.
- Unbounded static chi-square fits use `LsqFit`.
- Generic scalar objectives, bounds, constraints, priors, likelihood fits, and
  parameter-dependent covariance use `Optimization.jl`.

Avoid adding work inside objective functions unless it is mathematically
necessary. In particular, do not rebuild static covariance factors or allocate
large temporary matrices inside every objective call.

## Known Bottlenecks

- Dense covariance matrices scale as `O(n^2)` memory and `O(n^3)`
  factorization. Large correlated datasets need structured covariance
  operators, not denser micro-optimizations.
- Parameter-dependent x uncertainties require recomputing effective covariance
  terms. This is statistical work, not accidental overhead.
- Generic `Optimization.jl` objectives are flexible but slower than the
  specialized `LsqFit` path.
- Plotting startup is dominated by Makie/CairoMakie compilation on first use.
- Full gallery generation will become expensive once docs embed many high
  quality plots; split fast docs checks from full visual builds if needed.

## Numerical Rules

- Do not use explicit matrix inverse in production covariance calculations.
- Prefer factorization and linear solves.
- Do not silently repair invalid covariance matrices. If a covariance matrix is
  not finite, symmetric, and positive definite, fail early or make the repair an
  explicit documented user policy with a visible diagnostic.
- Reject non-finite data, parameters, uncertainties, and model output before
  optimizer internals can turn them into vague convergence failures.
- Keep diagnostics visible when covariance, Hessian, ndf, p-values, or bounds
  make local errors unreliable.
- Profile and contour scans default to marking failed refits as non-finite grid
  points instead of aborting the whole scan. `diagnose(profile_result)` and
  `diagnose(contour_result)` must surface those failures before users interpret
  intervals or contour topology.
- `profile(...; adaptive=true)` refines intervals that bracket the selected
  profile threshold. `contour(...; adaptive=true)` refines grid cells whose
  corner values bracket requested contour levels. These refinements improve
  diagnostic resolution without making the whole scan uniformly dense.
- `diagnostic_dashboard(...)` is a summary layer over `diagnose(...)`. It must
  not introduce independent statistical rules; it only converts structured
  findings into status, counts, and prioritized next actions for lab use.
- Add analytic reference tests before changing statistical semantics.
- Add benchmark coverage before changing a hot path.

## Hardening Policy

Smoke tests are not release evidence. They only show that a friendly path still
executes. The release gate for the core is a torture suite with deliberately
hostile but realistic cases:

- invalid and non-finite inputs,
- negative or zero uncertainties,
- singular, nonsymmetric, or ill-conditioned covariance matrices,
- badly scaled but identifiable models,
- active bounds and constraints,
- local minima and difficult initial values,
- large datasets with analytic Jacobians,
- parameter-dependent covariance and x-uncertainty paths,
- likelihood, profile, contour, and multi-dataset workflows.

When a torture test fails, prefer a small architectural fix over a narrow test
workaround. The intended outcome is not that every hostile fit succeeds; it is
that every result or failure is statistically interpretable.

## Extension Points

- Structured covariance API: banded, Toeplitz, low-rank-plus-diagonal, sparse
  precision, and custom whitening operators.
- In-place model evaluation for very large datasets.
- Analytic Jacobian hooks for likelihood models, not only XY models.
- Specialized optimizers for common likelihood families.
- ODE/PDE model adapters once the base API is stable.
- Literate gallery build system with light and dark plot exports.

## Release Checklist

- `julia --project=. -e 'using Pkg; Pkg.test()'` passes.
- `julia --project=docs docs/make.jl` builds without new warnings.
- `git diff --check` is clean.
- Public exports have docstrings and curated reference entries.
- Gallery examples run from a clean checkout.
- Benchmarks cover any modified hot path.
- The roadmap and maintenance notes mention any known limitation that users
  could reasonably hit.
