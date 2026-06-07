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
- Correlated Gaussian parameter constraints are prepared once as
  `PreparedParameterConstraint` objects. The same Cholesky factor and log
  determinant are reused for weighted residuals, chi-square, Gaussian NLL,
  likelihood objectives, likelihood Hessians, and likelihood covariance
  estimates.
- Unbounded static chi-square fits use `LsqFit`.
- Generic scalar objectives, bounds, constraints, priors, likelihood fits, and
  parameter-dependent covariance use `Optimization.jl`.

Avoid adding work inside objective functions unless it is mathematically
necessary. In particular, do not rebuild static covariance factors or allocate
large temporary matrices inside every objective call. If a quantity depends
only on the data, fixed uncertainty model, or fixed parameter constraints, it
belongs in a cache.

## Known Bottlenecks

- Dense covariance matrices scale as `O(n^2)` memory and `O(n^3)`
  factorization. Large correlated datasets need structured covariance
  operators, not denser micro-optimizations.
- Parameter-dependent x uncertainties require recomputing effective covariance
  terms. This is statistical work, not accidental overhead. For large datasets,
  prefer the public `x_derivative=(x, p) -> dy_dx` hook over the default
  pointwise AD derivative.
- Generic `Optimization.jl` objectives are flexible but slower than the
  specialized `LsqFit` path.
- Plotting startup is dominated by Makie/CairoMakie compilation on first use,
  but plotting is isolated in an optional CairoMakie package extension so the
  fitting/reporting core does not load Makie.
- Full gallery generation will become expensive once docs embed many high
  quality plots; split fast docs checks from full visual builds if needed.

## Plot Layout Rules

- When `show_legend=true` and `stats_position=:right`, `plot_fit` places the
  legend above the right-side report instead of consuming data-axis space.
- Plot styles are limited to three explicit output contracts: `:workbench`,
  `:showcase`, and `:publication`. Light/dark appearance and LaTeX conversion
  are independent options rather than additional style names.
- `:workbench` is the direct analysis/notebook contract: sans typography,
  restrained blue fit color, short error-bar whiskers, plain-text reporting
  unless the user explicitly passes LaTeX labels.
- `:showcase` is allowed to use a stronger but still scientific color
  hierarchy for documentation and talks. It must remain readable and must not
  become a decorative product-design theme.
- The right-side panel is a shared layout component, not a manually positioned
  axis. Legend entries, model text, parameters, and statistics are left-aligned
  and use one compact vertical rhythm. Custom Makie figures should use
  `plot_info_panel!`, `plot_theme`, and `plot_palette` instead of rebuilding
  this hierarchy.
- `plot_fit` preserves its declared `figure_size`. Layout compaction may remove
  unused grid tracks, but it must never resize the output around a short legend
  or report. This keeps data axes comparable across fits and prevents later
  Makie additions from changing the exported footprint.
- `plot_info_panel!` is top-aligned and does not report its compact content
  height to the parent layout. The data axis or compound diagnostic grid owns
  the row height; the report must never collapse or vertically center it.
- Standard fit, profile, contour, residual, and diagnostic plots all resolve
  color and light/dark appearance through the same style contract. A
  documentation-wide style selector must only reference plots that have
  rendered variants for every supported style and appearance.
- `plot_profile_matrix` is a composition layer over existing profile and
  contour refits. It must not introduce independent confidence-threshold
  semantics; threshold interpretation belongs to `profile`, `contour`, and the
  associated diagnostics.
- Post-fit annotation helpers (`fit_axis`, `add_curve!`, `add_points!`,
  `add_vline!`, `add_hline!`, `add_vband!`, `add_hband!`) are thin Makie
  wrappers. They should stay small, return Makie plot objects, validate obvious
  shape/range errors, and never trigger a refit or mutate `FitResult`.
- Default fit plots use larger axis, tick, legend, and report typography than
  the compact publication theme. Documentation assets use an additional
  web-specific scale because raster plots are reduced inside the Documenter
  content column.
- Documentation light/dark assets must share the
  `jufitter-plot-light`/`jufitter-plot-dark` classes. CSS guarantees exactly one
  visible variant for the active Documenter theme.
- Documentation plot-style switching is asset-based. A plot may declare
  `data-jufitter-plot-group` and `data-jufitter-plot-style`; the site selector
  then chooses the selected Makie-rendered style and falls back to `:showcase`
  when a page has not yet rendered all variants. Do not recolor or invert plots
  in CSS.
- Custom documentation plots must derive theme, palette, marker sizes,
  whisker widths, line widths, and right-side report typography from
  `plot_theme` and `plot_palette`. Only scientifically semantic elements, such
  as a threshold marker or a second model regime, should override individual
  Makie colors explicitly.
- Compound gallery figures with residual or pull panels are part of the same
  plotting contract. Poisson, histogram, damped-oscillator, constraints/profile,
  contour, and multi-dataset assets must be regenerated as
  `:workbench`, `:showcase`, and `:publication` variants for both light and
  dark appearances whenever their generator changes.
- Wide documentation plots expand from the article's left edge toward the
  available right side. They must not be viewport-centered because that can
  place the image underneath Documenter's fixed navigation sidebar.
- Gallery asset regeneration is a visual change and must be checked in both
  themes before commit.
- Contour diagnostics default to labeled confidence regions and an optional
  local-covariance overlay. Delta-cost heatmaps are opt-in because they are
  useful for surface inspection but slower to interpret as uncertainty regions.
- A gallery plot must state the uncertainty model, band type, and sigma level.
  Controlled demonstrations must be labeled as such; constructed data must
  never be presented as archival measurements.
- Documentation code cells that print a report, diagnostic dashboard, or
  derived scientific quantity should be followed by a compact
  `jufitter-cell-output` block. The block is not a second source of truth; it is
  a reader-facing checkpoint showing what the code is expected to produce. The
  output text must come from a real run of the displayed code, apart from
  explicit abridgement of repeated or verbose lines.

## Numerical Rules

- Do not use explicit matrix inverse in production covariance calculations.
- Prefer factorization and linear solves.
- Do not silently repair invalid covariance matrices. If a covariance matrix is
  not finite, symmetric, and positive definite, fail early or make the repair an
  explicit documented user policy with a visible diagnostic.
- Reject non-finite data, parameters, uncertainties, and model output before
  optimizer internals can turn them into vague convergence failures.
- Reject user-provided start values that are outside declared bounds. `p0` and
  explicit `initial_guesses` are part of the scientific input, so JuFitter must
  not silently clip them to the nearest allowed value. Generated multistart
  candidates may be constructed inside bounds, but user input should either be
  accepted as written or fail clearly.
- Reject invalid parameter-space metadata during problem normalization. Prior
  means and sigmas, fixed-parameter values and optional uncertainties, and
  correlated-constraint means/covariances must be finite. Correlated
  parameter-constraint covariance matrices must also be symmetric positive
  definite; do not rely on `Symmetric(cov)` to hide an asymmetric input matrix.
- Keep diagnostics visible when covariance, Hessian, ndf, p-values, or bounds
  make local errors unreliable.
- Fixed parameters are not allowed to bypass bounds. Profile and contour refits
  must apply the same invariant: scan points outside declared bounds are failed
  refits, not valid uncertainty samples.
- Profile and contour scans default to marking failed refits as non-finite grid
  points instead of aborting the whole scan. `diagnose(profile_result)` and
  `diagnose(contour_result)` must surface those failures before users interpret
  intervals or contour topology.
- `profile(...; adaptive=true)` refines intervals that bracket the selected
  profile threshold. `contour(...; adaptive=true)` refines grid cells whose
  corner values bracket requested contour levels. These refinements improve
  diagnostic resolution without making the whole scan uniformly dense.
- Contour levels are finite, positive delta-cost thresholds. JuFitter sorts and
  deduplicates them before scanning so adaptive refinement and filled-region
  plotting use one unambiguous ordering.
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
