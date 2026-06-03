# Documentation Plan

This page tracks the documentation work needed before JuFitter should be
advertised publicly. The target is a stable, modern, easy-to-host site with the
tone of Julia, Pluto, Makie, and SciML: computation made accessible, but without
weakening the mathematics.

## Publication Readiness

The current codebase is a solid v0 foundation: the core fit paths, statistical
semantics, plotting backend, examples, benchmarks, and tests are in place. It
is not ready for broad promotion until the documentation reaches the same
standard as the code.

The first public release should wait for:

- A complete quickstart path from installation to first plot.
- A polished gallery with generated figures.
- Function-level documentation for every public export.
- A theory section that explains which statistical model each helper uses.
- Troubleshooting pages for startup time, plotting backends, slow fits, and
  covariance pitfalls.
- A clean design layer with readable light and dark modes.

## Site Architecture

- Home: clear promise, minimal example, links to common entry points.
- Install: package install, project activation, first compile, plot backends,
  troubleshooting.
- Quickstart: one complete path with `fitplot`, `fit_model`, reports, and
  diagnostics.
- Gallery: example-driven documentation with generated figures.
- Concepts: statistics, backend design, performance, covariance semantics.
- API Reference: detailed docs for every exported type and function.
- Maintenance Notes: architecture, bottlenecks, known limitations, and future
  optimization points.

## Design Direction

The website should feel closer to a modern scientific product than a default
API dump. The practical target:

- Documenter.jl for stable Julia-native hosting and deployment.
- Evaluate DocumenterVitepress.jl for a more modern navigation and dark-mode
  experience.
- Literate examples for gallery pages, so code, text, and figures stay in sync.
- Light and dark plot themes generated from the same gallery scripts.
- Visual references: Makie/Beautiful Makie for plot galleries, SciML for
  large-scale technical navigation, Pluto for accessible language and immediacy.

## API Documentation Requirements

Every exported symbol must answer:

- What problem does this solve?
- What are the required arguments?
- What keyword arguments exist, with defaults and units/semantics where
  relevant?
- What statistical model is minimized?
- What does the result contain?
- What can go wrong, and what warning should the user expect?
- Minimal example and one realistic example where useful.

Priority public exports:

- `fitplot`, `plot_fit`, `plot_residuals`, `plot_diagnostics`,
  `plot_profile`, `plot_contour`.
- `fit`, `fit_model`, `fit_custom`, `fit_poisson_model`,
  `fit_histogram_model`, `fit_histogram_density`, `fit_unbinned_model`,
  `fit_extended_unbinned_model`, `fit_indexed_model`, `fit_multi_model`.
- `profile`, `profile_interval`, `contour`.
- `fit_report`, `report_text`.
- `FitProblem`, `FitResult`, `LikelihoodFitProblem`, `LikelihoodFitResult`,
  `FitStatistics`, `FitDiagnostics`, `ParameterPrior`, `FixedParameter`,
  `ParameterConstraint`, `ErrorComponent`, `ConstraintSpec`.

## Theory Requirements

The mathematics section must explain:

- Chi-square fits with diagonal and dense covariance.
- Full Gaussian negative log-likelihood and the log-determinant term.
- Effective variance for x uncertainties and when it is only an approximation.
- Poisson likelihood and Poisson deviance.
- Histogram, unbinned, and extended-unbinned likelihoods.
- Priors, parameter constraints, fixed parameters, bounds, and ndf semantics.
- Local covariance from Jacobian or Hessian.
- Profile likelihoods, contour levels, Wilks approximations, and limitations.
- AIC, BIC, p-values, and when not to trust them.

Each theory page should end with practical guidance: which JuFitter function to
use, which options matter, and which diagnostics to inspect.

## Gallery Build Plan

Examples should be executable and visually checked. The target workflow is:

- Write gallery scripts under `examples/gallery/`.
- Save outputs under `examples/output/`.
- Convert stable examples into docs pages via Literate or a small build helper.
- Use real data subsets where possible.
- Add a smoke test for each gallery script once it stabilizes.

Course-specific raw data can be used as source material, but public
documentation must only contain curated datasets with neutral names, minimal
columns, and provenance notes that explain the measurement type without local
paths or course-internal context.

## Current Documentation Debt

- `api.md` currently relies on raw autodocs. It needs curated pages organized by
  workflow.
- The gallery needs at least eight polished long-form examples with real or
  realistic datasets, visible uncertainty information, and light/dark plot
  exports.
- The math section needs an entry-level layer for engineers and beginners before
  the precise likelihood/covariance treatment.
