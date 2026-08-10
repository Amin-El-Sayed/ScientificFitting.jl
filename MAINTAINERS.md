# Maintainer Contract

This file records invariants that should survive individual implementations.
Current release evidence belongs in `RELEASE_AUDIT.md`, future work in
`ROADMAP.md`, commands in `RELEASE_CHECKLIST.md`, and contributor workflow in
`DEVELOPMENT.md`. The public editorial standard lives in `DOCUMENTATION.md`.
Do not duplicate those documents here.

## Source Boundaries

- `types.jl` normalizes public problem and result structures.
- `parameters.jl` maps full parameter vectors to optimizer-visible free
  vectors.
- `weights.jl` prepares covariance, whitening, residual, and Jacobian data.
- `costs.jl` defines chi-square and Gaussian likelihood semantics.
- `fit.jl` selects a compatible backend and constructs `FitResult`.
- `likelihood_fits.jl` implements count, histogram, unbinned, indexed,
  multi-dataset, and custom-objective workflows.
- `profile.jl` owns profile, contour, interval, and profile-matrix computation.
- `report.jl` owns Makie-free structured and text output.
- `plotting_api.jl` defines the optional public plotting boundary.
- `plotting.jl` is loaded only by `ext/JuFitterCairoMakieExt.jl`.

`Backend Design` in the documentation explains the full data flow. Keep this
list as an ownership map, not a second architecture chapter.

## Numerical Invariants

1. The statistical model is defined before backend selection. A requested fast
   backend must be rejected if it would discard bounds, constraints, priors,
   covariance terms, or likelihood normalization.
2. Static covariance factors, log determinants, and correlated parameter
   constraints are prepared once and reused. Work that depends only on fixed
   input does not belong inside every objective evaluation.
3. Production covariance calculations use factorizations and linear solves,
   never explicit matrix inverses.
4. Invalid scientific input is not silently repaired. Non-finite values,
   non-positive uncertainties, invalid bounds, and non-SPD covariance fail
   before optimizer internals obscure the cause.
5. User starting values are accepted as supplied or rejected. They are not
   clipped into bounds. Internally generated multistart candidates may be
   constructed inside the valid domain.
6. Parameter-dependent covariance keeps automatic-differentiation information
   through the actual objective factorization. Stripping dual values is allowed
   only for finite-value or symmetry validation.
7. Parameter covariance is a local quadratic approximation. Diagnostics must
   direct users to profiles or contours when bounds, conditioning, correlation,
   or scan geometry make symmetric errors unreliable.
8. Diagnostics interpret a result; they do not mutate it, repair it, or invent
   a second set of statistical rules.

## Large-Problem Contracts

- `inplace=true` models and Jacobians must fill every output and remain generic
  over scalar type so bounded and constrained AD paths can use them.
- `WhiteningOperator` represents the complete static observation covariance.
  The same operator whitens residuals and Jacobian columns; it is exclusive
  with other observation-uncertainty inputs.
- Dense covariance remains the exact small/medium path with `O(n^2)` storage and
  `O(n^3)` factorization. Structured large-data support belongs in operators or
  sparse factorizations, not dense micro-optimizations.
- Fixed parameters remain part of the scientific model and must satisfy bounds.
  Profile and contour refits obey the same parameter-domain rules as the main
  fit.

## Plotting Invariants

1. `using JuFitter` must not load Makie or CairoMakie. Fitting, reports,
   diagnostics, profiles, and profile-matrix computation remain usable without
   plotting dependencies.
2. Rendering never refits or mutates a result. `plot_profile_matrix(matrix)`
   must render a precomputed matrix without repeating profile scans.
3. Use Makie's `GridLayout`, `Auto`, and content observables rather than manual
   pixel guesses or a second private layout engine. Explicit fixed dimensions
   are reserved for a declared export constraint.
4. A standard right panel has one hierarchy: legend, model, parameters, then
   statistics. `plot_info_panel!`, `plot_theme`, and `plot_palette` are the
   shared composition hooks for compound figures.
5. Plot roles are `:screen` and `:article`; light/dark appearance and explicit
   element overrides are independent. The former `:lab`, `:workbench`,
   `:modern`, and publication names are compatibility aliases, not extra designs.
   Typography, markers, axes, grid, colors, bands, legends, and report panels
   come from one central preset; compound figures consume those tokens rather
   than branching on style or appearance locally.
6. Annotation helpers validate coordinates, return Makie plot objects, and do
   not trigger a fit. Non-finite residual, pull, ratio, curve, line, or band
   coordinates fail clearly instead of disappearing in a renderer. Full-axis
   reference spans use Makie's axis-relative primitives; they must not turn
   provisional pre-render limits into data coordinates.
7. Documentation style switching uses real light/dark Makie assets. It never
   recolors or inverts one raster image with CSS.
8. Text sizes are part of the plotting contract. Validate them in the rendered
   documentation at normal browser scale, not only in full-resolution PNGs;
   wide multi-panel figures otherwise hide unreadable labels behind apparently
   successful image generation. Figure dimensions describe the intended
   logical display footprint; use `px_per_unit` for raster density instead of
   enlarging the logical canvas and shrinking it again in HTML.

## Documentation Truth

- Displayed terminal output and derived values come from the displayed code.
  Snapshot tests are the source of truth; hand-transcribed plausible output is
  not acceptable.
- Tutorial code uses visible curated arrays or small public data files. Private
  paths, course-internal names, and hidden data generation do not enter public
  pages.
- Every uncertainty band, contour level, marker, and diagnostic status states
  what it represents.
- Generated site output and gallery scratch output remain ignored. Tracked
  documentation assets must have a reproducible generator and visual review.
- Gallery scripts must write normal user output only below `examples/output/`.
  A script may update tracked documentation assets only behind an explicit
  maintainer environment flag documented beside that script.
- Public API changes require a docstring, curated API-reference update, local
  and rendered link checks, and a focused contract test.

## Change Gate

For each behavioral change:

1. Add an analytic, reference, regression, or hostile-input test that would
   have failed before the change.
2. Run the focused suite while iterating and the broader gate required by the
   affected subsystem.
3. Update public documentation for user-visible behavior, this contract for a
   durable invariant, or `RELEASE_AUDIT.md` for evidence and limitations.
4. Run `git diff --check` and keep generated output out of the commit.

Friendly examples and smoke tests can aid development, but they are not release
evidence.
