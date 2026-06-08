# JuFitter Release Audit

Status: 2026-06-08

This document tracks what must be true before JuFitter should be advertised as a
serious scientific fitting library. Passing tests is necessary, but not
sufficient.

Publication policy: do not push, publish, register, deploy documentation, or
make the repository public without explicit manual approval from Amin_El_Sayed.
Local commits on `codex/*` branches are allowed only as reviewable checkpoints.

## Current Verification

- `git diff --check` passes for the current release-hardening branch.
- `julia --project=. --startup-file=no -e 'include("test/torture_runtests.jl")'`
  passes with 64 torture checks in about 65s on the local machine after adding
  validation that fixed/profiled parameter values cannot bypass declared bounds
  and that `p0` and user-provided `initial_guesses` are finite and inside
  declared bounds instead of being silently clipped before optimization. The
  same suite now checks that parameter priors, correlated parameter
  constraints, fixed-parameter metadata, and likelihood-fit parameter metadata
  reject non-finite or nonsymmetric inputs before solver dispatch. It also
  checks that likelihood observations, histogram edges, domains, quadrature
  tolerances, and expected-count output lengths fail before optimizer internals
  can turn invalid scientific input into vague convergence failures. Indexed
  and multi-dataset Gaussian wrappers now reject non-finite observations,
  non-positive or non-finite `sigma_y`, invalid indexed `cov_y`, and empty
  multi-fit inputs before objective evaluation.
- `julia --project=. -e 'using Test; include("test/statistics/profile_contour_reference.jl")'`
  passes with 47 profile/contour reference checks, including validation of
  finite positive ordered contour thresholds, finite positive profile
  thresholds, default scan controls, and explicit finite distinct scan grids.
- `julia --project=. -e 'using Test; include("test/statistics/likelihood_reference.jl")'`
  passes with 49 likelihood reference checks after the Poisson-and-histogram
  workflow rewrite.
- `julia --project=. -e 'using Test; include("test/statistics/covariance_semantics_reference.jl")'`
  passes with 23 covariance and constraint reference checks. This includes a
  finite-difference reference for value, gradient, and Hessian of a Gaussian NLL
  with full `cov_x` propagation, where the effective dense covariance depends on
  the fitted model parameter.
- `julia --project=. -e 'using Test; include("test/statistics/diagnostics_reference.jl")'`
  passes with 45 diagnostic reference checks in about 38s, including the
  local-covariance warning that recommends profile/contour intervals when
  active bounds or strong parameter correlations make symmetric errors suspect.
- `julia --project=. --startup-file=no -e 'include("test/core_runtests.jl")'`
  passes with 386 core checks in about 9m56s on the local machine after the
  CairoMakie plotting extension split. This is too slow for the default
  developer gate and must be split further before release CI is finalized.
- `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'` passes with
  439 checks in about 13m15s for the test phase after package test
  precompilation. This verifies that test extras load the CairoMakie extension
  correctly.
- `julia --project=docs --startup-file=no test/plots/fitplot.jl` passes with 64
  focused plot-regression checks in about 54s.
- `julia --project=docs --startup-file=no docs/make.jl` passes locally. The
  build output is ignored and must not be committed.
- `julia --project=. test/docs_gallery_gate.jl` passes with 225 checks. The
  gate enforces the current release standard for every public gallery page:
  scientific question, data context, model/cost explanation, complete Julia
  code, real diagnostic output, interpretation, failure modes, explicit
  one-sigma semantics, valid image assets, and complete
  `:workbench`/`:showcase`/`:publication` light/dark plot-style coverage.
- `julia --project=. --startup-file=no test/docs_public_release_gate.jl` passes
  with 183 checks. The gate first verifies that every page in the public
  Documenter navigation is covered, then scans those pages plus README for
  AI/placeholder wording, private local paths, author-handle leakage, and
  course-internal dataset language. It also rejects known stale public API
  identifiers such as `profile_curve` and `contour_grid`.
- `julia --project=. --startup-file=no test/docs_api_reference_gate.jl` passes
  locally. The gate verifies that every exported public binding except the
  module name has a REPL/Documenter-visible docstring, preventing `@autodocs`
  from exposing undocumented public names. It also requires each public export
  to appear on the curated API reference page.
- `README.md` has been rewritten from an API dump into a public landing page:
  concise project purpose, quickstart, installation status, documentation entry
  points, current scope, known limitations, development gates, and release
  policy.
- The README quickstart was executed in a temporary directory with
  `julia --project=docs --startup-file=no`; it produced the PDF output file,
  converged, and returned a diagnostic dashboard with `status = ok`.
- `julia --project=. --startup-file=no test/docs_link_gate.jl` passes locally.
  The gate validates 307 local Markdown links, HTML links, and image sources
  under `docs/src`, including `.html` links that should resolve to source
  `.md` pages before Documenter renders them.
- `julia --project=. --startup-file=no test/docs_html_link_gate.jl` passes
  after `docs/make.jl`. The gate checks `href` and `src` targets in the rendered
  HTML under `docs/build`, so navigation and asset references are validated in
  the actual static site layout rather than only in source Markdown.
- `julia --project=. --startup-file=no test/docs_visual_asset_gate.jl` passes
  locally. The gate validates documentation PNG headers, minimum dimensions,
  non-empty alt text, complete style/appearance coverage, consistent dimensions
  within each plot group, and rejects unreferenced gallery PNG leftovers. This
  is a visual-asset sanity gate; it does not replace future pixel-level
  regression testing.
- `julia --project=. --startup-file=no test/docs_visual_snapshot_gate.jl`
  passes locally. The gate checks SHA-256 snapshots for every committed
  documentation gallery PNG, so plot asset changes cannot enter unnoticed; a
  manifest update must be intentional after visual review.
- `.github/workflows/ci.yml` now separates the release checks into core tests on
  Linux and macOS, a full package test on Linux, and a documentation lane that
  runs source docs gates, plot regressions, `docs/make.jl`, rendered-link
  validation, and executable output snapshots. This still needs confirmation on
  GitHub Actions after the branch is pushed.
- `julia --project=. --startup-file=no test/docs_output_snapshots.jl` passes
  with 33 checks in about 2m02s. The gate executes the documented
  quickstart/gallery example scripts with snapshot markers and verifies that
  every documented `Real output` block is an ordered subset of real script
  output. This caught and fixed stale rounded numbers, missing statistic fields,
  and a manually summarized multi-dataset output block. The test uses
  `JUFITTER_DOC_SNAPSHOT_ONLY=1`, so it verifies computations and terminal
  output without re-rendering Makie assets.
- `test/performance_budget_gate.jl` now covers representative steady-state hot
  paths with deliberately broad budgets: 10k-point analytic linear
  least-squares, no-op bounds preserving the fast path, and 300-point dense
  covariance. This is wired into the core CI lane with a runner scale factor.
  It is a regression guard, not a publishable benchmark claim.
- `test/benchmark_contract_gate.jl` checks the benchmark release contract
  without measuring timings: the benchmark runner keeps the required hot-path
  cases, records the required summary fields, and README, performance docs,
  pre-release checklist, and CI all point to the same `--project=benchmarks`
  workflow.
- `RELEASE_CHECKLIST.md` now defines the local pre-release gate: clean
  repository state, core/package/statistical/torture tests, documentation gates,
  docs build and rendered-link validation, output snapshots, plot regressions,
  performance guard, benchmark evidence, optional Python interop, CI status,
  Git identity, and manual publication approval. `test/release_checklist_gate.jl`
  verifies that the checklist keeps these required commands and safeguards.
- The pre-release checklist now requires explicit review of three scientific
  limitation classes before public claims: dense covariance scaling versus
  structured covariance/custom whitening, local parameter-covariance validity
  versus profile/contour intervals, and finite-difference AD references for
  parameter-dependent dense `cov_x` propagation.
- `benchmarks/runbenchmarks.jl` can now write TOML baselines with `--save` and
  compare later runs with `--compare`. Saved summaries include Julia, JuFitter,
  OS, CPU, Julia-thread, BLAS-thread, git-commit, timing, memory, and allocation
  metadata. Baseline comparison fails if benchmark cases are missing from either
  side, so benchmark-set drift cannot pass as a valid performance comparison.
  Runner arguments are validated before dispatch: seconds must be finite and
  positive, tolerance finite and non-negative, and baseline paths non-empty.
  Local benchmark manifests and outputs are ignored because they are
  machine-specific.
- A real benchmark-runner smoke test passed locally with
  `julia --project=benchmarks --startup-file=no benchmarks/runbenchmarks.jl
  --seconds=0.01 --save=/tmp/jufitter-benchmark-smoke.toml
  --compare=/tmp/jufitter-benchmark-smoke.toml --tolerance=0.25`. This verified
  the full save/compare path for the current eight benchmark cases and metadata
  serialization. It is deliberately not release benchmark evidence because
  `--seconds=0.01` is too short for stable performance claims and no reference
  runner has been selected.
- `benchmarks/startup_probe.jl` now provides a dedicated startup smoke path for
  the core package. It starts a fresh Julia process, loads `JuFitter`, verifies
  that neither `Makie` nor `CairoMakie` was loaded, and can write a TOML summary
  of elapsed wall time and metadata. This supports the release claim that
  fitting/reporting stays Makie-free without turning startup timing into an
  unreviewed benchmark claim.
- `julia --project=. --startup-file=no benchmarks/startup_probe.jl
  --save=/tmp/jufitter-startup-probe.toml` passed locally. The fresh-process
  output included `loaded_plot_modules=` and `core_without_makie=true`, then the
  temporary TOML artifact was removed.
- `test/startup_probe_gate.jl` now executes the startup probe with a temporary
  TOML output file, verifies the fresh-process `core_without_makie=true`
  output, checks the saved metadata, and removes the temporary artifact. The
  gate is wired into CI and the pre-release checklist.
- `julia --project=. examples/gallery/10_multi_dataset_calibration.jl` prints
  the same diagnostic-dashboard sections that the Multi-Dataset gallery page
  shows. This fixed a documentation/example sync defect where the page
  contained dashboard output that was not produced by the executable example.
- `julia --project=. --startup-file=no -e 'using JuFitter; ...'` verifies that
  fitting and `report_text(...)` work without loading `CairoMakie` or `Makie`.
  The public plotting API now lives behind the optional CairoMakie package
  extension and reports a clear error if plot functions are called before
  `using CairoMakie`.
- `julia --project=docs --startup-file=no -e 'using JuFitter, CairoMakie; ...'`
  verifies that the CairoMakie extension loads and `plot_fit(...)` returns a
  Makie `Figure` with the expected default footprint.
- `julia --project=docs docs/make.jl` passes after the strict gallery gate was
  added and the gallery pages were brought to the gated structure. The generated
  `docs/build` directory was removed after verification.
- The documentation build passes after adding the `How JuFitter Works` page,
  reordering the gallery path, replacing beginner-facing generated data blocks
  on the first pages with explicit arrays, and adding the asset-based plot-style
  selector.
- Documentation examples that print fit reports or diagnostic dashboards now
  show compact notebook-style output blocks directly under the corresponding
  code cells. These output blocks are guarded by the executable snapshot gate,
  so stale or manually invented terminal output fails tests instead of relying
  on editorial review.
- The documentation toolbar now exposes theme and plot-style selectors directly
  instead of hiding light/dark mode behind Documenter's settings menu.
- The unused Documenter settings and article-toggle controls are removed from
  the documentation toolbar; light/dark and plot-style selection now use the
  JuFitter controls directly.
- Diagnostic dashboard status labels now use reader-facing actions:
  `ok - no immediate issue`, `review - inspect diagnostics`, and
  `critical - fix before use`. Critical summaries now say to fix the issue
  before using the result for conclusions instead of using publication-specific
  phrasing.
- The damped-oscillator gallery now uses a stated pointwise statistical
  uncertainty scale, keeps the constant-frequency model as a rejected
  model-criticism example, and shows a frequency-drift model whose real
  diagnostic dashboard reaches `ok`.
- The photoelectric gallery and executable example now present the threshold
  plot as a modular JuFitter/Makie extension: the fit results remain the source
  of truth, while the intersection band, marker, legend, and right-side report
  use the public style, palette, annotation, and information-panel APIs.
- Public documentation navigation now separates user documentation from
  engineering notes. Internal roadmap, documentation-plan, and research
  landscape pages remain in the repository but are no longer part of the
  public Documenter navigation.
- The XY-uncertainty and full-covariance gallery pages now use explicit
  notebook-style measured arrays in public code blocks, keep uncertainty-model
  construction separate from data generation, and provide real Makie-rendered
  `:workbench`, `:showcase`, and `:publication` plot assets for light and dark
  documentation appearances.
- The plot style contracts were tightened: `:workbench` is now a plain sans
  analysis style with shorter error-bar whiskers, `:showcase` has a distinct
  but restrained scientific color hierarchy, and the photoelectric custom plot
  now derives theme, palette, marker sizes, whisker widths, line widths, and
  report typography from the public plotting style API.
- The documentation style selector now covers the compound custom gallery
  figures as well as ordinary `plot_fit` figures: Poisson counts, histogram
  likelihoods, constraints/profile/contour plots, damped oscillator, and
  multi-dataset calibration all render real Makie assets for
  `:workbench`, `:showcase`, and `:publication` in light and dark appearances.
  The `How JuFitter Works` page now shows a flat left-to-right pipeline with
  explicit validation, cost-construction, solver-dispatch, post-fit-analysis,
  and output stages rather than a loose component inventory.
- Technical maintenance pages remain in the rendered documentation, but they
  are nested under `Reference > Technical Notes` rather than appearing as a
  top-level user path. The public documentation hygiene gate now prevents a
  top-level `Engineering Notes` navigation block from returning.
- Browser QA against a temporary local Documenter server verifies that the
  Poisson/histogram page loads, exposes 12 grouped plot assets, and switches
  the visible dark-mode image sources between `:publication` and `:workbench`
  correctly.
- Static gallery image reference validation passes for all Markdown sources
  under `docs/src`.
- The full core suite is now separated from plot rendering. It checks the
  optional plotting boundary by asserting that plotting calls without
  `using CairoMakie` fail with a clear extension-loading message.
- The exported plotting API stubs now carry public docstrings at the optional
  CairoMakie extension boundary. This prevents `@autodocs` from exposing bare
  placeholder functions for `plot_fit`, `fitplot`, annotation helpers, style
  helpers, and profile/contour plotting entry points.
- Focused plot regression tests pass with 64 checks after the plot-style
  architecture was reduced to three explicit contracts, light/dark appearance
  was separated from style, and the right-side report became a reusable,
  left-aligned layout component. The same gate also covers compatibility
  aliases, public theme/palette access, invalid style combinations, and
  style-aware profile, contour, residual, and diagnostic plots. It now also
  verifies that a right-side report cannot collapse the requested output
  footprint, that users can add Makie layout elements after `plot_fit` returns
  without destabilizing the figure, and that `plot_profile_matrix` exports a
  multi-parameter overview while rejecting invalid parameter selections. The
  same gate now covers post-fit annotation helpers for curves, points,
  vertical/horizontal reference lines, and vertical/horizontal bands, including
  non-finite annotation inputs and too-short curves. Residual, pull, and ratio
  diagnostic plots now reject non-finite displayed values, and ratio
  diagnostics fail clearly when a model prediction is zero.
- The damped-oscillator compound figure now uses the public shared theme,
  palette, and information-panel APIs instead of a private one-off layout.
- Documentation plots expand rightward from the article column rather than
  centering across the viewport, preventing wide figures from overlapping the
  fixed navigation sidebar.
- Contour plots changed
  from a heatmap-first default to labeled profile regions with an optional local
  covariance overlay. Profile plots now label the actual scan, local parabolic
  approximation, and interval threshold, and can focus the displayed delta-cost
  range without changing the scan. The plot path also covers failed-refit cells,
  invalid display limits, and fully non-finite contour surfaces.
- Phase 3 is reopened as the release gate. The project no longer treats
  friendly-path examples or smoke tests as evidence of production robustness.

## Release Blockers

- Public documentation outside the gated gallery still needs a page-by-page
  editorial pass before broad promotion. The legacy German mathematical audit
  was removed from the public
  documentation source, the main technical concept pages are now in polished
  English, and the public reference overview has been rewritten as a workflow
  map with current API names.
- The API reference is generated from docstrings, and every exported public
  binding now has baseline REPL/Documenter documentation. The public API page
  is curated by workflow area instead of relying on an unstructured
  `@autodocs` dump. The next release polish step is richer
  parameter-by-parameter reference entries for the highest-traffic workflows.
- Source and rendered documentation links are covered locally by
  `test/docs_link_gate.jl` and `test/docs_html_link_gate.jl`, and both are wired
  into `.github/workflows/ci.yml`. Remote CI execution still needs to be
  observed after pushing the branch.
- Documentation gallery PNGs now have byte-level snapshot coverage. This catches
  unintentional asset drift, but human visual review is still required for
  intentional plot style or renderer changes.
- There is no documentation deployment workflow yet.
- Browser screenshots can still time out on very large documentation pages.
  Use DOM/computed-style checks plus targeted static image inspection until a
  stable visual-regression workflow exists.

## Code And Numerical Blockers

- Plotting now uses a CairoMakie package extension, so fitting/reporting no
  longer requires Makie at load time. Plotting still has the expected
  CairoMakie first-use compilation cost, and the docs/examples environment must
  keep `CairoMakie` as an explicit dependency.
- Cold-start precompilation still needs a dedicated benchmark gate. The current
  evidence proves that Makie is not loaded by the fitting core, but it does not
  yet quantify the remaining Optimization/SciML first-use cost across clean
  environments.
- The torture suite is only starting. It must cover constraints, priors,
  likelihoods, profiles, contours, multi-dataset fits, bad scaling, local
  minima, invalid uncertainty models, and large datasets before release claims
  can use words like robust.
- The following items are not automatic v0 blockers if their limitations are
  documented honestly, but they are release-audit candidates for every serious
  public claim about scale, nonlinear uncertainty, or x-covariance derivatives:
  dense covariance asymptotics, local parameter-covariance validity, and
  parameter-dependent dense `cov_x` automatic differentiation.
- x-uncertainty propagation now has a public vectorized
  `x_derivative=(x, p) -> dy_dx` hook, which avoids the default point-by-point
  AD derivative when the model derivative is known. The default path remains
  pointwise for API simplicity, and very large workflows still need future
  in-place model/residual APIs.
- Dense covariance support is mathematically useful, but still `O(n^2)` memory
  and `O(n^3)` factorization. Large correlated datasets need structured
  covariance operators or custom whitening operators.
- The first release may keep dense covariance as the explicit exact small/medium
  path, but future work must add structured covariance or whitening operators
  for long time series, images, spectra, and detector arrays. That future path
  must define matrix-free whitening, log-determinant handling where needed, and
  diagnostics that tell users when dense covariance is the wrong model.
- Objective functions still allocate substantially through `ForwardDiff`
  Jacobians/Hessians and array-returning model calls. Very large datasets need
  an in-place residual/model API.
- The parameter-dependent dense covariance audit found and fixed a real AD
  defect: Cholesky validation stripped `ForwardDiff` dual information before
  factorization. Validation now strips duals only for finite/symmetry checks,
  while the factorization itself preserves AD values. A focused finite-
  difference regression test covers Gaussian-NLL value, gradient, and Hessian
  for dense `cov_x` propagation. This is not a substitute for future structured
  covariance operators, broader `cov_x` AD audits, or large-scale performance
  work.
- Full parameter-dependent dense `cov_x` remains an audit-sensitive path.
  Simple diagonal `sigma_x` propagation has clearer AD semantics; dense
  effective covariance can lose derivative information if future validation,
  conversion, or factorization code strips `ForwardDiff` duals too early.
  Release and refactor gates should keep finite-difference value/gradient/
  Hessian references for these paths.
- Invalid covariance matrices must not be silently repaired. If a future
  regularization/jitter policy is added, it must be explicit in the API and
  visible in diagnostics.
- Bounds and constraints currently rely on local Hessian/covariance
  approximations. Reports must be explicit about when errors are local and when
  profile intervals are required.
- Parameter covariance remains a local approximation. Nonlinear models, weak
  data, active bounds, and asymmetric likelihoods require profiles/contours for
  credible intervals. Diagnostics now add an explicit profile/contour
  recommendation for active bounds, ill-conditioned local covariance/Hessian
  estimates, and strong parameter correlations, but future work should detect
  nonlinearity and asymmetric likelihoods more directly and expose asymmetric
  interval summaries where the profile likelihood supports them.
- Profile and contour scans now expose failed refits through diagnostics and
  support adaptive refinement around profile thresholds and contour levels.
  Strongly curved/non-elliptic contours still need deeper diagnostic polish
  before release-grade claims.
- Multi-dataset fitting currently supports useful parameter mapping, but not the
  full uncertainty model expected from kafe2-level workflows.
- Python interoperability is not yet release-validated. The intended path is
  JuliaCall/PythonCall calling JuFitter from Python. A Makie-free example and an
  opt-in gate exist (`examples/python/fit_from_python.py`,
  `test/python_interop_gate.jl`), but this machine currently has Python
  3.12.13 without `juliacall`. Release claims still require running the gate in
  a clean Python environment and checking array conversion semantics. The
  install page documents the current policy: Python use is experimental or
  deferred unless `JUFITTER_RUN_PYTHON_INTEROP=1` passes with `juliacall`.
- The local performance-budget gate is wired into CI, but it has not yet been
  observed on GitHub Actions and does not replace saved `BenchmarkTools`
  baselines.
- The benchmark runner can save and compare baselines, but no release reference
  hardware or CI baseline has been selected yet.
- `diagnostic_dashboard(...)` now summarizes structured findings into status,
  severity counts, and deduplicated next actions. Profile contours now default
  to directly labeled 1-sigma/2-sigma regions with local-covariance comparison;
  residual/pull views and the combined diagnostic dashboard are still not yet
  at the level needed to compete with kafe2/Minuit-style workflows.

## Documentation Blockers

- The gallery now has structural, public-hygiene, output-snapshot, link,
  rendered-link, visual-asset sanity, and byte-level visual-snapshot gates.
  Human visual review is still required whenever plot assets, themes, or page
  CSS change.
- The executable output-snapshot gate is intentionally strict, but still slow
  enough to belong in a documentation/release CI lane rather than the default
  edit-test loop. It now avoids Makie asset rendering via
  `JUFITTER_DOC_SNAPSHOT_ONLY=1`; further speedups require reducing expensive
  fits or caching documented snapshot computations.
- The math section now starts with a beginner path, model-choice table, and
  minimal mental model before the formal likelihood/covariance reference. It
  still needs human subject-matter review before broad promotion.
- Technical maintenance notes remain accessible under `Reference > Technical
  Notes`; they no longer dominate the top-level user navigation.
- Private/local dataset language is now covered by the public documentation
  hygiene gate; keep extending that gate when new release-language risks are
  identified.
- Dark-mode plots must be checked visually whenever a gallery asset changes.

## CI And Packaging Blockers

- Add a docs-deploy job for GitHub Pages or the chosen static host.
- Confirm the new core/package/docs CI lanes on GitHub Actions after pushing.
- Run the Python interoperability release gate in CI or explicitly document why
  Python support is deferred from v0. The local opt-in gate exists, and the
  install page documents the deferred/experimental status, but the gate has not
  yet passed in a clean `juliacall` environment.
- Select release-reference hardware or CI runners and save benchmark baselines
  for hot paths with `benchmarks/runbenchmarks.jl --save=...`.
- Run `RELEASE_CHECKLIST.md` from a clean checkout before any public release,
  registration, documentation deployment, or announcement.

## Minimum Public v0 Criteria

- All tests pass on macOS and Linux CI.
- Docs build and link check pass in CI.
- README, install, quickstart, and gallery are polished enough for first-time
  users.
- At least eight gallery workflows are complete stories with formulas, full
  code, generated plots, and interpretation.
- Core statistical semantics are covered by reference tests.
- Known limitations are explicit and not hidden in optimistic marketing text.
- Performance claims are backed by reproducible benchmarks.
