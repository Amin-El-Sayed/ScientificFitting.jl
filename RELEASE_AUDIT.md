# JuFitter Release Audit

Status: 2026-08-04

This document tracks what must be true before JuFitter should be advertised as a
serious scientific fitting library. Passing tests is necessary, but not
sufficient.

Publication policy: do not push, publish, register, deploy documentation, or
make the repository public without explicit manual approval from Amin_El_Sayed.
Local commits on `codex/*` branches are allowed only as reviewable checkpoints.

- A documentation-driven plotting audit found that `add_vband!` and
  `add_hband!` materialized Makie's provisional pre-render axis limits as data
  coordinates. An annotation added before the first render could therefore
  expand the orthogonal axis to the default `0:10` range. Both helpers now use
  Makie's axis-relative `vspan!`/`hspan!` primitives, and a focused plot
  regression verifies that pre-render spans leave automatic data limits intact.

## Current Verification

- Julia 1.10.11 previously passed the Makie-free core gate with 432 checks and
  the full package test with 501 checks in an isolated environment resolved
  under Julia 1.10. The current in-place model/Jacobian slice additionally
  passes its 18-check focused reference, and the structured-whitening slice
  passes its current 35-check reference under Julia 1.10. On the current audit
  branch, a fresh Julia-1.10-resolved temporary environment also passes the
  82-check torture suite, the 56-check diagnostic suite, and the shared
  `StatsAPI.fit` binding check. The package
  `[compat]`, core CI matrix, and full-package Linux CI matrix support both
  Julia 1.10 and Julia 1.12; the complete current matrix still needs to be
  observed remotely after the branch is pushed.
- The last completed
  `julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'` release gate
  passed on the documentation-finalization branch with 597 checks in 53m52.3s
  after the likelihood-scale API correction. `test/runtests.jl` delegates to
  the single authoritative Makie-free core inventory before adding CairoMakie
  tests. This corrected a real gate defect: the previous 509-check package
  runner omitted the structured-whitening and in-place reference files even
  though their focused suites passed separately. The current full run followed
  the public documentation source cleanup and the focused API, likelihood,
  torture, plotting, gallery, output-snapshot, and documentation gates listed
  below.
- JuFitter's low-level methods extend `StatsAPI.fit` rather than creating a
  competing generic. A local namespace check with `JuFitter`, `Distributions`,
  and `StatsAPI` loaded together confirms the shared binding and executes a
  converged `fit(problem)` call without ambiguity.
- Solver iteration counts are no longer invented when a backend does not expose
  them. `FitResult`, `LikelihoodFitResult`, and `FitReport` use `missing` for an
  unavailable count, and text reports render `iterations = unavailable` instead
  of substituting the configured iteration limit. The focused public-API
  regression covers the LsqFit path.
- Likelihood objectives and reports now use one explicit ``-2\log L``
  convention. The former `nll_min` field and `*_nll` cost symbols were
  misleading because the implementation, covariance Hessian factor, profile
  thresholds, AIC, and BIC already used the twice-negative-log-likelihood
  scale. They are now named `minus2loglik_min` and `*_likelihood`; numerical
  values and inferential thresholds are unchanged. Analytic Gaussian,
  likelihood, parameter-dependent-covariance, and structured-whitening
  references pass with 141 checks, the public API regression passes with 139
  checks, and the complete statistical suite passes with 247 checks. The
  public documentation gate rejects the retired identifiers, while the
  33-check executable-output gate verifies the new report labels against real
  example runs.
- The plotting release slice passes locally with 79 focused API/layout tests,
  232 gallery-structure checks, 1201 visual-asset checks, and 83 intentional
  PNG snapshot checks. The complete Documenter build also passes, followed by
  2020 checks against links, assets, and the responsive architecture flow in
  the rendered HTML.
- The public gallery overview now contains only reader-facing learning paths,
  example summaries, and reproducible run instructions. Asset-generation and
  release-gate procedures remain in the maintainer documentation rather than
  leaking into the scientific handbook.
- The public architecture walkthrough now describes likelihood objectives and
  stored minima on the same explicit ``-2\log L`` scale as the implementation,
  Hessian covariance, profile thresholds, and information criteria.
- Linear Calibration now identifies its controlled, deliberately imperfect
  dataset, saves directly from the complete page code, reports parameter units,
  and no longer exposes the internal gallery asset generator.
- XY Uncertainties now uses an explicit dimensioned record whose x-resolution
  contribution is scientifically relevant without overwhelming the residual
  model. Its real fit has ``\chi^2/\mathrm{ndf}=0.954`` and an `ok` diagnostic,
  and the page includes the complete saved plot rather than referring to hidden
  asset-generation code.
- Full Covariance now separates local readout variance from a finite-time
  correlated component, uses a positive physical decay-rate parameter
  consistently in code, report, and plot, and compares the full matrix with a
  deliberately incomplete diagonal treatment. The displayed output and six
  style assets come from the same fit; the full model has
  ``\chi^2/\mathrm{ndf}=0.752`` and an `ok` diagnostic, while the diagonal
  approximation visibly understates the local parameter errors.
- Compound gallery figures for Poisson counts, histogram likelihoods, damped
  oscillation, and multi-dataset fits now use the same natural-width
  `plot_info_panel!` layout contract as ordinary fit plots. Ordinary layouts do
  not guess fixed side-panel widths or resize the figure around compact report
  content; explicit fixed widths remain an opt-in export control.
- `git diff --check` passes for the current release-hardening branch.
- `julia --project=. --startup-file=no -e 'include("test/torture_runtests.jl")'`
  passes with 79 torture checks in about 1m38s on the local machine. Three
  validation checks for the former `ci_level` storage option were removed with
  that inert public keyword; the remaining suite continues to cover the
  scientific input and solver failures below. It verifies that fixed/profiled
  parameter values cannot bypass declared bounds
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
  multi-fit inputs before objective evaluation. Solver controls now fail during
  construction, incompatible explicit `backend=:lsqfit` requests cannot discard
  statistical terms or constraints, and unsupported likelihood-fit keywords no
  longer disappear silently.
- `julia --project=. --startup-file=no -e 'using Test; include("test/statistics/profile_contour_reference.jl")'`
  passes with 85 profile/contour reference checks, including validation of
  finite positive ordered contour thresholds, finite positive profile
  thresholds, default scan controls, explicit finite distinct scan grids, and
  the Makie-free `profile_matrix` diagnostic object that combines profile,
  contour, per-panel diagnostic reports, and per-panel `:ok`/`:review`/`:stop`
  status before plotting. The gate now also verifies that non-finite,
  non-symmetric, or non-positive-definite local covariance overlays in contour
  diagnostics produce an explicit `:contour_local_covariance_unavailable`
  warning instead of crashing or silently omitting the local-covariance
  comparison. The same gate covers `profile_matrix_triage(...)`, which turns
  a profile-matrix diagnostic object into a severity-ordered list of panels to
  inspect.
- `julia --project=. -e 'using Test; include("test/statistics/likelihood_reference.jl")'`
  passes with 49 likelihood reference checks after the Poisson-and-histogram
  workflow rewrite.
- `julia --project=. -e 'using Test; include("test/statistics/covariance_semantics_reference.jl")'`
  passes with 27 covariance and constraint reference checks. This includes a
  finite-difference reference for value, gradient, and Hessian of Gaussian
  ``-2\log L`` with full `cov_x` propagation, where the effective dense
  covariance depends on the fitted model parameter. It also verifies that
  asymmetric Gaussian parameter information uses a continuous, normalized
  split-normal likelihood term rather than a side-dependent normalization that
  jumps at the central value.
- `julia --project=. --startup-file=no test/statistics/runtests.jl` passes with
  247 analytic and numerical reference checks after the split-normal
  correction.
- `julia --project=. --startup-file=no -e 'using Test; include("test/statistics/diagnostics_reference.jl")'`
  passes with 56 diagnostic reference checks in about 46s, including the
  local-covariance warning that recommends profile/contour intervals when
  active bounds or strong parameter correlations make symmetric errors suspect.
  The same gate now checks that a long same-sign pull run is reported with its
  concrete point and x interval, so structured residual warnings point to a
  data region a lab user can inspect. A converged stationary maximum with
  negative local curvature now produces a critical finding, dashboard status
  `:stop`, and `NaN` uncertainty/correlation values rather than fabricated zero
  errors.
- The complete core inventory is exercised by the 597-check package gate above.
  Its roughly 42-minute runtime is too slow for the default developer loop and
  remains a release/CI gate; focused reference and torture suites provide the
  edit-test loop.
- `julia --project=docs --startup-file=no test/plots/fitplot.jl` passes with 83
  focused plot-regression checks locally, including structured-whitening error
  bars, prediction-band semantics, and the no-marginals failure path.
- `julia --project=docs --startup-file=no docs/make.jl` passes locally. The
  build output is ignored and must not be committed.
- `julia --project=. --startup-file=no test/docs_gallery_gate.jl` passes with
  257 checks. The gate enforces the current release standard for every public
  gallery page:
  scientific question, data context, model/cost explanation, complete Julia
  code, real diagnostic output, interpretation, failure modes, explicit
  one-sigma semantics, valid image assets, and complete
  `:lab`/`:modern`/`:article` light/dark plot-style coverage.
  The Photoelectric Work Function checks additionally require the centered
  line parametrization and the baseline-corrected physical slope
  `m_emit - m_base` for both Planck's constant and the work function; the same
  contract is checked in the executable asset generator.
- `julia --project=. --startup-file=no test/docs_public_release_gate.jl` passes
  with 581 checks. The gate first verifies that every page in the public
  Documenter navigation is covered, then scans those pages plus README for
  AI/placeholder wording, private local paths, author-handle leakage, and
  course-internal dataset language. It also rejects known stale public API
  identifiers such as `profile_curve`, `contour_grid`, `nll_min`, and
  `gaussian_nll`, draft/tutorial residue phrases, public image tags without
  non-empty alt text, and the ungrouped TeX operator subscripts that previously
  caused browser-side KaTeX failures. The
  first-user checks additionally require a real rendered landing-page hero,
  visible style-switchable plot assets, an explicit CairoMakie import in the
  plotting quickstart, a single report emission, and no undeployed canonical
  URL. The Reference Map checks additionally preserve the distinction between
  fit result types, additive and complete covariance descriptions, parameter
  metadata and propagated information, profile computation, and the optional
  plotting boundary.
- `julia --project=. --startup-file=no test/docs_api_reference_gate.jl` passes
  with 15 checks. The gate verifies that every exported public binding except the
  module name has a REPL/Documenter-visible docstring, preventing `@autodocs`
  from exposing undocumented public names. It also requires each public export
  to appear on the curated API reference page and protects the entry-point,
  model-signature, uncertainty, failure, custom-objective, profile, and
  optional-plotting contracts.
- `README.md` has been rewritten from an API dump into a public landing page:
  concise project purpose, quickstart, installation status, documentation entry
  points, current scope, known limitations, development gates, and release
  policy.
- The README quickstart was executed in a temporary directory with
  `julia --project=docs --startup-file=no`; it produced the PDF output file,
  converged, and returned a diagnostic dashboard with `status = ok`.
- `julia --project=. --startup-file=no test/docs_link_gate.jl` passes locally.
  The gate validates 403 local Markdown links, HTML links, image sources, and
  the static-route contract under `docs/src`. Documentation now uses the same
  flat `.html` routes locally and on static hosts; CI no longer switches raw
  gallery-card links to an incompatible pretty-URL layout.
- `julia --project=. --startup-file=no test/docs_html_link_gate.jl` passes
  after `docs/make.jl` with 2020 checks. The gate checks `href` and `src`
  targets in the rendered HTML under `docs/build`, so navigation and asset
  references are validated in the actual static site layout rather than only in
  source Markdown. It also guards the bounded, top-to-bottom architecture flow
  against the fixed-width layout that previously overflowed the article.
- The public landing page, installation page, and quickstart now form one
  verified first-user path. Installation distinguishes the local pre-release
  workflow from the future registry command, fitting/reporting from optional
  plotting dependencies, and one-time compilation from steady-state runtime.
  The quickstart uses only the explicit arrays shown to the reader, renders a
  plot on the page, prints one real report plus one real diagnostic dashboard,
  and explains the prediction-band and `report` switches. Its documented output
  is checked by the executable snapshot gate.
- The landing page distinguishes least-squares `FitResult` from likelihood
  `LikelihoodFitResult` and no longer implies that goodness-of-fit information
  exists for every statistical cost.
- The Gallery overview now follows the same eight-step progression as the
  Documenter navigation: weighted Gaussian fitting, x uncertainty, correlated
  observations, nonlinear model criticism, derived quantities, count
  likelihoods, nonlinear profile geometry, and shared-parameter hypotheses.
  The gallery gate checks card and navigation order plus every named learning
  step.
- `julia --project=. --startup-file=no test/docs_visual_asset_gate.jl` passes
  locally with 1201 checks. The gate validates documentation PNG headers,
  minimum dimensions, non-empty alt text, complete style/appearance coverage,
  consistent dimensions within each plot group, and rejects unreferenced
  gallery PNG leftovers. This is a visual-asset sanity gate; it does not
  replace future pixel-level regression testing.
- `julia --project=. --startup-file=no test/docs_visual_snapshot_gate.jl`
  passes locally with 83 checks. The gate checks SHA-256 snapshots for every
  committed documentation gallery PNG, so plot asset changes cannot enter
  unnoticed; a manifest update must be intentional after visual review.
- `.github/workflows/ci.yml` now separates the release checks into core tests on
  Linux and macOS, a full package test on Linux, and a documentation lane that
  runs source docs gates, plot regressions, `docs/make.jl`, rendered-link
  validation, and executable output snapshots. This still needs confirmation on
  GitHub Actions after the branch is pushed.
- `julia --project=. --startup-file=no test/docs_output_snapshots.jl` passes
  with 33 checks. The gate executes the documented
  quickstart/gallery example scripts with snapshot markers and verifies that
  every documented `Real output` block is an ordered subset of real script
  output. This caught and fixed stale rounded numbers, missing statistic fields,
  and a manually summarized multi-dataset output block. The test uses
  `JUFITTER_DOC_SNAPSHOT_ONLY=1`, so it verifies computations and terminal
  output without re-rendering Makie assets.
- `test/performance_budget_gate.jl` now covers representative steady-state hot
  paths with deliberately broad budgets: out-of-place and in-place 10k-point
  analytic linear least-squares, no-op bounds preserving the fast path, and
  300-point dense covariance. This is wired into the core CI lane with a runner
  scale factor. It is a regression guard, not a publishable benchmark claim.
- `julia --project=. --startup-file=no test/performance_budget_gate.jl` passes
  locally with 14 checks. This verifies the current steady-state
  budget guard and the Makie-free core extension boundary, but it remains a
  regression guard rather than publishable benchmark evidence.
- `julia --project=. --startup-file=no test/numerics/inplace_model_reference.jl`
  passes with 18 checks under Julia 1.12 and Julia 1.10. It verifies numerical
  agreement between out-of-place and `model!(out, x, p)` fits for diagonal,
  dense, and sparse covariance, native in-place LsqFit evaluation, the general
  bounded optimizer path, and an in-place full Jacobian with a fixed parameter.
  Profile refits preserve the mutating contract. Invalid signatures and
  incomplete or non-finite mutating output buffers fail before solver dispatch.
  The matched 10k-point local benchmark records both variants' time, memory,
  and allocation counts without turning one machine's ratio into a public
  performance claim.
- `test/statistics/structured_whitening_reference.jl` passes with 35 checks.
  The public `WhiteningOperator` represents a complete static observation
  covariance through `whiten!(out, residual)`, an explicit covariance log
  determinant, and optional marginal standard deviations for plotting. Dense
  AR(1) references verify fitted parameters, parameter covariance, weighted
  residuals, chi-square, normalized ``-2\log L``, AIC/BIC, bounded
  Gaussian-likelihood AD, in-place model/Jacobian evaluation, and profile
  refits. Invalid signatures,
  incomplete or non-finite output, Vector-only methods, inconsistent marginal
  dimensions, and attempts to combine independent covariance descriptions fail
  before solver dispatch.
- The benchmark suite includes `fit/structured_whitening_100000`, and the
  performance gate compares 10k and 50k allocations to guard linear scaling.
  An exploratory local BenchmarkTools run of the 100k AR(1) case recorded about
  19 ms, 119 MiB, and 628 allocations after compilation. This is internal
  evidence only, not a portable performance claim or a replacement for a saved
  baseline on selected release hardware.
- `test/benchmark_contract_gate.jl` checks the benchmark release contract
  without measuring timings: the benchmark runner keeps the required hot-path
  cases including profile-matrix diagnostics, records the required summary
  fields, and README, performance docs, pre-release checklist, and CI all point
  to the same `--project=benchmarks` workflow.
- `JUFITTER_RUN_PYTHON_INTEROP=1 julia --project=. --startup-file=no
  test/python_interop_gate.jl` passes locally with 13 checks in an isolated
  `/tmp` Python virtual environment after installing `juliacall`. The example
  keeps JuliaCall's managed Julia project active, develops this checkout into
  that environment, verifies Python access to fit parameters, `report_text`,
  and `diagnostic_dashboard_text`, and checks that fitting/reporting does not
  load Makie or CairoMakie. This also fixed the earlier SciMLBase/PythonCall
  extension warning caused by activating the JuFitter project after JuliaCall
  had started.
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
  Release comparisons also fail on mismatched Julia version, OS, CPU, machine
  target, Julia threads, BLAS threads, or units unless
  `--allow-metadata-mismatch` is explicitly used for exploratory local
  comparisons. Runner arguments are validated before dispatch: seconds must be
  finite and positive, tolerance finite and non-negative, and baseline paths
  non-empty. Local benchmark manifests and outputs are ignored because they are
  machine-specific.
- A real benchmark-runner smoke test passed locally at the corresponding
  earlier revision with
  `julia --project=benchmarks --startup-file=no benchmarks/runbenchmarks.jl
  --seconds=0.01 --save=/tmp/jufitter-benchmark-smoke.toml
  --compare=/tmp/jufitter-benchmark-smoke.toml --tolerance=0.25`. This verified
  the full save/compare path for that revision's nine benchmark cases, including
  `diagnostics/saturation_profile_matrix`, and metadata serialization. It is
  deliberately not release benchmark evidence because `--seconds=0.01` is too
  short for stable performance claims and no reference runner has been
  selected.
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
- `julia --project=. --startup-file=no test/startup_probe_gate.jl` passes with
  12 checks. It executes the startup probe with a temporary TOML output file,
  verifies the fresh-process `core_without_makie=true` output, checks the saved
  metadata, and removes the temporary artifact. The gate is wired into CI and
  the pre-release checklist.
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
  `:lab`, `:modern`, and `:article` plot assets for light and dark
  documentation appearances.
- The plot style contracts were tightened: `:lab` is now a plain sans analysis
  style with shorter error-bar whiskers, `:modern` has a distinct but
  restrained scientific color hierarchy, and the photoelectric custom plot
  now derives theme, palette, marker sizes, whisker widths, line widths, and
  report typography from the public plotting style API.
- The documentation style selector now covers the compound custom gallery
  figures as well as ordinary `plot_fit` figures: Poisson counts, histogram
  likelihoods, constraints/profile/contour plots, damped oscillator, and
  multi-dataset calibration all render real Makie assets for
  `:lab`, `:modern`, and `:article` in light and dark appearances.
- The `How JuFitter Works` page now follows the implemented architecture in a
  bounded top-to-bottom flow: typed problem construction, validation, objective
  construction, compatible solver dispatch, result construction, and optional
  post-fit analysis. It distinguishes Gaussian and likelihood problem types,
  additive parameter information from bounds and fixed parameters, and primary
  fitting from profile/contour refits. Desktop, narrow mobile, and dark-mode
  browser checks show no horizontal overflow; the same mobile check caught and
  fixed an overflowing documentation toolbar.
- `Fitting for Practitioners` now begins with the measurement process and maps
  continuous observations, correlated Gaussian vectors, counts, histograms,
  and unbinned samples to the corresponding statistical interfaces. It uses
  explicit numerical examples to show why uncertainty scale and covariance
  change the cost, states the limits of first-order x-error propagation, and
  replaces the ndf-independent reduced-chi-square rule with the correct
  `sqrt(2 * ndf)` scale. Profile and contour guidance uses the implemented
  `ProfileInterval.profile_result` field and explains what asymmetric, clipped,
  open, or non-elliptic regions require in practice. Desktop light/dark browser
  checks caught and removed unreadable wide tables instead of masking them with
  page-specific CSS. A focused linear-fit run executes the documented
  `profile_interval`, profile diagnosis, contour, and contour-diagnosis access
  patterns successfully.
- Technical maintenance pages remain in the rendered documentation, but they
  are nested under `Reference > Technical Notes` rather than appearing as a
  top-level user path. The public documentation hygiene gate now prevents a
  top-level `Engineering Notes` navigation block from returning.
- Browser QA against a temporary local Documenter server verifies that the
  Poisson/histogram page loads, exposes 12 grouped plot assets, and switches
  the visible dark-mode image sources between `:article` and `:lab`
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
- Focused plot regression tests pass with 79 checks after the plot-style
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
- `ProfileMatrixResult` now retains the selected best-fit values and local
  covariance geometry together with its profile and contour scans.
  `plot_profile_matrix(matrix_result)` therefore renders a previously computed
  diagnostic without repeating any profile refits; the convenience method for
  a fit result computes once and delegates to the same renderer.
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
- The scoped v0 core is locally through its focused numerical gates. Friendly
  examples and smoke tests are not accepted as evidence; package, CI, and
  documentation gates below remain mandatory before publication.
- `Statistical Foundations` now forms a complete method chapter rather than a
  formula inventory. It derives Gaussian, correlated, x-error, Poisson,
  histogram, unbinned, and extended likelihood models; distinguishes local
  covariance from profile coverage; states Wilks and p-value limits; and
  constrains AIC/BIC comparisons to compatible normalized likelihoods. Light-
  and dark-mode browser checks confirm that equations and threshold tables
  render correctly.
- The API reference now provides an observation-model entry-point table,
  allocating and in-place model contracts, complete uncertainty and
  parameter-control semantics, solver and failure behavior, result-field
  definitions, likelihood assumptions, profile thresholds, report selection,
  and the optional CairoMakie boundary before the generated method docstrings.
  The page builds with method-specific `fit` documentation, passes light/dark
  browser review, and has no unresolved local or rendered links. The unused
  `ci_level` fit keyword was removed rather than presenting a scientifically
  inert option; profile and contour coverage remains explicit through
  `threshold` and `levels` where those choices actually affect computation.

## Release Blockers

- Every page in the configured public navigation has completed a page-level
  source and contract pass. This is not a substitute for final human visual and
  subject-matter review, but private planning pages, the stale hosted roadmap,
  and the generic research-landscape draft no longer leak into the generated
  site.
- The API reference and Reference Map have completed their page-level editorial
  and contract audits. The Reference Map now follows the actual public
  boundaries between observation models, additive and complete covariance
  descriptions, parameter information, `FitResult` and
  `LikelihoodFitResult`, local and profile uncertainty, and optional Makie
  output. The technical-note scope is also resolved: architecture and
  performance remain public, while the redundant public maintenance page was
  replaced by a compact root-level `MAINTAINERS.md` contract. Planning,
  evidence, release commands, and durable invariants now each have one
  authoritative home.
- `docs/src` now contains only Markdown pages listed in the public Documenter
  navigation. The durable editorial contract moved to root-level
  `DOCUMENTATION.md`; the root `ROADMAP.md` remains the only roadmap, and the
  redundant hidden research-landscape page was removed. The public hygiene gate
  enforces this one-to-one source/navigation boundary.
- Source and rendered documentation links are covered locally by
  `test/docs_link_gate.jl` and `test/docs_html_link_gate.jl`, and both are wired
  into `.github/workflows/ci.yml`. Remote CI execution still needs to be
  observed after pushing the branch.
- Documentation gallery PNGs now have byte-level snapshot coverage. This catches
  unintentional asset drift, but human visual review is still required for
  intentional plot style or renderer changes.
- There is no documentation deployment workflow yet.
- The repository still needs an explicitly chosen license and release citation
  metadata before publication. Add the approved `LICENSE` and `CITATION.cff`
  only after the maintainer has selected the license and reviewed the citation
  fields; do not infer that policy from the source code.
- Browser screenshots can still time out on very large documentation pages.
  Use DOM/computed-style checks plus targeted static image inspection until a
  stable visual-regression workflow exists.

## Known Numerical Limits And Post-v0 Work

- Plotting now uses a CairoMakie package extension, so fitting/reporting no
  longer requires Makie at load time. Plotting still has the expected
  CairoMakie first-use compilation cost, and the docs/examples environment must
  keep `CairoMakie` as an explicit dependency.
- Cold-start precompilation still needs a dedicated benchmark gate. The current
  evidence proves that Makie is not loaded by the fitting core, but it does not
  yet quantify the remaining Optimization/SciML first-use cost across clean
  environments.
- Hostile-input coverage is broad enough to support a scoped v0, but it is not
  evidence that every scientific data shape is handled. Continue extending it
  when new constraint, likelihood, profile/contour, multi-dataset, scaling, or
  covariance failure modes are found; public wording must stay scoped to the
  workflows represented by reference and torture tests.
- The following items are not automatic v0 blockers if their limitations are
  documented honestly, but they are release-audit candidates for every serious
  public claim about scale, nonlinear uncertainty, or x-covariance derivatives:
  dense covariance asymptotics, local parameter-covariance validity, and
  parameter-dependent dense `cov_x` automatic differentiation.
- x-uncertainty propagation now has a public vectorized
  `x_derivative=(x, p) -> dy_dx` hook, which avoids the default point-by-point
  AD derivative when the model derivative is known. Large least-squares models
  can additionally use the in-place model/Jacobian contract; the default x
  derivative remains pointwise for API simplicity.
- Dense covariance support is mathematically useful, but still `O(n^2)` memory
  and `O(n^3)` factorization. Large correlated datasets with a known complete
  static covariance can now use an application-specific `WhiteningOperator`.
- Static sparse `cov_y` now remains sparse for the unbounded least-squares
  backend. Validation and whitening use sparse CHOLMOD factors without
  materializing a dense covariance, and the benchmark runner includes
  `fit/sparse_covariance_5000`. Sparse covariance in constrained or
  Gaussian-likelihood optimizer paths remains a documented limitation because CHOLMOD
  solves do not propagate `ForwardDiff` dual numbers.
- `WhiteningOperator` supplies matrix-free static residual and Jacobian
  whitening plus the covariance log determinant required for normalized ``-2\log L`` and
  information criteria. It is deliberately a low-level scientific contract:
  JuFitter can validate its interface and finite output, but cannot prove that
  a custom operator satisfies `W'W = inv(C)` or that its determinant is correct.
  Built-in banded, Toeplitz, low-rank-plus-diagonal, sparse-precision,
  parameter-dependent, and structured x-covariance operators remain future
  work.
- General constrained and parameter-dependent-covariance objectives still
  allocate through `ForwardDiff` Jacobians/Hessians and dual-number output
  buffers. The in-place model contract removes avoidable user-model and LsqFit
  allocations, but matrix-free residual objectives and in-place likelihood
  evaluation remain future large-scale work.
- The parameter-dependent dense covariance audit found and fixed a real AD
  defect: Cholesky validation stripped `ForwardDiff` dual information before
  factorization. Validation now strips duals only for finite/symmetry checks,
  while the factorization itself preserves AD values. A focused finite-
  difference regression test covers Gaussian ``-2\log L`` value, gradient, and
  Hessian for dense `cov_x` propagation. This is not a substitute for future
  structured covariance operators, broader `cov_x` AD audits, or large-scale
  performance work.
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
  before release-grade claims. Invalid local covariance overlays are now
  surfaced explicitly as diagnostic findings, so a profile contour can remain
  usable even when the symmetric local ellipse is not.
- Multi-dataset fitting currently supports useful parameter mapping, but not the
  full uncertainty model expected from kafe2-level workflows.
- Python interoperability has a local clean-environment validation through
  JuliaCall/PythonCall, but it still needs CI confirmation before public release
  claims. The Makie-free example and opt-in gate
  (`examples/python/fit_from_python.py`, `test/python_interop_gate.jl`) check
  that Python can access fitted parameters, `report_text`,
  `diagnostic_dashboard_text`, and the fitting/reporting path without loading
  Makie or CairoMakie. Broader NumPy conversion and ownership semantics remain
  future documentation/test work beyond the minimal plain-array path.
- The local performance-budget gate is wired into CI, but it has not yet been
  observed on GitHub Actions and does not replace saved `BenchmarkTools`
  baselines.
- The benchmark runner can save and compare baselines, but no release reference
  hardware or CI baseline has been selected yet.
- `diagnostic_dashboard(...)` now summarizes structured findings into status,
  severity counts, and deduplicated next actions. Profile contours now default
  to directly labeled 1-sigma/2-sigma regions with local-covariance comparison;
  residual/pull diagnostics now identify long same-sign pull intervals, but the
  combined visual diagnostic dashboard is still not yet at the level needed to
  compete with kafe2/Minuit-style workflows.
- The optional CairoMakie implementation remains concentrated in a roughly
  2100-line rendering file. The public extension boundary and API are compact
  and tested, so a mechanical split is not required for numerical correctness,
  but future plot features should not add another private layout framework.
  Split style/layout, fit rendering, and diagnostic rendering only as a
  behavior-preserving maintenance change with the same visual gates.

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
- Public technical notes are limited to architecture and performance under
  `Reference > Technical Notes`. Durable contributor invariants live in
  `MAINTAINERS.md`; release status and commands remain in `RELEASE_AUDIT.md` and
  `RELEASE_CHECKLIST.md` rather than leaking a changing maintenance log into the
  user website.
- Private/local dataset language is now covered by the public documentation
  hygiene gate; keep extending that gate when new release-language risks are
  identified.
- Dark-mode plots must be checked visually whenever a gallery asset changes.

## CI And Packaging Blockers

- Choose and add the release license and `CITATION.cff`; verify package name,
  UUID, authorship, repository URL, and version metadata against the exact
  repository that will be registered.
- Julia 1.10 is the intentional support floor. Local Julia 1.10 evidence covers
  the earlier 432-check core gate and 501-check full package suite plus the
  current 18-check in-place and 35-check structured-whitening reference slices;
  confirm the complete current 1.10/1.12 core and package matrix on GitHub
  Actions after pushing.
- Add a docs-deploy job for GitHub Pages or the chosen static host.
- Confirm the new core/package/docs CI lanes on GitHub Actions after pushing.
- Run the Python interoperability release gate in CI. The local opt-in gate now
  passes in an isolated `juliacall` environment, but public release claims need
  the same path observed on the target CI or release machine.
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
