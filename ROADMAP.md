# JuFitter Roadmap

JuFitter is developed as a scientific fitting package with three non-negotiable
goals:

1. Beautiful, reliable plots from minimal input.
2. Statistically explicit models, costs, uncertainties, and intervals.
3. Fast, robust numerics with clear diagnostics when assumptions fail.

## Strategic Positioning

The external ecosystem has strong individual tools, but no single Julia package
currently combines a kafe2-like statistical model, Minuit-style profiling,
Makie-quality plotting, and a modern documentation experience.

JuFitter should compete by combining:

- the one-liner convenience expected from SciPy/lmfit workflows,
- the explicit covariance and constraint semantics expected from kafe2/Minuit,
- Julia-native performance and composability,
- a documented Python interoperability path through JuliaCall/PythonCall for
  users who want to call the fitting/reporting engine from Python without
  rewriting JuFitter,
- publication-grade plots as a default outcome, not a later styling step,
- documentation that teaches both usage and statistical meaning.

## Phase 0: Project Foundation

Status: baseline infrastructure complete; open polish items are tracked in later
phases.

Acceptance criteria:

- The repository has a clean baseline commit.
- Tests run through `Pkg.test()`.
- Documentation builds through Documenter.jl.
- Benchmarks have an executable starting point.
- The roadmap and plot design principles are written down before major refactors.
- Generated artifacts are ignored and not committed.

## Phase 1: Plot-First User Experience

Goal: a user can pass data, uncertainties, and a model and receive a polished
fit plot without manual layout repair.

Key deliverables:

- `fitplot(...)` convenience API for common workflows.
- Robust automatic margins for data, error bars, fit curves, bands, labels, and
  statistics panels.
- Clean default theme plus `:paper`, `:latex`, and future gallery themes.
- Optional residual or pull panel.
- Clear controls for parameter/report display: plot, console, both, or none.
- Multi-fit and multi-dataset plots as first-class workflows.

Acceptance criteria:

- Default plots render correctly for small, dense, heteroscedastic, and
  x/y-uncertain datasets.
- PNG, PDF, and SVG exports are smoke-tested.
- Plot customization follows Makie conventions where possible.

## Phase 2: Statistical Core

Status: complete for the current v1 core; future statistical features belong to
new scoped phases or explicit follow-up items.

Goal: every number in `FitResult` has an explicit statistical interpretation.

Key deliverables:

- Audited convention for chi-square, negative log-likelihood, deviance, AIC,
  BIC, p-values, and degrees of freedom.
- Analytic reference tests for linear Gaussian models.
- Reference tests for Poisson, histogram, unbinned, extended-unbinned, and
  multi-fit likelihoods.
- Profile-likelihood and contour validation.
- Component-based uncertainty models with clear covariance semantics.

Acceptance criteria:

- Standard Gaussian cases match analytic solutions.
- Likelihood minima match independent reference calculations.
- Warnings are emitted when reported uncertainties or p-values are not
  statistically reliable.

Completion evidence:

- Analytic Gaussian references cover diagonal, full-covariance, constraints,
  covariance scaling, and component covariance semantics.
- Likelihood references cover Poisson, histogram Poisson, unbinned, extended
  unbinned, and mapped multi-fit workflows.
- Profile and contour scans are validated against the covariance quadratic form.
- Diagnostics warn for unavailable p-values, non-positive degrees of freedom,
  active bounds, and ill-conditioned covariance/Hessian cases.

## Phase 3: Numerics, Performance, and Torture Testing

Status: reopened as the release gate. The previous benchmark baseline was useful
but not strict enough for a public scientific package.

Goal: robustness and speed are measured continuously.

Key deliverables:

- Torture tests for hostile inputs, bad scaling, singular or invalid
  uncertainty structures, active constraints, non-finite model output, and large
  datasets.
- Clear failure semantics: invalid scientific input must fail early with a
  useful error; numerical repair must never be silent.
- Benchmark suite for least squares, full covariance, likelihood fits, profiles,
  and plotting.
- No explicit matrix inverse in production covariance calculations.
- Stable factorization strategy with Cholesky/LDL/fallback diagnostics, not
  hidden jitter.
- Clear AD/Jacobian policy: analytic, automatic differentiation, finite
  differences as fallback.
- Parameter scaling and optimizer fallback policy.

Acceptance criteria:

- Benchmarks are reproducible from the repository and compared against saved
  baselines.
- Numerically bad inputs produce actionable errors or diagnostics.
- Common workflows stay fast as features are added.
- Pathological but meaningful fits are covered before gallery examples claim
  robustness.

Current evidence:

- `benchmarks/runbenchmarks.jl` covers least-squares, no-op bounds, dense
  covariance, bounded dense covariance, Poisson likelihood, plotting, and
  profile scans.
- Static diagonal and dense Gaussian covariance terms are factorized once per
  fit evaluation cache instead of being rebuilt on every objective call.
- Production covariance calculations avoid explicit matrix inverse calls; they
  use linear solves or stable eigen fallbacks for singular/ill-conditioned
  cases.
- The fast `LsqFit` path is preserved for unbounded static least-squares fits,
  including user-supplied no-op bounds such as `[-Inf, Inf]`.
- `LsqFit` backend Jacobians are reused when constructing `FitResult`, avoiding
  an avoidable AD pass for common least-squares workflows.
- The generic `Optimization.jl` path receives cached objective data, so static
  covariance factorizations and log determinants are not recomputed by AD.
- Dense parameter-dependent `cov_x` propagation is covered by a finite-
  difference reference test for Gaussian-NLL value, gradient, and Hessian. The
  Cholesky path now keeps AD information in the actual factorization and strips
  duals only for finite-value and symmetry validation.
- X-uncertainty propagation now accepts a vectorized `x_derivative=(x, p) ->
  dy_dx` hook, so large datasets can avoid the default pointwise AD derivative
  when the model derivative is known.
- Dense covariance remains intentionally documented as an expensive `O(n^2)`
  memory and `O(n^3)` factorization path; huge correlated datasets need future
  structured covariance operators rather than dense matrices.
- Parameter covariance is intentionally treated as a local approximation.
  Nonlinear models, weak data, active bounds, and asymmetric likelihoods can
  make `Cov(p̂)` optimistic or misleading; profile and contour workflows are
  the current mitigation and must remain visible in reports and documentation.
- Profile scans can adaptively refine threshold-crossing intervals, and contour
  scans can adaptively refine cells around requested contour levels. Reference
  tests cover both adaptive paths against the quadratic covariance benchmark.
- `diagnostic_dashboard(...)` summarizes structured `diagnose(...)` findings
  into status, severity counts, and deduplicated next actions for lab workflows.

Open hardening work:

- Apply small, local numerical optimizations before larger architecture work:
  cache static parameter-constraint factorizations, avoid repeated
  model-independent work inside objectives, and add tests for every such change.
- Improve adaptive contour refinement for strongly curved/non-elliptic regions
  beyond simple level-bracketing cells.
- Add a kafe2-inspired but JuFitter-native profile/contour matrix for quick
  multi-parameter diagnosis: diagonal profile scans, lower-triangle pairwise
  contours, local covariance overlays, compact labels, and machine-readable
  diagnostics. This should explain what the scientist should conclude rather
  than merely copying kafe2's visual grammar.
- Turn the diagnostic dashboard into a visual Makie report that combines fit
  quality, pulls, profiles, contours, and next actions.
- Add structured covariance/whitening operators for large correlated data.
  Target use cases include long time series, images, spectra, detector arrays,
  and other measurements where dense covariance storage or factorization is the
  wrong asymptotic model.
- Add in-place model and residual APIs for huge datasets.
- Add explicit optimizer fallback and parameter-scaling diagnostics.
- Expand torture tests until they cover constraints, priors, likelihoods,
  profiles, contours, and multi-dataset workflows under hostile conditions.

## Phase 3.5: Diagnostic Plots and Contours

Status: blocked on the Phase 3 hardening gate.

Goal: profile, contour, residual, pull, covariance, and likelihood diagnostics
must communicate statistical meaning as clearly as kafe2/Minuit-style tools,
with Makie-level visual quality.

Key deliverables:

- Profile and contour plots with correct likelihood thresholds, labels, legends,
  and clear distinction between local covariance and profile intervals.
- A profile/contour matrix that gives a fast overview of parameter uncertainty,
  correlation, non-parabolicity, active bounds, and failed refits across several
  parameters.
- Residual and pull diagnostics that explain outliers, heteroscedasticity, and
  systematic model failures.
- External legends/statistical reports that can sit beside the plot without
  wasting plot area.
- User-extensible plot annotation API for physical extrapolations, thresholds,
  derived quantities, and uncertainty markers.
- Modular plotting helpers for adding individual objects to existing fit
  figures: extra curves, reference bands, vertical/horizontal thresholds,
  derived-quantity markers, multi-fit legends, and consistent right-side
  reports without rerunning or rewriting the fit.

Acceptance criteria:

- Diagnostic plots are backed by numerical tests for profile/contour semantics.
- The profile/contour matrix is tested against a quadratic covariance reference
  and at least one non-elliptic profile example.
- Every diagnostic example explains what a scientist should conclude from the
  plot, not merely how to call the function.
- The gallery includes difficult, messy examples before simple aesthetic
  showcase examples.
- Plot objects can be added after `plot_fit` returns without changing the
  exported figure footprint or invalidating the layout.

## Phase 4: Documentation and Gallery

Status: postponed behind Phase 3 and Phase 3.5 quality gates. Documentation
work continues only where it records architecture decisions, known limitations,
and hardening results.

Goal: documentation is good enough to teach scientific fitting, not only list
function signatures.

Key deliverables:

- Documenter.jl site with a polished landing page.
- Quickstart, tutorials, theory guide, plot gallery, and API reference.
- Realistic examples from physics, engineering, and social sciences.
- Every major example includes generated plots.
- Function-level documentation for every exported type and function.
- Installation, troubleshooting, first-compile, plotting-backend, and
  performance guidance.
- A long-term maintenance section documenting architecture, hot paths,
  bottlenecks, and known limitations.
- Light and dark documentation design with matching plot themes.
- A reproducible gallery pipeline, ideally via Literate-style examples, so code,
  text, and figures stay synchronized.
- Gallery code should look like real notebook work. Prefer explicit curated
  data arrays or small CSV inputs over long synthetic data-generation blocks on
  beginner-facing pages. Synthetic generation belongs in hidden build scripts or
  clearly labeled controlled-method demonstrations.
- Reorder the gallery as a gradual progression: quick linear calibration,
  uncertainty variants, full covariance, nonlinear resonance/oscillation,
  constraints and profiles, likelihood/histogram workflows, multi-dataset
  workflows, then advanced publication/custom-plot workflows. Photoelectric
  threshold extraction should appear after the reader has seen linear fits and
  uncertainty propagation.
- Code examples must be commented like maintainable lab notebook code:
  concise comments for data meaning, uncertainty model, statistical choice,
  derived quantities, and plot annotations.
- Statistical concepts should be introduced before they are used in examples,
  or linked to a short explanation at first use. No unexplained chi-square,
  p-value, profile likelihood, or Wilks-threshold references in beginner pages.
- Notebook-style `Real output` blocks must be verified by
  `test/docs_output_snapshots.jl`, not edited by hand. If a page shows terminal
  output, the corresponding example script needs a snapshot marker or an
  equivalent executable source of truth.

Acceptance criteria:

- The docs build locally and in CI.
- Tutorials are executable or have a documented reason when they are not.
- Statistical concepts are explained with equations and practical settings.
- Gallery examples include simple, dense-data, XY-uncertainty, full-covariance,
  profile/contour, histogram, likelihood, and multi-fit workflows.

## Phase 4.5: Python Interoperability

Status: planned, not release-claimed.

Goal: Python users can call the mature JuFitter fitting/reporting engine without
the project becoming a second Python implementation.

Key deliverables:

- Validate the current package through JuliaCall/PythonCall from a clean Python
  environment.
- Provide a minimal Python example that fits arrays, reads parameter estimates,
  and prints `report_text(...)` without loading Makie.
- Document array conversion and ownership rules for NumPy inputs and Julia
  outputs.
- Keep plotting optional: Python interoperability must work for fitting and
  reports even when CairoMakie is not installed.
- Add a CI or release-gate job that runs the Python interoperability example
  when PythonCall/JuliaCall support is enabled.

Acceptance criteria:

- A documented `juliacall` example runs from Python on a clean checkout.
- Returned fit parameters, uncertainties, diagnostics, and text reports are
  accessible from Python.
- Limitations are explicit: Julia startup cost, package environment setup, and
  plotting-backend requirements are not hidden.
- The first public release has a clear path from `Pkg.add` to first successful
  plot and report.
- The docs explain when not to trust local errors, p-values, or chi-square/ndf.

## Local Browser QA

Browser-based documentation QA needs explicit permission for the local docs URL
used by the in-app browser, currently `http://localhost:8010`. Without that
permission, automated checks can still build the docs and inspect generated
files, but they cannot verify visible sidebar/plot geometry in the browser.
When browser access is available, each plot-layout change should measure the
visible plot, article column, sidebar, viewport width, and horizontal overflow
on at least one gallery page.

## Agent Operating Model

Agents are used only for scoped, non-overlapping work:

- Plot work owns `src/plotting.jl`, plot tests, and gallery examples.
- Statistical work owns cost semantics, reference tests, and theory docs.
- Performance work owns `benchmarks/` and isolated hot spots.
- Documentation work owns tutorials, page structure, and examples.

The main integration thread keeps architecture decisions and final review.
