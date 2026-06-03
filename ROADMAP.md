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
- Dense covariance remains intentionally documented as an expensive `O(n^2)`
  memory and `O(n^3)` factorization path; huge correlated datasets need future
  structured covariance operators rather than dense matrices.
- Profile scans can adaptively refine threshold-crossing intervals, and contour
  scans can adaptively refine cells around requested contour levels. Reference
  tests cover both adaptive paths against the quadratic covariance benchmark.
- `diagnostic_dashboard(...)` summarizes structured `diagnose(...)` findings
  into status, severity counts, and deduplicated next actions for lab workflows.

Open hardening work:

- Improve adaptive contour refinement for strongly curved/non-elliptic regions
  beyond simple level-bracketing cells.
- Turn the diagnostic dashboard into a visual Makie report that combines fit
  quality, pulls, profiles, contours, and next actions.
- Add structured covariance/whitening operators for large correlated data.
- Add in-place model and residual APIs for huge datasets.
- Add explicit optimizer fallback and parameter-scaling diagnostics.
- Add performance-budget CI for representative hot paths.
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
- Residual and pull diagnostics that explain outliers, heteroscedasticity, and
  systematic model failures.
- External legends/statistical reports that can sit beside the plot without
  wasting plot area.
- User-extensible plot annotation API for physical extrapolations, thresholds,
  derived quantities, and uncertainty markers.

Acceptance criteria:

- Diagnostic plots are backed by numerical tests for profile/contour semantics.
- Every diagnostic example explains what a scientist should conclude from the
  plot, not merely how to call the function.
- The gallery includes difficult, messy examples before simple aesthetic
  showcase examples.

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

Acceptance criteria:

- The docs build locally and in CI.
- Tutorials are executable or have a documented reason when they are not.
- Statistical concepts are explained with equations and practical settings.
- Gallery examples include simple, dense-data, XY-uncertainty, full-covariance,
  profile/contour, histogram, likelihood, and multi-fit workflows.
- The first public release has a clear path from `Pkg.add` to first successful
  plot and report.
- The docs explain when not to trust local errors, p-values, or chi-square/ndf.

## Agent Operating Model

Agents are used only for scoped, non-overlapping work:

- Plot work owns `src/plotting.jl`, plot tests, and gallery examples.
- Statistical work owns cost semantics, reference tests, and theory docs.
- Performance work owns `benchmarks/` and isolated hot spots.
- Documentation work owns tutorials, page structure, and examples.

The main integration thread keeps architecture decisions and final review.
