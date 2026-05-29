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

## Phase 3: Numerics and Performance

Goal: robustness and speed are measured continuously.

Key deliverables:

- Benchmark suite for least squares, full covariance, likelihood fits, profiles,
  and plotting.
- No explicit matrix inverse in production covariance calculations.
- Stable factorization strategy with Cholesky/LDL/fallback diagnostics.
- Clear AD/Jacobian policy: analytic, automatic differentiation, finite
  differences as fallback.
- Parameter scaling and optimizer fallback policy.

Acceptance criteria:

- Benchmarks are reproducible from the repository.
- Numerically bad inputs produce actionable diagnostics.
- Common workflows stay fast as features are added.

## Phase 4: Documentation and Gallery

Goal: documentation is good enough to teach scientific fitting, not only list
function signatures.

Key deliverables:

- Documenter.jl site with a polished landing page.
- Quickstart, tutorials, theory guide, plot gallery, and API reference.
- Realistic examples from physics, engineering, and social sciences.
- Every major example includes generated plots.

Acceptance criteria:

- The docs build locally and in CI.
- Tutorials are executable or have a documented reason when they are not.
- Statistical concepts are explained with equations and practical settings.

## Agent Operating Model

Agents are used only for scoped, non-overlapping work:

- Plot work owns `src/plotting.jl`, plot tests, and gallery examples.
- Statistical work owns cost semantics, reference tests, and theory docs.
- Performance work owns `benchmarks/` and isolated hot spots.
- Documentation work owns tutorials, page structure, and examples.

The main integration thread keeps architecture decisions and final review.
