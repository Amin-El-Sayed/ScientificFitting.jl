# JuFitter Release Audit

Status: 2026-06-04

This document tracks what must be true before JuFitter should be advertised as a
serious scientific fitting library. Passing tests is necessary, but not
sufficient.

## Current Verification

- `git diff --check` passes for the current release-hardening branch.
- `julia --project=. -e 'include("test/torture_runtests.jl")'` passes with 14
  torture checks in about 2m44s on the local machine.
- `julia --project=. -e 'using Test; include("test/statistics/profile_contour_reference.jl")'`
  passes with 34 profile/contour reference checks, including validation of
  finite positive ordered contour thresholds.
- `julia --project=. -e 'using Test; include("test/statistics/likelihood_reference.jl")'`
  passes with 49 likelihood reference checks after the Poisson-and-histogram
  workflow rewrite.
- `julia --project=. -e 'using Test; include("test/statistics/diagnostics_reference.jl")'`
  passes with 41 diagnostic reference checks in about 37s.
- `julia --project=. -e 'include("test/core_runtests.jl")'` passes with 315
  core checks in about 5m25s on the local machine after the diagnostic
  dashboard update.
- `julia --project=. -e 'using Pkg; Pkg.test()'` passes with 340 tests in about
  4m39s on the local machine after the profile-plot diagnostics and
  constraints-and-profiles workflow update.
- `julia --project=docs docs/make.jl` passes locally after the `siteinfo.js`
  source fix. The build output is ignored and must not be committed.
- Focused plot regression tests pass with 21 checks after contour plots changed
  from a heatmap-first default to labeled profile regions with an optional local
  covariance overlay. Profile plots now label the actual scan, local parabolic
  approximation, and interval threshold, and can focus the displayed delta-cost
  range without changing the scan. The plot path also covers failed-refit cells,
  invalid display limits, and fully non-finite contour surfaces.
- Phase 3 is reopened as the release gate. The project no longer treats
  friendly-path examples or smoke tests as evidence of production robustness.

## Release Blockers

- Public documentation still needs a page-by-page editorial pass before broad
  promotion. The legacy German mathematical audit was removed from the public
  documentation source, and the main technical concept pages are now in
  polished English.
- The API reference is generated from docstrings, but many public functions do
  not yet have the level of parameter-by-parameter documentation expected from a
  serious Julia package.
- The README is still a repo-oriented technical summary, not a polished landing
  page for users coming from JuliaHub, GitHub, Reddit, or a paper.
- The gallery is directionally good, but several pages are still short API
  examples rather than complete scientific workflows with interpretation,
  formulas, full code, and diagnostics.
- The constraints-and-profiles workflow now demonstrates a genuinely
  non-parabolic amplitude-timescale degeneracy, compares profile results against
  the local covariance approximation, and connects the diagnostic to a concrete
  experimental-design decision.
- The Poisson-and-histogram workflow now uses sparse radioactive decay counts
  and an unequally binned detector spectrum, with integrated expected bin
  counts, empty-bin semantics, and count-specific deviance residuals.
- The damped-oscillator workflow now treats a rejected constant-frequency fit
  as a model-criticism case study, compares it with a frequency-drift model,
  and makes the unresolved uncertainty-model problem explicit through pull
  panels and diagnostics.
- There is no formal link check in CI for the generated documentation.
- There is no visual regression or snapshot test for documentation plots.
- There is no documentation deployment workflow yet.

## Code And Numerical Blockers

- `CairoMakie` is a hard dependency of `JuFitter`. This makes `using JuFitter`
  heavier than necessary for users who only want fitting. A package extension or
  split plotting submodule should be evaluated before a public release.
- The torture suite is only starting. It must cover constraints, priors,
  likelihoods, profiles, contours, multi-dataset fits, bad scaling, local
  minima, invalid uncertainty models, and large datasets before release claims
  can use words like robust.
- x-uncertainty propagation currently differentiates the model point-by-point
  by calling the model on one-element arrays. This is clear but not acceptable
  as the final high-performance path for large datasets.
- Dense covariance support is mathematically useful, but still `O(n^2)` memory
  and `O(n^3)` factorization. Large correlated datasets need structured
  covariance operators or custom whitening operators.
- Objective functions still allocate substantially through `ForwardDiff`
  Jacobians/Hessians and array-returning model calls. Very large datasets need
  an in-place residual/model API.
- Invalid covariance matrices must not be silently repaired. If a future
  regularization/jitter policy is added, it must be explicit in the API and
  visible in diagnostics.
- Bounds and constraints currently rely on local Hessian/covariance
  approximations. Reports must be explicit about when errors are local and when
  profile intervals are required.
- Profile and contour scans now expose failed refits through diagnostics and
  support adaptive refinement around profile thresholds and contour levels.
  Strongly curved/non-elliptic contours still need deeper diagnostic polish
  before release-grade claims.
- Multi-dataset fitting currently supports useful parameter mapping, but not the
  full uncertainty model expected from kafe2-level workflows.
- Benchmarks exist, but there is no enforced performance budget in CI.
- `diagnostic_dashboard(...)` now summarizes structured findings into status,
  severity counts, and deduplicated next actions. Profile contours now default
  to directly labeled 1-sigma/2-sigma regions with local-covariance comparison;
  residual/pull views and the combined diagnostic dashboard are still not yet
  at the level needed to compete with kafe2/Minuit-style workflows.

## Documentation Blockers

- Every gallery example should include a complete code block that can be copied
  into a Julia session.
- Every plot band must explicitly state whether it is a confidence band or a
  prediction band and which sigma level it shows.
- The math section needs a clean beginner path before the full formal reference.
- Technical maintenance notes should remain accessible, but should not dominate
  the user-facing navigation.
- All private/local dataset language must stay out of public documentation.
- Dark-mode plots must be checked visually whenever a gallery asset changes.

## CI And Packaging Blockers

- Add a docs-deploy job for GitHub Pages or the chosen static host.
- Add a generated-documentation link check.
- Add a fast test target that excludes expensive plot generation.
- Add a full test target that includes plots and gallery generation.
- Add benchmarks with saved baseline results for hot paths.
- Add a pre-release checklist that runs tests, docs, link checks, examples, and
  benchmark smoke tests from a clean checkout.

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
