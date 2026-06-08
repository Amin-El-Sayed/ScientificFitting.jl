# Roadmap

The canonical development roadmap lives in the repository root as `ROADMAP.md`.
This documentation page mirrors its phase structure for the hosted docs.

## Phase 0

Project foundation: clean baseline, CI, test structure, benchmarks, and
Documenter.jl documentation.

## Phase 1

Plot-first user experience: beautiful defaults, robust layout, publication
exports, and Makie-compatible customization.

## Phase 2

Statistical core: audited likelihood conventions, covariance semantics,
constraints, profiles, contours, and reference tests.

Status: complete for the current v1 core. The statistical reference suite now
covers Gaussian, likelihood, covariance-component, profile/contour,
profile-matrix, and diagnostic-warning semantics.

## Phase 3

Numerics and performance: reproducible benchmarks, stable factorizations,
diagnostics, AD/Jacobian policy, and optimizer fallback strategy.

Status: reopened as the release gate. Static covariance factorizations are
cached, explicit production matrix inverses were removed, no-op bounds keep the
fast least-squares backend, and benchmark coverage now includes least-squares,
dense covariance, bounded covariance, Poisson likelihood, plotting, and
profiles. That evidence is useful, but it is not enough for public robustness
claims until the torture suite covers more hostile constraints, priors,
likelihood, profile, contour, and multi-dataset workflows.

Dense covariance fits remain the intentionally expensive path: they are correct
and tested, but very large correlated datasets need future structured
covariance operators instead of materialized dense matrices.

## Phase 4

Documentation and gallery: tutorials, theory guide, API reference, and
real-data examples with generated plots.

Status: in progress. The docs now have explicit pages for installation,
quickstart, gallery workflows, mathematics/statistics, plotting design,
performance, maintenance notes, and curated API reference. The remaining work
is the final page-by-page editorial and visual review before public promotion;
generated output blocks and plot assets are already guarded by release tests.
