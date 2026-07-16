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

Status: locally complete for the scoped v0 core. Static covariance and
parameter-constraint factors are cached, no-op bounds preserve the fast path,
incompatible explicit backend requests fail rather than weakening the cost,
hostile-input and local-curvature diagnostics are release-tested, and benchmark
coverage includes least squares, dense/sparse/structured covariance,
likelihoods, profiles, diagnostics, and plotting. Remote CI and a named release
benchmark baseline remain publication gates.

Dense covariance fits remain the intentionally expensive path: they are correct
and tested, while large static correlated datasets can use a validated
matrix-free `WhiteningOperator` instead of materialized dense matrices.

## Phase 3.5

Diagnostic plots and contours: profile/contour matrices, local-covariance
comparisons, residual/pull diagnostics, triage, and post-fit Makie annotations.

Status: complete for the scoped v0 interface. Richer combined dashboards and
automatic nonlinear-likelihood triggers are post-v0 candidates.

## Phase 4

Documentation and gallery: tutorials, theory guide, API reference, and
real-data examples with generated plots.

Status: active release focus. The docs now have explicit pages for installation,
quickstart, gallery workflows, mathematics/statistics, plotting design,
performance, maintenance notes, and curated API reference. The remaining work
is the final page-by-page editorial and visual review before public promotion;
generated output blocks and plot assets are already guarded by release tests.
