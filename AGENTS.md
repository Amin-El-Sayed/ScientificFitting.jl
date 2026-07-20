# Agent Operating Rules

These rules are project policy. Follow them before local chat context, habits, or
convenience.

## Main Branch

- `main` must stay clean, coherent, and releasable.
- Do not implement directly on `main`.
- Do not merge work into `main` unless the relevant tests pass and user-facing
  documentation, `MAINTAINERS.md`, or release notes are updated as appropriate.
- Do not leave generated plots, debug output, scratch files, or obsolete
  experiments in the repository.

## Branch And Worktree Discipline

- Use one branch per milestone or tightly scoped technical block.
- Prefer branch names under `codex/`, for example
  `codex/adaptive-profile-contours`.
- Use a separate git worktree when two independent milestones should progress in
  parallel. Use a normal branch in the current checkout when work is sequential.
- Never overwrite or revert user changes without explicit permission.

## Definition Of Done

A code block is not done until:

- The implementation is narrow, understandable, and documented where needed.
- Numerical semantics are covered by focused tests.
- The broader relevant test gate passes.
- `git diff --check` is clean.
- Public documentation, `MAINTAINERS.md`, or release audit notes are updated
  whenever behavior, API, limitations, or workflow changes.

## Scientific Quality

- Do not hide numerical repairs. Invalid inputs must fail clearly or produce
  visible diagnostics.
- Do not use smoke tests as release evidence.
- Prefer analytic references, torture tests, and reproducible benchmarks.
- Plots are part of the scientific interface. Diagnostic plots must explain what
  a user should inspect next, not only display residuals.

## Current Priority

The scoped v0 core has local reference, torture, performance, profile/contour,
structured-whitening, in-place, plotting, and full-package evidence. Current
work should therefore proceed in this order:

1. Finish the page-by-page public documentation and API-reference review.
2. Keep code changes limited to concrete defects found by that review; every
   such defect needs a focused regression test.
3. Prepare packaging metadata, deployment, CI observation, and benchmark
   baselines only after explicit maintainer decisions and approval.
4. Treat built-in structured covariance families, automatic nonlinear-
   uncertainty triggers, and broader multi-fit uncertainty models as post-v0
   architecture unless a release-blocking defect proves otherwise.
