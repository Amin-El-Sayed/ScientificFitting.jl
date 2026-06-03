# Agent Operating Rules

These rules are project policy. Follow them before local chat context, habits, or
convenience.

## Main Branch

- `main` must stay clean, coherent, and releasable.
- Do not implement directly on `main`.
- Do not merge work into `main` unless the relevant tests pass and user-facing
  documentation or maintenance notes are updated.
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
- Documentation, maintenance notes, or release audit notes are updated whenever
  behavior, API, limitations, or workflow changes.

## Scientific Quality

- Do not hide numerical repairs. Invalid inputs must fail clearly or produce
  visible diagnostics.
- Do not use smoke tests as release evidence.
- Prefer analytic references, torture tests, and reproducible benchmarks.
- Plots are part of the scientific interface. Diagnostic plots must explain what
  a user should inspect next, not only display residuals.

## Current Priority

Code robustness comes before documentation polish. The current technical backlog
is:

1. Robust and adaptive profile/contour scans.
2. Useful diagnosis dashboard and diagnostic plots.
3. Performance gates and benchmark baselines.
4. Structured covariance or custom whitening operators.
5. In-place model/residual APIs for very large datasets.
6. Full documentation pass only after the core behavior is stable.
