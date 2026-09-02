# Development Workflow

ScientificFitting is developed as a scientific package. The repository workflow is built
around keeping `main` clean and making every merged block testable,
understandable, and documented.

## Branches vs Worktrees

A **branch** is a named line of history. It can be checked out in the current
directory:

```bash
git switch -c codex/adaptive-profile-contours
```

A **worktree** is an additional checkout of the same repository. It lets two
branches exist in two directories at the same time:

```bash
git worktree add ../ScientificFitting-adaptive codex/adaptive-profile-contours
```

Use branches for sequential work. Use worktrees when two independent milestones
must progress in parallel without mixing files.

## Standard Milestone Flow

1. Start from clean `main`.
2. Create a scoped branch under `codex/`.
3. Implement one milestone or one coherent sub-block.
4. Run focused tests while iterating.
5. Run the broader gate before commit.
6. Update public docs, `MAINTAINERS.md`, or release notes as appropriate.
7. Commit with a clear message.
8. Merge to `main` only after the gate is green.

## Codex Operating Model

Use Codex as an implementation worker, not as an uncontrolled auto-merge bot.
The safe loop is:

1. Keep the current thread responsible for integration decisions.
2. Let Codex work on one scoped `codex/` branch or worktree at a time.
3. Ask for a concrete outcome, for example "finish adaptive contours" or
   "harden dense covariance errors".
4. Require tests, relevant documentation, and a clean diff before a commit.
5. Merge to `main` only from a reviewed, green branch.

Heartbeats and automations are useful for continuing the loop when nobody is
actively typing. They must be state-driven: each run reads the repository state,
the roadmap, and the audit notes before choosing the next task. They must not
blindly repeat an old checklist, and they must not merge to `main`.

Skills are local instruction bundles for specialized work. Use them when the
task matches the skill, for example browser verification for the documentation
site or spreadsheet tooling for CSV/XLSX data extraction. Skills do not replace
the repository rules in `AGENTS.md`; they only add task-specific procedure.

## Test Gates

Fast focused gates:

```bash
julia --project=. -e 'include("test/torture_runtests.jl")'
julia --project=. -e 'using Test; include("test/statistics/profile_contour_reference.jl")'
```

Core gate:

```bash
julia --project=. -e 'include("test/core_runtests.jl")'
```

Full gate:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Documentation gate:

```bash
julia --project=docs docs/make.jl
```

Always run:

```bash
git diff --check
```

## Merge Rule

Do not merge to `main` unless:

- The relevant focused tests pass.
- The core or full gate passes, depending on risk.
- Public documentation, `MAINTAINERS.md`, or release notes are updated when the
  corresponding contract changes.
- Generated artifacts and debug outputs are absent.
- The diff is reviewed for accidental unrelated changes.

## Generated Files

The following are local artifacts and should not be committed:

- `docs/build/`
- `examples/output/`
- `benchmarks/output/`
- `.DS_Store`
- temporary debug plots

If a generated figure is needed for documentation, place it intentionally under
`docs/src/assets/` and regenerate it from a tracked script.
