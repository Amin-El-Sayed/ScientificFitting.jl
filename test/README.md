# JuFitter Test Layout

The test suite is split by purpose, not by implementation file.

- `regression/`: current API and behavior coverage kept stable while the package is refactored.
- `statistics/`: analytic reference cases, likelihood conventions, coverage checks, and covariance semantics.
- `numerics/`: conditioning, fallback paths, optimizer behavior, and derivative accuracy.
- `api/`: user-facing constructor, keyword, validation, and error-message behavior.
- `plots/`: plot generation, layout robustness, export formats, and visual-regression checks.
- `docs_gallery_gate.jl`: fast structural gate for public gallery pages and
  referenced plot assets.
- `docs_public_release_gate.jl`: fast hygiene gate for public Documenter pages
  and README. It rejects AI/placeholder markers, private paths, and
  course-internal wording before broad promotion. It also blocks known stale
  public API identifiers that previously appeared in prose.
- `docs_api_reference_gate.jl`: fast API-reference gate that requires every
  exported public binding to have a REPL/Documenter-visible docstring.
- `docs_link_gate.jl`: fast local-link gate for Markdown links, HTML links, and
  image sources under `docs/src`.
- `docs_html_link_gate.jl`: rendered-site link gate for `docs/build` after
  `docs/make.jl` has run.
- `docs_visual_asset_gate.jl`: visual-asset sanity gate for documentation PNGs:
  valid PNG headers, minimum dimensions, consistent style-pair dimensions, and
  no unreferenced gallery PNG leftovers.
- `docs_output_snapshots.jl`: heavy release gate that executes the documented
  gallery/quickstart examples with snapshot markers and verifies that every
  documented `Real output` block is an ordered subset of real script output.
  It sets `JUFITTER_DOC_SNAPSHOT_ONLY=1`, so example scripts compute the same
  fit/report values but skip Makie asset rendering.
- `python_interop_gate.jl`: opt-in gate for calling JuFitter from Python through
  JuliaCall. The default run passes with an informational note; setting
  `JUFITTER_RUN_PYTHON_INTEROP=1` requires `python3` and `juliacall` and runs
  `examples/python/fit_from_python.py`.
- `performance_budget_gate.jl`: steady-state performance gate for representative
  hot paths. It warms compilation first and uses deliberately broad budgets to
  catch large regressions rather than benchmark-machine noise.

New tests should be narrow and deterministic. Expensive Monte Carlo or benchmark
checks belong in `benchmarks/` or a dedicated validation job, not in the default
unit-test path.

Run the current torture checks for the active hardening loop:

```julia
include("test/torture_runtests.jl")
```

Run the broader core checks without plot tests:

```julia
include("test/core_runtests.jl")
```

Run only the statistical reference checks with:

```julia
include("test/statistics/runtests.jl")
```

Run the documentation release gates with:

```julia
include("test/docs_gallery_gate.jl")
include("test/docs_public_release_gate.jl")
include("test/docs_api_reference_gate.jl")
include("test/docs_link_gate.jl")
include("test/docs_visual_asset_gate.jl")
# Run after `julia --project=docs docs/make.jl`.
include("test/docs_html_link_gate.jl")
include("test/docs_output_snapshots.jl")
```

Run the Python interoperability gate after installing JuliaCall with:

```bash
JUFITTER_RUN_PYTHON_INTEROP=1 julia --project=. --startup-file=no test/python_interop_gate.jl
```

If `juliacall` is not installed, leave the environment variable unset; the gate
then records that the external prerequisite is intentionally absent without
making default Julia tests depend on Python packaging state.

Run the performance budget gate with:

```bash
julia --project=. --startup-file=no test/performance_budget_gate.jl
```

On intentionally slow CI runners, set `JUFITTER_PERFORMANCE_BUDGET_SCALE` to a
larger value. Do not use this gate for marketing numbers; use
`benchmarks/runbenchmarks.jl` for measured benchmark reports.
