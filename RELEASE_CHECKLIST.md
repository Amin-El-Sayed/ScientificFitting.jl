# JuFitter Pre-Release Checklist

This checklist is the local release gate before JuFitter is advertised,
registered, deployed, pushed to a public repository, or announced. It is
deliberately stricter than the normal edit-test loop.

Publication policy: do not push, publish, register, deploy documentation, or
make the repository public without explicit manual approval from
Amin_El_Sayed.

## 1. Repository State

Start from a clean release branch or reviewed release candidate:

```bash
git status --short --branch
git diff --check
```

Required manual checks:

- The branch is not `main` unless this is the final reviewed merge point.
- The worktree contains no generated build output, debug plots, scratch files,
  or machine-local benchmark artifacts.
- `RELEASE_AUDIT.md`, `ROADMAP.md`, and public documentation describe known
  limitations honestly.

## 2. Core Correctness Gates

Run the core and package gates from a clean checkout:

```bash
julia --project=. --startup-file=no -e 'include("test/core_runtests.jl")'
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
```

Run focused scientific reference gates when touching the statistical core:

```bash
julia --project=. --startup-file=no -e 'using Test; include("test/statistics/profile_contour_reference.jl")'
julia --project=. --startup-file=no -e 'using Test; include("test/statistics/likelihood_reference.jl")'
julia --project=. --startup-file=no -e 'using Test; include("test/statistics/covariance_semantics_reference.jl")'
julia --project=. --startup-file=no -e 'using Test; include("test/statistics/diagnostics_reference.jl")'
julia --project=. --startup-file=no -e 'include("test/torture_runtests.jl")'
```

## 3. Documentation And Gallery Gates

Run the source documentation gates:

```bash
julia --project=. --startup-file=no test/docs_gallery_gate.jl
julia --project=. --startup-file=no test/docs_public_release_gate.jl
julia --project=. --startup-file=no test/docs_api_reference_gate.jl
julia --project=. --startup-file=no test/docs_link_gate.jl
julia --project=. --startup-file=no test/docs_visual_asset_gate.jl
julia --project=. --startup-file=no test/docs_visual_snapshot_gate.jl
```

Build the static site and then validate rendered output:

```bash
julia --project=docs --startup-file=no docs/make.jl
julia --project=. --startup-file=no test/docs_html_link_gate.jl
julia --project=. --startup-file=no test/docs_output_snapshots.jl
```

Run plot regression tests in the docs environment:

```bash
julia --project=docs --startup-file=no test/plots/fitplot.jl
```

Required manual checks:

- Review changed gallery PNGs visually before updating
  `test/docs_visual_snapshot_manifest.txt`.
- Check light and dark documentation modes.
- Remove `docs/build/` after local verification unless intentionally packaging
  the static site.

## 4. Performance And Benchmark Evidence

Run the startup probe to verify that the fitting/reporting core still starts
without loading Makie:

```bash
julia --project=. --startup-file=no benchmarks/startup_probe.jl --save=/tmp/jufitter-startup-probe.toml
julia --project=. --startup-file=no test/startup_probe_gate.jl
```

Run the benchmark contract gate:

```bash
julia --project=. --startup-file=no test/benchmark_contract_gate.jl
```

Run the performance regression guard:

```bash
julia --project=. --startup-file=no test/performance_budget_gate.jl
```

For public performance claims, use saved benchmark evidence instead of the
budget gate:

```bash
julia --project=benchmarks benchmarks/runbenchmarks.jl --save=benchmarks/output/local-baseline.toml
julia --project=benchmarks benchmarks/runbenchmarks.jl --compare=benchmarks/output/local-baseline.toml
```

Do not advertise speed claims until the reference hardware or CI runner is
named and the benchmark output is reviewed.

Do not use `--allow-metadata-mismatch` for release evidence. That flag is only
for exploratory comparisons when you knowingly compare different machines,
Julia versions, or thread configurations.

Required scientific limitation checks:

- Dense covariance matrices are the exact small/medium-data path, not a large
  data strategy. Before claiming support for long time series, images, spectra,
  or detector arrays, either provide structured covariance or custom whitening
  evidence, or state the current `O(n^2)` memory and `O(n^3)` factorization
  limitation clearly.
- Parameter covariance is a local approximation. For nonlinear models, weak
  data, active bounds, or asymmetric likelihoods, verify that documentation and
  diagnostics point users to profile/contour intervals instead of implying that
  symmetric covariance errors are final.
- Parameter-dependent dense `cov_x` propagation is an audit-sensitive AD path.
  Keep finite-difference value, gradient, and Hessian references in the release
  gate whenever validation, conversion, or factorization code touches dense
  effective covariance logic.

## 5. Python Interoperability

Python use is intended through JuliaCall/PythonCall, but it is not a default
Julia test dependency. If Python support is claimed for the release, run:

```bash
JUFITTER_RUN_PYTHON_INTEROP=1 julia --project=. --startup-file=no test/python_interop_gate.jl
```

If that gate is not run in a clean `juliacall` environment, document Python
support as experimental or deferred for the release.

## 6. CI, Deployment, And Publication

Before public promotion:

- GitHub Actions must pass on the target branch for core, package, docs, plots,
  and documentation gates.
- Documentation deployment must be reviewed separately; do not add or trigger
  public deployment without explicit approval.
- Confirm package metadata, version, license, README, release notes, and known
  limitations.
- Confirm Git identity before any public push:
  `Amin-El-Sayed <78275938+Amin-El-Sayed@users.noreply.github.com>`.
- Do not register the package, publish docs, or announce on Reddit/Discourse
  until Amin_El_Sayed has manually approved the exact release candidate.
