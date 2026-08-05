# Installation

JuFitter declares support for Julia 1.10 and later. The CI workflow is
configured to run the core and full-package gates on Julia 1.10 and Julia 1.12;
the complete current matrix must pass before public release. Julia 1.10 is the
compatibility floor.

!!! note "Pre-release package"
    JuFitter is not registered yet and the repository must not be published
    before maintainer review. Until that release decision, install from a local
    checkout. The registry command below documents the intended final path, not
    a command that works today.

## Current Local Checkout

From a terminal, instantiate the numerical core in the repository root:

```bash
cd /path/to/JuFitter
julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate()'
```

Start a Julia session in that environment:

```bash
julia --project=.
```

Then load the package:

```julia
using JuFitter
```

This loads fitting, likelihoods, diagnostics, profiles, contours, and text
reports. It does **not** load Makie.

The checked-out repository keeps plotting and documentation dependencies in a
separate environment. Instantiate it when running the gallery or building the
site:

```bash
julia --project=docs --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=docs examples/gallery/01_quickstart_linear.jl
```

The example writes its figure to the ignored `examples/output/` directory.

## Registry Installation After Release

Once JuFitter is registered, a normal Julia environment will use:

```julia
using Pkg
Pkg.add("JuFitter")
```

For static PNG, PDF, and SVG plots, install CairoMakie in the same environment:

```julia
using Pkg
Pkg.add(["JuFitter", "CairoMakie"])

using JuFitter
using CairoMakie
```

Keeping CairoMakie optional is deliberate. A batch analysis can fit data,
create reports, and run diagnostics without compiling a graphics stack.

## First Use And Compilation

The first `using JuFitter` in a new environment compiles the numerical core.
The first `using CairoMakie` and first rendered figure take longer because Julia
also compiles Makie's layout, text, and rendering methods. Later sessions reuse
the precompile cache unless Julia, package versions, preferences, or the target
environment change.

Do not use the full package test suite to check an installation; it is a slow
release gate. A core-only check is enough:

```bash
julia --project=. --startup-file=no -e 'using JuFitter; println("JuFitter core ready")'
```

For plotting, run the tracked quickstart example shown above. It exercises the
same API used by the first tutorial and confirms CairoMakie export.

## Python Interoperability (Experimental)

Python can call the Julia implementation through JuliaCall. This path reuses
JuFitter's fit results, reports, and diagnostics; it is not a separate Python
rewrite.

```bash
python3 -m pip install juliacall
python3 examples/python/fit_from_python.py
```

The tracked script develops the local checkout into JuliaCall's managed Julia
environment and keeps plotting out of the process. The release gate is opt-in:

```bash
JUFITTER_RUN_PYTHON_INTEROP=1 julia --project=. --startup-file=no test/python_interop_gate.jl
```

Python support remains experimental or deferred for public v0 claims until the
same path is observed on the selected release CI or release machine.

## Troubleshooting

| symptom | first check |
| --- | --- |
| `using JuFitter` is slow once | Let precompilation finish; this is not fit runtime. |
| Every fresh session recompiles | Reuse the same project and depot; check whether Julia or package versions keep changing. |
| `plot_fit` says the extension is unavailable | Add and load `CairoMakie` before calling plotting functions. |
| PDF or SVG export fails | Verify a minimal CairoMakie figure in the same environment; inspect backend and font errors first. |
| A fit is unexpectedly slow | Check for dense covariance, bounds, constraints, priors, parameter-dependent covariance, or pointwise x-derivatives. These select more general numerical paths. |
| Package versions will not resolve | Confirm Julia is at least 1.10 and instantiate a clean environment rather than mixing incompatible manifests. |

Continue with the [Quickstart](quickstart.md). For package internals and scaling
limits, see [Backend Design](backend_design.md) and [Performance](performance.md).
