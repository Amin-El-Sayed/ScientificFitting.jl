# Installation

ScientificFitting supports Julia 1.10 and later. CI targets Julia 1.10 and the
latest stable Julia for the core and full-package tests. The optional
NativeMinuit adapter requires Julia 1.11+; its integration tests also target
that minimum version. Documentation rendering stays pinned to Julia 1.12.

For Python, see [Python Interface](python.md).

!!! note "Early release"
    Version 0.2 is an early work-in-progress release. Please report unexpected
    behavior with a minimal reproducer and the complete diagnostic output.

## Install The Julia Package

Install the registered package from Julia's General registry:

```julia
using Pkg
Pkg.add("ScientificFitting")
```

Then load the numerical core:

```julia
using ScientificFitting
```

This loads fitting, likelihoods, diagnostics, profiles, contours, and text
reports. It does **not** load Makie.

## Work From A Checkout

From a terminal, instantiate the numerical core in the repository root:

```bash
cd /path/to/ScientificFitting
julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate()'
```

Start a Julia session in that environment:

```bash
julia --project=.
```

The checked-out repository keeps plotting and documentation dependencies in a
separate environment. Instantiate it when running the gallery or building the
site:

```bash
julia --project=docs --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=docs examples/gallery/01_quickstart_linear.jl
```

The example writes its figure to the ignored `examples/output/` directory.

For static PNG, PDF, and SVG plots, install CairoMakie in the same environment:

```julia
using Pkg
Pkg.add(["ScientificFitting", "CairoMakie"])

using ScientificFitting
using CairoMakie
```

Keeping CairoMakie optional is deliberate. A batch analysis can fit data,
create reports, and run diagnostics without compiling a graphics stack.

## First Use And Compilation

The first `using ScientificFitting` in a new environment compiles the numerical core.
The first `using CairoMakie` and first rendered figure take longer because Julia
also compiles Makie's layout, text, and rendering methods. Later sessions reuse
the precompile cache unless Julia, package versions, preferences, or the target
environment change.

Do not use the full package test suite to check an installation; it is a slow
release gate. A core-only check is enough:

```bash
julia --project=. --startup-file=no -e 'using ScientificFitting; println("ScientificFitting core ready")'
```

For plotting, run the tracked quickstart example shown above. It exercises the
same API used by the first tutorial and confirms CairoMakie export.

## Python Interface

In a Python 3.10+ virtual environment:

```bash
python -m pip install 'scientificfitting[plot]'
```

The Python interface uses NumPy models and optional native Matplotlib plots,
without Makie. JuliaCall provisions the registered Julia core automatically on
first use, which needs network access and compilation. Omit `[plot]` when no
plots are needed. Supported features, examples, and the conda-forge submission
status are on the [Python Interface](python.md) page.

## Troubleshooting

| symptom | first check |
| --- | --- |
| `using ScientificFitting` is slow once | Let precompilation finish; this is not fit runtime. |
| Every fresh session recompiles | Reuse the same project and depot; check whether Julia or package versions keep changing. |
| `plot_fit` says the extension is unavailable | Add and load `CairoMakie` before calling plotting functions. |
| PDF or SVG export fails | Verify a minimal CairoMakie figure in the same environment; inspect backend and font errors first. |
| A fit is unexpectedly slow | Check for dense covariance, bounds, constraints, priors, parameter-dependent covariance, or pointwise x-derivatives. These select more general numerical paths. |
| Package versions will not resolve | Confirm Julia is at least 1.10 and instantiate a clean environment rather than mixing incompatible manifests. |

Continue with the [Quickstart](quickstart.md). For package internals and scaling
limits, see [Backend Design](backend_design.md) and [Performance](performance.md).
