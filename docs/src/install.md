# Install

JuFitter is a Julia package. During active development, use it from a checked
out repository:

```julia
using Pkg
Pkg.activate("/path/to/JuFitter")
Pkg.instantiate()
```

From the repository root:

```bash
julia --project=.
```

Then:

```julia
using JuFitter
```

This loads the fitting, likelihood, diagnostics, profiling, and reporting core.
It does not load Makie.

When JuFitter is registered, the intended user-facing installation path is:

```julia
using Pkg
Pkg.add("JuFitter")
```

For plotting, install a Makie backend in the same environment and load it before
calling JuFitter plot functions:

```julia
using Pkg
Pkg.add(["JuFitter", "CairoMakie"])

using JuFitter
using CairoMakie
```

## First Compile

The first fitting-only run should compile JuFitter's numerical core without
Makie. The first plotting run can still take noticeably longer because Julia
compiles CairoMakie and its rendering dependencies. This is startup and
precompilation cost, not fit runtime. After the relevant methods are compiled,
repeated fits and plots in the same environment reuse the compiled code.

For a clean local check:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Plot Backends

JuFitter's plotting API is provided by an optional CairoMakie extension. The
core package remains usable for fitting and reporting without CairoMakie. Static
documentation figures use CairoMakie because it is the right default for PNG,
PDF, and SVG output. Later interactive examples can add GLMakie or WGLMakie
workflows, but static output stays the default for reproducible docs.

## Calling JuFitter From Python

The intended Python interoperability path is JuliaCall/PythonCall: Python starts
Julia, loads JuFitter, and calls the same fitting/reporting engine. This is
useful when a project is mostly Python but the fit should use JuFitter's
statistical model and diagnostics.

The minimal example is Makie-free:

```bash
pip install juliacall
python3 examples/python/fit_from_python.py
```

That script activates this repository, calls `fit_model`, and prints
`report_text`. It does not load CairoMakie, so it checks the numerical and text
reporting path without paying plotting compilation cost.

Release policy: Python support is not a public v0 claim until the opt-in gate
has passed in a clean Python environment:

```bash
JUFITTER_RUN_PYTHON_INTEROP=1 julia --project=. --startup-file=no test/python_interop_gate.jl
```

If that gate has not been run with `juliacall` installed, document Python use as
experimental or deferred for the release.

## Troubleshooting

If `using JuFitter` is slow the first time, wait for precompilation to finish.
If every new Julia session is slow, check that you are reusing the same project
environment and not creating a fresh temporary environment each time.

If PDF or SVG export fails, first verify that CairoMakie can render a minimal
figure in the same environment. Most export issues are backend or font related,
not fit related.

If a fit is unexpectedly slow, check whether the problem uses dense covariance
matrices, bounds, constraints, priors, or parameter-dependent x uncertainties.
Those features are supported, but they intentionally move the fit from the fast
least-squares path to the generic optimizer path.
