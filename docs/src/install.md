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
