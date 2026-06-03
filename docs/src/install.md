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

When JuFitter is registered, the intended user-facing installation path is:

```julia
using Pkg
Pkg.add("JuFitter")
```

## First Compile

The first run can take noticeably longer because Julia compiles JuFitter,
Makie, Optimization.jl, and their dependencies. This is startup and
precompilation cost, not fit runtime. After the package is compiled, repeated
fits in the same environment use the compiled methods.

For a clean local check:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Plot Backends

JuFitter currently uses CairoMakie for static publication output. It is the
right default for PNG, PDF, and SVG documentation figures. Later interactive
examples can add GLMakie or WGLMakie workflows, but static output stays the
default for reproducible docs.

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
