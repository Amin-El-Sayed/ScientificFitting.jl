"""Call JuFitter from Python through JuliaCall.

This example intentionally uses only the fitting, text-reporting, and diagnostic
API. It does not load CairoMakie/Makie, so it is the smallest realistic
interoperability check for Python users who want JuFitter's numerical engine
without Julia-side plotting.

Prerequisites:
    pip install juliacall

Run from the repository root:
    python3 examples/python/fit_from_python.py
"""

from pathlib import Path

from juliacall import Main as jl


ROOT = Path(__file__).resolve().parents[2]

jl.seval("using Pkg")
jl.Pkg.activate(str(ROOT))
jl.seval("using JuFitter")
jl.seval("linear_model(x, p) = @. p[1] * x + p[2]")
jl.seval(
    """
    function loaded_plot_modules_for_python()
        names = String[]
        for pkg in keys(Base.loaded_modules)
            pkg.name in ("Makie", "CairoMakie") && push!(names, pkg.name)
        end
        return sort(names)
    end
    """
)

x = [0.0, 1.0, 2.0, 3.0, 4.0]
y = [1.02, 2.95, 5.05, 7.02, 8.96]
sigma_y = [0.12, 0.12, 0.13, 0.13, 0.14]

result = jl.JuFitter.fit_model(
    jl.linear_model,
    x,
    y,
    p0=[1.0, 0.0],
    sigma_y=sigma_y,
)

params = list(result.params)
stderr = list(result.param_stderr)

print("JuFitter from Python")
print(f"slope = {params[0]:.8g} +/- {stderr[0]:.3g}")
print(f"offset = {params[1]:.8g} +/- {stderr[1]:.3g}")
print(jl.JuFitter.report_text(result, parameter_names=["slope", "offset"]))
print(jl.JuFitter.diagnostic_dashboard_text(result))

# This is a release contract: fitting and text diagnostics must stay usable from
# Python without pulling in the plotting backend.
plot_modules = list(jl.loaded_plot_modules_for_python())
print(f"plot modules loaded = {plot_modules}")
if plot_modules:
    raise RuntimeError(
        "Python fitting/reporting unexpectedly loaded plotting modules: "
        + ", ".join(plot_modules)
    )
