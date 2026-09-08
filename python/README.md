# ScientificFitting for Python

Fit NumPy models to measurement data with Gaussian errors, Poisson counts, or
your own likelihood. Bounds, shared parameters, profiles, and diagnostic
reports use the same numerical core as ScientificFitting.jl. Optional plots
are ordinary, editable Matplotlib figures, not Julia/Makie objects.

## Installation

In a Python 3.10+ environment:

```sh
python -m pip install 'scientificfitting[plot]'
```

JuliaCall installs a compatible Julia runtime and dependencies automatically.
First use requires internet access and compilation; subsequent fits reuse
that installation. The wheel is small, but Julia and its dependencies are
separate downloads and take additional disk space. No manual Julia setup
is required. Matplotlib and SciPy are optional (`plot` and `sparse` extras).

```python
import numpy as np
from scientificfitting import fit_model, plot_fit

def line(x, slope, offset):
    return slope * x + offset

x = np.array([0., 1., 2., 3.])
y = np.array([0.1, 1.2, 1.9, 3.2])
fit = fit_model(line, x, y, p0={"slope": 1., "offset": 0.}, sigma_y=0.2)
print(fit.report())
fig, ax = plot_fit(fit, xlabel="x / mm", ylabel="U / V")
ax.axvline(1.5, color="black", linestyle="--")
fig.savefig("calibration.pdf")
```

Measurement errors are inputs. Reported parameter errors are local covariance
approximations; profiles help examine asymmetry and non-quadratic behavior.
This is likelihood optimization, not posterior sampling. Python callbacks
use finite derivatives or supplied analytic Jacobians, not Julia dual numbers.

See the [Python guide](https://amin-el-sayed.github.io/ScientificFitting.jl/python.html)
for likelihoods, covariance, diagnostics, and native Matplotlib composition.
[Bug reports and scientific use cases](https://github.com/Amin-El-Sayed/ScientificFitting.jl/issues)
are welcome. MIT licensed, copyright Amin El Sayed.

## Development

To use a source checkout instead of the registered Julia core:

```sh
python -m pip install -e './python[plot,test]'
python python/develop.py
python -m pytest python/tests
```
