# JuFitter Example Gallery

Run examples from the repository root:

```bash
julia --project=docs examples/gallery/01_quickstart_linear.jl
```

Generated figures are written to `examples/output/`, which is ignored by git.

## Gallery Structure

- `01_quickstart_linear.jl`: minimal one-call `fitplot(x, y; sigma_y)` workflow.
- `02_xy_uncertainties_photoelectric.jl`: two-regime photoelectric threshold fit with x/y uncertainties and propagated line-intersection uncertainty.
- `03_plot_customization.jl`: theme, units, report placement, sigma band, export, and Makie keyword passthroughs.
- `04_covariance_and_effective_variance.jl`: full y-covariance and effective-variance x-uncertainty examples.
- `05_constraints_priors_profiles.jl`: bounds, inequality constraints, Gaussian priors, profile interval, and contour plot.
- `06_likelihood_workflows.jl`: Poisson, histogram, unbinned, extended-unbinned, indexed, custom, and multi-dataset likelihood fits.
- `07_plot_styles.jl`: controlled comparison of the `:workbench`,
  `:showcase`, and `:publication` contracts with identical scientific content.
- `08_damped_oscillator_decay.jl`: real damped-oscillator data, x/y
  uncertainties, constant-frequency versus frequency-drift model criticism,
  pull diagnostics, and light/dark docs export.
- `09_docs_gallery_suite.jl`: generates the public documentation gallery assets.
- `10_multi_dataset_calibration.jl`: tests full versus partial parameter
  sharing across three calibration channels with joint covariance propagation
  and per-dataset pull diagnostics.

The gallery is intentionally systematic. Avoid adding loose one-off scripts at
the top level; new examples should either extend an existing gallery file or add
a numbered workflow with clear scope.

## Python Interoperability

- `python/fit_from_python.py`: minimal JuliaCall example for Python users. It
  activates this Julia project, fits plain Python arrays, reads parameter
  estimates, and prints `report_text(...)` without loading CairoMakie.

This path is intentionally Julia-backed. JuFitter is not reimplemented in
Python; Python calls the Julia fitting/reporting engine through
`juliacall`.

Prerequisites:

```bash
python3 -m pip install juliacall
```

Run from the repository root:

```bash
python3 examples/python/fit_from_python.py
```

Current limitation: this is an interoperability example and release gate, not a
separate Python package. Startup time is Julia startup plus package loading, and
plotting from Python is deliberately not claimed yet.
