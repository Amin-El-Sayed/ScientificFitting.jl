# JuFitter Example Gallery

Run examples from the repository root:

```bash
julia --project=. examples/gallery/01_quickstart_linear.jl
```

Generated figures are written to `examples/output/`, which is ignored by git.

## Gallery Structure

- `01_quickstart_linear.jl`: minimal one-call `fitplot(x, y; sigma_y)` workflow.
- `02_xy_uncertainties_photoelectric.jl`: two-regime photoelectric threshold fit with x/y uncertainties and propagated line-intersection uncertainty.
- `03_plot_customization.jl`: theme, units, report placement, sigma band, export, and Makie keyword passthroughs.
- `04_covariance_and_effective_variance.jl`: full y-covariance and effective-variance x-uncertainty examples.
- `05_constraints_priors_profiles.jl`: bounds, inequality constraints, Gaussian priors, profile interval, and contour plot.
- `06_likelihood_workflows.jl`: Poisson, histogram, unbinned, extended-unbinned, indexed, custom, and multi-dataset likelihood fits.
- `07_plot_styles.jl`: controlled comparison of every public plot style with identical scientific content.
- `08_damped_oscillator_decay.jl`: real damped-oscillator data, x/y
  uncertainties, constant-frequency versus frequency-drift model criticism,
  pull diagnostics, and light/dark docs export.
- `09_docs_gallery_suite.jl`: generates the public documentation gallery assets.

The gallery is intentionally systematic. Avoid adding loose one-off scripts at
the top level; new examples should either extend an existing gallery file or add
a numbered workflow with clear scope.
