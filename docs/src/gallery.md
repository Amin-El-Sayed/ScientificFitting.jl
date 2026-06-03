# Gallery

The gallery is the public face of JuFitter. It should demonstrate scientific
workflows, not only API calls. Every major example must have a real question,
visible uncertainties, a fit result, and a plot that can stand on its own.

Run repository examples from the project root:

```bash
julia --project=. examples/gallery/01_quickstart_linear.jl
```

Generated figures are written to `examples/output/`.

## Current Gallery

```@raw html
<div class="jufitter-gallery-grid">
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" src="assets/gallery/linear_calibration_light.png" alt="Linear calibration fit">
<img class="jufitter-plot-dark" src="assets/gallery/linear_calibration_dark.png" alt="Linear calibration fit dark">
<div>
<span class="jufitter-tag">quickstart</span>
<span class="jufitter-tag">prediction band</span>
<h3><a href="gallery/linear_calibration.html">Linear calibration</a></h3>
<p>One-call Gaussian fit with heteroscedastic uncertainties, report panel, and visible prediction band.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" src="assets/gallery/photoelectric_threshold_light.png" alt="Photoelectric work-function fit">
<img class="jufitter-plot-dark" src="assets/gallery/photoelectric_threshold_dark.png" alt="Photoelectric work-function fit dark">
<div>
<span class="jufitter-tag">x/y errors</span>
<span class="jufitter-tag">extrapolation</span>
<h3><a href="gallery/photoelectric_threshold.html">Photoelectric work function</a></h3>
<p>Linear photoelectric fit with threshold marker, y-intercept work function, and propagated uncertainty.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" src="assets/gallery/damped_oscillator_decay_light.png" alt="Damped oscillator fit">
<img class="jufitter-plot-dark" src="assets/gallery/damped_oscillator_decay_dark.png" alt="Damped oscillator fit dark">
<div>
<span class="jufitter-tag">real data</span>
<span class="jufitter-tag">nonlinear</span>
<h3><a href="gallery/resonance_decay.html">Damped oscillator</a></h3>
<p>Mechanical oscillator decay with visible uncertainty band, nonlinear model,
and goodness-of-fit diagnostics.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" src="assets/gallery/full_covariance_decay_light.png" alt="Full covariance decay fit">
<img class="jufitter-plot-dark" src="assets/gallery/full_covariance_decay_dark.png" alt="Full covariance decay fit dark">
<div>
<span class="jufitter-tag">covariance</span>
<span class="jufitter-tag">correlations</span>
<h3><a href="gallery/full_covariance.html">Full covariance</a></h3>
<p>Exponential decay with a dense y-covariance matrix and correlated readout noise.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" src="assets/gallery/xy_uncertainties_light.png" alt="XY uncertainty fit">
<img class="jufitter-plot-dark" src="assets/gallery/xy_uncertainties_dark.png" alt="XY uncertainty fit dark">
<div>
<span class="jufitter-tag">effective variance</span>
<span class="jufitter-tag">x errors</span>
<h3><a href="gallery/xy_uncertainties.html">XY uncertainties</a></h3>
<p>Fit where x errors contribute through the local model slope, not just as cosmetic error bars.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" src="assets/gallery/constraints_priors_light.png" alt="Constrained fit">
<img class="jufitter-plot-dark" src="assets/gallery/constraints_priors_dark.png" alt="Constrained fit dark">
<div>
<span class="jufitter-tag">constraints</span>
<span class="jufitter-tag">profiles</span>
<h3><a href="gallery/constraints_profiles.html">Constraints and profiles</a></h3>
<p>Bounds, inequality constraints, Gaussian priors, profile scans, and two-parameter contours.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" src="assets/gallery/poisson_counts_light.png" alt="Poisson count fit">
<img class="jufitter-plot-dark" src="assets/gallery/poisson_counts_dark.png" alt="Poisson count fit dark">
<div>
<span class="jufitter-tag">likelihood</span>
<span class="jufitter-tag">counts</span>
<h3><a href="gallery/poisson_histogram.html">Poisson and histograms</a></h3>
<p>Count data and binned likelihood fits where Gaussian least squares is the wrong default.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" src="assets/gallery/multi_dataset_shared_slope_light.png" alt="Multi dataset fit">
<img class="jufitter-plot-dark" src="assets/gallery/multi_dataset_shared_slope_dark.png" alt="Multi dataset fit dark">
<div>
<span class="jufitter-tag">multi-fit</span>
<span class="jufitter-tag">shared parameters</span>
<h3><a href="gallery/multi_dataset.html">Multi-dataset fit</a></h3>
<p>Two datasets share a slope while each keeps its own offset and uncertainty model.</p>
</div>
</div>
</div>
```

The repository scripts are currently:

- `01_quickstart_linear.jl`: minimal linear Gaussian fit with uncertainty band.
- `02_xy_uncertainties_photoelectric.jl`: physics example with x/y
  uncertainties, bounds, multistart, LaTeX labels, and diagnostics.
- `03_plot_customization.jl`: themes, units, report placement, sigma bands,
  export, and Makie keyword passthrough.
- `04_covariance_and_effective_variance.jl`: full y covariance and
  effective-variance x uncertainty.
- `05_constraints_priors_profiles.jl`: bounds, inequality constraints,
  Gaussian priors, profile interval, and contour plot.
- `06_likelihood_workflows.jl`: Poisson, histogram, unbinned,
  extended-unbinned, indexed, custom, and multi-dataset likelihood fits.
- `07_plot_styles.jl`: clean default, dense-data minimal, and LaTeX-style paper
  plot layouts.
- `08_damped_oscillator_decay.jl`: real mechanical oscillator decay with x/y
  uncertainties, LaTeX labels, and light/dark documentation exports.
- `09_docs_gallery_suite.jl`: documentation asset generator for the public
  gallery pages.

## Planned Real-Data Gallery

The next gallery pass should convert curated public datasets into reproducible
examples and documentation pages. Each dataset must be renamed by physical
content, reduced to the columns needed for the example, and documented without
local paths, course-internal context, or unrelated acquisition notes.

## Required Example Types

- Simple linear fit with visible uncertainties and readable result panel.
- Pendulum period analysis: small-angle model, improved nonlinear correction,
  residual diagnostics, and model comparison.
- Driven or damped resonance: nonlinear line shape, damping, phase, and
  physically interpretable parameters.
- Interference or diffraction: oscillatory model with calibration and
  parameter correlations.
- Signal analysis: noisy time series, frequency extraction, and fit residuals.
- Photoelectric-effect workflow with multiple fits in one plot and a visible
  work-function marker.
- Circuit example: RC/RLC transient or frequency response with units and
  derived quantities.
- Histogram or count-data example with Poisson statistics and likelihood
  diagnostics.
- One non-physics example from engineering, biology, chemistry, or social data
  to show that JuFitter is not laboratory-course specific.
- Dense data fit with minimal markers, clear line, and no visual clutter.
- XY fit with x and y uncertainties.
- Full covariance fit with a visible correlation effect.
- Fit with priors, fixed parameters, bounds, and diagnostics.
- Profile interval and two-parameter contour.
- Poisson count fit and histogram fit.
- Unbinned and extended-unbinned likelihood examples.
- Multi-dataset fit with shared and dataset-specific parameters.
- Plot theme comparison including a dark-mode export.

## Visual Rules

- Data should not look like decoration. Marker size, alpha, line width, and
  bands must be chosen for the dataset density.
- Uncertainty bands must be visible enough to teach the concept.
- Result panels should help, not dominate. Large datasets need a quieter layout.
- Every gallery plot must work on light and dark documentation backgrounds.
- If a plot needs manual margins, that is a plotting bug or missing API option,
  not an acceptable tutorial workaround.
