# Gallery

The gallery is the fastest way to learn JuFitter. It is organized as a
progression from ordinary laboratory fits to correlated uncertainties,
nonlinear models, likelihoods, profiles, and multi-dataset hypotheses. Start
near the top, run the complete script, and treat each page as a worked
scientific analysis rather than a syntax sample.

Every finished gallery page follows the same pattern:

1. a concrete question,
2. measured quantities with units and uncertainties,
3. the model and its assumptions,
4. complete executable code,
5. a plot with visible uncertainty semantics,
6. fit diagnostics,
7. interpretation and realistic failure modes.

If you only want the shortest path from two arrays to a polished fit, start with
[Quickstart](@ref). If the fit looks suspicious, continue with
[Fitting for Practitioners](@ref). If you need to justify the method in a
report, read [Statistical Foundations](@ref).

The complete scripts can be run from the project root:

```bash
julia --project=docs examples/gallery/01_quickstart_linear.jl
```

Generated figures are written to `examples/output/`. Documentation assets under
`docs/src/assets/gallery/` are regenerated only when the public examples or plot
style contracts change.

## Recommended Path

```@raw html
<div class="jufitter-flow">
  <div class="jufitter-flow-step"><strong>First fit</strong><span>Linear calibration teaches the default data-model-diagnosis loop.</span></div>
  <div class="jufitter-flow-step"><strong>Uncertainty structure</strong><span>XY errors and full covariance explain when error bars are not independent decoration.</span></div>
  <div class="jufitter-flow-step"><strong>Real experiment</strong><span>Damped oscillation shows nonlinear parameters, units, and imperfect residuals.</span></div>
  <div class="jufitter-flow-step"><strong>Derived quantities</strong><span>Photoelectric threshold extraction adds regime choice and propagated intersections.</span></div>
  <div class="jufitter-flow-step"><strong>Likelihoods</strong><span>Poisson and histogram examples switch from least squares to count statistics.</span></div>
  <div class="jufitter-flow-step"><strong>Trust check</strong><span>Profiles, contours, and diagnostics test whether local symmetric errors are credible.</span></div>
</div>
```

## Current Gallery

```@raw html
<div class="jufitter-gallery-grid">
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="workbench" src="assets/gallery/linear_calibration_workbench_light.png" alt="Linear calibration in workbench style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="workbench" src="assets/gallery/linear_calibration_workbench_dark.png" alt="Linear calibration in workbench dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="showcase" src="assets/gallery/linear_calibration_showcase_light.png" alt="Linear calibration in showcase style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="showcase" src="assets/gallery/linear_calibration_showcase_dark.png" alt="Linear calibration in showcase dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="publication" src="assets/gallery/linear_calibration_publication_light.png" alt="Linear calibration in publication style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="publication" src="assets/gallery/linear_calibration_publication_dark.png" alt="Linear calibration in publication dark style">
<div>
<span class="jufitter-tag">first fit</span>
<span class="jufitter-tag">prediction band</span>
<h3><a href="gallery/linear_calibration.html">Linear calibration</a></h3>
<p>Estimate a calibration law from heteroscedastic measurements. This page is the controlled baseline for reading parameters, bands, residuals, and goodness of fit.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="workbench" src="assets/gallery/xy_uncertainties_workbench_light.png" alt="XY uncertainty fit in workbench style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="workbench" src="assets/gallery/xy_uncertainties_workbench_dark.png" alt="XY uncertainty fit in workbench dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="showcase" src="assets/gallery/xy_uncertainties_showcase_light.png" alt="XY uncertainty fit in showcase style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="showcase" src="assets/gallery/xy_uncertainties_showcase_dark.png" alt="XY uncertainty fit in showcase dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="publication" src="assets/gallery/xy_uncertainties_publication_light.png" alt="XY uncertainty fit in publication style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="publication" src="assets/gallery/xy_uncertainties_publication_dark.png" alt="XY uncertainty fit in publication dark style">
<div>
<span class="jufitter-tag">effective variance</span>
<span class="jufitter-tag">x errors</span>
<h3><a href="gallery/xy_uncertainties.html">XY uncertainties</a></h3>
<p>Include uncertainty in the independent variable through the local model slope. This example is useful when calibration, frequency, voltage, or position errors are not negligible.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="workbench" src="assets/gallery/full_covariance_decay_workbench_light.png" alt="Full covariance fit in workbench style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="workbench" src="assets/gallery/full_covariance_decay_workbench_dark.png" alt="Full covariance fit in workbench dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="showcase" src="assets/gallery/full_covariance_decay_showcase_light.png" alt="Full covariance fit in showcase style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="showcase" src="assets/gallery/full_covariance_decay_showcase_dark.png" alt="Full covariance fit in showcase dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="publication" src="assets/gallery/full_covariance_decay_publication_light.png" alt="Full covariance fit in publication style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="publication" src="assets/gallery/full_covariance_decay_publication_dark.png" alt="Full covariance fit in publication dark style">
<div>
<span class="jufitter-tag">covariance</span>
<span class="jufitter-tag">correlations</span>
<h3><a href="gallery/full_covariance.html">Full covariance</a></h3>
<p>Use a dense covariance matrix when measurements share readout noise. The point is not syntax; it is how correlations change parameter uncertainty and goodness-of-fit interpretation.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="workbench" src="assets/gallery/damped_oscillator_decay_workbench_light.png" alt="Damped oscillator fit in workbench style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="workbench" src="assets/gallery/damped_oscillator_decay_workbench_dark.png" alt="Damped oscillator fit in workbench dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="showcase" src="assets/gallery/damped_oscillator_decay_showcase_light.png" alt="Damped oscillator fit in showcase style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="showcase" src="assets/gallery/damped_oscillator_decay_showcase_dark.png" alt="Damped oscillator fit in showcase dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="publication" src="assets/gallery/damped_oscillator_decay_publication_light.png" alt="Damped oscillator fit in publication style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="publication" src="assets/gallery/damped_oscillator_decay_publication_dark.png" alt="Damped oscillator fit in publication dark style">
<div>
<span class="jufitter-tag">real data</span>
<span class="jufitter-tag">nonlinear</span>
<span class="jufitter-tag">model criticism</span>
<h3><a href="gallery/resonance_decay.html">Damped oscillator</a></h3>
<p>Discover why a visually convincing constant-frequency fit is statistically rejected, test a frequency-drift extension, and use pull structure to decide what must be investigated next.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="workbench" src="assets/gallery/photoelectric_threshold_workbench_light.png" alt="Photoelectric work-function fit in workbench style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="workbench" src="assets/gallery/photoelectric_threshold_workbench_dark.png" alt="Photoelectric work-function fit in workbench dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="showcase" src="assets/gallery/photoelectric_threshold_showcase_light.png" alt="Photoelectric work-function fit in showcase style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="showcase" src="assets/gallery/photoelectric_threshold_showcase_dark.png" alt="Photoelectric work-function fit in showcase dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="publication" src="assets/gallery/photoelectric_threshold_publication_light.png" alt="Photoelectric work-function fit in publication style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="publication" src="assets/gallery/photoelectric_threshold_publication_dark.png" alt="Photoelectric work-function fit in publication dark style">
<div>
<span class="jufitter-tag">x/y errors</span>
<span class="jufitter-tag">line intersection</span>
<h3><a href="gallery/photoelectric_threshold.html">Photoelectric work function</a></h3>
<p>Fit baseline and emission regimes separately, then propagate both covariance matrices into the threshold intersection and work function.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="workbench" src="assets/gallery/constraints_priors_workbench_light.png" alt="Constrained fit in workbench style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="workbench" src="assets/gallery/constraints_priors_workbench_dark.png" alt="Constrained fit in workbench dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="showcase" src="assets/gallery/constraints_priors_showcase_light.png" alt="Constrained fit in showcase style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="showcase" src="assets/gallery/constraints_priors_showcase_dark.png" alt="Constrained fit in showcase dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="publication" src="assets/gallery/constraints_priors_publication_light.png" alt="Constrained fit in publication style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="publication" src="assets/gallery/constraints_priors_publication_dark.png" alt="Constrained fit in publication dark style">
<div>
<span class="jufitter-tag">constraints</span>
<span class="jufitter-tag">profiles</span>
<h3><a href="gallery/constraints_profiles.html">Constraints and profiles</a></h3>
<p>An early saturation measurement leaves amplitude and time constant nonlinearly coupled. Profiles and two-parameter regions show exactly why the local covariance summary fails.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="workbench" src="assets/gallery/poisson_counts_workbench_light.png" alt="Poisson count fit in workbench style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="workbench" src="assets/gallery/poisson_counts_workbench_dark.png" alt="Poisson count fit in workbench dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="showcase" src="assets/gallery/poisson_counts_showcase_light.png" alt="Poisson count fit in showcase style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="showcase" src="assets/gallery/poisson_counts_showcase_dark.png" alt="Poisson count fit in showcase dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="publication" src="assets/gallery/poisson_counts_publication_light.png" alt="Poisson count fit in publication style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="publication" src="assets/gallery/poisson_counts_publication_dark.png" alt="Poisson count fit in publication dark style">
<div>
<span class="jufitter-tag">likelihood</span>
<span class="jufitter-tag">counts</span>
<h3><a href="gallery/poisson_histogram.html">Poisson and histograms</a></h3>
<p>Extract a radioactive half-life and a detector peak from sparse counts. Exact count semantics, integrated unequal bins, empty bins, and deviance residuals replace invented Gaussian error bars.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="workbench" src="assets/gallery/multi_dataset_shared_slope_workbench_light.png" alt="Multi-dataset fit in workbench style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="workbench" src="assets/gallery/multi_dataset_shared_slope_workbench_dark.png" alt="Multi-dataset fit in workbench dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="showcase" src="assets/gallery/multi_dataset_shared_slope_showcase_light.png" alt="Multi-dataset fit in showcase style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="showcase" src="assets/gallery/multi_dataset_shared_slope_showcase_dark.png" alt="Multi-dataset fit in showcase dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="publication" src="assets/gallery/multi_dataset_shared_slope_publication_light.png" alt="Multi-dataset fit in publication style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="publication" src="assets/gallery/multi_dataset_shared_slope_publication_dark.png" alt="Multi-dataset fit in publication dark style">
<div>
<span class="jufitter-tag">multi-fit</span>
<span class="jufitter-tag">shared parameters</span>
<span class="jufitter-tag">model comparison</span>
<h3><a href="gallery/multi_dataset.html">Multi-dataset fit</a></h3>
<p>Test whether three calibration channels may share one gain, identify the incompatible channel from its pulls, and propagate the gain difference from the joint covariance.</p>
</div>
</div>
</div>
```

## What Each Example Teaches

| Page | Scientific use case | Statistical focus | Diagnosis to inspect |
| --- | --- | --- | --- |
| Linear calibration | sensor or scale calibration | weighted Gaussian least squares | residual structure, prediction band, chi-square per degree of freedom |
| XY uncertainties | calibration with uncertain abscissa | effective-variance approximation | slope-dependent variance and approximation validity |
| Full covariance | repeated readout with shared noise | whitening with dense covariance | correlation effect on uncertainty and p-value |
| Damped oscillator | mechanical decay or resonance envelope | nonlinear least squares with correlated parameters | phase/frequency coupling, residual periodicity, profile asymmetry |
| Photoelectric work function | threshold from two fitted regimes | x/y uncertainty and derived quantity propagation | regime choice, intersection uncertainty, model range |
| Constraints and profiles | early saturation measurement | bounds, prior information, profile intervals, non-elliptic contours | unseen plateau and amplitude-timescale degeneracy |
| Poisson and histograms | counts, rates, binned events | Poisson deviance and likelihood fits | low-count bins, deviance residuals, empty-bin behavior |
| Multi-dataset fit | shared physics across runs | parameter mapping and joint costs | per-dataset residuals and shared-parameter tension |

## Executable Scripts

The pages are written for reading; the scripts under `examples/gallery/` are
written for running. They use the same models, uncertainty assumptions, and
diagnostics, but keep the surrounding explanation out of the source file.

- `01_quickstart_linear.jl`: minimal weighted Gaussian fit with visible
  uncertainty semantics.
- `02_xy_uncertainties_photoelectric.jl`: photoelectric threshold fit with x/y
  uncertainties, two fitted regimes, and propagated intersection uncertainty.
- `03_plot_customization.jl`: themes, units, reports, sigma bands, export, and
  Makie keyword passthrough.
- `04_covariance_and_effective_variance.jl`: dense y covariance and
  effective-variance x uncertainty.
- `05_constraints_priors_profiles.jl`: bounds, inequality constraints, Gaussian
  priors, profile intervals, and contours.
- `06_likelihood_workflows.jl`: radioactive decay counts, an integrated
  detector spectrum, unbinned, extended-unbinned, indexed, custom, and
  multi-dataset likelihood fits.
- `07_plot_styles.jl`: controlled comparison of the `:workbench`,
  `:showcase`, and `:publication` contracts using identical scientific
  content.
- `08_damped_oscillator_decay.jl`: real mechanical oscillator decay with x/y
  uncertainties and light/dark documentation exports.
- `09_docs_gallery_suite.jl`: reproducible documentation asset generator.
- `10_multi_dataset_calibration.jl`: full versus partial parameter sharing,
  per-dataset pulls, and joint-covariance propagation.

## Editorial Standard

A release-ready gallery page must let a reader copy the full code, reproduce
the plot, understand the statistical assumptions, and decide whether the fit
should be trusted. The repository enforces the structural part of that standard
with a gallery gate:

```bash
julia --project=. test/docs_gallery_gate.jl
julia --project=. test/docs_output_snapshots.jl
```

The gate checks that every public gallery page has a scientific question, data,
model or cost description, executable code, real output, diagnostics,
interpretation, failure modes, complete light/dark plot assets, and all three
public plot-style variants. Visual taste still requires human review, but the
site should no longer regress into missing pages, absent output blocks, broken
assets, or incomplete tutorial fragments. The output-snapshot gate is slower:
it executes the documented example scripts and verifies that notebook-style
`Real output` blocks are copied from real script output rather than edited by
hand.
