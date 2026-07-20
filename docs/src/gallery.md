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
  <div class="jufitter-flow-step"><strong>Measured x</strong><span>XY uncertainties show when horizontal error bars must enter the cost.</span></div>
  <div class="jufitter-flow-step"><strong>Shared noise</strong><span>Full covariance replaces independent errors when observations move together.</span></div>
  <div class="jufitter-flow-step"><strong>Model criticism</strong><span>Damped oscillation uses pull structure to expose missing nonlinear physics.</span></div>
  <div class="jufitter-flow-step"><strong>Derived quantity</strong><span>Photoelectric threshold extraction propagates two fits into one physical result.</span></div>
  <div class="jufitter-flow-step"><strong>Count data</strong><span>Poisson and histogram examples replace Gaussian residuals with count likelihoods.</span></div>
  <div class="jufitter-flow-step"><strong>Nonlinear uncertainty</strong><span>Profiles and contours test whether local symmetric errors remain credible.</span></div>
  <div class="jufitter-flow-step"><strong>Shared hypotheses</strong><span>A multi-dataset fit tests which parameters may be common across experiments.</span></div>
</div>
```

## Current Gallery

```@raw html
<div class="jufitter-gallery-grid">
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="lab" src="assets/gallery/linear_calibration_lab_light.png" alt="Linear calibration in lab style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="lab" src="assets/gallery/linear_calibration_lab_dark.png" alt="Linear calibration in lab dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="modern" src="assets/gallery/linear_calibration_modern_light.png" alt="Linear calibration in modern style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="modern" src="assets/gallery/linear_calibration_modern_dark.png" alt="Linear calibration in modern dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="article" src="assets/gallery/linear_calibration_article_light.png" alt="Linear calibration in article style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="article" src="assets/gallery/linear_calibration_article_dark.png" alt="Linear calibration in article dark style">
<div>
<span class="jufitter-tag">first fit</span>
<span class="jufitter-tag">prediction band</span>
<h3><a href="gallery/linear_calibration.html">Linear calibration</a></h3>
<p>Estimate a calibration law from heteroscedastic measurements. This page is the controlled baseline for reading parameters, bands, residuals, and goodness of fit.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="lab" src="assets/gallery/xy_uncertainties_lab_light.png" alt="XY uncertainty fit in lab style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="lab" src="assets/gallery/xy_uncertainties_lab_dark.png" alt="XY uncertainty fit in lab dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="modern" src="assets/gallery/xy_uncertainties_modern_light.png" alt="XY uncertainty fit in modern style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="modern" src="assets/gallery/xy_uncertainties_modern_dark.png" alt="XY uncertainty fit in modern dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="article" src="assets/gallery/xy_uncertainties_article_light.png" alt="XY uncertainty fit in article style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="article" src="assets/gallery/xy_uncertainties_article_dark.png" alt="XY uncertainty fit in article dark style">
<div>
<span class="jufitter-tag">effective variance</span>
<span class="jufitter-tag">x errors</span>
<h3><a href="gallery/xy_uncertainties.html">XY uncertainties</a></h3>
<p>Include uncertainty in the independent variable through the local model slope. This example is useful when calibration, frequency, voltage, or position errors are not negligible.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="lab" src="assets/gallery/full_covariance_decay_lab_light.png" alt="Full covariance fit in lab style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="lab" src="assets/gallery/full_covariance_decay_lab_dark.png" alt="Full covariance fit in lab dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="modern" src="assets/gallery/full_covariance_decay_modern_light.png" alt="Full covariance fit in modern style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="modern" src="assets/gallery/full_covariance_decay_modern_dark.png" alt="Full covariance fit in modern dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="article" src="assets/gallery/full_covariance_decay_article_light.png" alt="Full covariance fit in article style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="article" src="assets/gallery/full_covariance_decay_article_dark.png" alt="Full covariance fit in article dark style">
<div>
<span class="jufitter-tag">covariance</span>
<span class="jufitter-tag">correlations</span>
<h3><a href="gallery/full_covariance.html">Full covariance</a></h3>
<p>Use a dense covariance matrix when measurements share readout noise. The point is not syntax; it is how correlations change parameter uncertainty and goodness-of-fit interpretation.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="lab" src="assets/gallery/damped_oscillator_decay_lab_light.png" alt="Damped oscillator fit in lab style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="lab" src="assets/gallery/damped_oscillator_decay_lab_dark.png" alt="Damped oscillator fit in lab dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="modern" src="assets/gallery/damped_oscillator_decay_modern_light.png" alt="Damped oscillator fit in modern style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="modern" src="assets/gallery/damped_oscillator_decay_modern_dark.png" alt="Damped oscillator fit in modern dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="article" src="assets/gallery/damped_oscillator_decay_article_light.png" alt="Damped oscillator fit in article style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="article" src="assets/gallery/damped_oscillator_decay_article_dark.png" alt="Damped oscillator fit in article dark style">
<div>
<span class="jufitter-tag">real data</span>
<span class="jufitter-tag">nonlinear</span>
<span class="jufitter-tag">model criticism</span>
<h3><a href="gallery/resonance_decay.html">Damped oscillator</a></h3>
<p>Discover why a visually convincing constant-frequency fit is statistically rejected, test a frequency-drift extension, and use pull structure to decide what must be investigated next.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="lab" src="assets/gallery/photoelectric_threshold_lab_light.png" alt="Photoelectric work-function fit in lab style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="lab" src="assets/gallery/photoelectric_threshold_lab_dark.png" alt="Photoelectric work-function fit in lab dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="modern" src="assets/gallery/photoelectric_threshold_modern_light.png" alt="Photoelectric work-function fit in modern style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="modern" src="assets/gallery/photoelectric_threshold_modern_dark.png" alt="Photoelectric work-function fit in modern dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="article" src="assets/gallery/photoelectric_threshold_article_light.png" alt="Photoelectric work-function fit in article style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="article" src="assets/gallery/photoelectric_threshold_article_dark.png" alt="Photoelectric work-function fit in article dark style">
<div>
<span class="jufitter-tag">x/y errors</span>
<span class="jufitter-tag">line intersection</span>
<h3><a href="gallery/photoelectric_threshold.html">Photoelectric work function</a></h3>
<p>Fit baseline and emission regimes separately, then propagate both covariance matrices into the threshold intersection and work function.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="lab" src="assets/gallery/poisson_counts_lab_light.png" alt="Poisson count fit in lab style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="lab" src="assets/gallery/poisson_counts_lab_dark.png" alt="Poisson count fit in lab dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="modern" src="assets/gallery/poisson_counts_modern_light.png" alt="Poisson count fit in modern style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="modern" src="assets/gallery/poisson_counts_modern_dark.png" alt="Poisson count fit in modern dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="article" src="assets/gallery/poisson_counts_article_light.png" alt="Poisson count fit in article style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="article" src="assets/gallery/poisson_counts_article_dark.png" alt="Poisson count fit in article dark style">
<div>
<span class="jufitter-tag">likelihood</span>
<span class="jufitter-tag">counts</span>
<h3><a href="gallery/poisson_histogram.html">Poisson and histograms</a></h3>
<p>Extract a radioactive half-life and a detector peak from sparse counts. Exact count semantics, integrated unequal bins, empty bins, and deviance residuals replace invented Gaussian error bars.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="lab" src="assets/gallery/constraints_priors_lab_light.png" alt="Constrained fit in lab style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="lab" src="assets/gallery/constraints_priors_lab_dark.png" alt="Constrained fit in lab dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="modern" src="assets/gallery/constraints_priors_modern_light.png" alt="Constrained fit in modern style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="modern" src="assets/gallery/constraints_priors_modern_dark.png" alt="Constrained fit in modern dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="article" src="assets/gallery/constraints_priors_article_light.png" alt="Constrained fit in article style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="article" src="assets/gallery/constraints_priors_article_dark.png" alt="Constrained fit in article dark style">
<div>
<span class="jufitter-tag">constraints</span>
<span class="jufitter-tag">profiles</span>
<h3><a href="gallery/constraints_profiles.html">Constraints and profiles</a></h3>
<p>An early saturation measurement leaves amplitude and time constant nonlinearly coupled. Profiles and two-parameter regions show exactly why the local covariance summary fails.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="lab" src="assets/gallery/multi_dataset_shared_slope_lab_light.png" alt="Multi-dataset fit in lab style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="lab" src="assets/gallery/multi_dataset_shared_slope_lab_dark.png" alt="Multi-dataset fit in lab dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="modern" src="assets/gallery/multi_dataset_shared_slope_modern_light.png" alt="Multi-dataset fit in modern style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="modern" src="assets/gallery/multi_dataset_shared_slope_modern_dark.png" alt="Multi-dataset fit in modern dark style">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="article" src="assets/gallery/multi_dataset_shared_slope_article_light.png" alt="Multi-dataset fit in article style">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="article" src="assets/gallery/multi_dataset_shared_slope_article_dark.png" alt="Multi-dataset fit in article dark style">
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
| Poisson and histograms | counts, rates, binned events | Poisson deviance and likelihood fits | low-count bins, deviance residuals, empty-bin behavior |
| Constraints and profiles | early saturation measurement | bounds, prior information, profile intervals, non-elliptic contours | unseen plateau and amplitude-timescale degeneracy |
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
- `07_plot_styles.jl`: controlled comparison of the `:lab`,
  `:modern`, and `:article` contracts using identical scientific
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
