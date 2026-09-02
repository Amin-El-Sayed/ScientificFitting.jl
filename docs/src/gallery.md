# ScientificFitting

```@raw html
<section class="scientificfitting-hero">
  <div class="scientificfitting-kicker">Scientific model fitting for Julia</div>
  <p class="scientificfitting-lede">
    Simple fits stay simple. Difficult fits keep their statistics explicit.
  </p>
  <p>
    One Julia-native workflow for least squares and likelihood fits, x/y and
    correlated uncertainty, parameter constraints, profiles and contours,
    actionable diagnostics, and editable Makie figures.
  </p>
  <div class="scientificfitting-hero-actions">
    <a class="scientificfitting-button primary" href="gallery/linear_calibration.html">Start with a complete fit</a>
    <a class="scientificfitting-button" href="install.html">Install</a>
  </div>
</section>
```

Start with [Linear Calibration](gallery/linear_calibration.md) for the shortest
complete analysis, or choose the example closest to your data below. Every page
contains the measurements, assumptions, executable code, actual program output,
diagnostics, and scientific interpretation.

## Worked Examples

```@raw html
<div class="scientificfitting-gallery-grid">
<div class="scientificfitting-gallery-item">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-linear" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/linear_calibration_sans_panel_light.png" alt="Linear calibration in sans style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-linear" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/linear_calibration_sans_panel_dark.png" alt="Linear calibration in dark sans style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-linear" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/linear_calibration_sans_plot_light.png" alt="Linear calibration in sans style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-linear" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/linear_calibration_sans_plot_dark.png" alt="Linear calibration in dark sans style without result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-linear" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/linear_calibration_tex_panel_light.png" alt="Linear calibration in tex style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-linear" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/linear_calibration_tex_plot_light.png" alt="Linear calibration in tex style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-linear" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/linear_calibration_tex_panel_dark.png" alt="Linear calibration in dark tex style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-linear" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/linear_calibration_tex_plot_dark.png" alt="Linear calibration in dark tex style without result panel">
<div>
<span class="scientificfitting-tag">first fit</span>
<span class="scientificfitting-tag">prediction band</span>
<h3><a href="gallery/linear_calibration.html">Linear calibration</a></h3>
<p>Estimate a calibration law from heteroscedastic measurements. This page is the controlled baseline for reading parameters, bands, residuals, and goodness of fit.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-xy" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/xy_uncertainties_sans_panel_light.png" alt="XY uncertainty fit in sans style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-xy" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/xy_uncertainties_sans_panel_dark.png" alt="XY uncertainty fit in dark sans style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-xy" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/xy_uncertainties_sans_plot_light.png" alt="XY uncertainty fit in sans style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-xy" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/xy_uncertainties_sans_plot_dark.png" alt="XY uncertainty fit in dark sans style without result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-xy" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/xy_uncertainties_tex_panel_light.png" alt="XY uncertainty fit in tex style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-xy" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/xy_uncertainties_tex_plot_light.png" alt="XY uncertainty fit in tex style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-xy" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/xy_uncertainties_tex_panel_dark.png" alt="XY uncertainty fit in dark tex style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-xy" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/xy_uncertainties_tex_plot_dark.png" alt="XY uncertainty fit in dark tex style without result panel">
<div>
<span class="scientificfitting-tag">effective variance</span>
<span class="scientificfitting-tag">x errors</span>
<h3><a href="gallery/xy_uncertainties.html">XY uncertainties</a></h3>
<p>Include uncertainty in the independent variable through the local model slope. This example is useful when calibration, frequency, voltage, or position errors are not negligible.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-covariance" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/full_covariance_decay_sans_panel_light.png" alt="Full covariance fit in sans style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-covariance" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/full_covariance_decay_sans_panel_dark.png" alt="Full covariance fit in dark sans style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-covariance" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/full_covariance_decay_sans_plot_light.png" alt="Full covariance fit in sans style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-covariance" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/full_covariance_decay_sans_plot_dark.png" alt="Full covariance fit in dark sans style without result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-covariance" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/full_covariance_decay_tex_panel_light.png" alt="Full covariance fit in tex style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-covariance" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/full_covariance_decay_tex_plot_light.png" alt="Full covariance fit in tex style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-covariance" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/full_covariance_decay_tex_panel_dark.png" alt="Full covariance fit in dark tex style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-covariance" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/full_covariance_decay_tex_plot_dark.png" alt="Full covariance fit in dark tex style without result panel">
<div>
<span class="scientificfitting-tag">covariance</span>
<span class="scientificfitting-tag">correlations</span>
<h3><a href="gallery/full_covariance.html">Full covariance</a></h3>
<p>Use a dense covariance matrix when measurements share readout noise. The point is not syntax; it is how correlations change parameter uncertainty and goodness-of-fit interpretation.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-damped" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/damped_oscillator_decay_sans_panel_light.png" alt="Damped oscillator fit in sans style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-damped" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/damped_oscillator_decay_sans_panel_dark.png" alt="Damped oscillator fit in dark sans style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-damped" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/damped_oscillator_decay_sans_plot_light.png" alt="Damped oscillator fit in sans style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-damped" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/damped_oscillator_decay_sans_plot_dark.png" alt="Damped oscillator fit in dark sans style without result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-damped" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/damped_oscillator_decay_tex_panel_light.png" alt="Damped oscillator fit in tex style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-damped" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/damped_oscillator_decay_tex_plot_light.png" alt="Damped oscillator fit in tex style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-damped" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/damped_oscillator_decay_tex_panel_dark.png" alt="Damped oscillator fit in dark tex style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-damped" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/damped_oscillator_decay_tex_plot_dark.png" alt="Damped oscillator fit in dark tex style without result panel">
<div>
<span class="scientificfitting-tag">real data</span>
<span class="scientificfitting-tag">nonlinear</span>
<span class="scientificfitting-tag">model criticism</span>
<h3><a href="gallery/resonance_decay.html">Damped oscillator</a></h3>
<p>Discover why a visually convincing constant-frequency fit is statistically rejected, test a frequency-drift extension, and use pull structure to decide what must be investigated next.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-photoelectric" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/photoelectric_threshold_sans_panel_light.png" alt="Photoelectric work-function fit in sans style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-photoelectric" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/photoelectric_threshold_sans_panel_dark.png" alt="Photoelectric work-function fit in dark sans style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-photoelectric" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/photoelectric_threshold_sans_plot_light.png" alt="Photoelectric work-function fit in sans style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-photoelectric" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/photoelectric_threshold_sans_plot_dark.png" alt="Photoelectric work-function fit in dark sans style without result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-photoelectric" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/photoelectric_threshold_tex_panel_light.png" alt="Photoelectric work-function fit in tex style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-photoelectric" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/photoelectric_threshold_tex_plot_light.png" alt="Photoelectric work-function fit in tex style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-photoelectric" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/photoelectric_threshold_tex_panel_dark.png" alt="Photoelectric work-function fit in dark tex style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-photoelectric" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/photoelectric_threshold_tex_plot_dark.png" alt="Photoelectric work-function fit in dark tex style without result panel">
<div>
<span class="scientificfitting-tag">x/y errors</span>
<span class="scientificfitting-tag">line intersection</span>
<h3><a href="gallery/photoelectric_threshold.html">Photoelectric work function</a></h3>
<p>Fit baseline and emission regimes separately, then propagate both covariance matrices into the threshold intersection and work function.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-poisson" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/poisson_counts_sans_panel_light.png" alt="Poisson count fit in sans style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-poisson" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/poisson_counts_sans_panel_dark.png" alt="Poisson count fit in dark sans style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-poisson" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/poisson_counts_sans_plot_light.png" alt="Poisson count fit in sans style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-poisson" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/poisson_counts_sans_plot_dark.png" alt="Poisson count fit in dark sans style without result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-poisson" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/poisson_counts_tex_panel_light.png" alt="Poisson count fit in tex style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-poisson" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/poisson_counts_tex_plot_light.png" alt="Poisson count fit in tex style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-poisson" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/poisson_counts_tex_panel_dark.png" alt="Poisson count fit in dark tex style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-poisson" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/poisson_counts_tex_plot_dark.png" alt="Poisson count fit in dark tex style without result panel">
<div>
<span class="scientificfitting-tag">likelihood</span>
<span class="scientificfitting-tag">counts</span>
<h3><a href="gallery/poisson_histogram.html">Poisson and histograms</a></h3>
<p>Extract a radioactive half-life and a detector peak from sparse counts. Exact count semantics, integrated unequal bins, empty bins, and deviance residuals replace invented Gaussian error bars.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-constraints" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/constraints_priors_sans_panel_light.png" alt="Constrained fit in sans style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-constraints" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/constraints_priors_sans_panel_dark.png" alt="Constrained fit in dark sans style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-constraints" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/constraints_priors_sans_plot_light.png" alt="Constrained fit in sans style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-constraints" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/constraints_priors_sans_plot_dark.png" alt="Constrained fit in dark sans style without result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-constraints" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/constraints_priors_tex_panel_light.png" alt="Constrained fit in tex style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-constraints" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/constraints_priors_tex_plot_light.png" alt="Constrained fit in tex style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-constraints" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/constraints_priors_tex_panel_dark.png" alt="Constrained fit in dark tex style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-constraints" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/constraints_priors_tex_plot_dark.png" alt="Constrained fit in dark tex style without result panel">
<div>
<span class="scientificfitting-tag">constraints</span>
<span class="scientificfitting-tag">profiles</span>
<h3><a href="gallery/constraints_profiles.html">Constraints and profiles</a></h3>
<p>An early saturation measurement leaves amplitude and time constant nonlinearly coupled. Profiles and two-parameter regions show exactly why the local covariance summary fails.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-multi" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/multi_dataset_shared_slope_sans_panel_light.png" alt="Multi-dataset fit in sans style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-multi" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="show" src="assets/gallery/multi_dataset_shared_slope_sans_panel_dark.png" alt="Multi-dataset fit in dark sans style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-multi" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/multi_dataset_shared_slope_sans_plot_light.png" alt="Multi-dataset fit in sans style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-multi" data-scientificfitting-plot-style="sans" data-scientificfitting-plot-panel="hide" src="assets/gallery/multi_dataset_shared_slope_sans_plot_dark.png" alt="Multi-dataset fit in dark sans style without result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-multi" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/multi_dataset_shared_slope_tex_panel_light.png" alt="Multi-dataset fit in tex style with result panel">
<img class="scientificfitting-plot-light" data-scientificfitting-plot-group="gallery-multi" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/multi_dataset_shared_slope_tex_plot_light.png" alt="Multi-dataset fit in tex style without result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-multi" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="show" src="assets/gallery/multi_dataset_shared_slope_tex_panel_dark.png" alt="Multi-dataset fit in dark tex style with result panel">
<img class="scientificfitting-plot-dark" data-scientificfitting-plot-group="gallery-multi" data-scientificfitting-plot-style="tex" data-scientificfitting-plot-panel="hide" src="assets/gallery/multi_dataset_shared_slope_tex_plot_dark.png" alt="Multi-dataset fit in dark tex style without result panel">
<div>
<span class="scientificfitting-tag">multi-fit</span>
<span class="scientificfitting-tag">shared parameters</span>
<span class="scientificfitting-tag">model comparison</span>
<h3><a href="gallery/multi_dataset.html">Multi-dataset fit</a></h3>
<p>Test whether three calibration channels may share one gain, identify the incompatible channel from its pulls, and propagate the gain difference from the joint covariance.</p>
</div>
</div>
</div>
```

Need help judging a result? Continue with [Fitting for Practitioners](@ref).
For derivations, use [Statistical Foundations](@ref); for exact signatures and
defaults, use the [API Reference](@ref).
