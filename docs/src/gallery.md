# Gallery

The gallery is the fastest way to learn JuFitter. It is organized as a
progression from ordinary laboratory fits to correlated uncertainties,
nonlinear models, likelihoods, profiles, and multi-dataset hypotheses. Start
near the top, run the complete script, and treat each page as a worked
scientific analysis rather than a syntax sample.

Every gallery page follows the same pattern:

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

Each page contains the complete code used for its analysis and figure. From a
source checkout, the corresponding script can also be run from the project
root:

```bash
julia --project=docs examples/gallery/01_quickstart_linear.jl
```

The script writes its figure to `examples/output/` and prints the same fit
summary and diagnostic messages shown on the page.

## Recommended Path

```@raw html
<ol class="jufitter-flow">
  <li class="jufitter-flow-step"><strong>First fit</strong><span>Linear calibration teaches the default data-model-diagnosis loop.</span></li>
  <li class="jufitter-flow-step"><strong>Measured x</strong><span>XY uncertainties show when horizontal error bars must enter the cost.</span></li>
  <li class="jufitter-flow-step"><strong>Shared noise</strong><span>Full covariance replaces independent errors when observations move together.</span></li>
  <li class="jufitter-flow-step"><strong>Model criticism</strong><span>Damped oscillation uses pull structure to expose missing nonlinear physics.</span></li>
  <li class="jufitter-flow-step"><strong>Derived quantity</strong><span>Photoelectric threshold extraction propagates two fits into one physical result.</span></li>
  <li class="jufitter-flow-step"><strong>Count data</strong><span>Poisson and histogram examples replace Gaussian residuals with count likelihoods.</span></li>
  <li class="jufitter-flow-step"><strong>Nonlinear uncertainty</strong><span>Profiles and contours test whether local symmetric errors remain credible.</span></li>
  <li class="jufitter-flow-step"><strong>Shared hypotheses</strong><span>A multi-dataset fit tests which parameters may be common across experiments.</span></li>
</ol>
```

## Examples

```@raw html
<div class="jufitter-gallery-grid">
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/linear_calibration_sans_panel_light.png" alt="Linear calibration in sans style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/linear_calibration_sans_panel_dark.png" alt="Linear calibration in dark sans style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/linear_calibration_sans_plot_light.png" alt="Linear calibration in sans style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/linear_calibration_sans_plot_dark.png" alt="Linear calibration in dark sans style without result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/linear_calibration_tex_panel_light.png" alt="Linear calibration in tex style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/linear_calibration_tex_plot_light.png" alt="Linear calibration in tex style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/linear_calibration_tex_panel_dark.png" alt="Linear calibration in dark tex style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-linear" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/linear_calibration_tex_plot_dark.png" alt="Linear calibration in dark tex style without result panel">
<div>
<span class="jufitter-tag">first fit</span>
<span class="jufitter-tag">prediction band</span>
<h3><a href="gallery/linear_calibration.html">Linear calibration</a></h3>
<p>Estimate a calibration law from heteroscedastic measurements. This page is the controlled baseline for reading parameters, bands, residuals, and goodness of fit.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/xy_uncertainties_sans_panel_light.png" alt="XY uncertainty fit in sans style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/xy_uncertainties_sans_panel_dark.png" alt="XY uncertainty fit in dark sans style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/xy_uncertainties_sans_plot_light.png" alt="XY uncertainty fit in sans style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/xy_uncertainties_sans_plot_dark.png" alt="XY uncertainty fit in dark sans style without result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/xy_uncertainties_tex_panel_light.png" alt="XY uncertainty fit in tex style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/xy_uncertainties_tex_plot_light.png" alt="XY uncertainty fit in tex style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/xy_uncertainties_tex_panel_dark.png" alt="XY uncertainty fit in dark tex style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-xy" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/xy_uncertainties_tex_plot_dark.png" alt="XY uncertainty fit in dark tex style without result panel">
<div>
<span class="jufitter-tag">effective variance</span>
<span class="jufitter-tag">x errors</span>
<h3><a href="gallery/xy_uncertainties.html">XY uncertainties</a></h3>
<p>Include uncertainty in the independent variable through the local model slope. This example is useful when calibration, frequency, voltage, or position errors are not negligible.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/full_covariance_decay_sans_panel_light.png" alt="Full covariance fit in sans style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/full_covariance_decay_sans_panel_dark.png" alt="Full covariance fit in dark sans style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/full_covariance_decay_sans_plot_light.png" alt="Full covariance fit in sans style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/full_covariance_decay_sans_plot_dark.png" alt="Full covariance fit in dark sans style without result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/full_covariance_decay_tex_panel_light.png" alt="Full covariance fit in tex style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/full_covariance_decay_tex_plot_light.png" alt="Full covariance fit in tex style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/full_covariance_decay_tex_panel_dark.png" alt="Full covariance fit in dark tex style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-covariance" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/full_covariance_decay_tex_plot_dark.png" alt="Full covariance fit in dark tex style without result panel">
<div>
<span class="jufitter-tag">covariance</span>
<span class="jufitter-tag">correlations</span>
<h3><a href="gallery/full_covariance.html">Full covariance</a></h3>
<p>Use a dense covariance matrix when measurements share readout noise. The point is not syntax; it is how correlations change parameter uncertainty and goodness-of-fit interpretation.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/damped_oscillator_decay_sans_panel_light.png" alt="Damped oscillator fit in sans style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/damped_oscillator_decay_sans_panel_dark.png" alt="Damped oscillator fit in dark sans style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/damped_oscillator_decay_sans_plot_light.png" alt="Damped oscillator fit in sans style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/damped_oscillator_decay_sans_plot_dark.png" alt="Damped oscillator fit in dark sans style without result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/damped_oscillator_decay_tex_panel_light.png" alt="Damped oscillator fit in tex style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/damped_oscillator_decay_tex_plot_light.png" alt="Damped oscillator fit in tex style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/damped_oscillator_decay_tex_panel_dark.png" alt="Damped oscillator fit in dark tex style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-damped" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/damped_oscillator_decay_tex_plot_dark.png" alt="Damped oscillator fit in dark tex style without result panel">
<div>
<span class="jufitter-tag">real data</span>
<span class="jufitter-tag">nonlinear</span>
<span class="jufitter-tag">model criticism</span>
<h3><a href="gallery/resonance_decay.html">Damped oscillator</a></h3>
<p>Discover why a visually convincing constant-frequency fit is statistically rejected, test a frequency-drift extension, and use pull structure to decide what must be investigated next.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/photoelectric_threshold_sans_panel_light.png" alt="Photoelectric work-function fit in sans style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/photoelectric_threshold_sans_panel_dark.png" alt="Photoelectric work-function fit in dark sans style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/photoelectric_threshold_sans_plot_light.png" alt="Photoelectric work-function fit in sans style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/photoelectric_threshold_sans_plot_dark.png" alt="Photoelectric work-function fit in dark sans style without result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/photoelectric_threshold_tex_panel_light.png" alt="Photoelectric work-function fit in tex style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/photoelectric_threshold_tex_plot_light.png" alt="Photoelectric work-function fit in tex style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/photoelectric_threshold_tex_panel_dark.png" alt="Photoelectric work-function fit in dark tex style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-photoelectric" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/photoelectric_threshold_tex_plot_dark.png" alt="Photoelectric work-function fit in dark tex style without result panel">
<div>
<span class="jufitter-tag">x/y errors</span>
<span class="jufitter-tag">line intersection</span>
<h3><a href="gallery/photoelectric_threshold.html">Photoelectric work function</a></h3>
<p>Fit baseline and emission regimes separately, then propagate both covariance matrices into the threshold intersection and work function.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/poisson_counts_sans_panel_light.png" alt="Poisson count fit in sans style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/poisson_counts_sans_panel_dark.png" alt="Poisson count fit in dark sans style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/poisson_counts_sans_plot_light.png" alt="Poisson count fit in sans style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/poisson_counts_sans_plot_dark.png" alt="Poisson count fit in dark sans style without result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/poisson_counts_tex_panel_light.png" alt="Poisson count fit in tex style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/poisson_counts_tex_plot_light.png" alt="Poisson count fit in tex style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/poisson_counts_tex_panel_dark.png" alt="Poisson count fit in dark tex style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-poisson" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/poisson_counts_tex_plot_dark.png" alt="Poisson count fit in dark tex style without result panel">
<div>
<span class="jufitter-tag">likelihood</span>
<span class="jufitter-tag">counts</span>
<h3><a href="gallery/poisson_histogram.html">Poisson and histograms</a></h3>
<p>Extract a radioactive half-life and a detector peak from sparse counts. Exact count semantics, integrated unequal bins, empty bins, and deviance residuals replace invented Gaussian error bars.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/constraints_priors_sans_panel_light.png" alt="Constrained fit in sans style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/constraints_priors_sans_panel_dark.png" alt="Constrained fit in dark sans style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/constraints_priors_sans_plot_light.png" alt="Constrained fit in sans style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/constraints_priors_sans_plot_dark.png" alt="Constrained fit in dark sans style without result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/constraints_priors_tex_panel_light.png" alt="Constrained fit in tex style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/constraints_priors_tex_plot_light.png" alt="Constrained fit in tex style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/constraints_priors_tex_panel_dark.png" alt="Constrained fit in dark tex style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-constraints" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/constraints_priors_tex_plot_dark.png" alt="Constrained fit in dark tex style without result panel">
<div>
<span class="jufitter-tag">constraints</span>
<span class="jufitter-tag">profiles</span>
<h3><a href="gallery/constraints_profiles.html">Constraints and profiles</a></h3>
<p>An early saturation measurement leaves amplitude and time constant nonlinearly coupled. Profiles and two-parameter regions show exactly why the local covariance summary fails.</p>
</div>
</div>
<div class="jufitter-gallery-item">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/multi_dataset_shared_slope_sans_panel_light.png" alt="Multi-dataset fit in sans style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="sans" data-jufitter-plot-panel="show" src="assets/gallery/multi_dataset_shared_slope_sans_panel_dark.png" alt="Multi-dataset fit in dark sans style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/multi_dataset_shared_slope_sans_plot_light.png" alt="Multi-dataset fit in sans style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="sans" data-jufitter-plot-panel="hide" src="assets/gallery/multi_dataset_shared_slope_sans_plot_dark.png" alt="Multi-dataset fit in dark sans style without result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/multi_dataset_shared_slope_tex_panel_light.png" alt="Multi-dataset fit in tex style with result panel">
<img class="jufitter-plot-light" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/multi_dataset_shared_slope_tex_plot_light.png" alt="Multi-dataset fit in tex style without result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="tex" data-jufitter-plot-panel="show" src="assets/gallery/multi_dataset_shared_slope_tex_panel_dark.png" alt="Multi-dataset fit in dark tex style with result panel">
<img class="jufitter-plot-dark" data-jufitter-plot-group="gallery-multi" data-jufitter-plot-style="tex" data-jufitter-plot-panel="hide" src="assets/gallery/multi_dataset_shared_slope_tex_plot_dark.png" alt="Multi-dataset fit in dark tex style without result panel">
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

## Run An Example

The page code is arranged in execution order: load the measurements, define the
model, fit, inspect diagnostics, and construct the figure. It can be copied into
a Julia script or a Pluto notebook without depending on hidden setup cells.

If you are working from the JuFitter source tree, the matching scripts live in
`examples/gallery/`. For example:

```bash
julia --project=docs examples/gallery/08_damped_oscillator_decay.jl
```

The scripts use the same data files and statistical assumptions as the pages.
They are useful when you want to modify an analysis end to end; the code blocks
inside each page are better when you want to understand one step at a time.
