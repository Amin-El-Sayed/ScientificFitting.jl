# ScientificFitting

```@raw html
<section class="scientificfitting-hero">
  <div class="scientificfitting-kicker">Least squares and likelihood fits in Julia</div>
  <p class="scientificfitting-lede">
    Fit your model to measurements. Inspect the result, not just the curve.
  </p>
  <p>
    Supply standard deviations or covariances for Gaussian measurement errors,
    use Poisson counts, or define your own likelihood. Fit with parameter
    constraints, inspect residuals and profiles, and edit the resulting Makie
    figure without rewriting the fit.
  </p>
  <div class="scientificfitting-hero-actions">
    <a class="scientificfitting-button primary" href="quickstart.html">Start with a complete fit</a>
    <a class="scientificfitting-button" href="install.html">Install</a>
  </div>
</section>
```

**What do the uncertainties mean?** Measurement errors describe the assumed
observation distribution. Reported parameter errors are local covariance
approximations; profile-likelihood scans show where that approximation breaks
down. A band for the fitted model does not include the scatter of future
measurements unless explicitly stated. ScientificFitting uses optimization,
not posterior sampling. The [statistical chapters](statistics.md)
explain these distinctions and the assumptions behind interval coverage.

Start with [Linear Calibration](quickstart.md) for the shortest
complete analysis, or choose the example closest to your data below. Every page
contains the measurements, assumptions, executable code, actual program output,
diagnostics, and scientific interpretation.

## Worked Examples

New for v0.3: [LHCb open-data mass spectrum](gallery/lhcb_mass_spectrum.md).
Compose a named signal/background model, fit with NativeMinuit, and check yield
sensitivity to the peak shape. Includes source DOI and reproducible ROOT selection.

```@raw html
<div class="scientificfitting-gallery-grid">
<div class="scientificfitting-gallery-item">
<img src="assets/gallery/xy_uncertainties_sans_panel_light.png" alt="XY uncertainty fit in sans style with result panel">
<div>
<span class="scientificfitting-tag">effective variance</span>
<span class="scientificfitting-tag">x errors</span>
<h3><a href="gallery/xy_uncertainties.html">XY uncertainties</a></h3>
<p>Include uncertainty in the independent variable through the local model slope, when calibration, frequency, voltage, or position errors are not negligible.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img src="assets/gallery/full_covariance_decay_sans_panel_light.png" alt="Full covariance fit in sans style with result panel">
<div>
<span class="scientificfitting-tag">covariance</span>
<span class="scientificfitting-tag">correlations</span>
<h3><a href="gallery/full_covariance.html">Full covariance</a></h3>
<p>Use a dense covariance matrix when measurements share readout noise, and see how correlations change parameter uncertainty and goodness-of-fit interpretation.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img src="assets/gallery/damped_oscillator_decay_sans_panel_light.png" alt="Damped oscillator fit in sans style with result panel">
<div>
<span class="scientificfitting-tag">real data</span>
<span class="scientificfitting-tag">nonlinear</span>
<span class="scientificfitting-tag">model criticism</span>
<h3><a href="gallery/resonance_decay.html">Damped oscillator</a></h3>
<p>Discover why a visually convincing constant-frequency fit is statistically rejected, test a frequency-drift extension, and use pull structure to decide what must be investigated next.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img src="assets/gallery/photoelectric_threshold_sans_panel_light.png" alt="Photoelectric work-function fit in sans style with result panel">
<div>
<span class="scientificfitting-tag">x/y errors</span>
<span class="scientificfitting-tag">line intersection</span>
<h3><a href="gallery/photoelectric_threshold.html">Photoelectric work function</a></h3>
<p>Fit baseline and emission regimes separately, then propagate both covariance matrices into the threshold intersection and work function.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img src="assets/gallery/poisson_counts_sans_panel_light.png" alt="Poisson count fit in sans style with result panel">
<div>
<span class="scientificfitting-tag">likelihood</span>
<span class="scientificfitting-tag">counts</span>
<h3><a href="gallery/poisson_histogram.html">Poisson and histograms</a></h3>
<p>Extract a radioactive half-life and a detector peak from sparse counts. Exact count semantics, integrated unequal bins, empty bins, and deviance residuals replace invented Gaussian error bars.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img src="assets/gallery/constraints_priors_sans_panel_light.png" alt="Constrained fit in sans style with result panel">
<div>
<span class="scientificfitting-tag">constraints</span>
<span class="scientificfitting-tag">profiles</span>
<h3><a href="gallery/constraints_profiles.html">Constraints and profiles</a></h3>
<p>An early saturation measurement leaves amplitude and time constant nonlinearly coupled. Profiles and two-parameter regions show why the local covariance summary fails.</p>
</div>
</div>
<div class="scientificfitting-gallery-item">
<img src="assets/gallery/multi_dataset_shared_slope_sans_panel_light.png" alt="Multi-dataset fit in sans style with result panel">
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
For derivations, use [Statistical Foundations](statistics.md); for exact signatures and
defaults, use the [API Reference](@ref).
