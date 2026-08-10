# Statistical Foundations

Fitting is not primarily an optimization problem. The optimizer finds a
minimum; the probability model determines what that minimum means. This section
derives the statistical quantities used by JuFitter, states their assumptions,
and explains what must be checked before a result is reported.

## Choose A Reading Path

The chapters are organized by scientific question rather than package type.
You do not need to read all of them for every fit.

| Your data or question | Start here | Continue with |
| --- | --- | --- |
| measured values with known standard uncertainties | [Gaussian Fits and Covariance](gaussian_models.md) | [goodness of fit](parameter_inference.md#Goodness-Of-Fit) |
| correlated measurements or long time series | [correlated measurements](gaussian_models.md#Correlated-Measurements-And-Whitening) | [structured whitening](gaussian_models.md#Structured-Whitening) |
| uncertainty in both x and y | [uncertainty in x](gaussian_models.md#Uncertainty-In-X) | [parameter-dependent covariance](gaussian_models.md#Parameter-Dependent-Covariance) |
| external calibration, bounds, or fixed parameters | [Parameters and Fit Quality](parameter_inference.md) | [local covariance](parameter_inference.md#Local-Parameter-Covariance) |
| nonlinear or weakly constrained parameters | [Profiles and Contours](profiles_contours.md) | [Constraints and Profiles](gallery/constraints_profiles.md) |
| detector counts or histogram bins | [Likelihoods and Model Comparison](likelihood_models.md) | [Poisson deviance](likelihood_models.md#Poisson-Deviance) |
| individual, unbinned observations | [unbinned likelihoods](likelihood_models.md#Unbinned-And-Extended-Likelihoods) | [goodness-of-fit limits](likelihood_models.md#When-A-P-Value-Is-Not-Available) |
| competing models | [AIC and BIC](likelihood_models.md#Model-Comparison-With-AIC-And-BIC) | residual and scientific checks |

For operational advice, use [Fitting for Practitioners](fitting_for_practitioners.md).
For complete executable analyses, use the [Gallery](gallery.md).

## The Statistical Contract

Every fit combines five distinct objects:

```@raw html
<div class="jufitter-flow">
  <div class="jufitter-flow-step"><strong>Observations</strong><span>What was recorded, and under which sampling process?</span></div>
  <div class="jufitter-flow-step"><strong>Model</strong><span>What mean, density, or event rate does physics predict?</span></div>
  <div class="jufitter-flow-step"><strong>Uncertainty</strong><span>Which fluctuations are independent, correlated, or external?</span></div>
  <div class="jufitter-flow-step"><strong>Cost</strong><span>Which likelihood or quadratic form follows from those assumptions?</span></div>
  <div class="jufitter-flow-step"><strong>Inference</strong><span>Which estimates, intervals, tests, and diagnostics are justified?</span></div>
</div>
```

Changing the optimizer should not change the statistical model. Changing
``\sigma_y`` to a covariance matrix, or a Gaussian cost to a Poisson
likelihood, does.

JuFitter uses the convention

```math
C(\theta) = -2\log L(\theta)
```

for normalized likelihood costs. With this convention, likelihood-ratio
differences, Gaussian chi-square, local curvature, and profile thresholds share
the same scale.

For a static Gaussian covariance, the optimizer may minimize only ``\chi^2``
because the omitted normalization is constant in the parameters. JuFitter keeps
the two quantities separate: `stats.cost_min` is the objective that was
minimized, while `stats.minus2loglik_min` also includes the Gaussian
normalization ``n\log(2\pi)+\log\det V`` and normalized auxiliary terms. They
have the same minimizer only while the covariance is parameter-independent.

## Data, Model, And Residuals

For observations ``d`` and model predictions ``m(\theta)``, the raw residual is

```math
r(\theta)=d-m(\theta).
```

A residual is not yet statistically meaningful. A residual of ``0.2`` may be
tiny for a voltmeter with ``\sigma=1\,\mathrm{V}`` and enormous for one with
``\sigma=1\,\mathrm{mV}``. The uncertainty model supplies that scale.

For independent measurements, the standardized residual or pull is

```math
z_i=\frac{d_i-m_i(\theta)}{\sigma_i}.
```

If the model and uncertainties are correct, these pulls should fluctuate around
zero with a scale near one and no visible structure. At the fitted parameters
they are not independent ``\mathcal N(0,1)`` draws: estimating parameters
projects out fitted directions, and pointwise leverage changes their variance.
A pull plot is therefore a localization tool, not a second chi-square test. A
small total cost cannot replace the pattern check: alternating residuals, long
same-sign runs, or a frequency-dependent drift can reveal model failure even
when one summary number looks acceptable.

## What JuFitter Reports

The result fields follow these conventions:

| field | meaning |
| --- | --- |
| `stats.cost_min` | minimized selected objective: usually ``\chi^2`` for static Gaussian least squares and ``-2\log L`` for likelihood wrappers |
| `stats.minus2loglik_min` | normalized Gaussian ``-2\log L`` for an x-y `FitResult`; for likelihood/custom results it equals the stored objective and has likelihood meaning only when that objective follows the documented normalization |
| `stats.chi2` | Gaussian chi-square or Poisson deviance; `NaN` when no generic statistic exists |
| `stats.chi2_ndf` | chi-square-like statistic divided by positive ndf |
| `stats.pvalue` | upper-tail chi-square approximation where defined |
| `stats.aic`, `stats.bic` | information criteria from the normalized cost and free-parameter count |
| `param_covariance` | local Jacobian/Hessian approximation |
| `param_stderr` | square roots of valid covariance diagonal entries |

Use `diagnose(result)` or `diagnostic_dashboard(result)` before reporting. Use
`profile`, `profile_interval`, `contour`, or `profile_matrix` when local
quadratic errors are scientifically doubtful.

## Reporting Checklist

A reproducible result should state:

1. The observations, units, and selection or binning rules.
2. The model and physical meaning of every fitted parameter.
3. The uncertainty or likelihood model, including correlations and auxiliary
   parameter information.
4. The minimized cost, number of free parameters, and goodness-of-fit statistic
   only where its reference distribution is justified.
5. Whether covariance errors or profile intervals were reported.
6. Residual, pull, profile, or simulation checks relevant to the model.
7. Active bounds, fixed parameters, failed scans, or approximations that limit
   interpretation.

## Further Reading

- [Particle Data Group: Statistics review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-statistics.pdf) for likelihood, goodness-of-fit, confidence intervals, nuisance parameters, and asymptotic limits.
- [Baker and Cousins, *Clarification of the use of chi-square and likelihood functions in fits to histograms*](https://doi.org/10.1016/0167-5087(84)90016-4) for Poisson likelihood-ratio fitting.
- [Wilks, *The Large-Sample Distribution of the Likelihood Ratio*](https://doi.org/10.1214/aoms/1177732360) for the theorem behind common profile thresholds.
- [Akaike, *A new look at the statistical model identification*](https://doi.org/10.1109/TAC.1974.1100705) for the information criterion and its predictive interpretation.

Continue with [Gaussian Fits and Covariance](gaussian_models.md) for measured
values, [Likelihoods and Model Comparison](likelihood_models.md) for counts or
events, or [Fitting for Practitioners](fitting_for_practitioners.md) for a
decision workflow.
