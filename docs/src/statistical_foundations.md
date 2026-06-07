# Statistical Foundations

This page explains the statistical model behind JuFitter. The practical idea is
simple: the optimizer only finds the minimum. The scientific meaning comes from
the data model, uncertainty model, and cost function you choose.

```@raw html
<div class="jufitter-flow">
  <div class="jufitter-flow-step"><strong>Data</strong><span>Measurements, counts, histograms, or samples.</span></div>
  <div class="jufitter-flow-step"><strong>Model</strong><span>Prediction with named physical parameters.</span></div>
  <div class="jufitter-flow-step"><strong>Uncertainty</strong><span>Standard errors, covariance matrices, x errors, priors.</span></div>
  <div class="jufitter-flow-step"><strong>Cost</strong><span>χ², Gaussian NLL, Poisson NLL, or custom objective.</span></div>
  <div class="jufitter-flow-step"><strong>Optimizer</strong><span>Minimizes the selected cost without defining the statistics.</span></div>
  <div class="jufitter-flow-step"><strong>Diagnostics</strong><span>Parameters, covariance, profiles, contours, residuals, reports.</span></div>
</div>
```

## How To Use This Page

If you are fitting lab or engineering data for the first time, read this page in
three passes:

1. Pick the row in the decision table below that matches your data.
2. Read only the corresponding cost-function section.
3. Come back to profiles and contours when the dashboard warns that local
   covariance errors may not be enough.

If you already know likelihood theory, use the same page as a convention
reference: JuFitter uses the ``-2\log L`` convention, reports local covariance
from curvature, and treats profile/contour scans as the check on that local
approximation.

## Which Statistical Model Should I Use?

The fit is only as meaningful as the probability model behind it. The table
below is the practical starting point.

| data situation | use | why |
| --- | --- | --- |
| y measurements with independent standard uncertainties | `sigma_y` with Gaussian least squares | each residual is measured in its own standard deviations |
| x and y measurements both have relevant uncertainty | `sigma_x` plus `sigma_y` | x uncertainty is propagated through the local model slope |
| points share calibration, baseline, normalization, or readout noise | `cov_y`, `cov_x`, or uncertainty components | off-diagonal covariance prevents shared noise from being counted as independent evidence |
| counts in bins or detector channels | Poisson or histogram likelihood | count variance follows the expected count, especially when counts are small |
| unbinned samples from a distribution | unbinned likelihood | no information is lost by arbitrary binning |
| external knowledge about parameters | priors, fixed parameters, or parameter constraints | external information belongs in parameter space, not as fake data points |
| suspicious local errors, bounds, weak data, nonlinear model | profiles and contours | the likelihood shape, not the Hessian alone, defines credible intervals |

Two common mistakes are worth avoiding:

- Do not use ordinary least squares just because it is familiar. If the data
  are counts, sparse histograms, or strongly correlated measurements, the cost
  function should reflect that.
- Do not treat `param_stderr` as final when the likelihood is visibly
  non-parabolic. It is a local approximation; profiles and contours test when
  that approximation is safe.

## A Minimal Mental Model

For independent Gaussian measurements, each point contributes roughly

```math
\left(\frac{\text{observed}-\text{predicted}}{\text{standard uncertainty}}\right)^2.
```

A residual of one standard deviation contributes about one unit to ``\chi^2``.
Ten independent residuals of typical size one therefore give
``\chi^2 \approx 10``. This is why ``\chi^2/\mathrm{ndf}`` near one is a useful
first sanity check.

Correlations change this mental model. If two points share the same baseline
error, they are not two fully independent pieces of evidence. A covariance
matrix tells the fit which residual patterns are plausible together. That is
why the same visible scatter can imply different parameter errors depending on
whether the points are independent or correlated.

## Data, Model, Residual

For a vector of observations ``d`` and model predictions ``m(\theta)``, the
residual is

```math
r(\theta) = d - m(\theta).
```

The parameter vector ``\theta`` contains only quantities that are estimated by
the fit. Fixed parameters are inserted into the model, but they are not part of
the optimizer-visible free parameter vector.

## Gaussian Least Squares

For independent Gaussian uncertainties ``\sigma_i``, the familiar chi-square is

```math
\chi^2(\theta)
=
\sum_i
\left(
\frac{d_i - m_i(\theta)}{\sigma_i}
\right)^2.
```

This is the right default when every point has a known independent standard
uncertainty and the residuals should be approximately Gaussian.

For correlated Gaussian uncertainties, use the covariance matrix ``V``:

```math
\chi^2(\theta)=r(\theta)^T V^{-1} r(\theta).
```

The off-diagonal entries of ``V`` describe which residual patterns are plausible
together. For two points with identical variance ``\sigma^2`` and correlation
``\rho``,

```math
V=\sigma^2
\begin{pmatrix}
1 & \rho \\
\rho & 1
\end{pmatrix}.
```

Positive ``\rho`` means "both residuals high" is less surprising than one high
and one low. This is exactly what baseline drift, shared calibration constants,
or finite-memory electronics can do. A fully systematic scale error is often
better handled as a nuisance parameter or separate propagated uncertainty than
as a short-range covariance matrix.

JuFitter does not form ``V^{-1}`` explicitly in production calculations. It
uses factorization and linear solves, because those are more stable and faster
than materializing an inverse matrix.

## Full Gaussian Likelihood

If the data are distributed as

```math
d \sim \mathcal{N}(m(\theta), V(\theta)),
```

then JuFitter's Gaussian negative log-likelihood convention is

```math
\mathrm{NLL}(\theta) = -2\log L(\theta)
```

and

```math
\mathrm{NLL}(\theta)
=
n\log(2\pi)
+
\log\det V(\theta)
+
r(\theta)^T V(\theta)^{-1}r(\theta).
```

If ``V`` does not depend on the parameters, the first two terms are constant for
the location of the minimum. Minimizing ``\chi^2`` and minimizing the full
Gaussian NLL then give the same best-fit parameters.

If ``V`` depends on ``\theta``, the ``\log\det V(\theta)`` term is not constant
and must remain in the cost. JuFitter therefore switches `cost=:auto` to the
Gaussian NLL when the effective covariance is parameter-dependent.

## X Uncertainties

For a smooth one-dimensional model ``y=f(x,\theta)``, small x uncertainties can
be propagated into an effective y covariance:

```math
V_\mathrm{eff}(\theta)
=
V_y + J_x(\theta)V_xJ_x(\theta)^T.
```

For pointwise independent x errors this becomes

```math
\sigma_{\mathrm{eff},i}^2(\theta)
=
\sigma_{y,i}^2
+
\left(
\frac{\partial f(x_i,\theta)}{\partial x}
\right)^2
\sigma_{x,i}^2.
```

This first-order approximation is useful for small x uncertainties and smooth
models. If x errors are large, the model is strongly nonlinear, or the
measurement process constrains x and y jointly, a full errors-in-variables model
or orthogonal-distance formulation is statistically cleaner.

## Priors, Constraints, And Fixed Parameters

A Gaussian prior on parameter ``\theta_j`` adds

```math
\left(
\frac{\theta_j-\mu_j}{\tau_j}
\right)^2
```

to a chi-square cost. In the full NLL convention, its normalization is included:

```math
\log(2\pi\tau_j^2)
+
\left(
\frac{\theta_j-\mu_j}{\tau_j}
\right)^2.
```

Fixed parameters are different from tight priors. A fixed parameter is removed
from the optimizer problem entirely. If you attach an uncertainty to a fixed
parameter, that uncertainty is reportable context; it is not a force pulling the
fit toward the fixed value.

Bounds and constraints change the interpretation of local symmetric errors.
When a best-fit parameter is close to a bound, prefer profile intervals over a
single covariance-derived standard error.

## Parameter Covariance

For a locally linear least-squares problem with weighted Jacobian ``J_w``,

```math
\operatorname{Cov}(\hat\theta)
\approx
(J_w^T J_w)^{-1}.
```

If the measurement uncertainties are known from the experiment, this covariance
should not be automatically inflated to force ``\chi^2/\mathrm{ndf}`` to one.
JuFitter exposes this choice explicitly through
`scale_covariance=:auto | :never | :always`.

For general NLL fits, the covariance is a local curvature approximation:

```math
\operatorname{Cov}(\hat\theta)
\approx
2H^{-1},
```

where

```math
H=\nabla^2 \mathrm{NLL}(\hat\theta)
```

in the convention ``\mathrm{NLL}=-2\log L``.

This approximation is useful but not a promise that the likelihood is
parabolic. Nonlinear models, weak data, active bounds, and asymmetric
likelihoods can make the local covariance too optimistic or even misleading.
That is why JuFitter also exposes profile intervals and pairwise contours: they
inspect the cost surface away from the minimum instead of trusting only the
local Hessian.

## Goodness Of Fit

For Gaussian residuals with known uncertainties, the expected scale is

```math
\chi^2 \approx \mathrm{ndf},
```

with natural fluctuations of about

```math
\sqrt{2\,\mathrm{ndf}}.
```

The p-value ``P(\chi^2)`` is the probability, assuming the model and uncertainty
model are correct, of observing a chi-square at least as large as the measured
one. It is not the probability that the model is true.

Very small p-values suggest missing model structure, underestimated
uncertainties, outliers, wrong correlations, or optimizer failure. Very large
p-values can indicate overestimated uncertainties, hidden correlations, or
over-smoothed data.

## Profiles and Contours

Local covariance errors assume the cost is parabolic near the minimum:

```math
\Delta C(\theta_i)
\approx
\left(
\frac{\theta_i-\hat\theta_i}{\sigma_i}
\right)^2.
```

A profile scan checks that assumption. It fixes one parameter, refits all other
free parameters, and records

```math
C_\mathrm{prof}(a)
=
\min_{\theta_{-i}} C(\theta_i=a,\theta_{-i}).
```

The plotted quantity is

```math
\Delta C(a)=C_\mathrm{prof}(a)-C_\min.
```

For two parameters, a contour scan fixes both parameters on a grid and refits
the rest. In the ``-2\log L`` convention, common Wilks thresholds are
``\Delta C=1.0`` for an approximate one-parameter 1σ interval and
``\Delta C=2.30`` for an approximate joint two-parameter 1σ contour.

Profiles and contours are diagnostic tools. They show when symmetric local
errors are reliable and when the likelihood geometry is skewed, clipped,
banana-shaped, or otherwise non-Gaussian.

## Poisson Count Models

For observed counts ``n_i`` and expected counts ``\mu_i(\theta)``,

```math
P(n_i\mid\mu_i)=e^{-\mu_i}\frac{\mu_i^{n_i}}{n_i!}.
```

Ignoring constants that do not depend on the model parameters,

```math
-2\log L(\theta)
=
2\sum_i
\left(
\mu_i(\theta) - n_i\log\mu_i(\theta)
\right).
```

For goodness-of-fit, the Poisson deviance is

```math
D(\theta)
=
2\sum_i
\left[
\mu_i(\theta)-n_i+n_i\log\frac{n_i}{\mu_i(\theta)}
\right],
```

where ``n_i\log(n_i/\mu_i)`` is defined as zero when ``n_i=0``.

Use Poisson or histogram likelihoods for count data. Gaussian least squares can
be a poor approximation when counts are small or bins are sparse.
