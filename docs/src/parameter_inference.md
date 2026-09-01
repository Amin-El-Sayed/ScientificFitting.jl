# Parameters And Fit Quality

A fit often combines the primary dataset with external calibration information,
fixed values, or physical bounds. Those choices change both the objective and
the interpretation of uncertainty. This chapter separates those mechanisms,
then derives degrees of freedom, goodness-of-fit statistics, and the local
parameter covariance reported at the optimum.

## External Parameter Information

### Gaussian Parameter Terms

A symmetric Gaussian term centered at ``\mu_j`` with scale ``\tau_j`` adds

```math
\left(\frac{\theta_j-\mu_j}{\tau_j}\right)^2
```

to chi-square. In the normalized ``-2\log L`` convention it contributes

```math
\log(2\pi\tau_j^2)
+
\left(\frac{\theta_j-\mu_j}{\tau_j}\right)^2.
```

`parameter_priors` implements these scalar penalty terms, including a
piecewise scale for asymmetric quoted errors. `parameter_constraints`
implements the correlated multivariate analogue.

For example, an auxiliary calibration ``g=1.00\pm0.05`` contributes

```math
\left(\frac{1.10-1.00}{0.05}\right)^2=4
```

when the joint fit proposes ``g=1.10``. The primary data can still move the
estimate away from the calibration, but they must pay that likelihood cost.
The ``0.05`` is not appended to the final error afterward: ``g`` is fitted
jointly, so its correlation with every scientific parameter propagates through
the fit. If several calibration quantities share a reference, use one
correlated parameter constraint rather than independent scalar terms that count
the common information more than once.

For asymmetric scales ``\tau_-`` and ``\tau_+``, JuFitter uses the continuous,
normalized split-normal cost

```math
C_\mathrm{split}(\theta_j)
=
\log\frac{\pi}{2}
+2\log(\tau_-+\tau_+)
+
\left(\frac{\theta_j-\mu_j}{\tau_{\pm}}\right)^2,
```

where ``\tau_-`` applies below ``\mu_j`` and ``\tau_+`` above it. The shared
normalization is essential: using a different Gaussian normalization on each
side would make the cost discontinuous at the center.

There are two legitimate interpretations, and they should not be mixed
silently:

- **Auxiliary measurement:** a previous calibration measured the parameter.
  The term is part of a joint frequentist likelihood.
- **Prior information:** the term expresses prior belief. The optimum is then a
  penalized or MAP-like estimate, not a pure maximum-likelihood estimate.

JuFitter evaluates the same numerical term in either case; the scientific
interpretation belongs in the analysis. JuFitter does not turn that term into a
full Bayesian posterior or produce Bayesian posterior intervals.

### Fixed Parameters And Bounds

A fixed parameter is removed from the free optimization variables. It is not a
Gaussian term with an extremely small uncertainty. An optional uncertainty on
a `FixedParameter` is report metadata and does not propagate automatically into
the fitted covariance.

A bound defines the allowed parameter space but contributes no smooth penalty.
Local Hessian errors become unreliable near an active bound because the cost
surface is truncated. Use a profile interval and report the bound when it
affects the result.

Nonlinear equality and inequality constraints similarly change the accessible
parameter geometry. They can be scientifically necessary, but standard
unconstrained asymptotic formulas need not remain exact at their boundary.

## Degrees Of Freedom

For a regular Gaussian fit with ``n`` observations, known full-rank covariance,
and ``k`` free parameters,

```math
\mathrm{ndf}=n-k.
```

Fixed parameters do not count toward ``k``. JuFitter treats scalar Gaussian
parameter terms as one auxiliary observation and a correlated constraint on
``q`` parameters as ``q`` auxiliary observations. This is the natural counting
when those terms represent independent calibration measurements. If they are
subjective priors, interpreting the resulting `ndf` as a frequentist degree of
freedom is not justified.

The count ``n-k`` does not by itself prove a chi-square reference distribution.
For a nonlinear model it is the conventional dimension count; the calibration
of ``\chi^2_{\min}`` is generally asymptotic. Rank-deficient covariance,
parameters on boundaries, or a variance model estimated from the same residuals
require a separate derivation or simulation.

For likelihood fits with an explicit goodness-of-fit statistic, Gaussian
auxiliary terms contribute both their quadratic residuals to that statistic and
their dimensions to `ndf`. Counting the calibration observation without its
discrepancy, or vice versa, would produce an internally inconsistent p-value.

When ``\mathrm{ndf}\le 0``, reduced statistics and chi-square p-values are not
meaningful. JuFitter returns `NaN` for them and reports the problem rather than
inventing a number.

## Goodness Of Fit

For a linear Gaussian model with known, full-rank covariance,

```math
\chi^2_{\min} \sim \chi^2_{\mathrm{ndf}}.
```

This result is exact under the linear-model assumptions. For a regular
nonlinear Gaussian model it is an asymptotic approximation. It need not hold
when the covariance was tuned from the same residuals, a bound is active, the
model is weakly identified, or data-dependent filtering changed the sampling
process.

Its mean and standard deviation are

```math
E[\chi^2]=\mathrm{ndf},
\qquad
\operatorname{sd}(\chi^2)=\sqrt{2\,\mathrm{ndf}}.
```

With 28 degrees of freedom, values fluctuate naturally on a scale
``\sqrt{56}=7.5``. The corresponding reduced chi-square has a standard
deviation of about ``\sqrt{2/28}=0.27``. A universal rule such as
``0.9<\chi^2/\mathrm{ndf}<1.1`` would therefore be absurdly strict for a small
dataset and too weak for a very large one.

JuFitter reports the upper-tail p-value

```math
p=P\!\left(\chi^2_{\mathrm{ndf}}\ge\chi^2_\mathrm{observed}\right).
```

It is the probability of obtaining an equal or larger discrepancy under the
stated model and uncertainty assumptions. It is not:

- the probability that the model is true;
- the probability that a parameter lies in an interval;
- a measure of scientific importance;
- proof that residuals are structure-free.

A small p-value can result from a wrong model, underestimated uncertainties,
missing correlation, outliers, or optimizer failure. A very large p-value can
result from overestimated uncertainties, duplicated or correlated data counted
as independent, or data that have been smoothed before fitting. Always inspect
pulls or residuals alongside the scalar test.

## Local Parameter Covariance

The local error formula is a Taylor approximation, not an extra fit. For one
parameter and a smooth cost ``C=-2\log L``, expand around the minimum
``\hat\theta``:

```math
C(\hat\theta+\delta)
=
C(\hat\theta)
+C'(\hat\theta)\delta
+\frac{1}{2}C''(\hat\theta)\delta^2
+O(\delta^3).
```

At an interior minimum ``C'(\hat\theta)=0``. Matching the remaining quadratic
term to ``\Delta C\approx(\delta/\sigma_\theta)^2`` gives

```math
\sigma_\theta^2
\approx
\frac{2}{C''(\hat\theta)}.
```

For several Gaussian-fit parameters, let ``W`` whiten the observations and
define the weighted model Jacobian

```math
(J_w)_{ij}
=
\left[W\frac{\partial m}{\partial\theta_j}\right]_i.
```

Near the solution, the whitened residual vector obeys

```math
z(\hat\theta+\delta)
\approx
z(\hat\theta)-J_w\delta.
```

Substitution into ``\chi^2=z^Tz`` yields

```math
\Delta\chi^2
\approx
-2z(\hat\theta)^T J_w\delta
+\delta^T J_w^T J_w\delta.
```

The linear term vanishes at the least-squares optimum. Comparing the remaining
quadratic form with
``\Delta\chi^2\approx\delta^T\operatorname{Cov}(\hat\theta)^{-1}\delta``
gives

```math
\operatorname{Cov}(\hat\theta)
\approx
(J_w^T J_w)^{-1}.
```

For a general cost ``C=-2\log L`` with Hessian

```math
H=\nabla^2 C(\hat\theta),
```

the local covariance is

```math
\operatorname{Cov}(\hat\theta)\approx 2H^{-1}.
```

The factor two follows from the ``-2\log L`` convention. These formulas describe
the same local curvature in two equivalent forms: ``J_w^TJ_w`` for a locally
linear Gaussian residual problem, and the full objective Hessian otherwise.
They describe the curvature at one point. They are reliable when the estimator is well
identified, the minimum is interior, and the cost is approximately quadratic
over the reported uncertainty region.

They can fail for weak data, highly nonlinear parameterizations, degeneracies,
multiple minima, active bounds, or asymmetric likelihoods. A negative local
curvature is not an uncertainty. JuFitter marks materially invalid covariance
geometry as critical and returns `NaN` rather than fabricating a zero standard
error.

### Covariance Scaling

When measurement uncertainties are known externally, multiplying the parameter
covariance by ``\chi^2/\mathrm{ndf}`` changes their stated meaning and should not
be automatic. If instead the residual scale is unknown and estimated from the
same data, the classical least-squares estimate includes that factor.

JuFitter makes the policy explicit through
`scale_covariance=:auto | :never | :always`. Report which interpretation was
used; do not rescale merely to make reduced chi-square look closer to one.

Local covariance is deliberately a local summary. Continue with
[Profiles and Contours](profiles_contours.md) when bounds, nonlinearity,
asymmetry, or strong correlation make the quadratic approximation doubtful.
For a complete nonlinear workflow, see
[Constraints and Profiles](gallery/constraints_profiles.md).
