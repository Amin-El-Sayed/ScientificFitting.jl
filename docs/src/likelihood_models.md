# Likelihoods And Model Comparison

*Chapter 5 of 5 · Previous: [Profiles and Contours](profiles_contours.md)*

Counts, histogram bins, and individual events do not automatically carry
Gaussian error bars. Their sampling process defines a likelihood. This chapter
derives ScientificFitting's Poisson, histogram, unbinned, and extended likelihood costs,
then separates parameter estimation, goodness of fit, and model comparison.

## User-Defined Measurement Errors

For independent observations, specify a normalized density or mass
``q_i(y_i\mid f(x_i,p),p)``. The common likelihood cost is

```math
C(p)=-2\sum_i\log q_i(y_i\mid f(x_i,p),p).
```

`fit_likelihood_model` takes these log probabilities in one vectorized callback.
The distribution can differ between observations; known scales or trial counts
can be captured by the callback. Dependent observations instead need a **joint**
likelihood via `fit_custom`; multiplying their marginal densities is not enough.

For example, Laplace errors of known scale ``b`` have density
``q(y\mid\mu)=\exp(-|y-\mu|/b)/(2b)``. Thus
``C(\mu)=2\sum_i|y_i-\mu|/b+2n\log(2b)``: the optimum is a sample median,
not a mean, and the observed cost has corners rather than quadratic curvature.

```@example laplace_measurements
using ScientificFitting

y = [-1.2, -0.1, 0.2, 0.4, 0.8, 1.3, 5.0]
location(x, p) = fill(p[1], length(x))
# Known scale b = 1; retain the normalizing constant.
logprob(y, mu, p) = -abs.(y .- mu) .- log(2.0)
result = fit_likelihood_model(location, collect(eachindex(y)), y;
    logprob, p0=[0.1], optimizer=:nelder_mead)
(location=round(only(result.params); digits=6), local_error=only(result.param_stderr))
```

Nelder-Mead does not differentiate this cost. Its default
`parameter_covariance=:none` leaves the local error unavailable, not zero.
Profiles can still evaluate the actual cost over explicit `values` grids.
This makes the optimization usable without pretending that every distribution
has Hessian-based errors or the same confidence thresholds.

### A Moving Support Boundary

Suppose a trigger time ``\mu`` is followed by independent exponential delays of
known scale ``b``: ``Y_i=\mu+E_i``, ``E_i\sim\operatorname{Exp}(\text{scale}=b)``.
The likelihood is zero if ``\mu>\min_i y_i``. Within its support,

```math
\hat\mu=\min_i y_i,\qquad
\Delta C(\mu)=\frac{2n}{b}(\hat\mu-\mu),\quad\mu\le\hat\mu.
```

There is no parabolic minimum. A derivative-free optimizer can find the edge,
but changing optimizers cannot make a Hessian interval meaningful. Here an
exact calculation is simple: ``\hat\mu-\mu`` is exponential with rate ``n/b``.
An interval ``[\hat\mu-db/(2n),\hat\mu]`` therefore covers the true ``\mu`` with
probability ``1-e^{-d/2}``. The usual one-parameter threshold ``d=1`` covers only
**39.3%**, not 68.3%. For 68.3% coverage use ``d=-2\log(1-0.683)`` instead.
This is a property of the sampling model, not a numerical defect. For more
complicated non-regular likelihoods, interval calibration needs an appropriate
sampling calculation or simulation; ScientificFitting does not invent one.

## Poisson Counts And Histograms

For an observed non-negative integer count ``n_i`` with expected count
``\mu_i(\theta)>0``, the sampling model is

```math
P(n_i\mid\mu_i)
=
e^{-\mu_i}\frac{\mu_i^{n_i}}{n_i!}.
```

ScientificFitting minimizes the normalized cost

```math
C(\theta)
=
2\sum_i
\left[
\mu_i(\theta)-n_i\log\mu_i(\theta)+\log\Gamma(n_i+1)
\right].
```

The factorial term is constant for optimization but keeps likelihood values and
information criteria on a defined scale.

A zero-count bin is not a zero-uncertainty measurement. If ``n=0`` and
``\mu=0.5``, the bin contributes a finite likelihood cost and a Poisson deviance
of ``1``. The common Gaussian shortcut ``\sigma=\sqrt n`` would instead assign
zero uncertainty and break exactly where the Poisson model is most needed.

For histogram fits, the model must predict the expected count in each bin. A
density should be integrated across the bin:

```math
\mu_i(\theta)
=
N\int_{b_i}^{b_{i+1}} f(x\mid\theta)\,dx.
```

Evaluating a curved density only at the bin center can bias wide or uneven
bins. `fit_histogram_density` performs the bin integrations; use
`fit_histogram_model` when the supplied model already returns expected bin
counts.

### Poisson Deviance

The likelihood-ratio goodness-of-fit statistic against a saturated count model
is

```math
D(\theta)
=
2\sum_i
\left[
\mu_i-n_i+n_i\log\frac{n_i}{\mu_i}
\right],
```

with the logarithmic term defined as zero for ``n_i=0``. ScientificFitting stores this
deviance in the chi-square-like statistics fields for Poisson and histogram
fits.

Its chi-square calibration is asymptotic. With many very small expected counts,
the reported p-value can be inaccurate even though the Poisson likelihood used
for parameter estimation is correct. Combine sparse bins only when
scientifically defensible, or calibrate the deviance distribution with
simulation.

## Unbinned And Extended Likelihoods

For independent observations ``x_1,\ldots,x_n`` from a normalized density
``f(x\mid\theta)``, the unbinned cost is

```math
C(\theta)
=
-2\sum_{i=1}^{n}\log f(x_i\mid\theta).
```

The density must be positive and normalized over the sampling domain. An
unbinned fit preserves positions within bins, but it is not automatically more
correct: detector acceptance, truncation, resolution, and selection must still
appear in ``f``.

If the expected event count also carries information, describe the data as a
Poisson point process with intensity ``\lambda(x\mid\theta)`` and integrated
rate

```math
\Lambda(\theta)=\int_\Omega \lambda(x\mid\theta)\,dx.
```

The extended cost is

```math
C_\mathrm{ext}(\theta)
=
2\Lambda(\theta)
-2\sum_{i=1}^{n}\log\lambda(x_i\mid\theta).
```

`fit_unbinned_model` uses a normalized density. `fit_extended_unbinned_model`
uses an event intensity and numerically integrates it over the declared domain.

### When A P-Value Is Not Available

An unbinned likelihood value is not, by itself, a universal goodness-of-fit
statistic. Its absolute scale depends on the density and measurement units.
ScientificFitting therefore returns `NaN` rather than inventing chi-square, reduced
chi-square, or a p-value for generic unbinned and extended fits.

Goodness of fit then needs a separate statistic appropriate to the scientific
question: an empirical-CDF test in one dimension,
probability-integral-transform diagnostics, multidimensional simulation checks,
or a likelihood-ratio test against a specified alternative. Parameter
estimation and goodness of fit are related tasks, not the same calculation.

## Model Comparison With AIC And BIC

For ``k`` fitted parameters and maximized normalized likelihood ``L_{\max}``,

```math
\mathrm{AIC}=2k-2\log L_{\max},
```

```math
\mathrm{BIC}=k\log n-2\log L_{\max}.
```

Smaller is preferred within a valid comparison. AIC estimates relative
out-of-sample predictive information loss under its regularity assumptions.
BIC is a large-sample approximation with a stronger complexity penalty and
additional assumptions. Neither is a p-value, proof of a model, or a substitute
for residual inspection.

ScientificFitting uses its stored observation count `nobs` in the BIC arithmetic;
Gaussian auxiliary dimensions are included in that count. For strongly
correlated data there may be no unique iid-like effective sample size, so a BIC
number can be algebraically defined without having the usual model-selection
interpretation. State what ``n`` means before using BIC as evidence.

Only compare values when all candidates use:

- the same observations and event-selection domain;
- the same likelihood normalization and uncertainty model;
- the same treatment of auxiliary parameter information;
- comparable definitions of ``n`` and free-parameter count ``k``.

Do not compare AIC from a Gaussian fit to AIC from an unnormalized custom cost,
or models fitted to different subsets of data. When Gaussian parameter terms
are interpreted as priors rather than auxiliary observations, ordinary AIC/BIC
no longer have their standard maximum-likelihood interpretation. ScientificFitting still
reports the arithmetic values from its normalized joint cost; the analyst must
not overstate them.

A lower criterion means only that one candidate is preferred relative to the
others considered. A physically missing model is not rescued because it has the
best AIC in a poor candidate set.

See [Poisson and Histograms](gallery/poisson_histogram.md) for two executable
count workflows and [Damped Oscillator](gallery/resonance_decay.md) for a model
comparison whose conclusion is supported by residual structure and a nested
cost difference, not AIC alone.

---

**Previous chapter:** [Profiles and Contours](profiles_contours.md) · **End of
the sequence.** Return to [Statistical Foundations](statistical_foundations.md)
for the common reporting contract, or apply the methods in the
[Gallery](gallery.md).
