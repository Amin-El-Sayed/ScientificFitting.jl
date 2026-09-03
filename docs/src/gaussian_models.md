# Gaussian Fits And Covariance

*Chapter 2 of 5 · Previous:
[Statistical Foundations](statistical_foundations.md) · Next:
[Parameters and Fit Quality](parameter_inference.md)*

Use a Gaussian measurement model when the recorded quantity is continuous and
its uncertainty describes a standard deviation or covariance. This chapter
starts with independent errors, then shows what changes when observations move
together, the covariance is structured, or uncertainty in x must be propagated
into the response.

## Gaussian Least Squares

Start from the measurement equation, not from a least-squares algorithm:

```math
d_i=m_i(\theta)+\epsilon_i,
\qquad
\epsilon_i\sim\mathcal N(0,\sigma_i^2).
```

For independent observations, the probability density of the complete dataset
is the product of the one-point densities,

```math
L(\theta)
=
\prod_{i=1}^{n}
\frac{1}{\sqrt{2\pi}\,\sigma_i}
\exp\!\left[
-\frac{1}{2}
\left(\frac{d_i-m_i(\theta)}{\sigma_i}\right)^2
\right].
```

Taking ``-2\log`` turns the product into a sum:

```math
-2\log L(\theta)
=
\sum_i \log(2\pi\sigma_i^2)
+
\sum_i
\left(\frac{d_i-m_i(\theta)}{\sigma_i}\right)^2.
```

When the quoted ``\sigma_i`` are known and do not depend on ``\theta``, the
first sum is constant. Maximizing the Gaussian likelihood is therefore
equivalent to minimizing

```math
\chi^2(\theta)
=
\sum_{i=1}^{n}
\left(\frac{d_i-m_i(\theta)}{\sigma_i}\right)^2.
```

Each one-sigma residual contributes one unit. For example, pulls
``(1,-1,0.5,-0.5)`` contribute

```math
\chi^2 = 1^2+(-1)^2+0.5^2+(-0.5)^2=2.5.
```

This cost is the right default when:

- the response variable is continuous;
- the quoted uncertainties describe repeated-measurement scatter;
- the Gaussian approximation is reasonable;
- correlations are absent or have already been modeled.

It is not the right default merely because a least-squares solver is
convenient. Sparse counts, censored observations, and strongly non-Gaussian
measurements require a likelihood matching their sampling process.

## Correlated Measurements And Whitening

Let ``V=\operatorname{Cov}(d)`` be the observation covariance. Generalized least
squares uses

```math
\chi^2(\theta)=r(\theta)^T V^{-1}r(\theta).
```

The diagonal entries are variances. Off-diagonal entries describe residual
patterns that tend to move together.

### A Two-Point Example

Take two residuals with equal standard uncertainty ``\sigma`` and correlation
``\rho=0.8``:

```math
V=\sigma^2
\begin{pmatrix}
1 & 0.8\\
0.8 & 1
\end{pmatrix}.
```

Two one-sigma residual patterns then have very different costs:

```math
r_\mathrm{common}=\sigma(1,1)^T,
\qquad
\chi^2_\mathrm{common}=\frac{2}{1+\rho}=1.11,
```

```math
r_\mathrm{opposite}=\sigma(1,-1)^T,
\qquad
\chi^2_\mathrm{opposite}=\frac{2}{1-\rho}=10.
```

A common shift is plausible because the measurements share noise. Opposite
shifts are not. Treating the points as independent would assign ``\chi^2=2`` to
both patterns and therefore misstate both goodness of fit and parameter
uncertainty.

Typical sources of correlation are shared calibration constants, common
background subtraction, baseline drift, and finite-memory electronics. A
global systematic scale uncertainty is often clearer as a nuisance parameter
or external parameter constraint than as an arbitrary short-range covariance.

The representation should follow the mechanism:

| uncertainty source | useful representation |
| --- | --- |
| independent readout scatter | pointwise ``\sigma_i`` |
| repeated samples sharing noise | covariance matrix or whitening operator |
| uncertain calibration constant | fitted nuisance parameter with auxiliary information |
| plausible but unquantified bias | sensitivity analysis, not an invented Gaussian error |

A quantified systematic effect is therefore not automatically an extra number
to attach after the fit. If it moves several observations coherently, encode
that mechanism in the joint model before parameter uncertainties and goodness
of fit are computed. If its size is not probabilistically quantified, vary
defensible alternatives and report the sensitivity instead of hiding a guessed
distribution inside the covariance.

ScientificFitting never forms ``V^{-1}`` explicitly. Dense covariance is factorized and
applied through triangular solves, which is both more stable and faster than
materializing an inverse.

For a Cholesky factorization ``V=LL^T``, define

```math
z=L^{-1}r.
```

Then

```math
r^T V^{-1}r
=r^T L^{-T}L^{-1}r
=z^Tz.
```

This transformation is **whitening**: it rotates and rescales the correlated
residual vector so that, under a correct Gaussian model, ``z`` has covariance
``I``. It does not smooth, average, or discard observations. It expresses the
same covariance model in coordinates where ordinary squared pulls can be
summed.

## Structured Whitening

A whitening operation ``W`` satisfies

```math
W^TW=V^{-1},
\qquad
\chi^2=\lVert Wr\rVert^2.
```

This is the same statistical model in a data structure that can exploit
problem-specific sparsity or recurrences. For an AR(1) residual process with

```math
V_{ij}=\sigma^2\rho^{|i-j|},
```

the interior innovations are proportional to

```math
(Wr)_i=\frac{r_i-\rho r_{i-1}}{\sigma\sqrt{1-\rho^2}}.
```

The operation is ``O(n)`` although the equivalent dense matrix contains
``O(n^2)`` values and its generic factorization costs ``O(n^3)``.

`WhiteningOperator` accepts the whitening operation and ``\log\det V``. The
log determinant does not change a static chi-square minimum, but it is required
for a normalized Gaussian likelihood and comparable information criteria. The
operator represents one complete, static observation covariance. A
parameter-dependent covariance needs a parameter-dependent likelihood path;
one fixed whitening operation and one fixed log determinant cannot describe it.
Before using a custom operator at scale, verify ``\lVert Wr\rVert^2`` and
``\log\det V`` against a small dense reference problem.

## Full Gaussian Likelihood

If

```math
d\sim\mathcal N\!\left(m(\theta),V(\theta)\right),
```

then ScientificFitting evaluates

```math
C(\theta)
=
n\log(2\pi)
+\log\det V(\theta)
+r(\theta)^T V(\theta)^{-1}r(\theta).
```

If ``V`` is independent of ``\theta``, the first two terms are constant in the
optimization. The full Gaussian ``-2\log L`` cost and chi-square therefore have
the same minimizer, although the normalized likelihood value is needed for
AIC/BIC.

### Parameter-Dependent Covariance

If ``V`` changes with ``\theta``, the determinant term changes the optimum and
must not be dropped. A simple example is a scale model whose uncertainty is
specified as a fraction of its prediction: increasing the prediction increases
both the residual scale and ``\det V``. Keeping only the quadratic residual
term would reward arbitrarily inflated uncertainties.

ScientificFitting selects `cost=:gaussian_likelihood` for parameter-dependent covariance
when `cost=:auto`. This occurs for effective x-uncertainty propagation and
active model-relative uncertainty components.

## Uncertainty In X

For ``y=f(x,\theta)``, a small perturbation in x changes y by

```math
\delta y \approx \frac{\partial f}{\partial x}\,\delta x.
```

First-order propagation therefore gives

```math
V_\mathrm{eff}(\theta)
=
V_y+J_x(\theta)V_xJ_x(\theta)^T,
```

where ``J_x`` contains ``\partial f(x_i,\theta)/\partial x_i``. For independent
x and y errors,

```math
\sigma_{\mathrm{eff},i}^2(\theta)
=
\sigma_{y,i}^2
+
\left(\frac{\partial f}{\partial x}\right)^2
\sigma_{x,i}^2.
```

ScientificFitting's pointwise propagation uses the diagonal matrix
``J_x=\operatorname{diag}(\partial f(x_i,\theta)/\partial x_i)``. Correlations
within ``V_x`` are retained, but cross-covariance between measured x and y is
not represented by this approximation.

For a local slope of ``3``, ``\sigma_x=0.2``, and ``\sigma_y=0.4``, the x error
alone contributes ``0.6`` in y units and

```math
\sigma_\mathrm{eff}=\sqrt{0.4^2+0.6^2}=0.72.
```

Because the slope can depend on fitted parameters, so can
``V_\mathrm{eff}``. The full Gaussian likelihood cost is then required.

This method is a local linearization, not a general errors-in-variables model.
It is appropriate for small x errors and smooth, single-valued models. Large x
errors, strong curvature across an error bar, latent true x values, selection
effects, or correlated x-y measurement errors require a more explicit
measurement model.

See [Full Covariance](gallery/full_covariance.md) and
[XY Uncertainties](gallery/xy_uncertainties.md) for executable examples.

---

**Previous chapter:** [Statistical Foundations](statistical_foundations.md) ·
**Next chapter:** [Parameters and Fit Quality](parameter_inference.md), for
auxiliary measurements, degrees of freedom, goodness of fit, and local
parameter uncertainty.
