# Multi-Dataset Fit: Which Parameters May Be Shared?

Several datasets often measure the same physics without sharing every
instrument parameter. A simultaneous fit can use that structure, but sharing a
parameter is a scientific hypothesis, not a numerical convenience.

This controlled calibration-transfer example asks whether three readout
channels may use one common gain. The answer is no: channels A and B are
compatible, while channel C requires a separate gain.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/multi_dataset_shared_slope_light.png" alt="Three-channel calibration transfer with all-shared and partial-sharing models, local fit bands, and pull panels">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/multi_dataset_shared_slope_dark.png" alt="Three-channel calibration transfer with all-shared and partial-sharing models, local fit bands, and pull panels in dark mode">
```

The solid lines show the accepted partial-sharing model. The dashed lines show
what happens when all channels are forced to share one gain. In the main panel
the difference is easy to underestimate; the first pull panel makes the failure
of the all-shared hypothesis clear.

## The Scientific Question

Three nominally identical readout channels are calibrated against the same
reference input. Their electronic zero points may differ, but the proposed
transfer procedure assumes one common sensitivity:

```math
y_i(x)=g x+b_i,
\qquad i\in\{A,B,C\}.
```

If this assumption is valid, one accurately calibrated gain can be transferred
between channels after correcting their offsets. If one channel has a different
gain, that transfer introduces a systematic error that grows with input.

This page uses a controlled, reproducible benchmark dataset designed to expose
that decision. The three channels have different x sampling, heteroscedastic
absolute y uncertainties, independent offsets, imperfect scatter, and a
deliberately incompatible gain for channel C.

The controlled construction is useful here because the correct sharing
structure is known. It tests whether the analysis finds the hidden gain
mismatch rather than merely producing three attractive lines.

## The Multi-Dataset Cost

For independent Gaussian measurements, the simultaneous cost is

```math
\chi^2(p)
= \sum_i \sum_j
  \left[
    \frac{y_{ij}-f_i(x_{ij};p_i)}{\sigma_{ij}}
  \right]^2.
```

Each local model receives only the global parameters listed by its
`parameter_map`. The first hypothesis uses

```math
p_\mathrm{all}=(g,b_A,b_B,b_C),
```

with maps

```math
p_A=(p_1,p_2),\qquad
p_B=(p_1,p_3),\qquad
p_C=(p_1,p_4).
```

All three channels therefore use exactly the same gain ``g``.

The second hypothesis uses

```math
p_\mathrm{partial}=(g_{AB},b_A,b_B,g_C,b_C),
```

with maps

```math
p_A=(p_1,p_2),\qquad
p_B=(p_1,p_3),\qquad
p_C=(p_4,p_5).
```

Channels A and B share ``g_{AB}``; channel C has its own ``g_C``.

## Complete Fit Code

```julia
using JuFitter
using LinearAlgebra

linear_channel(x, p) = @. p[1] * x + p[2]

x_a = collect(0.0:1.0:10.0)
x_b = collect(0.5:1.0:9.5)
x_c = collect(0.0:1.25:10.0)

sigma_a = @. 0.075 + 0.008 * x_a
sigma_b = @. 0.085 + 0.006 * x_b
sigma_c = @. 0.080 + 0.009 * x_c

pattern_a = 1.65 .* [0.2, -0.7, 0.4, -0.5, 0.8, -0.3, 0.1, -0.8, 0.6, -0.2, 0.5]
pattern_b = 1.65 .* [-0.4, 0.7, -0.5, 0.1, 0.8, -0.6, 0.3, -0.2, 0.6, -0.3]
pattern_c = 1.65 .* [0.3, -0.6, 0.7, -0.4, 0.1, 0.8, -0.5, 0.4, -0.7]

y_a = linear_channel(x_a, [1.82, 0.72]) .+ sigma_a .* pattern_a
y_b = linear_channel(x_b, [1.82, -0.46]) .+ sigma_b .* pattern_b
y_c = linear_channel(x_c, [1.91, 0.12]) .+ sigma_c .* pattern_c

models = [linear_channel, linear_channel, linear_channel]
x_sets = [x_a, x_b, x_c]
y_sets = [y_a, y_b, y_c]
sigma_sets = [sigma_a, sigma_b, sigma_c]

all_shared_result = fit_multi_model(
    models,
    x_sets,
    y_sets;
    p0=[1.85, 0.7, -0.4, 0.1],
    sigma_y=sigma_sets,
    parameter_map=[[1, 2], [1, 3], [1, 4]],
    parameter_names=["shared gain", "offset A", "offset B", "offset C"],
)

partial_shared_result = fit_multi_model(
    models,
    x_sets,
    y_sets;
    p0=[1.82, 0.7, -0.4, 1.90, 0.1],
    sigma_y=sigma_sets,
    parameter_map=[[1, 2], [1, 3], [4, 5]],
    parameter_names=["gain A/B", "offset A", "offset B", "gain C", "offset C"],
)

gain_gap_gradient = [-1.0, 0.0, 0.0, 1.0, 0.0]
gain_gap = partial_shared_result.params[4] - partial_shared_result.params[1]
sigma_gain_gap = sqrt(dot(
    gain_gap_gradient,
    partial_shared_result.param_covariance * gain_gap_gradient,
))

println(report_text(all_shared_result))
println(report_text(partial_shared_result))
println("gain C - gain A/B = ", gain_gap, " +/- ", sigma_gain_gap)
```

`parameter_map` is the essential part of the interface. The same local function
`linear_channel(x, p)` is reused for every channel, while the maps define which
global parameters each call receives.

## Diagnose The All-Shared Hypothesis

The fully shared model converges. Its fitted gain is a compromise between
incompatible channels:

```math
g_\mathrm{all}=1.8478\pm0.0068.
```

Convergence does not make that compromise physically valid. The fit gives

```math
\frac{\chi^2_\mathrm{all}}{\mathrm{ndf}}=2.035,
\qquad
P(\chi^2_\mathrm{all})=0.00139.
```

The first pull panel shows the specific failure. Channel C is systematically
below the common-gain model at low input and above it at high input. That sign
change is the residual signature of a slope mismatch. Channels A and B are
pulled in the opposite direction because the shared gain is forced to
compromise.

This is why inspecting only the global ``\chi^2`` is insufficient. The global
number says that something is wrong; per-dataset pulls show which sharing
assumption is wrong.

## Fit Only The Defensible Sharing Structure

Allowing channel C to have an independent gain gives

```math
g_{AB}=1.8232\pm0.0081,
\qquad
g_C=1.9060\pm0.0124.
```

The gain difference must be propagated from the **joint** covariance matrix:

```math
\Delta g = g_C-g_{AB},
\qquad
\sigma_{\Delta g}^2
= \nabla\Delta g^\mathsf{T}
  \operatorname{Cov}(p)
  \nabla\Delta g.
```

Here,

```math
\Delta g=0.0828\pm0.0148,
```

which differs from zero by approximately ``5.6\sigma`` under the local Gaussian
approximation.

The revised fit gives

```math
\frac{\chi^2_\mathrm{partial}}{\mathrm{ndf}}=0.870,
\qquad
P(\chi^2_\mathrm{partial})=0.651.
```

The second pull panel no longer contains a channel-dependent trend. Sharing the
gain between A and B is supported by this dataset; transferring channel C's
gain is not.

## Why Compare AIC Here?

The partial-sharing model adds one parameter. A lower ``\chi^2`` is therefore
expected even if the extra freedom is unnecessary. AIC adds a parameter-count
penalty:

```math
\mathrm{AIC}=2k-2\log L_\max.
```

Both models use the same data and Gaussian cost, so their AIC values may be
compared:

```math
\Delta\mathrm{AIC}
= \mathrm{AIC}_\mathrm{all}-\mathrm{AIC}_\mathrm{partial}
\approx 29.2.
```

That improvement is much larger than the penalty for one extra gain parameter.
It reinforces the residual and gain-difference evidence against the
all-shared-gain hypothesis.

AIC is not proof that the partial-sharing model is physically true. It only
compares the candidate models supplied here. A nonlinear response, correlated
calibration errors, or a shared reference-standard uncertainty could still
require another model.

## What The Bands Mean

The main panel shows one **local 1σ fit band** for each channel under the
partial-sharing model. For a local linear model,

```math
\sigma_\mathrm{fit}^2(x)
= J(x)\,\operatorname{Cov}(p_i)\,J^\mathsf{T}(x),
\qquad
J(x)=(x,1).
```

These bands describe uncertainty in each fitted mean response. They are not
prediction bands for a future observation and therefore do not add
``\sigma_{ij}^2``. They also do not include uncertainty in the reference x
values or correlations caused by a shared calibration standard.

## What Can Go Wrong

- **Sharing by naming convention:** parameters should be shared because the
  physical model requires it, not because two columns have the same label.
- **Reading only the aggregate fit statistic:** one failing dataset can be
  hidden inside an acceptable-looking global curve.
- **Fitting every dataset independently:** this discards real shared
  information and prevents direct propagation through the joint covariance.
- **Ignoring shared systematic uncertainty:** `fit_multi_model` currently
  supports dataset-specific diagonal `sigma_y`; a common reference-standard
  error requires a richer covariance treatment than this example.
- **Using AIC across different data or likelihoods:** the comparison is
  meaningful here because both candidates use exactly the same observations
  and cost definition.

## Reproduce The Figure

The tracked script contains the complete fits, covariance propagation, pull
construction, local bands, and light/dark Makie figure:

```bash
julia --project=. examples/gallery/10_multi_dataset_calibration.jl
```

The practical conclusion is precise: **transfer the gain between channels A and
B, but calibrate channel C separately.**

Next, use [Full Covariance](full_covariance.md) when datasets share systematic
uncertainty, or [Constraints And Profiles](constraints_profiles.md) when the
shared-parameter geometry is nonlinear and local errors may fail.
