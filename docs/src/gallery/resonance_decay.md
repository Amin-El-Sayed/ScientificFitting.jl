# Damped Oscillator

This example fits a freely decaying mechanical oscillator. The dataset is dense
enough to stress marker styling and uncertainty bands, but still compact enough
to run quickly as a documentation example.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/damped_oscillator_decay_light.png" alt="Damped oscillator fit in the light documentation theme">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/damped_oscillator_decay_dark.png" alt="Damped oscillator fit in the dark documentation theme">
```

## Question

The measured angle decreases slowly while the system keeps oscillating. We want
the damping constant and angular frequency, including propagated parameter
uncertainties and a visible confidence band.

## Model

```math
\phi(t) = A e^{-\lambda t}\cos(\omega t + \phi_0).
```

The parameters are the initial amplitude ``A``, damped angular frequency
``\omega``, phase ``\phi_0``, and damping constant ``\lambda``. The fit uses
both angle uncertainty and a small timing uncertainty, so the cost function is
not just an unweighted least-squares curve through dense data.

## Workflow

1. Load the curated CSV with columns `time_s`, `phi_rad`, and `sigma_phi_rad`.
2. Define the four-parameter underdamped oscillator model.
3. Provide several initial guesses because phase-periodic nonlinear fits can
   have local minima.
4. Fit with bounds, angle uncertainties, and timing uncertainties.
5. Inspect ``\chi^2/\mathrm{ndf}``, ``P(\chi^2)``, parameter errors, and the
   residual structure before trusting the result.

The plot uses a **1σ prediction band**. That band includes propagated parameter
uncertainty plus the observation uncertainty, so it is the right band when the
reader wants to see the expected scale of new measurements. For a pure
parameter-confidence band, use `band=:confidence`; for a wider visual envelope,
increase `nsigma` explicitly and label it in the legend.

## Run It

From the repository root:

```bash
julia --project=. examples/gallery/08_damped_oscillator_decay.jl
```

The curated data is distributed with the examples:

```text
examples/data/damped_oscillator/pohl_wheel_free_decay.csv
```
