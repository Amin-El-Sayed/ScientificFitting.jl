# Full Covariance

Many real measurements have correlations between points. This example uses a
dense y-covariance matrix instead of independent error bars.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/full_covariance_decay_light.png" alt="Full covariance exponential fit">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/full_covariance_decay_dark.png" alt="Full covariance exponential fit in dark mode">
```

## Model

```math
y(t) = A e^{-\lambda t} + C
```

Use `cov_y` when neighboring points share readout noise, baseline drift, a
normalization uncertainty, or any systematic source that makes residuals
correlated.

## Core Call

```julia
result = fit_model(model, x, y; p0=[1.5, -0.7, 0.0], cov_y=cov_y)
```
