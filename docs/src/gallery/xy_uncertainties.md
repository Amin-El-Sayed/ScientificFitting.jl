# XY Uncertainties

This example shows effective-variance fitting: uncertainty in x changes the
vertical cost through the local model slope.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/xy_uncertainties_light.png" alt="XY uncertainty fit">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/xy_uncertainties_dark.png" alt="XY uncertainty fit in dark mode">
```

## Model

```math
y = m x + b
```

For a smooth model, an x uncertainty ``\sigma_x`` contributes approximately
``(f'(x)\sigma_x)^2`` to the effective y variance. That is statistically
different from only drawing horizontal error bars.

## Use This When

- calibration points have uncertainty in both axes,
- the model slope is steep enough that x uncertainty matters,
- ignoring x errors would visibly underestimate parameter errors.
