# Constraints And Profiles

This workflow combines bounded parameters, an inequality constraint, a Gaussian
prior, and likelihood-profile diagnostics.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/constraints_priors_light.png" alt="Constrained quadratic fit">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/constraints_priors_dark.png" alt="Constrained quadratic fit in dark mode">
```

## Model

```math
y = a x^2 + b x + c,\qquad a \ge 0
```

The prior on ``c`` encodes external information, while the profile and contour
show whether the local covariance approximation is enough.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/curvature_profile_light.png" alt="Curvature profile">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/curvature_profile_dark.png" alt="Curvature profile in dark mode">
```

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/curvature_slope_contour_light.png" alt="Curvature slope contour">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/curvature_slope_contour_dark.png" alt="Curvature slope contour in dark mode">
```

## What It Shows

- `bounds` for hard physical ranges.
- `constraints` for inequalities.
- `parameter_priors` for external Gaussian information.
- `profile` and `contour` for nonlocal uncertainty diagnostics.
