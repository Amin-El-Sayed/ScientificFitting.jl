# Diagnostic Plotting

These functions visualize an existing fit result or a stored profile scan.
Only `plot_profile_matrix(result)` performs new fits; passing a
`ProfileMatrixResult` renders stored numerical results.

```julia
using ScientificFitting
using CairoMakie
```

## Choose A Figure

| Question | Function | Performs optimization? |
|---|---|---:|
| Where does the model miss the data? | [`plot_residuals`](@ref) | no |
| Are residual, pull, and ratio views mutually consistent? | [`plot_diagnostics`](@ref) | no |
| Is one parameter locally parabolic? | [`plot_profile`](@ref) | no |
| Is a parameter pair described by a covariance ellipse? | [`plot_contour`](@ref) | no |
| Which parameters need closer inspection? | [`plot_profile_matrix`](@ref) | yes for a fit result; no for a stored matrix |

All functions accept `theme=:sans` or `:tex`, `appearance=:light` or
`:dark`, an optional `theme_override`, and standard file output keywords.

## Residuals, Pulls, And Ratios

```julia
plot_residuals(result; kind=:pull)
plot_diagnostics(result)
```

| `kind` | Displayed value | Reference line |
|---|---|---:|
| `:residual` | ``y_i-f_i`` with available y errors | 0 |
| `:pull` | weighted or whitened residual coordinate | 0 |
| `:ratio` | ``y_i/f_i`` with propagated y-error ratio | 1 |

With dense covariance, whitened coordinates are not pointwise pulls in the
original measurement order. Ratios are rejected when a fitted value is zero or
non-finite.

Shared keywords are `filename`, `format`, `theme`, `appearance`,
`theme_override`, `figure_size`, `xlabel`, `color`, `reference_color`,
`marker`, `markersize`, `error_whiskerwidth`, `axis_kwargs`,
`scatter_kwargs`, and `errorbars_kwargs`. `plot_diagnostics` also accepts
`reference_line_kwargs`.

## One-Parameter Profiles

```julia
plot_profile(profile_result; local_sigma=result.param_stderr[i])
```

| Concern | Keywords |
|---|---|
| Output and style | `filename`, `format`, `theme`, `appearance`, `theme_override`, `figure_size` |
| Labels | `title`, `xlabel`, `ylabel` |
| Profile | `line_color`, `line_width`, `profile_label`, `line_kwargs` |
| Local parabola | `local_sigma`, `local_color`, `local_linewidth`, `local_linestyle`, `local_label`, `local_line_kwargs` |
| Threshold | `threshold_color`, `threshold_label`, `threshold_kwargs` |
| Layout | `show_legend`, `legend_position`, `delta_max`, `axis_kwargs`, `legend_kwargs` |

`local_sigma` and `delta_max` must be positive. `delta_max` changes only the
displayed range. The default legend occupies a row below the data axis;
`legend_position=:right` selects a bounded side column.

## Two-Parameter Contours

```julia
plot_contour(
    contour_result;
    local_covariance=result.param_covariance,
    local_center=result.params[[i, j]],
)
```

| Concern | Keywords |
|---|---|
| Output and style | `filename`, `format`, `theme`, `appearance`, `theme_override`, `figure_size` |
| Labels | `title`, `xlabel`, `ylabel`, `axis_kwargs` |
| Profile regions | `show_regions`, `show_profile_lines`, `level_colors`, `region_colors`, `line_color`, `contour_kwargs` |
| Optional heatmap | `show_heatmap`, `colormap`, `heatmap_kwargs` |
| Local covariance | `local_covariance`, `local_center`, `local_line_color`, `local_linewidth`, `local_linestyle`, `local_contour_kwargs` |
| Legend | `show_legend`, `legend_position`, `legend_kwargs` |

Filled regions are the default; the common two-parameter thresholds 2.30 and
6.18 are labeled as one- and two-sigma regions. The covariance approximation
is drawn separately. A heatmap is opt-in.

## Profile Matrices

```julia
matrix = profile_matrix(result; parameters=[1, 2, 3])
plot_profile_matrix(matrix; parameter_names=["A", "lambda", "offset"])
```

`plot_profile_matrix(result)` computes the required profile and contour refits.
Its scan controls are `parameters`, `parameter_names`, `npoints_profile`,
`npoints_contour`, `nsigma`, `profile_threshold`, `contour_levels`, `adaptive`,
`max_refinements`, and `max_points`.

Both methods accept `filename`, `format`, `theme`, `appearance`,
`theme_override`, `panel_status_mode`, `delta_max`, and `figure_size`.
`panel_status_mode` is `:issues`, `:all`, or `:none` and controls labels
independently of visual style.

## Failure Contract

| Invalid request | Result |
|---|---|
| CairoMakie extension not loaded | `ArgumentError` naming CairoMakie |
| Unsupported residual kind or non-finite coordinates | `ArgumentError` |
| Ratio with a zero/non-finite model prediction | `ArgumentError` |
| Non-positive profile display scale | `ArgumentError` |
| Incompatible local covariance or contour geometry | `ArgumentError` |
| Invalid profile-matrix status mode or display names | `ArgumentError` |

## API Documentation

```@docs
ScientificFitting.plot_residuals
ScientificFitting.plot_diagnostics
ScientificFitting.plot_profile
ScientificFitting.plot_contour
ScientificFitting.plot_profile_matrix
```
