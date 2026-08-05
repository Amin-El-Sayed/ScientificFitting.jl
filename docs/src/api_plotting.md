# Plotting

Fitting, profiles, diagnostics, and text reports do not require Makie. Plotting
is activated only after loading CairoMakie:

```julia
using JuFitter
using CairoMakie
```

[`plot_fit`](@ref) returns a `Figure`. [`fitplot`](@ref) fits and returns the
named tuple `(result, figure)`. Use `fit_axis(figure)` as the stable extension
point for native Makie calls or JuFitter's `add_*!` helpers.

For `fitplot`, `report=:plot`, `:console`, `:both`, or `:none` controls the
information panel and terminal output. For `plot_fit`, use `show_stats` to
control the panel. `fit_range=:axis` extrapolates the fitted curve over the
padded visible x range; use `fit_range=:data` or an explicit `xgrid` at a
physical domain boundary.

The maintained output roles are `:lab`, `:screen`, and `:article`. Explicit
Makie keywords override defaults only for the element receiving them. See
[Plotting And Customization](plotting_design.md) for layout and extension
examples. Visual role design is still subject to maintainer acceptance before
v0; the API contracts on this page do not imply final approval of the presets.

```@docs
JuFitter.fitplot
JuFitter.plot_fit
JuFitter.fit_axis
JuFitter.add_curve!
JuFitter.add_points!
JuFitter.add_vline!
JuFitter.add_hline!
JuFitter.add_vband!
JuFitter.add_hband!
JuFitter.plot_theme
JuFitter.plot_palette
JuFitter.plot_info_panel!
JuFitter.plot_residuals
JuFitter.plot_diagnostics
JuFitter.plot_profile
JuFitter.plot_contour
JuFitter.plot_profile_matrix
```
