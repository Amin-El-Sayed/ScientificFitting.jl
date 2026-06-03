# Linear Calibration

This is the smallest useful JuFitter workflow: measured calibration points,
point-by-point uncertainties, a weighted fit, and a plot that states exactly
what its uncertainty band means.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/linear_calibration_light.png" alt="Linear calibration fit">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/linear_calibration_dark.png" alt="Linear calibration fit in dark mode">
```

## Model

```math
U(x) = m x + b
```

The uncertainty is heteroscedastic: points at larger ``x`` are measured with a
slightly larger standard uncertainty. JuFitter therefore minimizes

```math
\chi^2(m,b)=\sum_i
\left(\frac{U_i-(m x_i+b)}{\sigma_{U,i}}\right)^2.
```

The plotted band is a **1σ prediction band**. It combines the local parameter
uncertainty of the fitted line with the expected scatter of one new
measurement. This is intentionally wider than a pure confidence band for the
mean curve.

## Complete Code

```julia
using CairoMakie
using JuFitter
using LaTeXStrings

x = collect(range(0.0, 10.0; length=28))
scatter_scale = @. 0.22 + 0.025 * x
sigma_y = @. 0.10 + 0.012 * x
y = @. 0.82 + 1.72 * x +
         scatter_scale * (0.55 * sin(1.25 * x) + 0.18 * cos(3.7 * x))

model(x, p) = @. p[1] * x + p[2]

result = fit_model(model, x, y; p0=[1.5, 0.5], sigma_y=sigma_y)

plot_fit(
    result;
    title=L"\mathrm{Sensor\ calibration}",
    model_label=L"U(x)=m x + b",
    xlabel=L"x",
    xunit=L"\mathrm{mm}",
    ylabel=L"U",
    yunit=L"\mathrm{V}",
    parameter_names=[L"m", L"b"],
    latex_labels=true,
    latex_stats=true,
    band=:prediction,
    nsigma=1,
    band_label=L"1\sigma\ \mathrm{prediction\ band}",
    show_legend=true,
    stats_position=:right,
    stats_mode=:full,
)
```

## Run It

```bash
julia --project=. examples/gallery/09_docs_gallery_suite.jl
```
