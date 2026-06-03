# Photoelectric Work Function

This example uses the photoelectric effect as an actual analysis workflow, not
only as a line-fit demo. The measured stopping voltage is fitted as a function
of light frequency, then the fit is extrapolated to extract the work function
from the y-axis intercept.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/photoelectric_threshold_light.png" alt="Photoelectric work-function fit">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/photoelectric_threshold_dark.png" alt="Photoelectric work-function fit in dark mode">
```

## Physics

```math
eU_0 = h\nu - \Phi
```

Equivalently, with frequency in THz,

```math
U_0(\nu_\mathrm{THz}) = m\nu_\mathrm{THz} + b,
```

where ``m = 10^{12}h/e`` and ``b=-\Phi/e``. Numerically, the work function in
electronvolts is simply

```math
\Phi[\mathrm{eV}] = -b[\mathrm{V}].
```

The threshold frequency is the x-intercept:

```math
\nu_0 = -\frac{b}{m}.
```

## Error Propagation

The fit returns the covariance matrix of ``m`` and ``b``:

```math
V =
\begin{pmatrix}
\sigma_m^2 & \mathrm{cov}(m,b) \\
\mathrm{cov}(m,b) & \sigma_b^2
\end{pmatrix}.
```

The work-function uncertainty is direct:

```math
\sigma_\Phi[\mathrm{eV}] = \sigma_b[\mathrm{V}].
```

For the threshold frequency, use first-order propagation:

```math
\sigma_{\nu_0}^2 =
\begin{pmatrix} b/m^2 & -1/m \end{pmatrix}
V
\begin{pmatrix} b/m^2 \\ -1/m \end{pmatrix}.
```

## Workflow

1. Convert wavelengths to frequencies.
2. Fit only the points above threshold with a linear model.
3. Keep the below-threshold points in the plot as visual evidence for the
   physical cutoff.
4. Extrapolate the fitted line to ``\nu=0`` and ``U_0=0``.
5. Mark ``\Phi`` and ``\nu_0`` directly in the plot with propagated
   uncertainties.

## Complete Code

```julia
using CairoMakie
using JuFitter
using LinearAlgebra

const c = 299_792_458.0
const e = 1.602176634e-19

wavelength_nm = [150.0, 200.0, 250.0, 300.0, 350.0,
                 400.0, 450.0, 500.0, 550.0, 600.0]
voltage = [5.99, 3.87, 2.69, 1.78, 1.28, 0.77, 0.50, 0.15, 0.0, 0.0]
sigma_wavelength_nm = fill(0.01, length(wavelength_nm))
sigma_voltage = fill(0.04, length(voltage))

frequency_THz = @. c / (wavelength_nm * 1e-9) / 1e12
sigma_frequency_THz = @. c * (sigma_wavelength_nm * 1e-9) /
                         (wavelength_nm * 1e-9)^2 / 1e12

fit_mask = voltage .> 0.0
model(nu_THz, p) = @. p[1] * nu_THz + p[2]

result = fit_model(
    model,
    frequency_THz[fit_mask],
    voltage[fit_mask];
    p0=[0.004, -2.2],
    sigma_y=sigma_voltage[fit_mask],
    sigma_x=sigma_frequency_THz[fit_mask],
    bounds=([0.0, -20.0], [0.02, 5.0]),
    initial_guesses=[[0.004, -2.2], [0.0042, -2.6], [0.0038, -1.8]],
)

slope, intercept = result.params
cov = result.param_covariance

h_fit = slope * e / 1e12
sigma_h = sqrt(cov[1, 1]) * e / 1e12

work_function_eV = -intercept
sigma_work_function_eV = sqrt(cov[2, 2])

threshold_THz = -intercept / slope
threshold_gradient = [intercept / slope^2, -1 / slope]
sigma_threshold_THz = sqrt(dot(threshold_gradient, cov * threshold_gradient))

println("h = ", h_fit, " +/- ", sigma_h, " J s")
println("Phi = ", work_function_eV, " +/- ", sigma_work_function_eV, " eV")
println("nu0 = ", threshold_THz, " +/- ", sigma_threshold_THz, " THz")
```

The documentation asset generator adds the Makie annotations for the intercept,
threshold marker, uncertainty band, and report panel.
