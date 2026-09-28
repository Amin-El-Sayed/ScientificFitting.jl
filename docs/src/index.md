# ScientificFitting.jl

Nonlinear curve and likelihood fitting with explicit uncertainties: y errors,
x errors, full covariance matrices, matrix-free whitening, Poisson and
unbinned likelihoods, parameter priors and constraints, profile-likelihood
intervals, diagnostics, and Makie figures — one workflow from data to
reported result, in Julia and [from Python](python.md).

```@example landing
using ScientificFitting

x = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
U = [1.31, 2.13, 3.09, 3.90, 5.15, 5.95, 7.11, 7.86]
sigma_U = fill(0.12, length(U))

result = fit_model((x, p) -> p[1] .* x .+ p[2], x, U;
                   p0=[1.0, 0.0], sigma_y=sigma_U)
println(report_text(result; parameter_names=["m", "b"]))
```

Where to go next:

- **First complete fit**, including the plot and the diagnosis:
  [Quickstart](quickstart.md).
- **Which function fits your data**: the
  [entry-point table](api.md#Choose-An-Entry-Point).
- **Complete executable analyses**, from x-y uncertainties to an LHCb mass
  spectrum: [Gallery](gallery.md).
- **What the numbers mean**, derived once and precisely:
  [Statistics Reference](statistics.md).
- **What this package does not do**: [Scope and Alternatives](limits.md).
