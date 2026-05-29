# JuFitter

JuFitter is a Julia package for scientific fitting with a focus on beautiful
default plots, explicit statistical semantics, robust numerics, and readable
workflows for scientists and engineers.

```julia
using JuFitter

x = collect(range(0.0, 10.0; length=200))
model(x, p) = @. p[1] * x + p[2]
sigma_y = fill(0.2, length(x))
y = model(x, [2.0, 1.0]) .+ sigma_y .* sin.(1.8 .* x)

result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y)
plot_fit(result; filename="fit.pdf")
```

## Development Priorities

- One-line fit plots that look publication-ready without manual repairs.
- Full control for users who need LaTeX, custom fonts, colors, units, panels,
  reports, and export settings.
- Statistically explicit Gaussian, Poisson, histogram, unbinned, extended, and
  multi-fit workflows.
- Robust covariance handling and profile-likelihood uncertainty estimates.
- Documentation that explains both the API and the statistics behind it.

See the roadmap and plotting design pages for the active development plan.
