# JuFitter

<section class="jufitter-hero">
<div class="jufitter-kicker">Scientific fitting for Julia</div>

JuFitter is a Julia package for scientific fitting. Its goal is simple:
give scientists and engineers one clear path from data and uncertainties to a
statistically meaningful fit, a readable report, and a plot that is already good
enough to show.

</section>

```julia
using JuFitter

x = collect(range(0.0, 10.0; length=200))
model(x, p) = @. p[1] * x + p[2]
sigma_y = fill(0.2, length(x))
y = model(x, [2.0, 1.0]) .+ sigma_y .* sin.(1.8 .* x)

fit = fitplot(
    model,
    x,
    y;
    p0=[1.0, 0.0],
    sigma_y=sigma_y,
    parameter_names=["m", "b"],
    filename="fit.pdf",
)

result = fit.result
```

## What JuFitter Optimizes For

- One-line fit plots that look publication-ready without manual repairs.
- Explicit Gaussian, Poisson, histogram, unbinned, extended, indexed, custom,
  and multi-dataset workflows.
- Covariance-aware uncertainty propagation, profile likelihoods, contours, and
  diagnostics when local errors are not trustworthy.
- Makie-based output with clean defaults, publication export, and full
  customization when needed.
- A documentation style that makes statistical fitting accessible without
  hiding the mathematics.

## Where To Start

<div class="jufitter-card-grid">
<div class="jufitter-card">
<h3>Start</h3>
<p>Install the package, understand first compile time, and produce the first fit plot.</p>
<p><a href="quickstart.html">Run the quickstart</a></p>
</div>
<div class="jufitter-card">
<h3>Gallery</h3>
<p>Start from complete scientific workflows: data, model, fit, diagnostics, interpretation.</p>
<p><a href="gallery.html">Open the gallery</a></p>
</div>
<div class="jufitter-card">
<h3>Guides</h3>
<p>Choose uncertainty inputs, diagnose bad fits, and decide what to inspect next.</p>
<p><a href="fitting_for_practitioners.html">Read the practitioner guide</a></p>
</div>
<div class="jufitter-card">
<h3>Statistics</h3>
<p>Understand and justify the methods behind chi-square, likelihoods, covariance, profiles, and contours.</p>
<p><a href="statistical_foundations.html">Read the theory</a></p>
</div>
</div>
