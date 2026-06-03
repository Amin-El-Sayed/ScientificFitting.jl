# Poisson And Histogram Fits

Count data should usually be fitted with a likelihood for counts, not by
pretending that every bin is Gaussian.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/poisson_counts_light.png" alt="Poisson count fit">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/poisson_counts_dark.png" alt="Poisson count fit in dark mode">
```

## Poisson Model

```math
\mu(x) = \exp(a + b x)
```

The model returns expected counts. The fit minimizes a Poisson likelihood and
the plot shows the natural ``\sqrt{\mu}`` count scale.

```@raw html
<img class="jufitter-plot jufitter-plot-light" src="../assets/gallery/histogram_likelihood_light.png" alt="Histogram likelihood fit">
<img class="jufitter-plot jufitter-plot-dark" src="../assets/gallery/histogram_likelihood_dark.png" alt="Histogram likelihood fit in dark mode">
```

## Histogram Model

For histogram fits, the model returns one expected count per bin. This is the
right interface when the binning is part of the measurement.

```julia
result = fit_histogram_model(expected_counts, edges, counts; p0, bounds)
```
