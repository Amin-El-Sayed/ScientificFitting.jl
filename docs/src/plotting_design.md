# Plotting Design

Plotting is a primary feature, not decoration after fitting.

## Default User Story

The common workflow should be short and reliable:

```julia
fitplot(x, y; sigma_y)
fitplot(model, x, y; p0, sigma_y)
plot_fit(result)
```

The output should have sensible margins, readable labels, visible uncertainties,
a clear fit curve, and an optional compact summary without requiring manual
layout tuning.

## Design Principles

- Beautiful defaults first; customization second.
- Automatic layout must account for error bars, bands, labels, legends, and
  statistics panels.
- All high-level options should map cleanly to Makie concepts.
- Export quality must be consistent across PNG, PDF, and SVG.
- Plot functions should return the Makie `Figure` for further user control.

## Planned High-Level Controls

- `theme=:clean | :paper | :latex | custom`
- `xlabel`, `ylabel`, `xunit`, `yunit`, `title`
- `parameter_names`
- `report=:plot | :console | :both | :none`
- `band=:none | :confidence | :prediction`
- `confidence_level` or `nsigma`
- `show_grid`, `show_legend`, `show_residuals`, `show_pulls`
- `axis_kwargs`, `line_kwargs`, `scatter_kwargs`, `band_kwargs`,
  `legend_kwargs`

## Acceptance Tests

- Small datasets do not produce cramped plots.
- Dense datasets do not obscure the fit curve.
- Error bars and confidence bands are fully inside the visible range.
- Long labels, LaTeX labels, and parameter summaries are not clipped.
- Multi-dataset and multi-fit plots use clear visual hierarchy.
