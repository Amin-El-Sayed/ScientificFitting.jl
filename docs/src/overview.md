# JuFitter Overview

JuFitter provides:

- `LsqFit` default backend for unconstrained/static-weight fits
- `Optimization.jl` backend for constrained fits and dynamic effective-variance weighting
- covariance-aware uncertainty propagation to parameter errors and confidence bands
- CairoMakie publication-style output
- extractable `FitReport` objects via `fit_report`
- readable plain-text summaries via `report_text`

See `examples/basic_fit.jl` for a complete run.

For the planned statistical architecture and implementation rules, see
`docs/statistical_foundations.md`.

For the current cost-function backend, see `docs/backend_design.md`.
