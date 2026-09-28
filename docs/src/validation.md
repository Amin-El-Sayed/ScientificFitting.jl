# Validation

This page collects the external references the test suite checks on every
run and the measured performance behind the package's scaling claims. All
timings below were taken on an Apple M3 Pro with Julia 1.13.0; rerun the
quoted scripts for your hardware.

## NIST StRD Certified Problems

The [NIST Statistical Reference Datasets](https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml)
publish nonlinear regression problems with certified parameter estimates,
standard deviations, and residual sums of squares to 11 significant digits.
`test/statistics/nist_strd_reference.jl` fits eight of them, unweighted, from
Start 2 — the second, closer of NIST's two published sets of starting values —
on every core-suite run and checks all three certified quantities:

| Problem | n | Parameters | Difficulty (NIST) |
| --- | ---: | ---: | --- |
| Misra1a | 14 | 2 | lower |
| Chwirut2 | 54 | 3 | lower |
| Gauss1 | 250 | 8 | lower |
| MGH17 | 33 | 5 | average |
| Rat42 | 9 | 3 | higher |
| Thurber | 37 | 7 | higher |
| BoxBOD | 6 | 2 | higher |
| Eckerle4 | 35 | 3 | higher |

Agreement: parameters within a relative tolerance of ``10^{-6}``
(``10^{-4}`` for MGH17, ``10^{-5}`` for Thurber), certified residual sums of
squares within ``10^{-9}`` (``10^{-8}`` for MGH17), and certified standard
deviations within ``10^{-3}``–``10^{-4}``. The standard-deviation check is an external
validation of the ``\chi^2/\mathrm{ndf}`` covariance scaling for fits without
supplied uncertainties
([Covariance Scaling](statistics.md#Covariance-Scaling)).

## Independent Cross-Checks

- The [LHCb mass spectrum](gallery/lhcb_mass_spectrum.md) example ships an
  independent NumPy/SciPy/iminuit implementation (`lhcb_reference.py`,
  including a MINOS interval) whose results the documented fit reproduces.
- Results for fits with x uncertainties differ reproducibly but slightly from
  ODR-convention tools (`scipy.odr`, kafe2, York); the difference and its
  Monte-Carlo quantification are derived in the
  [Statistics Reference](statistics.md#Relation-To-ODR-And-The-York-Method).
- The residual-structure diagnostics are calibrated by construction:
  `test/statistics/scaling_consistency_reference.jl` runs 300 correct fits at
  each of ``n \in \{15, 25, 50, 200, 500\}`` and requires a per-code
  false-positive rate at or below 8%.

## Matrix-Free Whitening At Large n

A structured covariance does not require its dense matrix. For a fit whose
residual covariance is modeled as an AR(1) process (first-order
autoregressive: each point correlated with its neighbor) through a
`WhiteningOperator` (`benchmarks/whitening_scaling.jl`), after warm-up:

| n | operator fit | dense-covariance fit | dense memory |
| ---: | ---: | ---: | ---: |
| 10⁴ | 0.005 s | 17.4 s | 0.8 GB |
| 10⁵ | 0.021 s | — | 80 GB (infeasible) |
| 10⁶ | 0.26 s | — | 8 TB (infeasible) |

Operator and dense path agree in the fitted parameters where both run; the
dense column stops where the covariance no longer fits in memory. The
dense-memory column is not measured but computed: storage for the covariance
matrix alone is ``8n^2`` bytes of `Float64`. The operator path is ``O(n)`` in
time and memory for this structure.

## Time To First Fit

In Julia, every named function or closure is its own type, so each distinct
model function normally triggers a fresh compilation of the fitting
pipeline. With `derivatives=:finite`, all models pass through one fixed
callback signature and the pipeline compiles once
([Backend Design](backend_design.md)). `benchmarks/ttfx_probe.jl` times
every cell below in its own fresh Julia process, so each number is what a
new session pays; the script prints the table rows directly:

| Measurement | `derivatives=:auto` | `derivatives=:finite` |
| --- | ---: | ---: |
| first fit after `using` | 5.0 s | 0.5 s |
| each additional model type | 5.1 s | 0.06 s |
| Poisson fit, additional rate model | 9.5 s | 0.11 s |

Use `derivatives=:finite` when exploring model variants interactively; use
the default `:auto` for production fits, where automatic differentiation
specializes once per model and then runs at full speed. The finite-difference
path is checked against closed-form least-squares solutions in
`test/numerics/finite_derivatives_reference.jl`.
