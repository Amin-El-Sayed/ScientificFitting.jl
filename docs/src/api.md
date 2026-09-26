# API Reference

Start with [Quickstart](quickstart.md) for a first fit or
[Statistical Foundations](statistics.md) for the derivations behind each
method.

## Choose An Entry Point

| Data and sampling model | Entry point | Model contract | Result |
|---|---|---|---|
| Numeric ``x`` and ``y`` with Gaussian uncertainties | [`fit_model`](@ref) | `model(x, p) -> y_hat` | [`FitResult`](@ref) |
| Independent observations with a custom distribution | [`fit_likelihood_model`](@ref) | `model(x, p)` and `logprob(y, y_hat, p)` | [`LikelihoodFitResult`](@ref) |
| Events described by an upstream distribution | [`fit_distribution`](@ref) | `make_distribution(p) -> Distribution`, or an optional BuildConstructors model | [`LikelihoodFitResult`](@ref) |
| Binned events described by an upstream distribution | [`fit_distribution`](@ref) | Same factory, `edges, counts`; known `total_count` or extended mixture yields | [`LikelihoodFitResult`](@ref) |
| Independent counts | [`fit_poisson_model`](@ref) | `model(x, p) -> expected_counts` | [`LikelihoodFitResult`](@ref) |
| Histogram with expected bin counts | [`fit_histogram_model`](@ref) | `expected_counts(edges, p) -> mu` | [`LikelihoodFitResult`](@ref) |
| Histogram from a normalized density | [`fit_histogram_density`](@ref) | `pdf(x, p) -> density` | [`LikelihoodFitResult`](@ref) |
| Independent unbinned observations | [`fit_unbinned_model`](@ref) | `pdf(x, p) -> density` | [`LikelihoodFitResult`](@ref) |
| Unbinned events with a parameter-dependent rate | [`fit_extended_unbinned_model`](@ref) | `rate(x, p) -> event_rate` | [`LikelihoodFitResult`](@ref) |
| Observations addressed by non-numeric indices | [`fit_indexed_model`](@ref) | `model(indices, p) -> y_hat` | [`LikelihoodFitResult`](@ref) |
| Several datasets sharing parameters | [`fit_multi_model`](@ref) | one `model_i(x_i, p_i)` per dataset | [`LikelihoodFitResult`](@ref) |
| A custom scalar objective | [`fit_custom`](@ref) | `objective(p) -> scalar` | [`LikelihoodFitResult`](@ref) |

Minimal call shapes, with required keywords shown:

| Entry point | Minimal call |
|---|---|
| `fit_model` | `fit_model(model, x, y; p0=[...], sigma_y=[...])` |
| `fit_likelihood_model` | `fit_likelihood_model(model, x, y; p0=[...], logprob=logprob)` |
| `fit_distribution` | `fit_distribution(make_distribution, observations; p0=[...])` |
| `fit_poisson_model` | `fit_poisson_model(expected_counts, x, counts; p0=[...])` |
| `fit_histogram_model` | `fit_histogram_model(expected_per_bin, edges, counts; p0=[...])` |
| `fit_histogram_density` | `fit_histogram_density(pdf, edges, counts; p0=[...], total_count=sum(counts))` |
| `fit_unbinned_model` | `fit_unbinned_model(pdf, observations; p0=[...])` |
| `fit_extended_unbinned_model` | `fit_extended_unbinned_model(rate, observations, (a, b); p0=[...])` |
| `fit_indexed_model` | `fit_indexed_model(model, indices, y; p0=[...], cov_y=C)` |
| `fit_multi_model` | `fit_multi_model(models, xs, ys; p0=[...], sigma_y=sigma_sets)` |
| `fit_custom` | `fit_custom(cost; p0=[...], nobs=n)` |

Complete analyses are in the [Gallery](gallery.md).

For reusable low-level workflows, construct [`FitProblem`](@ref) or
[`LikelihoodFitProblem`](@ref) and call [`fit`](@ref), which extends the
`StatsAPI.fit` generic.

## Common Conventions

### Parameters And Model Functions

`p0` fixes the parameter order. Indices in bounds, fixed values, constraints,
profiles, and contours are one-based into that vector.

The allocating model contract is:

```text
model(x, p) -> vector with length(y)
```

For allocation-sensitive fits, use:

```text
model!(out, x, p)
jacobian!(J, x, p)  # optional
```

and pass `inplace=true`; ScientificFitting validates that the callbacks fill every
output. On automatic-differentiation paths, `p`, `out`, and `J` may contain
non-`Float64` scalar types; mutating functions must not hard-code `Float64`
buffers.

Input observations and starting values are copied to `Float64` storage; models
must return finite values wherever the solver evaluates them.

### Parameter Control

`fit_model` and the likelihood wrappers accept:

| Keyword | Meaning |
|---|---|
| `p0` | Required complete starting vector. |
| `bounds=(lower, upper)` | Componentwise closed bounds; use `+/-Inf` for an open side. |
| `fixed_parameters` | Remove parameters from the optimizer with `FixedParameter`, `i => value`, or equivalent named tuples. |
| `parameter_priors` | Independent normalized Gaussian or split-normal terms. |
| `parameter_constraints` | Correlated Gaussian terms on selected parameters. |
| `constraints` | General nonlinear constraints; `ineq(p) <= 0` and `eq(p) == 0`. |

A `FixedParameter` uncertainty is report metadata, never an objective term; see
[`FixedParameter`](@ref) and
[Fixed Parameters And Bounds](statistics.md#Fixed-Parameters-And-Bounds).

### Solver Control

| Keyword | Default | Contract |
|---|---:|---|
| `solver` | `nothing` | Optional [`OptimizationSolver`](@ref) or [`NativeMinuitSolver`](@ref); leave legacy `backend`/`optimizer` at `:auto`. |
| `maxiters` | `500` for `fit_model`, `1000` for likelihood wrappers | Positive per-candidate budget: iterations for LsqFit/Optim, objective calls for NLopt/NativeMinuit. |
| `tol` | `1e-10`, or `1e-6` with `derivatives=:finite` | Positive solver-specific stopping tolerance, not a statistical error; see [Solver Adapters](api_fitting.md#Solver-Adapters). |
| `initial_guesses` | `nothing` | Additional complete starting vectors. |
| `multistart` | `1` | Total candidate budget including `p0`. Additional `initial_guesses` are tried next, then deterministic candidates from bounds or scaled `p0` where available. Set this above 1 to use additional starts. |

`fit_model` additionally accepts:

| Keyword | Default | Contract |
|---|---:|---|
| `backend` | `:auto` | `:auto`, `:lsqfit`, or `:optimization`. |
| `cost` | `:auto` | `:chi2` or full `:gaussian_likelihood` on the ``-2\log L`` scale; `:auto` uses the latter for parameter-dependent covariance. |
| `scale_covariance` | `:auto` | `:auto`, `:always`, or `:never`; see [Parameter Covariance](@ref parameter-covariance-reference). |
| `jacobian` | `nothing` | Analytic model Jacobian, allocating or in-place according to `inplace`. |
| `x_derivative` | `nothing` | Vector ``\partial f/\partial x`` for efficient x-uncertainty propagation. |

Automatic backend routing is specified in
[Backend Design](backend_design.md); `result.backend` records the choice.

## Reference Sections

| Need | Reference |
|---|---|
| Fit inputs, uncertainty objects, constraints, likelihoods, and solver behavior | [Fitting](api_fitting.md) |
| Result fields, covariance, profiles, contours, diagnostics, and reports | [Results And Diagnostics](api_results.md) |
| Optional Makie boundary, fit figures, annotations, and visual contracts | [Fit Plotting](api_plotting.md) |
| Residual, pull, profile, contour, and profile-matrix figures | [Diagnostic Plotting](api_plotting_diagnostics.md) |
