# API Reference

This reference defines JuFitter's public contracts. Start with
[Quickstart](quickstart.md) for a first fit or
[Statistical Foundations](statistical_foundations.md) to choose the method
chapter that matches the observation process.

## Choose An Entry Point

The fitting function is selected by the observation model, not by the plotting
style or optimizer.

| Data and sampling model | Entry point | Model contract | Result |
|---|---|---|---|
| Numeric ``x`` and ``y`` with Gaussian uncertainties | [`fit_model`](@ref) | `model(x, p) -> y_hat` | [`FitResult`](@ref) |
| Independent counts | [`fit_poisson_model`](@ref) | `model(x, p) -> expected_counts` | [`LikelihoodFitResult`](@ref) |
| Histogram with expected bin counts | [`fit_histogram_model`](@ref) | `expected_counts(edges, p) -> mu` | [`LikelihoodFitResult`](@ref) |
| Histogram from a normalized density | [`fit_histogram_density`](@ref) | `pdf(x, p) -> density` | [`LikelihoodFitResult`](@ref) |
| Independent unbinned observations | [`fit_unbinned_model`](@ref) | `pdf(x, p) -> density` | [`LikelihoodFitResult`](@ref) |
| Unbinned events with a parameter-dependent rate | [`fit_extended_unbinned_model`](@ref) | `rate(x, p) -> event_rate` | [`LikelihoodFitResult`](@ref) |
| Observations addressed by non-numeric indices | [`fit_indexed_model`](@ref) | `model(indices, p) -> y_hat` | [`LikelihoodFitResult`](@ref) |
| Several datasets sharing parameters | [`fit_multi_model`](@ref) | one `model_i(x_i, p_i)` per dataset | [`LikelihoodFitResult`](@ref) |
| A custom scalar objective | [`fit_custom`](@ref) | `objective(p) -> scalar` | [`LikelihoodFitResult`](@ref) |

Minimal call shapes, with required keywords shown explicitly:

| Entry point | Minimal call |
|---|---|
| `fit_model` | `fit_model(model, x, y; p0=[...], sigma_y=[...])` |
| `fit_poisson_model` | `fit_poisson_model(expected_counts, x, counts; p0=[...])` |
| `fit_histogram_model` | `fit_histogram_model(expected_per_bin, edges, counts; p0=[...])` |
| `fit_histogram_density` | `fit_histogram_density(pdf, edges, counts; p0=[...], total_count=sum(counts))` |
| `fit_unbinned_model` | `fit_unbinned_model(pdf, observations; p0=[...])` |
| `fit_extended_unbinned_model` | `fit_extended_unbinned_model(rate, observations, (a, b); p0=[...])` |
| `fit_indexed_model` | `fit_indexed_model(model, indices, y; p0=[...], cov_y=C)` |
| `fit_multi_model` | `fit_multi_model(models, xs, ys; p0=[...], sigma_y=sigma_sets)` |
| `fit_custom` | `fit_custom(cost; p0=[...], nobs=n)` |

These are signatures, not one shared executable example. Each model must obey
the contract in the first table, and every placeholder must be replaced by data
of matching dimensions. Complete analyses are in the [Gallery](gallery.md).

For reusable low-level workflows, construct [`FitProblem`](@ref) or
[`LikelihoodFitProblem`](@ref) and call [`fit`](@ref). `JuFitter.fit` extends
the `StatsAPI.fit` generic, so it coexists with StatsBase and Distributions.

## Common Conventions

### Parameters And Model Functions

`p0` fixes the parameter order. Parameter indices in bounds, fixed values,
constraints, profiles, and contours are Julia's one-based indices into that
vector.

The allocating model contract is:

```julia
model(x, p) -> vector with length(y)
```

For allocation-sensitive fits, use:

```julia
model!(out, x, p)
jacobian!(J, x, p)  # optional
```

and pass `inplace=true`. JuFitter validates that the callbacks fill every
output. On automatic-differentiation paths, `p`, `out`, and `J` may contain
non-`Float64` scalar types; mutating functions must not hard-code `Float64`
buffers internally.

Input observations and starting values are copied to `Float64` storage. A model
must return finite values at every parameter point used by the solver.

### Parameter Control

The following controls are accepted by `fit_model` and the likelihood wrappers.

| Keyword | Meaning |
|---|---|
| `p0` | Required complete starting vector. |
| `bounds=(lower, upper)` | Componentwise closed bounds; use `+/-Inf` for an open side. |
| `fixed_parameters` | Remove parameters from the optimizer with `FixedParameter`, `i => value`, or equivalent named tuples. |
| `parameter_priors` | Independent normalized Gaussian or split-normal terms. |
| `parameter_constraints` | Correlated Gaussian terms on selected parameters. |
| `constraints` | General nonlinear constraints; `ineq(p) <= 0` and `eq(p) == 0`. |

Constraint callbacks receive the complete parameter vector in `p0` order,
including fixed entries. JuFitter maps it to reduced free-parameter coordinates
internally.

`FixedParameter(..., sigma)` stores that uncertainty in reports and in the fixed
parameter's covariance diagonal, with zero fitted cross-covariances. It does not
change the objective or propagate uncertainty into free fitted parameters. Use
a `ParameterPrior` or `ParameterConstraint` when an external measurement must
participate in the fit.

### Solver Control

| Keyword | Default | Contract |
|---|---:|---|
| `maxiters` | `500` for `fit_model`, `1000` for likelihood wrappers | Positive iteration limit for each candidate. |
| `tol` | `1e-10` | Positive absolute and relative solver tolerance. |
| `initial_guesses` | `nothing` | Additional complete starting vectors. |
| `multistart` | `1` | Deterministic candidates generated from finite bounds or scaled versions of `p0`. |

`fit_model` additionally accepts:

| Keyword | Default | Contract |
|---|---:|---|
| `backend` | `:auto` | `:auto`, `:lsqfit`, or `:optimization`. |
| `cost` | `:auto` | `:chi2` or full `:gaussian_likelihood` on the ``-2\log L`` scale; `:auto` uses the latter for parameter-dependent covariance. |
| `scale_covariance` | `:auto` | `:auto`, `:always`, or `:never`; see [Parameter Covariance](@ref parameter-covariance-reference). |
| `jacobian` | `nothing` | Analytic model Jacobian, allocating or in-place according to `inplace`. |
| `x_derivative` | `nothing` | Vector ``\partial f/\partial x`` for efficient x-uncertainty propagation. |

`backend=:auto` uses LsqFit only when static chi-square least squares represents
the complete problem. An incompatible explicit `backend=:lsqfit` request raises
an error rather than dropping bounds, constraints, priors, or
parameter-dependent covariance.

## Reference Sections

| Need | Reference |
|---|---|
| Fit inputs, uncertainty objects, constraints, likelihoods, and solver behavior | [Fitting](api_fitting.md) |
| Result fields, covariance, profiles, contours, diagnostics, and reports | [Results And Diagnostics](api_results.md) |
| Optional Makie boundary, fit figures, annotations, and diagnostic plots | [Plotting](api_plotting.md) |
