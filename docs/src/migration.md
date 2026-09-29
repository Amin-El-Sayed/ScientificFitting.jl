# Migrating From Other Tools

Concept-for-concept translations for the most common starting points.

## From scipy.optimize.curve_fit

| scipy | ScientificFitting |
| --- | --- |
| `curve_fit(f, x, y, p0)` | `fit_model(f, x, y; p0)` — the model takes the parameter vector: `f(x, p)` |
| `sigma=σ` (1-D) | `sigma_y=σ` |
| `sigma=Σ` (2-D) | `cov_y=Σ` |
| `absolute_sigma=True` | `scale_covariance=:never` |
| `absolute_sigma=False` (default) | `scale_covariance=:auto` scales by ``\chi^2/\mathrm{ndf}`` exactly when no uncertainties were supplied; with supplied `sigma_y` it does **not** rescale — request `:always` explicitly for scipy's default behavior |
| `bounds=(lo, hi)` | `bounds=(lo, hi)` — vectors with one entry per parameter; scalars do not broadcast |
| `pcov` (second return value, `popt, pcov = curve_fit(…)`) | `result.param_covariance`; `np.sqrt(np.diag(pcov))` corresponds to `result.param_stderr` |
| — | goodness of fit, p-value, diagnostics, and profile intervals come with every fit |

## From scipy.odr

`scipy.odr` is deprecated since SciPy 1.17 and scheduled for removal in 1.19.

| scipy.odr | ScientificFitting |
| --- | --- |
| `Model(f)` with `f(beta, x)` | the model takes `f(x, p)` — argument order is flipped |
| `RealData(x, y, sx=…, sy=…)` | `sigma_x=…, sigma_y=…` (standard deviations, as in odr) |
| `ODR(data, model, beta0=p0).run()` | `fit_model(f, x, y; p0, sigma_x, sigma_y)` |
| `output.beta` | `result.params` |
| `output.sd_beta` | `result.param_stderr` — odr always rescales by ``\chi^2/\mathrm{ndf}`` (`res_var`); pass `scale_covariance=:always` to reproduce that convention |

Fitted values differ slightly by construction: odr omits the log-determinant
term of the effective-variance objective; see
[Cross-Checking Against ODR-Convention Tools](@ref) below.

## From lmfit

| lmfit | ScientificFitting |
| --- | --- |
| `Parameters()` with `vary=False` | `fixed_parameters=[FixedParameter(i, value)]` |
| `Parameter(min=…, max=…)` | `bounds` |
| algebraic parameter constraints | nonlinear `constraints=(eq=…, ineq=…)` |
| `conf_interval` (F-test) | `profile_interval(result, i)` — likelihood-ratio crossing on the same scale as `param_stderr` ([Covariance Scaling](statistics.md#Covariance-Scaling)) |
| `fit_report()` | `report_text(result)`; structured access via `fit_report(result)` and `diagnose(result)` |
| model composition with prefixes | not built in; compose Julia functions, or see [BuildConstructors](interfaces.md) |

## From iminuit

| iminuit | ScientificFitting |
| --- | --- |
| `Minuit(cost, …)` with a builtin cost | the matching entry point: `fit_model`, `fit_poisson_model`, `fit_unbinned_model`, … |
| custom cost function | `fit_custom(objective; p0, nobs)` on the ``-2\log L`` scale |
| `m.migrad()` | minimization runs inside every `fit_*` call with the default solver; for the MIGRAD algorithm itself, `solver=NativeMinuitSolver()` (Julia ≥ 1.11, [extension](interfaces.md)) |
| HESSE errors | `result.param_stderr` / `param_covariance` |
| MINOS asymmetric errors | `profile_interval(result, i)`; unbracketed sides come back `NaN` with a finding, and `diagnose` reports failed refits |
| `NormalConstraint` | `parameter_priors` (scalar, also asymmetric) and `parameter_constraints` (correlated) |
| `fixed` parameters | `fixed_parameters` |

## From LsqFit.jl

| LsqFit | ScientificFitting |
| --- | --- |
| `curve_fit(model, x, y, p0)` | `fit_model(model, x, y; p0)` — same model signature; LsqFit remains the automatic fast path for least squares with parameter-independent weights ([Backend Design](backend_design.md)) |
| `wt` weights | `sigma_y` (standard deviations, not weights) |
| `stderror(fit)` | `result.param_stderr` |
| `margin_error`, `confint` | `profile_interval` for likelihood-ratio intervals |

## Cross-Checking Against ODR-Convention Tools

`scipy.odr`, kafe2, and York-style line fits use the effective-variance
objective **without** its log-determinant term; ScientificFitting keeps the
term because the effective covariance is parameter-dependent. Expect small,
reproducible differences in x-error fits — neither convention is
systematically better in the checked configurations, and with `sigma_x` alone
the determinant term is required for a well-posed problem. Derivation and
Monte-Carlo numbers:
[Statistics Reference](statistics.md#Relation-To-ODR-And-The-York-Method).
