
"""
    BootstrapSamples

Parameter vectors returned by [`parametric_bootstrap`](@ref). `params` has one
row per successful refit and one column per parameter in `p0` order; fixed
parameters keep their fixed value in every row. `n_requested` is the number of
simulated data sets and `n_failed` the number of refits that failed or did not
converge. Failed refits are not part of `params`, so
`size(params, 1) == n_requested - n_failed`.
"""
struct BootstrapSamples
    params::Matrix{Float64}
    n_requested::Int
    n_failed::Int
end

const _BOOTSTRAP_SUPPORTED_CASES =
    "Supported are Gaussian fits with sigma_y, cov_y, or no observation uncertainty and " *
    "without x uncertainty, error components, or WhiteningOperator; fit_poisson_model, " *
    "fit_histogram_model, and fit_histogram_density; and a LikelihoodFitProblem constructed " *
    "with bootstrap_resampler."

_bootstrap_unsupported(reason::AbstractString) = ArgumentError(
    "parametric resampling is not defined for $reason. " * _BOOTSTRAP_SUPPORTED_CASES)

_covariance_root(factor::Cholesky) = factor.L

# CHOLMOD factorizes P*C*P' = L*L' with a fill-reducing permutation P, so the
# root of C itself is P'*L: the rows of L returned to the original order.
_covariance_root(factor::SparseArrays.CHOLMOD.Factor) = sparse(factor.L)[invperm(factor.p), :]

function _bootstrap_resampler(result::FitResult)
    problem = result.problem
    reason = problem.whitening !== nothing ?
        "a WhiteningOperator, which does not expose a covariance root to simulate noise" :
        has_x_uncertainty(problem) ?
        "x uncertainties: resampling the latent true x values is a separate model" :
        any(component -> component.active, problem.error_components) ?
        "error components, whose covariance may depend on the observed y or on the model" : nothing
    reason === nothing || throw(_bootstrap_unsupported(reason))

    base_root = if problem.sigma_y !== nothing
        Diagonal(problem.sigma_y)
    elseif problem.cov_y !== nothing
        _covariance_root(_stable_cholesky(problem.cov_y))
    else
        I
    end
    # Resample with the observation model the fit actually reports: the same
    # chi2/ndf factor that scales param_covariance rescales the noise here, so
    # all three interval methods share one scale (unweighted fits resample at
    # the estimated residual standard deviation, :never keeps sigma = 1).
    scale = if base_root === I && result.stats.cost == :gaussian_likelihood
        # The profiled-sigma likelihood of an unweighted fit (its information
        # criteria section) estimates sigma^2 = RSS/n; the reported curvature
        # lives on that scale, so the resampled noise must too.
        rss = result.stats.chi2
        isfinite(rss) && rss > 0 || throw(_bootstrap_unsupported(
            "an unweighted likelihood fit whose residual sum of squares is " *
            "not positive (RSS = $rss)"))
        rss / length(problem.y)
    else
        if result.stats.cost == :chi2 &&
                _should_scale_covariance(problem, result.options.scale_covariance) &&
                !(result.stats.ndf > 0 && isfinite(result.stats.chi2) && result.stats.chi2 > 0)
            throw(_bootstrap_unsupported(
                "a fit whose residual scale sqrt(chi2/ndf) is required but undefined " *
                "(chi2 = $(result.stats.chi2), ndf = $(result.stats.ndf))"))
        end
        _covariance_scale(result)
    end
    root = sqrt(scale) * base_root
    n = length(problem.y)
    return (rng, p) -> _with_data(problem, _model_values(problem, p) .+ root * randn(rng, n), p)
end

function _bootstrap_resampler(result::LikelihoodFitResult)
    resampler = result.problem.bootstrap_resampler
    resampler === nothing && throw(_bootstrap_unsupported(
        "an objective without bootstrap_resampler, which does not define how to simulate its data"))
    return resampler
end

_refit(result::FitResult, problem::FitProblem) = fit(
    problem;
    cost=result.options.cost,
    maxiters=result.options.maxiters,
    tol=result.options.tol,
    scale_covariance=result.options.scale_covariance,
    solver=result.options.solver,
)

_refit(result::LikelihoodFitResult, problem::LikelihoodFitProblem) = fit(
    problem;
    maxiters=result.options.maxiters,
    tol=result.options.tol,
    parameter_covariance=result.options.parameter_covariance,
    solver=result.options.solver,
)

"""
    parametric_bootstrap(result; nsamples=999, rng=Random.default_rng())
        -> BootstrapSamples

Simulate `nsamples` data sets from the fitted model at `p̂ = result.params`,
refit each one with the options of the original fit starting at `p̂`, and
collect the fitted parameter vectors in a [`BootstrapSamples`](@ref). Only the
observations are simulated; the model, uncertainties, bounds, constraints,
fixed parameters, and inputs such as `total_count` of `fit_histogram_density`
keep their values from `result.problem`. `multistart` and `initial_guesses` of
the original fit are not used: every refit has the single start `p̂`.

Simulated data by problem type:
- Gaussian `FitProblem` with a static observation covariance from `sigma_y`
  or `cov_y`: `y* = model(x, p̂) + L*z` with independent standard normal `z`,
  where `L = Diagonal(sigma_y)` or the Cholesky factor of `cov_y`, multiplied
  by the square root of the same `chi2/ndf` factor that scales
  `param_covariance` ([Covariance Scaling](@ref)). Unweighted fits therefore
  resample at the estimated residual standard deviation, and
  `scale_covariance=:never` keeps the stated `sigma = 1`; when the factor is
  required but undefined (`ndf <= 0`, or `chi2` not finite and positive), the
  call raises `ArgumentError`.
- [`fit_poisson_model`](@ref), [`fit_histogram_model`](@ref), and
  [`fit_histogram_density`](@ref): independent counts `n*_i ~ Poisson(mu_i(p̂))`
  per observation or bin, with `mu_i` computed exactly as in the objective.
- A `LikelihoodFitProblem` constructed with `bootstrap_resampler`.

Every other problem raises `ArgumentError`: x uncertainties (resampling the
latent true x values is a separate model), error components, a
`WhiteningOperator`, parameter priors and parameter constraints (their
auxiliary measurements would stay fixed while the data vary, which narrows the
distribution), and objectives without a resampler such as `fit_custom` or the
unbinned fits.

A refit that throws or returns `converged == false` is dropped and counted in
`n_failed`. As soon as more than 10% of the `nsamples` refits have failed, an
`ErrorException` reports the counts instead of returning a distribution biased
toward data sets that are easy to fit. `nsamples < 1` raises `DomainError`. The
samples are reproducible for a given state of `rng`.
"""
function parametric_bootstrap(
    result::Union{FitResult, LikelihoodFitResult};
    nsamples::Integer=999,
    rng::AbstractRNG=Random.default_rng(),
)
    nsamples >= 1 || throw(DomainError(nsamples, "nsamples must be >= 1"))
    problem = result.problem
    isempty(problem.parameter_priors) && isempty(problem.parameter_constraints) ||
        throw(_bootstrap_unsupported("parameter priors or parameter constraints, whose " *
            "auxiliary measurements would stay fixed while the data vary"))
    resample = _bootstrap_resampler(result)

    phat = result.params
    params = Matrix{Float64}(undef, nsamples, length(phat))
    kept = 0
    n_failed = 0
    for _ in 1:nsamples
        resampled = resample(rng, phat)
        refit = try
            _refit(result, resampled)
        catch err
            err isa InterruptException && rethrow()
            err
        end
        if refit isa Exception || !refit.converged
            n_failed += 1
            10 * n_failed > nsamples && throw(ErrorException(
                "at least $n_failed of $nsamples bootstrap refits failed or did not converge, " *
                "more than the allowed 10%; the remaining samples would favor data sets that " *
                "are easy to fit. Last failure: " *
                _bootstrap_failure_text(refit)))
        else
            kept += 1
            params[kept, :] .= refit.params
        end
    end
    return BootstrapSamples(params[1:kept, :], nsamples, n_failed)
end

# Failure summaries must stay readable: a solver DomainError can carry whole
# dual-number vectors through showerror.
function _bootstrap_failure_text(refit)
    text = refit isa Exception ? sprint(showerror, refit) : refit.message
    line = first(split(text, '\n'))
    return length(line) > 200 ? first(line, 200) * "…" : line
end
