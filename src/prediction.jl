"""
    predict(result::FitResult, x=result.problem.x; uncertainty=false)

Evaluate the fitted model on finite 1D coordinates `x`, without a plotting
backend. With `uncertainty=true`, return `(mean, sigma)`, where `sigma` is the
local standard uncertainty of the fitted mean: `sqrt(diag(J * Cov(p) * J'))`.
It excludes new-observation noise and is not a predictive interval for future
measurements. Nonlinear/asymmetric uncertainty may require profile inference.

The model Jacobian follows `result.problem.derivatives`; an explicitly supplied
analytic Jacobian takes precedence. Invalid local covariance raises an error
instead of silently drawing a zero-width band.
"""
function predict(result::FitResult, x::AbstractVector=result.problem.x; uncertainty::Bool=false)
    grid = _float_vector(x)
    _assert_finite_vector("prediction coordinates", grid)
    mean = _model_values(result.problem, result.params; x=grid)
    return uncertainty ? (mean=mean, sigma=_prediction_band_sigma(result, grid)) : mean
end

"""Propagate parameter covariance for every renderer using one numerical implementation."""
function _prediction_band_sigma(result::FitResult, xgrid::AbstractVector)
    J = _parameter_jacobian(result.problem, result.params; x=xgrid)
    cov = result.param_covariance
    all(isfinite, cov) || throw(ArgumentError("prediction uncertainty requires finite parameter covariance"))
    terms = (J * cov) .* J
    variances = vec(sum(terms; dims=2))
    # Only absorb rounding-scale cancellation, not negative inferred variances.
    tolerance = 100eps(Float64) .* vec(sum(abs, terms; dims=2))
    all(isfinite, variances) && all(variances .>= -tolerance) || throw(ArgumentError(
        "prediction uncertainty requires a valid local parameter covariance",
    ))
    return sqrt.(max.(variances, 0.0))
end
