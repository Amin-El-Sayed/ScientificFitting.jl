"""
Typed boundary for foreign callbacks evaluated with finite differences.

The function field is deliberately abstract: the solver specializes on its
declared return type, not each runtime-created Python closure. One dynamic
dispatch per callback lets these fits share precompiled numerical code. Native
Julia models bypass this adapter and retain automatic differentiation.
"""
struct _TypedCallback{R}
    f::Function
end
(callback::_TypedCallback{R})(args...) where {R} = callback.f(args...)::R

"""Validate the differentiation policy stored with a problem and its refits."""
function _validate_derivatives(mode::Symbol)
    mode in (:auto, :finite) || throw(ArgumentError("derivatives must be :auto or :finite"))
    return mode
end

"""Default stopping tolerance; differenced gradients have a higher numerical noise floor."""
_default_fit_tolerance(mode::Symbol) = mode == :finite ? 1e-6 : 1e-10

"""Select SciML's derivative provider without passing dual numbers to foreign callbacks."""
function _optimization_ad(problem; second_order::Bool=false)
    ad = _derivative_mode(problem) == :finite ? AutoFiniteDiff(fdtype=Val(:central)) : AutoForwardDiff()
    return second_order ? DifferentiationInterface.SecondOrder(ad, ad) : ad
end

"""Differentiate vector predictions/residuals with the same policy as the objective."""
function _derivative_jacobian(problem, f, p)
    return _derivative_mode(problem) == :finite ?
           FiniteDiff.finite_difference_jacobian(f, p, Val(:central)) :
           ForwardDiff.jacobian(f, p)
end

"""Differentiate the complete cost, including covariance normalization and constraints."""
function _derivative_hessian(problem, f, p)
    return _derivative_mode(problem) == :finite ?
           FiniteDiff.finite_difference_hessian(f, p) : ForwardDiff.hessian(f, p)
end
