using ScientificFitting
using LsqFit
using Test

@testset "Solver limits and convergence are not fabricated" begin
    x = collect(1.0:8.0)
    y = [1.9, 2.1, 2.0, 1.8, 2.2, 2.1, 1.9, 2.0]
    model(x, p) = fill(exp(p[1]), length(x))
    jacobian(x, p) = fill(exp(p[1]), length(x), 1)
    model!(out, x, p) = (out .= exp(p[1]); nothing)
    jacobian!(out, x, p) = (out .= exp(p[1]); nothing)

    for inplace in (false, true), analytic in (false, true)
        function_used = inplace ? model! : model
        jac = analytic ? (inplace ? jacobian! : jacobian) : nothing
        reference = analytic ?
            LsqFit.curve_fit(function_used, jac, x, y, [-2.0];
                inplace, maxIter=1, x_tol=1e-10, g_tol=1e-10) :
            LsqFit.curve_fit(function_used, x, y, [-2.0];
                inplace, maxIter=1, x_tol=1e-10, g_tol=1e-10)
        short = fit_model(function_used, x, y; p0=[-2.0], sigma_y=ones(8),
            inplace, jacobian=jac, maxiters=1, tol=1e-10)
        @test !LsqFit.isconverged(reference)
        @test !short.converged
        @test short.params ≈ LsqFit.coef(reference) atol=1e-12
        @test occursin("did not converge", short.message)
        @test any(f -> f.code == :optimizer_not_converged, diagnose(short).findings)
        complete = fit_model(function_used, x, y; p0=[-2.0], sigma_y=ones(8),
            inplace, jacobian=jac, maxiters=100, tol=1e-10)
        @test complete.converged
        @test complete.params[1] ≈ log(2) atol=1e-8
    end

    # A finite objective at a stopped nuisance fit is not a profile minimum.
    objective(p) = p[1]^2 + (p[2]^2 - p[1] - 1)^2
    result = fit_custom(objective; p0=[0.0, 1.0], nobs=20, maxiters=1)
    @test result.converged
    nuisance = ScientificFitting._refit_with_fixed(result, [FixedParameter(1, 0.5)])
    @test !nuisance.converged
    @test isfinite(nuisance.stats.cost_min)
    scan = profile(result, 1; values=[-0.5, 0.0, 0.5])
    @test isinf(scan.delta_cost[3])
    @test any(f -> f.code == :profile_refit_failed, diagnose(scan).findings)
    @test_throws ErrorException profile(result, 1; values=[-0.5, 0.0, 0.5], on_failure=:throw)
end
