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

@testset "Multistart ranks status before cost and keeps the best fallback" begin
    record(converged, cost) = (; converged, stats=(cost_min=cost,))
    prefer = ScientificFitting._prefer_fit
    @test prefer(record(true, 1.), record(false, 0.))
    @test !prefer(record(false, 0.), record(true, 1.))
    for converged in (false, true)
        @test prefer(record(converged, 1.), nothing)
        @test prefer(record(converged, 1.), record(converged, 2.))
        @test !prefer(record(converged, 2.), record(converged, 1.))
        @test !prefer(record(converged, 1.), record(converged, 1.))
        for bad in (Inf, NaN)
            @test !prefer(record(converged, bad), nothing)
            @test !prefer(record(converged, bad), record(false, 1.))
        end
    end

    # Neither run may move with a one-evaluation budget; the second is better.
    custom = fit_custom(p -> (p[1]-2)^2; p0=[0.], initial_guesses=[[1.]],
        multistart=2, nobs=10, optimizer=:nelder_mead, maxiters=1)
    @test !custom.converged
    @test custom.params == [1.]
    @test custom.stats.cost_min == 1.

    # Exercise the Gaussian loop as well, using independent stopped-run results.
    x, y = collect(1.:8.), fill(2., 8)
    model(x, p) = fill(exp(p[1]), length(x))
    runs = [fit_model(model, x, y; p0=[start], sigma_y=ones(8), maxiters=1)
            for start in (-2., -1.)]
    @test all(r -> !r.converged, runs)
    @test runs[2].stats.cost_min < runs[1].stats.cost_min
    combined = fit_model(model, x, y; p0=[-2.], sigma_y=ones(8), maxiters=1,
                         initial_guesses=[[-1.]], multistart=2)
    @test !combined.converged
    @test combined.params == runs[2].params
end
