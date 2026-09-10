using ScientificFitting
using LinearAlgebra
using Test
import NativeMinuit

@testset "NativeMinuit objective, covariance and free coordinates" begin
    covariance = [0.04 0.006; 0.006 0.09]
    precision = inv(covariance)
    center = [1.5, -0.7]
    objective(p) = dot(p .- center, precision * (p .- center))
    for derivatives in (:auto, :finite)
        solver = NativeMinuitSolver(steps=[0.2, 0.3], strategy=2)
        r = fit_custom(objective; p0=[0., 0.], nobs=10, derivatives, solver, tol=1e-6)
        @test r.converged
        @test r.backend == :native_minuit
        @test r.params ≈ center atol=1e-5
        @test r.param_covariance ≈ covariance rtol=1e-5
        @test r.options.solver === solver
        @test ismissing(r.iterations)
        native = r.solver_result.raw
        @test native.errordef == 1
        @test native.valid
        @test native.fval ≈ r.stats.cost_min atol=1e-12
        NativeMinuit.hesse!(native)
        @test native.covariance ≈ r.param_covariance rtol=1e-4
    end

    cost(p) = (p[1]-2)^2 + (p[3]-p[1])^2 + p[2]^2
    r = fit_custom(cost; p0=[0.5, 0., 0.5], nobs=10,
        fixed_parameters=[FixedParameter(2, 0.4)],
        parameter_priors=[ParameterPrior(3, 2., 1.)],
        bounds=([0., 0., 0.], [4., 1., 4.]),
        solver=NativeMinuitSolver(steps=[0.1, 0.2, 0.3]), tol=1e-6)
    @test r.converged
    @test r.params ≈ [2., 0.4, 2.] atol=1e-4
    @test r.solver_result.parameter_indices == [1, 3]
    @test r.param_stderr[2] == 0
    refit = ScientificFitting._refit_with_fixed(r, [FixedParameter(1, 1.)])
    @test refit.converged
    @test refit.options.solver === r.options.solver
    @test refit.params ≈ [1., 0.4, 1.5] atol=1e-4
    @test refit.solver_result.parameter_indices == [3]
    scan = profile(r, 1; values=[1., 2., 3.], on_failure=:throw)
    @test scan.delta_cost ≈ [1.5, 0., 1.5] atol=1e-5

    x, y = [0., 1., 2., 3.], [1.1, 2.9, 5.2, 6.8]
    line(x, p) = @. p[1]*x + p[2]
    g = fit_model(line, x, y; p0=[1., 0.], sigma_y=fill(0.2, 4),
                  solver=NativeMinuitSolver(), tol=1e-6)
    baseline = fit_model(line, x, y; p0=[1., 0.], sigma_y=fill(0.2, 4))
    @test g.converged
    @test g.params ≈ baseline.params atol=1e-5
    @test g.param_covariance ≈ baseline.param_covariance rtol=1e-5
    @test g.stats.cost_min ≈ baseline.stats.cost_min atol=1e-7
    fixed = ScientificFitting._refit_with_fixed(g, [FixedParameter(1, 2.)])
    @test fixed.backend == :native_minuit
    @test fixed.converged
    @test fixed.options.solver === g.options.solver
end

@testset "NativeMinuit controls never silently change the problem" begin
    solver = NativeMinuitSolver()
    cost(p) = (p[1]-1)^2 + 100*(p[2]-p[1]^2)^2
    r = fit_custom(cost; p0=[-1.2, 1.], nobs=10, maxiters=1,
                   solver, parameter_covariance=:none)
    @test !r.converged
    @test r.solver_result.raw.fmin.internal.reached_call_limit
    @test occursin("call_limit=true", r.message)
    @test all(isnan, r.param_stderr)

    r = fit_custom(p -> (p[1]-2)^2; p0=[0.5], nobs=10,
                   bounds=([0.], [1.]), solver, tol=1e-6)
    @test r.converged
    @test 1-1e-5 <= r.params[1] <= 1
    @test_throws ArgumentError fit_custom(cost; p0=[0., 0.], nobs=10,
        constraints=(ineq=p -> [sum(p)-2],), solver)
    @test_throws ArgumentError NativeMinuitSolver(errordef=0.5)
    @test_throws ArgumentError NativeMinuitSolver(grad=identity)
    @test_throws ArgumentError NativeMinuitSolver(steps=[1., 0.])
    @test_throws ArgumentError fit_custom(cost; p0=[0., 0.], nobs=10,
        solver=NativeMinuitSolver(steps=[1.]))
end

@testset "Asymmetric profile interval agrees with native MINOS" begin
    # Four events constrain a positive mean; the coupled nuisance must be refit.
    cost(p) = 2*(p[1] - 4log(p[1])) + (p[2] - p[1]/4)^2 / 0.04
    r = fit_custom(cost; p0=[3., 0.5], nobs=10,
        bounds=([0.05, -Inf], [20., Inf]), solver=NativeMinuitSolver(), tol=1e-6)
    @test r.converged
    @test r.params ≈ [4., 1.] atol=1e-4
    native = r.solver_result.raw
    NativeMinuit.minos!(native, "p1"; tol=1e-6)
    err = native.merrors["p1"]
    @test err.lower_valid && err.upper_valid
    @test !err.lower_par_limit && !err.upper_par_limit
    @test !err.lower_fcn_limit && !err.upper_fcn_limit
    interval = profile_interval(r, 1; values=collect(range(1., 8.; length=81)))
    @test interval.lower ≈ r.params[1] + err.lower atol=3e-3
    @test interval.upper ≈ r.params[1] + err.upper atol=3e-3
    @test interval.uncertainty_plus > interval.uncertainty_minus + 0.5
    for crossing in (interval.lower, interval.upper)
        @test cost([crossing, crossing/4]) - r.stats.cost_min ≈ 1 atol=2e-3
    end
end
