using JuFitter
using LinearAlgebra
using SparseArrays
using Test

@testset "In-place model and Jacobian contracts" begin
    x = collect(range(-2.0, 3.0; length=200))
    sigma_y = @. 0.08 + 0.01 * (x + 2.0)
    y = @. 1.35 * x - 0.42 + 0.03 * sin(2.1 * x)

    model(x, p) = @. p[1] * x + p[2]
    model_calls = Ref(0)
    function model!(out, x, p)
        model_calls[] += 1
        @. out = p[1] * x + p[2]
        return nothing
    end
    jacobian_calls = Ref(0)
    function jacobian!(J, x, p)
        jacobian_calls[] += 1
        J[:, 1] .= x
        J[:, 2] .= 1
        return nothing
    end

    reference = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y)
    inplace = fit_model(model!, x, y; p0=[1.0, 0.0], sigma_y=sigma_y, inplace=true)
    analytic = fit_model(
        model!,
        x,
        y;
        p0=[1.0, 0.0],
        sigma_y=sigma_y,
        jacobian=jacobian!,
        inplace=true,
    )
    bounded = fit_model(
        model!,
        x,
        y;
        p0=[1.0, 0.0],
        sigma_y=sigma_y,
        bounds=([-10.0, -10.0], [10.0, 10.0]),
        inplace=true,
    )
    cov_y = Matrix(Diagonal(sigma_y .^ 2))
    covariance_reference = fit_model(model, x, y; p0=[1.0, 0.0], cov_y=cov_y)
    covariance_inplace = fit_model(model!, x, y; p0=[1.0, 0.0], cov_y=cov_y, inplace=true)
    sparse_cov_y = spdiagm(
        -1 => fill(0.0005, length(x) - 1),
        0 => sigma_y .^ 2,
        1 => fill(0.0005, length(x) - 1),
    )
    sparse_reference = fit_model(model, x, y; p0=[1.0, 0.0], cov_y=sparse_cov_y)
    sparse_inplace = fit_model(model!, x, y; p0=[1.0, 0.0], cov_y=sparse_cov_y, inplace=true)
    fixed = fit_model(
        model!,
        x,
        y;
        p0=[1.0, reference.params[2]],
        sigma_y=sigma_y,
        fixed_parameters=2 => reference.params[2],
        jacobian=jacobian!,
        inplace=true,
    )
    scan = profile(analytic, 1; npoints=7, nsigma=1.5, adaptive=false)

    @test inplace.backend == :lsqfit
    @test analytic.backend == :lsqfit
    @test bounded.backend == :optimization
    @test model_calls[] > 0
    @test jacobian_calls[] > 0
    @test isapprox(inplace.params, reference.params; rtol=1e-9, atol=1e-10)
    @test isapprox(analytic.params, reference.params; rtol=1e-9, atol=1e-10)
    @test isapprox(bounded.params, reference.params; rtol=1e-7, atol=1e-8)
    @test isapprox(analytic.param_covariance, reference.param_covariance; rtol=1e-8)
    @test isapprox(covariance_inplace.params, covariance_reference.params; rtol=1e-9, atol=1e-10)
    @test isapprox(sparse_inplace.params, sparse_reference.params; rtol=1e-9, atol=1e-10)
    @test isapprox(fixed.params[1], reference.params[1]; rtol=1e-9, atol=1e-10)
    @test scan isa ProfileResult
    @test all(isfinite, scan.delta_cost)
    @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, 0.0], inplace=true)
    @test_throws ArgumentError fit_model(
        model!,
        x,
        y;
        p0=[1.0, 0.0],
        jacobian=(x, p) -> hcat(x, ones(length(x))),
        inplace=true,
    )
    @test_throws ArgumentError fit_model(
        (out, x, p) -> (out[1] = p[1] * x[1] + p[2]),
        x,
        y;
        p0=[1.0, 0.0],
        inplace=true,
    )
    @test_throws ArgumentError fit_model(
        model!,
        x,
        y;
        p0=[1.0, 0.0],
        jacobian=(J, x, p) -> (J[:, 1] .= x),
        inplace=true,
    )
end
