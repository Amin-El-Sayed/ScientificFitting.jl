using JuFitter
using LinearAlgebra
using Test

@testset "Profile and contour statistical references" begin
    x = collect(range(-2.0, 2.0; length=21))
    sigma_y = @. 0.12 + 0.02 * (x + 2.0)
    model(x, p) = @. p[1] * x + p[2]
    y = model(x, [1.7, -0.4])

    result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y, scale_covariance=:never)
    covariance = result.param_covariance
    precision = Symmetric(covariance) \ Matrix{Float64}(I, size(covariance))
    slope_sigma = sqrt(covariance[1, 1])
    offset_sigma = sqrt(covariance[2, 2])

    @testset "One-dimensional profile follows covariance parabola" begin
        slope_values = result.params[1] .+ slope_sigma .* [-2.0, -1.0, 0.0, 1.0, 2.0]
        prof = profile(result, 1; values=slope_values, threshold=1.0)
        expected_delta = @. abs2(slope_values - result.params[1]) / covariance[1, 1]

        @test prof.parameter_index == 1
        @test prof.best_value == result.params[1]
        @test isapprox(prof.delta_cost, expected_delta; atol=2e-8, rtol=2e-8)
        @test argmin(prof.delta_cost) == 3
        @test !any(f -> f.code == :profile_not_parabolic, diagnose(prof; local_sigma=slope_sigma).findings)
    end

    @testset "Profile interval reproduces one-sigma covariance errors" begin
        slope_values = result.params[1] .+ slope_sigma .* [-2.0, -1.0, 0.0, 1.0, 2.0]
        interval = profile_interval(result, 1; values=slope_values, threshold=1.0)

        @test isapprox(interval.lower, result.params[1] - slope_sigma; atol=2e-8, rtol=2e-8)
        @test isapprox(interval.upper, result.params[1] + slope_sigma; atol=2e-8, rtol=2e-8)
        @test isapprox(interval.uncertainty_minus, slope_sigma; atol=2e-8, rtol=2e-8)
        @test isapprox(interval.uncertainty_plus, slope_sigma; atol=2e-8, rtol=2e-8)
    end

    @testset "Adaptive profile refines threshold crossings" begin
        slope_values = result.params[1] .+ slope_sigma .* [-2.0, 0.0, 2.0]
        coarse = profile(result, 1; values=slope_values, threshold=1.0)
        refined = profile(result, 1; values=slope_values, threshold=1.0, adaptive=true, max_refinements=3)
        lower, upper = JuFitter._profile_crossings(refined)

        @test length(refined.values) > length(coarse.values)
        @test issorted(refined.values)
        @test isapprox(lower, result.params[1] - slope_sigma; atol=5e-3, rtol=5e-3)
        @test isapprox(upper, result.params[1] + slope_sigma; atol=5e-3, rtol=5e-3)
    end

    @testset "Profile inputs fail clearly before refits" begin
        slope_values = result.params[1] .+ slope_sigma .* [-1.0, 0.0, 1.0]

        @test profile(result, 1; values=slope_values, npoints=1).values == sort(slope_values)
        @test_throws ArgumentError profile(result, 1; npoints=2)
        @test_throws ArgumentError profile(result, 1; nsigma=Inf)
        @test_throws ArgumentError profile(result, 1; threshold=0.0)
        @test_throws ArgumentError profile(result, 1; threshold=NaN)
        @test_throws ArgumentError profile(result, 1; values=[slope_values[1], NaN, slope_values[3]])
        @test_throws ArgumentError profile(result, 1; values=[slope_values[1], slope_values[1], slope_values[3]])
    end

    @testset "Two-dimensional contour follows covariance quadratic form" begin
        slope_values = result.params[1] .+ slope_sigma .* [-1.0, 0.0, 1.0]
        offset_values = result.params[2] .+ offset_sigma .* [-1.0, 0.0, 1.0]
        cont = contour(result, 1, 2; xvalues=slope_values, yvalues=offset_values, levels=[2.30])

        expected_delta = Matrix{Float64}(undef, length(slope_values), length(offset_values))
        for ix in eachindex(slope_values), iy in eachindex(offset_values)
            delta = [slope_values[ix] - result.params[1], offset_values[iy] - result.params[2]]
            expected_delta[ix, iy] = dot(delta, precision * delta)
        end

        @test cont.parameter_indices == (1, 2)
        @test cont.levels == [2.30]
        @test isapprox(cont.delta_cost, expected_delta; atol=3e-8, rtol=3e-8)
        @test argmin(vec(cont.delta_cost)) == 5
        @test !any(
            f -> f.code == :contour_not_elliptic,
            diagnose(cont; local_covariance=covariance, local_center=result.params[[1, 2]]).findings,
        )

        sorted_cont = contour(
            result,
            1,
            2;
            xvalues=slope_values,
            yvalues=offset_values,
            levels=[6.18, 2.30, 2.30],
        )
        @test sorted_cont.levels == [2.30, 6.18]
        @test_throws ArgumentError contour(result, 1, 2; levels=Float64[])
        @test_throws ArgumentError contour(result, 1, 2; levels=[0.0, 2.30])
        @test_throws ArgumentError contour(result, 1, 2; levels=[NaN, 2.30])
        @test contour(result, 1, 2; xvalues=slope_values, yvalues=offset_values, npoints=1).levels == [2.30, 6.18]
        @test_throws ArgumentError contour(result, 1, 2; npoints=1)
        @test_throws ArgumentError contour(result, 1, 2; nsigma=NaN)
        @test_throws ArgumentError contour(result, 1, 2; xvalues=[slope_values[1], NaN], yvalues=offset_values)
        @test_throws ArgumentError contour(result, 1, 2; xvalues=[slope_values[1]], yvalues=offset_values)
        @test_throws ArgumentError contour(result, 1, 2; xvalues=slope_values, yvalues=[offset_values[1], offset_values[1]])
    end

    @testset "Adaptive contour refines level-crossing cells" begin
        slope_values = result.params[1] .+ slope_sigma .* [-2.0, 0.0, 2.0]
        offset_values = result.params[2] .+ offset_sigma .* [-2.0, 0.0, 2.0]
        coarse = contour(result, 1, 2; xvalues=slope_values, yvalues=offset_values, levels=[2.30])
        refined = contour(
            result,
            1,
            2;
            xvalues=slope_values,
            yvalues=offset_values,
            levels=[2.30],
            adaptive=true,
            max_refinements=2,
        )

        @test length(refined.x_values) > length(coarse.x_values)
        @test length(refined.y_values) > length(coarse.y_values)
        @test issorted(refined.x_values)
        @test issorted(refined.y_values)

        finite_delta = refined.delta_cost[isfinite.(refined.delta_cost)]
        @test minimum(finite_delta) <= 2.30 <= maximum(finite_delta)
        @test !any(
            f -> f.code == :contour_not_elliptic,
            diagnose(refined; local_covariance=covariance, local_center=result.params[[1, 2]]).findings,
        )
    end

    @testset "Profile diagnosis catches non-parabolic and unbracketed scans" begin
        values = collect(range(-2.0, 2.0; length=9))
        skewed_delta = [v < 0 ? v^2 : 0.35 * v^2 for v in values]
        skewed = ProfileResult(1, values, skewed_delta, skewed_delta, 1.0, 0.0)
        skewed_report = diagnose(skewed; local_sigma=1.0, tolerance=0.25)

        @test any(f -> f.code == :profile_not_parabolic, skewed_report.findings)

        narrow_values = [-0.4, 0.0, 0.4]
        narrow_delta = narrow_values .^ 2
        narrow = ProfileResult(1, narrow_values, narrow_delta, narrow_delta, 1.0, 0.0)
        narrow_report = diagnose(narrow; local_sigma=1.0)

        @test any(f -> f.code == :profile_threshold_not_bracketed, narrow_report.findings)

        failed = ProfileResult(1, values, [0.0, 0.5, Inf, 1.2, 2.0, Inf, 2.0, 1.2, 0.5], [0.0, 0.5, Inf, 1.2, 2.0, Inf, 2.0, 1.2, 0.5], 1.0, 0.0)
        failed_report = diagnose(failed; local_sigma=1.0)

        @test any(f -> f.code == :profile_refit_failed, failed_report.findings)
    end

    @testset "Contour diagnosis catches non-elliptic and unbracketed scans" begin
        values = collect(range(-2.0, 2.0; length=9))
        elliptic_delta = Matrix{Float64}(undef, length(values), length(values))
        warped_delta = Matrix{Float64}(undef, length(values), length(values))
        for ix in eachindex(values), iy in eachindex(values)
            xvalue = values[ix]
            yvalue = values[iy]
            elliptic_delta[ix, iy] = xvalue^2 + yvalue^2
            warped_delta[ix, iy] = (xvalue < 0 ? xvalue^2 : 0.25 * xvalue^2) + yvalue^2
        end

        warped = ContourResult((1, 2), values, values, warped_delta, warped_delta, [2.30])
        warped_report = diagnose(warped; local_covariance=Matrix{Float64}(I, 2, 2), local_center=[0.0, 0.0], tolerance=0.5)

        @test any(f -> f.code == :contour_not_elliptic, warped_report.findings)

        narrow_values = [-0.4, 0.0, 0.4]
        narrow_delta = Matrix{Float64}(undef, length(narrow_values), length(narrow_values))
        for ix in eachindex(narrow_values), iy in eachindex(narrow_values)
            narrow_delta[ix, iy] = narrow_values[ix]^2 + narrow_values[iy]^2
        end
        narrow = ContourResult((1, 2), narrow_values, narrow_values, narrow_delta, narrow_delta, [2.30])
        narrow_report = diagnose(narrow; local_covariance=Matrix{Float64}(I, 2, 2), local_center=[0.0, 0.0])

        @test any(f -> f.code == :contour_levels_not_bracketed, narrow_report.findings)

        failed_delta = copy(elliptic_delta)
        failed_delta[2, 3] = Inf
        failed_delta[4, 5] = NaN
        failed = ContourResult((1, 2), values, values, failed_delta, failed_delta, [2.30])
        failed_report = diagnose(failed; local_covariance=Matrix{Float64}(I, 2, 2), local_center=[0.0, 0.0])

        @test any(f -> f.code == :contour_refit_failed, failed_report.findings)
    end
end
