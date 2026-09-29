using ScientificFitting
using Test
using LinearAlgebra
using Random
using Statistics

# References for the covariance-scale contract (audit blockers B1-B3) and the
# sample-size calibration of the residual-structure diagnostics (D1). Every
# numerical comparison is against a closed-form ordinary-least-squares result,
# a hand-written likelihood normalization, or a rate calibrated by
# construction - never against the package's own output.

@testset "Covariance scaling consistency references" begin
    linear_model(x, p) = @. p[1] + p[2] * x

    @testset "Unweighted OLS matches the closed-form scaled solution" begin
        rng = Xoshiro(7)
        n = 25
        x = collect(range(0.0, 5.0; length=n))
        y = 1.0 .+ 2.0 .* x .+ 0.3 .* randn(rng, n)

        result = fit_model(linear_model, x, y; p0=[0.0, 1.0])

        # Closed-form OLS: beta = (X'X)^-1 X'y, sigma_hat^2 = RSS/(n-2),
        # stderr_j = sigma_hat * sqrt([(X'X)^-1]_jj).
        X = hcat(ones(n), x)
        beta = (X' * X) \ (X' * y)
        rss = sum(abs2, y .- X * beta)
        sigma2_hat = rss / (n - 2)
        cov_ref = sigma2_hat .* inv(X' * X)

        @test isapprox(result.params, beta; rtol=1e-8)
        @test isapprox(result.stats.chi2, rss; rtol=1e-8)
        @test isapprox(result.param_stderr, sqrt.(diag(cov_ref)); rtol=1e-6)
        @test isapprox(result.param_covariance, cov_ref; rtol=1e-6)

        # B1: the profile interval must agree with the reported standard error.
        # Both uncertainty fields are magnitudes relative to the best value.
        for index in 1:2
            interval = profile_interval(result, index)
            @test isapprox(interval.uncertainty_plus, result.param_stderr[index]; rtol=2e-4)
            @test isapprox(interval.uncertainty_minus, result.param_stderr[index]; rtol=2e-4)
        end

        # B1 consequence: an exactly linear model is exactly parabolic on the
        # reported scale - no false :profile_not_parabolic.
        prof = ScientificFitting.profile(result, 1)
        report = diagnose(prof; local_sigma=result.param_stderr[1])
        @test !any(f -> f.code == :profile_not_parabolic, report.findings)

        # Contours share the scale: the delta surface must match the local
        # covariance ellipse of the reported (scaled) covariance.
        pair = ScientificFitting.contour(result, 1, 2; levels=[2.30])
        pair_report = diagnose(
            pair;
            local_covariance=result.param_covariance,
            local_center=result.params,
        )
        @test !any(f -> f.code == :contour_not_elliptic, pair_report.findings)

        matrix = profile_matrix(result; npoints_profile=9, npoints_contour=5)
        @test all(status -> status == :ok, values(matrix.panel_status))
    end

    @testset "Profiled-sigma likelihood normalization is unit invariant (B2)" begin
        rng = Xoshiro(11)
        n = 40
        x = collect(range(0.0, 10.0; length=n))
        y = 1.0 .+ 2.0 .* x .+ 0.02 .* x .^ 2 .+ 0.5 .* randn(rng, n)
        quadratic_model(x, p) = @. p[1] + p[2] * x + p[3] * x^2

        fit_lin = fit_model(linear_model, x, y; p0=[1.0, 2.0])

        # Hand formula: -2 log L with sigma profiled out, sigma counted in AIC.
        rss = fit_lin.stats.chi2
        m2ll_ref = n * (log(2pi * rss / n) + 1)
        @test isapprox(fit_lin.stats.minus2loglik_min, m2ll_ref; rtol=1e-10)
        @test isapprox(fit_lin.stats.aic, m2ll_ref + 2 * (2 + 1); rtol=1e-10)
        @test isapprox(fit_lin.stats.bic, m2ll_ref + log(n) * (2 + 1); rtol=1e-10)

        # The model comparison must not depend on the unit of y.
        fit_quad = fit_model(quadratic_model, x, y; p0=[1.0, 2.0, 0.0])
        daic = fit_quad.stats.aic - fit_lin.stats.aic
        fit_lin_mv = fit_model(linear_model, x, 1000.0 .* y; p0=[1000.0, 2000.0])
        fit_quad_mv = fit_model(quadratic_model, x, 1000.0 .* y; p0=[1000.0, 2000.0, 0.0])
        daic_mv = fit_quad_mv.stats.aic - fit_lin_mv.stats.aic
        @test isapprox(daic, daic_mv; atol=1e-6)

        # Weighted control: with known sigma the normalization keeps its
        # explicit log-determinant and sigma is not counted as a parameter.
        sigma = fill(0.5, n)
        weighted = fit_model(linear_model, x, y; p0=[1.0, 2.0], sigma_y=sigma)
        m2ll_weighted = n * log(2pi) + sum(log.(sigma .^ 2)) + weighted.stats.chi2
        @test isapprox(weighted.stats.minus2loglik_min, m2ll_weighted; rtol=1e-10)
        @test isapprox(weighted.stats.aic, m2ll_weighted + 2 * 2; rtol=1e-10)
    end

    @testset "scale_covariance is honest in the likelihood branch (B3)" begin
        rng = Xoshiro(13)
        n = 30
        x = collect(range(0.0, 5.0; length=n))
        y = 1.0 .+ 2.0 .* x .+ 0.1 .* randn(rng, n)
        sigma_y = fill(0.1, n)
        sigma_x = fill(0.05, n)

        @test_throws ArgumentError fit_model(
            linear_model, x, y;
            p0=[1.0, 2.0], sigma_y=sigma_y, sigma_x=sigma_x,
            scale_covariance=:always,
        )

        auto = fit_model(linear_model, x, y; p0=[1.0, 2.0], sigma_y=sigma_y, sigma_x=sigma_x)
        never = fit_model(linear_model, x, y; p0=[1.0, 2.0], sigma_y=sigma_y, sigma_x=sigma_x,
            scale_covariance=:never)
        @test auto.stats.cost == :gaussian_likelihood
        @test auto.param_stderr == never.param_stderr
    end

    @testset "Fixed parameters keep a unit correlation diagonal (D5)" begin
        rng = Xoshiro(17)
        n = 20
        x = collect(range(0.0, 5.0; length=n))
        y = 1.0 .+ 2.0 .* x .+ 0.1 .* randn(rng, n)
        result = fit_model(linear_model, x, y; p0=[1.0, 2.0], sigma_y=fill(0.1, n),
            fixed_parameters=[FixedParameter(1, 1.0)])
        @test result.param_correlation[1, 1] == 1.0
        @test result.param_correlation[2, 2] == 1.0
        @test result.param_correlation[1, 2] == 0.0
    end

    @testset "Overparametrized fits report NaN information criteria (D7)" begin
        result = fit_model((x, p) -> @.(p[1] + p[2] * x + p[3] * x^2),
            [0.0, 1.0], [1.0, 3.0]; p0=[1.0, 1.0, 0.0])
        @test result.stats.ndf < 0
        @test isnan(result.stats.aic)
        @test isnan(result.stats.bic)
    end

    @testset "Explicit initial guesses are always tried (D3)" begin
        problem = FitProblem((x, p) -> @.(p[1] * x + p[2]),
            [1.0, 2.0, 3.0], [2.0, 4.0, 6.0]; p0=[1.0, 0.0])
        guess = [5.0, -1.0]
        candidates = ScientificFitting._initial_candidates(problem, [guess], 1)
        @test length(candidates) == 2
        @test candidates[1] == [1.0, 0.0]
        @test candidates[2] == guess
    end

    @testset "Residual-structure diagnostics are calibrated across n (D1)" begin
        rng = Xoshiro(20250925)
        structure_codes = (:structured_residual_signs, :autocorrelated_pulls, :large_pull, :extreme_pull)
        repetitions = 300
        for n in (15, 25, 50, 200, 500)
            counts = Dict(code => 0 for code in structure_codes)
            x = collect(range(0.0, 5.0; length=n))
            sigma = fill(0.2, n)
            for _ in 1:repetitions
                y = 1.0 .+ 2.0 .* x .+ 0.2 .* randn(rng, n)
                result = fit_model(linear_model, x, y; p0=[1.0, 2.0], sigma_y=sigma)
                # The residual-structure checks run in diagnose(result), not in
                # the stored numerical diagnostics.
                for finding in diagnose(result).findings
                    finding.code in structure_codes || continue
                    counts[finding.code] += 1
                end
            end
            for code in structure_codes
                @test counts[code] / repetitions <= 0.08
            end
        end

        # Power control: genuine structure must still be flagged.
        x = collect(range(-3.0, 3.0; length=41))
        curved = @. 0.4 * x^2 + 0.8 * x - 1.0
        biased = fit_model(linear_model, x, curved; p0=[0.0, 0.0], sigma_y=fill(0.05, 41))
        @test any(f -> f.code == :structured_residual_signs, diagnose(biased).findings)
    end
end
