using JuFitter
using LinearAlgebra
using Test

@testset "Torture input validation and pathological fits" begin
    x = collect(range(0.0, 1.0; length=8))
    y = @. 2.0 * x + 1.0
    model(x, p) = @. p[1] * x + p[2]

    @testset "Unphysical data and uncertainties fail before optimization" begin
        @test_throws ArgumentError fit_model(model, Float64[], Float64[]; p0=[1.0, 0.0])
        @test_throws ArgumentError fit_model(model, [0.0, NaN], [1.0, 2.0]; p0=[1.0, 0.0])
        @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, Inf])
        @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=fill(0.0, length(x)))
        @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, 0.0], sigma_x=[0.1, -0.1, fill(0.1, length(x) - 2)...])
    end

    @testset "Bad covariance matrices are not silently repaired" begin
        nonsymmetric = Matrix{Float64}(I, length(x), length(x))
        nonsymmetric[1, 2] = 0.2
        singular = ones(length(x), length(x))

        @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, 0.0], cov_y=nonsymmetric)
        @test_throws ArgumentError fit_model(model, x, y; p0=[1.0, 0.0], cov_y=singular)
    end

    @testset "Fixed and profiled parameters cannot bypass bounds" begin
        bounded = fit_model(
            model,
            x,
            y;
            p0=[1.5, 1.0],
            sigma_y=fill(0.1, length(x)),
            bounds=([0.0, -Inf], [2.0, Inf]),
        )

        @test bounded.converged
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.5, 1.0],
            sigma_y=fill(0.1, length(x)),
            bounds=([0.0, -Inf], [2.0, Inf]),
            fixed_parameters=(index=1, value=3.0),
        )

        prof = profile(bounded, 1; values=[1.0, 2.0, 3.0], on_failure=:inf)
        @test isfinite(prof.cost_values[1])
        @test isfinite(prof.cost_values[2])
        @test isinf(prof.cost_values[3])
        @test any(f -> f.code == :profile_refit_failed, diagnose(prof).findings)
        @test_throws ArgumentError profile(bounded, 1; values=[1.0, 2.0, 3.0], on_failure=:throw)
    end

    @testset "User start values are validated, not silently repaired" begin
        bounds = ([0.0, -Inf], [2.0, Inf])

        @test fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            bounds=bounds,
            initial_guesses=[[1.5, 0.0]],
            multistart=2,
        ).converged

        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[3.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            bounds=bounds,
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            bounds=bounds,
            initial_guesses=[[NaN, 0.0]],
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            bounds=bounds,
            initial_guesses=[[3.0, 0.0]],
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            bounds=bounds,
            initial_guesses=[[1.0]],
        )
    end

    @testset "Parameter priors and constraints reject invalid metadata" begin
        valid_constraint = (
            indices=[1, 2],
            mean=[2.0, 1.0],
            covariance=[0.25 0.02; 0.02 0.09],
        )
        @test fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            parameter_priors=(index=1, mean=2.0, sigma=0.5),
            parameter_constraints=valid_constraint,
        ).converged

        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            parameter_priors=(index=1, mean=NaN, sigma=0.5),
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            parameter_priors=(index=1, mean=2.0, sigma=Inf),
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            parameter_constraints=(indices=[1, 2], mean=[2.0, Inf], covariance=[0.25 0.02; 0.02 0.09]),
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            parameter_constraints=(indices=[1, 2], mean=[2.0, 1.0], covariance=[1.0 0.0; 0.4 1.0]),
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            parameter_constraints=(indices=[1, 2], mean=[2.0, 1.0], covariance=[1.0 NaN; NaN 1.0]),
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            fixed_parameters=(index=1, value=Inf),
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            fixed_parameters=(index=1, value=1.0, sigma=Inf),
        )
        @test_throws ArgumentError fit_custom(
            p -> sum(abs2, p);
            p0=[1.0],
            nobs=3,
            parameter_priors=(index=1, mean=0.0, sigma=Inf),
        )
    end

    @testset "Likelihood observations fail before optimizer internals" begin
        count_x = collect(1.0:4.0)
        counts = [5.0, 4.0, 3.0, 2.0]
        count_model(x, p) = fill(exp(p[1]), length(x))
        expected_counts(edges, p) = fill(exp(p[1]), length(edges) - 1)
        pdf(x, p) = exp(-0.5 * (x - p[1])^2) / sqrt(2pi)
        rate(x, p) = exp(p[1])

        @test fit_poisson_model(count_model, count_x, counts; p0=[log(3.5)]).converged
        @test fit_histogram_model(expected_counts, 0.0:1.0:4.0, counts; p0=[log(3.5)]).converged
        @test fit_unbinned_model(pdf, [-0.2, 0.0, 0.4]; p0=[0.1]).converged
        @test fit_extended_unbinned_model(rate, [0.2, 0.4, 0.8], (0.0, 1.0); p0=[0.0]).converged

        @test_throws ArgumentError fit_poisson_model(count_model, [1.0, NaN], [1.0, 2.0]; p0=[0.0])
        @test_throws ArgumentError fit_poisson_model(count_model, count_x, [5.0, NaN, 3.0, 2.0]; p0=[0.0])
        @test_throws ArgumentError fit_poisson_model(count_model, count_x, [5.0, 4.5, 3.0, 2.0]; p0=[0.0])
        @test_throws ArgumentError fit_poisson_model((x, p) -> fill(exp(p[1]), length(x) - 1), count_x, counts; p0=[0.0])
        @test_throws ArgumentError fit_histogram_model(expected_counts, [0.0, 1.0, NaN, 4.0, 5.0], counts; p0=[0.0])
        @test_throws ArgumentError fit_histogram_model(expected_counts, [0.0, 1.0, 1.0, 3.0, 4.0], counts; p0=[0.0])
        @test_throws ArgumentError fit_histogram_model(expected_counts, 0.0:1.0:4.0, [5.0, 4.0, 3.5, 2.0]; p0=[0.0])
        @test_throws ArgumentError fit_histogram_density(pdf, 0.0:1.0:4.0, counts; p0=[0.0], total_count=Inf)
        @test_throws ArgumentError fit_histogram_density(pdf, 0.0:1.0:4.0, counts; p0=[0.0], rtol=NaN)
        @test_throws ArgumentError fit_unbinned_model(pdf, [-0.2, NaN, 0.4]; p0=[0.0])
        @test_throws ArgumentError fit_extended_unbinned_model(rate, [0.2, NaN], (0.0, 1.0); p0=[0.0])
        @test_throws ArgumentError fit_extended_unbinned_model(rate, [0.2, 1.4], (0.0, 1.0); p0=[0.0])
        @test_throws ArgumentError fit_extended_unbinned_model(rate, [0.2, 0.4], (0.0, Inf); p0=[0.0])
        @test_throws ArgumentError fit_extended_unbinned_model(rate, [0.2, 0.4], (0.0, 1.0); p0=[0.0], rtol=0.0)
    end

    @testset "Non-finite model output is a model error, not optimizer noise" begin
        bad_model(x, p) = [xi == 0.0 ? NaN : p[1] / xi + p[2] for xi in x]
        @test_throws ArgumentError fit_model(bad_model, x, y; p0=[1.0, 0.0], sigma_y=fill(0.1, length(x)))
    end

    @testset "Large badly scaled but identifiable linear fit remains accurate" begin
        n = 20_000
        xs = collect(range(-1.0e6, 1.0e6; length=n))
        true_p = [2.5e-6, -3.0e3]
        ys = @. true_p[1] * xs + true_p[2] + 0.02 * sin(xs / 1.0e5)
        sigma = fill(0.05, n)
        scaled_model(x, p) = @. p[1] * x + p[2]
        jacobian(x, p) = hcat(x, ones(length(x)))
        design = hcat(xs, ones(n))
        reference = (design' * design) \ (design' * ys)

        result = fit_model(
            scaled_model,
            xs,
            ys;
            p0=[1.0e-6, -2.5e3],
            sigma_y=sigma,
            jacobian=jacobian,
            scale_covariance=:never,
            maxiters=200,
        )

        @test result.converged
        @test isapprox(result.params[1], reference[1]; atol=2e-12, rtol=2e-8)
        @test isapprox(result.params[2], reference[2]; atol=2e-8, rtol=2e-11)
        @test result.stats.ndf == n - 2
        @test isfinite(result.stats.chi2)
        @test result.diagnostics.covariance_condition < 1e25
    end
end
