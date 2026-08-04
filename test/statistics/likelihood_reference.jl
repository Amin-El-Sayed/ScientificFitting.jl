using Distributions
using JuFitter
using LinearAlgebra
using SpecialFunctions
using Statistics
using Test

function _poisson_minus2loglik(counts, mu)
    return 2.0 * sum(n -> mu - n * log(mu) + loggamma(n + 1.0), counts)
end

function _poisson_deviance_reference(counts, mu)
    total = 0.0
    for n in counts
        total += n == 0 ? 2.0 * mu : 2.0 * (mu - n + n * log(n / mu))
    end
    return total
end

@testset "Likelihood statistical references" begin
    @testset "Constant Poisson model matches analytic MLE" begin
        x = collect(1.0:8.0)
        counts = Float64[2, 4, 3, 7, 6, 5, 9, 4]
        model(x, p) = fill(p[1], length(x))

        result = fit_poisson_model(model, x, counts; p0=[3.0], bounds=([1e-9], [100.0]))
        muhat = mean(counts)
        expected_cost = _poisson_minus2loglik(counts, muhat)
        expected_gof = _poisson_deviance_reference(counts, muhat)
        ndf = length(counts) - 1

        @test result.converged
        @test result.stats.cost == :poisson_likelihood
        @test isapprox(result.params[1], muhat; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_covariance[1, 1], muhat / length(counts); atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_stderr[1], sqrt(muhat / length(counts)); atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.cost_min, expected_cost; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.minus2loglik_min, expected_cost; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.chi2, expected_gof; atol=2e-8, rtol=2e-8)
        @test result.stats.ndf == ndf
        @test isapprox(result.stats.chi2_ndf, expected_gof / ndf; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.pvalue, ccdf(Chisq(ndf), expected_gof); atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.aic, expected_cost + 2.0; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.bic, expected_cost + log(length(counts)); atol=2e-8, rtol=2e-8)
    end

    @testset "Poisson goodness of fit includes auxiliary Gaussian data" begin
        x = collect(1.0:6.0)
        counts = Float64[2, 3, 4, 2, 5, 3]
        prior_mean = 5.0
        prior_sigma = 0.8
        model(x, p) = fill(p[1], length(x))

        result = fit_poisson_model(
            model,
            x,
            counts;
            p0=[3.0],
            bounds=([1e-9], [100.0]),
            parameter_priors=(index=1, mean=prior_mean, sigma=prior_sigma),
        )
        muhat = result.params[1]
        expected_gof = _poisson_deviance_reference(counts, muhat) +
                       abs2((muhat - prior_mean) / prior_sigma)
        ndf = length(counts) + 1 - 1

        @test result.converged
        @test isapprox(result.stats.chi2, expected_gof; atol=2e-8, rtol=2e-8)
        @test result.stats.ndf == ndf
        @test isapprox(result.stats.chi2_ndf, expected_gof / ndf; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.pvalue, ccdf(Chisq(ndf), expected_gof); atol=2e-8, rtol=2e-8)
    end

    @testset "Poisson goodness of fit includes correlated parameter data" begin
        x = collect(range(-1.0, 1.5; length=7))
        counts = Float64[2, 3, 4, 6, 8, 11, 15]
        constraint_mean = [log(4.5), 0.55]
        constraint_cov = [0.09 0.015; 0.015 0.04]
        model(x, p) = @. exp(p[1] + p[2] * x)

        result = fit_poisson_model(
            model,
            x,
            counts;
            p0=[log(4.0), 0.5],
            parameter_constraints=(
                indices=[1, 2],
                mean=constraint_mean,
                covariance=constraint_cov,
            ),
        )
        expected_counts = model(x, result.params)
        data_gof = sum(zip(counts, expected_counts)) do (n, mu)
            n == 0 ? 2.0 * mu : 2.0 * (mu - n + n * log(n / mu))
        end
        delta = result.params .- constraint_mean
        expected_gof = data_gof + dot(delta, constraint_cov \ delta)
        ndf = length(counts) + length(constraint_mean) - length(result.params)

        @test result.converged
        @test isapprox(result.stats.chi2, expected_gof; atol=2e-8, rtol=2e-8)
        @test result.stats.ndf == ndf
        @test isapprox(result.stats.pvalue, ccdf(Chisq(ndf), expected_gof); atol=2e-8, rtol=2e-8)
    end

    @testset "Unbinned normal location fit matches analytic MLE" begin
        data = Float64[-1.1, -0.4, 0.2, 0.5, 0.9, 1.4]
        sigma = 0.7
        normal_location_pdf(x, p) = exp(-0.5 * abs2((x - p[1]) / sigma)) / (sigma * sqrt(2pi))

        result = fit_unbinned_model(normal_location_pdf, data; p0=[0.0])
        muhat = mean(data)
        expected_cost = length(data) * log(2pi * sigma^2) + sum(abs2, (data .- muhat) ./ sigma)
        expected_cov = sigma^2 / length(data)

        @test result.converged
        @test result.stats.cost == :unbinned_likelihood
        @test isapprox(result.params[1], muhat; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_covariance[1, 1], expected_cov; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_stderr[1], sqrt(expected_cov); atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.cost_min, expected_cost; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.minus2loglik_min, expected_cost; atol=2e-8, rtol=2e-8)
        @test isnan(result.stats.chi2)
        @test isnan(result.stats.chi2_ndf)
        @test isnan(result.stats.pvalue)
        @test isapprox(result.stats.aic, expected_cost + 2.0; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.bic, expected_cost + log(length(data)); atol=2e-8, rtol=2e-8)
    end

    @testset "Constant histogram Poisson model matches analytic MLE" begin
        edges = collect(0.0:1.0:6.0)
        counts = Float64[3, 5, 2, 7, 4, 6]
        expected_counts(edges, p) = fill(p[1], length(edges) - 1)

        result = fit_histogram_model(expected_counts, edges, counts; p0=[4.0], bounds=([1e-9], [100.0]))
        muhat = mean(counts)
        expected_cost = _poisson_minus2loglik(counts, muhat)
        expected_gof = _poisson_deviance_reference(counts, muhat)
        ndf = length(counts) - 1

        @test result.converged
        @test result.stats.cost == :histogram_poisson_likelihood
        @test isapprox(result.params[1], muhat; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_covariance[1, 1], muhat / length(counts); atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.cost_min, expected_cost; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.chi2, expected_gof; atol=2e-8, rtol=2e-8)
        @test result.stats.ndf == ndf
        @test isapprox(result.stats.pvalue, ccdf(Chisq(ndf), expected_gof); atol=2e-8, rtol=2e-8)
    end

    @testset "Extended unbinned homogeneous rate matches analytic MLE" begin
        data = Float64[0.1, 0.4, 1.2, 1.8, 2.1, 2.7]
        domain = (0.0, 3.0)
        rate(x, p) = exp(p[1])

        result = fit_extended_unbinned_model(rate, data, domain; p0=[0.0])
        n = length(data)
        length_domain = domain[2] - domain[1]
        lambda_hat = n / length_domain
        expected_param = log(lambda_hat)
        expected_cost = 2.0 * n - 2.0 * n * log(lambda_hat)
        expected_cov = 1.0 / n

        @test result.converged
        @test result.stats.cost == :extended_unbinned_likelihood
        @test isapprox(result.params[1], expected_param; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_covariance[1, 1], expected_cov; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.cost_min, expected_cost; atol=2e-8, rtol=2e-8)
        @test isnan(result.stats.chi2)
        @test isnan(result.stats.pvalue)
        @test isapprox(result.stats.aic, expected_cost + 2.0; atol=2e-8, rtol=2e-8)
        @test isapprox(result.stats.bic, expected_cost + log(n); atol=2e-8, rtol=2e-8)
    end

    @testset "Mapped multi-fit matches stacked weighted least squares" begin
        x1 = collect(range(0.0, 3.0; length=8))
        x2 = collect(range(0.5, 3.5; length=7))
        sigma1 = fill(0.06, length(x1))
        sigma2 = fill(0.08, length(x2))
        y1 = @. 1.8 * x1 + 0.4 + 0.02 * sin(x1)
        y2 = @. 1.8 * x2 - 0.7 + 0.02 * cos(x2)
        local_linear(x, p) = @. p[1] * x + p[2]

        X = vcat(
            hcat(x1, ones(length(x1)), zeros(length(x1))),
            hcat(x2, zeros(length(x2)), ones(length(x2))),
        )
        y = vcat(y1, y2)
        sigma = vcat(sigma1, sigma2)
        w = 1.0 ./ abs2.(sigma)
        fisher = X' * (w .* X)
        rhs = X' * (w .* y)
        params = Symmetric(fisher) \ rhs
        cov = Symmetric(fisher) \ Matrix{Float64}(I, size(fisher))
        residuals = y .- X * params
        chi2 = sum(abs2, residuals ./ sigma)
        ndf = length(y) - length(params)

        result = fit_multi_model(
            [local_linear, local_linear],
            [x1, x2],
            [y1, y2];
            p0=[1.0, 0.0, 0.0],
            sigma_y=[sigma1, sigma2],
            parameter_map=[[1, 2], [1, 3]],
        )

        @test result.converged
        @test result.stats.cost == :multi_chi2
        @test isapprox(result.params, params; atol=2e-8, rtol=2e-8)
        @test isapprox(result.param_covariance, cov; atol=2e-7, rtol=2e-7)
        @test isapprox(result.stats.chi2, chi2; atol=2e-8, rtol=2e-8)
        @test result.stats.ndf == ndf
        @test isapprox(result.stats.pvalue, ccdf(Chisq(ndf), chi2); atol=2e-8, rtol=2e-8)
    end
end
