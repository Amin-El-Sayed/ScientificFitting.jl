module TuringLikelihoodReference

using ScientificFitting, Distributions, ForwardDiff, Random, Statistics, Test
using Turing
using FlexiChains: ess, rhat, mcse, Extra

"""Attach a single explicit prior to an SF data-only -2log(L) callback."""
@model function rate_posterior(data_cost, prior)
    rate ~ prior
    @addlogprob! -data_cost([rate])/2
end

@testset "Turing external likelihood contract" begin
    counts = [0, 3, 1, 4, 2, 5, 0, 2]
    x = collect(eachindex(counts))
    rate_model(x, p) = fill(p[1], length(x))
    options = (; p0=[2.], bounds=([0.], [Inf]))
    ordinary = fit_poisson_model(rate_model, x, counts; options...)
    distribution = fit_distribution(p -> Poisson(p[1]), counts; options...)
    calibrated = fit_poisson_model(rate_model, x, counts; options...,
        parameter_priors=[ParameterPrior(1, 5., .7)])
    @test ordinary.converged && distribution.converged && calibrated.converged
    @test only(ordinary.params) ≈ mean(counts) atol=1e-7
    @test only(distribution.params) ≈ mean(counts) atol=1e-7
    @test only(calibrated.params) > only(ordinary.params)

    prior = Gamma(2., 3.)  # Distributions uses shape and scale, not rate.
    exact = Gamma(2 + sum(counts), inv(1/3 + length(counts)))
    # The SF auxiliary term must not leak into a separately specified Turing prior.
    for result in (ordinary, distribution, calibrated)
        posterior = rate_posterior(result.problem.objective, prior)
        offset = logjoint(posterior, (rate=2.,)) - logpdf(exact, 2.)
        for rate in (.03, .5, 2., 10.)
            target(r) = logjoint(posterior, (rate=r,))
            @test loglikelihood(posterior, (; rate)) ≈ loglikelihood(Poisson(rate), counts)
            @test logprior(posterior, (; rate)) ≈ logpdf(prior, rate)
            @test target(rate) - logpdf(exact, rate) ≈ offset
            @test ForwardDiff.derivative(target, rate) ≈
                (1 + sum(counts))/rate - (1/3 + length(counts))
            @test ForwardDiff.derivative(r -> ForwardDiff.derivative(target, r), rate) ≈
                -(1 + sum(counts))/rate^2
        end
    end

    # Four chains exercise Turing's positive-parameter transform and actual NUTS.
    posterior = rate_posterior(distribution.problem.objective, prior)
    chain = sample(Xoshiro(20260912), posterior, NUTS(500, .85; adtype=AutoForwardDiff()),
        MCMCSerial(), 2000, 4; progress=false, verbose=false)
    draws = chain[@varname(rate)]
    @test size(draws) == (2000, 4)
    @test all(>(0), draws)
    @test rhat(chain)[@varname(rate)] < 1.01
    @test ess(chain; kind=:bulk)[@varname(rate)] > 1000
    @test ess(chain; kind=:tail)[@varname(rate)] > 1000
    @test !any(chain[Extra(:numerical_error)])
    @test abs(mean(draws) - mean(exact)) < 6mcse(chain)[@varname(rate)]
    @test std(draws) ≈ std(exact) rtol=.07
    for q in (.16, .5, .84)
        @test abs(quantile(vec(Array(draws)), q) - quantile(exact, q)) < .1std(exact)
    end
end

end
