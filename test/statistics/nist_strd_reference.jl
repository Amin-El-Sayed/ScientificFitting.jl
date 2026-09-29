using ScientificFitting
using Test

# NIST StRD nonlinear regression suite: certified parameter estimates,
# standard deviations, and residual sums of squares from an external national
# metrology reference. The fits are unweighted, so they also validate the
# profiled-sigma covariance scaling (chi2/ndf) against certified values.
include(joinpath(@__DIR__, "..", "data", "nist_strd.jl"))

# Certified values carry 11 significant digits; solver tolerances and problem
# conditioning limit what a general-purpose fit reproduces. The bounds below
# are per-problem observations, not targets: tighten only with evidence.
const NIST_RTOL = Dict(
    "Misra1a" => (params=1e-6, sd=1e-4, rss=1e-9),
    "Chwirut2" => (params=1e-6, sd=1e-4, rss=1e-9),
    "Gauss1" => (params=1e-6, sd=1e-4, rss=1e-9),
    "MGH17" => (params=1e-4, sd=1e-3, rss=1e-8),
    "Rat42" => (params=1e-6, sd=1e-4, rss=1e-9),
    "Thurber" => (params=1e-5, sd=1e-3, rss=1e-9),
    "BoxBOD" => (params=1e-6, sd=1e-4, rss=1e-9),
    "Eckerle4" => (params=1e-6, sd=1e-4, rss=1e-9),
)

@testset "NIST StRD certified references" begin
    for problem in NIST_PROBLEMS
        tol = NIST_RTOL[problem.name]
        @testset "$(problem.name)" begin
            result = fit_model(problem.model, problem.x, problem.y;
                               p0=problem.start2, maxiters=2000)
            @test result.converged
            @test isapprox(result.params, problem.certified; rtol=tol.params)
            @test isapprox(result.stats.chi2, problem.certified_rss; rtol=tol.rss)
            # NIST standard deviations use sigma_hat^2 = RSS/(n-p); with no
            # supplied uncertainties, scale_covariance=:auto applies exactly
            # that factor.
            @test isapprox(result.param_stderr, problem.certified_sd; rtol=tol.sd)
        end
    end
end
