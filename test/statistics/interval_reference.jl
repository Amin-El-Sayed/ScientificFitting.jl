using Distributions
using LinearAlgebra
using Random
using ScientificFitting
using Statistics
using Test

_interval_line(x, p) = @. p[1] * x + p[2]
_interval_constant(x, p) = fill(p[1], length(x))

# Weighted least squares for the straight line in closed form:
# p̂ = (XᵀWX)⁻¹XᵀWy and Cov(p̂) = (XᵀWX)⁻¹ with W = diag(1/sigma²).
function _interval_line_reference(x, y, sigma)
    X = hcat(x, ones(length(x)))
    w = 1.0 ./ abs2.(sigma)
    fisher = Symmetric(X' * (w .* X))
    params = fisher \ (X' * (w .* y))
    cov = inv(fisher)
    rss = sum(abs2, y .- X * params)
    return (; params, cov, rss)
end

function _interval_argument_error(f, fragment::AbstractString)
    err = try
        f()
        nothing
    catch caught
        caught
    end
    return err isa ArgumentError && occursin(fragment, err.msg)
end

# Sample quantiles of n draws from a density f are asymptotically normal
# (Bahadur): Var(ξ̂_q) = q(1 - q) / (n f(ξ_q)²) and, for q₁ < q₂,
# Cov(ξ̂_q₁, ξ̂_q₂) = q₁(1 - q₂) / (n f(ξ_q₁) f(ξ_q₂)). For N(μ, sd²) the density at
# μ ± z·sd is φ(z)/sd. With q = (1 - level)/2 and z = Φ⁻¹((1 + level)/2):
#   SE(bound) = sd * sqrt(q(1 - q) / n) / φ(z)
#   SE(width) = sd * sqrt(level(1 - level) / n) / φ(z)
# (the width variance is [2q(1 - q) - 2q²] / (n f²) = level(1 - level) / (n f²)).
# Assertions use four standard errors.
_interval_bound_tolerance(sd, level, n) =
    4 * sd * sqrt((1 - level) / 2 * (1 + level) / 2 / n) / pdf(Normal(), quantile(Normal(), (1 + level) / 2))
_interval_width_tolerance(sd, level, n) =
    4 * sd * sqrt(level * (1 - level) / n) / pdf(Normal(), quantile(Normal(), (1 + level) / 2))

_interval_width(interval::ConfidenceInterval) = interval.upper - interval.lower

@testset "Confidence interval references" begin
    x = collect(range(-2.0, 3.0; length=15))
    sigma_y = @. 0.2 + 0.03 * (x + 2.0)
    y = @. 2.0 * x + 1.5 + sigma_y * sin(3.1 * x + 0.4)
    # scale_covariance=:never states the given sigma_y as known, so the exact
    # sampling distribution is p̂ ~ N(p, (XᵀWX)⁻¹).
    known = fit_model(_interval_line, x, y; p0=[1.0, 0.0], sigma_y=sigma_y, scale_covariance=:never)
    reference = _interval_line_reference(x, y, sigma_y)
    default_level = 0.6827

    @testset "Linear Gaussian with known sigma: hessian and profile are exact" begin
        @test known.converged
        # LsqFit stops once its step falls below x_tol·‖p̂‖ with the default
        # x_tol = 1e-10 (default_fit_tolerance); rtol 1e-9 allows an order of
        # magnitude for the error remaining after that step.
        @test isapprox(known.params, reference.params; rtol=1e-9)
        # The 2-sigma level 2Φ(2) - 1 has the exact half-width 2 sd, written
        # independently of the quantile function used by the implementation.
        levels = ((default_level, quantile(Normal(), (1 + default_level) / 2)), (2 * cdf(Normal(), 2.0) - 1, 2.0))
        # :hessian is centered on the fitted p̂ by definition, so its limits are
        # compared around known.params; the half-width z·sd is exact for a
        # linear model up to rounding.
        # Profile limits come from linear interpolation of delta_cost = t² (t in
        # units of sd) on the final grid. confidence_interval scans ±nsigma = 5 sd
        # with npoints = 121 and halves every bracketing step max_refinements = 3
        # times, so the bracket width is h = (10 / 120) / 2³ = 1/96 sd. The chord
        # exceeds the parabola by at most h²/4 and the slope at the crossing is
        # 2z, so the limit moves by at most h²/(8z) sd; the refit costs add at
        # most the solver tolerance 1e-10 relative.
        h = (2 * 5 / (121 - 1)) / 2^3
        for (level, z) in levels, index in 1:2
            sd = sqrt(reference.cov[index, index])
            lower = reference.params[index] - z * sd
            upper = reference.params[index] + z * sd

            hessian = confidence_interval(known, index; method=:hessian, level=level)
            @test hessian.method == :hessian
            @test hessian.parameter_index == index
            @test hessian.level == level
            @test isapprox(hessian.lower, known.params[index] - z * sd; rtol=1e-10)
            @test isapprox(hessian.upper, known.params[index] + z * sd; rtol=1e-10)

            profiled = confidence_interval(known, index; method=:profile, level=level)
            @test profiled.method == :profile
            @test profiled.level == level
            @test abs(profiled.lower - lower) <= sd * h^2 / (8 * z) + 1e-10 * abs(lower)
            @test abs(profiled.upper - upper) <= sd * h^2 / (8 * z) + 1e-10 * abs(upper)
        end
        @test confidence_interval(known, 1).method == :profile
        @test confidence_interval(known, 1).level == default_level
    end

    @testset "Bootstrap of a linear Gaussian fit matches the exact normal quantiles" begin
        # The refit is linear in y* = ŷ + Diagonal(sigma_y)·z, so
        # p̂* = p̂ + (XᵀWX)⁻¹XᵀW·Diagonal(sigma_y)·z ~ N(p̂, (XᵀWX)⁻¹) exactly; only
        # the percentile estimate carries Monte Carlo error.
        nsamples = 20_000
        z = quantile(Normal(), (1 + default_level) / 2)
        for index in 1:2
            sd = sqrt(reference.cov[index, index])
            interval = confidence_interval(known, index; method=:bootstrap, nsamples=nsamples, rng=Xoshiro(41 + index))
            tolerance = _interval_bound_tolerance(sd, default_level, nsamples)
            @test interval.method == :bootstrap
            @test interval.parameter_index == index
            @test abs(interval.lower - (reference.params[index] - z * sd)) <= tolerance
            @test abs(interval.upper - (reference.params[index] + z * sd)) <= tolerance
        end
    end

    @testset "Bootstrap resamples on the reported covariance scale" begin
        # Unweighted line: :auto scales the covariance by s² = RSS/(n - 2), so
        # the hessian width is 2z·s·sd₀ with sd₀ = sqrt((XᵀX)⁻¹ᵢᵢ). Resampling
        # y* = ŷ + s·z gives p̂* ~ N(p̂, s²(XᵀX)⁻¹), the same width. With
        # scale_covariance=:never the stated sigma = 1 applies to both.
        y_unweighted = @. 2.0 * x + 1.5 + 0.15 * sin(3.1 * x + 0.4)
        scaled = fit_model(_interval_line, x, y_unweighted; p0=[1.0, 0.0])
        unscaled = fit_model(_interval_line, x, y_unweighted; p0=[1.0, 0.0], scale_covariance=:never)
        unit = _interval_line_reference(x, y_unweighted, ones(length(x)))
        s = sqrt(unit.rss / (length(x) - 2))
        sd_unscaled = sqrt(unit.cov[1, 1])
        sd_scaled = s * sd_unscaled
        z = quantile(Normal(), (1 + default_level) / 2)
        nsamples = 20_000
        tolerance_scaled = _interval_width_tolerance(sd_scaled, default_level, nsamples)
        tolerance_unscaled = _interval_width_tolerance(sd_unscaled, default_level, nsamples)
        # The data make the two scale hypotheses distinguishable at this nsamples.
        @test 2z * (sd_unscaled - sd_scaled) > tolerance_scaled + tolerance_unscaled

        hessian_scaled = confidence_interval(scaled, 1; method=:hessian)
        @test isapprox(_interval_width(hessian_scaled), 2z * sd_scaled; rtol=1e-10)
        bootstrap_scaled = confidence_interval(scaled, 1; method=:bootstrap, nsamples=nsamples, rng=Xoshiro(51))
        @test abs(_interval_width(bootstrap_scaled) - _interval_width(hessian_scaled)) <= tolerance_scaled

        hessian_unscaled = confidence_interval(unscaled, 1; method=:hessian)
        @test isapprox(_interval_width(hessian_unscaled), 2z * sd_unscaled; rtol=1e-10)
        bootstrap_unscaled = confidence_interval(unscaled, 1; method=:bootstrap, nsamples=nsamples, rng=Xoshiro(52))
        @test abs(_interval_width(bootstrap_unscaled) - 2z * sd_unscaled) <= tolerance_unscaled
    end

    @testset "Poisson rate bootstrap is the Poisson distribution of the count" begin
        # A constant rate on k observations has the MLE λ̂ = Σnᵢ/k. The resample
        # draws nᵢ* ~ Poisson(λ̂), so the refit MLE is Σnᵢ*/k with
        # Σnᵢ* ~ Poisson(kλ̂); for k = 1 the refit returns the drawn count itself.
        # The bounds keep the rate physical; the MLE stays interior because
        # P(Σnᵢ* = 0) = exp(-kλ̂) < 1e-21, so no refit has a reason to fail.
        #
        # Discreteness: with probability ≈ 1 - 6e-5 the empirical cdf at the
        # relevant counts lies within m = 4·sqrt(q(1 - q)/n) of the exact cdf, and
        # the order statistics bracketing the type-7 sample quantile sit within
        # 1/n of rank q. The sample quantile therefore lies between the exact
        # quantiles at q ∓ (m + 1/n); the first @test of each pair checks that
        # this bracket is within ±1 count of quantile(Poisson(kλ̂), q).
        # The refit values equal integers up to the solver tolerance; 1e-6 of a
        # count separates "refit equals the drawn count" from a wrong estimate.
        nsamples = 10_000
        for (counts, seed) in (([50], 61), ([47, 55, 50, 52], 62))
            k = length(counts)
            fitted = fit_poisson_model(_interval_constant, zeros(k), counts; p0=[40.0], bounds=([0.0], [Inf]))
            @test fitted.converged
            @test abs(k * fitted.params[1] - sum(counts)) <= 1e-6

            samples = parametric_bootstrap(fitted; nsamples=nsamples, rng=Xoshiro(seed))
            @test samples.n_requested == nsamples
            @test samples.n_failed == 0
            @test size(samples.params) == (nsamples, 1)
            totals = k .* samples.params[:, 1]
            @test maximum(abs.(totals .- round.(totals))) <= 1e-6

            exact = Poisson(k * fitted.params[1])
            interval = confidence_interval(fitted, 1; method=:bootstrap, nsamples=nsamples, rng=Xoshiro(seed))
            for (q, limit) in (((1 - default_level) / 2, interval.lower), ((1 + default_level) / 2, interval.upper))
                expected = quantile(exact, q)
                margin = 4 * sqrt(q * (1 - q) / nsamples) + 1 / nsamples
                @test expected - 1 <= quantile(exact, q - margin) && quantile(exact, q + margin) <= expected + 1
                @test abs(k * limit - expected) <= 1 + 1e-6
            end
        end
    end

    @testset "Contracts" begin
        for level in (0.0, 1.0, -0.1, 1.5, NaN)
            @test_throws DomainError confidence_interval(known, 1; method=:hessian, level=level)
        end
        @test_throws DomainError confidence_interval(known, 1; method=:bootstrap, nsamples=0)
        @test_throws DomainError parametric_bootstrap(known; nsamples=0)

        @test _interval_argument_error(() -> confidence_interval(known, 0), "outside 1:2")
        @test _interval_argument_error(() -> confidence_interval(known, 3), "outside 1:2")
        @test _interval_argument_error(() -> confidence_interval(known, :slope), "stores no parameter names")
        @test _interval_argument_error(() -> confidence_interval(known, 1; method=:nonsense), "method must be")
        for keywords in ((; nsamples=100), (; rng=Xoshiro(1)), (; npoints=31))
            @test _interval_argument_error(
                () -> confidence_interval(known, 1; method=:hessian, keywords...), "does not use the keyword")
        end
        @test _interval_argument_error(
            () -> confidence_interval(known, 1; method=:profile, nsamples=100), "does not use the keyword")

        named = fit_poisson_model(_interval_line, [0.0, 1.0], [40, 60]; p0=[10.0, 40.0],
                                  parameter_names=["slope", "offset"])
        @test _interval_argument_error(() -> confidence_interval(named, :amplitude), "unknown parameter name")
        without_covariance = fit_poisson_model(_interval_line, [0.0, 1.0], [40, 60]; p0=[10.0, 40.0],
                                               parameter_covariance=:none)
        @test _interval_argument_error(
            () -> confidence_interval(without_covariance, 1; method=:hessian), "is not finite")

        fixed = fit_model(_interval_line, x, y; p0=[1.0, 1.5], sigma_y=sigma_y, fixed_parameters=[FixedParameter(2, 1.5)])
        @test _interval_argument_error(() -> confidence_interval(fixed, 2), "is fixed")
        @test all(==(1.5), parametric_bootstrap(fixed; nsamples=50, rng=Xoshiro(3)).params[:, 2])

        eiv = fit_model(_interval_line, x, y; p0=[1.0, 0.0], sigma_y=sigma_y, sigma_x=fill(0.05, length(x)))
        @test _interval_argument_error(() -> parametric_bootstrap(eiv; nsamples=10), "x uncertainties")
        @test _interval_argument_error(
            () -> confidence_interval(eiv, 1; method=:bootstrap, nsamples=10), "x uncertainties")
        prior = fit_model(_interval_line, x, y; p0=[1.0, 0.0], sigma_y=sigma_y,
                          parameter_priors=[ParameterPrior(1, 2.0, 0.1)])
        @test _interval_argument_error(() -> parametric_bootstrap(prior; nsamples=10), "parameter priors")
        custom = fit_custom(p -> sum(abs2, (y .- (p[1] .* x .+ p[2])) ./ sigma_y); p0=[1.0, 0.0], nobs=length(x))
        @test _interval_argument_error(() -> parametric_bootstrap(custom; nsamples=10), "bootstrap_resampler")
    end

    @testset "Failed refits are counted and reported" begin
        # maxiters = 1 is inherited by every refit, so no refit can meet the
        # 1e-10 stopping tolerance. The guard fires at the first failure count
        # with 10·n_failed > nsamples, i.e. after fld(nsamples, 10) + 1 refits.
        stalled = fit_model(_interval_line, x, y; p0=[1.0, 0.0], sigma_y=sigma_y, maxiters=1)
        @test !stalled.converged
        nsamples = 20
        counts = "$(fld(nsamples, 10) + 1) of $nsamples"
        for attempt in (() -> parametric_bootstrap(stalled; nsamples=nsamples, rng=Xoshiro(71)),
                        () -> confidence_interval(stalled, 1; method=:bootstrap, nsamples=nsamples, rng=Xoshiro(71)))
            err = try
                attempt()
                nothing
            catch caught
                caught
            end
            @test err isa ErrorException
            @test err isa ErrorException && occursin(counts, err.msg)
        end
    end

    @testset "Parameter names resolve to the same interval as the index" begin
        named = fit_poisson_model(_interval_line, [0.0, 0.25, 0.5, 0.75, 1.0], [98, 105, 117, 121, 134];
                                  p0=[30.0, 90.0], parameter_names=["slope", "offset"])
        for (name, index) in ((:slope, 1), ("offset", 2))
            by_name = confidence_interval(named, name)
            by_index = confidence_interval(named, index)
            @test by_name.parameter_index == index
            @test isfinite(by_name.lower) && isfinite(by_name.upper)
            @test by_name == by_index
        end
    end
end

@testset "Unbounded Poisson refits stay feasible" begin
    # Regression: a line search stepping to a negative rate must see infinite
    # cost, not an exception — a throwing objective killed exactly the small
    # resamples (loss rate P(N <= 40 | Poisson(50)) = 8.6%), biasing the
    # percentile interval high while staying under the 10% abort.
    fitres = fit_poisson_model((x, p) -> fill(p[1], length(x)), [0.0], [50.0];
        p0=[50.000001])
    samples = parametric_bootstrap(fitres; nsamples=2000, rng=Xoshiro(5))
    @test samples.n_failed == 0
    # One observation: each refit equals its resampled count, so the bounds
    # are Poisson(50) sample quantiles. Quantile SE = sqrt(q(1-q)/n) = 0.0327
    # at 4 SE; divided by the pmf near the bounds (~0.023) that is at most
    # 1.5 counts, so atol = 2 counts.
    ci = confidence_interval(fitres, 1; method=:bootstrap, nsamples=2000, rng=Xoshiro(5))
    @test isapprox(ci.lower, quantile(Poisson(50.0), (1 - ci.level) / 2); atol=2)
    @test isapprox(ci.upper, quantile(Poisson(50.0), (1 + ci.level) / 2); atol=2)
end

@testset "Profiled-sigma likelihood resamples at RSS/n" begin
    # Regression: an unweighted fit with cost=:gaussian_likelihood reports
    # curvature from the sigma-profiled likelihood (sigma_hat^2 = RSS/n);
    # resampling at sigma = 1 produced ~5x wider bootstrap intervals.
    x = collect(1.0:10.0)
    y = 2 .* x .+ 1 .+ 0.3 .* sin.(3 .* x)
    res = fit_model((x, p) -> p[1] .* x .+ p[2], x, y; p0=[1.0, 0.0],
        cost=:gaussian_likelihood)
    h = confidence_interval(res, 1; method=:hessian)
    b = confidence_interval(res, 1; method=:bootstrap, nsamples=4000, rng=Xoshiro(9))
    width_h = h.upper - h.lower
    # Width SE (Bahadur): sd*sqrt(level(1-level)/n)/phi(z); with width ~ 2z*sd
    # the relative 4-SE tolerance is 4*sqrt(level(1-level)/n)/(2z*phi(z)) = 0.053.
    @test isapprox(b.upper - b.lower, width_h; rtol=0.053)
end
