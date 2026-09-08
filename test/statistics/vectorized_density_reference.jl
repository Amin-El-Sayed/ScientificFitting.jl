using Test
using ScientificFitting
using ScientificFitting: ForwardDiff
using ScientificFitting.QuadGK: quadgk
using Statistics

@testset "Vectorized density references" begin
    data = [0.12, 0.28, 0.51, 0.62, 0.75, 1.3, 1.8]
    scalar(x, p) = exp(-x/p[1])/p[1]
    batch(x, p) = @. exp(-x/p[1])/p[1]
    starts = (; p0=[0.5], bounds=([0.01], [5.0]))
    for derivatives in (:auto, :finite)
        @testset "$derivatives event fits and refits" begin
            result = fit_unbinned_model(batch, data; starts..., vectorized=true, derivatives)
            reference = fit_unbinned_model(scalar, data; starts..., derivatives)
            @test result.converged
            @test result.params ≈ [mean(data)] atol=2e-6
            @test result.param_covariance[1, 1] ≈ mean(data)^2/length(data) rtol=2e-5
            @test result.stats.cost_min ≈ reference.stats.cost_min atol=1e-8
            values = [0.5, mean(data), 1.0]
            scan = profile(result, 1; values)
            expected = 2*length(data) .* (log.(values ./ mean(data)) .+ mean(data) ./ values .- 1)
            @test scan.delta_cost ≈ expected atol=1e-8

            intensity(x, p) = fill(p[1], length(x))
            extended = fit_extended_unbinned_model(intensity, data, (0, 2);
                p0=[3.0], bounds=([0.01], [20.0]), vectorized=true, derivatives)
            @test extended.converged
            @test extended.params ≈ [length(data)/2] atol=2e-5
            @test extended.param_covariance[1, 1] ≈ length(data)/4 rtol=2e-5
            @test extended.stats.cost_min ≈ 2*length(data)*(1-log(length(data)/2)) atol=1e-8
        end

        @testset "$derivatives unequal-bin quadrature" begin
            edges = [0.0, 0.3, 0.9, 2.0, 5.0]
            counts = [24, 35, 28, 12]
            expected_counts(edges, p) = -100 .* diff(exp.(-edges ./ p[1]))
            reference = fit_histogram_model(expected_counts, edges, counts; starts..., derivatives)
            result = fit_histogram_density(batch, edges, counts; starts..., derivatives,
                total_count=100, vectorized=true, rtol=1e-10)
            @test result.converged
            @test result.params ≈ reference.params atol=3e-6
            @test result.param_covariance ≈ reference.param_covariance rtol=3e-5
            @test result.stats.cost_min ≈ reference.stats.cost_min atol=1e-8
        end
    end

    @testset "Adaptive batches preserve values and parameter derivatives" begin
        sizes = Int[]
        oscillatory(x, p) = (push!(sizes, length(x)); @. exp(p[1])*(1+0.2*cos(100*x)))
        integrate(p) = first(quadgk(ScientificFitting._density_integrand(oscillatory, p, true),
                                    0.0, 1.0; rtol=1e-10))
        exact(p) = exp(p[1])*(1+0.2*sin(100)/100)
        @test integrate([0.4]) ≈ exact([0.4]) rtol=1e-10
        @test length(sizes) > 1 # This integrand actually requires adaptive refinement.
        @test all(>(1), sizes)
        @test ForwardDiff.gradient(integrate, [0.4]) ≈ [exact([0.4])] rtol=1e-9
        @test ForwardDiff.hessian(integrate, [0.4]) ≈ reshape([exact([0.4])], 1, 1) rtol=1e-8
    end

    @testset "Callback count and malformed batches" begin
        many = collect(range(0.01, 4.0; length=2000))
        calls = Ref(0)
        counted(x, p) = (calls[] += 1; batch(x, p))
        result = fit_unbinned_model(counted, many; starts..., vectorized=true, derivatives=:finite)
        @test result.converged
        @test result.params ≈ [mean(many)] atol=2e-6
        @test calls[] < 1000 # One call per objective, not per event.
        for bad in ((x,p) -> 1.0, (x,p) -> ones(length(x)+1),
                    (x,p) -> fill(NaN, length(x)), (x,p) -> fill(-1.0, length(x)))
            @test_throws ArgumentError fit_unbinned_model(bad, data; starts..., vectorized=true)
        end
        @test_throws ArgumentError fit_histogram_density((x,p) -> 1.0, [0., 1., 2.], [2, 3];
            starts..., vectorized=true)
    end
end
