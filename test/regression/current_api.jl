using JuFitter
using LinearAlgebra
using StatsAPI
using Test

@testset "JuFitter" begin
    @test JuFitter.fit === StatsAPI.fit

    @testset "Linear fit with sigma_y" begin
        x = collect(range(0.0, 10.0; length=120))
        p_true = [2.0, -1.5]
        model(x, p) = @. p[1] * x + p[2]
        sigma_y = fill(0.2, length(x))
        y = model(x, p_true) .+ sigma_y .* sin.(1.3 .* x)

        result = fit_model(model, x, y; p0=[1.2, 0.0], sigma_y=sigma_y)

        @test result.backend == :lsqfit
        @test result.converged
        @test isapprox(result.params[1], p_true[1]; atol=0.08)
        @test isapprox(result.params[2], p_true[2]; atol=0.2)
        @test result.stats.ndf == length(x) - length(p_true)
        @test isfinite(result.stats.chi2)
        @test isfinite(result.stats.chi2_ndf)
        @test isfinite(result.stats.pvalue)
        @test size(result.param_covariance) == (2, 2)
        @test size(result.param_correlation) == (2, 2)
    end

    @testset "Nonlinear fit" begin
        x = collect(range(0.0, 2.0; length=100))
        p_true = [1.8, -0.9, 0.2]
        model(x, p) = @. p[1] * exp(p[2] * x) + p[3]
        y = model(x, p_true) .+ 0.03 .* cos.(2.1 .* x)

        result = fit_model(model, x, y; p0=[1.3, -0.5, 0.0], sigma_y=fill(0.03, length(x)))

        @test result.converged
        @test isapprox(result.params[1], p_true[1]; atol=0.2)
        @test isapprox(result.params[2], p_true[2]; atol=0.2)
        @test isapprox(result.params[3], p_true[3]; atol=0.2)
    end

    @testset "Diagonal and full covariance consistency" begin
        x = collect(range(-1.0, 1.0; length=80))
        p_true = [0.7, 0.9]
        model(x, p) = @. p[1] * x + p[2]
        sigma_y = fill(0.15, length(x))
        y = model(x, p_true) .+ sigma_y .* sin.(2.7 .* x)

        cov_y = Diagonal(sigma_y .^ 2)

        result_sigma = fit_model(model, x, y; p0=[0.3, 0.3], sigma_y=sigma_y)
        result_cov = fit_model(model, x, y; p0=[0.3, 0.3], cov_y=cov_y)

        @test isapprox(result_sigma.params[1], result_cov.params[1]; atol=1e-3)
        @test isapprox(result_sigma.params[2], result_cov.params[2]; atol=1e-3)
    end

    @testset "Effective-variance path with x-uncertainty" begin
        x = collect(range(0.0, 6.0; length=100))
        p_true = [1.2, 0.4]
        model(x, p) = @. p[1] * x + p[2]

        sigma_y = fill(0.1, length(x))
        sigma_x = fill(0.15, length(x))
        x_obs = x .+ sigma_x .* cos.(2.0 .* x)
        y = model(x_obs, p_true) .+ sigma_y .* sin.(2.9 .* x_obs)

        result = fit_model(model, x_obs, y; p0=[0.8, 0.0], sigma_y=sigma_y, sigma_x=sigma_x)

        @test result.backend == :optimization
        @test result.stats.cost == :gaussian_nll
        @test result.converged
        @test isfinite(result.stats.chi2)
        @test isfinite(result.stats.cost_min)
        @test result.stats.cost_min != result.stats.chi2
    end

    @testset "Vectorized x derivative avoids pointwise AD for x uncertainty" begin
        x = collect(range(0.0, 5.0; length=80))
        p_true = [1.4, -0.3]
        sigma_y = fill(0.08, length(x))
        sigma_x = fill(0.12, length(x))

        scalar_model_calls = Ref(0)
        function guarded_model(x, p)
            if length(x) == 1
                scalar_model_calls[] += 1
                error("pointwise x derivative path was used")
            end
            return @. p[1] * x + p[2]
        end

        derivative_calls = Ref(0)
        function x_derivative(x, p)
            derivative_calls[] += 1
            return fill(p[1], length(x))
        end

        y = guarded_model(x, p_true) .+ sigma_y .* sin.(1.9 .* x)
        result = fit_model(
            guarded_model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=sigma_y,
            sigma_x=sigma_x,
            x_derivative=x_derivative,
        )

        @test result.backend == :optimization
        @test result.converged
        @test derivative_calls[] > 0
        @test scalar_model_calls[] == 0
        @test isfinite(result.stats.nll_min)
    end

    @testset "Invalid x derivative fails clearly" begin
        x = collect(range(0.0, 2.0; length=20))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [1.0, 0.2])
        sigma_y = fill(0.1, length(x))
        sigma_x = fill(0.1, length(x))

        short_derivative(x, p) = fill(p[1], length(x) - 1)
        bad_derivative(x, p) = fill(NaN, length(x))

        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=sigma_y,
            sigma_x=sigma_x,
            x_derivative=short_derivative,
        )
        @test_throws ArgumentError fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=sigma_y,
            sigma_x=sigma_x,
            x_derivative=bad_derivative,
        )
    end

    @testset "Gaussian NLL for static diagonal covariance" begin
        x = collect(range(0.0, 2.0; length=60))
        model(x, p) = @. p[1] * x + p[2]
        sigma_y = fill(0.2, length(x))
        y = model(x, [1.1, 0.3]) .+ sigma_y .* sin.(x)

        result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=sigma_y, cost=:gaussian_nll)
        expected_nll = result.stats.chi2 + length(x) * log(2 * pi) + sum(log.(sigma_y .^ 2))

        @test result.backend == :optimization
        @test result.stats.cost == :gaussian_nll
        @test isapprox(result.stats.nll_min, expected_nll; atol=1e-8)
        @test isapprox(result.stats.cost_min, result.stats.nll_min; atol=1e-8)
        @test isapprox(result.stats.aic, result.stats.nll_min + 2 * length(result.params); atol=1e-8)
    end

    @testset "Analytic Jacobian is used by LsqFit backend" begin
        x = collect(range(-2.0, 2.0; length=40))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [0.7, -0.4])
        jacobian_calls = Ref(0)
        function jacobian(x, p)
            jacobian_calls[] += 1
            return hcat(x, ones(length(x)))
        end

        result = fit_model(model, x, y; p0=[0.0, 0.0], sigma_y=fill(0.1, length(x)), jacobian=jacobian)

        @test result.backend == :lsqfit
        @test jacobian_calls[] > 0
        @test isapprox(result.params[1], 0.7; atol=1e-10)
        @test isapprox(result.params[2], -0.4; atol=1e-10)
    end

    @testset "General constraints path" begin
        x = collect(range(0.0, 5.0; length=90))
        p_true = [2.2, 1.0]
        model(x, p) = @. p[1] * x + p[2]
        sigma_y = fill(0.12, length(x))
        y = model(x, p_true) .+ sigma_y .* sin.(1.7 .* x)

        constraints = (
            eq = p -> [p[2] - 1.0],
            ineq = p -> [p[1] - 3.0],
        )

        result = fit_model(
            model,
            x,
            y;
            p0=[1.0, 1.0],
            sigma_y=sigma_y,
            constraints=constraints,
        )

        @test result.backend == :optimization
        @test result.converged
        @test result.params[1] <= 3.0 + 1e-6
        @test isapprox(result.params[2], 1.0; atol=5e-3)
    end

    @testset "Parameter prior penalty" begin
        x = collect(range(0.0, 2.0; length=80))
        p_true = [3.0, 0.2]
        model(x, p) = @. p[1] * p[2] * x
        sigma_y = fill(0.05, length(x))
        y = model(x, p_true) .+ sigma_y .* sin.(2.5 .* x)

        result = fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.4],
            sigma_y=sigma_y,
            parameter_priors=(index=2, mean=0.2, sigma=0.005),
        )

        @test result.backend == :optimization
        @test result.converged
        @test isapprox(result.params[2], 0.2; atol=0.02)
        @test result.stats.ndf == length(x) + 1 - 2
    end

    @testset "Fixed parameters and asymmetric priors" begin
        x = collect(range(0.0, 4.0; length=80))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [2.0, 1.0])

        fixed_result = fit_model(
            model,
            x,
            y;
            p0=[0.0, 1.0],
            sigma_y=fill(0.1, length(x)),
            fixed_parameters=(index=2, value=1.0, sigma_minus=0.2, sigma_plus=0.4),
        )

        @test fixed_result.converged
        @test isapprox(fixed_result.params[1], 2.0; atol=1e-10)
        @test fixed_result.params[2] == 1.0
        @test fixed_result.stats.ndf == length(x) - 1
        @test isapprox(fixed_result.param_stderr[2], 0.4; atol=1e-12)
        fixed_report = fit_report(fixed_result; parameter_names=["m", "b"])
        @test fixed_report.parameters[2].fixed
        @test fixed_report.parameters[2].uncertainty_minus == 0.2
        @test fixed_report.parameters[2].uncertainty_plus == 0.4

        prior_result = fit_model(
            model,
            x,
            y;
            p0=[0.0, 2.0],
            sigma_y=fill(0.1, length(x)),
            parameter_priors=(index=2, mean=1.0, sigma_minus=0.1, sigma_plus=0.5),
        )

        @test prior_result.converged
        @test isapprox(prior_result.params[2], 1.0; atol=0.05)
    end

    @testset "Profile and contour grids" begin
        x = collect(range(-2.0, 2.0; length=60))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [1.3, -0.5]) .+ 0.03 .* sin.(x)
        result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=fill(0.05, length(x)))

        prof = profile(result, 1; npoints=9, nsigma=2)
        min_index = argmin(prof.cost_values)

        @test length(prof.values) == 9
        @test prof.delta_cost[min_index] <= 1e-6
        @test isapprox(prof.values[min_index], result.params[1]; atol=prof.values[2] - prof.values[1])

        cont = contour(result, 1, 2; npoints=5, nsigma=1)

        @test size(cont.cost_values) == (5, 5)
        @test size(cont.delta_cost) == (5, 5)
        @test minimum(cont.delta_cost) <= 1e-5
        @test cont.parameter_indices == (1, 2)

        interval = profile_interval(result, 1; npoints=9, nsigma=2)
        @test interval.parameter_index == 1
        @test isfinite(interval.uncertainty_minus)
        @test isfinite(interval.uncertainty_plus)
        profile_report = fit_report(result; parameter_names=["m", "b"], errors=:profile, profile_npoints=9, profile_nsigma=2)
        @test profile_report.parameters[1].uncertainty_minus > 0
        @test profile_report.parameters[1].uncertainty_plus > 0
    end

    @testset "Correlated parameter constraints" begin
        x = collect(range(0.0, 3.0; length=80))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [1.5, 0.5]) .+ 0.02 .* sin.(x)
        constraint = (
            indices=[1, 2],
            mean=[1.5, 0.5],
            covariance=[0.05^2 0.5 * 0.05 * 0.1; 0.5 * 0.05 * 0.1 0.1^2],
        )

        result = fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.05, length(x)),
            parameter_constraints=constraint,
        )

        @test result.converged
        @test result.backend == :optimization
        @test result.stats.ndf == length(x) + 2 - 2
        @test isapprox(result.params[1], 1.5; atol=0.03)
        @test isapprox(result.params[2], 0.5; atol=0.03)
    end

    @testset "Poisson, histogram, unbinned, and multi fits" begin
        x = collect(1.0:10.0)
        counts = [2, 4, 5, 6, 9, 11, 12, 14, 15, 19]
        poisson_model(x, p) = @. exp(p[1] + p[2] * x)
        poisson_result = fit_poisson_model(
            poisson_model,
            x,
            counts;
            p0=[0.0, 0.1],
            bounds=([-10.0, -10.0], [10.0, 10.0]),
        )

        @test poisson_result.converged
        @test poisson_result.stats.cost == :poisson_nll
        @test isfinite(poisson_result.stats.chi2)
        @test poisson_result.stats.ndf == length(counts) - 2

        edges = collect(0.0:1.0:5.0)
        hist_counts = [4, 9, 13, 20, 31]
        expected_counts(edges, p) = [p[1] * (edges[i + 1] - edges[i]) * exp(p[2] * (edges[i] + edges[i + 1]) / 2) for i in 1:(length(edges) - 1)]
        hist_result = fit_histogram_model(
            expected_counts,
            edges,
            hist_counts;
            p0=[3.0, 0.3],
            bounds=([1e-6, -5.0], [100.0, 5.0]),
        )

        @test hist_result.converged
        @test hist_result.stats.cost == :histogram_poisson_nll
        @test hist_result.stats.ndf == length(hist_counts) - 2

        data = [-1.1, -0.2, 0.1, 0.3, 0.9, 1.2]
        normal_pdf(x, p) = exp(-0.5 * ((x - p[1]) / p[2])^2) / (p[2] * sqrt(2 * pi))
        unbinned_result = fit_unbinned_model(
            normal_pdf,
            data;
            p0=[0.0, 1.0],
            bounds=([-5.0, 0.1], [5.0, 5.0]),
        )

        @test unbinned_result.converged
        @test unbinned_result.stats.cost == :unbinned_nll
        @test isapprox(unbinned_result.params[1], sum(data) / length(data); atol=1e-6)
        @test isnan(unbinned_result.stats.pvalue)

        x1 = collect(0.0:1.0:5.0)
        x2 = collect(0.0:1.0:5.0)
        model1(x, p) = @. p[1] * x + p[2]
        model2(x, p) = @. p[1] * x + p[3]
        y1 = model1(x1, [2.0, 1.0, -1.0])
        y2 = model2(x2, [2.0, 1.0, -1.0])
        multi_result = fit_multi_model(
            [model1, model2],
            [x1, x2],
            [y1, y2];
            p0=[0.0, 0.0, 0.0],
            sigma_y=[fill(0.1, length(x1)), fill(0.1, length(x2))],
            parameter_names=["slope", "offset1", "offset2"],
        )

        @test multi_result.converged
        @test multi_result.stats.cost == :multi_chi2
        @test isapprox(multi_result.params[1], 2.0; atol=1e-8)
        @test isapprox(multi_result.params[2], 1.0; atol=1e-8)
        @test isapprox(multi_result.params[3], -1.0; atol=1e-8)
        @test fit_report(multi_result).parameters[1].name == "slope"

        prof = profile(poisson_result, 1; npoints=5)
        @test length(prof.values) == 5
        @test minimum(prof.delta_cost) <= 1e-5
    end

    @testset "Custom, indexed, extended unbinned, mapped MultiFit, and diagnostics" begin
        custom = fit_custom(
            p -> sum(abs2, p .- [1.0, 2.0]);
            p0=[0.0, 0.0],
            gof=p -> sum(abs2, p .- [1.0, 2.0]),
            nobs=5,
            parameter_names=["a", "b"],
        )
        @test custom.converged
        @test custom.stats.cost == :custom
        @test isapprox(custom.params, [1.0, 2.0]; atol=1e-10)
        @test fit_report(custom).parameters[2].name == "b"

        indices = [:a, :b, :a, :c]
        y = [1.0, 2.0, 1.1, 3.0]
        indexed_model(indices, p) = [p[idx == :a ? 1 : idx == :b ? 2 : 3] for idx in indices]
        indexed = fit_indexed_model(
            indexed_model,
            indices,
            y;
            p0=[0.0, 0.0, 0.0],
            sigma_y=fill(0.1, length(y)),
        )
        @test indexed.converged
        @test indexed.stats.cost == :indexed_chi2
        @test isapprox(indexed.params, [1.05, 2.0, 3.0]; atol=1e-8)

        rate(x, p) = exp(p[1])
        extended = fit_extended_unbinned_model(rate, [0.1, 0.2, 0.8, 0.9], (0.0, 1.0); p0=[0.0])
        @test extended.converged
        @test extended.stats.cost == :extended_unbinned_nll
        @test isapprox(exp(extended.params[1]), 4.0; atol=1e-6)

        x = collect(0.0:1.0:3.0)
        local_linear(x, p) = @. p[1] * x + p[2]
        y1 = local_linear(x, [2.0, 1.0])
        y2 = local_linear(x, [2.0, -1.0])
        mapped = fit_multi_model(
            [local_linear, local_linear],
            [x, x],
            [y1, y2];
            p0=[0.0, 0.0, 0.0],
            parameter_map=[[1, 2], [1, 3]],
            sigma_y=[fill(0.1, length(x)), fill(0.1, length(x))],
            parameter_names=["m", "b1", "b2"],
        )
        @test mapped.converged
        @test isapprox(mapped.params, [2.0, 1.0, -1.0]; atol=1e-8)
        @test fit_report(mapped).parameters[3].name == "b2"

        bounded = fit_custom(
            p -> abs2(p[1] - 1.0);
            p0=[0.25],
            nobs=2,
            bounds=([0.0], [0.5]),
        )
        @test bounded.converged
        @test !isempty(bounded.diagnostics.active_bounds)
        @test any(contains("active bounds"), bounded.diagnostics.warnings)
    end

    @testset "Covariance scaling policy and multistart" begin
        x = collect(range(0.0, 5.0; length=80))
        model(x, p) = @. p[1] * exp(-p[2] * x)
        y = model(x, [2.0, 0.8])
        sigma_y = fill(0.05, length(x))

        scaled = fit_model(model, x, y; p0=[1.0, 0.2], sigma_y=sigma_y, scale_covariance=:always)
        unscaled = fit_model(model, x, y; p0=[1.0, 0.2], sigma_y=sigma_y, scale_covariance=:never)

        @test scaled.converged
        @test unscaled.converged
        @test maximum(unscaled.param_stderr) > maximum(scaled.param_stderr)

        hard = fit_model(
            model,
            x,
            y;
            p0=[0.05, 4.0],
            sigma_y=sigma_y,
            bounds=([0.01, 0.01], [5.0, 5.0]),
            initial_guesses=[[0.05, 4.0], [2.0, 0.8]],
            multistart=2,
        )

        @test hard.converged
        @test isapprox(hard.params[1], 2.0; atol=1e-2)
        @test isapprox(hard.params[2], 0.8; atol=1e-2)
    end

    @testset "No-op bounds keep the fast least-squares backend" begin
        x = collect(range(-1.0, 1.0; length=30))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [1.3, -0.2])

        result = fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            sigma_y=fill(0.05, length(x)),
            bounds=([-Inf, -Inf], [Inf, Inf]),
        )

        @test result.backend == :lsqfit
        @test isapprox(result.params, [1.3, -0.2]; atol=1e-8)
    end

    @testset "Named error components" begin
        x = collect(range(0.0, 4.0; length=50))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [1.5, 0.4]) .+ 0.02 .* sin.(x)

        result = fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            error_components=[
                (name=:stat, target=:y, mode=:absolute, values=fill(0.05, length(x))),
                (name=:scale, target=:y, mode=:relative, values=0.01),
                (name=:disabled, target=:y, mode=:absolute, values=10.0, active=false),
            ],
        )

        @test result.converged
        @test result.backend == :optimization
        @test result.stats.cost == :chi2
        @test isapprox(result.params[1], 1.5; atol=0.03)
        @test length(result.problem.error_components) == 3

        model_relative = fit_model(
            model,
            x,
            y;
            p0=[1.0, 0.0],
            error_components=[
                (name=:stat, target=:y, mode=:absolute, values=0.05),
                (name=:model_scale, target=:y, mode=:model_relative, values=0.03),
            ],
        )

        @test model_relative.converged
        @test model_relative.stats.cost == :gaussian_nll
    end

    @testset "Fit report" begin
        x = collect(range(0.0, 3.0; length=40))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [1.4, -0.2])

        result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=fill(0.1, length(x)))
        report = fit_report(result; parameter_names=["m", "b"])
        text = report_text(report)

        @test length(report.parameters) == 2
        @test report.parameters[1].name == "m"
        @test isapprox(report.parameters[1].value, result.params[1])
        @test report.statistics.ndf == result.stats.ndf
        @test report.statistics.cost == :chi2
        @test isapprox(report.statistics.cost_min, result.stats.cost_min)
        @test size(report.covariance) == (2, 2)
        @test report.diagnostics isa FitDiagnostics
        @test occursin("chi2/ndf", text)
        @test occursin("cost_min", text)
    end

    @testset "Plotting extension boundary" begin
        x = collect(range(0.0, 4.0; length=50))
        p_true = [1.5, 0.6]
        model(x, p) = @. p[1] * x + p[2]
        sigma_y = fill(0.08, length(x))
        y = model(x, p_true) .+ sigma_y .* cos.(3.3 .* x)

        result = fit_model(model, x, y; p0=[1.0, 0.2], sigma_y=sigma_y)
        err = try
            plot_fit(result; title="Linear Smoke Fit")
            nothing
        catch caught
            caught
        end

        @test err isa ArgumentError
        @test occursin("optional CairoMakie plotting extension", sprint(showerror, err))
        @test occursin("using CairoMakie", sprint(showerror, err))
    end
end
