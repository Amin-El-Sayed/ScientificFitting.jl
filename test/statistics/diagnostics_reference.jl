using JuFitter
using Test

@testset "Statistical diagnostic warnings" begin
    @testset "Unbinned likelihood warns about unavailable p-values" begin
        data = Float64[-0.8, -0.1, 0.2, 0.7]
        sigma = 0.5
        normal_location_pdf(x, p) = exp(-0.5 * abs2((x - p[1]) / sigma)) / (sigma * sqrt(2pi))

        result = fit_unbinned_model(normal_location_pdf, data; p0=[0.0])

        @test isnan(result.stats.chi2)
        @test isnan(result.stats.pvalue)
        @test any(contains("goodness-of-fit statistic is unavailable"), result.diagnostics.warnings)
        report = diagnose(result)
        @test report isa DiagnosticReport
        @test any(f -> f.code == :gof_unavailable, report.findings)
    end

    @testset "Infinite covariance and Hessian condition numbers are warned" begin
        flat = fit_custom(_ -> 0.0; p0=[1.0], gof=_ -> 0.0, nobs=2)

        @test isinf(flat.diagnostics.covariance_condition)
        @test isinf(flat.diagnostics.hessian_condition)
        @test any(contains("parameter covariance is ill-conditioned"), flat.diagnostics.warnings)
        @test any(contains("cost Hessian is ill-conditioned"), flat.diagnostics.warnings)
        report = diagnose(flat)
        @test any(f -> f.code == :ill_conditioned_covariance, report.findings)
        @test any(f -> f.code == :ill_conditioned_hessian, report.findings)
    end

    @testset "Non-positive degrees of freedom warns about reduced statistics" begin
        x = [0.0, 1.0]
        y = [1.0, 3.0]
        model(x, p) = @. p[1] * x + p[2]

        result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=fill(0.1, length(x)))

        @test result.stats.ndf == 0
        @test isnan(result.stats.chi2_ndf)
        @test isnan(result.stats.pvalue)
        @test any(contains("non-positive degrees of freedom"), result.diagnostics.warnings)
        @test any(f -> f.code == :nonpositive_ndf, diagnose(result).findings)
    end

    @testset "Actionable diagnosis flags model and uncertainty pathologies" begin
        x = collect(range(-3.0, 3.0; length=41))
        quadratic = @. 0.4 * x^2 + 0.8 * x - 1.0
        linear_model(x, p) = @. p[1] * x + p[2]

        bad = fit_model(linear_model, x, quadratic; p0=[0.0, 0.0], sigma_y=fill(0.05, length(x)))
        bad_report = diagnose(bad)

        @test any(f -> f.code == :very_large_reduced_chi2, bad_report.findings)
        @test any(f -> f.code == :tiny_pvalue, bad_report.findings)
        @test any(f -> f.code == :structured_residual_signs, bad_report.findings)
        pull_finding = only(filter(f -> f.code == :extreme_pull, bad_report.findings))
        @test contains(pull_finding.evidence, "point")
        @test contains(pull_finding.evidence, "x =")
        @test contains(diagnose_text(bad_report), "action:")

        too_good = fit_model(linear_model, x, linear_model(x, [1.2, -0.3]); p0=[1.0, 0.0], sigma_y=fill(10.0, length(x)))
        @test any(f -> f.code == :very_small_reduced_chi2, diagnose(too_good).findings)

        bounded = fit_model(
            linear_model,
            x,
            linear_model(x, [2.0, 1.0]);
            p0=[0.0, 0.0],
            sigma_y=fill(0.1, length(x)),
            bounds=([-Inf, -Inf], [1.0, Inf]),
        )
        @test any(f -> f.code == :active_bounds, diagnose(bounded).findings)
    end

    @testset "Diagnosis flags strong parameter correlations" begin
        x = collect(range(99.0, 101.0; length=50))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [0.8, -20.0])

        result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=fill(0.05, length(x)))
        report = diagnose(result)

        @test any(f -> f.code == :strong_parameter_correlation, report.findings)
        finding = only(filter(f -> f.code == :strong_parameter_correlation, report.findings))
        @test finding.severity in (:warning, :critical)
        @test contains(finding.recommendation, "Re-center")
    end
end
