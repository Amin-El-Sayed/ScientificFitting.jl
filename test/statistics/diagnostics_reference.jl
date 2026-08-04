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

    @testset "Negative local curvature stops reporting" begin
        stationary_maximum = fit_custom(
            p -> -abs2(p[1]);
            p0=[0.0],
            nobs=2,
            bounds=([-1.0], [1.0]),
        )

        @test stationary_maximum.converged
        @test stationary_maximum.param_covariance[1, 1] < 0
        @test isnan(stationary_maximum.param_stderr[1])
        @test isnan(stationary_maximum.param_correlation[1, 1])
        @test any(
            f -> f.code == :nonpositive_parameter_covariance && f.severity == :critical,
            diagnose(stationary_maximum).findings,
        )
        @test any(contains("symmetric parameter errors must not be reported"), stationary_maximum.diagnostics.warnings)
        @test diagnostic_dashboard(stationary_maximum).status == :stop
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

        block_x = collect(1.0:30.0)
        block_y = vcat(fill(1.0, 10), fill(-1.0, 10), fill(1.0, 10))
        constant_model(x, p) = fill(p[1], length(x))
        block_result = fit_model(constant_model, block_x, block_y; p0=[0.0], sigma_y=fill(0.2, length(block_x)))
        block_report = diagnose(block_result)

        run_finding = only(filter(f -> f.code == :long_same_sign_pull_run, block_report.findings))
        @test run_finding.severity == :warning
        @test contains(run_finding.evidence, "point 1 to 10")
        @test contains(run_finding.evidence, "x = 1.0 to 10.0")
        @test contains(run_finding.recommendation, "acquisition interval")

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
        bounded_report = diagnose(bounded)
        @test any(f -> f.code == :active_bounds, bounded_report.findings)
        @test any(f -> f.code == :local_covariance_requires_profile_check, bounded_report.findings)
    end

    @testset "Likelihood results use goodness-of-fit diagnostics" begin
        x = collect(0.0:1.0:5.0)
        local_linear(x, p) = @. p[1] * x + p[2]
        sigma = fill(0.1, length(x))
        y_a = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0]
        y_b = [0.0, 2.4, 4.8, 7.2, 9.6, 12.0]

        result = fit_multi_model(
            [local_linear, local_linear],
            [x, x],
            [y_a, y_b];
            p0=[2.0, 1.0, 0.0],
            sigma_y=[sigma, sigma],
            parameter_map=[[1, 2], [1, 3]],
        )

        @test result.converged
        @test result.stats.pvalue < 0.01
        @test any(f -> f.code in (:small_pvalue, :tiny_pvalue), diagnose(result).findings)
        @test diagnostic_dashboard(result).status != :ok
    end

    @testset "Diagnostic dashboard prioritizes lab next actions" begin
        clean_dashboard = diagnostic_dashboard(DiagnosticReport(DiagnosticFinding[], "Synthetic clean report."))

        @test clean_dashboard isa DiagnosticDashboard
        @test clean_dashboard.status == :ok
        @test clean_dashboard.severity_counts[:critical] == 0
        @test clean_dashboard.severity_counts[:warning] == 0
        @test isempty(clean_dashboard.next_actions)
        @test contains(diagnostic_dashboard_text(clean_dashboard), "status = ok")

        clean_model(x, p) = @. p[1] * x + p[2]
        bad_x = collect(range(-3.0, 3.0; length=41))
        bad_y = @. 0.4 * bad_x^2 + 0.8 * bad_x - 1.0
        bad = fit_model(clean_model, bad_x, bad_y; p0=[0.0, 0.0], sigma_y=fill(0.05, length(bad_x)))
        bad_dashboard = diagnostic_dashboard(bad; max_actions=3)

        @test bad_dashboard.status == :stop
        @test bad_dashboard.severity_counts[:critical] > 0
        @test length(bad_dashboard.next_actions) == 3
        @test allunique(lowercase.(strip.(bad_dashboard.next_actions)))
        @test contains(diagnostic_dashboard_text(bad_dashboard), "Next actions:")

        warning_report = DiagnosticReport(
            DiagnosticFinding[
                DiagnosticFinding(:warning, :a, "A", "evidence", "Inspect residuals."),
                DiagnosticFinding(:warning, :b, "B", "evidence", "Inspect residuals."),
                DiagnosticFinding(:info, :c, "C", "evidence", "Record context."),
            ],
            "Synthetic warning report.",
        )
        warning_dashboard = diagnostic_dashboard(warning_report)

        @test warning_dashboard.status == :review
        @test warning_dashboard.severity_counts[:warning] == 2
        @test warning_dashboard.next_actions == ["Inspect residuals.", "Record context."]
    end

    @testset "Diagnosis flags strong parameter correlations" begin
        x = collect(range(99.0, 101.0; length=50))
        model(x, p) = @. p[1] * x + p[2]
        y = model(x, [0.8, -20.0])

        result = fit_model(model, x, y; p0=[1.0, 0.0], sigma_y=fill(0.05, length(x)))
        report = diagnose(result)

        @test any(f -> f.code == :strong_parameter_correlation, report.findings)
        @test any(f -> f.code == :local_covariance_requires_profile_check, report.findings)
        finding = only(filter(f -> f.code == :strong_parameter_correlation, report.findings))
        @test finding.severity in (:warning, :critical)
        @test contains(finding.recommendation, "Re-center")
        local_covariance_finding = only(filter(f -> f.code == :local_covariance_requires_profile_check, report.findings))
        @test contains(local_covariance_finding.evidence, "strong parameter correlation")
        @test contains(local_covariance_finding.recommendation, "profile or contour")
    end
end
