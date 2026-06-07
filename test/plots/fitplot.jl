using JuFitter
using CairoMakie: Label, save
using Test

@testset "fitplot convenience API" begin
    x = collect(range(0.0, 3.0; length=20))
    sigma_y = fill(0.1, length(x))
    y = @. 1.5 * x - 0.2 + sigma_y * sin(2.0 * x)

    quick = fitplot(
        x,
        y;
        sigma_y=sigma_y,
        report=:none,
        parameter_names=["m", "b"],
        title="Quick linear fit",
        xlabel="x",
        xunit="s",
        ylabel="y",
        yunit="V",
        band=:none,
    )

    @test quick.result isa FitResult
    @test quick.figure !== nothing
    @test quick.result.converged
    @test length(quick.result.params) == 2

    right_panel_out = joinpath(mktempdir(), "fitplot_right_panel.png")
    right_panel = fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=right_panel_out,
        format=:png,
        report=:plot,
        show_legend=true,
        stats_position=:right,
        stats_mode=:full,
        band=:prediction,
        band_label="1-sigma prediction band",
        parameter_names=["m", "b"],
    )

    @test right_panel.figure !== nothing
    @test isfile(right_panel_out)
    @test size(right_panel.figure.scene) == (1220, 720)

    custom_out = joinpath(mktempdir(), "fitplot_customized.png")
    Label(
        right_panel.figure[2, 1:2],
        "User annotation";
        tellwidth=false,
        halign=:left,
    )
    save(custom_out, right_panel.figure)
    @test size(right_panel.figure.scene) == (1220, 720)
    @test isfile(custom_out)

    extension_out = joinpath(mktempdir(), "fitplot_extensions.png")
    ax = fit_axis(right_panel.figure)
    @test ax !== nothing
    @test add_curve!(ax, x -> 1.5 * x - 0.2; color=:black, linestyle=:dash, label="reference") !== nothing
    @test add_curve!(ax, [0.5, 2.5], [0.7, 3.4]; color=:gray40) !== nothing
    @test add_points!(ax, 1.2, 1.6; marker=:star5, color=:gold) !== nothing
    @test add_vline!(ax, 1.0; color=:gray50, linestyle=:dot) !== nothing
    @test add_hline!(ax, 2.0; color=:gray50, linestyle=:dot) !== nothing
    @test add_vband!(ax, 1.4, 1.7; color=(:gray70, 0.18)) !== nothing
    @test add_hband!(ax, 2.2, 2.6; color=(:gray70, 0.12)) !== nothing
    save(extension_out, right_panel.figure)
    @test size(right_panel.figure.scene) == (1220, 720)
    @test isfile(extension_out)
    @test_throws ArgumentError fit_axis(right_panel.figure; index=0)
    @test_throws ArgumentError add_curve!(ax, [1.0, 2.0], [1.0])
    @test_throws ArgumentError add_curve!(ax, x -> x; n=1)
    @test_throws ArgumentError add_vband!(ax, 2.0, 1.0)
    @test_throws ArgumentError add_hband!(ax, 2.0, 1.0)

    model(x, p) = @. p[1] * exp(-p[2] * x) + p[3]
    y2 = model(x, [2.0, 0.7, 0.1]) .+ sigma_y .* cos.(1.4 .* x)
    out = joinpath(mktempdir(), "fitplot.png")

    nonlinear = fitplot(
        model,
        x,
        y2;
        p0=[1.0, 0.3, 0.0],
        sigma_y=sigma_y,
        filename=out,
        format=:png,
        report=:none,
        theme=:clean,
        nsigma=2,
        band_label="2-sigma band",
        parameter_names=["A", "lambda", "C"],
    )

    @test nonlinear.result isa FitResult
    @test nonlinear.result.converged
    @test isfile(out)

    workbench_out = joinpath(mktempdir(), "fitplot_workbench.png")
    showcase_out = joinpath(mktempdir(), "fitplot_showcase_dark.png")
    publication_out = joinpath(mktempdir(), "fitplot_publication.png")

    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=workbench_out,
        format=:png,
        report=:none,
        theme=:workbench,
        parameter_names=["m", "b"],
    )
    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=showcase_out,
        format=:png,
        report=:none,
        theme=:showcase,
        appearance=:dark,
        parameter_names=["m", "b"],
    )
    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=publication_out,
        format=:png,
        report=:none,
        theme=:publication,
        parameter_names=["m", "b"],
        model_label="y = m x + b",
    )

    @test isfile(workbench_out)
    @test isfile(showcase_out)
    @test isfile(publication_out)
    @test plot_theme(:showcase; appearance=:dark) !== nothing
    @test plot_palette(:publication).fit_color == "#000000"
    @test_throws ArgumentError plot_fit(quick.result; theme=:unknown)
    @test_throws ArgumentError plot_fit(quick.result; theme=:dark, appearance=:light)

    contour_values = collect(range(-2.0, 2.0; length=17))
    contour_delta = [xv^2 + 0.5 * xv * yv + yv^2 for xv in contour_values, yv in contour_values]
    contour_result = ContourResult((1, 2), contour_values, contour_values, contour_delta, contour_delta, [2.30, 6.18])
    contour_out = joinpath(mktempdir(), "contour_levels.png")
    heatmap_out = joinpath(mktempdir(), "contour_heatmap.png")
    failed_out = joinpath(mktempdir(), "contour_failed_refit.png")

    @test plot_contour(
        contour_result;
        filename=contour_out,
        format=:png,
        theme=:showcase,
        appearance=:dark,
        local_covariance=[1.0 -0.25; -0.25 1.0],
        local_center=[0.0, 0.0],
    ) !== nothing
    @test plot_contour(contour_result; filename=heatmap_out, format=:png, show_heatmap=true) !== nothing
    profile_result = ProfileResult(1, contour_values, contour_values .^ 2, contour_values .^ 2, 1.0, 0.0)
    profile_out = joinpath(mktempdir(), "profile_legend.png")
    @test plot_profile(
        profile_result;
        filename=profile_out,
        format=:png,
        theme=:workbench,
        appearance=:dark,
        local_sigma=1.0,
        delta_max=4.0,
    ) !== nothing
    @test_throws ArgumentError plot_profile(profile_result; delta_max=0.0)
    failed_delta = copy(contour_delta)
    failed_delta[2, 3] = Inf
    failed_contour = ContourResult((1, 2), contour_values, contour_values, failed_delta, failed_delta, [2.30, 6.18])
    @test plot_contour(failed_contour; filename=failed_out, format=:png) !== nothing
    all_failed = fill(Inf, size(contour_delta))
    @test_throws ArgumentError plot_contour(
        ContourResult((1, 2), contour_values, contour_values, all_failed, all_failed, [2.30]),
    )
    @test isfile(contour_out)
    @test isfile(heatmap_out)
    @test isfile(profile_out)
    @test isfile(failed_out)

    matrix_out = joinpath(mktempdir(), "profile_matrix.png")
    matrix_fig = plot_profile_matrix(
        quick.result;
        parameters=[1, 2],
        parameter_names=["m", "b"],
        filename=matrix_out,
        format=:png,
        theme=:workbench,
    )
    @test matrix_fig !== nothing
    @test isfile(matrix_out)
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=Int[])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=[1, 1])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=[1, 3])

    residual_out = joinpath(mktempdir(), "residual_showcase_dark.png")
    diagnostics_out = joinpath(mktempdir(), "diagnostics_workbench_dark.png")
    @test plot_residuals(
        quick.result;
        filename=residual_out,
        format=:png,
        theme=:showcase,
        appearance=:dark,
    ) !== nothing
    @test plot_diagnostics(
        quick.result;
        filename=diagnostics_out,
        format=:png,
        theme=:workbench,
        appearance=:dark,
    ) !== nothing
    @test isfile(residual_out)
    @test isfile(diagnostics_out)
end
