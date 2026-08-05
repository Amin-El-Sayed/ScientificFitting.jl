using JuFitter
using CairoMakie: Axis, Figure, Label, save, scatter!, with_theme
using LaTeXStrings
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
    @test size(right_panel.figure.scene) == (1200, 720)

    custom_out = joinpath(mktempdir(), "fitplot_customized.png")
    Label(
        right_panel.figure[2, 1:2],
        "User annotation";
        tellwidth=false,
        halign=:left,
    )
    save(custom_out, right_panel.figure)
    @test size(right_panel.figure.scene) == (1200, 720)
    @test isfile(custom_out)

    extension_out = joinpath(mktempdir(), "fitplot_extensions.png")
    ax = fit_axis(right_panel.figure)
    @test ax !== nothing
    @test add_curve!(ax, x -> 1.5 * x - 0.2; color=:black, linestyle=:dash, label="reference") !== nothing
    @test add_curve!(ax, [0.5, 2.5], [0.7, 3.4]; color=:gray40) !== nothing
    @test add_points!(ax, 1.2, 1.6; marker=:star5, color=:gray25) !== nothing
    @test add_vline!(ax, 1.0; color=:gray50, linestyle=:dot) !== nothing
    @test add_hline!(ax, 2.0; color=:gray50, linestyle=:dot) !== nothing
    @test add_vband!(ax, 1.4, 1.7; color=(:gray70, 0.18)) !== nothing
    @test add_hband!(ax, 2.2, 2.6; color=(:gray70, 0.12)) !== nothing
    save(extension_out, right_panel.figure)
    @test size(right_panel.figure.scene) == (1200, 720)
    @test isfile(extension_out)
    @test_throws ArgumentError fit_axis(right_panel.figure; index=0)
    @test_throws ArgumentError add_curve!(ax, [1.0, 2.0], [1.0])
    @test_throws ArgumentError add_curve!(ax, [1.0], [2.0])
    @test_throws ArgumentError add_curve!(ax, [1.0, NaN], [2.0, 3.0])
    @test_throws ArgumentError add_curve!(ax, x -> x == 0 ? Inf : x; xspan=(0.0, 1.0))
    @test_throws ArgumentError add_curve!(ax, x -> x; n=1)
    @test_throws ArgumentError add_curve!(ax, x -> x; xspan=(0.0, Inf))
    @test_throws ArgumentError add_points!(ax, [1.0, Inf], [2.0, 3.0])
    @test_throws ArgumentError add_vline!(ax, NaN)
    @test_throws ArgumentError add_hline!(ax, Inf)
    @test_throws ArgumentError add_vband!(ax, 2.0, 1.0)
    @test_throws ArgumentError add_vband!(ax, NaN, 1.0)
    @test_throws ArgumentError add_hband!(ax, 2.0, 1.0)
    @test_throws ArgumentError add_hband!(ax, 0.0, Inf)

    # Axis-relative spans must not turn Makie's provisional 0:10 limits into
    # data limits when annotations are added before the first render.
    span_figure = Figure()
    span_axis = Axis(span_figure[1, 1])
    scatter!(span_axis, [0.0, 1.0], [2.0, 3.0])
    add_vband!(span_axis, 0.25, 0.5; color=(:gray70, 0.18))
    add_hband!(span_axis, 2.3, 2.7; color=(:gray70, 0.12))
    save(joinpath(mktempdir(), "axis_relative_spans.png"), span_figure)
    span_rect = span_axis.finallimits[]
    span_xlimits = (span_rect.origin[1], span_rect.origin[1] + span_rect.widths[1])
    span_ylimits = (span_rect.origin[2], span_rect.origin[2] + span_rect.widths[2])
    @test span_xlimits[1] < 0.0 < 1.0 < span_xlimits[2]
    @test span_ylimits[1] < 2.0 < 3.0 < span_ylimits[2]
    @test span_xlimits[2] - span_xlimits[1] < 1.5
    @test span_ylimits[2] - span_ylimits[1] < 1.5

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
        theme=:lab,
        nsigma=2,
        band_label="2-sigma band",
        parameter_names=["A", "lambda", "C"],
    )

    @test nonlinear.result isa FitResult
    @test nonlinear.result.converged
    @test isfile(out)

    function model!(out, x, p)
        @. out = p[1] * exp(-p[2] * x) + p[3]
        return nothing
    end
    inplace_plot = fitplot(
        model!,
        x,
        y2;
        p0=[1.0, 0.3, 0.0],
        sigma_y=sigma_y,
        inplace=true,
        report=:none,
        band=:none,
    )
    @test inplace_plot.result.backend == :lsqfit
    @test inplace_plot.result.converged
    @test inplace_plot.figure !== nothing

    rho = 0.45
    whitening_sigma = 0.10
    innovation_sigma = whitening_sigma * sqrt(1 - rho^2)
    function whiten!(out, residual)
        out[1] = residual[1] / whitening_sigma
        @inbounds for i in 2:length(residual)
            out[i] = (residual[i] - rho * residual[i - 1]) / innovation_sigma
        end
        return nothing
    end
    whitening = WhiteningOperator(
        whiten!;
        logdet_covariance=2 * length(x) * log(whitening_sigma) +
                              (length(x) - 1) * log1p(-rho^2),
        marginal_sigma=whitening_sigma,
    )
    linear_model(x, p) = @. p[1] * x + p[2]
    whitening_plot = fitplot(
        linear_model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening,
        report=:none,
    )
    @test whitening_plot.result.backend == :lsqfit
    @test whitening_plot.result.problem.whitening === whitening
    @test whitening_plot.figure !== nothing

    whitening_without_marginals = WhiteningOperator(
        whiten!;
        logdet_covariance=whitening.logdet_covariance,
    )
    result_without_marginals = fit_model(
        linear_model,
        x,
        y;
        p0=[1.0, 0.0],
        whitening=whitening_without_marginals,
    )
    @test_throws ArgumentError plot_fit(
        result_without_marginals;
        show_stats=false,
        show_legend=false,
        band=:prediction,
    )
    @test plot_fit(
        result_without_marginals;
        show_stats=false,
        show_legend=false,
        band=:confidence,
    ) !== nothing

    screen_out = joinpath(mktempdir(), "fitplot_screen.png")
    lab_out = joinpath(mktempdir(), "fitplot_lab.png")
    legacy_out = joinpath(mktempdir(), "fitplot_legacy_dark.png")
    article_out = joinpath(mktempdir(), "fitplot_article.png")

    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=lab_out,
        format=:png,
        report=:none,
        theme=:lab,
        parameter_names=["m", "b"],
    )
    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=screen_out,
        format=:png,
        report=:none,
        theme=:screen,
        parameter_names=["m", "b"],
    )
    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=legacy_out,
        format=:png,
        report=:none,
        theme=:modern,
        appearance=:dark,
        parameter_names=["m", "b"],
    )
    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=article_out,
        format=:png,
        report=:none,
        theme=:article,
        parameter_names=["m", "b"],
        model_label="y = m x + b",
    )

    @test isfile(screen_out)
    @test isfile(lab_out)
    @test isfile(legacy_out)
    @test isfile(article_out)
    @test plot_theme(:modern; appearance=:dark) !== nothing
    screen_style = plot_palette(:screen)
    lab_style = plot_palette(:lab)
    modern_style = plot_palette(:modern)
    workbench_style = plot_palette(:workbench)
    article_style = plot_palette(:article)
    @test lab_style != screen_style
    @test modern_style == screen_style
    @test workbench_style == lab_style
    @test lab_style.data_marker == :cross
    @test screen_style.data_marker == :circle
    @test lab_style.xgridvisible && lab_style.ygridvisible
    @test lab_style.topspinevisible && lab_style.rightspinevisible
    @test screen_style.xgridvisible && screen_style.ygridvisible
    @test !screen_style.topspinevisible && !screen_style.rightspinevisible
    @test !article_style.xgridvisible && !article_style.ygridvisible
    @test article_style.topspinevisible && article_style.rightspinevisible
    @test screen_style.titlealign == :left
    @test article_style.tickalign == 1.0
    @test screen_style.ticklabelsize >= 22
    @test screen_style.xlabelsize >= 28
    @test screen_style.ylabelsize >= 28
    @test screen_style.stats_fontsize >= 24
    @test screen_style.legend_labelsize >= 22
    @test lab_style.ticklabelsize >= 20
    @test lab_style.xlabelsize >= 25
    @test lab_style.ylabelsize >= 25
    @test lab_style.stats_fontsize >= 21
    @test lab_style.legend_labelsize >= 20
    @test article_style.ticklabelsize >= 26
    @test article_style.xlabelsize >= 32
    @test article_style.ylabelsize >= 32
    @test article_style.stats_fontsize >= 30
    @test article_style.legend_labelsize >= 26
    @test minimum(style.spinewidth for style in (lab_style, screen_style, article_style)) >= 1.5
    @test maximum(style.error_whiskerwidth for style in (lab_style, screen_style, article_style)) <= 6
    @test all(length(unique(style.series_colors)) == length(style.series_colors) for
        style in (lab_style, screen_style, article_style))
    @test lab_style.fit_color == "#005a8d"
    @test screen_style.fit_color == :dodgerblue
    @test screen_style.band_color == screen_style.fit_color
    @test screen_style.secondary_color == "#b83280"
    @test article_style.fit_color == "#0072b2"
    @test article_style.secondary_color == "#cc79a7"
    @test article_style.data_color == article_style.background_color
    @test article_style.data_strokecolor == article_style.axis_color
    @test article_style.data_strokewidth > screen_style.data_strokewidth

    style_signature(style) = (
        style.data_marker,
        style.fit_color,
        style.xgridvisible,
        style.ygridvisible,
        style.topspinevisible,
        style.rightspinevisible,
        style.titlealign,
        style.tickalign,
    )
    @test length(unique(style_signature.((lab_style, screen_style, article_style)))) == 3

    screen_axis = fit_axis(plot_fit(quick.result; theme=:screen, show_stats=false))
    lab_axis = fit_axis(plot_fit(quick.result; theme=:lab, show_stats=false))
    modern_axis = fit_axis(plot_fit(quick.result; theme=:modern, show_stats=false))
    article_axis = fit_axis(plot_fit(quick.result; theme=:article, show_stats=false))
    @test screen_axis.xgridvisible[] && screen_axis.ygridvisible[]
    @test lab_axis.xgridvisible[] && lab_axis.ygridvisible[]
    @test lab_axis.topspinevisible[] && lab_axis.rightspinevisible[]
    @test modern_axis.xgridvisible[] && modern_axis.ygridvisible[]
    @test article_axis.topspinevisible[] && article_axis.rightspinevisible[]

    panel_figure = with_theme(plot_theme(:screen)) do
        Figure(size=(720, 420))
    end
    @test plot_info_panel!(
        panel_figure[1, 1];
        theme=:screen,
        model_label="y = m x + b",
        parameter_lines=["m = 1.0 +/- 0.1"],
        statistic_lines=["chi2/ndf = 1.0"],
    ) !== nothing
    latex_figure = plot_fit(
        quick.result;
        theme=:article,
        latex_labels=true,
        latex_stats=true,
        model_label=L"U_0(\nu)=h\nu/e-\Phi/e",
        xlabel=L"\nu",
        xunit=L"\mathrm{THz}",
    )
    @test latex_figure !== nothing
    plotting_extension = Base.get_extension(JuFitter, :JuFitterCairoMakieExt)
    @test plotting_extension !== nothing
    article_dark_diagnostics = plotting_extension._diagnostic_colors(:article, :dark)
    @test article_dark_diagnostics.levels[1] == plot_palette(:article; appearance=:dark).fit_color
    @test all(color -> first(color) != :black, article_dark_diagnostics.regions)
    @test plotting_extension._panel_status_color(:review, :light) == "#4b5560"
    @test plotting_extension._panel_status_color(:review, :dark) == "#d3dae0"
    @test_throws ArgumentError plot_fit(quick.result; theme=:unknown)
    @test_throws ArgumentError plot_fit(quick.result; theme=:dark, appearance=:light)
    @test_throws ArgumentError plot_fit(quick.result; fit_range=:unknown)

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
        theme=:modern,
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
        theme=:lab,
        appearance=:dark,
        local_sigma=1.0,
        delta_max=4.0,
        threshold_kwargs=(linewidth=2.5,),
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
        theme=:lab,
        panel_status_mode=:all,
    )
    @test matrix_fig !== nothing
    @test isfile(matrix_out)
    precomputed_matrix = profile_matrix(
        quick.result;
        parameters=[1, 2],
        parameter_names=["m", "b"],
        npoints_profile=5,
        npoints_contour=3,
        nsigma=1,
        contour_levels=[2.30],
    )
    precomputed_matrix_out = joinpath(mktempdir(), "precomputed_profile_matrix.png")
    @test plot_profile_matrix(
        precomputed_matrix;
        parameter_names=[L"m", L"b"],
        filename=precomputed_matrix_out,
        format=:png,
        theme=:modern,
    ) !== nothing
    @test isfile(precomputed_matrix_out)
    @test_throws ArgumentError plot_profile_matrix(precomputed_matrix; parameter_names=["m"])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=Int[])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=[1, 1])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=[1, 3])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=[1, 2], panel_status_mode=:bad)

    residual_out = joinpath(mktempdir(), "residual_modern_dark.png")
    diagnostics_out = joinpath(mktempdir(), "diagnostics_lab_dark.png")
    @test plot_residuals(
        quick.result;
        filename=residual_out,
        format=:png,
        theme=:modern,
        appearance=:dark,
    ) !== nothing
    @test plot_diagnostics(
        quick.result;
        filename=diagnostics_out,
        format=:png,
        theme=:lab,
        appearance=:dark,
    ) !== nothing
    @test plot_residuals(
        quick.result;
        theme=:lab,
        marker=:rect,
        markersize=12,
        error_whiskerwidth=4,
        scatter_kwargs=(color=:red,),
        errorbars_kwargs=(whiskerwidth=2,),
    ) !== nothing
    @test plot_diagnostics(
        quick.result;
        theme=:modern,
        marker=:diamond,
        markersize=10,
        error_whiskerwidth=4,
        scatter_kwargs=(color=:red,),
        errorbars_kwargs=(whiskerwidth=2,),
        reference_line_kwargs=(linewidth=2.0,),
    ) !== nothing
    @test isfile(residual_out)
    @test isfile(diagnostics_out)

    zero_model(x, p) = p[1] .* x
    zero_result = fit_model(zero_model, x, zero.(x); p0=[0.0], sigma_y=sigma_y)
    @test_throws ArgumentError plot_residuals(zero_result; kind=:ratio)
    @test_throws ArgumentError plot_diagnostics(zero_result)
end
