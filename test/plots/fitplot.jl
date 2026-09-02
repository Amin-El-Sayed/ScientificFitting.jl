using ScientificFitting
using CairoMakie: Axis, Figure, Label, errorbars!, save, scatter!, with_theme
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
        show_panel=false,
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
        show_panel=true,
        show_legend=true,
        stats_position=:right,
        stats_mode=:full,
        band=:prediction,
        band_label="1-sigma prediction band",
        parameter_names=["m", "b"],
    )

    @test right_panel.figure !== nothing
    @test isfile(right_panel_out)
    @test size(right_panel.figure.scene) == (1040, 640)

    custom_out = joinpath(mktempdir(), "fitplot_customized.png")
    Label(
        right_panel.figure[2, 1:2],
        "User annotation";
        tellwidth=false,
        halign=:left,
    )
    save(custom_out, right_panel.figure)
    @test size(right_panel.figure.scene) == (1040, 640)
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
    @test size(right_panel.figure.scene) == (1040, 640)
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
        show_panel=false,
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
        show_panel=false,
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
        show_panel=false,
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
        show_panel=false,
        show_legend=false,
        band=:prediction,
    )
    @test plot_fit(
        result_without_marginals;
        show_panel=false,
        show_legend=false,
        band=:confidence,
    ) !== nothing

    analysis_out = joinpath(mktempdir(), "fitplot_analysis.png")
    lab_out = joinpath(mktempdir(), "fitplot_lab.png")
    legacy_out = joinpath(mktempdir(), "fitplot_legacy_dark.png")
    article_out = joinpath(mktempdir(), "fitplot_article.png")

    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=lab_out,
        format=:png,
        show_panel=false,
        theme=:lab,
        parameter_names=["m", "b"],
    )
    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=analysis_out,
        format=:png,
        show_panel=false,
        theme=:analysis,
        parameter_names=["m", "b"],
    )
    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=legacy_out,
        format=:png,
        show_panel=false,
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
        show_panel=false,
        theme=:article,
        parameter_names=["m", "b"],
        model_label="y = m x + b",
    )

    @test isfile(analysis_out)
    @test isfile(lab_out)
    @test isfile(legacy_out)
    @test isfile(article_out)
    @test plot_theme(:modern; appearance=:dark) !== nothing
    sans_style = plot_palette(:sans)
    tex_style = plot_palette(:tex)
    analysis_style = plot_palette(:analysis)
    screen_style = plot_palette(:screen)
    lab_style = plot_palette(:lab)
    modern_style = plot_palette(:modern)
    workbench_style = plot_palette(:workbench)
    presentation_style = plot_palette(:presentation)
    showcase_style = plot_palette(:showcase)
    article_style = plot_palette(:article)
    @test analysis_style == sans_style
    @test screen_style == analysis_style
    @test lab_style == analysis_style
    @test workbench_style == analysis_style
    @test modern_style == sans_style
    @test presentation_style == sans_style
    @test showcase_style == presentation_style
    @test article_style == tex_style
    @test sans_style.name == :sans
    @test tex_style.name == :tex
    @test tex_style.diagnostic_scale > sans_style.diagnostic_scale
    @test sans_style.data_marker == :circle
    @test sans_style.xgridvisible && sans_style.ygridvisible
    @test !sans_style.topspinevisible && !sans_style.rightspinevisible
    @test !tex_style.xgridvisible && !tex_style.ygridvisible
    @test tex_style.topspinevisible && tex_style.rightspinevisible
    @test sans_style.titlealign == :left
    @test tex_style.tickalign == 1.0
    @test sans_style.ticklabelsize >= 22
    @test sans_style.xlabelsize >= 28
    @test sans_style.ylabelsize >= 28
    @test sans_style.stats_fontsize >= 24
    @test sans_style.legend_labelsize >= 22
    @test sans_style.subplot_titlesize >= sans_style.ticklabelsize
    @test tex_style.ticklabelsize >= 26
    @test tex_style.xlabelsize >= 32
    @test tex_style.ylabelsize >= 32
    @test tex_style.stats_fontsize >= 30
    @test tex_style.legend_labelsize >= 26
    @test tex_style.subplot_titlesize >= tex_style.ticklabelsize
    @test minimum(style.spinewidth for style in
        (sans_style, tex_style)) >= 1.5
    @test maximum(style.error_whiskerwidth for style in
        (sans_style, tex_style)) <= 6
    @test all(style.data_markersize <= 9 for style in (sans_style, tex_style))
    @test sans_style.data_markersize == tex_style.data_markersize
    @test sans_style.error_linewidth == tex_style.error_linewidth
    @test sans_style.error_linewidth <= 1.1
    @test sans_style.xerr_color == sans_style.axis_color
    @test sans_style.yerr_color == sans_style.axis_color
    @test tex_style.xerr_color == tex_style.axis_color
    @test tex_style.yerr_color == tex_style.axis_color
    @test all(length(unique(style.series_colors)) == length(style.series_colors) for
        style in (sans_style, tex_style))
    @test sans_style.fit_color == :dodgerblue
    @test sans_style.band_color == sans_style.fit_color
    @test sans_style.secondary_color == :red
    @test tex_style.fit_color == "#0072b2"
    @test tex_style.secondary_color == "#d55e00"
    @test tex_style.data_color == tex_style.background_color
    @test tex_style.data_strokecolor == tex_style.axis_color
    @test tex_style.data_strokewidth > sans_style.data_strokewidth

    # Direct Makie composition must inherit the same error-bar contract as
    # `plot_fit`; disabling Makie's color cycle keeps errors tied to the axes.
    themed_errorbars = with_theme(plot_theme(:sans)) do
        figure = Figure()
        axis = Axis(figure[1, 1])
        errorbars!(axis, [1.0], [2.0], [0.2])
    end
    @test themed_errorbars.color[] == sans_style.axis_color
    @test themed_errorbars.linewidth[] == sans_style.error_linewidth
    @test themed_errorbars.whiskerwidth[] == sans_style.error_whiskerwidth

    style_signature(style) = (
        style.data_marker,
        style.fit_color,
        style.xgridvisible,
        style.ygridvisible,
        style.topspinevisible,
        style.rightspinevisible,
        style.titlealign,
        style.tickalign,
        style.ticklabelsize,
        style.data_markersize,
    )
    @test length(unique(style_signature.((sans_style, tex_style)))) == 2

    sans_axis = fit_axis(plot_fit(quick.result; theme=:sans, show_panel=false))
    lab_axis = fit_axis(plot_fit(quick.result; theme=:lab, show_panel=false))
    modern_axis = fit_axis(plot_fit(quick.result; theme=:modern, show_panel=false))
    tex_axis = fit_axis(plot_fit(quick.result; theme=:tex, show_panel=false))
    @test sans_axis.xgridvisible[] && sans_axis.ygridvisible[]
    @test lab_axis.xgridvisible[] && lab_axis.ygridvisible[]
    @test !lab_axis.topspinevisible[] && !lab_axis.rightspinevisible[]
    @test modern_axis.xgridvisible[] && modern_axis.ygridvisible[]
    @test tex_axis.topspinevisible[] && tex_axis.rightspinevisible[]

    # Style and information density are orthogonal: both styles support both
    # panel states, and both default to a self-contained figure.
    sans_default = plot_fit(quick.result; theme=:sans)
    tex_default = plot_fit(quick.result; theme=:tex)
    sans_without_panel = plot_fit(quick.result; theme=:sans, show_panel=false)
    tex_without_panel = plot_fit(quick.result; theme=:tex, show_panel=false)
    @test size(sans_default.scene) == sans_style.figure_size_with_panel
    @test size(tex_default.scene) == tex_style.figure_size_with_panel
    @test size(sans_without_panel.scene) == sans_style.figure_size_without_panel
    @test size(tex_without_panel.scene) == tex_style.figure_size_without_panel

    panel_figure = with_theme(plot_theme(:sans)) do
        Figure(size=(720, 420))
    end
    panel = plot_info_panel!(
        panel_figure[1, 1];
        theme=:sans,
        model_label="y = m x + b",
        parameter_lines=["m = 1.0 +/- 0.1"],
        statistic_lines=["chi2/ndf = 1.0"],
    )
    @test panel !== nothing
    @test panel.width[] isa Auto
    bounded_panel = plot_info_panel!(
        panel_figure[1, 2];
        width=320,
        parameter_lines=["a deliberately long parameter description that must wrap"],
    )
    @test bounded_panel.width[] == 320
    plotting_extension = Base.get_extension(ScientificFitting, :ScientificFittingCairoMakieExt)
    @test plotting_extension !== nothing
    @test occursin('\n', plotting_extension._wrap_panel_text(
        "a deliberately long parameter description that must wrap",
        100,
        20,
    ))
    @test_throws ArgumentError plot_info_panel!(panel_figure[1, 1]; width=0)
    latex_figure = plot_fit(
        quick.result;
        theme=:tex,
        latex_labels=true,
        latex_stats=true,
        model_label=L"U_0(\nu)=h\nu/e-\Phi/e",
        xlabel=L"\nu",
        xunit=L"\mathrm{THz}",
    )
    @test latex_figure !== nothing
    @test fit_axis(quick.figure).xlabel[] == "x / s"
    @test ScientificFitting._strip_math_delims(String(fit_axis(latex_figure).xlabel[])) ==
          "\\nu\\,/\\,\\mathrm{THz}"

    mktemp() do _, report_stream
        redirect_stdout(report_stream) do
            fitplot(quick.result; show_panel=false, print_report=true, band=:none)
        end
        flush(report_stream)
        seekstart(report_stream)
        @test occursin("Fit report", read(report_stream, String))
    end

    @test isnothing(plotting_extension._validate_panel_status_mode(:issues))
    @test isnothing(plotting_extension._validate_panel_status_mode(:all))
    tex_dark_diagnostics = plotting_extension._diagnostic_colors(:tex, :dark)
    @test tex_dark_diagnostics.levels[1] == plot_palette(:tex; appearance=:dark).fit_color
    @test all(color -> first(color) != :black, tex_dark_diagnostics.regions)
    @test plotting_extension._panel_status_color(:review, :light) == "#4b5560"
    @test plotting_extension._panel_status_color(:review, :dark) == "#d3dae0"
    @test_throws ArgumentError plot_fit(quick.result; theme=:unknown)
    @test_throws ArgumentError plot_fit(quick.result; theme=:custom)
    @test_throws ArgumentError plot_fit(quick.result; theme=:dark)
    @test_throws ArgumentError plot_fit(quick.result; fit_range=:unknown)
    @test_throws ArgumentError plot_fit(quick.result; nsigma=0)
    @test_throws ArgumentError plot_fit(quick.result; nsigma=Inf)
    @test_throws ArgumentError plot_fit(quick.result; limit_padding=-0.01)

    contour_values = collect(range(-2.0, 2.0; length=17))
    contour_delta = [xv^2 + 0.5 * xv * yv + yv^2 for xv in contour_values, yv in contour_values]
    contour_result = ContourResult((1, 2), contour_values, contour_values, contour_delta, contour_delta, [2.30, 6.18])
    @test plotting_extension._contour_level_label(2.30) ==
          "profile 1σ region (2 params, Δcost = 2.3)"
    @test plotting_extension._local_contour_label([2.30, 6.18]) ==
          "local covariance 1σ/2σ contours (parabolic approximation)"
    @test plotting_extension._local_contour_label([1.0, 4.0]) ==
          "local covariance contours (parabolic approximation)"
    contour_out = joinpath(mktempdir(), "contour_levels.png")
    heatmap_out = joinpath(mktempdir(), "contour_heatmap.png")
    failed_out = joinpath(mktempdir(), "contour_failed_refit.png")

    contour_figure = plot_contour(
        contour_result;
        filename=contour_out,
        format=:png,
        theme=:modern,
        appearance=:dark,
        local_covariance=[1.0 -0.25; -0.25 1.0],
        local_center=[0.0, 0.0],
    )
    @test contour_figure !== nothing
    contour_axis = only(filter(content -> content isa Axis, contour_figure.content))
    @test contour_axis.scene.viewport[].widths[1] >= 0.75 * size(contour_figure.scene)[1]

    # Descriptive confidence labels must not reduce the data-column width.
    contour_without_legend = plot_contour(contour_result; show_legend=false)
    save(joinpath(mktempdir(), "contour_without_legend.png"), contour_without_legend)
    contour_axis_without_legend = only(filter(content -> content isa Axis, contour_without_legend.content))
    @test contour_axis.scene.viewport[].widths[1] ==
          contour_axis_without_legend.scene.viewport[].widths[1]

    contour_right = plot_contour(
        contour_result;
        legend_position=:right,
        local_covariance=[1.0 -0.25; -0.25 1.0],
        local_center=[0.0, 0.0],
    )
    save(joinpath(mktempdir(), "contour_right_legend.png"), contour_right)
    contour_right_axis = only(filter(content -> content isa Axis, contour_right.content))
    @test contour_right_axis.scene.viewport[].widths[1] >= 0.5 * size(contour_right.scene)[1]
    @test plot_contour(contour_result; filename=heatmap_out, format=:png, show_heatmap=true) !== nothing
    profile_result = ProfileResult(1, contour_values, contour_values .^ 2, contour_values .^ 2, 1.0, 0.0)
    profile_out = joinpath(mktempdir(), "profile_legend.png")
    profile_figure = plot_profile(
        profile_result;
        filename=profile_out,
        format=:png,
        theme=:lab,
        appearance=:dark,
        local_sigma=1.0,
        delta_max=4.0,
        threshold_kwargs=(linewidth=2.5,),
    )
    @test profile_figure !== nothing
    profile_axis = only(filter(content -> content isa Axis, profile_figure.content))
    @test profile_axis.scene.viewport[].widths[1] >= 0.75 * size(profile_figure.scene)[1]
    @test_throws ArgumentError plot_profile(profile_result; delta_max=0.0)
    @test_throws ArgumentError plot_profile(profile_result; legend_position=:inside)
    @test_throws ArgumentError plot_contour(contour_result; legend_position=:inside)
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
    sans_matrix = plot_profile_matrix(precomputed_matrix; theme=:sans)
    tex_matrix = plot_profile_matrix(precomputed_matrix; theme=:tex)
    sans_matrix_axes = filter(content -> content isa Axis, sans_matrix.content)
    tex_matrix_axes = filter(content -> content isa Axis, tex_matrix.content)
    @test length(sans_matrix_axes) == 4
    @test length(tex_matrix_axes) == 4
    @test first(tex_matrix_axes).titlesize[] > first(sans_matrix_axes).titlesize[]
    @test minimum(axis.xticklabelsize[] for axis in sans_matrix_axes) >= 20
    @test minimum(axis.xticklabelsize[] for axis in tex_matrix_axes) >= 22
    @test_throws ArgumentError plot_profile_matrix(precomputed_matrix; parameter_names=["m"])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=Int[])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=[1, 1])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=[1, 3])
    @test_throws ArgumentError plot_profile_matrix(quick.result; parameters=[1, 2], panel_status_mode=:bad)

    residual_out = joinpath(mktempdir(), "residual_modern_dark.png")
    diagnostics_out = joinpath(mktempdir(), "diagnostics_analysis_dark.png")
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
