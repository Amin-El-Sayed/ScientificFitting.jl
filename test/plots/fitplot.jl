using JuFitter
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

    minimal_out = joinpath(mktempdir(), "fitplot_minimal.png")
    paper_out = joinpath(mktempdir(), "fitplot_paper.png")

    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=minimal_out,
        format=:png,
        report=:none,
        theme=:minimal,
        parameter_names=["m", "b"],
    )
    fitplot(
        x,
        y;
        sigma_y=sigma_y,
        filename=paper_out,
        format=:png,
        report=:none,
        theme=:paper,
        parameter_names=["m", "b"],
        model_label="y = m x + b",
    )

    @test isfile(minimal_out)
    @test isfile(paper_out)
end
