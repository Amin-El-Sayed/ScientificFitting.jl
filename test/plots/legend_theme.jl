using Test, ScientificFitting, CairoMakie, LaTeXStrings

@testset "Titled Makie legends inherit appearance" begin
    for style in (:sans, :tex), appearance in (:light, :dark)
        palette = plot_palette(style; appearance)
        with_theme(plot_theme(style; appearance)) do
            fig = Figure()
            for (i, title) in enumerate(("Group", L"\eta"))
                legend = Legend(fig[i,1], [LineElement()], ["component"], title)
                @test legend.titlecolor[] == palette.axis_color
                @test legend.labelcolor[] == palette.axis_color
            end
            custom = Legend(fig[3,1], [LineElement()], ["component"], "Group";
                titlecolor=:red)
            @test custom.titlecolor[] == :red
        end
    end
end

@testset "Diagnostic legends fit without squeezing the data" begin
    xs = collect(range(-3., 3.; length=17))
    delta = [x^2 + y^2 for x in xs, y in xs]
    result = ContourResult((1,2), xs, xs, delta, delta, [2.30,6.18])
    scan = ProfileResult(1, xs, xs.^2, xs.^2, 1., 0.)
    mktempdir() do directory
        for style in (:sans, :tex), appearance in (:light, :dark), position in (:below, :right)
            fig = plot_contour(result; theme=style, appearance,
                legend_position=position, local_covariance=[1. 0.; 0. 1.],
                local_center=[0.,0.], style=FitPlotStyle(figure_size=(820,700)))
            save(joinpath(directory, "$(style)_$(appearance)_$(position).svg"), fig)
            axis = only(filter(item -> item isa Axis, fig.content))
            legend = only(filter(item -> item isa Legend, fig.content))
            box = legend.layoutobservables.computedbbox[]
            @test all(box.origin .>= 0)
            @test all(box.origin .+ box.widths .<= collect(size(fig.scene)) .+ 1)
            @test axis.scene.viewport[].widths[1] >= 420
            @test axis.scene.viewport[].widths[2] >= 300
            @test all(collect(size(fig.scene)) .>= [820,700])
        end

        # Oversized custom typography must grow the same shared layout.
        fig = plot_profile(scan; theme=:tex, local_sigma=1.,
            style=FitPlotStyle(figure_size=(600,400)), legend_kwargs=(labelsize=44,))
        save(joinpath(directory, "large_profile_legend.svg"), fig)
        legend = only(filter(item -> item isa Legend, fig.content))
        box = legend.layoutobservables.computedbbox[]
        @test all(box.origin .+ box.widths .<= collect(size(fig.scene)) .+ 1)
        @test only(filter(item -> item isa Axis, fig.content)).scene.viewport[].widths[1] >= 420
        @test only(filter(item -> item isa Axis, fig.content)).scene.viewport[].widths[2] >= 300
    end
end
