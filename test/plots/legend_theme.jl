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
