using CairoMakie, LaTeXStrings, Printf

"""Compare fitted CMS bin means and residuals; this function never runs a fit."""
function cms_mass_plot(results, edges, counts; theme=:sans, appearance=:light, show_panel=true)
    pal = plot_palette(theme; appearance)
    colors = pal.series_colors[1:3]
    centers = (edges[1:end-1] + edges[2:end])/2
    step_x = repeat(edges; inner=2)[2:end-1]
    means = [cms_binmeans(result, edges) for result in results]
    return with_theme(plot_theme(theme; appearance)) do
        fig = Figure(size=(1060, show_panel ? 1100 : 930), backgroundcolor=pal.background_color)
        ax = Axis(fig[1,1]; title="CMS dimuons: model checks", ylabel="candidates / bin")
        # Counts exceed 1000 in every bin; sqrt(n) bars are the large-count approximation.
        errorbars!(ax, centers, counts, sqrt.(counts); color=pal.yerr_color,
            whiskerwidth=pal.error_whiskerwidth)
        scatter!(ax, centers, counts; color=pal.data_color,
            markersize=pal.data_markersize, label="data")
        for (mean, color, label) in zip(means, colors,
                ("A: exponential", "B: quadratic", "C: two widths"))
            lines!(ax, step_x, repeat(mean; inner=2); color, label)
        end
        Legend(fig[0,1], ax; nbanks=2, tellwidth=false)
        axes = [ax]
        hidexdecorations!(ax; grid=false)
        for (i, mean) in enumerate(means)
            label = theme==:tex ? (L"r_A", L"r_B", L"r_C")[i] : "r ($(('A','B','C')[i]))"
            rx = Axis(fig[i+1,1]; ylabel=label,
                xlabel=theme==:tex ? L"m_{\mu\mu}\,/\,(\mathrm{GeV}/c^2)" : "dimuon mass / (GeV/c^2)")
            # +/-2 is a pointwise visual reference, not a simultaneous acceptance band.
            hspan!(rx, -2., 2.; color=(pal.reference_color, .10))
            barplot!(rx, centers, cms_residuals(counts, mean); width=.0017, color=colors[i])
            hlines!(rx, [0.]; color=pal.stats_color)
            ylims!(rx, -5.8, 5.8)
            i < 3 && hidexdecorations!(rx; grid=false)
            rowsize!(fig.layout, i+1, Auto(1))
            push!(axes, rx)
        end
        linkxaxes!(axes...)
        xlims!(ax, first(edges), last(edges))
        ylims!(ax, 0, 1.08maximum(counts))
        rowsize!(fig.layout, 1, Auto(3))
        if show_panel
            stats = [@sprintf("%s: D / ndf = %.1f / %d", label, r.stats.chi2, r.stats.ndf)
                for (label,r) in zip(("A", "B", "C"), results)]
            plot_info_panel!(fig[5,1]; theme, appearance,
                title="2,551,454 candidates; all three models leave structure",
                statistic_lines=stats)
        end
        resize_plot_to_layout!(fig; minimum_axis_size=(700, nothing))
        return fig
    end
end

"""Show data shapes by muon geometry and the conditional core-scale estimates."""
function cms_kinematic_plot(results, edges, groups; theme=:sans, appearance=:light, show_panel=true)
    pal = plot_palette(theme; appearance)
    labels = ["[0, 0.9)", "[0.9, 1.4)", "[1.4, infinity)"]
    eta_label = theme==:tex ? L"\max(|\eta_1|,|\eta_2|)" : "max(|eta_1|, |eta_2|)"
    step_x = repeat(edges; inner=2)[2:end-1]
    return with_theme(plot_theme(theme; appearance)) do
        fig = Figure(size=(1020, show_panel ? 900 : 770), backgroundcolor=pal.background_color)
        ax = Axis(fig[1,1]; title="Different kinematics, different mass shapes",
            ylabel="fraction / bin", xlabel=theme==:tex ?
                L"m_{\mu\mu}\,/\,(\mathrm{GeV}/c^2)" : "dimuon mass / (GeV/c^2)")
        for i in eachindex(results)
            lines!(ax, step_x, repeat(groups[:,i]/sum(groups[:,i]); inner=2);
                color=pal.series_colors[i], label=labels[i])
        end
        xlims!(ax, 2.95, 3.25)  # zoom only the display; every fit uses [2.8, 3.4)
        Legend(fig[0,1], ax, eta_label; nbanks=3, tellwidth=false)
        rx = Axis(fig[2,1]; xticks=(1:3, labels), xlabel=eta_label,
            ylabel=theme==:tex ? L"\sigma_{\mathrm{core}}\,/\,(\mathrm{MeV}/c^2)" :
                "core scale / (MeV/c^2)")
        scales = [1000r.params[2] for r in results]
        errors = [1000r.param_stderr[2] for r in results]
        errorbars!(rx, 1:3, scales, errors; color=pal.yerr_color,
            whiskerwidth=pal.error_whiskerwidth)
        scatter!(rx, 1:3, scales; color=collect(pal.series_colors[1:3]), markersize=pal.data_markersize)
        xlims!(rx, .6, 3.4)
        rowsize!(fig.layout, 1, Auto(2))
        rowsize!(fig.layout, 2, Auto(1))
        if show_panel
            plot_info_panel!(fig[3,1]; theme, appearance,
                title="Model A in each group; local errors only",
                statistic_lines=[@sprintf("%s: D / ndf = %.1f / %d", label, r.stats.chi2, r.stats.ndf)
                    for (label,r) in zip(labels, results)])
        end
        resize_plot_to_layout!(fig; minimum_axis_size=(700, nothing))
        return fig
    end
end
