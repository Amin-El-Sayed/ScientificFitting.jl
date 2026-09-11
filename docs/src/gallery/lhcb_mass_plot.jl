using CairoMakie, LaTeXStrings, ForwardDiff, Printf

"""Plot already-fitted bin expectations and deviance residuals; never refit."""
function lhcb_mass_plot(result, constructor, edges, counts;
        theme=:sans, appearance=:light, show_panel=true)
    names = keys(parameter_values(result))
    build(p) = build_model(constructor, NamedTuple{names}(Tuple(p)))
    binmeans(model) = sum(n .* diff(cdf.(d, edges))
        for (d, n) in zip(components(model), DistributionsHEP.yields(model)))
    means = binmeans(fitted_model(result))
    jac = ForwardDiff.jacobian(p -> binmeans(build(p)), result.params)
    sigma = sqrt.(max.(vec(sum((jac * result.param_covariance) .* jac; dims=2)), 0))
    residual = [sign(n-m)*sqrt(max(2*(m-n+(n==0 ? 0 : n*log(n/m))), 0))
        for (n,m) in zip(counts, means)]
    centers = (edges[1:end-1] + edges[2:end])/2
    # Garwood intervals for bin means are asymmetric, including at small counts.
    lower = [n==0 ? 0. : quantile(Chisq(2n), 0.158655254)/2 for n in counts]
    upper = [quantile(Chisq(2(n+1)), 0.841344746)/2 for n in counts]
    step_x = repeat(edges; inner=2)[2:end-1]
    pal = plot_palette(theme; appearance)
    return with_theme(plot_theme(theme; appearance)) do
        fig = Figure(size=(show_panel ? 1180 : 900, 740), backgroundcolor=pal.background_color)
        ax = Axis(fig[1,1]; title="LHCb open data: three-kaon mass",
            ylabel="candidates / bin")
        band!(ax, step_x, repeat(means-sigma; inner=2), repeat(means+sigma; inner=2);
            color=(pal.band_color, 0.28), label="local 1-sigma mean band")
        lines!(ax, step_x, repeat(means; inner=2);
            color=pal.fit_color, label="two-width peak + background")
        background = last(DistributionsHEP.yields(fitted_model(result))) .*
            diff(cdf.(last(components(fitted_model(result))), edges))
        lines!(ax, step_x, repeat(background; inner=2);
            color=pal.reference_color, linestyle=:dash, label="background")
        errorbars!(ax, centers, counts, counts-lower, upper-counts;
            color=pal.yerr_color, whiskerwidth=pal.error_whiskerwidth)
        scatter!(ax, centers, counts; color=pal.data_color,
            markersize=pal.data_markersize, label="data (68% Poisson intervals)")
        hidexdecorations!(ax; grid=false)
        rx = Axis(fig[2,1]; xlabel=theme==:tex ? L"m_{KKK}\,/\,\mathrm{MeV}\,c^{-2}" :
            "three-kaon mass / (MeV/c^2)", ylabel="deviance residual")
        barplot!(rx, centers, residual; width=4., color=pal.fit_color)
        hlines!(rx, [0.]; color=pal.stats_color)
        hlines!(rx, [-2., 2.]; color=pal.reference_color, linestyle=:dash)
        linkxaxes!(ax, rx)
        xlims!(ax, first(edges), last(edges))
        ylims!(ax, 0, 1.08maximum(upper))
        rowsize!(fig.layout, 1, Auto(3))
        rowsize!(fig.layout, 2, Auto(1))
        if show_panel
            p, e = result.params, result.param_stderr
            lines = theme==:tex ? Any[
                LaTeXString(@sprintf("N_s = %.0f \\pm %.0f", p[6], e[6])),
                LaTeXString(@sprintf("N_b = %.0f \\pm %.0f", p[7], e[7])),
                LaTeXString(@sprintf("\\mu = %.2f \\pm %.2f\\;\\mathrm{MeV}/c^2", p[1], e[1])),
                LaTeXString(@sprintf("\\sigma = %.2f \\pm %.2f\\;\\mathrm{MeV}/c^2", p[2], e[2])),
            ] : Any[@sprintf("signal = %.0f +/- %.0f candidates", p[6], e[6]),
                @sprintf("background = %.0f +/- %.0f candidates", p[7], e[7]),
                @sprintf("centroid = %.2f +/- %.2f MeV/c^2", p[1], e[1]),
                @sprintf("core width = %.2f +/- %.2f MeV/c^2", p[2], e[2])]
            plot_info_panel!(fig[1:2,2]; theme, appearance, legend_source=ax,
                title="Window-conditional fit", parameter_lines=lines,
                statistic_lines=[@sprintf("D / ndf = %.2f / %d", result.stats.chi2, result.stats.ndf),
                    @sprintf("asymptotic p = %.3f", result.stats.pvalue)])
        else
            Legend(fig[3,1], ax; nbanks=2, tellwidth=false)
        end
        resize_plot_to_layout!(fig; minimum_axis_size=(600, nothing))
        fig
    end
end
