using JuFitter
using LaTeXStrings
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

# A controlled style comparison: scientific content stays fixed while the three
# production-oriented visual contracts change.
x = collect(range(0.0, 10.0; length=90))
sigma_y = 0.12 .+ 0.01 .* x
y = @. 1.85 * x + 0.7 + sigma_y * sin(1.6 * x)

for style in (:workbench, :showcase, :publication)
    typography = if style == :publication
        (
            title=L"\mathrm{Controlled\ style\ comparison}",
            model_label=L"U(t)=m t+b",
            xlabel=L"t",
            xunit=L"\mathrm{s}",
            ylabel=L"U",
            yunit=L"\mathrm{V}",
            parameter_names=[L"m", L"b"],
            latex_labels=true,
            latex_stats=true,
            band_label=L"1\sigma\ \mathrm{prediction\ band}",
        )
    else
        (
            title="Controlled style comparison",
            model_label="U(t) = m t + b",
            xlabel="t",
            xunit="s",
            ylabel="U",
            yunit="V",
            parameter_names=["m", "b"],
            latex_labels=false,
            latex_stats=false,
            band_label="1σ prediction band",
        )
    end
    fitplot(
        x,
        y;
        typography...,
        sigma_y=sigma_y,
        filename=example_output("07_style_$(style).pdf"),
        theme=style,
        band=:prediction,
        nsigma=1,
        show_legend=true,
        stats_position=:right,
        stats_mode=:full,
        report=:plot,
    )
end

println("Saved controlled style comparison to ", example_output("07_style_*.pdf"))
