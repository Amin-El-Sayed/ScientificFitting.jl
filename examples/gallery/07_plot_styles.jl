using JuFitter
using LaTeXStrings
include(joinpath(@__DIR__, "..", "_example_utils.jl"))

# The same fit rendered in the three main visual modes:
# clean default, dense-data minimal, and LaTeX-style paper output.
x = collect(range(0.0, 10.0; length=90))
sigma_y = 0.12 .+ 0.01 .* x
y = @. 1.85 * x + 0.7 + sigma_y * sin(1.6 * x)

fitplot(
    x,
    y;
    sigma_y=sigma_y,
    filename=example_output("07_style_clean.pdf"),
    theme=:clean,
    xlabel="time",
    xunit="s",
    ylabel="signal",
    yunit="V",
    parameter_names=["slope", "offset"],
    report=:plot,
)

fitplot(
    x,
    y;
    sigma_y=sigma_y,
    filename=example_output("07_style_minimal.pdf"),
    theme=:minimal,
    xlabel="time",
    xunit="s",
    ylabel="signal",
    yunit="V",
    parameter_names=["slope", "offset"],
    report=:plot,
)

fitplot(
    x,
    y;
    sigma_y=sigma_y,
    filename=example_output("07_style_paper.pdf"),
    theme=:paper,
    latex_stats=true,
    xlabel=L"t",
    xunit=L"\mathrm{s}",
    ylabel=L"U",
    yunit=L"\mathrm{V}",
    parameter_names=[L"m", L"b"],
    report=:plot,
)

println("Saved style examples to ", example_output("07_style_*.pdf"))
