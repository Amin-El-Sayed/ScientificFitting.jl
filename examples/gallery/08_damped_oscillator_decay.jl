using JuFitter
using LaTeXStrings

const DATA_FILE = joinpath(@__DIR__, "..", "data", "damped_oscillator", "pohl_wheel_free_decay.csv")
const OUTPUT_DIR = joinpath(@__DIR__, "..", "output")
const DOC_ASSET_DIR = joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "gallery")

function load_damped_oscillator(path)
    rows = readlines(path)[2:end]
    time = Float64[]
    phi = Float64[]
    sigma_phi = Float64[]

    for row in rows
        fields = split(row, ",")
        length(fields) == 3 || continue
        push!(time, parse(Float64, fields[1]))
        push!(phi, parse(Float64, fields[2]))
        push!(sigma_phi, parse(Float64, fields[3]))
    end

    return time, phi, sigma_phi
end

time, phi, sigma_phi = load_damped_oscillator(DATA_FILE)
sigma_t = fill(0.0005, length(time))

# Underdamped oscillator solution with amplitude, frequency, phase, and damping.
model(t, p) = @. p[1] * exp(-p[4] * t) * cos(p[2] * t + p[3])

initial_guesses = [
    [1.9, 3.26, -0.7, 0.0035],
    [1.6, 3.20, 0.5, 0.0020],
    [2.2, 3.35, -2.0, 0.0060],
]

result = fit_model(
    model,
    time,
    phi;
    p0=initial_guesses[1],
    sigma_y=sigma_phi,
    sigma_x=sigma_t,
    bounds=([0.0, 2.0, -20.0, 0.0], [5.0, 5.0, 20.0, 0.05]),
    initial_guesses=initial_guesses,
    maxiters=2000,
)

mkpath(OUTPUT_DIR)
mkpath(DOC_ASSET_DIR)

parameter_names = [L"A", L"\omega", L"\phi_0", L"\lambda"]
model_label = L"\phi(t)=A e^{-\lambda t}\cos(\omega t+\phi_0)"

common_plot_kwargs = (
    parameter_names=parameter_names,
    title=L"\mathrm{Damped\ oscillator\ decay}",
    model_label=model_label,
    xlabel=L"t",
    xunit=L"\mathrm{s}",
    ylabel=L"\phi",
    yunit=L"\mathrm{rad}",
    latex_labels=true,
    latex_stats=true,
    stats_position=:right,
    stats_mode=:full,
    stats_sigdigits=5,
    stats_fontsize=13,
    figure_size=(1320, 760),
    band=:prediction,
    nsigma=1,
    band_label=L"1\sigma\ \mathrm{prediction\ band}",
    show_legend=true,
    legend_position=:lt,
    data_label=L"\mathrm{measurement}",
    fit_label=L"\mathrm{fit}",
    data_markersize=2.8,
    fit_linewidth=1.45,
    band_alpha=0.24,
)

plot_fit(
    result;
    common_plot_kwargs...,
    theme=:minimal,
    filename=joinpath(OUTPUT_DIR, "08_damped_oscillator_decay_light.png"),
    format=:png,
)

plot_fit(
    result;
    common_plot_kwargs...,
    theme=:dark,
    band_alpha=0.30,
    filename=joinpath(OUTPUT_DIR, "08_damped_oscillator_decay_dark.png"),
    format=:png,
)

plot_fit(
    result;
    common_plot_kwargs...,
    theme=:minimal,
    filename=joinpath(DOC_ASSET_DIR, "damped_oscillator_decay_light.png"),
    format=:png,
)

plot_fit(
    result;
    common_plot_kwargs...,
    theme=:dark,
    band_alpha=0.30,
    filename=joinpath(DOC_ASSET_DIR, "damped_oscillator_decay_dark.png"),
    format=:png,
)

println(report_text(result; parameter_names=["A", "omega", "phi0", "lambda"]))
