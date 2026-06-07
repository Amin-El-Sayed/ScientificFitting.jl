const MULTI_SNAPSHOT_ONLY = get(ENV, "JUFITTER_DOC_SNAPSHOT_ONLY", "0") == "1"
const MULTI_RENDER_DOC_ASSETS = !MULTI_SNAPSHOT_ONLY

if MULTI_RENDER_DOC_ASSETS
    using CairoMakie
end

using JuFitter
using LinearAlgebra
using Printf

const MULTI_OUTPUT_DIR = joinpath(@__DIR__, "..", "output")
const MULTI_DOC_ASSET_DIR = joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "gallery")
const MULTI_EMIT_DOC_OUTPUT_SNAPSHOTS = get(ENV, "JUFITTER_DOC_OUTPUT_SNAPSHOTS", "0") == "1"

function emit_multi_doc_output_snapshot(body::Function, id::AbstractString)
    MULTI_EMIT_DOC_OUTPUT_SNAPSHOTS || return nothing

    println("=== JUFITTER_DOC_OUTPUT_BEGIN ", id, " ===")
    body()
    println("=== JUFITTER_DOC_OUTPUT_END ", id, " ===")
    return nothing
end

linear_channel(x, p) = @. p[1] * x + p[2]

# Controlled three-channel calibration record with realistic heteroscedastic scatter.
x_a = collect(0.0:1.0:10.0)
x_b = collect(0.5:1.0:9.5)
x_c = collect(0.0:1.25:10.0)

sigma_a = @. 0.075 + 0.008 * x_a
sigma_b = @. 0.085 + 0.006 * x_b
sigma_c = @. 0.080 + 0.009 * x_c

pattern_a = 1.65 .* [0.2, -0.7, 0.4, -0.5, 0.8, -0.3, 0.1, -0.8, 0.6, -0.2, 0.5]
pattern_b = 1.65 .* [-0.4, 0.7, -0.5, 0.1, 0.8, -0.6, 0.3, -0.2, 0.6, -0.3]
pattern_c = 1.65 .* [0.3, -0.6, 0.7, -0.4, 0.1, 0.8, -0.5, 0.4, -0.7]

y_a = linear_channel(x_a, [1.82, 0.72]) .+ sigma_a .* pattern_a
y_b = linear_channel(x_b, [1.82, -0.46]) .+ sigma_b .* pattern_b
y_c = linear_channel(x_c, [1.91, 0.12]) .+ sigma_c .* pattern_c

models = [linear_channel, linear_channel, linear_channel]
x_sets = [x_a, x_b, x_c]
y_sets = [y_a, y_b, y_c]
sigma_sets = [sigma_a, sigma_b, sigma_c]

# Hypothesis one: all three channels share one physical gain.
all_shared_result = fit_multi_model(
    models,
    x_sets,
    y_sets;
    p0=[1.85, 0.7, -0.4, 0.1],
    sigma_y=sigma_sets,
    parameter_map=[[1, 2], [1, 3], [1, 4]],
    parameter_names=["shared gain", "offset A", "offset B", "offset C"],
)

# Hypothesis two: A and B share a gain, while channel C has its own gain.
partial_shared_result = fit_multi_model(
    models,
    x_sets,
    y_sets;
    p0=[1.82, 0.7, -0.4, 1.90, 0.1],
    sigma_y=sigma_sets,
    parameter_map=[[1, 2], [1, 3], [4, 5]],
    parameter_names=["gain A/B", "offset A", "offset B", "gain C", "offset C"],
)

gain_gap_gradient = [-1.0, 0.0, 0.0, 1.0, 0.0]
gain_gap = partial_shared_result.params[4] - partial_shared_result.params[1]
sigma_gain_gap = sqrt(dot(
    gain_gap_gradient,
    partial_shared_result.param_covariance * gain_gap_gradient,
))

function multi_theme(dark::Bool)
    foreground = dark ? "#edf2f4" : "#14151a"
    muted = dark ? "#b8c1ca" : "#5b6270"
    grid = dark ? ("#2a313a", 0.85) : ("#e9eef4", 0.95)
    background = dark ? "#111318" : "#ffffff"
    return Theme(
        fontsize=20,
        font="TeX Gyre Heros",
        Figure=(backgroundcolor=background,),
        Axis=(
            backgroundcolor=background,
            xlabelsize=27,
            ylabelsize=27,
            titlesize=26,
            xticklabelsize=20,
            yticklabelsize=20,
            xlabelcolor=foreground,
            ylabelcolor=foreground,
            titlecolor=foreground,
            xticklabelcolor=muted,
            yticklabelcolor=muted,
            xtickcolor=muted,
            ytickcolor=muted,
            xgridcolor=grid,
            ygridcolor=grid,
            topspinevisible=false,
            rightspinevisible=false,
            leftspinecolor=muted,
            bottomspinecolor=muted,
        ),
        Legend=(framevisible=false, labelcolor=foreground, labelsize=19, patchsize=(30, 16)),
    )
end

function confidence_band(result, parameter_indices, x)
    local_covariance = result.param_covariance[parameter_indices, parameter_indices]
    jacobian = hcat(x, ones(length(x)))
    variance = vec(sum((jacobian * local_covariance) .* jacobian; dims=2))
    return sqrt.(clamp.(variance, 0.0, Inf))
end

function dataset_pulls(result, parameter_maps)
    return [
        (y_sets[i] .- linear_channel(x_sets[i], result.params[parameter_maps[i]])) ./ sigma_sets[i]
        for i in eachindex(x_sets)
    ]
end

function add_pull_reference!(axis, xmin, xmax, color_1sigma, color_2sigma, zero_color)
    band!(axis, [xmin, xmax], [-2.0, -2.0], [2.0, 2.0]; color=color_2sigma)
    band!(axis, [xmin, xmax], [-1.0, -1.0], [1.0, 1.0]; color=color_1sigma)
    hlines!(axis, [0.0]; color=zero_color, linewidth=1.4)
end

fmt(x, digits=4) = @sprintf("%.*g", digits, x)

function save_multi_dataset_calibration(
    filename;
    dark::Union{Nothing, Bool}=nothing,
    style::Symbol=:showcase,
    appearance::Symbol=dark === nothing ? :light : (dark ? :dark : :light),
)
    MULTI_RENDER_DOC_ASSETS || return nothing

    dark_mode = appearance == :dark
    palette = plot_palette(style; appearance=appearance)
    foreground = palette.stats_color
    muted = palette.stats_muted_color
    background = dark_mode ? "#111318" : (style == :showcase ? "#fbfcfd" : "#ffffff")
    colors = if style == :publication
        dark_mode ? ["#edf2f4", "#b8c1ca", "#8d96a3"] : ["#101216", "#606874", "#8a929c"]
    elseif style == :showcase
        dark_mode ? [palette.fit_color, "#f0b35f", "#d686bd"] : [palette.fit_color, "#b05a36", "#9b4d86"]
    else
        dark_mode ? [palette.fit_color, "#c8a04d", "#af7ac5"] : [palette.fit_color, "#8a6f22", "#6e5aae"]
    end
    bands = [(colors[i], style == :publication ? 0.10 : (dark_mode ? 0.14 : 0.24)) for i in eachindex(colors)]
    pull_1sigma = (palette.band_color, dark_mode ? 0.13 : 0.24)
    pull_2sigma = (palette.band_color, dark_mode ? 0.06 : 0.12)
    markers = [:circle, :rect, :diamond]
    labels = ["channel A", "channel B", "channel C"]
    partial_maps = [[1, 2], [1, 3], [4, 5]]
    all_shared_maps = [[1, 2], [1, 3], [1, 4]]

    figure = with_theme(plot_theme(style; appearance=appearance)) do
        Figure(size=(1680, 1040), backgroundcolor=background)
    end
    fit_axis = Axis(
        figure[1, 1];
        title="Three-channel calibration transfer",
        ylabel="channel response y (V)",
    )
    shared_pull_axis = Axis(
        figure[2, 1];
        title="Pulls: all channels forced to share one gain",
        titlealign=:left,
        titlesize=19,
        ylabel="pull rᵢ",
        ylabelsize=22,
    )
    partial_pull_axis = Axis(
        figure[3, 1];
        title="Pulls: channels A/B share a gain; C is independent",
        titlealign=:left,
        titlesize=19,
        xlabel="reference input x",
        ylabel="pull rᵢ",
        ylabelsize=22,
    )

    x_grid = collect(range(0.0, 10.0; length=300))
    for i in eachindex(x_sets)
        partial_parameters = partial_shared_result.params[partial_maps[i]]
        partial_prediction = linear_channel(x_grid, partial_parameters)
        partial_sigma = confidence_band(partial_shared_result, partial_maps[i], x_grid)
        all_shared_prediction = linear_channel(x_grid, all_shared_result.params[all_shared_maps[i]])

        band!(
            fit_axis,
            x_grid,
            partial_prediction .- partial_sigma,
            partial_prediction .+ partial_sigma;
            color=bands[i],
        )
        lines!(
            fit_axis,
            x_grid,
            all_shared_prediction;
            color=(colors[i], 0.48),
            linestyle=:dash,
            linewidth=max(1.5, palette.fit_linewidth - 0.5),
        )
        lines!(
            fit_axis,
            x_grid,
            partial_prediction;
            color=colors[i],
            linewidth=palette.fit_linewidth,
        )
        errorbars!(
            fit_axis,
            x_sets[i],
            y_sets[i],
            sigma_sets[i];
            color=(colors[i], 0.50),
            whiskerwidth=palette.error_whiskerwidth,
        )
        scatter!(
            fit_axis,
            x_sets[i],
        y_sets[i];
            color=colors[i],
            marker=markers[i],
            markersize=palette.data_markersize + 1.2,
        )
    end
    hidexdecorations!(fit_axis; grid=false)

    shared_pulls = dataset_pulls(all_shared_result, all_shared_maps)
    partial_pulls = dataset_pulls(partial_shared_result, partial_maps)
    add_pull_reference!(shared_pull_axis, 0.0, 10.0, pull_1sigma, pull_2sigma, (foreground, 0.55))
    add_pull_reference!(partial_pull_axis, 0.0, 10.0, pull_1sigma, pull_2sigma, (foreground, 0.55))
    for i in eachindex(x_sets)
        lines!(shared_pull_axis, x_sets[i], shared_pulls[i]; color=(colors[i], 0.42), linewidth=1.2)
        scatter!(shared_pull_axis, x_sets[i], shared_pulls[i]; color=colors[i], marker=markers[i], markersize=6)
        lines!(partial_pull_axis, x_sets[i], partial_pulls[i]; color=(colors[i], 0.42), linewidth=1.2)
        scatter!(partial_pull_axis, x_sets[i], partial_pulls[i]; color=colors[i], marker=markers[i], markersize=6)
    end
    hidexdecorations!(shared_pull_axis; grid=false)
    linkxaxes!(fit_axis, shared_pull_axis, partial_pull_axis)
    ylims!(shared_pull_axis, -3.4, 3.4)
    ylims!(partial_pull_axis, -3.4, 3.4)

    side = GridLayout()
    figure[1:3, 2] = side
    legend_elements = LegendElement[
        MarkerElement(color=colors[i], marker=markers[i], markersize=12)
        for i in eachindex(labels)
    ]
    append!(
        legend_elements,
        [
            LineElement(color=foreground, linewidth=palette.fit_linewidth),
            LineElement(color=(foreground, 0.55), linestyle=:dash, linewidth=max(1.5, palette.fit_linewidth - 0.5)),
            PolyElement(color=(palette.band_color, dark_mode ? 0.18 : 0.28)),
        ],
    )
    Legend(
        side[1, 1],
        legend_elements,
        [
            labels...,
            "partial-sharing model",
            "all-shared-gain hypothesis",
            "local 1σ fit band",
        ];
        framevisible=false,
        tellheight=true,
        nbanks=2,
        labelsize=palette.stats_fontsize + 1,
    )
    report_axis = Axis(side[2, 1]; backgroundcolor=:transparent)
    hidedecorations!(report_axis)
    hidespines!(report_axis)
    text!(
        report_axis,
        0,
        1;
        text="partial-sharing result\n" *
             "  gain A/B = $(fmt(partial_shared_result.params[1], 6)) ± $(fmt(partial_shared_result.param_stderr[1], 2))\n" *
             "  gain C = $(fmt(partial_shared_result.params[4], 6)) ± $(fmt(partial_shared_result.param_stderr[4], 2))\n" *
             "  gain C − gain A/B\n" *
             "    $(fmt(gain_gap, 5)) ± $(fmt(sigma_gain_gap, 2))\n" *
             "    $(fmt(gain_gap / sigma_gain_gap, 3))σ from zero\n\n" *
             "all-shared-gain hypothesis\n" *
             "  χ²/ndf = $(fmt(all_shared_result.stats.chi2_ndf, 4))\n" *
             "  P(χ²) = $(fmt(all_shared_result.stats.pvalue, 4))\n" *
             "  AIC = $(fmt(all_shared_result.stats.aic, 5))\n\n" *
             "partial-sharing model\n" *
             "  χ²/ndf = $(fmt(partial_shared_result.stats.chi2_ndf, 4))\n" *
             "  P(χ²) = $(fmt(partial_shared_result.stats.pvalue, 4))\n" *
             "  AIC = $(fmt(partial_shared_result.stats.aic, 5))\n\n" *
             "ΔAIC = $(fmt(all_shared_result.stats.aic - partial_shared_result.stats.aic, 5))\n" *
             "Do not transfer channel C's gain.",
        space=:relative,
        align=(:left, :top),
        color=foreground,
        fontsize=palette.stats_fontsize + 2,
        lineheight=1.08,
    )

    rowsize!(side, 1, Auto())
    rowsize!(side, 2, Relative(1))
    rowsize!(figure.layout, 1, Relative(0.58))
    rowsize!(figure.layout, 2, Relative(0.21))
    rowsize!(figure.layout, 3, Relative(0.21))
    colsize!(figure.layout, 2, Fixed(600))
    save(filename, figure)
end

if MULTI_RENDER_DOC_ASSETS
    mkpath(MULTI_OUTPUT_DIR)
    mkpath(MULTI_DOC_ASSET_DIR)

    for (dark, suffix) in ((false, "light"), (true, "dark"))
        save_multi_dataset_calibration(
            joinpath(MULTI_OUTPUT_DIR, "10_multi_dataset_calibration_$(suffix).png");
            dark=dark,
        )
        save_multi_dataset_calibration(
            joinpath(MULTI_DOC_ASSET_DIR, "multi_dataset_shared_slope_$(suffix).png");
            dark=dark,
        )
    end

    for style in (:workbench, :showcase, :publication), appearance in (:light, :dark)
        save_multi_dataset_calibration(
            joinpath(MULTI_DOC_ASSET_DIR, "multi_dataset_shared_slope_$(style)_$(appearance).png");
            style=style,
            appearance=appearance,
        )
    end
end

println("All-shared-gain hypothesis")
println(report_text(all_shared_result))
println()
println("Partial-sharing model")
println(report_text(partial_shared_result))
println("gain C - gain A/B = ", gain_gap, " +/- ", sigma_gain_gap)
println()
println("All-shared diagnostic dashboard")
println(diagnostic_dashboard_text(all_shared_result))
println("Partial-sharing diagnostic dashboard")
println(diagnostic_dashboard_text(partial_shared_result))
emit_multi_doc_output_snapshot("multi_dataset") do
    println("All-shared-gain hypothesis")
    println(report_text(all_shared_result))
    println()
    println("Partial-sharing model")
    println(report_text(partial_shared_result))
    println("gain C - gain A/B = ", gain_gap, " +/- ", sigma_gain_gap)
    println()
    println("All-shared diagnostic dashboard")
    println(diagnostic_dashboard_text(all_shared_result))
    println("Partial-sharing diagnostic dashboard")
    println(diagnostic_dashboard_text(partial_shared_result))
end
