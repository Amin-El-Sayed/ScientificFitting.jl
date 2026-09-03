using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")
const GALLERY_SRC = joinpath(DOCS_SRC, "gallery")

const GALLERY_PAGES = [
    "linear_calibration.md",
    "xy_uncertainties.md",
    "full_covariance.md",
    "resonance_decay.md",
    "photoelectric_threshold.md",
    "constraints_profiles.md",
    "poisson_histogram.md",
    "multi_dataset.md",
]

function read_gallery_page(name)
    path = joinpath(GALLERY_SRC, name)
    return path, read(path, String)
end

function has_any(text, patterns)
    return any(pattern -> occursin(pattern, text), patterns)
end

function image_sources(text)
    return [match.captures[1] for match in eachmatch(r"<img[^>]+src=\"([^\"]+)\"", text)]
end

function style_groups(text)
    groups = Dict{String,Set{Tuple{String,String,String}}}()
    for img_match in eachmatch(r"<img([^>]+)>", text)
        tag = img_match.captures[1]
        group_match = match(r"data-scientificfitting-plot-group=\"([^\"]+)\"", tag)
        style_match = match(r"data-scientificfitting-plot-style=\"([^\"]+)\"", tag)
        src_match = match(r"src=\"([^\"]+)\"", tag)
        any(isnothing, (group_match, style_match, src_match)) && continue

        group = group_match.captures[1]
        style = style_match.captures[1]
        src = src_match.captures[1]
        panel_match = match(r"data-scientificfitting-plot-panel=\"([^\"]+)\"", tag)
        panel = isnothing(panel_match) ? "" : panel_match.captures[1]
        appearance = occursin("_dark", src) ? "dark" : occursin("_light", src) ? "light" : ""
        get!(groups, group, Set{Tuple{String,String,String}}())
        push!(groups[group], (style, panel, appearance))
    end
    return groups
end

function resolve_doc_asset(page_path, src)
    return normpath(joinpath(dirname(page_path), src))
end

@testset "Documentation gallery release gate" begin
    required_style_pairs = Set(
        (style, appearance)
        for style in ("sans", "tex")
        for appearance in ("light", "dark")
    )
    required_panel_variants = Set(
        (style, panel, appearance)
        for style in ("sans", "tex")
        for panel in ("show", "hide")
        for appearance in ("light", "dark")
    )

    for page in GALLERY_PAGES
        path, text = read_gallery_page(page)

        @testset "$page" begin
            @test occursin(r"^# "m, text)
            @test has_any(text, [r"^## .*Question"m, r"^## Scientific Question"m])
            @test has_any(text, [r"^## Data"m, r"^## The Measurement"m])
            @test has_any(text, [r"^## .*Model"m, r"^## .*Cost"m, r"^## Poisson Likelihood"m])
            @test has_any(text, [r"```julia"m, r"^## Complete .*Code"m, r"^## Complete .*Fit"m])
            @test occursin("scientificfitting-cell-output", text)
            @test occursin("Fit diagnostic dashboard", text)
            @test has_any(text, [r"^## Diagnostics"m, r"^## What To Inspect"m, r"^## Diagnose"m])
            @test has_any(text, [r"^## Interpretation"m, r"^## Decision"m, r"^## Read "m, r"^## Reading "m, r"^## Why Compare"m])
            @test has_any(text, [r"^## What Can Go Wrong"m, r"^## Failure Modes"m, r"^## What To Do Before"m])
            @test occursin("1σ", text) || occursin("1-sigma", text) || occursin("one-sigma", text)
            @test !occursin("P1", text)
            @test !occursin("Praktikum", text)
            @test !occursin("placeholder prose", lowercase(text))
            @test !occursin("being rewritten", text)
            @test !occursin("not all of them are finished", text)

            sources = image_sources(text)
            @test !isempty(sources)
            for src in sources
                @test isfile(resolve_doc_asset(path, src))
            end

            groups = style_groups(text)
            @test !isempty(groups)
            for (_, variants) in groups
                style_pairs = Set((style, appearance) for (style, _, appearance) in variants)
                @test required_style_pairs ⊆ style_pairs

                panels = Set(panel for (_, panel, _) in variants if !isempty(panel))
                if !isempty(panels)
                    @test panels == Set(("show", "hide"))
                    @test required_panel_variants ⊆ variants
                end
            end
        end
    end

    overview = read(joinpath(DOCS_SRC, "gallery.md"), String)
    @test !occursin("being rewritten", overview)
    @test !occursin("remaining documentation passes", overview)
    @test !occursin("not all of them are finished", overview)
    @test !occursin("toy-like synthetic examples", overview)

    expected_path = [
        "Linear calibration",
        "XY uncertainties",
        "Full covariance",
        "Damped oscillator",
        "Photoelectric work function",
        "Poisson and histograms",
        "Constraints and profiles",
        "Multi-dataset fit",
    ]
    card_positions = [findfirst(">$(title)</a></h3>", overview) for title in expected_path]
    @test all(!isnothing, card_positions)
    @test issorted(first.(card_positions))

    make_source = read(joinpath(ROOT, "docs", "make.jl"), String)
    navigation_pages = [
        "linear_calibration",
        "xy_uncertainties",
        "full_covariance",
        "resonance_decay",
        "photoelectric_threshold",
        "poisson_histogram",
        "constraints_profiles",
        "multi_dataset",
    ]
    navigation_positions = [
        findfirst("gallery/$(page).md", make_source) for page in navigation_pages
    ]
    @test all(!isnothing, navigation_positions)
    @test issorted(first.(navigation_positions))

    @test occursin("Simple fits stay simple", overview)
    @test occursin("actual program output", overview)
    @test !occursin("## Recommended Path", overview)
    @test !occursin("## What Each Example Teaches", overview)

    full_covariance = read(joinpath(GALLERY_SRC, "full_covariance.md"), String)
    @test occursin("exp(-p[2] * t)", full_covariance)
    @test !occursin("exp(p[2] * t)", full_covariance)
    @test occursin("sigma_stat^2 * (i == j)", full_covariance)
    @test occursin("sigma_corr^2 * exp", full_covariance)
    @test occursin("filename=\"full_covariance_decay.pdf\"", full_covariance)

    photoelectric = read(joinpath(GALLERY_SRC, "photoelectric_threshold.md"), String)
    @test occursin("photoelectric_slope = me - mb", photoelectric)
    @test occursin("(nu_THz - reference_frequency_THz)", photoelectric)
    @test occursin("h_fit = photoelectric_slope * elementary_charge", photoelectric)
    @test occursin("work_function_eV = photoelectric_slope * threshold_THz", photoelectric)

    resonance = read(joinpath(GALLERY_SRC, "resonance_decay.md"), String)
    resonance_script = read(
        joinpath(ROOT, "examples", "gallery", "08_damped_oscillator_decay.jl"),
        String,
    )
    @test occursin("sigma_angle = data.sigma_angle", resonance)
    @test occursin("sigma_angle = data.sigma_angle", resonance_script)
    @test !occursin("0.5 .* data.sigma_angle", resonance)
    @test !occursin("0.5 .* data.sigma_angle", resonance_script)
    @test occursin("1\\ \\mathrm{mm}", resonance)
    @test occursin("chi2/ndf = 0.248081", resonance)
    @test occursin("status = review - inspect diagnostics", resonance)
    @test occursin("SCIENTIFICFITTING_RENDER_DOC_ASSETS", resonance_script)
    @test occursin("SCIENTIFICFITTING_DOC_ASSET_DIR", resonance_script)
    @test occursin("width=panel_width", resonance_script)
    @test !occursin("h_fit = me *", photoelectric)
    @test !occursin("examples/gallery/09_docs_gallery_suite.jl", photoelectric)

    constraints = read(joinpath(GALLERY_SRC, "constraints_profiles.md"), String)
    @test occursin("using CairoMakie", constraints)
    @test occursin("response = [", constraints)
    @test occursin("ScientificFitting.contour(", constraints)
    @test occursin("plot_profile_matrix(profile_overview)", constraints)
    @test occursin("amplitude_profile = amplitude_interval.profile_result", constraints)
    @test !occursin("residual_pattern", constraints)

    gallery_generator = read(
        joinpath(ROOT, "examples", "gallery", "09_docs_gallery_suite.jl"),
        String,
    )
    @test occursin(
        "photoelectric_slope = emission_slope - baseline_slope",
        gallery_generator,
    )
    @test !occursin("h_fit = emission_slope *", gallery_generator)
    @test occursin("SCIENTIFICFITTING_RENDER_DOC_ASSETS", gallery_generator)
    @test occursin("SCIENTIFICFITTING_DOC_ASSET_GROUP", gallery_generator)
    @test occursin("observed_density = counts ./ widths", gallery_generator)
    @test occursin("expected_density = expected ./ widths", gallery_generator)
    @test !occursin("background = result.params[4] .* widths", gallery_generator)
    @test occursin("profile_overview = profile_matrix(", gallery_generator)
    @test occursin(r"plot_profile_matrix\(\s*profile_overview;", gallery_generator)

    multi_dataset = read(joinpath(GALLERY_SRC, "multi_dataset.md"), String)
    @test occursin("y_a = [", multi_dataset)
    @test occursin("nested_pvalue = ccdf(Chisq(1), delta_chi2)", multi_dataset)
    @test occursin("status = review - inspect diagnostics", multi_dataset)
    @test !occursin("pattern_a", multi_dataset)

    multi_dataset_script = read(
        joinpath(ROOT, "examples", "gallery", "10_multi_dataset_calibration.jl"),
        String,
    )
    @test occursin("const MULTI_RENDER_DOC_ASSETS", multi_dataset_script)
    @test occursin("get(ENV, \"SCIENTIFICFITTING_RENDER_DOC_ASSETS\", \"0\") == \"1\"", multi_dataset_script)
    @test occursin("nested_pvalue = ccdf(Chisq(1), delta_chi2)", multi_dataset_script)
    @test !occursin("pattern_a", multi_dataset_script)

    for script_name in (
        "05_constraints_priors_profiles.jl",
        "06_likelihood_workflows.jl",
    )
        script = read(joinpath(ROOT, "examples", "gallery", script_name), String)
        @test occursin("using CairoMakie", script)
    end

    constraints_script = read(
        joinpath(ROOT, "examples", "gallery", "05_constraints_priors_profiles.jl"),
        String,
    )
    @test occursin("ScientificFitting.contour(", constraints_script)
end
