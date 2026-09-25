using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")
const DOCS_MAKE = joinpath(ROOT, "docs", "make.jl")
const PUBLIC_DOC_PAGES = [
    "index.md",
    "install.md",
    "python.md",
    "quickstart.md",
    "how_scientificfitting_works.md",
    "gallery.md",
    "gallery/linear_calibration.md",
    "gallery/xy_uncertainties.md",
    "gallery/full_covariance.md",
    "gallery/resonance_decay.md",
    "gallery/photoelectric_threshold.md",
    "gallery/constraints_profiles.md",
    "gallery/poisson_histogram.md",
    "gallery/multi_dataset.md",
    "gallery/lhcb_mass_spectrum.md",
    "fitting_for_practitioners.md",
    "interfaces.md",
    "plotting_design.md",
    "statistical_foundations.md",
    "gaussian_models.md",
    "parameter_inference.md",
    "profiles_contours.md",
    "likelihood_models.md",
    "api.md",
    "api_fitting.md",
    "api_results.md",
    "api_plotting.md",
    "api_plotting_diagnostics.md",
    "citation.md",
    "backend_design.md",
    "performance.md",
]

const PUBLIC_TEXT_FILES = vcat(joinpath.(DOCS_SRC, PUBLIC_DOC_PAGES), [joinpath(ROOT, "README.md")])

# Content classes that must never reach the published site: leaked private
# context, unfinished-work markers, and tooling disclosure. These forbid kinds
# of content, not specific documentation wording.
const FORBIDDEN_PUBLIC_PATTERNS = Pair{String, Regex}[
    "AI/LLM disclosure text" => r"(?i)\b(as an ai|chatgpt|large language model|ai[- ]?generated|ai slop)\b",
    "placeholder marker" => r"(?i)\b(todo|fixme|lorem ipsum|placeholder prose|being rewritten|work in progress|coming soon|to be written|to be added)\b",
    "private local path" => r"(?i)((?<![\w/])/Users/|file:///Users/|Documents/Projekte|private P1|P1-Praktikum|Praktikum)",
    "private author handle in public prose" => r"(?i)\bAmin_El_Sayed\b",
    "course-internal wording" => r"(?i)\b(course[- ]internal|lab-course-internal|private dataset)\b",
]

function public_file_text(path)
    isfile(path) || error("public documentation file missing: $(relpath(path, ROOT))")
    return read(path, String)
end

function documenter_navigation_pages()
    text = read(DOCS_MAKE, String)
    pages = String[]
    for match in eachmatch(r"\"([^\"]+\.md)\"", text)
        push!(pages, match.captures[1])
    end
    return sort(unique(pages))
end

function docs_source_markdown_pages()
    pages = String[]
    for (directory, _, filenames) in walkdir(DOCS_SRC)
        for filename in filenames
            endswith(filename, ".md") || continue
            push!(pages, relpath(joinpath(directory, filename), DOCS_SRC))
        end
    end
    return sort(pages)
end

function markdown_outside_docs()
    # Audit publishable source, not ignored pytest/virtual-environment caches.
    paths = split(read(`git -C $ROOT ls-files --cached --others --exclude-standard -z`, String), '\0')
    return sort(unique(filter(path -> endswith(path, ".md") && !startswith(path, "docs/"), paths)))
end

function markdown_image_alt_texts(text::AbstractString)
    return [match.captures[1] for match in eachmatch(r"!\[([^\]]*)\]\([^)]+\)", text)]
end

function html_image_tags(text::AbstractString)
    return [match.captures[1] for match in eachmatch(r"<img\b([^>]*)>", text)]
end

function html_image_alt_text(tag::AbstractString)
    match_result = match(r"\balt=\"([^\"]*)\"", tag)
    return isnothing(match_result) ? nothing : match_result.captures[1]
end

@testset "Public documentation release hygiene" begin
    @testset "Private-path pattern matches machine-local paths only" begin
        pattern = Dict(FORBIDDEN_PUBLIC_PATTERNS)["private local path"]
        @test occursin(pattern, "`/Users/example/project`")
        @test occursin(pattern, "file:///Users/example/project")
        @test !occursin(pattern, "https://matplotlib.org/stable/users/explain/axes/constrainedlayout_guide.html")
    end

    @testset "Documenter navigation covers every page" begin
        @test documenter_navigation_pages() == sort(setdiff(PUBLIC_DOC_PAGES, ["index.md"]))
        @test docs_source_markdown_pages() == sort(PUBLIC_DOC_PAGES)
    end

    @testset "Repository root stays package-facing" begin
        root_markdown = sort(filter(name -> endswith(name, ".md"), readdir(ROOT)))
        # Both distributable packages need a public README, not internal notes.
    end

    @testset "Configured public files exist" begin
        for path in PUBLIC_TEXT_FILES
            @test isfile(path)
        end
    end

    for path in PUBLIC_TEXT_FILES
        text = public_file_text(path)
        rel = relpath(path, ROOT)
        @testset "$rel" begin
            @test !isempty(strip(text))
            for (label, pattern) in FORBIDDEN_PUBLIC_PATTERNS
                @testset "$label" begin
                    @test !occursin(pattern, text)
                end
            end
            @testset "image alt text" begin
                for alt in markdown_image_alt_texts(text)
                    @test !isempty(strip(alt))
                end
                for tag in html_image_tags(text)
                    alt = html_image_alt_text(tag)
                    @test !isnothing(alt)
                    @test !isempty(strip(alt))
                end
            end
        end
    end
end
