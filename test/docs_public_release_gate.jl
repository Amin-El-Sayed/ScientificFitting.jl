using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")
const DOCS_MAKE = joinpath(ROOT, "docs", "make.jl")

const PUBLIC_DOC_PAGES = [
    "index.md",
    "install.md",
    "quickstart.md",
    "how_jufitter_works.md",
    "gallery.md",
    "gallery/linear_calibration.md",
    "gallery/xy_uncertainties.md",
    "gallery/full_covariance.md",
    "gallery/resonance_decay.md",
    "gallery/photoelectric_threshold.md",
    "gallery/constraints_profiles.md",
    "gallery/poisson_histogram.md",
    "gallery/multi_dataset.md",
    "fitting_for_practitioners.md",
    "plotting_design.md",
    "statistical_foundations.md",
    "api.md",
    "overview.md",
    "backend_design.md",
    "performance.md",
    "maintenance.md",
]

const PUBLIC_TEXT_FILES = vcat(joinpath.(DOCS_SRC, PUBLIC_DOC_PAGES), [joinpath(ROOT, "README.md")])

const FORBIDDEN_PUBLIC_PATTERNS = Pair{String, Regex}[
    "AI/LLM disclosure text" => r"(?i)\b(as an ai|chatgpt|large language model|ai[- ]?generated|ai slop)\b",
    "placeholder marker" => r"(?i)\b(todo|fixme|lorem ipsum|placeholder prose|being rewritten|not all of them are finished|work in progress|coming soon|to be written|to be added)\b",
    "draft/tutorial residue" => r"(?i)\b(draft-only|toy example|left as an exercise|why this example matters|synthetic perfect-data)\b",
    "private local path" => r"(?i)(/Users/|Documents/Projekte|private P1|P1-Praktikum|Praktikum)",
    "private author handle in public prose" => r"(?i)\bAmin_El_Sayed\b",
    "course-internal wording" => r"(?i)\b(course[- ]internal|lab-course-internal|private dataset)\b",
    "stale public API identifier" => r"\b(profile_curve|contour_grid)\b",
    "ungrouped TeX operator subscript" => r"_(\\min|\\max)\b",
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

function documenter_make_text()
    return read(DOCS_MAKE, String)
end

function install_page_text()
    return public_file_text(joinpath(DOCS_SRC, "install.md"))
end

function quickstart_page_text()
    return public_file_text(joinpath(DOCS_SRC, "quickstart.md"))
end

function home_page_text()
    return public_file_text(joinpath(DOCS_SRC, "index.md"))
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
    @testset "Documenter navigation coverage" begin
        @test documenter_navigation_pages() == sort(PUBLIC_DOC_PAGES)
        @test !occursin(r"(?m)^\s{8}\"Engineering Notes\"\s*=>", documenter_make_text())
        @test occursin(r"(?m)^\s{12}\"Technical Notes\"\s*=>", documenter_make_text())
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

    @testset "Python interoperability documentation" begin
        text = install_page_text()
        @test occursin("juliacall", text)
        @test occursin("examples/python/fit_from_python.py", text)
        @test occursin("JUFITTER_RUN_PYTHON_INTEROP=1", text)
        @test occursin("experimental or deferred", text)
    end

    @testset "First-user path is executable and honest" begin
        home = home_page_text()
        quickstart = quickstart_page_text()
        install = install_page_text()

        @test occursin("```@raw html\n<section class=\"jufitter-hero\">", home)
        @test occursin("data-jufitter-plot-group=\"home-first-fit\"", home)
        @test !occursin("collect(range", home)
        @test occursin("using CairoMakie", quickstart)
        @test !occursin("collect(range", quickstart)
        @test occursin("report=:both", quickstart)
        @test !occursin("println(report_text", quickstart)
        @test !occursin("Pkg.test()", install)
        @test !occursin("canonical=", documenter_make_text())
    end
end
