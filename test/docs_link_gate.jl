using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")
const DOCS_MAKE = joinpath(ROOT, "docs", "make.jl")

function docs_markdown_files()
    files = String[]
    for (dir, _, names) in walkdir(DOCS_SRC)
        for name in names
            endswith(name, ".md") || continue
            push!(files, joinpath(dir, name))
        end
    end
    return sort(files)
end

function is_external_link(target::AbstractString)
    return startswith(target, "http://") ||
           startswith(target, "https://") ||
           startswith(target, "mailto:") ||
           startswith(target, "#") ||
           startswith(target, "@ref") ||
           isempty(strip(target))
end

function strip_fragment(target::AbstractString)
    return first(split(target, '#'; limit=2))
end

function resolve_docs_target(page::AbstractString, raw_target::AbstractString)
    target = strip_fragment(strip(raw_target))
    isempty(target) && return nothing

    if endswith(target, ".html")
        target = target[1:end - length(".html")] * ".md"
    end

    return normpath(joinpath(dirname(page), target))
end

function html_targets(text::AbstractString)
    targets = Pair{String,String}[]
    for match in eachmatch(r"\b(href|src)=\"([^\"]+)\"", text)
        push!(targets, match.captures[1] => match.captures[2])
    end
    return targets
end

function markdown_targets(text::AbstractString)
    targets = String[]
    for match in eachmatch(r"!?\[[^\]]*\]\(([^\)\s]+)(?:\s+\"[^\"]*\")?\)", text)
        push!(targets, match.captures[1])
    end
    return targets
end

@testset "Documentation local links" begin
    @testset "Portable static routes" begin
        make_source = read(DOCS_MAKE, String)
        @test occursin(r"prettyurls\s*=\s*false", make_source)
        @test !occursin(r"prettyurls\s*=\s*get\(ENV", make_source)
    end

    for page in docs_markdown_files()
        text = read(page, String)

        @testset "$(relpath(page, ROOT))" begin
            for target in markdown_targets(text)
                is_external_link(target) && continue
                resolved = resolve_docs_target(page, target)
                @test resolved !== nothing
                @test isfile(resolved)
            end

            for (attr, target) in html_targets(text)
                is_external_link(target) && continue
                resolved = resolve_docs_target(page, target)
                @test resolved !== nothing
                @test isfile(resolved)

                if attr == "href" && endswith(strip_fragment(target), ".html")
                    @test endswith(resolved, ".md")
                end
            end
        end
    end
end
