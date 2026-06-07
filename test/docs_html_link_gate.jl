using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_BUILD = joinpath(ROOT, "docs", "build")

function html_files()
    files = String[]
    isdir(DOCS_BUILD) || return files
    for (dir, _, names) in walkdir(DOCS_BUILD)
        for name in names
            endswith(name, ".html") || continue
            push!(files, joinpath(dir, name))
        end
    end
    return sort(files)
end

function strip_fragment(target::AbstractString)
    return first(split(target, '#'; limit=2))
end

function is_ignored_rendered_target(target::AbstractString)
    stripped = strip(target)
    cleaned = strip_fragment(stripped)
    return isempty(stripped) ||
           startswith(stripped, "#") ||
           startswith(stripped, "http://") ||
           startswith(stripped, "https://") ||
           startswith(stripped, "mailto:") ||
           startswith(stripped, "javascript:") ||
           basename(cleaned) == "versions.js"
end

function rendered_targets(text::AbstractString)
    targets = Pair{String,String}[]
    for match in eachmatch(r"\b(href|src)=\"([^\"]+)\"", text)
        push!(targets, match.captures[1] => match.captures[2])
    end
    return targets
end

function resolve_rendered_target(page::AbstractString, target::AbstractString)
    cleaned = strip_fragment(strip(target))
    isempty(cleaned) && return nothing
    startswith(cleaned, "/") && return normpath(joinpath(DOCS_BUILD, cleaned[2:end]))
    return normpath(joinpath(dirname(page), cleaned))
end

@testset "Rendered documentation links" begin
    @test isdir(DOCS_BUILD)
    pages = html_files()
    @test !isempty(pages)

    for page in pages
        text = read(page, String)
        @testset "$(relpath(page, DOCS_BUILD))" begin
            for (attr, target) in rendered_targets(text)
                is_ignored_rendered_target(target) && continue
                resolved = resolve_rendered_target(page, target)
                @test resolved !== nothing
                @test isfile(resolved) || isdir(resolved)

                if attr == "href" && !isnothing(resolved) && isdir(resolved)
                    @test isfile(joinpath(resolved, "index.html"))
                end
            end
        end
    end
end
