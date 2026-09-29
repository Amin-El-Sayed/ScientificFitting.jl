using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")

# Fenced code that Documenter either executes (@example/@setup/@repl) or
# presents as Julia source. julia-repl transcripts are excluded: they contain
# prompt characters and are covered by the output snapshot gate.
const CODE_FENCE = r"```(?:julia|@example[^\n]*|@setup[^\n]*|@repl[^\n]*)\n(.*?)\n```"s

function markdown_pages()
    pages = String[]
    for (dir, _, names) in walkdir(DOCS_SRC)
        for name in names
            endswith(name, ".md") && push!(pages, joinpath(dir, name))
        end
    end
    return sort(pages)
end

function code_blocks(text::AbstractString)
    return [match.captures[1] for match in eachmatch(CODE_FENCE, text)]
end

function has_parse_error(ex)
    ex isa Expr || return false
    ex.head in (:error, :incomplete) && return true
    return any(has_parse_error, ex.args)
end

function parses_cleanly(code::AbstractString)
    parsed = try
        Meta.parseall(String(code))
    catch
        return false
    end
    return !has_parse_error(parsed)
end

function image_sources(text::AbstractString)
    sources = String[]
    for match in eachmatch(r"<img[^>]+src=\"([^\"]+)\"", text)
        push!(sources, match.captures[1])
    end
    for match in eachmatch(r"!\[[^\]]*\]\(([^)\s]+)", text)
        push!(sources, match.captures[1])
    end
    return sources
end

function is_external(src::AbstractString)
    return startswith(src, "http://") || startswith(src, "https://") || startswith(src, "data:")
end

@testset "Documentation code parses and assets exist" begin
    for page in markdown_pages()
        text = read(page, String)
        rel = relpath(page, ROOT)
        @testset "$rel" begin
            for (i, code) in enumerate(code_blocks(text))
                @testset "code block $i" begin
                    @test parses_cleanly(code)
                end
            end
            for src in image_sources(text)
                is_external(src) && continue
                # Figures rendered by the page's own @setup blocks only exist
                # after the Documenter build; the rendered link gate checks them.
                generated = occursin("```@setup", text) && endswith(src, ".svg")
                generated && continue
                @test isfile(normpath(joinpath(dirname(page), src)))
            end
        end
    end
end
