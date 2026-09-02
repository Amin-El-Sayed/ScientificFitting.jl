using Test

const ROOT = abspath(joinpath(@__DIR__, ".."))
const DOCS_SRC = joinpath(ROOT, "docs", "src")

const OUTPUT_EXPECTATIONS = [
    ("quickstart.md", "quickstart", 1),
    (joinpath("gallery", "linear_calibration.md"), "linear_calibration", 1),
    (joinpath("gallery", "xy_uncertainties.md"), "xy_uncertainties", 1),
    (joinpath("gallery", "full_covariance.md"), "full_covariance", 1),
    (joinpath("gallery", "resonance_decay.md"), "resonance_constant", 1),
    (joinpath("gallery", "resonance_decay.md"), "resonance_drift", 2),
    (joinpath("gallery", "photoelectric_threshold.md"), "photoelectric_threshold", 1),
    (joinpath("gallery", "constraints_profiles.md"), "constraints_profiles", 1),
    (joinpath("gallery", "poisson_histogram.md"), "poisson_decay", 1),
    (joinpath("gallery", "poisson_histogram.md"), "histogram_likelihood", 2),
    (joinpath("gallery", "multi_dataset.md"), "multi_dataset", 1),
]

const EXECUTABLE_PAGES = unique(first.(OUTPUT_EXPECTATIONS))

function normalize_output(text::AbstractString)
    cleaned = replace(text, "\r\n" => "\n", "\r" => "\n")
    lines = String.(split(cleaned, '\n'))
    while !isempty(lines) && isempty(strip(first(lines)))
        popfirst!(lines)
    end
    while !isempty(lines) && isempty(strip(last(lines)))
        pop!(lines)
    end
    return lines
end

function documented_output_blocks(relative_page::AbstractString)
    page = read(joinpath(DOCS_SRC, relative_page), String)
    pattern = r"<div class=\"scientificfitting-cell-output\">.*?<pre>(.*?)</pre>.*?</div>"s
    return [match.captures[1] for match in eachmatch(pattern, page)]
end

function documented_julia_blocks(relative_page::AbstractString)
    page = read(joinpath(DOCS_SRC, relative_page), String)
    return [match.captures[1] for match in eachmatch(r"^```julia\n(.*?)\n```"ms, page)]
end

function markdown_pages()
    pages = String[]
    for (directory, _, filenames) in walkdir(DOCS_SRC)
        for filename in filenames
            endswith(filename, ".md") || continue
            push!(pages, relpath(joinpath(directory, filename), DOCS_SRC))
        end
    end
    return sort(pages)
end

function marker_outputs(text::AbstractString)
    pattern = r"=== SCIENTIFICFITTING_DOC_OUTPUT_BEGIN ([^\n]+?) ===\n(.*?)\n=== SCIENTIFICFITTING_DOC_OUTPUT_END \1 ==="s
    outputs = Dict{String,String}()
    for match in eachmatch(pattern, text)
        id, body = match.captures
        outputs[id] = body
    end
    return outputs
end

function run_snapshot_script(script::AbstractString)
    cmd = `$(Base.julia_cmd()) --project=docs --startup-file=no $script`
    return read(
        setenv(
            cmd,
            "SCIENTIFICFITTING_DOC_OUTPUT_SNAPSHOTS" => "1",
            "SCIENTIFICFITTING_DOC_SNAPSHOT_ONLY" => "1",
        ),
        String,
    )
end

function documented_page_outputs(text::AbstractString)
    pattern = r"=== SCIENTIFICFITTING_DOC_PAGE_BEGIN ([^\n]+?) ===\n(.*?)\n=== SCIENTIFICFITTING_DOC_PAGE_END \1 ==="s
    outputs = Dict{String,String}()
    for match in eachmatch(pattern, text)
        page, body = match.captures
        outputs[page] = body
    end
    return outputs
end

function run_documented_pages(relative_pages)
    return mktempdir() do root
        entries = NamedTuple[]
        source_data = joinpath(ROOT, "examples", "data", "damped_oscillator")

        for (index, relative_page) in enumerate(relative_pages)
            blocks = documented_julia_blocks(relative_page)
            isempty(blocks) && error("$relative_page has no executable Julia blocks")

            directory = joinpath(root, "page_$index")
            mkpath(directory)

            # The oscillator page intentionally uses a public repository-relative
            # data path. Recreate it so every page runs from a clean directory.
            target_data = joinpath(directory, "examples", "data", "damped_oscillator")
            mkpath(target_data)
            for filename in readdir(source_data)
                cp(joinpath(source_data, filename), joinpath(target_data, filename))
            end

            script = joinpath(directory, "documented_example.jl")
            write(script, join(blocks, "\n\n"))
            push!(entries, (; index, relative_page, directory, script))
        end

        # One Julia process loads ScientificFitting and Makie once. A fresh module and
        # working directory still isolate names and generated files per page.
        driver = joinpath(root, "run_documented_pages.jl")
        open(driver, "w") do io
            println(io, "function run_page(index, page, directory, script)")
            println(io, "    println(\"=== SCIENTIFICFITTING_DOC_PAGE_BEGIN \", page, \" ===\")")
            println(io, "    page_module = Module(Symbol(\"ScientificFittingDocsPage\", index))")
            println(io, "    cd(directory) do")
            println(io, "        Base.include(page_module, script)")
            println(io, "    end")
            println(io, "    println(\"=== SCIENTIFICFITTING_DOC_PAGE_END \", page, \" ===\")")
            println(io, "end")
            for entry in entries
                println(
                    io,
                    "run_page(", entry.index, ", ", repr(entry.relative_page), ", ",
                    repr(entry.directory), ", ", repr(entry.script), ")",
                )
            end
        end

        command = `$(Base.julia_cmd()) --project=$(joinpath(ROOT, "docs")) --startup-file=no $driver`
        stdout = IOBuffer()
        stderr = IOBuffer()
        try
            run(pipeline(command; stdout, stderr))
        catch
            message = String(take!(stderr))
            isempty(strip(message)) || @info "Documented examples emitted stderr" message
            rethrow(error)
        end
        warnings = String(take!(stderr))
        return documented_page_outputs(String(take!(stdout))), warnings
    end
end

@testset "Documentation output snapshots" begin
    @testset "Every Julia fence parses" begin
        for relative_page in markdown_pages()
            for (index, block) in enumerate(documented_julia_blocks(relative_page))
                @testset "$relative_page block $index" begin
                    @test Meta.parseall("begin\n$block\nend"; filename=relative_page) isa Expr
                end
            end
        end
    end

    page_outputs, page_stderr = run_documented_pages(EXECUTABLE_PAGES)
    gallery_output = run_snapshot_script(joinpath(ROOT, "examples", "gallery", "09_docs_gallery_suite.jl"))
    resonance_output = run_snapshot_script(joinpath(ROOT, "examples", "gallery", "08_damped_oscillator_decay.jl"))
    snapshots = merge(marker_outputs(gallery_output), marker_outputs(resonance_output))
    snapshots["quickstart"] = page_outputs["quickstart.md"]

    @test isempty(strip(page_stderr))

    for (_, id, _) in OUTPUT_EXPECTATIONS
        @test haskey(snapshots, id)
    end

    @testset "Every result cell has one real generator snapshot" begin
        expected_ids = Set(snapshot_id for (_, snapshot_id, _) in OUTPUT_EXPECTATIONS)
        @test Set(keys(snapshots)) == expected_ids

        for relative_page in EXECUTABLE_PAGES
            expected_count = count(expectation -> expectation[1] == relative_page, OUTPUT_EXPECTATIONS)
            @test length(documented_output_blocks(relative_page)) == expected_count
        end
    end

    for (relative_page, snapshot_id, block_index) in OUTPUT_EXPECTATIONS
        blocks = documented_output_blocks(relative_page)
        @test block_index <= length(blocks)

        documented = normalize_output(blocks[block_index])
        actual = normalize_output(snapshots[snapshot_id])
        @test documented == actual
        if documented != actual
            @info "Documented output differs from its real generator snapshot" relative_page snapshot_id
        end
    end

    @testset "Complete documentation workflows execute exactly as shown" begin
        for relative_page in EXECUTABLE_PAGES
            @testset "$relative_page" begin
                documented_blocks = documented_output_blocks(relative_page)
                documented_output = reduce(
                    vcat,
                    normalize_output.(documented_blocks);
                    init=String[],
                )

                @test haskey(page_outputs, relative_page)
                @test normalize_output(page_outputs[relative_page]) == documented_output
            end
        end
    end
end
