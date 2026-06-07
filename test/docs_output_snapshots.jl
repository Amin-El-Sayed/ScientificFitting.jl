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

function normalize_output(text::AbstractString)
    cleaned = replace(text, "\r\n" => "\n", "\r" => "\n")
    lines = split(cleaned, '\n')
    return strip.(lines)
end

function nonempty_lines(text::AbstractString)
    return filter(!isempty, normalize_output(text))
end

function documented_output_blocks(relative_page::AbstractString)
    page = read(joinpath(DOCS_SRC, relative_page), String)
    pattern = r"<div class=\"jufitter-cell-output\">.*?<pre>(.*?)</pre>.*?</div>"s
    return [match.captures[1] for match in eachmatch(pattern, page)]
end

function marker_outputs(text::AbstractString)
    pattern = r"=== JUFITTER_DOC_OUTPUT_BEGIN ([^\n]+?) ===\n(.*?)\n=== JUFITTER_DOC_OUTPUT_END \1 ==="s
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
            "JUFITTER_DOC_OUTPUT_SNAPSHOTS" => "1",
            "JUFITTER_DOC_SNAPSHOT_ONLY" => "1",
        ),
        String,
    )
end

function is_ordered_line_subset(expected::Vector{SubString{String}}, actual::Vector{SubString{String}})
    actual_index = 1
    for line in expected
        found = false
        while actual_index <= length(actual)
            if actual[actual_index] == line
                found = true
                actual_index += 1
                break
            end
            actual_index += 1
        end
        found || return false, line
    end
    return true, ""
end

@testset "Documentation output snapshots" begin
    gallery_output = run_snapshot_script(joinpath(ROOT, "examples", "gallery", "09_docs_gallery_suite.jl"))
    resonance_output = run_snapshot_script(joinpath(ROOT, "examples", "gallery", "08_damped_oscillator_decay.jl"))
    snapshots = merge(marker_outputs(gallery_output), marker_outputs(resonance_output))

    for (_, id, _) in OUTPUT_EXPECTATIONS
        @test haskey(snapshots, id)
    end

    for (relative_page, snapshot_id, block_index) in OUTPUT_EXPECTATIONS
        blocks = documented_output_blocks(relative_page)
        @test block_index <= length(blocks)

        documented = nonempty_lines(blocks[block_index])
        actual = nonempty_lines(snapshots[snapshot_id])
        ok, missing = is_ordered_line_subset(documented, actual)
        @test ok
        if !ok
            @info "Documented output line is not present in real snapshot" relative_page snapshot_id missing
        end
    end
end
