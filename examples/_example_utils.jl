const EXAMPLE_OUTPUT_DIR = joinpath(@__DIR__, "output")
mkpath(EXAMPLE_OUTPUT_DIR)

example_output(filename::AbstractString) = joinpath(EXAMPLE_OUTPUT_DIR, filename)

function print_result_summary(title::AbstractString, result; parameter_names=nothing)
    println()
    println(title)
    println(repeat("-", length(title)))
    println(report_text(result; parameter_names=parameter_names))
end
