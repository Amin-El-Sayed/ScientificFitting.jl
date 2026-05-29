const EXAMPLE_OUTPUT_DIR = joinpath(@__DIR__, "output")
mkpath(EXAMPLE_OUTPUT_DIR)

example_output(filename::AbstractString) = joinpath(EXAMPLE_OUTPUT_DIR, filename)
