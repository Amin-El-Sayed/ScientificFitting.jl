# Time-to-first-fit probe behind the table in docs/src/validation.md.
# Every cell is timed in its own fresh Julia process, so each number includes
# exactly the compilation a user pays in a new session; rerun this script to
# reproduce the table on your hardware.
using Printf

const PROJECT = abspath(joinpath(@__DIR__, ".."))

function timed_cell(code::AbstractString)
    cmd = `$(Base.julia_cmd()) --project=$PROJECT --startup-file=no -e $code`
    return parse(Float64, strip(read(cmd, String)))
end

first_fit(mode) = """
using ScientificFitting
x = collect(1.0:30.0); y = 2 .* x .+ 1 .+ 0.01 .* sin.(x)
mA(x, p) = @. p[1] * x + p[2]
t = @elapsed fit_model(mA, x, y; p0=[1.0, 0.0], derivatives=$mode)
println(t)
"""

# One model type is already fitted; the cell times the next, distinct type.
second_model(mode) = """
using ScientificFitting
x = collect(1.0:30.0); y = 2 .* x .+ 1 .+ 0.01 .* sin.(x)
mA(x, p) = @. p[1] * x + p[2]
fit_model(mA, x, y; p0=[1.0, 0.0], derivatives=$mode)
mB(x, p) = @. p[1] * exp(-x / p[2])
yb = 5 .* exp.(-x ./ 7.0)
t = @elapsed fit_model(mB, x, yb; p0=[4.0, 5.0], derivatives=$mode)
println(t)
"""

# A Gaussian fit and one Poisson rate model are warm; the cell times a second,
# distinct rate model through the Poisson pipeline.
poisson_extra(mode) = """
using ScientificFitting
x = collect(1.0:30.0); y = 2 .* x .+ 1 .+ 0.01 .* sin.(x)
mA(x, p) = @. p[1] * x + p[2]
fit_model(mA, x, y; p0=[1.0, 0.0], derivatives=$mode)
cx = collect(1.0:8.0); counts = round.(Int, 40 .* exp.(-cx ./ 4))
r1(x, p) = @. p[1] * exp(-x / p[2])
fit_poisson_model(r1, cx, counts; p0=[45.0, 3.5], derivatives=$mode)
r2(x, p) = @. p[1] * exp(-x / p[2]) + 0 * p[1]
t = @elapsed fit_poisson_model(r2, cx, counts; p0=[45.0, 3.5], derivatives=$mode)
println(t)
"""

for (label, cell) in (
    ("first fit after `using`", first_fit),
    ("each additional model type", second_model),
    ("Poisson fit, additional rate model", poisson_extra),
)
    t_auto = timed_cell(cell(":auto"))
    t_finite = timed_cell(cell(":finite"))
    @printf "%-36s :auto %6.2f s   :finite %6.2f s\n" label t_auto t_finite
end
