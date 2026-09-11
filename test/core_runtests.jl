using Test

@testset "ScientificFitting core test suite" begin
    for file in (
        "regression/current_api.jl",
        "statistics/covariance_semantics_reference.jl",
        "statistics/diagnostics_reference.jl",
        "statistics/linear_gaussian_reference.jl",
        "statistics/structured_whitening_reference.jl",
        "statistics/likelihood_reference.jl",
        "statistics/vectorized_density_reference.jl",
        "statistics/observation_likelihood_reference.jl",
        "statistics/distribution_objects.jl",
        "statistics/profile_contour_reference.jl",
        "numerics/inplace_model_reference.jl",
        "numerics/finite_derivatives_reference.jl",
        "numerics/nonsmooth_likelihood_reference.jl",
        "numerics/solver_status_reference.jl",
        "numerics/solver_contract.jl",
        "numerics/torture_inputs.jl",
    )
        # Show progress before compilation; Julia's test summary appears only at the end.
        println("Core test: ", file)
        flush(stdout)
        @time include(file)
    end
end
