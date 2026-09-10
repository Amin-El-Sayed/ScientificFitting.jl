module ScientificFitting

using ADTypes
import StatsAPI: fit, predict
using Distributions
using DifferentiationInterface
using FiniteDiff
using ForwardDiff
using LaTeXStrings
using LinearAlgebra
using LsqFit
using Optimization
import OptimizationNLopt
using OptimizationOptimJL
using PrecompileTools: @setup_workload, @compile_workload
using QuadGK
using SpecialFunctions
using SparseArrays
using Statistics

include("formatting.jl")
include("solvers.jl")
include("types.jl")
include("derivatives.jl")
include("parameters.jl")
include("diagnostics.jl")
include("weights.jl")
include("costs.jl")
include("fit.jl")
include("likelihood_fits.jl")
include("profile.jl")
include("prediction.jl")
include("plotting_api.jl")
include("report.jl")
include("precompile.jl")

export ConstraintSpec
export ParameterPrior
export FixedParameter
export ParameterConstraint
export ErrorComponent
export WhiteningOperator
export FitOptions
export AbstractFitSolver, OptimizationSolver, NativeMinuitSolver, FitSolverResult
export solver_capabilities, solve_fit
export FitProblem
export FitResult
export LikelihoodFitProblem
export LikelihoodFitResult
export FitStatistics
export FitDiagnostics
export DiagnosticFinding
export DiagnosticReport
export DiagnosticDashboard
export FitReport
export ParameterEstimate
export ProfileResult
export ContourResult
export ProfileInterval
export ProfileMatrixResult
export ProfileMatrixPanelTriage
export fit
export predict
export fit_model
export fit_custom
export fit_likelihood_model
export fit_poisson_model
export fit_histogram_model
export fit_histogram_density
export fit_unbinned_model
export fit_extended_unbinned_model
export fit_indexed_model
export fit_multi_model
export fit_report
export profile
export profile_interval
export contour
export profile_matrix
export profile_matrix_triage
export fitplot
export plot_fit
export fit_axis
export add_curve!
export add_points!
export add_vline!
export add_hline!
export add_vband!
export add_hband!
export plot_theme
export plot_palette
export plot_info_panel!
export resize_plot_to_layout!
export plot_residuals
export plot_diagnostics
export plot_profile
export plot_contour
export plot_profile_matrix
export diagnose
export diagnose_text
export diagnostic_dashboard
export diagnostic_dashboard_text
export report_text

end # module ScientificFitting
