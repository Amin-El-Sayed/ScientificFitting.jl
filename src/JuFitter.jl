module JuFitter

using CairoMakie
using ADTypes
using Distributions
using DifferentiationInterface
using ForwardDiff
using LaTeXStrings
using LinearAlgebra
using LsqFit
using Optimization
using OptimizationOptimJL
using QuadGK
using SpecialFunctions
using SparseArrays
using Statistics

include("types.jl")
include("parameters.jl")
include("diagnostics.jl")
include("weights.jl")
include("costs.jl")
include("fit.jl")
include("likelihood_fits.jl")
include("profile.jl")
include("plotting.jl")
include("report.jl")

export ConstraintSpec
export ParameterPrior
export FixedParameter
export ParameterConstraint
export ErrorComponent
export FitOptions
export FitProblem
export FitResult
export LikelihoodFitProblem
export LikelihoodFitResult
export FitStatistics
export FitDiagnostics
export DiagnosticFinding
export DiagnosticReport
export FitReport
export ParameterEstimate
export ProfileResult
export ContourResult
export ProfileInterval
export fit
export fit_model
export fit_custom
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
export fitplot
export plot_fit
export plot_residuals
export plot_diagnostics
export plot_profile
export plot_contour
export diagnose
export diagnose_text
export report_text

end # module JuFitter
