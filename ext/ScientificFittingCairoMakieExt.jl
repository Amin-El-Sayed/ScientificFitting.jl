module ScientificFittingCairoMakieExt

using CairoMakie
using LaTeXStrings
using LinearAlgebra
using Statistics

using ScientificFitting

import ScientificFitting:
    add_curve!,
    add_hband!,
    add_hline!,
    add_points!,
    add_vband!,
    add_vline!,
    contour,
    diagnostic_dashboard_text,
    fit_axis,
    fit_model,
    fitplot,
    plot_contour,
    plot_diagnostics,
    plot_fit,
    plot_info_panel!,
    plot_palette,
    plot_profile,
    plot_profile_matrix,
    plot_residuals,
    plot_theme,
    profile,
    report_text,
    _fmt_value,
    _model_dydx,
    _model_values,
    _parameter_jacobian,
    _strip_math_delims,
    _weighted_data_residual,
    _xerror_for_plot,
    _yerror_for_plot

include("../src/plotting.jl")

end
