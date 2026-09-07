"""Native Python fitting with a shared Julia numerical core and optional Matplotlib.

Development API: models receive `(x, **parameters)` as NumPy arrays and floats.
Matplotlib is imported only when plotting; Julia/Makie figures are never exposed.
"""

from ._core import (
    Result, fit_custom, fit_extended_unbinned_model, fit_histogram_density,
    fit_histogram_model, fit_indexed_model, fit_likelihood_model, fit_model,
    fit_multi_model, fit_poisson_model, fit_unbinned_model,
)
from ._inputs import ErrorComponent, WhiteningOperator
from ._results import (
    ContourResult, DiagnosticFinding, DiagnosticReport, FitReport, ParameterEstimate,
    ProfileInterval, ProfileMatrixPanelTriage, ProfileMatrixResult, ProfileResult,
)
# Renderers import Matplotlib inside calls, preserving real signatures/docstrings
# on the public API without making plotting dependencies mandatory at import.
from ._plotting import add_report, plot_fit, plot_style
from ._diagnostic_plots import (
    plot_contour, plot_diagnostics, plot_profile, plot_profile_matrix, plot_residuals,
)


__all__ = [
    "Result", "ErrorComponent", "WhiteningOperator", "fit_model", "fit_custom", "fit_likelihood_model", "fit_poisson_model",
    "fit_histogram_model", "fit_histogram_density", "fit_unbinned_model",
    "fit_extended_unbinned_model", "fit_indexed_model", "fit_multi_model", "plot_fit",
    "DiagnosticFinding", "DiagnosticReport", "FitReport", "ParameterEstimate",
    "ProfileResult", "ProfileInterval", "ContourResult", "ProfileMatrixResult",
    "ProfileMatrixPanelTriage",
    "plot_style", "add_report", "plot_profile", "plot_contour", "plot_profile_matrix",
    "plot_residuals", "plot_diagnostics",
]
