"""Native Python fitting with a shared Julia numerical core and optional Matplotlib.

Development API: models receive `(x, *parameters)` as NumPy arrays and floats.
Matplotlib is imported only by `plot_fit`; Julia/Makie figures are never exposed.
"""

from ._core import (
    Result, fit_custom, fit_extended_unbinned_model, fit_histogram_density,
    fit_histogram_model, fit_indexed_model, fit_likelihood_model, fit_model,
    fit_multi_model, fit_poisson_model, fit_unbinned_model,
)


def plot_fit(result, **kwargs):
    """Draw an existing Gaussian fit; return a native Matplotlib `(Figure, Axes)`.

    See `scientificfitting._plotting.plot_fit` for display options. Importing or
    fitting with this package does not import a plotting backend.
    """
    from ._plotting import plot_fit as draw

    return draw(result, **kwargs)


__all__ = [
    "Result", "fit_model", "fit_custom", "fit_likelihood_model", "fit_poisson_model",
    "fit_histogram_model", "fit_histogram_density", "fit_unbinned_model",
    "fit_extended_unbinned_model", "fit_indexed_model", "fit_multi_model", "plot_fit",
]
