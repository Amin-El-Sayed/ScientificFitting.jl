"""Numeric ownership and uncertainty specifications at the Python boundary."""

from dataclasses import dataclass

import numpy as np


def _real_array(values):
    """Convert real inputs without silently dropping imaginary components."""
    if np.iscomplexobj(values):
        raise ValueError("fitting inputs must be real; complex values cannot be discarded")
    return np.asarray(values, dtype=np.float64)


def _array(values, *, ndim=1):
    array = _real_array(values).copy()
    if array.ndim != ndim:
        raise ValueError(f"expected a {ndim}-dimensional numeric array")
    return array


def _readonly(values):
    """Borrow a numeric input for one callback, without changing its owner's flags."""
    view = np.asarray(values).view()
    view.flags.writeable = False
    return view


def _snapshot(values):
    """Own result arrays independently of Julia; never change backend storage."""
    array = np.array(values, dtype=np.float64, copy=True)
    array.flags.writeable = False
    return array


def _covariance(values):
    """Copy SciPy sparse input as canonical CSC buffers, never as a dense matrix."""
    if hasattr(values, "tocsc"):
        from scipy.sparse import issparse

        if issparse(values):
            if np.issubdtype(values.dtype, np.complexfloating):
                raise ValueError("covariance must be real")
            csc = values.tocsc(copy=True).astype(np.float64, copy=False)
            # Validate before handing index buffers to Julia's sparse constructor.
            csc.check_format(full_check=True)
            csc.sum_duplicates()
            csc.sort_indices()
            return {"shape": csc.shape, "data": csc.data,
                    "indices": csc.indices, "indptr": csc.indptr}
    return _array(values, ndim=2)


def _error_values(values):
    if hasattr(values, "tocsc"):
        return _covariance(values)
    array = _real_array(values)
    if array.ndim == 0:
        return float(array)
    return _array(array, ndim=1) if array.ndim == 1 else _covariance(array)


@dataclass(frozen=True)
class ErrorComponent:
    """Named Gaussian uncertainty source, combined in covariance space.

    `target` is "x" or "y". `mode` is "absolute" (standard deviations),
    "relative" (fractions of measured values), "model_relative" (fractions
    of fitted y), or "covariance" (a covariance matrix or vector of standard
    deviations). `values` may be scalar, pointwise, dense, or SciPy sparse as
    appropriate for the mode. `active=False` keeps a source out of the fit.
    Inputs are copied when fitting; use dataclasses.replace to toggle sources
    for a new fit. Parameter-dependent sources retain their log determinant.
    """

    name: str
    target: str
    mode: str
    values: object
    active: bool = True

    def _payload(self):
        return (self.name, self.target, self.mode, _error_values(self.values), self.active)


@dataclass(frozen=True)
class WhiteningOperator:
    """Complete static observation covariance through a matrix-free operator.

    `apply(residual)` returns W @ residual, with W.T @ W = inv(C).
    With `inplace=True`, use `apply(out, residual)` and fill every output entry,
    returning None. Inputs are read-only NumPy views; output views write into
    the Julia buffer and must not be retained. The operator must be linear and
    remain unchanged while its fit result is used.

    `logdet_covariance` is log(det(C)), not log(det(W)). `marginal_sigma`
    optionally supplies marginal standard deviations for error bars, not extra
    fit weights. This represents the whole covariance and cannot be combined
    with other observation errors. No Julia dual numbers reach the callback.
    """

    apply: object
    logdet_covariance: float
    marginal_sigma: object = None
    inplace: bool = False

    def _payload(self):
        def whiten(out, residual):
            target, values = np.asarray(out), _readonly(residual)
            if self.inplace:
                if self.apply(target, values) is not None:
                    raise TypeError("in-place whitening must fill out and return None")
            else:
                result = np.asarray(self.apply(values))
                if result.shape != target.shape or np.iscomplexobj(result):
                    raise ValueError("whitening must return one real value per residual")
                np.copyto(target, result)

        marginal = None if self.marginal_sigma is None else _error_values(self.marginal_sigma)
        return (whiten, float(self.logdet_covariance), marginal)
