"""Independent gallery check with SciPy and C++ Minuit2 (no Julia/SF import).

Run from a checkout with numpy, scipy and iminuit installed. Reads the literal
counts from the executable page; its linked ROOT extraction retains provenance.
"""
from pathlib import Path
import re
import numpy as np
from scipy.special import gammaln, ndtr
from scipy.stats import chi2
from iminuit import Minuit


def expected(edges, center, width, width_ratio, core_fraction,
             background_scale, signal_yield, background_yield):
    """Closed-form, window-normalized bin integrals, including the shared mean."""
    core = ndtr((edges - center) / width)
    tail = ndtr((edges - center) / (width * width_ratio))
    signal_cdf = core_fraction * core + (1-core_fraction) * tail
    signal = np.diff(signal_cdf) / (signal_cdf[-1] - signal_cdf[0])
    background_cdf = -np.expm1(-(edges-edges[0]) / background_scale)
    background = np.diff(background_cdf) / background_cdf[-1]
    return signal_yield*signal + background_yield*background


def compare():
    page = Path(__file__).resolve().parents[1] / "docs/src/gallery/lhcb_mass_spectrum.md"
    text = re.search(r"counts = \[\n(.*?)\n\]", page.read_text(), flags=re.S).group(1)
    counts = np.fromstring(text.replace("\n", ""), sep=",")
    edges = np.arange(5200., 5601., 5.)
    assert len(counts) == 80 and counts.sum() == 7368

    def cost(center, width, width_ratio, core_fraction, background_scale,
             signal_yield, background_yield):
        means = expected(edges, center, width, width_ratio, core_fraction,
                         background_scale, signal_yield, background_yield)
        return 2*np.sum(means-counts*np.log(means)+gammaln(counts+1))

    for single in (True, False):
        fit = Minuit(cost, 5284., 14., 2., 1. if single else .8, 400., 6500., 1000.)
        fit.limits = [(5250,5310), (3,30), (1.05,5), (.05,1),
                      (20,5000), (0,15000), (0,15000)]
        fit.fixed["core_fraction"] = fit.fixed["width_ratio"] = single
        fit.errordef, fit.tol = 1., 1e-5  # Cost convention is -2 log L.
        fit.migrad().hesse()
        assert fit.valid and fit.accurate and not fit.fmin.has_parameters_at_limit
        means = expected(edges, *fit.values)
        deviance = 2*np.sum(means-counts+counts*np.log(counts/means))
        target = 101.9637663 if single else 69.43671155
        assert abs(deviance-target) < 1e-5
        print("\none width" if single else "\ntwo widths")
        print(f"D/ndf = {deviance:.8f}/{80-fit.nfit}; p = {chi2.sf(deviance, 80-fit.nfit):.8f}")
        # Local curvature errors, not Minuit's bound-transformed errors.
        for name in fit.parameters:
            print(f"{name:18s} = {fit.values[name]:.10g} +/- {np.sqrt(fit.covariance[name,name]):.7g}")
        if not single:
            fit.minos("signal_yield")
            interval = fit.merrors["signal_yield"]
            assert interval.lower_valid and interval.upper_valid
            print("Ns MINOS interval:", fit.values["signal_yield"]+interval.lower,
                  fit.values["signal_yield"]+interval.upper)


if __name__ == "__main__":
    compare()
