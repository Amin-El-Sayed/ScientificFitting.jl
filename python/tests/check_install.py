"""Check an installed base wheel, with real fits and no optional Python packages.

Run in a fresh environment after `pip install <wheel>`. `--source` stages an
unreleased Julia checkout via JuliaPkg; omitting it requires a registry install
and rejects persisted development or repository overrides.
It never supplies or installs a Julia executable itself.
"""

import argparse
import importlib.util
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np


def validate_core_source(info, source):
    """Check Pkg's resolved provenance, not just the requested installation mode."""
    if source is None:
        assert info["registry"] and not info["path"] and not info["repo"] and info["tree_hash"], (
            "expected a registry-installed core; use a new JuliaPkg environment "
            "without persisted source overrides, or pass --source explicitly"
        )
    else:
        assert info["path"] and Path(info["source"]).resolve() == Path(source).resolve(), (
            "loaded core does not match the requested development checkout"
        )


def check(source=None):
    import scientificfitting as sf

    assert "site-packages" in Path(sf.__file__).parts, "test the installed wheel, not an editable checkout"
    assert "juliacall" not in sys.modules
    assert importlib.util.find_spec("matplotlib") is None
    assert importlib.util.find_spec("scipy") is None
    if source:
        import juliapkg
        # Only the source selection differs from the eventual registry install.
        juliapkg.add("ScientificFitting", uuid="4ec0e560-7423-4ff9-937f-0d4da6f5f8f8",
                     path=str(Path(source).resolve()), dev=True)

    x = np.array([-1., 0., 1., 2., 3., 4.])
    y = np.array([-1.1, 0.8, 2.9, 4.7, 6.8, 8.7])
    design = np.column_stack([x, np.ones_like(x)])
    expected = np.linalg.lstsq(design, y, rcond=None)[0]
    covariance = np.linalg.inv(design.T @ design) * 0.2**2
    def fit():
        return sf.fit_model(lambda x, slope, offset: slope*x + offset, x, y,
                            p0={"slope": 1., "offset": 0.}, sigma_y=0.2)

    start = perf_counter()
    result = fit()
    first_fit = perf_counter() - start
    assert result.converged
    np.testing.assert_allclose(result.params, expected, atol=1e-6)
    np.testing.assert_allclose(result.covariance, covariance, rtol=1e-5)
    assert "slope" in result.report()
    assert "Fit diagnostic dashboard" in result.diagnose()
    start = perf_counter()
    fit()
    repeat_fit = perf_counter() - start
    scan = result.profile("slope", values=result.params[0] + result.stderr[0]*np.array([-1., 0., 1.]))
    np.testing.assert_allclose(scan.delta_cost, [1., 0., 1.], atol=3e-5)
    counts = sf.fit_poisson_model(lambda x, rate: np.full_like(x, rate), [0, 1, 2], [4, 6, 8],
                                  p0={"rate": 5.}, bounds={"rate": (0.1, 20.)})
    assert counts.converged
    np.testing.assert_allclose(counts.params, [6.], atol=1e-5)
    np.testing.assert_allclose(counts.covariance, [[2.]], rtol=2e-5)
    from scientificfitting._runtime import _backend
    backend = _backend()
    core_install = dict(backend.seval('''
        using Pkg
        let info = Pkg.dependencies()[Base.PkgId(ScientificFitting).uuid]
            pydict(registry=info.is_tracking_registry, path=info.is_tracking_path,
                   repo=info.is_tracking_repo, source=info.source, tree_hash=info.tree_hash)
        end
    '''))
    validate_core_source(core_install, source)
    assert not backend.seval('any(m -> nameof(m) in (:Makie, :CairoMakie), values(Base.loaded_modules))')
    assert "matplotlib" not in sys.modules and "scipy" not in sys.modules
    return dict(platform=platform.platform(), python=platform.python_version(),
                julia=backend.seval("string(VERSION)"), core=backend.seval("string(pkgversion(ScientificFitting))"),
                package=str(Path(sf.__file__).resolve()), source_argument=bool(source),
                core_source=backend.seval("pathof(ScientificFitting)"),
                core_install=core_install,
                first_fit_seconds=first_fit, repeat_fit_seconds=repeat_fit,
                maximum_parameter_error=float(np.max(np.abs(result.params-expected))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = json.dumps(check(args.source), indent=2)
    print(report)
    if args.output:
        args.output.write_text(report + "\n", encoding="utf-8")
