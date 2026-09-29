"""Inspect real build artifacts, including the wheel rebuilt from its sdist."""

from email import message_from_bytes
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import zipfile

import pytest
from juliapkg.compat import Compat
from semver import Version
from check_install import validate_core_source
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def distributions(tmp_path_factory):
    directory = tmp_path_factory.mktemp("distributions")
    build = subprocess.run(
        [sys.executable, "-m", "build", "--no-isolation", "--outdir", str(directory), str(ROOT / "python")],
        capture_output=True, text=True, timeout=120,
    )
    assert build.returncode == 0, build.stdout + build.stderr
    return directory


def test_installer_pin_allows_compatible_core_patches():
    core = tomllib.loads((ROOT / "Project.toml").read_text())
    pin = json.loads((ROOT / "python/scientificfitting/juliapkg.json").read_text())["packages"]["ScientificFitting"]
    assert set(pin) == {"uuid", "version"}  # Never ship a local development path.
    assert pin["uuid"] == core["uuid"]
    current, compatible = Version.parse(core["version"]), Compat.parse(pin["version"])
    wrapper = Version.parse(tomllib.loads((ROOT / "python/pyproject.toml").read_text())["project"]["version"])
    assert (wrapper.major, wrapper.minor) == (current.major, current.minor)
    assert current in compatible
    # Core bugfixes must not force a new Python release with an identical wrapper.
    assert current.bump_patch() in compatible
    assert current.bump_minor() not in compatible


def test_wheel_metadata_license_and_bridge(distributions):
    license_text = (ROOT / "LICENSE").read_bytes()
    with zipfile.ZipFile(next(distributions.glob("*.whl"))) as wheel:
        paths = wheel.namelist()
        metadata = message_from_bytes(wheel.read(next(p for p in paths if p.endswith("/METADATA"))))
        assert metadata["Name"] == "scientificfitting"
        assert metadata["Version"] == tomllib.loads((ROOT / "python/pyproject.toml").read_text())["project"]["version"]
        assert metadata["License-Expression"] == "MIT"
        assert metadata["License-File"] == "LICENSE"
        assert metadata["Description-Content-Type"] == "text/markdown"
        assert "NumPy models" in metadata.get_payload()
        assert wheel.read(next(p for p in paths if p.endswith("/licenses/LICENSE"))) == license_text
        runtime = ROOT / "python/scientificfitting"
        expected = {"scientificfitting/" + p.relative_to(runtime).as_posix(): p
                    for p in runtime.rglob("*") if p.is_file() and p.suffix in (".py", ".jl", ".json")}
        assert {p for p in paths if p.startswith("scientificfitting/") and not p.endswith("/")} == set(expected)
        for name, source in expected.items():
            assert wheel.read(name) == source.read_bytes()
        assert not any(p.startswith(("tests/", "docs/", "src/")) for p in paths)
        assert not any(p.endswith((".dylib", ".so", ".dll", ".png")) for p in paths)


def test_source_archive_is_self_contained(distributions):
    with tarfile.open(next(distributions.glob("*.tar.gz"))) as archive:
        prefix = archive.getnames()[0].split("/")[0] + "/"
        for name in ("LICENSE", "README.md", "pyproject.toml", "scientificfitting/_bridge.jl"):
            assert archive.extractfile(prefix + name).read() == (ROOT / "python" / name).read_bytes()


@pytest.mark.parametrize("mode", ["path", "repo"])
def test_registry_check_rejects_retained_overrides(mode):
    info = dict(registry=False, path=False, repo=False, source="checkout", tree_hash="a"*40)
    info[mode] = True
    with pytest.raises(AssertionError, match="registry-installed core"):
        validate_core_source(info, None)


def test_registry_check_requires_a_resolved_tree():
    info = dict(registry=True, path=False, repo=False, source="package", tree_hash=None)
    with pytest.raises(AssertionError, match="registry-installed core"):
        validate_core_source(info, None)
    info["tree_hash"] = "a"*40
    validate_core_source(info, None)


def test_development_check_requires_the_requested_source(tmp_path):
    info = dict(registry=False, path=True, repo=False, source=str(tmp_path), tree_hash=None)
    validate_core_source(info, tmp_path)
    with pytest.raises(AssertionError, match="requested development checkout"):
        validate_core_source(info, tmp_path / "different")
    info["path"] = False
    with pytest.raises(AssertionError, match="requested development checkout"):
        validate_core_source(info, tmp_path)
