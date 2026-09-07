"""Inspect real build artifacts, including the wheel rebuilt from its sdist."""

from email import message_from_bytes
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import zipfile

import pytest
from packaging.version import Version
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


def test_installer_pin_matches_core_version():
    core = tomllib.loads((ROOT / "Project.toml").read_text())
    python = tomllib.loads((ROOT / "python/pyproject.toml").read_text())["project"]
    pin = json.loads((ROOT / "python/scientificfitting/juliapkg.json").read_text())["packages"]["ScientificFitting"]
    assert pin == {"uuid": core["uuid"], "version": "~" + core["version"]}
    assert Version(python["version"]).release == Version(core["version"]).release


def test_wheel_metadata_license_and_bridge(distributions):
    license_text = (ROOT / "LICENSE").read_bytes()
    with zipfile.ZipFile(next(distributions.glob("*.whl"))) as wheel:
        paths = wheel.namelist()
        metadata = message_from_bytes(wheel.read(next(p for p in paths if p.endswith("/METADATA"))))
        assert metadata["Name"] == "scientificfitting"
        assert metadata["License-Expression"] == "MIT"
        assert metadata["License-File"] == "LICENSE"
        assert metadata["Description-Content-Type"] == "text/markdown"
        assert "NumPy models" in metadata.get_payload()
        assert wheel.read(next(p for p in paths if p.endswith("/licenses/LICENSE"))) == license_text
        for name in ("_bridge.jl", "juliapkg.json"):
            assert wheel.read("scientificfitting/" + name) == (ROOT / "python/scientificfitting" / name).read_bytes()
        assert not any(p.startswith(("tests/", "docs/", "src/")) for p in paths)
        assert not any(p.endswith((".dylib", ".so", ".dll", ".png")) for p in paths)


def test_source_archive_is_self_contained(distributions):
    with tarfile.open(next(distributions.glob("*.tar.gz"))) as archive:
        prefix = archive.getnames()[0].split("/")[0] + "/"
        for name in ("LICENSE", "README.md", "pyproject.toml", "scientificfitting/_bridge.jl"):
            assert archive.extractfile(prefix + name).read() == (ROOT / "python" / name).read_bytes()
