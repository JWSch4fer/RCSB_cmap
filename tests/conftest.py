
import os
import matplotlib
import pytest

# Use non-interactive backend for tests
matplotlib.use("Agg")

@pytest.fixture(autouse=True)
def _set_cwd_tmp(tmp_path, monkeypatch):
    # run each test in its own temp directory to avoid file collisions
    monkeypatch.chdir(tmp_path)
    # avoid writing cache or huge files in repo root
    os.environ.setdefault("MPLCONFIGDIR", str(tmp_path / ".mplconfig"))
