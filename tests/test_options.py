"""Unit tests for potpyri.utils.options."""
import os

import pytest

from potpyri.utils import options
from potpyri.instruments import instrument_getter


def test_init_options():
    """init_options returns parser with instrument choices."""
    parser = options.init_options()
    assert parser is not None
    for action in parser._get_positional_actions():
        if getattr(action, "dest", None) == "instrument":
            assert hasattr(action, "choices")
            assert "GMOS" in action.choices
            assert "LRIS" in action.choices
            return
    pytest.fail("instrument positional not found")


def test_add_paths_raises_when_data_path_missing():
    """add_paths raises when data_path does not exist."""
    tel = instrument_getter("GMOS")
    with pytest.raises(Exception, match="does not exist"):
        options.add_paths("/nonexistent/path/12345", "files.txt", tel)


def test_add_paths_creates_dirs(tmp_path):
    """add_paths creates raw, bad, red, log, cal, work and returns paths dict."""
    tel = instrument_getter("GMOS")
    paths = options.add_paths(str(tmp_path), "files.txt", tel)
    assert "data" in paths
    assert paths["data"] == os.path.abspath(tmp_path)
    assert "raw" in paths
    assert "red" in paths
    assert "log" in paths
    assert "cal" in paths
    assert "work" in paths
    assert "filelist" in paths
    assert "code" in paths
    assert os.path.exists(paths["raw"])
    assert os.path.exists(paths["red"])
    assert os.path.exists(paths["log"])
