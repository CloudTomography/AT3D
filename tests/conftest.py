"""pytest configuration for AT3D tests.

Tests use relative paths (e.g. ``../default_config.json``, ``./data/...``)
that assume the working directory is the ``tests/`` folder. This fixture
changes CWD to ``tests/`` at the start of the session and restores it
afterwards so that the paths resolve correctly regardless of where pytest
is invoked from.
"""
import os
import pathlib

import pytest

_TESTS_DIR = pathlib.Path(__file__).resolve().parent


@pytest.fixture(autouse=True, scope="session")
def _chdir_to_tests_dir():
    original = os.getcwd()
    os.chdir(_TESTS_DIR)
    yield
    os.chdir(original)
