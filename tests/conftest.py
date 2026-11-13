"""Shared pytest fixture registration for CLI integration scenarios."""

import pytest

from mdscope_test_support import tiny_cli_project


@pytest.fixture
def project(tmp_path):
    """Build an isolated, six-frame CLI project only for tests requesting it."""
    return tiny_cli_project(tmp_path)
