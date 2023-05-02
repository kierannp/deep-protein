"""
Unit and regression test for the deep_protein package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import deep_protein


def test_deep_protein_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "deep_protein" in sys.modules
