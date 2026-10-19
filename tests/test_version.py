"""Unit tests for __version__.py."""

import libdiffpy  # noqa


def test_package_version():
    """Ensure the package version is defined and not set to the initial
    placeholder."""
    assert hasattr(libdiffpy, "__version__")
    assert libdiffpy.__version__ != "0.0.0"
