"""
maass-form-core: Shared infrastructure for Maass forms computation.

Provides common database patterns, coefficient management, mathematical space
abstractions, and utility functions used by domain-specific packages
(maass-forms-hilbert, maass-forms-klein).
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("maass_form_core")
except PackageNotFoundError:
    __version__ = "0.0.0"
