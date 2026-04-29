"""
Mathematical space base classes for Maass forms.

Provides abstract base classes for mathematical spaces and form elements
that domain-specific packages inherit from.
"""

from maass_form_core.spaces.base_element import MaassFormElement
from maass_form_core.spaces.base_space import MaassFormSpace

__all__ = [
    "MaassFormSpace",
    "MaassFormElement",
]
