"""
Test that comp_manager can serialize SageMath types used by maass_forms_hilbert.

These tests document known serialization failures in comp_manager's JSONEncoder,
which lacks support for SageMath types like NumberFieldFractionalIdeal.
"""

import pytest
from comp_manager.utils.serialization import serialize
from sage.rings.number_field.number_field import QuadraticField


def test_serialize_number_field_fractional_ideal():
    """comp_manager must be able to serialize NumberFieldFractionalIdeal objects."""
    K = QuadraticField(2)
    ideal = K.ideal(1)
    serialize((ideal,))


def test_serialize_number_field_element():
    """comp_manager must be able to serialize NumberFieldElement objects."""
    K = QuadraticField(2)
    a = K.gen()
    serialize((a,))


def test_serialize_number_field_order():
    """comp_manager must be able to serialize NumberField ring of integers."""
    K = QuadraticField(2)
    OK = K.ring_of_integers()
    serialize((OK,))