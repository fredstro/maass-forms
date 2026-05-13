"""
Eigenvalue search routines for Kleinian Maass forms.

This subpackage provides parallel, database-backed eigenvalue search
over a one-dimensional grid of spectral parameters s = 1/2 + i*r.

It is the Kleinian (1D) analogue of :mod:`maass_forms_hilbert.search.search`,
where the spectral parameter is multi-dimensional. For the Kleinian case the
spectral parameter is a single complex number, so the grid is one-dimensional
and Broyden's method reduces to the secant method.
"""

from maass_forms_klein.search.search import (
    ARITHMETIC_LOCATORS,
    GENERAL_LOCATORS,
    check_coefficients_of_computed_object,
    coeff_diff_fun,
    compute_on_interval,
    compute_one_spectral_parameter,
    create_grid,
    distance_between_points,
    is_arithmetic_group,
    search_eigenvalues_via_relation,
    secant_iteration,
    sorting,
)

__all__ = [
    "ARITHMETIC_LOCATORS",
    "GENERAL_LOCATORS",
    "check_coefficients_of_computed_object",
    "coeff_diff_fun",
    "compute_on_interval",
    "compute_one_spectral_parameter",
    "create_grid",
    "distance_between_points",
    "is_arithmetic_group",
    "search_eigenvalues_via_relation",
    "secant_iteration",
    "sorting",
]
