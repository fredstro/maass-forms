r"""Exact Ford-floor height in the survey runner's meridian-normalised frame.

The package's holonomy matrices and the pure-NumPy survey runner have the same
cusp lattice up to a complex similarity, but the package meridian need not have
length one.  Heights transform by the modulus of that similarity.  This module
uses only the scalar floor height; it does not identify individual sphere
centres across frames, which would require the translational offset ``delta``.
"""

from math import isfinite

from maass_forms_klein.hyperbolic_space.face_pairing import _exact_floor, floor_height
from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup


def certified_floor_height(name: str) -> float:
    r"""Return exact power-vertex ``Y0`` with meridian translation normalised to one.

    The faces come from the horoball-completed Ford-domain path, not the
    word-length-limited fallback.  The package-frame height is divided by the
    modulus of its meridian translation; this is invariant under the unknown
    horizontal frame offset.  Numerical exactness is limited by the holonomy's
    double precision.

    INPUT:

    - ``name`` -- SnapPy manifold identifier, e.g. ``'4_1'``.

    OUTPUT:

    - float; the minimal Ford-floor height in the survey frame.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.modform.survey_floor import certified_floor_height
        sage: bool(abs(certified_floor_height('4_1')**2 - 2/3) < 1e-10)  # long time
        True
    """
    group = KleinianGroup(name)
    meridian = group.named_generators()["M"]
    meridian_translation = complex(meridian[0, 0] * meridian[0, 1])
    scale = abs(meridian_translation)
    if not isfinite(scale) or scale <= 0:
        raise ValueError(f"invalid meridian translation scale for {name}: {scale}")
    faces = group.ford_faces(certified=True)
    package_height = float(floor_height(group, faces=faces, exact=True))
    if not isfinite(package_height) or package_height <= 0:
        raise ValueError(f"invalid Ford floor height for {name}: {package_height}")
    return package_height / scale


def exact_cover_floor(cover, v1: complex, v2: complex) -> float:
    r"""Power-vertex floor of the emitted, centred-cell covering matrices.

    This is exact for the quadtree-certified *selected cover*, which can have a
    slightly lower floor than the complete Ford sphere family.  It is a safe
    fallback when the horoball-frame face enumeration cannot be aligned; it is
    not presented as the exact floor of the full Ford domain.

    INPUT:

    - ``cover`` -- matrices emitted by :func:`hejhal.covering_data`.
    - ``v1``, ``v2`` -- reduced cusp translation basis in the same frame.

    OUTPUT:

    - float; exact power-vertex floor of the selected cover.
    """
    circles = []
    for matrix in cover:
        c, d = complex(matrix[1, 0]), complex(matrix[1, 1])
        if abs(c) > 1e-12:
            circles.append((None, -d / c, 1.0 / abs(c)))
    if not circles:
        raise ValueError("cover has no isometric hemispheres")
    height = _exact_floor(circles, v1, v2)[0]
    if not isfinite(height) or height <= 0:
        raise ValueError("cover has no positive power-vertex floor")
    return height
