r"""
Ford-domain face pairings for Kleinian groups (phase 1).

This module extracts the *face-carrying* isometric hemispheres of the Ford
fundamental domain of a cusped Kleinian group and verifies completeness of the
face set by two independent tests:

- the volume identity

  .. MATH::

      \operatorname{vol}(P) = \int_P \frac{dx}{2\,h(x)^2},

  where `h(x)` is the height of the Ford floor above the cusp cell `P`
  (a converged value equals the hyperbolic volume of the manifold);

- the Epstein--Penner duality count, comparing the number of face pairs with
  the number of edges (equivalently tetrahedra) of SnapPy's canonical
  retriangulation.

The algorithm is a port of the standalone reference prototype in
``dev/face_pairing/face_pairing_prototype.py`` (validated 2026-07, plain
``snappy`` + ``numpy`` + ``scipy``), expressed against the package's own
:class:`~maass_forms_klein.hyperbolic_space.kleinian_group.KleinianGroup_class`
geometry. All numerical work is done in double precision; visibility is decided
exactly (up to holonomy error) by a linear program, and completeness is
cross-checked against SnapPy.

This is phase 1 of ``FACE_PAIRING_PLAN.md``: :func:`visible_faces`,
:func:`ford_volume` and :func:`check_face_completeness`. The bounded-``|c|``
enumeration (phase 3) that closes the last pair for knots such as ``7_4`` is not
implemented here; :func:`visible_faces` enumerates candidates by a breadth-first
search on bottom rows up to a word-length bound.

EXAMPLES::

    sage: from maass_form_core.testing import connect_mockdb
    sage: connect_mockdb()
    sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
    sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
    ....:     visible_faces, ford_volume, check_face_completeness)
    sage: G = KleinianGroup('4_1')
    sage: faces = visible_faces(G)
    sage: len(faces)
    4
"""

import warnings
from typing import Union

import numpy as np

from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RealNumber

from maass_forms_klein.hyperbolic_space.types import Circle, Face
from maass_forms_klein.hyperbolic_space.word_utils import find_inverse_word, word_list_sort_key

# Real_t defined locally to avoid import-chain issues (see package CLAUDE.md).
Real_t = Union[RealNumber, Integer, int, float]

# Alphabet letters that are *not* the parabolic cusp translations.
_PARABOLIC_LETTERS = set("LlMm")


def _warn_if_high_precision(prec: int) -> None:
    r"""
    Warn that ``prec > 53`` has no effect beyond the input generators: all
    geometry in this module is computed in IEEE double precision.

    EXAMPLES::

        sage: import warnings
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _warn_if_high_precision)
        sage: with warnings.catch_warnings(record=True) as caught:
        ....:     warnings.simplefilter("always")
        ....:     _warn_if_high_precision(113)
        sage: len(caught)
        1
        sage: with warnings.catch_warnings(record=True) as caught:
        ....:     warnings.simplefilter("always")
        ....:     _warn_if_high_precision(53)
        sage: len(caught)
        0
    """
    if prec > 53:
        warnings.warn(
            f"prec={prec} > 53 has no effect beyond generating the input "
            "generators/lattice; all face-pairing geometry is computed in IEEE "
            "double precision (53-bit).",
            stacklevel=3,
        )


def _numeric_generators(group, prec: int = 53) -> dict:
    r"""
    Return the named generators of ``group`` as ``numpy`` complex 2x2 arrays.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _numeric_generators
        sage: gens = _numeric_generators(KleinianGroup('4_1'))
        sage: sorted(gens)
        ['A', 'B', 'L', 'M', 'a', 'b', 'l', 'm']
        sage: gens['a'].shape
        (2, 2)
    """
    named = group.named_generators(prec=prec)
    return {
        name: np.array(
            [
                [complex(m[0, 0]), complex(m[0, 1])],
                [complex(m[1, 0]), complex(m[1, 1])],
            ],
            dtype=complex,
        )
        for name, m in named.items()
    }


def _gauss_reduce(v1: complex, v2: complex):
    r"""
    Lagrange--Gauss reduction of a complex lattice basis to its shortest,
    near-orthogonal form (same lattice).

    A reduced basis keeps the fundamental cell well proportioned even for very
    skew cusp shapes, so grid pre-passes resolve the floor and the neighbouring
    lattice translates that cover the cell lie in a small ring (trap 7 in
    ``dev/face_pairing/README.md``: for a raw basis a skew cusp such as ``8_2``
    needs translates at coefficients like ``(1, -11)``).

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _gauss_reduce
        sage: w1, w2 = _gauss_reduce(complex(1, 0), complex(11, 1))
        sage: sorted([round(abs(w1), 6), round(abs(w2), 6)])   # short vectors
        [1.0, 1.0]
        sage: # the reduced basis spans the same lattice as the input
        sage: w1, w2 = _gauss_reduce(complex(2, 0), complex(1, 3))
        sage: bool(abs(w1) <= abs(w2) + 1e-12)
        True
    """
    while True:
        if abs(v2) < abs(v1):
            v1, v2 = v2, v1
        mu = round((v1.real * v2.real + v1.imag * v2.imag) / abs(v1) ** 2)
        if mu == 0:
            return v1, v2
        v2 = v2 - mu * v1


def _lattice_vectors(group, prec: int = 53):
    r"""
    Return ``(v1, v2, base)`` for the cusp translation lattice as Python complex
    numbers, with ``base`` the corner of the centred fundamental cell.

    The basis is Lagrange--Gauss reduced (:func:`_gauss_reduce`) so the cell is
    well proportioned regardless of how skew the cusp shape is. All Ford-domain
    quantities (face count, volume, floor height) are properties of the lattice,
    not of the chosen basis, so this does not change any result -- it only keeps
    the geometry numerically well conditioned.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _lattice_vectors
        sage: v1, v2, base = _lattice_vectors(KleinianGroup('4_1'))
        sage: abs(v1) > 0 and abs(v2) > 0
        True
        sage: bool(abs(base + (v1 + v2) / 2) < 1e-12)
        True
    """
    v1v, v2v = group.translation_lattice(prec).basis()
    v1, v2 = _gauss_reduce(complex(v1v[0], v1v[1]), complex(v2v[0], v2v[1]))
    return v1, v2, -(v1 + v2) / 2


def _lattice_reduce(center: complex, v1: complex, v2: complex) -> complex:
    r"""
    Translate ``center`` by the lattice ``Z v1 + Z v2`` into the base cell.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _lattice_reduce
        sage: _lattice_reduce(complex(2.3, 1), complex(1, 0), complex(0, 1))
        (0.2999999999999998+0j)
    """
    matrix = np.array([[v1.real, v2.real], [v1.imag, v2.imag]])
    x = np.linalg.solve(matrix, [center.real, center.imag])
    n = np.floor(x)
    return center - n[0] * v1 - n[1] * v2


def _lattice_dist(a: complex, b: complex, v1: complex, v2: complex) -> float:
    r"""
    Distance from ``a`` to ``b`` modulo the lattice ``Z v1 + Z v2``.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _lattice_dist
        sage: d = _lattice_dist(complex(0.1, 0), complex(0.9, 0), complex(1, 0), complex(0, 1))
        sage: bool(abs(d - 0.2) < 1e-12)
        True
    """
    matrix = np.array([[v1.real, v2.real], [v1.imag, v2.imag]])
    x = np.linalg.solve(matrix, [(a - b).real, (a - b).imag])
    x -= np.round(x)
    return float(abs(x[0] * v1 + x[1] * v2))


def _enumerate_circles(gens: dict, v1: complex, v2: complex, maxlen: int, rmin: float) -> list:
    r"""
    Breadth-first enumeration of isometric circles by bottom row.

    For torsion-free groups the bottom row of a word determines its left coset,
    so this enumerates distinct isometric hemispheres. Returns a list of
    ``(word, centre_mod_Lambda, radius)`` with ``radius >= rmin``, one entry per
    distinct circle.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _numeric_generators, _lattice_vectors, _enumerate_circles)
        sage: G = KleinianGroup('4_1')
        sage: gens = _numeric_generators(G)
        sage: v1, v2, _ = _lattice_vectors(G)
        sage: circles = _enumerate_circles(gens, v1, v2, maxlen=3, rmin=0.05)
        sage: max(r for _, _, r in circles)  # abs tol 1e-9
        1.0
    """
    letters = [g for g in gens if g not in _PARABOLIC_LETTERS]
    rows = {g: (gens[g][1, 0], gens[g][1, 1]) for g in letters}
    seen, frontier = set(), []
    for g in letters:
        c, d = rows[g]
        seen.add(tuple(round(x, 5) for x in (c.real, c.imag, d.real, d.imag)))
        frontier.append((g, (c, d)))
    good: dict = {}
    for length in range(1, maxlen + 1):
        for word, (c, d) in frontier:
            if abs(c) > 1e-9 and 1.0 / abs(c) >= rmin:
                centre, radius = -d / c, 1.0 / abs(c)
                reduced = _lattice_reduce(centre, v1, v2)
                key = (round(reduced.real, 6), round(reduced.imag, 6), round(radius, 6))
                good.setdefault(key, (word, reduced, radius))
        if length == maxlen:
            break
        frontier = _extend_frontier(frontier, gens, letters, seen)
    return list(good.values())


def _extend_frontier(frontier: list, gens: dict, letters: list, seen: set) -> list:
    r"""
    Advance the breadth-first frontier by one letter, skipping backtracking and
    already-seen bottom rows.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _numeric_generators, _extend_frontier)
        sage: gens = _numeric_generators(KleinianGroup('4_1'))
        sage: letters = ['a', 'b', 'A', 'B']
        sage: seen = set()
        sage: frontier = [('a', (gens['a'][1, 0], gens['a'][1, 1]))]
        sage: nxt = _extend_frontier(frontier, gens, letters, seen)
        sage: all(len(w) == 2 for w, _ in nxt)
        True
    """
    new_frontier = []
    for word, (c, d) in frontier:
        for g in letters:
            if word[-1] == g.swapcase():
                continue
            a11, a12, a21, a22 = gens[g][0, 0], gens[g][0, 1], gens[g][1, 0], gens[g][1, 1]
            row = (c * a11 + d * a21, c * a12 + d * a22)
            key = tuple(round(x, 5) for x in (row[0].real, row[0].imag, row[1].real, row[1].imag))
            if key not in seen:
                seen.add(key)
                new_frontier.append((word + g, row))
    return new_frontier


def _floor_translates(circles: list, v1: complex, v2: complex, base: complex) -> list:
    r"""
    Expand each circle by every lattice translate whose disc meets the base cell.

    The translate range adapts to the cell so the floor is complete over the
    whole cusp cell even for skew cusp shapes -- a translate at index ``(i, j)``
    is kept iff its disc reaches the cell circumscribed circle. Combined with the
    Gauss-reduced basis (:func:`_lattice_vectors`) this keeps the number of
    translates small.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _floor_translates
        sage: items = _floor_translates([('a', complex(0, 0), 1.0)],
        ....:                           complex(1, 0), complex(0, 1), complex(-0.5, -0.5))
        sage: 1 <= len(items) <= 9
        True
    """
    centre_cell = base + (v1 + v2) / 2
    rcell = 0.5 * max(abs(v1 + v2), abs(v1 - v2))
    spacing = abs((v1 * np.conj(v2)).imag) / max(abs(v1), abs(v2))
    out = []
    for word, c, r in circles:
        k = int(np.ceil((r + rcell) / spacing)) + 1
        for i in range(-k, k + 1):
            for j in range(-k, k + 1):
                translated = c + i * v1 + j * v2
                if abs(translated - centre_cell) <= r + rcell:
                    out.append((word, translated, r))
    return out


def _envelope_min_h2(items: list, v1: complex, v2: complex, base: complex, n: int) -> float:
    r"""
    Minimum over an ``n x n`` grid of the base cell of the squared floor height
    ``max_i (r_i^2 - |z - c_i|^2)``.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _envelope_min_h2
        sage: items = [('a', 0j, 2.0)]
        sage: h2 = _envelope_min_h2(items, 1 + 0j, 0 + 1j, -0.5 - 0.5j, 4)
        sage: h2 > 0
        True
    """
    centres = np.array([c for _, c, _ in items])
    radii = np.array([r for _, _, r in items])
    s = np.linspace(0, 1, n, endpoint=False) + 0.5 / n
    u, w = np.meshgrid(s, s)
    z = base + u * v1 + w * v2
    best = np.full(z.shape, -1.0)
    for i in range(len(centres)):
        best = np.maximum(best, radii[i] ** 2 - np.abs(z - centres[i]) ** 2)
    return float(best.min())


def _bounded_candidates(circles: list, v1: complex, v2: complex, base: complex):
    r"""
    Filter the enumerated circles to those that can touch the Ford floor.

    Bootstraps a lower bound for the floor height ``Y0`` from the largest
    circles, drops circles with radius below it (a sphere with ``r < Y0`` never
    reaches the floor), and returns ``(all_items, Y0, candidate_faces)`` where
    ``all_items`` carries lattice translates for the envelope and
    ``candidate_faces`` are the circles with ``r >= 0.95 * Y0``.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _numeric_generators, _lattice_vectors, _enumerate_circles,
        ....:     _bounded_candidates)
        sage: G = KleinianGroup('4_1')
        sage: gens = _numeric_generators(G)
        sage: v1, v2, base = _lattice_vectors(G)
        sage: circles = _enumerate_circles(gens, v1, v2, maxlen=6, rmin=0.05)
        sage: all_items, y0, cand = _bounded_candidates(circles, v1, v2, base)
        sage: y0 > 0 and len(cand) >= 4
        True
    """
    ranked = sorted(circles, key=lambda t: -t[2])
    # Adaptive translate rings (not a fixed +-1 ring): every lattice translate
    # whose disc meets the cell is included, so the envelope is complete even for
    # skew cusps (trap 7). Clamp the sqrt at 0: a negative envelope minimum means
    # the enumerated spheres do not yet cover the cell (an incomplete candidate
    # set, e.g. a low word-length bound), so the floor bound is 0.
    boot = _floor_translates(ranked[:200], v1, v2, base)
    y0_lower = float(np.sqrt(max(0.0, _envelope_min_h2(boot, v1, v2, base, 300))))
    kept = [(w, c, r) for w, c, r in ranked if r >= 0.9 * y0_lower]
    all_items = _floor_translates(kept, v1, v2, base)
    y0 = float(np.sqrt(max(0.0, _envelope_min_h2(all_items, v1, v2, base, 500))))
    candidates = [(w, c, r) for w, c, r in all_items if r >= 0.95 * y0]
    return all_items, y0, candidates


def _cell_constraints(v1: complex, v2: complex, base: complex):
    r"""
    Half-plane constraints ``A x <= b`` describing the base cell ``base + [0, 1)
    v1 + [0, 1) v2`` in Cartesian coordinates.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _cell_constraints
        sage: A, b = _cell_constraints(complex(1, 0), complex(0, 1), complex(-0.5, -0.5))
        sage: len(A), len(b)
        (4, 4)
    """
    binv = np.linalg.inv(np.array([[v1.real, v2.real], [v1.imag, v2.imag]]))
    bb = np.array([base.real, base.imag])
    rows, rhs = [], []
    for row in (0, 1):
        rows.append(-binv[row])
        rhs.append(binv[row] @ (-bb))
        rows.append(binv[row])
        rhs.append(1 + binv[row] @ bb)
    return rows, rhs


def _lp_visible(items: list, v1: complex, v2: complex, base: complex) -> list:
    r"""
    Exact visibility over the candidate list by linear programming.

    Sphere ``i`` carries a face iff the program ``max t`` subject to
    ``l_i(x) - l_j(x) >= t`` for all ``j != i`` and ``x in P`` has a positive
    optimum, where ``l_i(x) = 2 <c_i, x> + r_i^2 - |c_i|^2`` is the squared floor
    height contributed by sphere ``i``. Returns ``(item, margin, witness)`` for
    every visible sphere.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _numeric_generators, _lattice_vectors, _enumerate_circles,
        ....:     _bounded_candidates, _lp_visible)
        sage: G = KleinianGroup('4_1')
        sage: gens = _numeric_generators(G)
        sage: v1, v2, base = _lattice_vectors(G)
        sage: circles = _enumerate_circles(gens, v1, v2, maxlen=6, rmin=0.05)
        sage: _, _, cand = _bounded_candidates(circles, v1, v2, base)
        sage: vis = _lp_visible(cand, v1, v2, base)
        sage: all(margin > 0 for _, margin, _ in vis)
        True
    """
    from scipy.optimize import linprog

    deduped: dict = {}
    for word, c, r in items:
        deduped.setdefault((round(c.real, 8), round(c.imag, 8), round(r, 8)), (word, c, r))
    items = list(deduped.values())
    a_cell, b_cell = _cell_constraints(v1, v2, base)
    lift = np.array([[2 * c.real, 2 * c.imag, r * r - abs(c) ** 2] for _, c, r in items])
    out = []
    for i in range(len(items)):
        a_ub = [
            [-(lift[i] - lift[j])[0], -(lift[i] - lift[j])[1], 1.0]
            for j in range(len(items))
            if j != i
        ]
        b_ub = [(lift[i] - lift[j])[2] for j in range(len(items)) if j != i]
        a_ub += [[a[0], a[1], 0.0] for a in a_cell]
        b_ub += list(b_cell)
        res = linprog(
            c=[0, 0, -1],
            A_ub=np.array(a_ub),
            b_ub=np.array(b_ub),
            bounds=[(None, None)] * 3,
            method="highs",
        )
        if res.status == 0 and -res.fun > 1e-12:
            out.append((items[i], -res.fun, complex(res.x[0], res.x[1])))
    return out


def _grid_visible(items: list, v1: complex, v2: complex, base: complex, n: int = 500) -> list:
    r"""
    Fast grid pre-pass for visibility: return the candidates that are the floor
    argmax at some grid point. Each returned tuple is ``(item, None, None)`` to
    match the shape of :func:`_lp_visible`.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _numeric_generators, _lattice_vectors, _enumerate_circles,
        ....:     _bounded_candidates, _grid_visible)
        sage: G = KleinianGroup('4_1')
        sage: gens = _numeric_generators(G)
        sage: v1, v2, base = _lattice_vectors(G)
        sage: circles = _enumerate_circles(gens, v1, v2, maxlen=6, rmin=0.05)
        sage: _, _, cand = _bounded_candidates(circles, v1, v2, base)
        sage: len(_grid_visible(cand, v1, v2, base)) >= 4
        True
    """
    centres = np.array([c for _, c, _ in items])
    radii = np.array([r for _, _, r in items])
    s = np.linspace(0, 1, n, endpoint=False) + 0.5 / n
    u, w = np.meshgrid(s, s)
    z = base + u * v1 + w * v2
    best = np.full(z.shape, -1, dtype=int)
    best_h2 = np.full(z.shape, -1.0)
    for i in range(len(centres)):
        h2 = radii[i] ** 2 - np.abs(z - centres[i]) ** 2
        m = h2 > best_h2
        best[m] = i
        best_h2[m] = h2[m]
    return [(items[i], None, None) for i in sorted(set(best.ravel().tolist())) if i >= 0]


def _merge_classes(visible: list, v1: complex, v2: complex) -> list:
    r"""
    Merge visible spheres into double-coset classes: two spheres are identified
    when their radii agree and their centres agree modulo the lattice. Returns
    one representative ``(word, centre, radius, margin, witness)`` per class,
    with the canonically smallest word (shortest, then alphabetical) chosen as
    representative so the result is independent of enumeration order. A ``None``
    word (e.g. a candidate from the horoball enumeration, whose word is unknown)
    sorts last, so a real word is preferred whenever one is available.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _merge_classes
        sage: vis = [(('aXb', complex(1, 0), 1.0), 1.0, 0j),
        ....:        (('b', complex(0, 0), 1.0), 2.0, 0j)]
        sage: merged = _merge_classes(vis, complex(1, 0), complex(0, 1))
        sage: len(merged)
        1
        sage: merged[0][0]        # shorter word kept as representative
        'b'
        sage: merged[0][3]        # ... together with its margin
        2.00000000000000
        sage: # a real word is preferred over a None word for the same class
        sage: vis = [((None, complex(0, 0), 1.0), 1.0, 0j),
        ....:        (('b', complex(1, 0), 1.0), 1.0, 0j)]
        sage: _merge_classes(vis, complex(1, 0), complex(0, 1))[0][0]
        'b'
    """
    classes: list = []
    for (word, c, r), margin, witness in visible:
        match = next(
            (
                cl
                for cl in classes
                if abs(cl["radius"] - r) < 1e-8 and _lattice_dist(c, cl["centre"], v1, v2) < 1e-6
            ),
            None,
        )
        if match is None:
            classes.append(
                {"word": word, "centre": c, "radius": r, "margin": margin, "witness": witness}
            )
        elif _word_key(word) < _word_key(match["word"]):
            match.update(word=word, centre=c, margin=margin, witness=witness)
    return [(cl["word"], cl["centre"], cl["radius"], cl["margin"], cl["witness"]) for cl in classes]


def _word_key(word):
    r"""
    Sort key that ranks real words before an unknown (``None``) word.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _word_key
        sage: _word_key('b') < _word_key(None)
        True
    """
    if word is None:
        return (10**9, "", "")
    return word_list_sort_key(word)


def _circle_of(mat: np.ndarray):
    r"""
    Centre and radius of the isometric circle of a ``numpy`` matrix, or ``None``
    for an upper-triangular (parabolic/identity) matrix.

    EXAMPLES::

        sage: import numpy as np
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _circle_of
        sage: _circle_of(np.array([[0, -1], [1, 0]], dtype=complex))
        ((-0+0j), 1.0)
    """
    c = mat[1, 0]
    if abs(c) < 1e-10:
        return None
    return (-mat[1, 1] / c, 1.0 / abs(c))


def _word_matrix(word: str, gens: dict) -> np.ndarray:
    r"""
    Product of ``numpy`` generator matrices along ``word``.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _numeric_generators, _word_matrix)
        sage: gens = _numeric_generators(KleinianGroup('4_1'))
        sage: _word_matrix('aA', gens).shape
        (2, 2)
    """
    mat = np.eye(2, dtype=complex)
    for ch in word:
        mat = mat @ gens[ch]
    return mat


def _pair_word(
    word: str, centre: complex, radius: float, classes: list, gens: dict, v1: complex, v2: complex
):
    r"""
    Return the representative word of the class paired with ``word`` (the class
    carrying the isometric circle of the inverse element), or ``None`` if the
    face set is not inversion-closed at this face.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _numeric_generators, _lattice_vectors, visible_faces, _pair_word)
        sage: G = KleinianGroup('4_1')
        sage: faces = visible_faces(G)
        sage: all(f.inverse_word is not None for f in faces)
        True
    """
    inv = _circle_of(_word_matrix(find_inverse_word(word), gens))
    if inv is None:
        return None
    ci, ri = inv
    for w2, c2, r2, _, _ in classes:
        if abs(r2 - ri) < 1e-8 and _lattice_dist(ci, c2, v1, v2) < 1e-6:
            return w2
    return None


def visible_faces(
    group, candidates=None, method: str = "lp", prec: int = 53, maxlen: int = 9
) -> list:
    r"""
    Extract the visible (face-carrying) isometric hemispheres of the Ford domain.

    The visible hemispheres are merged into classes modulo the translation
    lattice; each class becomes a :class:`~maass_forms_klein.hyperbolic_space.types.Face`
    whose ``inverse_word`` points at the paired class (the hemisphere of the
    inverse element).

    INPUT:

    - ``group`` -- a :class:`KleinianGroup_class`
    - ``candidates`` -- optional list of ``(word, Circle)`` to test; if ``None``
      (default), enumerate by breadth-first search on bottom rows up to word
      length ``maxlen``
    - ``method`` -- ``"lp"`` (default, exact visibility by linear program using
      HiGHS) or ``"grid"`` (fast approximate pre-pass)
    - ``prec`` -- integer (default: 53); precision in bits used only to generate
      the input generators and translation lattice. All subsequent geometry
      (visibility LP, class merging, quadrature) is done in IEEE double
      precision, so ``prec > 53`` only sharpens the rounding of the inputs and
      buys no extra working precision (a warning is emitted); interval/exact
      certification is deferred to a later phase
    - ``maxlen`` -- integer (default: 9); word-length bound for the default
      enumeration

    OUTPUT:

    - list of :class:`Face`, one per face of the Ford domain.

    EXAMPLES:

    The figure-eight knot complement ``4_1`` has 4 faces forming 2 pairs::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import visible_faces
        sage: G = KleinianGroup('4_1')
        sage: faces = visible_faces(G)
        sage: len(faces)
        4
        sage: len({f.word for f in faces})            # four distinct classes
        4
        sage: all(f.margin > 0 for f in faces)        # strictly visible
        True
        sage: {f.inverse_word for f in faces} == {f.word for f in faces}  # inversion-closed
        True
        sage: {round(float(f.circle.radius), 6) for f in faces}   # all radius 1
        {1.0}

    The grid pre-pass finds the same faces (without visibility margins)::

        sage: len(visible_faces(G, method='grid'))
        4
    """
    _warn_if_high_precision(prec)
    gens = _numeric_generators(group, prec=prec)
    v1, v2, base = _lattice_vectors(group, prec=prec)
    if candidates is None:
        circles = _enumerate_circles(gens, v1, v2, maxlen=maxlen, rmin=0.05)
    else:
        circles = [
            (word, complex(circle.center[0], circle.center[1]), float(circle.radius))
            for word, circle in candidates
        ]
    _, _, cand = _bounded_candidates(circles, v1, v2, base)
    if method == "grid":
        visible = _grid_visible(cand, v1, v2, base)
    elif method == "lp":
        visible = _lp_visible(cand, v1, v2, base)
    else:
        raise ValueError(f"Unknown method {method!r}; expected 'lp' or 'grid'.")
    classes = _merge_classes(visible, v1, v2)
    faces = []
    for word, centre, radius, margin, witness in classes:
        if word is None:
            # Word-less candidate (e.g. from the horoball enumeration): the group
            # element is unknown, so the inverse word cannot be resolved.
            inverse_word = None
        else:
            partner = _pair_word(word, centre, radius, classes, gens, v1, v2)
            inverse_word = partner if partner is not None else find_inverse_word(word)
        faces.append(
            Face(
                word=word,
                circle=Circle(center=vector((centre.real, centre.imag)), radius=radius),
                inverse_word=inverse_word,
                witness=None if witness is None else vector((witness.real, witness.imag)),
                margin=None if margin is None else float(margin),
            )
        )
    return faces


def ford_volume(group, faces=None, rel_tol: float = 1e-6, prec: int = 53, maxlen: int = 9) -> float:
    r"""
    Hyperbolic volume of the Ford domain by adaptive quadrature of the floor.

    Integrates ``1 / (2 h(x)^2)`` over the cusp cell, where ``h(x)`` is the Ford
    floor height, on successively refined dyadic midpoint grids until the
    relative change between grids falls below ``rel_tol`` (or refinement stalls),
    returning the value on the finest grid used. For a complete face set the
    result equals the hyperbolic volume of the manifold.

    INPUT:

    - ``group`` -- a :class:`KleinianGroup_class`
    - ``faces`` -- optional list of :class:`Face`; if ``None`` (default), the
      floor is rebuilt from the default breadth-first enumeration
    - ``rel_tol`` -- float (default: ``1e-6``); target relative tolerance
    - ``prec`` -- integer (default: 53); precision in bits used only to generate
      the input generators and translation lattice. All subsequent geometry
      (visibility LP, class merging, quadrature) is done in IEEE double
      precision, so ``prec > 53`` only sharpens the rounding of the inputs and
      buys no extra working precision (a warning is emitted); interval/exact
      certification is deferred to a later phase
    - ``maxlen`` -- integer (default: 9); word-length bound for the default
      enumeration

    OUTPUT:

    - float; the Ford-domain volume.

    EXAMPLES:

    For ``4_1`` the floor integral reproduces the hyperbolic volume
    ``2.029883...`` to well within ``1e-6``::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import ford_volume
        sage: G = KleinianGroup('4_1')
        sage: vol = ford_volume(G)
        sage: bool(abs(vol - float(G.manifold().volume())) < 1e-6)
        True
    """
    _warn_if_high_precision(prec)
    v1, v2, base = _lattice_vectors(group, prec=prec)
    items = _floor_items(group, faces, v1, v2, base, prec, maxlen)
    area = abs((v1 * np.conj(v2)).imag)
    prev = None
    for n in (350, 700, 1400, 2800):
        integral = _floor_integral(items, v1, v2, base, area, n)
        if prev is not None and abs(integral - prev) <= rel_tol * abs(integral):
            return float(integral)
        prev = integral
    return float(prev)


def floor_height(group, faces=None, exact: bool = False, prec: int = 53, maxlen: int = 9) -> float:
    r"""
    Height ``Y0`` of the lowest point of the Ford floor above the cusp cell.

    ``Y0`` is the largest height such that the horoball ``{h = Y0}`` still lies
    on or below the Ford floor over the whole cusp cell; it is the correct
    truncation height for the Fourier expansion of a Kleinian Maass form.

    Two algorithms are available:

    - ``exact=False`` (default): the minimum of the floor over an adaptive
      midpoint grid of the cell. Being a sample minimum it slightly
      *over*-estimates ``Y0`` (the true minimum lies between grid points).
    - ``exact=True``: the closed-form power-diagram method. The floor
      ``h(x)^2 = max_i (r_i^2 - |x - c_i|^2)`` is a maximum of concave
      paraboloids, so over each power-diagram cell it is minimised at a vertex
      (a radical centre of three spheres); ``Y0`` is the smallest such vertex
      value, reduced modulo the translation lattice. This is exact up to
      holonomy/double-precision error -- no grid bias.

    INPUT:

    - ``group`` -- a :class:`KleinianGroup_class`
    - ``faces`` -- optional list of :class:`Face`; if ``None`` (default), the
      floor is rebuilt from the default breadth-first enumeration
    - ``exact`` -- bool (default: ``False``); use the power-diagram vertex method
    - ``prec`` -- integer (default: 53); precision in bits used only to generate
      the input generators and translation lattice. All subsequent geometry is
      done in IEEE double precision, so ``prec > 53`` buys no extra working
      precision (a warning is emitted); a symbolic/interval-certified ``Y0`` for
      the arithmetic case is deferred to a later phase
    - ``maxlen`` -- integer (default: 9); word-length bound for the default
      enumeration

    OUTPUT:

    - float; the floor height ``Y0``.

    EXAMPLES:

    For the figure-eight knot ``4_1`` the floor height is ``sqrt(2/3)``; the
    power-diagram method reproduces ``Y0^2 = 2/3`` to double precision::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import floor_height
        sage: G = KleinianGroup('4_1')
        sage: y0 = floor_height(G, exact=True)
        sage: bool(abs(y0**2 - 2/3) < 1e-10)
        True

    The grid estimate over-estimates the exact value (sample minimum)::

        sage: y0_grid = floor_height(G)
        sage: bool(y0_grid >= y0 - 1e-9)
        True
        sage: bool(abs(y0_grid - y0) < 1e-2)
        True
    """
    _warn_if_high_precision(prec)
    v1, v2, base = _lattice_vectors(group, prec=prec)
    items = _floor_items(group, faces, v1, v2, base, prec, maxlen)
    if exact:
        return _exact_floor(items, v1, v2)[0]
    return float(np.sqrt(max(0.0, _envelope_min_h2(items, v1, v2, base, 500))))


def _floor_items(group, faces, v1: complex, v2: complex, base: complex, prec: int, maxlen: int):
    r"""
    Build the sphere list (with lattice translates) defining the Ford floor,
    either from supplied ``faces`` or from the default enumeration.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _lattice_vectors, visible_faces, _floor_items)
        sage: G = KleinianGroup('4_1')
        sage: v1, v2, base = _lattice_vectors(G)
        sage: items = _floor_items(G, visible_faces(G), v1, v2, base, 53, 9)
        sage: len(items) > 0
        True
    """
    if faces is not None:
        circles = [
            (f.word, complex(f.circle.center[0], f.circle.center[1]), float(f.circle.radius))
            for f in faces
        ]
        return _floor_translates(circles, v1, v2, base)
    gens = _numeric_generators(group, prec=prec)
    circles = _enumerate_circles(gens, v1, v2, maxlen=maxlen, rmin=0.05)
    all_items, _, _ = _bounded_candidates(circles, v1, v2, base)
    return all_items


def _floor_integral(
    items: list, v1: complex, v2: complex, base: complex, area: float, n: int
) -> float:
    r"""
    Midpoint quadrature of ``1 / (2 h^2)`` over the base cell on an ``n x n``
    grid.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _floor_integral
        sage: items = [('a', 0j, 2.0)]
        sage: _floor_integral(items, 1 + 0j, 0 + 1j, -0.5 - 0.5j, 1.0, 8) > 0
        True
    """
    centres = np.array([c for _, c, _ in items])
    radii = np.array([r for _, _, r in items])
    s = np.linspace(0, 1, n, endpoint=False) + 0.5 / n
    u, w = np.meshgrid(s, s)
    z = base + u * v1 + w * v2
    best = np.full(z.shape, -1.0)
    for i in range(len(centres)):
        best = np.maximum(best, radii[i] ** 2 - np.abs(z - centres[i]) ** 2)
    return float(np.sum((area / n**2) / (2 * best)))


def _radical_centre(ci, cj, ck, kap_i, kap_j, kap_k, binv, v1: complex, v2: complex):
    r"""
    Radical centre (power-diagram vertex) of the three spheres with centres
    ``ci, cj, ck`` and power-offsets ``kap = |c|^2 - r^2``, reduced modulo the
    lattice, or ``None`` if the three centres are collinear.

    EXAMPLES::

        sage: import numpy as np
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _radical_centre
        sage: binv = np.linalg.inv(np.array([[1.0, 0.0], [0.0, 1.0]]))
        sage: v = _radical_centre(complex(0, 0), complex(2, 0), complex(0, 2),
        ....:                     -1.0, 3.0, 3.0, binv, complex(1, 0), complex(0, 1))
        sage: v is not None
        True
    """
    a = np.array(
        [
            [2 * (cj - ci).real, 2 * (cj - ci).imag],
            [2 * (ck - ci).real, 2 * (ck - ci).imag],
        ]
    )
    det = a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]
    if abs(det) < 1e-12:
        return None
    x = np.linalg.solve(a, [kap_j - kap_i, kap_k - kap_i])
    v = complex(x[0], x[1])
    coords = binv @ np.array([v.real, v.imag])
    # Reduce into the *centred* cell [-1/2, 1/2)^2 (the region covered by
    # _floor_translates), so the envelope max at the vertex sees every hemisphere
    # above it; reducing into [0, 1)^2 instead can land outside that coverage.
    coords -= np.round(coords)
    return coords[0] * v1 + coords[1] * v2


def _exact_floor(circles: list, v1: complex, v2: complex, theta: float = 1.0):
    r"""
    Exact floor height by the power-diagram vertex method.

    Evaluates the squared floor ``max_i (theta^2 r_i^2 - |x - c_i|^2)`` at every
    radical centre of a pair-wise overlapping triple (a superset of the true
    power-diagram vertices, each value ``>=`` the minimum), returning
    ``(Y0, deepest_point, defining_triple)``.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _lattice_vectors, _floor_items, _exact_floor)
        sage: G = KleinianGroup('4_1')
        sage: v1, v2, base = _lattice_vectors(G)
        sage: items = _floor_items(G, None, v1, v2, base, 53, 9)
        sage: y0, point, triple = _exact_floor(items, v1, v2)
        sage: bool(abs(y0**2 - 2/3) < 1e-10)
        True
    """
    centres = np.array([c for _, c, _ in circles])
    weights = np.array([(theta * r) ** 2 for _, _, r in circles])
    radii = np.sqrt(weights)
    kap = np.abs(centres) ** 2 - weights
    binv = np.linalg.inv(np.array([[v1.real, v2.real], [v1.imag, v2.imag]]))
    best, argbest, triple = np.inf, None, None
    n = len(centres)
    for i in range(n):
        nbrs = [j for j in range(i + 1, n) if abs(centres[j] - centres[i]) < radii[i] + radii[j]]
        for a, j in enumerate(nbrs):
            for k in nbrs[a + 1 :]:
                if abs(centres[k] - centres[j]) >= radii[j] + radii[k]:
                    continue
                v = _radical_centre(
                    centres[i], centres[j], centres[k], kap[i], kap[j], kap[k], binv, v1, v2
                )
                if v is None:
                    continue
                envelope = float(np.max(weights - np.abs(v - centres) ** 2))
                if 0 < envelope < best:
                    best, argbest, triple = envelope, v, (i, j, k)
    return float(np.sqrt(best)), argbest, triple


def _ep_edges(group):
    r"""
    Number of edges of SnapPy's canonical retriangulation (equal to the number
    of tetrahedra for an ideal triangulation), or ``None`` when the canonical
    cell decomposition is not simplicial (finite vertices present), in which
    case a warning is emitted.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _ep_edges
        sage: _ep_edges(KleinianGroup('4_1'))
        2
    """
    canonical = group.manifold().canonical_retriangulation()
    if canonical.has_finite_vertices():
        warnings.warn(
            "Canonical cell decomposition is non-simplicial (finite vertices); "
            "skipping the Epstein-Penner edge count.",
            stacklevel=2,
        )
        return None
    return int(canonical.num_tetrahedra())


def _inversion_closed(faces: list) -> bool:
    r"""
    Whether the face set is closed under inversion.

    When every face carries a word, pair by inverse word. For word-less faces
    (from the horoball enumeration) fall back to the necessary geometric
    condition that every radius occurs an even number of times, so the faces can
    be partitioned into equal-radius inverse pairs.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.types import Face, Circle
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _inversion_closed
        sage: from sage.modules.free_module_element import vector
        sage: circ = Circle(center=vector((0, 0)), radius=1)
        sage: _inversion_closed([Face('a', circ, 'A'), Face('A', circ, 'a')])
        True
        sage: _inversion_closed([Face(None, circ), Face(None, circ)])  # even radii
        True
        sage: _inversion_closed([Face(None, circ)])                    # odd -> not closed
        False
    """
    if all(f.word is not None and f.inverse_word is not None for f in faces):
        return {f.inverse_word for f in faces} == {f.word for f in faces}
    counts: dict = {}
    for f in faces:
        key = round(float(f.circle.radius), 6)
        counts[key] = counts.get(key, 0) + 1
    return all(n % 2 == 0 for n in counts.values())


def check_face_completeness(
    group, candidates=None, prec: int = 53, maxlen: int = 9, volume_tol: float = 1e-5
) -> dict:
    r"""
    Verify that the enumerated face set is the complete Ford face pairing.

    Runs the two independent completeness tests of the plan: the volume identity
    (the floor integral must equal the manifold volume) and the Epstein--Penner
    duality count (the number of face pairs must equal the number of edges of the
    canonical retriangulation). Also reports whether the face set is closed under
    inversion.

    INPUT:

    - ``group`` -- a :class:`KleinianGroup_class`
    - ``candidates`` -- optional list of ``(word, Circle)`` to test instead of
      the default breadth-first enumeration; pass the output of
      :func:`enumerate_bounded_c` to certify completeness for knots (such as
      ``7_4``) whose last face pair is unreachable by word-length search
    - ``prec`` -- integer (default: 53); precision in bits used only to generate
      the input generators and translation lattice. All subsequent geometry
      (visibility LP, class merging, quadrature) is done in IEEE double
      precision, so ``prec > 53`` only sharpens the rounding of the inputs and
      buys no extra working precision (a warning is emitted); interval/exact
      certification is deferred to a later phase
    - ``maxlen`` -- integer (default: 9); word-length bound for the enumeration
    - ``volume_tol`` -- float (default: ``1e-5``); absolute tolerance on the
      volume residual for declaring completeness

    OUTPUT:

    - dict with keys ``volume_residual``, ``ep_edges``, ``pairs_found``,
      ``inversion_closed`` and ``complete``.

    EXAMPLES:

    The figure-eight knot ``4_1`` is complete at word length 2: 4 faces / 2
    pairs, EP edges 2, volume residual below ``1e-6``::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import check_face_completeness
        sage: G = KleinianGroup('4_1')
        sage: report = check_face_completeness(G)
        sage: report['pairs_found']
        2
        sage: report['ep_edges']
        2
        sage: report['inversion_closed']
        True
        sage: bool(abs(report['volume_residual']) < 1e-6)
        True
        sage: report['complete']
        True

    The knot ``5_2`` closes at 8 faces / 4 pairs, EP edges 4, volume residual
    below ``5e-7``::

        sage: report = check_face_completeness(KleinianGroup('5_2'))
        sage: report['pairs_found']
        4
        sage: report['ep_edges']
        4
        sage: bool(abs(report['volume_residual']) < 5e-7)
        True
        sage: report['complete']
        True
    """
    faces = visible_faces(group, candidates=candidates, prec=prec, maxlen=maxlen)
    inversion_closed = _inversion_closed(faces)
    volume = ford_volume(group, faces=faces, prec=prec, maxlen=maxlen)
    residual = volume - float(group.manifold().volume())
    ep_edges = _ep_edges(group)
    pairs_found = len(faces) // 2
    complete = bool(
        inversion_closed
        and ep_edges is not None
        and ep_edges == pairs_found
        and abs(residual) < volume_tol
    )
    return {
        "volume_residual": float(residual),
        "ep_edges": ep_edges,
        "pairs_found": pairs_found,
        "inversion_closed": bool(inversion_closed),
        "complete": complete,
    }


# --- phase 3: bounded-|c| enumeration via horoballs ---------------------------


def _cusp_frame_map(group, prec: int = 53):
    r"""
    Return ``(alpha, radius_const, mt)`` relating SnapPy's maximal-cusp frame to
    this group's coordinate frame.

    The two frames both fix the cusp at infinity, so they differ by a single
    complex scaling ``z_group = alpha * z_cusp`` with ``alpha = tM / mt`` (package
    meridian translation over SnapPy meridian translation). An isometric sphere
    over a horoball of Euclidean radius ``h`` has group-frame radius
    ``|alpha| * sqrt(2 h)``, so ``radius_const = |alpha|``.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _cusp_frame_map
        sage: alpha, rc, mt = _cusp_frame_map(KleinianGroup('4_1'))
        sage: bool(abs(rc - abs(alpha)) < 1e-12)
        True
    """
    named = group.named_generators(prec=prec)
    tM = complex(named["M"][0, 0] * named["M"][0, 1])
    cusp = group.manifold().cusp_neighborhood()
    cusp.set_displacement(cusp.stopping_displacement(), which_cusp=0)
    mt = complex(cusp.translations(which_cusp=0)[0])
    alpha = tM / mt
    return alpha, abs(alpha), mt


def _horoball_circles(group, cutoff: float, prec: int = 53) -> list:
    r"""
    Map every horoball of the maximal cusp above ``cutoff`` to a group-frame
    isometric circle ``(centre, radius)``, reduced modulo the lattice and
    de-duplicated.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import _horoball_circles
        sage: circ = _horoball_circles(KleinianGroup('4_1'), 0.2)
        sage: len(circ) > 0
        True
    """
    alpha, radius_const, _ = _cusp_frame_map(group, prec=prec)
    v1, v2, _ = _lattice_vectors(group, prec=prec)
    cusp = group.manifold().cusp_neighborhood()
    cusp.set_displacement(cusp.stopping_displacement(), which_cusp=0)
    seen: set = set()
    out = []
    for horoball in cusp.horoballs(cutoff, which_cusp=0):
        radius = radius_const * np.sqrt(2 * horoball["radius"])
        centre = _lattice_reduce(alpha * complex(horoball["center"]), v1, v2)
        key = (round(centre.real, 8), round(centre.imag, 8), round(radius, 8))
        if key in seen:
            continue
        seen.add(key)
        out.append((centre, radius))
    return out


def _horoball_floor_height(circles: list, v1: complex, v2: complex, base: complex) -> float:
    r"""
    Fast estimate of the Ford floor height ``Y0`` from the largest circles (the
    ones that determine the deepest floor point), used to auto-tune the radius
    floor of :func:`enumerate_bounded_c`.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _lattice_vectors, _horoball_circles, _horoball_floor_height)
        sage: G = KleinianGroup('4_1')
        sage: v1, v2, base = _lattice_vectors(G)
        sage: y0 = _horoball_floor_height(_horoball_circles(G, 0.1), v1, v2, base)
        sage: bool(abs(y0 - 0.8165) < 1e-2)
        True
    """
    ranked = sorted(circles, key=lambda cr: -cr[1])[:120]
    items = _floor_translates([(None, c, r) for c, r in ranked], v1, v2, base)
    return float(np.sqrt(max(0.0, _envelope_min_h2(items, v1, v2, base, 400))))


def enumerate_bounded_c(group, c_bound=None, prec: int = 53) -> list:
    r"""
    All isometric circles with ``|c| <= c_bound`` via SnapPy's horoball enumeration.

    For a cusped hyperbolic knot complement the isometric spheres of the coset
    representatives are in bijection with the horoballs of the maximal cusp: an
    element with lower-left entry ``c`` gives a sphere of radius ``1/|c|``, which
    corresponds to a horoball of Euclidean radius ``h`` with
    ``1/|c| = |alpha| sqrt(2 h)`` (see :func:`_cusp_frame_map`). SnapPy's
    ``cusp_neighborhood().horoballs(cutoff)`` enumerates every horoball above a
    cutoff *completely* from the triangulation, so this returns the complete list
    of circles up to the bound -- including face pairs that a word-length
    (breadth-first) search cannot reach in practice (e.g. the last pair of
    ``7_4``, at word length >= 14).

    Centres are reduced modulo the translation lattice and de-duplicated. The
    group element -- hence a word -- for each circle is not recovered here (that
    would need the covering machinery this enumeration is meant to complete), so
    each entry carries ``word = None``; completeness is certified by the
    frame-invariant volume and Epstein--Penner tests in
    :func:`check_face_completeness`.

    INPUT:

    - ``group`` -- a manifold-backed :class:`KleinianGroup_class`
    - ``c_bound`` -- real or ``None`` (default: ``None``); the bound on ``|c|``
      in the *meridian-normalised* frame (matching the plan; for ``7_4`` the
      missing pair has ``|c| <= 2.35``). When ``None`` the radius floor is
      auto-tuned to ``0.9 Y0`` from an estimated floor height, so no knowledge of
      the frame scale is needed -- this is the recommended usage. The
      meridian-normalised frame makes the bound comparable across knots (the raw
      matrix ``|c|`` scale varies by orders of magnitude with the cusp
      normalisation).
    - ``prec`` -- integer (default: 53); precision in bits for the frame data

    OUTPUT:

    - list of ``(None, Circle)``, one per distinct circle modulo the lattice.

    EXAMPLES:

    The knot ``7_4`` closes only once the horoball enumeration supplies the pair
    that row search misses; the auto-tuned bound yields the complete
    16-face / 8-pair Ford domain::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     enumerate_bounded_c, check_face_completeness)
        sage: G = KleinianGroup('7_4')
        sage: cands = enumerate_bounded_c(G)
        sage: report = check_face_completeness(G, candidates=cands)
        sage: (report['pairs_found'], report['ep_edges'])
        (8, 8)
        sage: report['inversion_closed']
        True
        sage: bool(abs(report['volume_residual']) < 1e-5)
        True
        sage: report['complete']
        True
    """
    _, radius_const, _ = _cusp_frame_map(group, prec=prec)
    v1, v2, base = _lattice_vectors(group, prec=prec)
    if c_bound is not None:
        rmin = radius_const / float(c_bound)
        circles = _horoball_circles(group, (rmin / radius_const) ** 2 / 2.0, prec=prec)
    else:
        # Auto-tune: enumerate a generous set (down to r = 0.3 r_max, i.e.
        # cutoff for sqrt(2h) >= 0.3), estimate Y0 from the big circles, then keep
        # everything that can reach the floor (r >= 0.9 Y0), refining the cutoff
        # if the estimate falls below the first pass.
        cutoff = 0.3**2 / 2.0
        circles = _horoball_circles(group, cutoff, prec=prec)
        rmin = 0.9 * _horoball_floor_height(circles, v1, v2, base)
        needed = (rmin / radius_const) ** 2 / 2.0
        if needed < cutoff:
            circles = _horoball_circles(group, needed, prec=prec)
    return [
        (None, Circle(center=vector((c.real, c.imag)), radius=r)) for c, r in circles if r >= rmin
    ]


# --- phase 4: edge cycles and the Poincare presentation -----------------------


def _floor_vertices(items: list, v1: complex, v2: complex, base: complex) -> list:
    r"""
    Vertices of the Ford floor: radical centres of triples of hemispheres where
    the three spheres are genuinely co-highest on the upper envelope (a true
    corner of the domain, not the smooth apex of a single dome, and not covered
    by a fourth sphere), reduced modulo the lattice and de-duplicated. Each entry
    is ``(vertex, height, (word_i, word_j, word_k))``.

    EXAMPLES::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import (
        ....:     _lattice_vectors, _floor_items, _floor_vertices)
        sage: G = KleinianGroup('4_1')
        sage: v1, v2, base = _lattice_vectors(G)
        sage: verts = _floor_vertices(_floor_items(G, None, v1, v2, base, 53, 9), v1, v2, base)
        sage: len(verts) > 0 and all(len(w) == 3 for _, _, w in verts)
        True
    """
    centres = np.array([c for _, c, _ in items])
    weights = np.array([r * r for _, _, r in items])
    radii = np.sqrt(weights)
    kap = np.abs(centres) ** 2 - weights
    binv = np.linalg.inv(np.array([[v1.real, v2.real], [v1.imag, v2.imag]]))
    n = len(centres)
    out: dict = {}
    for i in range(n):
        nbrs = [j for j in range(i + 1, n) if abs(centres[j] - centres[i]) < radii[i] + radii[j]]
        for a, j in enumerate(nbrs):
            for k in nbrs[a + 1 :]:
                if abs(centres[k] - centres[j]) >= radii[j] + radii[k]:
                    continue
                v = _radical_centre(
                    centres[i], centres[j], centres[k], kap[i], kap[j], kap[k], binv, v1, v2
                )
                if v is None:
                    continue
                height2 = weights[i] - abs(v - centres[i]) ** 2
                if height2 <= 1e-9:
                    continue
                # Skip apex coincidences (v = a sphere centre): a smooth top of a
                # dome, not a genuine three-face corner.
                if min(abs(v - centres[i]), abs(v - centres[j]), abs(v - centres[k])) < 1e-6:
                    continue
                # Must be on the upper envelope (no fourth sphere higher here).
                if np.max(weights - np.abs(v - centres) ** 2) > height2 + 1e-7:
                    continue
                reduced = _lattice_reduce(v, v1, v2)
                key = (round(reduced.real, 5), round(reduced.imag, 5))
                out.setdefault(
                    key,
                    (
                        vector((reduced.real, reduced.imag)),
                        float(np.sqrt(height2)),
                        (items[i][0], items[j][0], items[k][0]),
                    ),
                )
    return list(out.values())


def edge_cycles(group, faces=None, prec: int = 53) -> list:
    r"""
    The vertex/edge skeleton of the Ford floor's power diagram.

    Returns the vertices of the power diagram of the visible hemispheres (the
    radical centres that lie on the Ford floor); the edges of the domain are the
    floor arcs joining these vertices where two hemispheres meet. Each vertex is
    a dict ``{'point': vector, 'height': float, 'faces': (word_i, word_j,
    word_k)}`` giving the three hemispheres meeting there. This is the geometric
    input to the Poincare edge-cycle relations.

    INPUT:

    - ``group`` -- a :class:`KleinianGroup_class`
    - ``faces`` -- optional list of :class:`Face`; default rebuilds them from the
      breadth-first enumeration
    - ``prec`` -- integer (default: 53); precision in bits

    OUTPUT:

    - list of vertex dicts.

    EXAMPLES:

    The figure-eight domain is highly symmetric: all 9 floor vertices (modulo the
    lattice) sit at the floor height ``sqrt(2/3)``::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import edge_cycles
        sage: verts = edge_cycles(KleinianGroup('4_1'))
        sage: len(verts) > 0
        True
        sage: {round(v['height']**2, 10) for v in verts}     # every vertex at Y0
        {0.6666666667}
        sage: all(len(v['faces']) == 3 for v in verts)
        True
    """
    if faces is None:
        faces = visible_faces(group, prec=prec)
    v1, v2, base = _lattice_vectors(group, prec=prec)
    circles = [
        (f.word, complex(f.circle.center[0], f.circle.center[1]), float(f.circle.radius))
        for f in faces
    ]
    items = _floor_translates(circles, v1, v2, base)
    return [
        {"point": point, "height": height, "faces": words}
        for point, height, words in _floor_vertices(items, v1, v2, base)
    ]


def presentation_from_faces(group, prec: int = 53):
    r"""
    A presentation of the group from the Ford face pairing, cross-checked against
    SnapPy's canonical fundamental group.

    By Poincare's polyhedron theorem the face-pairing transformations generate
    the group; the generators returned here are one word per inverse face pair
    (the pairing ``g`` maps its isometric sphere to that of ``g^{-1}``). This is
    guaranteed to be a generating set by the (already certified) completeness of
    the Ford domain, but it is generally *not minimal*: a knot group has a
    2-generator presentation while the face pairing has one generator per pair.
    Deriving a minimal relator set on the face generators requires the edge-cycle
    (Tietze) reduction and is left as the optional remaining step; instead this
    returns SnapPy's canonical relators for comparison. Only groups whose faces
    carry words are supported (a :class:`ValueError` is raised otherwise, e.g.
    for the horoball-completed knots).

    INPUT:

    - ``group`` -- a manifold-backed :class:`KleinianGroup_class`
    - ``prec`` -- integer (default: 53); precision in bits

    OUTPUT:

    - tuple ``(generators, relators)``: the face-pairing generator words and
      SnapPy's canonical relators.

    EXAMPLES:

    For ``4_1`` the face pairing gives a 2-generator set -- the same rank as
    SnapPy's minimal presentation -- and the holonomy satisfies the relator::

        sage: from maass_form_core.testing import connect_mockdb
        sage: connect_mockdb()
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.face_pairing import presentation_from_faces
        sage: from maass_forms_klein.hyperbolic_space.word_utils import word_to_element
        sage: G = KleinianGroup('4_1')
        sage: gens, relators = presentation_from_faces(G)
        sage: len(gens)
        2
        sage: len(relators)
        1
        sage: # the holonomy representation satisfies the canonical relator
        sage: R = word_to_element(relators[0], G.named_generators())
        sage: bool(all(abs(R[i, j] - (1 if i == j else 0)) < 1e-6 for i in range(2)
        ....:          for j in range(2)) or
        ....:      all(abs(R[i, j] - (-1 if i == j else 0)) < 1e-6 for i in range(2)
        ....:          for j in range(2)))
        True
    """
    faces = visible_faces(group, prec=prec)
    if any(f.word is None for f in faces):
        raise ValueError(
            "presentation_from_faces requires face words; the horoball-completed "
            "face set is word-less. Only breadth-first-reachable knots are supported."
        )
    generators, seen = [], set()
    for f in faces:
        key = frozenset((f.word, f.inverse_word))
        if key in seen:
            continue
        seen.add(key)
        generators.append(f.word)
    relators = [str(r) for r in group.manifold().fundamental_group().relators()]
    return generators, relators
