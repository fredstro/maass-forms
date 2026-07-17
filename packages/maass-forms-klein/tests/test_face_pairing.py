"""
Regression tests for Ford-domain face pairings (phase 1).

Covers the regression table of ``FACE_PAIRING_PLAN.md`` for the knots that close
in phase 1 (4_1, 5_2, 6_1), the inversion-closedness and LP-margin-positivity
invariants, the Epstein--Penner non-simplicial guard path, and a property test
that an incomplete face set is detected by the volume identity. The 7_4 row is
deferred to phase 3 (bounded-|c| horoball enumeration).

Numerical tolerances live in the ``REGRESSION`` table next to the expected
counts, not as inline magic numbers.
"""

import random

import pytest

from maass_form_core.testing import connect_mockdb

connect_mockdb()

from maass_forms_klein.hyperbolic_space.face_pairing import (  # noqa: E402
    _ep_edges,
    check_face_completeness,
    edge_cycles,
    enumerate_bounded_c,
    floor_height,
    ford_volume,
    presentation_from_faces,
    visible_faces,
)
from maass_forms_klein.hyperbolic_space.word_utils import word_to_element  # noqa: E402
from maass_forms_klein.hyperbolic_space.kleinian_group import (  # noqa: E402
    KleinianGroup,
    _face_to_facedb,
    _facedb_to_face,
)

# Expected face data with absolute tolerances stored next to the regression data.
# ``volume_tol`` bounds the floor-volume residual; ``margin_min`` is a lower bound
# on the LP visibility margins (a stored floor -- the margins are strictly
# positive and dwarf the ~1e-12 double-precision/holonomy error). Margins are
# frame-dependent, so in the package frame they are smaller than the plan's
# meridian-frame ">= 1e-2"; the measured minima are 4_1:1.0, 5_2:2.1e-4, 6_1:1.0e-4.
REGRESSION = {
    "4_1": {"faces": 4, "pairs": 2, "ep_edges": 2, "volume_tol": 1e-6, "margin_min": 0.5},
    "5_2": {"faces": 8, "pairs": 4, "ep_edges": 4, "volume_tol": 5e-7, "margin_min": 1e-4},
    "6_1": {"faces": 12, "pairs": 6, "ep_edges": 6, "volume_tol": 4e-7, "margin_min": 5e-5},
}


@pytest.mark.parametrize(
    "knot",
    [
        "4_1",
        pytest.param("5_2", marks=pytest.mark.slow),
        pytest.param("6_1", marks=pytest.mark.slow),
    ],
)
def test_regression_row(knot):
    """The face count, pairing, EP edges and volume residual match the table."""
    expected = REGRESSION[knot]
    group = KleinianGroup(knot)

    faces = visible_faces(group)
    assert len(faces) == expected["faces"]
    # LP visibility margins are strictly positive, above the stored floor.
    assert min(float(face.margin) for face in faces) >= expected["margin_min"]
    # The face set is closed under inversion.
    assert {face.word for face in faces} == {face.inverse_word for face in faces}

    report = check_face_completeness(group)
    assert report["pairs_found"] == expected["pairs"]
    assert report["ep_edges"] == expected["ep_edges"]
    assert report["inversion_closed"] is True
    assert abs(report["volume_residual"]) < expected["volume_tol"]
    assert report["complete"] is True


@pytest.mark.slow
def test_incomplete_face_set_is_detected():
    """Random proper subsets of the 6_1 face set fail the volume test."""
    knot = "6_1"
    tol = REGRESSION[knot]["volume_tol"]
    group = KleinianGroup(knot)
    faces = visible_faces(group)
    true_volume = float(group.manifold().volume())

    # The complete face set reproduces the hyperbolic volume.
    assert abs(ford_volume(group, faces=faces) - true_volume) < tol

    # Dropping any genuine face opens a gap in the floor: the exposed region
    # drops onto lower hemispheres, so the volume integral grows past tolerance.
    # Several seeded random subsets (deterministic) exercise this detector.
    rng = random.Random(20240716)
    for _ in range(8):
        drop = set(rng.sample(range(len(faces)), rng.randint(1, 3)))
        subset = [face for i, face in enumerate(faces) if i not in drop]
        residual = ford_volume(group, faces=subset) - true_volume
        assert residual > tol, f"undetected incomplete subset dropping {sorted(drop)}"


class _FakeCanonical:
    """Stand-in for a SnapPy canonical retriangulation."""

    def __init__(self, finite_vertices: bool, num_tetrahedra: int):
        self._finite = finite_vertices
        self._num_tetrahedra = num_tetrahedra

    def has_finite_vertices(self) -> bool:
        return self._finite

    def num_tetrahedra(self) -> int:
        return self._num_tetrahedra


class _FakeManifold:
    def __init__(self, canonical: _FakeCanonical):
        self._canonical = canonical

    def canonical_retriangulation(self) -> _FakeCanonical:
        return self._canonical


class _FakeGroup:
    def __init__(self, manifold: _FakeManifold):
        self._manifold = manifold

    def manifold(self) -> _FakeManifold:
        return self._manifold


def test_ep_edges_counts_tetrahedra_when_simplicial():
    """A simplicial canonical decomposition reports its tetrahedron count."""
    group = _FakeGroup(_FakeManifold(_FakeCanonical(finite_vertices=False, num_tetrahedra=6)))
    assert _ep_edges(group) == 6


def test_ep_edges_guards_non_simplicial_decomposition():
    """A non-simplicial (finite-vertex) decomposition returns None and warns."""
    group = _FakeGroup(_FakeManifold(_FakeCanonical(finite_vertices=True, num_tetrahedra=99)))
    with pytest.warns(UserWarning, match="non-simplicial"):
        assert _ep_edges(group) is None


# --- phase 2: exact floor height, group methods, database round-trip ----------


def test_exact_floor_height_is_sqrt_two_thirds():
    """4_1 has floor height sqrt(2/3); the power-diagram method is exact."""
    group = KleinianGroup("4_1")
    y0 = floor_height(group, exact=True)
    assert abs(y0**2 - 2 / 3) < 1e-10
    # the grid estimate is a sample minimum, so it over-estimates Y0.
    assert floor_height(group) >= y0 - 1e-9


def test_group_ford_methods():
    """The cached KleinianGroup Ford methods agree with the regression row."""
    group = KleinianGroup("4_1")
    assert abs(group.ford_floor_height() ** 2 - 2 / 3) < 1e-10
    assert abs(group.ford_volume() - float(group.manifold().volume())) < 1e-6
    report = group.check_ford_complete()
    assert (report["pairs_found"], report["ep_edges"], report["complete"]) == (2, 2, True)


@pytest.mark.slow
def test_find_max_y_uses_certified_floor_height():
    """find_max_y starts from the certified Ford floor height for a knot space."""
    from maass_forms_klein.modform.compute_coefficients import _certified_starting_y, find_max_y
    from maass_forms_klein.modform.kmaass_space import KleinianMaassFormSpace

    space = KleinianMaassFormSpace("4_1")
    y0 = float(space.group().ford_floor_height())

    # the certified start is just below the exact floor height
    assert abs(_certified_starting_y(space) - 0.995 * y0) < 1e-9

    # with no explicit Y, find_max_y uses it: a valid height (pullback raises
    # every point) that beats the historical 0.73 default for this knot.
    max_y = float(find_max_y(space, M=1, Q=2))
    assert 0.73 < max_y <= y0 + 1e-9


def test_face_facedb_roundtrip():
    """A face survives conversion to a FaceDB document and back."""
    face = visible_faces(KleinianGroup("4_1"))[0]
    restored = _facedb_to_face(_face_to_facedb(face))
    assert restored.word == face.word
    assert restored.inverse_word == face.inverse_word
    assert abs(float(restored.circle.radius) - float(face.circle.radius)) < 1e-12
    assert abs(float(restored.circle.center[0]) - float(face.circle.center[0])) < 1e-12


def test_ford_domain_db_roundtrip():
    """ford_faces seeds a FordDomainDB document that reconstructs the faces."""
    group = KleinianGroup("4_1")
    faces = group.ford_faces()
    doc = group._load_ford_domain()
    assert doc is not None
    assert len(doc.faces) == len(faces)
    assert {f.word for f in doc.faces} == {f.word for f in faces}
    assert {f.partner_word for f in doc.faces} == {f.inverse_word for f in faces}


def test_face_persists_through_mock_database():
    """A single Face -> FaceDB survives a real save/query cycle in the mock DB."""
    from maass_forms_klein.database.models import FordDomainDB

    face = KleinianGroup("4_1").ford_faces()[0]
    label = "test_face_pairing_single"
    FordDomainDB.objects(label=label).delete()
    try:
        FordDomainDB(label=label, faces=[_face_to_facedb(face)], y0=0.816, complete=True).save()
        loaded = FordDomainDB.objects(label=label).first()
        assert loaded is not None
        assert len(loaded.faces) == 1
        # reconstruct the Face from the document actually read back from the DB
        restored = _facedb_to_face(loaded.faces[0])
        assert restored.word == face.word
        assert restored.inverse_word == face.inverse_word
        assert abs(float(restored.circle.radius) - float(face.circle.radius)) < 1e-12
        assert abs(float(restored.circle.center[0]) - float(face.circle.center[0])) < 1e-12
        assert abs(float(restored.circle.center[1]) - float(face.circle.center[1])) < 1e-12
        assert loaded.y0 == 0.816
        assert loaded.complete is True
    finally:
        FordDomainDB.objects(label=label).delete()


def test_full_face_pairing_persists_through_mock_database():
    """The whole 4_1 face pairing round-trips through a saved FordDomainDB."""
    from maass_forms_klein.database.models import FordDomainDB

    faces = KleinianGroup("4_1").ford_faces()
    label = "test_face_pairing_full"
    FordDomainDB.objects(label=label).delete()
    try:
        FordDomainDB(label=label, faces=[_face_to_facedb(f) for f in faces]).save()
        loaded = FordDomainDB.objects(label=label).first()
        restored = [_facedb_to_face(f) for f in loaded.faces]
        assert len(restored) == len(faces)
        assert {f.word for f in restored} == {f.word for f in faces}
        assert {f.inverse_word for f in restored} == {f.inverse_word for f in faces}
        # radii and centres survive serialisation for every face
        by_word = {f.word: f for f in faces}
        for r in restored:
            original = by_word[r.word]
            assert abs(float(r.circle.radius) - float(original.circle.radius)) < 1e-12
            assert abs(float(r.circle.center[0]) - float(original.circle.center[0])) < 1e-12
            assert abs(float(r.circle.center[1]) - float(original.circle.center[1])) < 1e-12
    finally:
        FordDomainDB.objects(label=label).delete()


def test_ford_domain_db_unique_label_and_delete():
    """FordDomainDB enforces one document per label and delete removes it."""
    from maass_forms_klein.database.models import FordDomainDB

    label = "test_face_pairing_unique"
    FordDomainDB.objects(label=label).delete()
    try:
        FordDomainDB(label=label, faces=[]).save()
        assert FordDomainDB.objects(label=label).count() == 1
    finally:
        FordDomainDB.objects(label=label).delete()
    assert FordDomainDB.objects(label=label).count() == 0


@pytest.mark.slow
def test_seven_four_incomplete_in_phase_one():
    """
    7_4 needs the phase-3 horoball enumeration to close its last pair (the
    missing pair only appears at word length >= 14); with the breadth-first
    search it is (correctly) reported incomplete. A short enumeration bound is
    used since no word length reachable in phase 1 closes the domain.
    """
    report = check_face_completeness(KleinianGroup("7_4"), maxlen=4)
    assert report["complete"] is False
    # some faces are found, but fewer than the complete 8 pairs (16 faces).
    assert 0 < report["pairs_found"] < 8


# --- phase 3: bounded-|c| horoball enumeration --------------------------------


def test_enumerate_bounded_c_closes_seven_four():
    """The horoball enumeration completes 7_4 to 16 faces / 8 pairs."""
    group = KleinianGroup("7_4")
    candidates = enumerate_bounded_c(group)  # auto-tuned bound
    # every candidate is word-less (the group element is not recovered)
    assert all(word is None for word, _ in candidates)
    report = check_face_completeness(group, candidates=candidates)
    assert report["pairs_found"] == 8
    assert report["ep_edges"] == 8
    assert report["inversion_closed"] is True
    assert abs(report["volume_residual"]) < 1e-5
    assert report["complete"] is True


def test_group_ford_faces_certified_closes_seven_four():
    """KleinianGroup.ford_faces / check_ford_complete close 7_4 via horoballs."""
    group = KleinianGroup("7_4")
    faces = group.ford_faces(certified=True)
    assert len(faces) == 16
    report = group.check_ford_complete(certified=True)
    assert (report["pairs_found"], report["ep_edges"], report["complete"]) == (8, 8, True)


@pytest.mark.parametrize("bad", [0, -1.0])
def test_enumerate_bounded_c_rejects_nonpositive_bound(bad):
    """A non-positive c_bound is rejected (would divide by zero / keep all)."""
    with pytest.raises(ValueError, match="c_bound must be positive"):
        enumerate_bounded_c(KleinianGroup("4_1"), c_bound=bad)


# 8-crossing knots: EP edge count and completeness via the horoball pipeline.
# 8_1 stresses trap 8 (sliver faces -> LP counts), 8_2 stresses trap 7 (skew
# cusp tau = 10.98 - 2.99i -> Gauss-reduced basis + adaptive translate rings).
EIGHT_CROSSING = {"8_1": 10, "8_2": 10, "8_3": 10}


@pytest.mark.slow
@pytest.mark.parametrize("knot", sorted(EIGHT_CROSSING))
def test_horoball_closes_eight_crossing(knot):
    """Skew and sliver 8-crossing knots close via the robustness-pass pipeline."""
    group = KleinianGroup(knot)
    report = group.check_ford_complete(certified=True)
    assert report["ep_edges"] == EIGHT_CROSSING[knot]
    assert report["pairs_found"] == EIGHT_CROSSING[knot]
    assert report["inversion_closed"] is True
    assert abs(report["volume_residual"]) < 2e-5
    assert report["complete"] is True


# --- phase 4: edge cycles and presentation ------------------------------------


def _is_pm_identity(matrix) -> bool:
    entries = [complex(matrix[i, j]) for i in range(2) for j in range(2)]
    identity = [1, 0, 0, 1]
    return all(abs(e - v) < 1e-6 for e, v in zip(entries, identity, strict=True)) or all(
        abs(e + v) < 1e-6 for e, v in zip(entries, identity, strict=True)
    )


def test_edge_cycles_figure_eight():
    """4_1's floor vertices are all at Y0 = sqrt(2/3) (its symmetry)."""
    verts = edge_cycles(KleinianGroup("4_1"))
    assert len(verts) > 0
    assert all(len(v["faces"]) == 3 for v in verts)
    assert {round(v["height"] ** 2, 8) for v in verts} == {round(2 / 3, 8)}


@pytest.mark.parametrize(
    "knot,ngens",
    [("4_1", 2), pytest.param("5_2", 4, marks=pytest.mark.slow)],
)
def test_presentation_from_faces(knot, ngens):
    """Face-pairing generators + the holonomy satisfies the canonical relators."""
    group = KleinianGroup(knot)
    generators, relators = presentation_from_faces(group)
    assert len(generators) == ngens
    named = group.named_generators()
    for relator in relators:
        assert _is_pm_identity(word_to_element(relator, named))
    # the face-pairing rank matches SnapPy for 4_1 (both minimal)
    if knot == "4_1":
        assert group.manifold().fundamental_group().num_generators() == ngens


def test_horoball_faces_are_wordless():
    """The horoball-completed face set carries no words (why a presentation of
    those knots needs the deferred word-recovery step)."""
    group = KleinianGroup("7_4")
    faces = visible_faces(group, candidates=enumerate_bounded_c(group))
    assert all(f.word is None for f in faces)
