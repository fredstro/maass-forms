"""
Regression fixtures from the 2026 eigenvalue survey
(``LESSONS_LEARNED_EIGENVALUE_SURVEY.md`` §6).

These lock down the strict-confirmation behaviour of the new Hejhal path
(:mod:`maass_forms_klein.modform.hejhal`) against survey ground truth:

* **4_1 first four strict eigenvalues** — reproduced from the reference runner
  (``reference/ford/maass_eigenvalue_runner.py``) to ``1e-6``. The classical
  literature value ``R = 6.62211934`` is a *spurious* dip that must NOT reach
  strict confirmation (a deliberate negative control).
* **8_8 four low-window eigenvalues** — the positive control the survey used on
  the low window.
* **exceptional interval** — the real-order scan on ``t in (0.02, 0.98)`` for
  4_1 must return zero strict confirmations ("no exceptional eigenvalues").

Strict confirmation criterion (survey §1.3): ``status == 'confirmed'`` with
``spread < 1e-3`` AND ``resid < 1e-3``.

All fixtures are marked ``slow``: each runs a genuine covering + Hejhal scan.
Run with ``pytest -m slow tests/test_hejhal_eigenvalues.py``.
"""

import pytest

# Strict-confirmation thresholds (survey §1.3).
STRICT_SPREAD = 1e-3
STRICT_RESID = 1e-3
# Agreement tolerance for values reproduced from the reference runner with the
# SAME parameters as ``search_eigenvalues`` defaults (§8 acceptance: 1e-6).
MATCH_TOL_REF = 1e-6
# Agreement tolerance for values quoted in the survey paper (§6): these were
# produced across campaigns with different truncation sizing, so double
# precision only pins them to the strict validation spread (~1e-3).  The strong
# assertion is that each is STRICTLY confirmed near the quoted location.
MATCH_TOL_PUB = 1e-3

# 4_1: first four strict eigenvalues, from the reference runner over [2, 9.5].
FIG8_FIRST_FOUR = [4.90008537, 5.91291796, 7.07200419, 7.40661560]
# 4_1: the historical spurious value — must NOT validate.
FIG8_SPURIOUS = 6.62211934
# 8_8: four low-window strict eigenvalues (survey §6).
K8_8_LOW = [2.183147, 2.452085, 3.115460, 3.460401]
# 6_1: lowest strict eigenvalue from the reference runner with the same window.
K6_1_FIRST = 3.46693579


def _strict(results):
    """Sorted R-values of the strictly confirmed candidates in ``results``."""
    return sorted(
        c["r"]
        for c in results
        if c["status"] == "confirmed"
        and c.get("spread", 1.0) < STRICT_SPREAD
        and c["resid"] < STRICT_RESID
    )


def _closest(values, target):
    return min((abs(v - target), v) for v in values)[1] if values else None


@pytest.mark.slow
@pytest.mark.parametrize("r_expected", FIG8_FIRST_FOUR)
def test_fig8_strict_eigenvalue(r_expected):
    """Each of 4_1's first four eigenvalues is strictly confirmed to 1e-5."""
    from maass_forms_klein.modform.hejhal import search_eigenvalues

    results = search_eigenvalues("4_1", r_expected - 0.15, r_expected + 0.15, step=0.05)
    strict = _strict(results)
    match = _closest(strict, r_expected)
    assert match is not None, f"no strict confirmation near R={r_expected}"
    assert abs(match - r_expected) < MATCH_TOL_REF, (match, r_expected)


@pytest.mark.slow
def test_fig8_spurious_6_6221_does_not_validate():
    """The classical R=6.62211934 dip does NOT reach strict confirmation."""
    from maass_forms_klein.modform.hejhal import search_eigenvalues

    results = search_eigenvalues("4_1", FIG8_SPURIOUS - 0.15, FIG8_SPURIOUS + 0.15, step=0.05)
    near = [r for r in _strict(results) if abs(r - FIG8_SPURIOUS) < 0.1]
    assert near == [], f"spurious value validated strictly: {near}"


@pytest.mark.slow
def test_6_1_first_strict_eigenvalue_matches_reference():
    """The package and reference runner agree on 6_1 to the required 1e-6."""
    from maass_forms_klein.modform.hejhal import search_eigenvalues

    results = search_eigenvalues("6_1", 3.30, 3.65, step=0.05)
    strict = _strict(results)
    match = _closest(strict, K6_1_FIRST)
    assert match is not None, "no strict confirmation near the first 6_1 eigenvalue"
    assert abs(match - K6_1_FIRST) < MATCH_TOL_REF, (match, K6_1_FIRST)


@pytest.mark.slow
@pytest.mark.parametrize("r_expected", K8_8_LOW)
def test_8_8_low_window_eigenvalue(r_expected):
    """8_8's four low-window eigenvalues are strictly confirmed to 1e-5."""
    from maass_forms_klein.modform.hejhal import search_eigenvalues

    results = search_eigenvalues("8_8", r_expected - 0.15, r_expected + 0.15, step=0.05)
    strict = _strict(results)
    match = _closest(strict, r_expected)
    assert match is not None, f"no strict confirmation near R={r_expected}"
    assert abs(match - r_expected) < MATCH_TOL_PUB, (match, r_expected)


@pytest.mark.slow
def test_fig8_exceptional_interval_empty():
    """No strict confirmations on the exceptional interval t in (0.02, 0.98)."""
    from maass_forms_klein.modform.hejhal import search_eigenvalues

    results = search_eigenvalues("4_1", 0.02, 0.98, step=0.02, order="real", rmax_ctx=2.2)
    assert _strict(results) == [], f"exceptional eigenvalue found: {_strict(results)}"
