"""Exact Ford-floor integration and pull-back regression checks (survey section 3)."""

import numpy as np
import pytest

from maass_form_core.testing import connect_mockdb

connect_mockdb()

from maass_forms_klein.modform.hejhal import HejhalContext  # noqa: E402


# Heights in the meridian-normalised frame of the pure-NumPy survey runner.
FLOOR_HEIGHTS = {
    "4_1": 0.816496580927726,
    "5_2": 0.5595282,
    "6_1": 0.4019750,
    "7_4": 0.4437142,
}


@pytest.mark.slow
@pytest.mark.parametrize("knot, expected", FLOOR_HEIGHTS.items())
def test_certified_floor_matches_survey_frame(knot, expected):
    """Wrong meridian scale or an incomplete face set changes this height."""
    from maass_forms_klein.modform.survey_floor import certified_floor_height

    assert abs(certified_floor_height(knot) - expected) < 1e-6


@pytest.mark.slow
@pytest.mark.parametrize("knot", ("4_1", "6_1"))
def test_hejhal_uses_exact_floor_and_pullback_raises_all_samples(knot):
    """Grid-biased Y0 or a misplaced cover violates the survey smoke test."""
    context = HejhalContext(knot, rmax=3.7)
    assert abs(context.Y0 - FLOOR_HEIGHTS[knot]) < 1e-6
    for height in (context.Y1, context.Y2):
        assert min(context.pre[height]["ystar"]) >= 0.98 * context.Y0


def test_pullback_smoke_guard_rejects_a_low_or_nonfinite_sample():
    """A partly misplaced cover must not silently pass validation."""
    from maass_forms_klein.modform.hejhal import _require_pullback_above_floor

    _require_pullback_above_floor(np.array([0.98, 1.1]), 1.0)
    with pytest.raises(RuntimeError, match=r"below 0\.98"):
        _require_pullback_above_floor(np.array([1.1, 0.979]), 1.0)
    with pytest.raises(RuntimeError, match="non-finite"):
        _require_pullback_above_floor(np.array([1.1, np.nan]), 1.0)


@pytest.mark.slow
def test_8_8_uses_exact_certified_cover_floor_when_horoball_frame_fails():
    """The failed horoball auto-bound must not block a covered Hejhal run."""
    context = HejhalContext("8_8", rmax=2.3)
    assert context.floor_source == "certified_cover"
    assert abs(context.Y0 - 0.3950617391458645) < 1e-6
    assert context.Y0 < context.cover_floor_estimate
    for height in (context.Y1, context.Y2):
        assert min(context.pre[height]["ystar"]) >= 0.98 * context.Y0


def test_pullback_reduces_cell_boundary_to_centred_half_open_interval():
    """A half-lattice point must agree with the group's centred-cell convention."""
    from maass_forms_klein.modform.hejhal import pullback_points

    points, heights = pullback_points(np.array([0.5 + 0j, -0.5 + 0j]), 0.4, [], 2j)
    assert np.allclose(points, [-0.5 + 0j, -0.5 + 0j])
    assert np.allclose(heights, [0.4, 0.4])
