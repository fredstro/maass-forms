# 2026 eigenvalue survey: integration notes

The source survey spec is `knot-maass/LESSONS_LEARNED_EIGENVALUE_SURVEY.md`
in the sibling checkout. This package is `maass_forms_klein` in the monorepo.

## Indicator discrepancy (survey §1.1–1.2)

The prose spec calls for the smallest singular value of one row-equilibrated,
column-normalised system stacked at two heights. The named, battle-tested
`reference/ford/maass_eigenvalue_runner.py` instead solves the two single-height
systems separately and uses the aligned distance between their right singular
vectors. The current package path follows that executable reference.

A literal stacked-σ-min experiment on 2026-09-21 found the true 4_1 and 6_1
eigenvalues to better than `1e-6`, but it also strictly confirmed a false dip
at `R = 6.632802` inside the mandatory `R = 6.6221` negative-control window.
At that false dip the stacked σ-min was `1.9e-7`, while the two-height
alignment residual was `0.300`; at the true 4_1 value `R = 4.90008537`, the
alignment residual was `7.8e-4`. The stacked experiment also moved 8_8's
first low-window value by about `0.0038`. Do not switch the scan indicator to
bare stacked σ-min without an additional criterion that passes **all** survey
fixtures. The experimental code was reverted; see the regression tests.

## Exact floor and frame conventions (survey §3 and §4.2)

The package's `face_pairing.floor_height(..., exact=True)` already uses power
vertices. Its holonomy frame does **not** always have a unit meridian:

| Knot | Package meridian length | Package-frame floor | Unit-meridian floor |
| --- | ---: | ---: | ---: |
| 4_1 | 1.000000 | 0.8164966 | 0.8164966 |
| 5_2 | 0.138141 | 0.0772938 | 0.5595282 |
| 6_1 | 0.308210 | 0.1238928 | 0.4019750 |
| 7_4 | 2.559502 | 1.1356876 | 0.4437142 |

The Hejhal runner works in a unit-meridian frame, so the scalar height is
divided by the magnitude of the package meridian translation. A horizontal
origin offset does not change the minimum height over a complete lattice
cell. It **does** change individual sphere centres, so package horoball
centres must not be mixed with the runner's holonomy matrices without
recovering and validating the frame offset `δ`.

For 8_8, horoball-face auto-tuning cannot establish a positive floor in its
current group frame; the centred-cell NumPy cover remains quadtree-certified.
The fallback computes the exact power-vertex floor of that **selected cover**
(`0.3950617391`), not of the complete Ford sphere family, and records
`floor_source = "certified_cover"`. All sampling heights are guarded by
`min(y*) >= 0.98 Y0`; a failure is an error, not a reason to lower `Y0`.
