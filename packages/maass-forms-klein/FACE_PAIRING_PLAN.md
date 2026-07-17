# Implementation plan: Ford-domain face pairings in `maass_forms_klein`

Status: proposal (2026-07-12). Companion to the paper section "Computing covering
generators for the Ford domain", subsection "From covering sets to face pairings".
Prototype code: `KnotMaass/scripts/face_pairing_prototype.py` (standalone
snappy + numpy/scipy; validated on 4_1, 5_2, 6_1, 7_4 — see "Validated facts" below).

## Handoff notes (for a local Claude Code session)

Working directory: `packages/maass-forms-klein/` in the maass-forms-mono
repository. Read the package `CLAUDE.md` first and follow it: venv_test +
`make install` after every file change, Sage doctests with `EXAMPLES::` on all
public functions (use `maass_form_core.testing.connect_mockdb` where the DB is
touched), ruff clean, line length 100.

Everything in this plan has a validated, standalone reference implementation
in `dev/face_pairing/` (plain python3, no Sage) — start by RUNNING those
scripts and reproducing their printed numbers, then port them into
`src/maass_forms_klein/hyperbolic_space/face_pairing.py` phase by phase.
`dev/face_pairing/README.md` lists six traps already discovered the hard way
(sphere translation side, centre conventions, pull-back sanity check, grid vs
exact Y0, mod-Lambda class merging, skew dual boxes) — read it before writing
geometry code. The acceptance tests are the regression table below; do not
mark a phase done until the corresponding rows reproduce inside the package
(as doctests or pytest in `tests/test_face_pairing.py`).

Suggested order: phase 1 (visibility + volume + EP checks, port of
`face_pairing_prototype.py`), phase 2 (exact Y0, port of
`exact_floor_height.py`, wire into `find_max_y`), phase 3 (horoball
enumeration — recipe in the dev README, ~1 day, closes 7_4), phase 4 last and
optional. Keep commits per phase.

## Goal

Extend the covering-generator machinery (`hyperbolic_space/utils.find_covering_generators`,
`KleinianGroup_class.covering_generators`) so the package can compute, verify and cache
the *face-pairing structure* of the Ford domain:

1. extract the visible (face-carrying) hemispheres from any candidate list
   (power-diagram / LP visibility);
2. enumerate candidates by |c| with the a priori bound |c| <= 1/Y0 instead of by
   word length (correct termination criterion for face completeness);
3. verify completeness via the volume identity vol = ∫_P dx / (2 h(x)^2) and via
   the Epstein–Penner duality count (face pairs = edges of the canonical
   decomposition, from snappy's `canonical_retriangulation`);
4. derive exact floor height Y0, inversion pairing, and (later) edge cycles /
   Poincaré presentation.

## Validated facts the implementation must reproduce (regression targets)

| knot | faces | pairs | EP edges | max word len | Y0 (face-complete) | vol residual |
|------|-------|-------|----------|--------------|--------------------|--------------|
| 4_1  | 4     | 2     | 2        | 2            | sqrt(2/3) = 0.816497 | < 1e-6 |
| 5_2  | 8     | 4     | 4        | 7            | 0.55969            | < 5e-7 |
| 6_1  | 12    | 6     | 6        | 6            | 0.40228            | < 4e-7 |
| 7_4  | 16    | 8     | 8        | >= 14 (last pair via horoballs) | 0.44377 | < 2e-6 |

(*) 7_4 is the acceptance test for the |c|-bounded enumeration (phase 3): both
completeness tests certify one face pair with |c| <= 2.35 and word length >= 14
that breadth-first search cannot reach in practice (2e6 distinct bottom rows at
length 13). Note: covering-stage Y0 estimates in the existing table
(e.g. 0.397 for 5_2) are lower bounds; face-complete values supersede them.

Additional validation (2026-07-13): the horoball pipeline closes ALL 20
hyperbolic 8-crossing knots (19 with pairs = EP edges and |vol residual|
<= 2e-5; 8_18 volume-complete with a non-simplicial canonical decomposition,
EP guard engaged). Requires traps 7-8 in dev/face_pairing/README.md
(Gauss-reduced lattice basis; LP counts for sliver faces - 8_1 has a face
pair of area share 2.5e-6).

## Module layout

New module: `src/maass_forms_klein/hyperbolic_space/face_pairing.py`
(pure Python; no new Cython initially). Reuses:

- `types.py`: `Circle`, `Parallelogram`, `Rectangle` (extend with `PowerCell`,
  `Face` dataclasses);
- `word_utils.py`: `word_to_element`, `word_to_circle`, `find_inverse_word`,
  `reduce_word`;
- `geometry_utils.py`: `matrix_to_circle`, dedup helpers;
- `parallelogram.py`: cell geometry;
- `kleinian_group.py`: `KleinianGroup_class` (new cached methods, below);
- `maass_form_core.testing.connect_mockdb` for doctests.

## API (phase-tagged)

```python
# --- phase 1: visibility + volume test (grid + LP) ------------------------
def visible_faces(
    group: KleinianGroup_class,
    candidates: list[tuple[str, Circle]] | None = None,   # default: BFS list
    method: str = "lp",            # "lp" (exact, scipy.optimize.linprog/HiGHS)
                                   # or "grid" (fast pre-pass)
    prec: int = 53,
) -> list[Face]:
    """Face = (word, circle, inverse_word, witness, margin).
    Classes merged modulo the translation lattice (centre equal mod Lambda,
    equal radius); pairing via find_inverse_word + class lookup."""

def floor_height(group, faces=None, exact: bool = False) -> Real_t:
    """Y0. exact=True: min over power-diagram vertices/edges of the envelope
    (phase 2); else adaptive grid lower bound. Replaces the trial-and-error
    loop in compute_coefficients.find_max_y."""

def ford_volume(group, faces=None, rel_tol: float = 1e-6) -> Real_t:
    """Adaptive quadrature of 1/(2 h^2) over the cusp cell; Richardson
    extrapolation on dyadic grids until rel_tol or non-convergence."""

def check_face_completeness(group) -> dict:
    """{'volume_residual': ..., 'ep_edges': ..., 'pairs_found': ...,
        'inversion_closed': bool, 'complete': bool}.
    ep_edges via snappy canonical_retriangulation (guard: has_finite_vertices
    => canonical decomposition non-simplicial; then skip the count, warn)."""

# --- phase 3: bounded-|c| enumeration via horoballs -------------------------
def enumerate_bounded_c(
    group, c_bound: Real_t, prec: int = 53,
) -> list[tuple[str | None, Circle]]:
    """All double cosets with |c| <= c_bound, via SnapPy's cusp-neighbourhood
    horoball enumeration (Weeks): elements with |c| <= X are exactly the
    horoballs of Euclidean diameter >= const/X^2, and
    manifold.cusp_neighborhood().horoballs(cutoff) enumerates them
    COMPLETELY from the triangulation. Steps: (1) maximal cusp,
    translations (m_t, l_t); (2) rescale centres by 1/m_t into the
    meridian-normalised frame; (3) circle radius r = A*sqrt(diam), constant
    A fixed by the cusp normalisation (calibrate once against a known |c|);
    (4) reduce mod Lambda, dedupe. Word representatives, when needed,
    recovered by pull-back of a point below the horoball centre.
    VALIDATED (2026-07-13, in session): for 7_4, 656 horoballs give the
    complete 20-sphere list with r >= 0.40, including the inverse pair of
    radius 0.45340 invisible to row-BFS at word length 13; the face set
    closes to 16 faces / 8 pairs (= EP edges), volume residual < 2e-6,
    corrected Y0 = 0.44377. The row-BFS is retained only as an independent
    cross-check; a self-contained (SnapPy-independent, interval-certified)
    enumeration is a long-term nice-to-have, no longer blocking."""

# --- phase 4: pairing/presentation ------------------------------------------
def edge_cycles(group, faces) -> list[list]:
    """Edges of the power diagram + cycle relations (Poincaré)."""
def presentation_from_faces(group) -> tuple[list, list]:
    """Generators (face words + L, M) and relations; compare with
    group.manifold().fundamental_group() as a doctest check."""
```

`KleinianGroup_class` additions (all `@cached_method`, DB-backed like
`covering_generators`):

```python
def ford_faces(self, certified: bool = False) -> list[Face]: ...
def ford_floor_height(self) -> Real_t: ...      # replaces ad hoc find_max_y
def ford_volume(self) -> Real_t: ...
def check_ford_complete(self) -> dict: ...
```

## Database (extends `database/models.py`)

New document `FordDomainDB` (pattern: `KleinianGroupDB` / `Word`):

- `label` (manifold), `faces` (list of {word, centre, radius, partner_word}),
- `y0`, `y0_exact` (bool), `volume_residual`, `ep_edges`, `complete` (bool),
- `c_bound_enumerated` (how far the |c| enumeration was pushed),
- provenance: package version, precision, date.

Seeding logic mirrors `covering_generators`: on construction, try
`FordDomainDB.objects(label=...)`; recertify cheaply (LP margins + volume)
rather than recompute.

## Integration points

1. `compute_coefficients.get_pb_pts_set_params` / `find_max_y`: replace the
   "decrease Y by 2% until pullback works" loop with `ford_floor_height()`;
   larger certified Y improves truncation (M0 ~ (R + D ln 10)/(2 pi delta Y))
   and conditioning. Keep the old path as fallback when faces are unknown.
2. `search/search.py`: no change needed (locators are Y-agnostic), but expose
   Y0 in `compute_one_spectral_parameter` metadata for the DB record.
3. Paper appendix table ("number-field data for the knots considered"):
   `ford_faces` + `ford_floor_height` fill the columns delta, max sphere
   radius, number of side pairings, Y0 directly.
4. `covering_generators`: after phase 3, seed the covering search with the
   face set (it *is* a covering set — Lemma in the paper), making the
   covering step a pure verification.

## Testing

- Doctests (Sage style, `connect_mockdb`): 4_1 end-to-end in every public
  function (faces=4, pairs=2, Y0^2 = 2/3 to 1e-10, volume to 1e-6).
- `tests/test_face_pairing.py` (pytest, marker `slow` for 5_2/6_1/7_4):
  regression table above; inversion-closedness; LP margin positivity;
  EP guard path (mock a non-simplicial canonical decomposition).
- Property test: random subsets of the 6_1 face set must FAIL the volume
  test (residual > 0) — the detector detects.
- Numerical hygiene: all LP margins and volume residuals asserted against
  absolute tolerances stored next to the regression data, not inline magic
  numbers.

## Phases and effort

| phase | content | acceptance | est. effort |
|-------|---------|------------|-------------|
| 1 | `visible_faces` (grid+LP), `ford_volume`, `check_face_completeness`, grid `floor_height`; doctests 4_1 | 4_1/5_2/6_1 rows of the regression table | 2–3 days |
| 2 | exact `floor_height` + exact Y0 via power-diagram vertices; `FordDomainDB`; `find_max_y` integration | Y0(4_1)^2 == 2/3 exactly (symbolic for arithmetic case), DB round-trip | 2 days |
| 3 | `enumerate_bounded_c` via horoballs; 7_4 completion | 7_4: 8 pairs, volume residual < 1e-5, new pair's word via pull-back | ~1 day (validated in session) |
| 4 | edge cycles, Poincaré presentation, comparison with SnapPy relators | presentation isomorphic for 4_1, 5_2 | 1 week, optional for the paper |

## Risks / notes

- Phase-3 risk resolved: the horoball route is complete by construction
  (Weeks' algorithm). Residual dependency: SnapPy's cusp machinery as an
  oracle - acceptable (it already supplies the holonomy), and the volume/EP
  tests keep it honest. 7_4 regression target updated: 16 faces / 8 pairs,
  Y0 = 0.44377, residual < 2e-6.
- The LP path adds a scipy dependency for HiGHS; scipy is already transitively
  present in Sage environments — confirm in tox before relying on it.
- Everything is double precision; the LP margins (>= 1e-2 in all observed
  cases) dwarf holonomy error, but store margins in the DB so a future
  interval-arithmetic pass (Arb) can re-certify.
- `canonical_retriangulation` can be expensive for larger knots; cache
  `ep_edges` in `FordDomainDB`.
