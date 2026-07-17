# Phase 4 — Edge-cycle Poincaré relations: hand-over for future work

Status: **incomplete / blocked**, precisely diagnosed. Written 2026-07-13.
Companion to `../../FACE_PAIRING_PLAN.md` (phase 4) and `README.md` in this
directory. Read those first.

This document is a self-contained hand-over so a future session can resume the
edge-cycle **relator** computation from the exact point of failure instead of
re-deriving everything. It records: what is already shipped, the two hard
constraints, the full edge-cycle algorithm that was attempted, every verified
sub-result, the precise failure and its root cause, and the recommended correct
construction.

---

## 1. Goal

Complete `presentation_from_faces` so it returns a genuine **Poincaré
presentation** derived from the Ford face pairing:

- **Generators**: the face-pairing transformations (one word per inverse face
  pair). *Already done and verified* (see §3).
- **Relators**: the **edge-cycle relations** — for each edge of the Ford domain,
  the product of face-pairings encountered while rotating once around the edge
  equals the identity (for a manifold; cone angle 2π). **This is the missing
  piece.**

Acceptance target from the plan: "presentation isomorphic for 4_1, 5_2".

---

## 2. Two hard constraints (context, not bugs)

1. **Word recovery is required and only available for BFS-reachable knots.** The
   edge-cycle product is a word in the face-pairing elements, so every face must
   carry a word. This holds for `4_1, 5_2, 6_1` (breadth-first reachable) but
   **not** for the horoball-completed knots (`7_4`, all 8-crossing), whose faces
   are word-less (see `README.md` / phase 3). So a presentation is only
   attainable for the word-bearing knots; `presentation_from_faces` already
   raises `ValueError` on word-less face sets.

2. **Group isomorphism is undecidable in general.** "Presentation isomorphic to
   SnapPy's" cannot be checked directly. The strongest *rigorous* verifications
   available are:
   - each relator word evaluates to ±I as a matrix (exact, since entries are
     algebraic — checked numerically at double precision, LP margins dwarf the
     error);
   - the abelianization of `<generators | relators>` equals `H_1` of the knot
     complement = `Z` (Smith normal form of the exponent-sum matrix);
   - relator count matches the deficiency (a knot group has deficiency 1: `k`
     generators need `k-1` independent relators).

---

## 3. What is already shipped (do not redo)

In `src/maass_forms_klein/hyperbolic_space/face_pairing.py`:

- **`edge_cycles(group, faces=None)`** — the power-diagram vertex/edge skeleton
  of the Ford floor. Returns each floor vertex as
  `{'point': vector, 'height': float, 'faces': (word_i, word_j, word_k)}` — the
  radical centres where three visible hemispheres are genuinely co-highest
  (apex points excluded via `|v - centre| > eps`; deduped modulo the lattice).
  For `4_1` it returns 12 vertices, all at height² = 2/3 = Y0². Helper:
  `_floor_vertices`, which reuses `_radical_centre`.

- **`presentation_from_faces(group)`** — returns `(generators, relators)` where
  `generators` = one word per inverse face pair (the Poincaré generators;
  correct by Poincaré + the certified Ford completeness), and `relators` =
  SnapPy's canonical relators (for the cross-check). For `4_1`: 2 generators
  (matches SnapPy's minimal rank); `5_2`, `6_1`: 4, 6 (non-minimal, expected).
  Verified: the holonomy satisfies the relators (±I).

Tests in `tests/test_face_pairing.py`: `test_edge_cycles_figure_eight`,
`test_presentation_from_faces` (4_1, 5_2). All green.

**The remaining work is purely: replace the placeholder `relators` (currently
SnapPy's) with the true edge-cycle relators computed from the geometry.**

---

## 4. Verified machinery (all in scratchpad prototypes; reusable)

### 4.1 Isometric circle of a group element

For a `numpy` `2×2` complex matrix `M = [[a,b],[c,d]]` (det 1), the isometric
sphere projects to the circle:

```
centre = -d/c,   radius = 1/|c|      (None if |c| ~ 0: parabolic/identity)
```

This is `_circle_of` in the module.

### 4.2 Isometry action on upper half-space H³ (VERIFIED)

Point of H³ as `(z, h)` with `z ∈ ℂ`, `h > 0` (the quaternion `z + h·j`). For
`M = [[a,b],[c,d]]`:

```python
def act(M, (z, h)):
    a, b, c, d = M[0,0], M[0,1], M[1,0], M[1,1]
    den = abs(c*z + d)**2 + abs(c)**2 * h**2
    z2  = ((a*z + b) * conj(c*z + d) + a * conj(c) * h**2) / den
    h2  = h / den
    return (z2, h2)
```

**Verified**: `g` maps the apex of `S_g` (at `(centre_g, radius_g)`) exactly to
the apex of `S_{g⁻¹}`. This is the correctness anchor — trust `act`.

### 4.3 Cusp translation lattice and the L/M representation

Common misread to avoid: for `4_1`, `L = [[-1, -3.4641],[0,-1]]` — the (0,1)
entry is `-3.4641` **real**. So (via `act` on a high point, empirically):

```
tauL = 3.4641   (real)     tauM = i   (imaginary)   -> rank 2, independent.
```

`L` translates `z -> z + tauL`, `M` translates `z -> z + tauM` (both act as pure
translations in PSL; `L`'s `-1` diagonal is `-I` up to sign, irrelevant in
PSL(2,ℂ)). The package `translation_lattice` basis is `[[3.4641, 0],[0, 1]]`
(as `[real, imag]` rows); `_lattice_vectors` returns these Gauss-reduced.

**Important:** `visible_faces` centres are reduced with the *Gauss-reduced*
basis, while the L/M translates use the *raw* `tauL, tauM` basis. These span the
same lattice Λ but different fundamental cells — do **not** mix a raw-basis
reduction with Gauss-basis face centres (that was a real bug: a mod-Λ
"reconstruction" gave 16/30 wrong on-top faces). Safe approach used: a single
large fixed arrangement, no reduction (§4.4).

### 4.4 The sphere arrangement (translates with exact element words)

Build the arrangement of all visible faces and their L/M translates, each tagged
with an **exact element word** (so the group element is exact — no reconstructing
from circles):

```python
def tword(n1, n2):
    return ('L'*n1 if n1>=0 else 'l'*(-n1)) + ('M'*n2 if n2>=0 else 'm'*(-n2))

arr = []                       # (word, matrix, centre, radius)
for f in visible_faces(G):
    for n1 in range(-8, 9):
        for n2 in range(-8, 9):
            w = f.word + tword(n1, n2)     # element = faceword · L^n1 · M^n2
            M = word_to_element(w) as numpy
            centre, radius = circle_of(M)
            if abs(centre) < 16:
                arr.append((w, M, centre, radius))
```

Key fact (why the translate word is `faceword · L^n1 M^n2`): the sphere of `g`
shifted by `+λ` is the sphere of `g · P_{-λ}` where `P_{-λ}` is translation by
`-λ` (RIGHT multiplication shifts the centre; LEFT multiplication by a parabolic
leaves the circle fixed — trap 1 in `README.md`). With `L, M` translating by
`+tauL, +tauM`, a translate by `λ = n1·tauL + n2·tauM` is `faceword · L^{-?}`…
in practice `faceword + tword(n1,n2)` sweeps the whole neighbourhood, so exact
signs are handled by enumerating all `(n1,n2)`.

`top2(z)` = the two arrangement spheres maximising `radius² - |z - centre|²`.
A point is on an **edge** iff the top two heights are equal (and a strict third
is lower); on a **vertex** iff the top three are equal.

---

## 5. The edge-cycle algorithm that was attempted

Theory (standard Poincaré, cusp at ∞):

- Ford domain `D` = exterior of all isometric spheres ∩ cusp fundamental cell.
- Face on sphere `I(g)` is paired to `I(g⁻¹)` by `g`.
- **Neighbour relation (derived, believed correct):** the copy of `D` adjacent
  across face `I(g)` is `g⁻¹·D`. (Because `g` maps face `I(g)` of `D` onto face
  `I(g⁻¹)` of `D`, so `D` and `g·D` are adjacent across `I(g⁻¹)`; substituting
  gives: across `I(f)` the neighbour is `f⁻¹·D`.)
- **Edge cycle:** fix a point `p ∈ H³` on an edge `E`. Copies around `E` are
  `D = W_0 D, W_1 D, …`. In copy `i`, `p` sits at `q_i = W_i⁻¹(p)`, which lies on
  `E` in `D`-coordinates ⇒ two equal-top spheres. Crossing the outgoing sphere
  (element `g_out`, in `D`-coords) advances `W_{i+1} = W_i · g_out⁻¹`. The
  incoming face in the next copy is `I(g_out⁻¹)`. Terminate when `W_m = ±I`; the
  relator is the accumulated word (should evaluate to ±I).

Pseudocode of the walk (the version tried):

```python
W = I; incoming = None
loop:
    q = act(inv(W), p)                       # p in D-coords of current copy
    i1, i2, h1, h2 = top2(q.z)               # two edge faces (must be equal height)
    if incoming is None: cross = i1          # arbitrary first choice
    else:                                    # identify incoming by circle, cross the other
        cross = the index whose circle is FARTHER from `incoming`
    S = arr[cross].matrix
    W = W @ inv(S)                           # neighbour across I(g_out)
    incoming = circle_of(inv(S))             # = I(g_out^{-1})
    if step >= 1 and W ≈ ±I: return relator
```

Start point: from a floor vertex `v` with faces `(a,b,c)`, step a small `ε`
along the radical axis of a pair `(a,b)` (direction `1j·(c_b - c_a)/|c_b - c_a|`)
into the edge interior; keep the direction where the two spheres stay co-top.

---

## 6. The failure (precisely characterised)

**Symptom:** the walk provably stays on edges (top-two heights equal every step),
but **no edge ever closes to ±I**. Instead the accumulated `W` **drifts by a
parabolic** (lattice translation).

Concrete `4_1` trace (start edge, faces `AAABBBlmmm`, `abM`):

```
s0: q=(0.433+0.25j)  cross=AAABBBlmmm   |W±I|=0
s1: q=(-1.732+0.5j)  cross=BAbABBa      |W±I|=1
s2: q=(3.031+1.75j)  cross=BAlm         |W±I|=4
s3: q=(0.433+1.25j)  cross=AAABBBlmmmm  |W±I|=1     <- s3.z = s0.z + i = s0 + tauM
s4: q=(-1.732+0.5j)  ... repeats s1 ...
```

After 3 steps the walk returns to the **start edge shifted by `tauM`**, then
loops forever, so `W` grows and never reaches ±I.

**Convention sweep (rules out simple sign errors):** all 8 combinations of
{cross the incoming-match vs the other} × {left vs right `W` accumulation} ×
{`act(W)` vs `act(inv(W))`}, over all 24 edge points of `4_1` → **0 closed**.
So the failure is *not* a direction/sign convention.

**Root cause (confirmed): coset-representative ambiguity.**
An isometric **circle** determines the group element only **up to left
multiplication by the parabolic stabiliser `Γ_∞`** (`g` and `P·g` share a circle
but act differently). The walk identifies the outgoing face by its circle, then
crosses using the arrangement's representative (`faceword·tword`). That
representative is a *valid* face-pairing (both its sphere and inverse-sphere are
visible faces), **but not necessarily the one whose neighbour copy `g_out⁻¹·D`
continues the rotation around the fixed edge.** Choosing the wrong coset rep
injects a fixed parabolic each cycle — exactly the observed `tauM` drift. No
choice among circle-shadow representatives fixes this, because the information
needed (which lift of the edge geodesic) is not in the 2-D shadow.

---

## 7. Recommended correct construction (the fix)

Stop resolving the outgoing element from the **shadow circle**. Instead track the
edge as a **geodesic in H³** and require the pairing to fix it:

1. **Represent an edge by its geodesic**, i.e. the two endpoints of the H³
   geodesic containing it. A finite edge runs between two floor vertices
   `v1 = (z1, h1)`, `v2 = (z2, h2)` (both computable — they are the radical
   centres in `_floor_vertices`). The geodesic's ideal endpoints on `∂H³ = ℂ∪∞`
   are computable from `v1, v2` (the semicircle/vertical line through them).

2. **Rotate around the edge by its exact stabiliser step.** The isometries that
   move `D` to the next wedge around `E` must **map the edge geodesic to itself**
   (they fix the axis of the edge). Among all elements whose isometric sphere is
   the outgoing face (the `Γ_∞`-coset), select the unique `g_out` such that
   `g_out` maps the edge geodesic onto the corresponding edge geodesic of the
   neighbour — equivalently, whose ideal endpoints transform consistently.

3. **Terminate** when the accumulated `W` maps the edge geodesic back to itself
   with the same orientation; for a manifold this forces `W = ±I` and yields the
   relator. Verify: relator matrix ±I, then abelianisation `= Z`.

Equivalently, this is the "algorithmic Poincaré theorem" with **exact edge
tracking** (see e.g. Epstein–Petronio, or Riley's / Marden's treatments; SnapPy's
own `fundamental_group` builds relators from the triangulation this way). The
2-D power-diagram walk is not enough on its own precisely because of the coset
ambiguity; the edge geodesic supplies the missing lift.

A cheaper practical alternative worth trying first: at each step, instead of
picking the coset rep `faceword·tword`, **choose the outgoing element by
look-ahead** — among `{ P · (candidate) : P ∈ small Λ }`, pick the one whose
`g_out⁻¹·D` places `q_{next}` back on the *same* edge geodesic (not a translate).
If that single change makes `4_1` close to ±I over the 24 edges, it confirms the
diagnosis and likely generalises.

---

## 8. Verification checklist for a completed implementation

- [ ] `4_1`: an edge cycle closes; relator matrix ≈ ±I (atol 1e-5).
- [ ] Collect relators over all edge orbits; reduce to independent set.
- [ ] Abelianisation (Smith normal form of exponent-sum matrix) `= Z`.
- [ ] `5_2`: same; 4 generators, 3 independent relators, `H_1 = Z`.
- [ ] Compare with `group.manifold().fundamental_group()`: same abelianisation;
      ideally each SnapPy relator is a consequence (spot-check numerically).
- [ ] Wire into `presentation_from_faces` (replace the SnapPy-relator
      placeholder); add pytest for `4_1`/`5_2` relators → ±I and `H_1 = Z`.

---

## 9. Environment / reproduction notes

- Working env: `/Users/fredrik/.virtualenvs/maass-forms` (passagemath 10.5.22,
  snappy 3.3.2, python 3.13); `sage` on its `bin/`. Run mono src via
  `PYTHONPATH=<core/src>:<hilbert/src>:<klein/src>`; needs `mongomock`,
  `freezegun` pip-installed. (See the `klein-working-sage-env` memory.)
- For pure numeric prototyping (no Sage): snappy+numpy+scipy in an isolated venv
  works, but the edge-cycle code needs `word_to_element` from the package, so use
  the sage env above.
- Reference numbers used above are for `KleinianGroup('4_1')` at `prec=53`:
  `tauL = 3.4641016…` (real), `tauM = i`; floor height `Y0 = sqrt(2/3) =
  0.816497`; 4 faces / 2 pairs; face generators `['AAABBB', 'ab']`.
- All the phase-4 edge-cycle experiments lived in scratchpad files `ec*.py`,
  `step1.py`, `verify_recon.py`, `search.py` — none were committed; the module is
  unchanged and green.
