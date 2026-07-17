# dev/face_pairing — validated prototypes for FACE_PAIRING_PLAN.md

Standalone reference implementations (plain python3: snappy + numpy + scipy,
NO Sage required) of every algorithm in `../../FACE_PAIRING_PLAN.md`, exactly
as validated in the July 2026 research sessions. Port these into
`src/maass_forms_klein/hyperbolic_space/face_pairing.py` following the
package conventions (Sage doctests, ruff, `make install` after every change —
see the package `CLAUDE.md`). The numbers these scripts print are the
regression targets in the plan; the port must reproduce them.

## Files

- `compute_ford_cover.py` — cusp-normalised holonomy (meridian scaled to 1),
  covering-generator BFS with certified quadtree cover test, greedy
  reduction, grid floor height. Run: `python3 compute_ford_cover.py 4_1 5_2`.
- `face_pairing_prototype.py` — bounded-|c| circle enumeration by row-BFS,
  LP visibility (faces), class merging mod Lambda, inversion pairing check,
  volume test, EP duality test. Run:
  `python3 face_pairing_prototype.py 4_1 5_2 --maxlen 9`
  (expect: 4_1 COMPLETE 2 pairs, 5_2 COMPLETE 4 pairs).
- `exact_floor_height.py` — exact Y0^theta via power-diagram vertices
  (radical centres). Run: `python3 exact_floor_height.py`
  (expect: 4_1 -> 0.816496580928 = sqrt(2/3); 7_4 -> 0.426277 on the
  bundled BFS-13 sphere list; the horoball-complete list gives 0.44377).
- `fig_data.json` — reference sphere lists (r >= 0.30) for 4_1 and 7_4 in
  the meridian-normalised frame, used by `exact_floor_height.py`'s demo.

## The horoball enumeration (phase 3) — validated recipe, not yet a script

```python
M = snappy.Manifold('7_4')
cn = M.cusp_neighborhood()
cn.set_displacement(cn.stopping_displacement())      # maximal cusp
mt, lt = (complex(x) for x in cn.translations())
H = cn.horoballs(cutoff)                             # COMPLETE above cutoff
# centre -> centre/mt  (meridian-normalised frame)
# sphere radius r = A*sqrt(2*h['radius']);  A calibrated once against a
# known sphere (derive exactly from the cusp normalisation in the port)
# reduce centres mod Lambda, dedupe -> full bounded-|c| circle list
```

Validated: 656 horoballs -> the complete 20-sphere list (r >= 0.40) for
7_4, including the radius-0.45340 inverse pair that row-BFS misses through
word length 13 (2e6 rows); faces close to 16/8 = EP edges, volume residual
< 2e-6, Y0 = 0.44377.

## 8-crossing validation (2026-07-13)

The horoball pipeline was run on all 20 hyperbolic 8-crossing knots
(8_19 is a torus knot). Result: ALL close — 19 with face pairs equal to the
EP edge count and volume residual <= 2e-5; 8_18 is volume-complete
(residual 7e-7) but its canonical decomposition is non-simplicial, so the
EP count guard correctly declines the comparison. Two additional traps were
found (7 and 8 below). Calibration is exact: A = 1/|meridian translation|
(at the maximal cusp the largest horoballs have diameter exactly 1).
The horoball frame differs from the holonomy frame by a global translation
(harmless for torus computations; account for it when recovering words).

## Traps discovered the hard way (do not rediscover)

1. Isometric spheres translate under RIGHT multiplication:
   circle(gamma T) is circle(gamma) shifted by -t for T = [[1,t],[0,1]];
   LEFT multiplication leaves the sphere fixed (bottom row unchanged).
2. `compute_ford_cover.covering_search` cover entries store centres
   ALREADY translated; the (n1, n2) bookkeeping must not be applied again.
3. After any pull-back construction, assert that (almost) no sample point
   at height Y < Y0 is outside every hemisphere — this catches both traps.
4. Grid minima OVERestimate Y0 (min over sample points); the power-vertex
   algorithm is exact and cheaper.
5. Dedup/class-merging must compare centres modulo Lambda with a tolerance
   (boundary rounding splits classes; cost us a wrong face count once).
6. The dual-lattice index box must be adaptive in k1 for skew cusp shapes
   (short vectors at large |k1|) — relevant if you touch the eigenvalue side.
7. GAUSS-REDUCE the cusp lattice basis (Lagrange-Gauss on (tau, 1)) before
   building cells and translate rings: for skew cusps (8_2 has
   tau = 10.98 - 2.99i) the geometric neighbour translates have coefficient
   vectors like (1, -11), far outside a +-1 ring, silently truncating the
   envelope (volume too big, spurious extra "faces").
8. Grid-based face COUNTS miss sliver faces (8_1 has a face pair of area
   share 2.5e-6 - a dozen pixels at n = 2200). Use the LP visibility test
   for counts; grids only for the volume integral and Y0 pre-passes.
