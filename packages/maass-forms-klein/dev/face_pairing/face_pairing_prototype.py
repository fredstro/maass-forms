"""
Prototype: extract Ford-domain face pairings from the covering enumeration.

Validated results (2026-07-12), reproduced by running this file:
    4_1 : 4 faces / 2 pairs,  EP edges 2, volume residual < 1e-6
    5_2 : 8 faces / 4 pairs,  EP edges 4, volume residual < 5e-7, Y0 = 0.55969
    6_1 : 12 faces / 6 pairs, EP edges 6, volume residual < 4e-7
    7_4 : 14 faces / 7 pairs found up to word length 13; EP edges 8 and
          volume residual +6.0e-4 certify ONE missing pair at length >= 14
          (radius >= Y0 = 0.427, i.e. |c| <= 2.35).

Companion to FACE_PAIRING_PLAN.md; standalone (snappy + numpy + scipy),
independent of the maass_forms_klein package. Depends on
compute_ford_cover.py (same directory) for the cusp-normalised holonomy.

Usage:  python3 face_pairing_prototype.py 4_1 [5_2 ...] [--maxlen N]
"""
import argparse
import numpy as np
from scipy.optimize import linprog

from compute_ford_cover import normalise_group, lattice_reduce


def enumerate_circles(gens, v1, v2, maxlen=10, rmin=0.05):
    """BFS on bottom rows (row determines the left coset for torsion-free
    groups); returns dict {(centre mod Lambda, radius) -> word}."""
    letters = ['a', 'b', 'A', 'B']
    gm = {g: (complex(gens[g][0, 0]), complex(gens[g][0, 1]),
              complex(gens[g][1, 0]), complex(gens[g][1, 1])) for g in letters}
    seen, frontier = set(), []
    for g in letters:
        row = (gm[g][2], gm[g][3])
        key = tuple(round(x, 5) for x in (row[0].real, row[0].imag,
                                          row[1].real, row[1].imag))
        seen.add(key)
        frontier.append((g, row))
    good = {}
    for length in range(1, maxlen + 1):
        for w, (c, d) in frontier:
            if abs(c) > 1e-9 and 1.0/abs(c) >= rmin:
                cen, r = -d/c, 1.0/abs(c)
                cred, _ = lattice_reduce(cen, v1, v2)
                k = (round(cred.real, 6), round(cred.imag, 6), round(r, 6))
                good.setdefault(k, (w, cred, r))
        if length == maxlen:
            break
        nf = []
        for w, (c, d) in frontier:
            for g in letters:
                if w[-1] != g.swapcase():
                    a11, a12, a21, a22 = gm[g]
                    row = (c*a11 + d*a21, c*a12 + d*a22)
                    key = tuple(round(x, 5) for x in (row[0].real, row[0].imag,
                                                      row[1].real, row[1].imag))
                    if key not in seen:
                        seen.add(key)
                        nf.append((w + g, row))
        frontier = nf
    return good


def lattice_dist(a, b, v1, v2):
    B = np.array([[v1.real, v2.real], [v1.imag, v2.imag]])
    x = np.linalg.solve(B, [(a - b).real, (a - b).imag])
    x -= np.round(x)
    return abs(x[0]*v1 + x[1]*v2)


def envelope(items, v1, v2, base, n):
    C = np.array([c for _, c, _ in items])
    R = np.array([r for _, _, r in items])
    s = np.linspace(0, 1, n, endpoint=False) + 0.5/n
    u, w = np.meshgrid(s, s)
    Z = base + u*v1 + w*v2
    best = np.full(Z.shape, -1)
    besth = np.full(Z.shape, -1.0)
    for i in range(len(C)):
        h2 = R[i]**2 - np.abs(Z - C[i])**2
        m = h2 > besth
        best[m] = i
        besth[m] = h2[m]
    return best, besth


def lp_visible(items, v1, v2, base):
    """Exact visibility over the candidate list: sphere i carries a face iff
    the LP  max t  s.t.  l_i(x) - l_j(x) >= t (j != i),  x in P  has a
    positive optimum, where l_i(x) = 2<c_i,x> + r_i^2 - |c_i|^2."""
    ded = {}
    for w, c, r in items:
        ded.setdefault((round(c.real, 8), round(c.imag, 8), round(r, 8)),
                       (w, c, r))
    items = list(ded.values())
    Binv = np.linalg.inv(np.array([[v1.real, v2.real], [v1.imag, v2.imag]]))
    bb = np.array([base.real, base.imag])
    A_cell, b_cell = [], []
    for row in (0, 1):
        A_cell.append(-Binv[row]); b_cell.append(Binv[row] @ (-bb))
        A_cell.append(Binv[row]);  b_cell.append(1 + Binv[row] @ bb)
    L = np.array([[2*c.real, 2*c.imag, r*r - abs(c)**2] for _, c, r in items])
    out = []
    for i in range(len(items)):
        A_ub = [[-(L[i]-L[j])[0], -(L[i]-L[j])[1], 1.0]
                for j in range(len(items)) if j != i]
        b_ub = [(L[i]-L[j])[2] for j in range(len(items)) if j != i]
        A_ub += [[a[0], a[1], 0.0] for a in A_cell]
        b_ub += list(b_cell)
        res = linprog(c=[0, 0, -1], A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                      bounds=[(None, None)]*3, method='highs')
        if res.status == 0 and -res.fun > 1e-12:
            out.append((items[i], -res.fun))
    return out


def inv_word(w):
    return ''.join(ch.swapcase() for ch in reversed(w))


def analyse(name, maxlen=10):
    import snappy
    M, G, gens, v1, v2 = normalise_group(name)
    base = -(v1 + v2)/2
    circles = enumerate_circles(gens, v1, v2, maxlen=maxlen, rmin=0.05)
    # bootstrap Y0 from the largest circles, then filter: a sphere with
    # r < Y0 can never touch the Ford floor (h_i <= r < Y0 <= envelope),
    # so it is irrelevant to faces, floor height and volume alike.
    ranked = sorted(circles.values(), key=lambda t: -t[2])
    boot = [(w, c + i*v1 + j*v2, r) for w, c, r in ranked[:200]
            for i in (-1, 0, 1) for j in (-1, 0, 1)]
    _, bh0 = envelope(boot, v1, v2, base, 300)
    y0_lower = float(np.sqrt(bh0.min()))          # lower bound for true Y0
    allitems = [(w, c + i*v1 + j*v2, r) for w, c, r in ranked
                if r >= 0.9*y0_lower
                for i in (-1, 0, 1) for j in (-1, 0, 1)]
    _, besth = envelope(allitems, v1, v2, base, 500)
    Y0 = float(np.sqrt(besth.min()))
    cand = [(w, c, r) for w, c, r in allitems if r >= 0.95*Y0]
    faces = lp_visible(cand, v1, v2, base)
    # merge into double-coset classes
    classes = []
    for (w, c, r), margin in faces:
        for cl in classes:
            if abs(cl[2] - r) < 1e-8 and lattice_dist(c, cl[1], v1, v2) < 1e-6:
                break
        else:
            classes.append((w, c, r, margin))
    # inversion pairing (necessary condition for a face set)
    words = [cl[0] for cl in classes]
    from compute_ford_cover import circle_of
    def word_mat(wd):
        m = np.eye(2, dtype=complex)
        for ch in wd:
            m = m @ gens[ch]
        return m
    paired = 0
    for w, c, r, _ in classes:
        ci, ri = circle_of(word_mat(inv_word(w)))
        if any(abs(r2 - ri) < 1e-8 and lattice_dist(ci, c2, v1, v2) < 1e-6
               for _, c2, r2, _ in classes):
            paired += 1
    # volume test (converged grid)
    vol = None
    for n in (700, 1400):
        _, bh = envelope(allitems, v1, v2, base, n)
        area = abs((v1*np.conj(v2)).imag)/n**2
        vol = float(np.sum(area/(2*bh)))
    vol_true = float(snappy.Manifold(name).volume())
    # EP duality test
    C = snappy.Manifold(name).canonical_retriangulation()
    simplicial = not C.has_finite_vertices()
    ep_edges = C.num_tetrahedra() if simplicial else None
    print(f"=== {name} (enumeration to word length {maxlen}) ===")
    print(f"  faces: {len(classes)}  (pairs {len(classes)/2}), "
          f"inversion-closed: {paired == len(classes)}")
    print(f"  Y0 = {Y0:.5f}")
    print(f"  volume: {vol:.7f}  vs true {vol_true:.7f}  "
          f"residual {vol - vol_true:+.2e}")
    print(f"  EP edges: {ep_edges if simplicial else 'non-simplicial'}"
          f"  -> {'COMPLETE' if simplicial and ep_edges == len(classes)/2 and abs(vol-vol_true) < 1e-5 else 'INCOMPLETE or unverified'}")
    for w, c, r, m in sorted(classes, key=lambda t: (-t[2], t[0])):
        print(f"    {w:<14s} r={r:.6f}  margin={m:.2e}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('knots', nargs='*', default=['4_1'])
    ap.add_argument('--maxlen', type=int, default=10)
    args = ap.parse_args()
    for k in args.knots:
        analyse(k, maxlen=args.maxlen)
