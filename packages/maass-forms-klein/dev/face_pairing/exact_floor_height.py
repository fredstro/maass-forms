"""
Exact floor height Y_0^theta of a finite hemisphere family, via power-diagram
vertices (radical centres).  Companion to the proposition "Exact computation
of the floor height" in the Ford-domain section.

The floor  h_theta(x)^2 = max_i (theta^2 r_i^2 - |x - c_i|^2)  is a maximum
of concave paraboloids whose argmax regions are the cells of the power
diagram with weights w_i = theta^2 r_i^2; a concave function on a convex
cell is minimised at a vertex, so

    (Y_0^theta)^2 = min over radical centres v_ijk (reduced mod Lambda)
                    of  max_i (w_i - |v - c_i|^2).

Evaluating over ALL triples is safe (superset of the true vertices, each
value >= the minimum). Certification: the horosphere at height Y is below
the theta-floor iff the discs of radii sqrt(theta^2 r_i^2 - Y^2) cover P,
so the quadtree covering test brackets the result rigorously.

Validation:
    4_1 : Y_0 = 0.816496580928 = sqrt(2/3) (exact), x* = -i/sqrt(3)
    7_4 : Y_0 = 0.426277 at x* = 0.15897 - 1.75177i
"""
import numpy as np


def exact_floor(circles, v1, v2, theta=1.0):
    """circles: [(centre complex, radius float)] including all
    Lambda-translates meeting a neighbourhood of the base cell.
    Returns (Y0_theta, deepest_point, defining_triple)."""
    C = np.array([c for c, _ in circles])
    W = np.array([(theta * r) ** 2 for _, r in circles])
    R = np.sqrt(W)
    kap = np.abs(C) ** 2 - W
    B = np.array([[v1.real, v2.real], [v1.imag, v2.imag]])
    Binv = np.linalg.inv(B)
    n = len(C)
    best, argbest, tri = np.inf, None, None
    for i in range(n):
        nbrs = [j for j in range(n)
                if j > i and abs(C[j] - C[i]) < R[i] + R[j]]
        for a, j in enumerate(nbrs):
            for k in nbrs[a + 1:]:
                if abs(C[k] - C[j]) >= R[j] + R[k]:
                    continue
                A = np.array([[2 * (C[j] - C[i]).real, 2 * (C[j] - C[i]).imag],
                              [2 * (C[k] - C[i]).real, 2 * (C[k] - C[i]).imag]])
                det = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
                if abs(det) < 1e-12:
                    continue
                x = np.linalg.solve(A, [kap[j] - kap[i], kap[k] - kap[i]])
                v = complex(x[0], x[1])
                co = Binv @ np.array([v.real, v.imag])
                co -= np.floor(co)                    # reduce mod Lambda
                v = co[0] * v1 + co[1] * v2
                E = float(np.max(W - np.abs(v - C) ** 2))
                if 0 < E < best:
                    best, argbest, tri = E, v, (i, j, k)
    return float(np.sqrt(best)), argbest, tri


if __name__ == '__main__':
    import json
    data = json.load(open('fig_data.json'))
    for key in ('4_1', '7_4'):
        d = data[key]
        v1, v2 = complex(*d['v1']), complex(*d['v2'])
        base = [(complex(*c['center']), c['radius']) for c in d['cand']
                if c['radius'] >= 0.9 * d['y0']]
        circ = [(c + i * v1 + j * v2, r) for c, r in base
                for i in (-1, 0, 1) for j in (-1, 0, 1)]
        ded = {}
        for c, r in circ:
            ded.setdefault((round(c.real, 8), round(c.imag, 8),
                            round(r, 8)), (c, r))
        circ = list(ded.values())
        for th in (1.0, 0.999):
            y0, x, _ = exact_floor(circ, v1, v2, theta=th)
            print(f"{key} theta={th}: Y0 = {y0:.12f}  x* = {x:.8f}")
