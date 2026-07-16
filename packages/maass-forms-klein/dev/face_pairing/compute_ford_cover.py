"""
Recompute covering generators for Ford-type fundamental domains of knot
complement groups, replicating the algorithm in
maass_forms_klein.hyperbolic_space.utils.find_covering_generators
(numpy double precision).
"""
import numpy as np
import snappy
import itertools, json, sys

EPS = 1e-10

def sl2c(G, word):
    m = G.SL2C(word)
    return np.array([[complex(m[0,0]), complex(m[0,1])],
                     [complex(m[1,0]), complex(m[1,1])]], dtype=complex)

def normalise_group(name, prec_hp=False):
    M = snappy.ManifoldHP(name) if prec_hp else snappy.Manifold(name)
    G = M.fundamental_group()
    gens = {}
    for g in G.generators():
        gens[g] = sl2c(G, g)
        gens[g.upper()] = np.linalg.inv(gens[g])
    L = sl2c(G, G.longitude())
    Mm = sl2c(G, G.meridian())
    # fixed point of the (rank-2) cusp: use whichever parabolic has |c| away from 0
    P = L if abs(L[1,0]) >= abs(Mm[1,0]) else Mm
    a, b, c, d = P[0,0], P[0,1], P[1,0], P[1,1]
    if abs(c) < EPS:
        Q = np.eye(2, dtype=complex)      # cusp already at infinity
    else:
        z0 = (a - d) / (2*c)
        # Moebius map sending z0 -> infinity, det = 1
        Q = np.array([[0, -1], [1, -z0]], dtype=complex)
    Qinv = np.linalg.inv(Q)
    ngens = {k: Q @ g @ Qinv for k, g in gens.items()}
    Ln = Q @ L @ Qinv
    Mn = Q @ Mm @ Qinv
    # translations: parabolic upper triangular [[e, alpha],[0, 1/e]], z -> z + e*alpha
    assert abs(Ln[1,0]) < 1e-8 and abs(Mn[1,0]) < 1e-8, (Ln, Mn)
    tL = Ln[0,0] * Ln[0,1]
    tM = Mn[0,0] * Mn[0,1]
    # rescale so that the meridian translation equals 1:  z -> z / tM,
    # i.e. conjugate by diag(mu, 1/mu) with mu^2 = 1/tM.
    mu = 1/np.sqrt(tM)
    D = np.array([[mu, 0], [0, 1/mu]], dtype=complex)
    Dinv = np.linalg.inv(D)
    ngens = {k: D @ g @ Dinv for k, g in ngens.items()}
    Ln = D @ Ln @ Dinv
    Mn = D @ Mn @ Dinv
    tL, tM = tL/tM, 1.0 + 0j          # tL is now the cusp shape
    ngens['L'] = Ln; ngens['l'] = np.linalg.inv(Ln)
    ngens['M'] = Mn; ngens['m'] = np.linalg.inv(Mn)
    return M, G, ngens, tL, tM

def circle_of(mat):
    c = mat[1,0]
    if abs(c) < EPS:
        return None
    return (-mat[1,1]/c, 1.0/abs(c))   # center, radius

# ---- lattice reduction of circle centres --------------------------------
def lattice_reduce(center, v1, v2):
    """Translate center by lattice Z v1 + Z v2 to representative in the base cell
    (coords in [0,1)); return (new_center, (n1, n2))."""
    B = np.array([[v1.real, v2.real], [v1.imag, v2.imag]])
    x = np.linalg.solve(B, np.array([center.real, center.imag]))
    n = np.floor(x)
    nc = center - n[0]*v1 - n[1]*v2
    return nc, (-int(n[0]), -int(n[1]))

# ---- covering test: exact quadtree (disc is convex => cell in disc iff corners in) ----
def cell_in_some_circle(corners, circles, scale):
    for (cc, r) in circles:
        rr = (r*scale)**2
        if all(abs(z-cc)**2 <= rr for z in corners):
            return True
    return False

def parallelogram_covered(base, v1, v2, circles, scale=0.999, min_frac=1/512):
    """Quadtree check that base + [0,1)v1 + [0,1)v2 is covered by scaled discs.
    Circles are pruned per node (disc must meet the cell bounding box)."""
    carr = np.array([c for c, _ in circles])
    rarr = np.array([r*scale for _, r in circles])
    stack = [(base, v1, v2, 1.0, np.arange(len(circles)))]
    while stack:
        b, w1, w2, frac, idx = stack.pop()
        if len(idx) == 0:
            return False, b + 0.5*(w1+w2)
        corners = np.array([b, b+w1, b+w2, b+w1+w2])
        # cell covered by a single disc iff all 4 corners inside it (convexity)
        d2 = np.abs(corners[None, :] - carr[idx][:, None])**2
        if np.any(np.all(d2 <= (rarr[idx]**2)[:, None], axis=1)):
            continue
        if frac < min_frac:
            return False, b + 0.5*(w1+w2)
        h1, h2 = 0.5*w1, 0.5*w2
        for s in (b, b+h1, b+h2, b+h1+h2):
            # prune: keep discs whose center is within radius + cell circumradius
            mid = s + 0.5*(h1+h2)
            circum = 0.5*abs(h1+h2) + 0.5*abs(h1-h2)
            keep = idx[np.abs(carr[idx] - mid) <= rarr[idx] + circum]
            stack.append((s, h1, h2, frac*0.5, keep))
    return True, None

# ---- main search --------------------------------------------------------
def covering_search(name, maxlen=8, scale=0.999, min_frac=1/512, max_words=200000,
                    verbose=True):
    M, G, gens, tL, tM = normalise_group(name)
    v1, v2 = tL, tM
    # base cell centred at 0 like translation_fundamental_domain (base at -(v1+v2)/2)
    base = -(v1+v2)/2
    letters = [g for g in gens if g not in 'LlMm']
    # BFS over reduced words
    seen_circ = {}   # rounded (center, radius) -> word  (after lattice reduction)
    frontier = [(g,) for g in letters]
    mats = {(g,): gens[g] for g in letters}
    all_circles = []   # (word_str_with_translation, center, radius)
    lat_diam = max(abs(v1), abs(v2))
    for length in range(1, maxlen+1):
        new_frontier = []
        for w in frontier:
            mat = mats[w]
            circ = circle_of(mat)
            if circ is not None:
                c0, r = circ
                if r > 1e-3:  # ignore microscopic spheres
                    cred, (n1, n2) = lattice_reduce(c0, v1, v2)
                    key = (round(cred.real, 6), round(cred.imag, 6), round(r, 6))
                    if key not in seen_circ:
                        word = ''.join(w)
                        seen_circ[key] = (word, (n1, n2))
                        all_circles.append((word, (n1, n2), cred, r))
            if length < maxlen:
                for g in letters:
                    if w[-1] != g.swapcase():
                        w2 = w + (g,)
                        mats[w2] = mat @ gens[g]
                        new_frontier.append(w2)
        # free memory of old level
        for w in frontier:
            mats.pop(w, None)
        frontier = new_frontier
        if len(seen_circ) > max_words:
            break
        # try to certify coverage with circles found so far (incl. neighbour translates)
        cand = expand_translates(all_circles, v1, v2)
        cand = [x for x in cand if circle_meets_cell(x[2], x[3], base, v1, v2)]
        cand.sort(key=lambda x: -x[3])
        ok, _ = parallelogram_covered(base, v1, v2, [(x[2], x[3]) for x in cand],
                                      scale=scale, min_frac=min_frac)
        if verbose:
            print(f"  len<={length}: {len(seen_circ)} distinct circles, "
                  f"{len(cand)} near cell, covered={ok}", flush=True)
        if ok:
            small = build_cover_greedy(cand, base, v1, v2, scale, min_frac)
            reduced = greedy_reduce(small, base, v1, v2, scale, min_frac)
            y0 = floor_height(base, v1, v2, [(x[2], x[3]) for x in cand])
            y0r = floor_height(base, v1, v2, [(x[2], x[3]) for x in reduced])
            return dict(name=name, v1=v1, v2=v2, base=base, maxlen=length,
                        n_circles=len(seen_circ), cover=reduced, group=G,
                        y0=y0, y0_reduced=y0r, cand=cand)
    return dict(name=name, v1=v1, v2=v2, base=base, maxlen=maxlen,
                n_circles=len(seen_circ), cover=None, group=G)

def expand_translates(all_circles, v1, v2, k=1):
    out = []
    for word, (n1, n2), c, r in all_circles:
        for i in range(-k, k+1):
            for j in range(-k, k+1):
                out.append((word, (n1+i, n2+j), c + i*v1 + j*v2, r))
    return out

def circle_meets_cell(c, r, base, v1, v2, margin=1.0):
    # cheap bounding-box test
    xs = [base.real, (base+v1).real, (base+v2).real, (base+v1+v2).real]
    ys = [base.imag, (base+v1).imag, (base+v2).imag, (base+v1+v2).imag]
    return (min(xs)-r <= c.real <= max(xs)+r) and (min(ys)-r <= c.imag <= max(ys)+r)

def build_cover_greedy(cand, base, v1, v2, scale, min_frac):
    """Build a small cover: repeatedly take an uncovered witness point and add
    the largest disc that contains it comfortably."""
    cand = sorted(cand, key=lambda x: -x[3])
    selected = []
    used = set()
    while True:
        ok, witness = parallelogram_covered(
            base, v1, v2, [(y[2], y[3]) for y in selected] or [(base, 1e-12)],
            scale=scale, min_frac=min_frac)
        if ok:
            return selected
        best, best_score = None, -1
        for i, (w, t, c, r) in enumerate(cand):
            if i in used:
                continue
            margin = r*scale - abs(witness - c)
            if margin > 0:
                score = margin + r        # prefer big discs that contain it well
                if score > best_score:
                    best, best_score = i, score
        if best is None:
            # witness not inside any candidate: give up, return everything
            return cand
        used.add(best)
        selected.append(cand[best])

def greedy_reduce(cand, base, v1, v2, scale, min_frac):
    circles = list(cand)
    circles.sort(key=lambda x: x[3])   # remove small first
    keep = list(circles)
    for x in circles:
        if len(keep) == 1:
            break
        trial = [y for y in keep if y is not x]
        ok, _ = parallelogram_covered(base, v1, v2, [(y[2], y[3]) for y in trial],
                                      scale=scale, min_frac=min_frac)
        if ok:
            keep = trial
    keep.sort(key=lambda x: (-x[3], len(x[0])))
    return keep

def floor_height(base, v1, v2, circles, n=400):
    """min over the cell of  max_k sqrt(r_k^2 - |z - c_k|^2)_+ :
    the height of the lowest point of the Ford domain floor (approx)."""
    s = np.linspace(0, 1, n, endpoint=False) + 0.5/n
    u, w = np.meshgrid(s, s)
    z = base + u*v1 + w*v2
    h2 = np.zeros(z.shape)
    for c, r in circles:
        h2 = np.maximum(h2, r*r - np.abs(z - c)**2)
    return float(np.sqrt(np.maximum(h2, 0).min()))

def translation_word(n1, n2):
    w = ''
    w += ('L' if n1 > 0 else 'l') * abs(n1)
    w += ('M' if n2 > 0 else 'm') * abs(n2)
    return w

if __name__ == '__main__':
    knots = sys.argv[1:] or ['4_1']
    results = {}
    for k in knots:
        print(f"=== {k} ===", flush=True)
        try:
            import os
            res = covering_search(k, maxlen=int(os.environ.get('MAXLEN', 8)))
        except Exception as e:
            print(f"  FAILED: {e}")
            continue
        if res['cover'] is None:
            print(f"  no cover up to length {res['maxlen']} ({res['n_circles']} circles)")
            results[k] = dict(v1=str(res['v1']), v2=str(res['v2']), cover=None,
                              n_circles=res['n_circles'], maxlen=res['maxlen'])
            continue
        cov = res['cover']
        print(f"  v1={res['v1']:.6f} v2={res['v2']:.6f}  "
              f"Y0={res['y0']:.6f} (reduced: {res['y0_reduced']:.6f})")
        print(f"  covering set: {len(cov)} spheres, max word length "
              f"{max(len(w) for w,_,_,_ in cov)}, radii "
              f"{sorted(set(round(r,6) for _,_,_,r in cov), reverse=True)}")
        for w, (n1, n2), c, r in cov:
            print(f"    {w+translation_word(n1,n2):<20s} centre {c:.6f}  radius {r:.6f}")
        results[k] = dict(v1=[res['v1'].real, res['v1'].imag],
                          v2=[res['v2'].real, res['v2'].imag],
                          base=[res['base'].real, res['base'].imag],
                          y0=res['y0'], y0_reduced=res['y0_reduced'],
                          n_circles=res['n_circles'], maxlen=res['maxlen'],
                          cover=[dict(word=w, trans=[n1,n2],
                                      center=[c.real, c.imag], radius=r)
                                 for w, (n1, n2), c, r in cov])
    json.dump(results, open('/home/claude/ford/results.json', 'a'), indent=1)
