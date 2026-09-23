r"""
Hejhal two-height eigenvalue search via aligned SVD null vectors.

This module folds the eigenvalue-location machinery of the 2026 survey
(``LESSONS_LEARNED_EIGENVALUE_SURVEY.md`` §1) into the package.  It is lifted
from the battle-tested pure-numpy reference
``reference/ford/maass_eigenvalue_runner.py`` (run end-to-end on all 28 knots
up to 8 crossings and cross-validated between machines to ``1e-7``) and its
exceptional-interval driver ``reference/ford/exceptional_scan.py``.

Why this replaces the ``modform/search.py`` sign-change locator (§1.1):

* The old path normalises ``c(1) = 1``, deletes a row/column, solves the square
  system at two heights and looks for sign changes of ``Re(C0[i] - C1[i])``.
  For Kleinian groups this is fragile: which coefficient is reliably nonzero
  varies by knot and eigenvalue, sign changes are lost when the indicator dips
  without crossing, and symmetric forms with ``c(1) = 0`` are missed entirely.
* The reference-tested path builds a matrix at each height ``Y1 > Y2``, then
  row-equilibrates and column-normalises each before taking its SVD. The
  indicator is the aligned distance of the two smallest right singular
  vectors; the null vectors provide coefficient vectors without pinning
  ``c(1)``. Dips are located by golden-section refinement.

The survey prose requests bare ``sigma_min`` of a stacked system, but the
named reference runner does not use it. A literal stacked implementation
failed the 4_1 negative control; see ``dev/EIGENVALUE_SURVEY_INTEGRATION_NOTES.md``.

Every dip is re-computed with independent parameters (``digits + 0.75``,
different height factors) and triaged by first-pass residual (§1.3):
``resid >= 0.10`` rejected, ``0.02 <= resid < 0.10`` reported as a candidate,
``resid < 0.02`` fully re-scanned.  Strict confirmation is
``spread < 1e-3`` AND ``resid < 1e-3``.

The imaginary/real order switch (``order='imaginary'|'real'``) selects the main
window (``s = 1 + iR``) or the exceptional interval (``s = 1 + t`` real); it is
threaded straight into :class:`maass_forms_klein.functions.besselk_quad.KBessel`
(§2.2).

.. NOTE::

    Covering matrices remain lifted from the self-contained NumPy reference.
    Sampling heights now use the package's power-vertex Ford floor when the
    horoball-face path succeeds, or the exact floor of the certified emitted
    cover otherwise. The origin offset between frames prevents mixing sphere
    centres without an explicit alignment.
"""

import logging

import numpy as np

from maass_forms_klein.functions.besselk_quad import KBessel
from maass_forms_klein.modform.survey_floor import certified_floor_height, exact_cover_floor

log = logging.getLogger(__name__)

EPS = 1e-10


# ---------------------------------------------------------------------------
# group setup (cusp-normalised holonomy, meridian scaled to 1)
# ---------------------------------------------------------------------------
def normalise_group(name, simplify=True):
    r"""
    Cusp-normalised holonomy of a snappy manifold.

    Returns ``(gens, tau, volume)`` where ``gens`` maps generator letters (and
    their inverses) to ``2x2`` complex matrices with the meridian translation
    scaled to ``1`` and the longitude to the cusp shape ``tau``.

    EXAMPLES::

        sage: from maass_forms_klein.modform.hejhal import normalise_group
        sage: gens, tau, vol = normalise_group('4_1')
        sage: bool(abs(tau.imag) - 3.4641016 < 1e-5)
        True
        sage: bool(abs(vol - 2.0298832) < 1e-6)
        True
    """
    import snappy

    M = snappy.Manifold(name)
    G = M.fundamental_group(simplify_presentation=simplify)

    def sl2c(word):
        m = G.SL2C(word)
        return np.array(
            [
                [complex(m[0, 0]), complex(m[0, 1])],
                [complex(m[1, 0]), complex(m[1, 1])],
            ],
            dtype=complex,
        )

    gens = {}
    for g in G.generators():
        gens[g] = sl2c(g)
        gens[g.upper()] = np.linalg.inv(gens[g])
    L = sl2c(G.longitude())
    Mm = sl2c(G.meridian())
    P = L if abs(L[1, 0]) >= abs(Mm[1, 0]) else Mm
    a, _b, c, d = P[0, 0], P[0, 1], P[1, 0], P[1, 1]
    if abs(c) < EPS:
        Q = np.eye(2, dtype=complex)
    else:
        z0 = (a - d) / (2 * c)
        Q = np.array([[0, -1], [1, -z0]], dtype=complex)
    Qi = np.linalg.inv(Q)
    ngens = {k: Q @ g @ Qi for k, g in gens.items()}
    Ln, Mn = Q @ L @ Qi, Q @ Mm @ Qi
    assert abs(Ln[1, 0]) < 1e-8 and abs(Mn[1, 0]) < 1e-8
    tL = Ln[0, 0] * Ln[0, 1]
    tM = Mn[0, 0] * Mn[0, 1]
    mu = 1 / np.sqrt(tM)
    D = np.array([[mu, 0], [0, 1 / mu]], dtype=complex)
    Di = np.linalg.inv(D)
    ngens = {k: D @ g @ Di for k, g in ngens.items()}
    tau = tL / tM
    return ngens, tau, float(M.volume())


def reduced_basis(tau):
    r"""
    Lagrange--Gauss reduced basis of the cusp lattice ``Z*tau + Z``.

    Spans the same lattice as ``(tau, 1)`` by a unimodular change of basis, so
    the Fourier expansion, pull-back and covering cell are unchanged; skew cusp
    shapes get dramatically smaller dual boxes (survey §2.3, §3.2).

    EXAMPLES::

        sage: from maass_forms_klein.modform.hejhal import reduced_basis
        sage: v2, v1 = reduced_basis(3.4641016151j)
        sage: bool(abs(v2) > 3.0 and abs(v1 - 1.0) < 1e-9)
        True
    """
    v1, v2 = complex(tau), 1.0 + 0j
    for _ in range(100):
        if abs(v2) < abs(v1):
            v1, v2 = v2, v1
        mu = round((v2.real * v1.real + v2.imag * v1.imag) / abs(v1) ** 2)
        if mu == 0:
            break
        v2 = v2 - mu * v1
    return v2, v1  # longer vector first, matching (tau, 1) ordering


# ---------------------------------------------------------------------------
# covering generators (BFS + certified cover; as in the paper section)
# ---------------------------------------------------------------------------
def _circle(mat):
    c = mat[1, 0]
    if abs(c) < EPS:
        return None
    return (-mat[1, 1] / c, 1.0 / abs(c))


def _lattice_reduce(center, v1, v2):
    B = np.array([[v1.real, v2.real], [v1.imag, v2.imag]])
    x = np.linalg.solve(B, [center.real, center.imag])
    n = np.floor(x)
    return center - n[0] * v1 - n[1] * v2


def _covered(base, v1, v2, circles, scale=0.999, min_frac=1 / 512):
    """Certified quadtree cover test; returns ``(bool, witness_or_None)``."""
    if not circles:
        return False, base + 0.5 * (v1 + v2)
    carr = np.array([c for c, _ in circles])
    rarr = np.array([r * scale for _, r in circles])
    stack = [(base, v1, v2, 1.0, np.arange(len(circles)))]
    while stack:
        b, w1, w2, frac, idx = stack.pop()
        if len(idx) == 0:
            return False, b + 0.5 * (w1 + w2)
        corners = np.array([b, b + w1, b + w2, b + w1 + w2])
        d2 = np.abs(corners[None, :] - carr[idx][:, None]) ** 2
        if np.any(np.all(d2 <= (rarr[idx] ** 2)[:, None], axis=1)):
            continue
        if frac < min_frac:
            return False, b + 0.5 * (w1 + w2)
        h1, h2 = 0.5 * w1, 0.5 * w2
        for s in (b, b + h1, b + h2, b + h1 + h2):
            mid = s + 0.5 * (h1 + h2)
            circum = 0.5 * abs(h1 + h2) + 0.5 * abs(h1 - h2)
            keep = idx[np.abs(carr[idx] - mid) <= rarr[idx] + circum]
            stack.append((s, h1, h2, frac * 0.5, keep))
    return True, None


def covering_data(gens, tau, maxlen=14, verbose=False, extra_levels=3):
    r"""
    Certified covering matrices and a sampled floor estimate for the cusp cell.

    Returns ``(cover_mats, floor_estimate)``: matrices whose isometric
    hemispheres (with lattice translations applied) cover the centred cusp
    cell, and a grid estimate of the floor of the *selected* subset.

    See the reference runner for the deepening/greedy-reduction rationale.
    HejhalContext does not use the grid estimate for sampling heights.
    """
    v1, v2 = reduced_basis(tau)
    best = None
    first_len = None
    flat = 0
    base = -(v1 + v2) / 2
    letters = list(gens)

    def _row_key(m):
        t = (
            round(m[1, 0].real, 5) + 0.0,
            round(m[1, 0].imag, 5) + 0.0,
            round(m[1, 1].real, 5) + 0.0,
            round(m[1, 1].imag, 5) + 0.0,
        )
        return min(t, tuple(-x + 0.0 for x in t))

    frontier = [(g,) for g in letters]
    mats = {(g,): gens[g] for g in letters}
    seen_rows = {_row_key(gens[g]) for g in letters}
    circ_words = {}
    for length in range(1, maxlen + 1):
        for w in frontier:
            circ = _circle(mats[w])
            if circ is not None and circ[1] > 5e-3:
                cred = _lattice_reduce(circ[0], v1, v2)
                key = (round(cred.real, 6), round(cred.imag, 6), round(circ[1], 6))
                if key not in circ_words:
                    t = circ[0] - cred
                    T = np.array([[1, t], [0, 1]], dtype=complex)
                    circ_words[key] = (w, mats[w] @ T, cred, circ[1])
        cand = []
        for _w, m, c, rr in circ_words.values():
            for i in (-1, 0, 1):
                for j in (-1, 0, 1):
                    cc = c + i * v1 + j * v2
                    T = np.array([[1, -(i * v1 + j * v2)], [0, 1]], dtype=complex)
                    cand.append((m @ T, cc, rr))
        xs = [p.real for p in (base, base + v1, base + v2, base + v1 + v2)]
        ys = [p.imag for p in (base, base + v1, base + v2, base + v1 + v2)]
        cand = [
            t
            for t in cand
            if (
                min(xs) - t[2] <= t[1].real <= max(xs) + t[2]
                and min(ys) - t[2] <= t[1].imag <= max(ys) + t[2]
            )
        ]
        cand.sort(key=lambda t: -t[2])
        ok, _ = _covered(base, v1, v2, [(c, r) for _, c, r in cand])
        if ok:
            n = 400
            s = np.linspace(0, 1, n, endpoint=False) + 0.5 / n
            u, w2m = np.meshgrid(s, s)
            Z = base + u * v1 + w2m * v2
            h2 = np.zeros(Z.shape)
            big = [t for t in cand if t[2] >= 0.2 * cand[0][2]]
            for _, c, r in big:
                h2 = np.maximum(h2, r * r - np.abs(Z - c) ** 2)
            Y0 = float(np.sqrt(max(h2.min(), 0.0)))
            selected = []
            while True:
                okk, wit = _covered(base, v1, v2, [(c, r) for _, c, r in selected])
                if okk:
                    break
                bt, score = None, -1.0
                for t in cand:
                    margin = t[2] * 0.999 - abs(wit - t[1])
                    if margin > 0 and margin + t[2] > score:
                        bt, score = t, margin + t[2]
                if bt is None:
                    selected = cand
                    break
                selected.append(bt)
            h2s = np.zeros(Z.shape)
            for _, c, r in selected:
                h2s = np.maximum(h2s, r * r - np.abs(Z - c) ** 2)
            Y0s = float(np.sqrt(max(h2s.min(), 0.0)))
            rest = [t for t in cand if not any(t is sel for sel in selected)]
            while Y0s < 0.98 * Y0 and rest and len(selected) < 40:
                imin = np.unravel_index(np.argmin(h2s), h2s.shape)
                zmin = Z[imin]
                depth = [(t[2] ** 2 - abs(zmin - t[1]) ** 2, k) for k, t in enumerate(rest)]
                dmax, kmax = max(depth)
                if dmax <= h2s[imin] + 1e-12:
                    break
                t = rest.pop(kmax)
                selected.append(t)
                h2s = np.maximum(h2s, t[2] ** 2 - np.abs(Z - t[1]) ** 2)
                Y0s = float(np.sqrt(max(h2s.min(), 0.0)))
            if verbose:
                log.debug(
                    "  cover at word length %d: %d spheres, reduced to %d;"
                    " Y0 = %.4f (full-list floor %.4f)",
                    length,
                    len(cand),
                    len(selected),
                    Y0s,
                    Y0,
                )
            prev = best[1] if best is not None else 0.0
            if best is None or Y0s > best[1]:
                best = ([m for m, _, _ in selected], Y0s)
            if first_len is None:
                first_len = length
            flat = flat + 1 if best[1] < 1.01 * prev else 0
            if (
                length >= first_len + extra_levels
                or length >= maxlen
                or flat >= 2
                or len(frontier) > 400000
            ):
                return best
        if length < maxlen:
            nf = []
            for w in frontier:
                mat = mats[w]
                for g in letters:
                    if w[-1] != g.swapcase():
                        m2 = mat @ gens[g]
                        key = _row_key(m2)
                        if key in seen_rows:
                            continue
                        seen_rows.add(key)
                        w2t = (*w, g)
                        mats[w2t] = m2
                        nf.append(w2t)
            for w in frontier:
                mats.pop(w, None)
            frontier = nf
            if len(nf) > 1200000 or len(seen_rows) > 6000000:
                break
    if best is not None:
        return best
    raise RuntimeError("no covering set found up to word length %d" % maxlen)


# ---------------------------------------------------------------------------
# pull-back into the Ford domain (centred cell, §3.3)
# ---------------------------------------------------------------------------
def pullback_points(xs, Y, cover, tau, max_iter=200):
    r"""
    Pull each ``(x, Y)`` back into the Ford domain; return ``(x*, y*)`` arrays.

    Reduction is into the CENTRED cell (``x`` in ``[-1/2, 1/2)`` in lattice
    coordinates), the same convention the covering certification uses (§3.3).
    """
    v1, v2 = reduced_basis(tau)
    B = np.array([[v1.real, v2.real], [v1.imag, v2.imag]])
    Binv = np.linalg.inv(B)
    xs_out = np.empty(len(xs), dtype=complex)
    ys_out = np.empty(len(xs), dtype=float)
    for k, x0 in enumerate(xs):
        x, y = complex(x0), float(Y)
        for _ in range(max_iter):
            co = Binv @ np.array([x.real, x.imag])
            co -= np.floor(co + 0.5)
            x = co[0] * v1 + co[1] * v2
            moved = False
            for m in cover:
                a, b = m[0, 0], m[0, 1]
                c, d = m[1, 0], m[1, 1]
                w = c * x + d
                D = (w.real**2 + w.imag**2) + abs(c) ** 2 * y * y
                if D < 1.0 - 1e-13:
                    x = ((a * x + b) * w.conjugate() + a * c.conjugate() * y * y) / D
                    y = y / D
                    moved = True
                    break
            if not moved:
                break
        else:
            raise RuntimeError("pull-back did not terminate")
        xs_out[k] = x
        ys_out[k] = y
    return xs_out, ys_out


# ---------------------------------------------------------------------------
# linear system, indicator and scan
# ---------------------------------------------------------------------------
def _require_pullback_above_floor(heights, floor_height):
    """Reject incomplete or misplaced covers before building the linear system."""
    heights = np.asarray(heights, dtype=float)
    if heights.size == 0 or not np.all(np.isfinite(heights)):
        raise RuntimeError("pull-back produced empty or non-finite heights")
    if float(np.min(heights)) < 0.98 * floor_height:
        raise RuntimeError(
            f"pull-back minimum height {float(np.min(heights)):.6g} is below 0.98 "
            f"times the exact floor {floor_height:.6g}; covering set or frame misplaced"
        )


class HejhalContext:
    r"""
    Two-height Hejhal system for a knot complement, with the aligned-null-vector
    indicator from the validated reference runner.

    INPUT:

    - ``name`` -- a snappy manifold identifier (e.g. ``'4_1'``).
    - ``digits`` -- target truncation digits driving the matrix size
      (default: ``2.5``).
    - ``rmax`` -- largest spectral parameter used to size the truncation
      (default: ``9.0``).
    - ``yfacs`` -- ``(Y1, Y2)`` sampling heights as fractions of the floor
      ``Y0`` (default: ``(0.95, 0.88)``).
    - ``order`` -- ``'imaginary'`` (main window) or ``'real'`` (exceptional
      interval); threaded into every :class:`KBessel` (default:
      ``'imaginary'``).
    - ``verbose`` -- print covering / lattice diagnostics (default: ``False``).

    EXAMPLES:

    Build the figure-eight context and check the indicator dips at the first
    strict eigenvalue ``R = 4.90008537`` (survey §6)::

        sage: from maass_forms_klein.modform.hejhal import HejhalContext
        sage: ctx = HejhalContext('4_1', rmax=6.5)          # long time
        sage: bool(ctx.indicator(4.90008537) < 1e-3)        # long time
        True
        sage: bool(ctx.indicator(4.6) > 1e-2)               # long time
        True
    """

    @staticmethod
    def _compute_cover(name, gens, tau, verbose):
        try:
            return covering_data(gens, tau, verbose=verbose)
        except RuntimeError:
            pass
        if verbose:
            log.debug("  covering BFS failed; retrying with geometric presentation")
        gens_g, tau_g, _ = normalise_group(name, simplify=False)
        if abs(tau_g - tau) < 1e-4:
            try:
                return covering_data(gens_g, tau_g, maxlen=10, extra_levels=2, verbose=verbose)
            except RuntimeError:
                pass
        elif verbose:
            log.debug("  geometric frame mismatch (tau %s vs %s)", tau_g, tau)
        for _ in range(2):
            gens2, tau2, _ = normalise_group(name)
            if abs(tau2 - tau) < 1e-9:
                try:
                    return covering_data(gens2, tau2, maxlen=12, verbose=verbose)
                except RuntimeError:
                    continue
        raise RuntimeError(f"covering failed for {name} in all presentations")

    def __init__(
        self,
        name,
        digits=2.5,
        rmax=9.0,
        yfacs=(0.95, 0.88),
        mextra=0,
        order="imaginary",
        verbose=False,
    ):
        self.name = name
        gens, tau, vol = normalise_group(name)
        self.tau, self.vol = tau, vol
        cover, cover_floor_estimate = self._compute_cover(name, gens, tau, verbose)
        self.cover, self.cover_floor_estimate = cover, cover_floor_estimate
        v1, v2 = reduced_basis(tau)
        try:
            self.Y0 = certified_floor_height(name)
            self.floor_source = "full_ford"
        except ValueError as exc:
            # A horoball-frame origin offset can make the package's face
            # auto-bound fail even when this centred-cell cover is certified.
            # Do not mix its centres with holonomy matrices: take only the
            # scalar power-vertex floor of the emitted cover instead.
            self.Y0 = exact_cover_floor(cover, v1, v2)
            self.floor_source = "certified_cover"
            log.warning("%s: full Ford floor unavailable (%s); using exact cover floor", name, exc)
        B = np.array([[v1.real, v1.imag], [v2.real, v2.imag]])
        Wd = np.linalg.inv(B)
        self.w1 = complex(Wd[0, 0], Wd[1, 0])
        self.w2 = complex(Wd[0, 1], Wd[1, 1])
        self.Y1, self.Y2 = yfacs[0] * self.Y0, yfacs[1] * self.Y0
        # truncation from the true K-Bessel decay exponent eta (survey §2.3)
        R = rmax
        target = digits * np.log(10)
        lo, hi = R * 1.0001, 8.0 * R
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            eta = np.sqrt(mid**2 - R**2) - R * np.arccosh(mid / R)
            if eta < target:
                lo = mid
            else:
                hi = mid
        Mnorm = hi / (2 * np.pi * self.Y2)
        cutoff = Mnorm + mextra * min(abs(self.w1), abs(self.w2))
        vecs, idx = [], []
        w1a2 = abs(self.w1) ** 2
        Wspan = int(np.ceil(cutoff / abs(self.w1))) + 1
        k2 = 0
        while True:
            found = False
            for s in [0] if k2 == 0 else [k2, -k2]:
                sv = s * self.w2
                k1c = -(self.w1.real * sv.real + self.w1.imag * sv.imag) / w1a2
                for k1 in range(int(np.floor(k1c)) - Wspan, int(np.ceil(k1c)) + Wspan + 1):
                    v = k1 * self.w1 + sv
                    if abs(v) <= cutoff:
                        vecs.append(v)
                        idx.append((k1, s))
                        found = True
            if not found and k2 > 0:
                break
            k2 += 1
        self.vecs = np.array(vecs)
        self.index = idx
        self.M1 = max(abs(k1) for k1, _ in idx)
        self.M2 = max(abs(k2v) for _, k2v in idx)
        self.zero = idx.index((0, 0))
        order_ = np.argsort(np.abs(self.vecs))
        nz = [i for i in order_ if i != self.zero]
        self.ref = nz[0]
        self.tests = nz[1:5]
        Q1, Q2 = self.M1 + 1, self.M2 + 1
        m1 = np.arange(2 * Q1) / (2.0 * Q1)
        m2 = np.arange(2 * Q2) / (2.0 * Q2)
        g1, g2 = np.meshgrid(m1, m2, indexing="ij")
        self.xs = (g1 * v1 + g2 * v2).ravel()
        self.N = len(self.xs)
        if verbose:
            log.debug(
                "  lattice tau = %s; dual box %dx%d -> %d vectors,"
                " %d sample points; Y1=%.4f Y2=%.4f",
                tau,
                self.M1,
                self.M2,
                len(self.vecs),
                self.N,
                self.Y1,
                self.Y2,
            )
        self.pre = {}
        for Y in (self.Y1, self.Y2):
            xstar, ystar = pullback_points(self.xs, Y, cover, tau)
            _require_pullback_above_floor(ystar, self.Y0)
            unmoved = int(np.sum(np.abs(ystar - Y) < 1e-12))
            if unmoved > max(2, 0.01 * len(ystar)):
                raise RuntimeError(
                    f"pull-back failed: {unmoved}/{len(ystar)} sample points at "
                    f"height Y={Y:.4f} < Y0={self.Y0:.4f} are outside every "
                    f"covering hemisphere - covering set misplaced?"
                )
            absv = np.abs(self.vecs)
            args_pb = 2 * np.pi * np.outer(absv, ystar)
            args_diag = 2 * np.pi * absv * Y

            def bpair(a, b):
                return a.real * b.real + a.imag * b.imag

            ph_star = np.exp(2j * np.pi * bpair(self.vecs[:, None], xstar[None, :]))
            ph_samp = np.exp(-2j * np.pi * bpair(self.vecs[:, None], self.xs[None, :]))
            kb = KBessel(np.concatenate([args_pb.ravel(), args_diag]), order=order)
            self.pre[Y] = {
                "kb": kb,
                "ystar": ystar,
                "ph_star": ph_star,
                "ph_samp": ph_samp,
                "nargs_pb": args_pb.size,
            }
        self.order = order

    def set_order(self, order):
        r"""
        Switch every :class:`KBessel` between imaginary and real order, so a
        context built for the main window can drive the exceptional scan.

        EXAMPLES::

            sage: from maass_forms_klein.modform.hejhal import HejhalContext
            sage: ctx = HejhalContext('4_1', rmax=2.2)      # long time
            sage: ctx.set_order('real'); ctx.order          # long time
            'real'
        """
        for p in self.pre.values():
            p["kb"].set_order(order)
        self.order = order

    def coefficients(self, R, Y):
        r"""
        Coefficient vector as the smallest right singular vector of the
        row-equilibrated + column-normalised cuspidal system (survey §1.1).

        Robust to eigenfunctions with vanishing low-order coefficients: no
        forced normalisation is made up front.
        """
        p = self.pre[Y]
        D, N = len(self.vecs), self.N
        kall = p["kb"].eval(R)
        K_pb = kall[: p["nargs_pb"]].reshape(D, N)
        K_diag = kall[p["nargs_pb"] :]
        P = (p["ystar"][None, :] * K_pb) * p["ph_star"]
        V = (p["ph_samp"] @ P.T) / N
        V[np.arange(D), np.arange(D)] -= Y * K_diag
        keep = [i for i in range(D) if i != self.zero]  # c(0)=0
        A = V[np.ix_(keep, keep)]
        rsc = np.abs(A).max(axis=1)
        rsc[rsc == 0] = 1.0
        A = A / rsc[:, None]
        csc = np.linalg.norm(A, axis=0)
        csc[csc == 0] = 1.0
        A = A / csc[None, :]
        _, _s, Vh = np.linalg.svd(A)
        chat = Vh[-1].conj()
        c = np.zeros(D, dtype=complex)
        c[keep] = chat / csc
        sel = np.argsort(np.abs(self.vecs[keep]))[:40]
        piv = keep[int(sel[np.argmax(np.abs(c[np.array(keep)[sel]]))])]
        c /= c[piv]
        return c

    def indicator(self, R):
        r"""
        Aligned two-height distance of the null coefficient vectors, restricted
        to the smallest-``|v|`` components.  Dips to truncation level at a
        genuine eigenvalue.
        """
        c1 = self.coefficients(R, self.Y1)
        c2 = self.coefficients(R, self.Y2)
        sel = np.argsort(np.abs(self.vecs))[: min(40, len(self.vecs))]
        a, b = c1[sel], c2[sel]
        beta = np.vdot(b, a) / np.vdot(b, b)
        return float(np.linalg.norm(a - beta * b) / np.linalg.norm(a))


def scan_and_refine(ctx, rmin, rmax, step, tol=1e-8, verbose=False):
    r"""
    Grid scan of ``ctx.indicator`` plus golden-section refinement of every dip
    (local minimum below ``0.6``).  Returns ``(grid, vals, candidates)`` where
    ``candidates`` is a list of ``(r0, residual)`` pairs (survey §1.1).
    """
    grid = np.arange(rmin, rmax + step / 2, step)
    vals = [float(ctx.indicator(float(r))) for r in grid]
    vals = np.array(vals)
    candidates = []
    gr = (np.sqrt(5) - 1) / 2
    for i in range(1, len(grid) - 1):
        if vals[i] < vals[i - 1] and vals[i] < vals[i + 1] and vals[i] < 0.6:
            a, b = grid[i - 1], grid[i + 1]
            x1, x2 = b - gr * (b - a), a + gr * (b - a)
            f1, f2 = ctx.indicator(x1), ctx.indicator(x2)
            for _ in range(36):
                if f1 < f2:
                    b, x2, f2 = x2, x1, f1
                    x1 = b - gr * (b - a)
                    f1 = ctx.indicator(x1)
                else:
                    a, x1, f1 = x1, x2, f2
                    x2 = a + gr * (b - a)
                    f2 = ctx.indicator(x2)
                if b - a < tol:
                    break
            r0 = 0.5 * (a + b)
            d0 = float(ctx.indicator(r0))
            candidates.append((r0, d0))
            if verbose:
                log.debug("    dip near %.3f -> r = %.8f, resid = %.2e", grid[i], r0, d0)
    return grid, vals, candidates


def _validate(ctx2, cands, step, order="imaginary"):
    """Independent-parameter validation with residual triage (survey §1.3)."""
    ctx2.set_order(order)
    out = []
    for r0, resid in cands:
        if resid >= 0.10:
            out.append({"r": r0, "resid": resid, "status": "rejected"})
        elif resid >= 0.02:
            out.append({"r": r0, "resid": resid, "status": "candidate"})
        else:
            lo = max(r0 - 2 * step, 1e-3)
            _, _, c2 = scan_and_refine(ctx2, lo, r0 + 2 * step, step / 2, tol=1e-6)
            match = [x for x in c2 if abs(x[0] - r0) < 5 * step]
            if match:
                r1, res1 = min(match, key=lambda t: abs(t[0] - r0))
                spread = abs(r1 - r0)
                combined_resid = max(resid, res1)
                out.append(
                    {
                        "r": 0.5 * (r0 + r1),
                        "spread": spread,
                        "resid": combined_resid,
                        "status": (
                            "confirmed"
                            if spread < 1e-3 and combined_resid < 1e-3
                            else "unconfirmed"
                        ),
                    }
                )
            else:
                out.append({"r": r0, "resid": resid, "status": "unconfirmed"})
    return out


def search_eigenvalues(
    name,
    rmin,
    rmax,
    step=0.05,
    digits=2.5,
    order="imaginary",
    yfacs=(0.95, 0.88),
    rmax_ctx=None,
    validate=True,
    verbose=False,
):
    r"""
    Locate Maass-form eigenvalues with the reference-tested two-height
    aligned-null-vector indicator, independent-parameter validation and triage.

    INPUT:

    - ``name`` -- a snappy manifold identifier (e.g. ``'4_1'``).
    - ``rmin``, ``rmax`` -- scan range for the spectral parameter (``R`` for
      imaginary order, ``t`` for real order).
    - ``step`` -- grid step (default: ``0.05``).
    - ``digits`` -- truncation digits (default: ``2.5``).
    - ``order`` -- ``'imaginary'`` or ``'real'`` (default: ``'imaginary'``).
    - ``yfacs`` -- sampling height fractions of ``Y0`` (default: ``(0.95, 0.88)``).
    - ``rmax_ctx`` -- spectral parameter used to size the truncation; defaults
      to ``rmax`` (pass a fixed value, e.g. ``2.2``, for the exceptional scan).
    - ``validate`` -- run the independent-parameter validation pass
      (default: ``True``).

    OUTPUT: a list of dicts, one per dip, with keys ``r``, ``resid``,
    ``status`` (``'confirmed'``/``'candidate'``/``'rejected'``/
    ``'unconfirmed'``) and, when confirmed, ``spread``.  Strict confirmation is
    ``status == 'confirmed'`` with ``spread < 1e-3`` and ``resid < 1e-3``.

    EXAMPLES:

    The figure-eight knot's first strict eigenvalue near ``R = 4.9001``, found
    and confirmed from a narrow window (survey §6)::

        sage: from maass_forms_klein.modform.hejhal import search_eigenvalues
        sage: res = search_eigenvalues('4_1', 4.75, 5.05, step=0.05)   # long time
        sage: strict = [c for c in res if c['status'] == 'confirmed'   # long time
        ....:           and c.get('spread', 1) < 1e-3 and c['resid'] < 1e-3]
        sage: bool(strict and abs(strict[0]['r'] - 4.90008537) < 1e-5)  # long time
        True
    """
    rc = rmax_ctx if rmax_ctx is not None else rmax
    ctx = HejhalContext(name, digits=digits, rmax=rc, yfacs=yfacs, order=order, verbose=verbose)
    _, _, cands = scan_and_refine(ctx, rmin, rmax, step, verbose=verbose)
    if not validate:
        return [{"r": r0, "resid": d0, "status": "candidate"} for r0, d0 in cands]
    ctx2 = HejhalContext(name, digits=digits + 0.75, rmax=rc, yfacs=(0.92, 0.84), order=order)
    return _validate(ctx2, cands, step, order=order)
