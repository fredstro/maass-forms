r"""
Vectorised Gauss--Legendre evaluation of the K-Bessel function.

This module implements the single integral-representation evaluator that was
used throughout the 2026 eigenvalue survey (see
``LESSONS_LEARNED_EIGENVALUE_SURVEY.md`` §2).  Rather than a scalar call per
argument (as in :func:`maass_forms_klein.modform.utils.bessel_function`), a
single :class:`KBessel` instance is built once from the whole argument set of a
Hejhal system and evaluated vectorised in ``float64`` chunks -- the form the
survey's matrix builder relies on.

Two orders are supported through the ``order`` constructor flag (or
:meth:`KBessel.set_order`):

- ``order='imaginary'`` (default) evaluates the *scaled* Bessel function of
  imaginary order,

  .. MATH::

      \widetilde{K}_{iR}(x) = e^{\pi R/2}\, K_{iR}(x)
      = e^{\pi R/2}\int_0^\infty e^{-x\cosh t}\cos(Rt)\,dt,

  which is the kernel of the main survey window (spectral parameter
  `s = 1 + iR`, `\lambda = 1 + R^2 \geq 1`).

- ``order='real'`` evaluates the real-order Bessel function

  .. MATH::

      K_t(x) = \int_0^\infty e^{-x\cosh u}\cosh(tu)\,du,

  the kernel of the *exceptional* interval (`s = 1 + t` real,
  `\lambda = 1 - t^2 < 1`).  As noted in §2.2 of the survey this is a two-line
  change from the imaginary-order kernel: ``cos(R t) -> cosh(t u)`` and drop
  the `e^{\pi R/2}` prefactor.  Everything downstream (pull-back, linear
  system, scaling, indicator, refinement, validation) is identical.

The nodes are chosen once from the largest argument via the truncation
analysis, so a single instance is reused across every spectral parameter in a
scan.  Accuracy is ``~1e-8`` for the imaginary order with 160 nodes over
`R \le 10`; the scaled result suffers intrinsic cancellation `~e^{-\pi R/2}`,
which is why ``float64`` (double precision) is required and sufficient for
Hejhal's method.
"""

import numpy as np


class KBessel:
    r"""
    Gauss--Legendre quadrature evaluator for `K_{iR}(x)` / `K_t(x)`.

    INPUT:

    - ``args`` -- iterable of nonnegative real arguments `x` at which the
      Bessel function is evaluated.  All arguments share a single set of
      quadrature nodes.
    - ``nodes`` -- number of Gauss--Legendre nodes (default: ``160``).
    - ``tmax`` -- upper truncation of the `t`-integral; if ``None`` it is
      chosen from the smallest positive argument so that the tail is below
      ``e^{-46}`` (default: ``None``).
    - ``chunk`` -- arguments are evaluated in blocks of this size to bound
      the temporary outer-product memory (default: ``200000``).
    - ``order`` -- ``'imaginary'`` for `\widetilde{K}_{iR}` (default) or
      ``'real'`` for `K_t`.

    EXAMPLES:

    The scaled imaginary-order Bessel function, validated against ``mpmath``
    (`\widetilde{K}_{iR}(x) = e^{\pi R/2} K_{iR}(x)`)::

        sage: from maass_forms_klein.functions.besselk_quad import KBessel
        sage: kb = KBessel([1.0, 5.0, 2.0, 0.5])
        sage: v = kb.eval(4.0)
        sage: bool(abs(v[0] - (-1.1570441018)) < 1e-6)      # R=4, x=1
        True
        sage: bool(abs(v[1] - 0.4296469267) < 1e-6)         # R=4, x=5
        True
        sage: bool(abs(kb.eval(6.5)[2] - 0.3051448917) < 1e-6)   # R=6.5, x=2
        True
        sage: bool(abs(kb.eval(2.0)[3] - 0.3818681483) < 1e-6)   # R=2, x=0.5
        True

    The real-order Bessel function `K_t(x)`, validated against
    ``scipy.special.kv``::

        sage: kbr = KBessel([1.0, 2.0, 0.5], order='real')
        sage: bool(abs(kbr.eval(0.5)[0] - 0.4610685044) < 1e-6)   # t=0.5, x=1
        True
        sage: bool(abs(kbr.eval(0.3)[1] - 0.1160369743) < 1e-6)   # t=0.3, x=2
        True

    The order can be switched at runtime, which is how the exceptional-interval
    scan reuses a context built for imaginary order::

        sage: kb.set_order('real')
        sage: bool(abs(kb.eval(0.5)[0] - 0.4610685044) < 1e-6)
        True
        sage: kb.set_order('imaginary')
        sage: bool(abs(kb.eval(4.0)[0] - (-1.1570441018)) < 1e-6)
        True

    An unknown order is rejected::

        sage: KBessel([1.0], order='complex')
        Traceback (most recent call last):
        ...
        ValueError: order must be 'imaginary' or 'real', got 'complex'
    """

    def __init__(self, args, nodes=160, tmax=None, chunk=200000, order="imaginary"):
        self.args = np.asarray(args, dtype=float)
        self.chunk = chunk
        pos = self.args[self.args > 1e-8]
        xmin = max(pos.min() if len(pos) else 1e-3, 1e-3)
        if tmax is None:
            tmax = float(np.arccosh(1.0 + 46.0 / xmin))
        t, w = np.polynomial.legendre.leggauss(nodes)
        self.t = 0.5 * tmax * (t + 1.0)
        self.w = 0.5 * tmax * w
        self.cosh_t = np.cosh(self.t)
        self.set_order(order)

    def set_order(self, order):
        r"""
        Set the Bessel order kind to ``'imaginary'`` or ``'real'``.

        EXAMPLES::

            sage: from maass_forms_klein.functions.besselk_quad import KBessel
            sage: kb = KBessel([1.0])
            sage: kb.order
            'imaginary'
            sage: kb.set_order('real'); kb.order
            'real'
        """
        if order not in ("imaginary", "real"):
            raise ValueError(f"order must be 'imaginary' or 'real', got {order!r}")
        self.order = order

    def eval(self, param):
        r"""
        Evaluate the Bessel function at every stored argument.

        INPUT:

        - ``param`` -- the spectral parameter: `R` for imaginary order
          (returns `\widetilde{K}_{iR}(x) = e^{\pi R/2} K_{iR}(x)`), or `t`
          for real order (returns `K_t(x)`).

        OUTPUT: a ``numpy`` ``float64`` array aligned with ``self.args``.

        EXAMPLES::

            sage: from maass_forms_klein.functions.besselk_quad import KBessel
            sage: kb = KBessel([1.0, 5.0])
            sage: len(kb.eval(4.0))
            2
        """
        if self.order == "imaginary":
            cw = np.cos(param * self.t) * self.w
        else:
            cw = np.cosh(param * self.t) * self.w
        out = np.empty(len(self.args))
        for i in range(0, len(self.args), self.chunk):
            blk = self.args[i : i + self.chunk]
            out[i : i + self.chunk] = np.exp(-np.outer(blk, self.cosh_t)) @ cw
        if self.order == "imaginary":
            return np.exp(0.5 * np.pi * param) * out
        return out
