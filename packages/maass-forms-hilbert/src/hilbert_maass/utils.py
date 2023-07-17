from hilbert_modgroup.pullback import HilbertPullback
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.rings.infinity import Infinity


def get_Q_from_bounds(P: HilbertPullback, b: tuple) -> tuple:
    """
    Find the bounds for the box [-C,C]^n for the integrer coordinates correponding to the box [-b,b]^n in the lattice.
    """
    C = 1
    for ida in P.group().ideal_cusp_representatives():
        t = matrix(P.basis_matrix_ideal(ida)).transpose().norm(Infinity)
        if t > C:
            C = t
    C = C*max(b)
    return (C,) * len(b)


@cached_function
def map_tuple_to_int(index_tuple: tuple, min_tix: tuple, max_tix: int) -> int:
    """
    Map a tuple (a,b) with min_tix < a, b < max_tix to
    an integer (max_tix - min_tx +1)*(a - min_tix) + b - min_tx
    and generalise this to longer tuples.

    EXAMPLES::

        sage: from

    """
    if not isinstance(index_tuple, tuple):
        raise ValueError("Call with tuple!")
    n = len(index_tuple)
    range_tix = (max_tix - min_tix + 1)
    # Check if any tuple elements are out of bounds.
    if max(index_tuple) - min_tix > range_tix or min(index_tuple) - min_tix < 0:
        raise ValueError(f"The out of bounds value(s) in {index_tuple}!")
    return sum(range_tix**(n-i-1)*(index_tuple[i] - min_tix) for i in range(n))


@cached_function
def map_int_to_tuple(index: int, min_tix: int, max_tix: int, len_tuple: int) -> tuple:
    r"""
    Map integer to tuple (the inverse of map_tuple_to_int)
    :param index:
    :param min_tix:
    :param max_tix:
    :param len_tuple:
    :return:
    """
    range_tix = (max_tix - min_tix + 1)
    ix_t = []
    for i in range(len_tuple-1, -1, -1):
        t = index % range_tix
        ix_t.append(t + min_tix)
        index = (index - t)/range_tix
    ix_t.reverse()
    return tuple(ix_t)
