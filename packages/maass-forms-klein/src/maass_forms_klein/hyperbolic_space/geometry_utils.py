from copy import deepcopy
from typing import ParamSpec, Union

from sage.functions.other import imag, real
from sage.misc.cachefunc import cached_function
from sage.misc.functional import sqrt
from sage.modules.free_module_element import vector
from sage.plot.circle import circle
from sage.plot.polygon import polygon
from sage.rings.complex_mpfr import ComplexField
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.real_mpfr import RR, RealNumber
from sage.structure.element import Matrix, Vector

from maass_forms_klein.hyperbolic_space.parallelogram import (
    parallelogram_covered_by_circles,
    parallelogram_in_circle,
    parallelogram_intersect_circle,
    split_parallelogram,
)
from maass_forms_klein.hyperbolic_space.types import Circle, Parallelogram, Rectangle
from maass_forms_klein.hyperbolic_space.word_utils import word_list_sort_key, word_to_circle

# Define types locally to avoid import chain issues
Real_t = Union[RealNumber, Integer, Rational, int, float]
Integer_t = Union[Integer, int]


# Define get_epsilon locally to avoid import issues
def get_epsilon(element: Real_t):
    """
    Get machine epsilon for given precision.

    INPUT:

    - ``element`` -- a numerical element whose precision determines epsilon

    OUTPUT:

    A small positive number representing machine epsilon x 16 for the given precision.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import get_epsilon
        sage: get_epsilon(1.0)
        1.77...e-15
        sage: from sage.rings.real_mpfr import RealField
        sage: RF = RealField(100)
        sage: get_epsilon(RF(1.0))
        1.26...e-29

    """
    if not hasattr(element, "base_ring"):
        dprec = 53
    elif element.base_ring().epsilon() == 0:
        return 0
    else:
        dprec = element.base_ring().prec()
    return 2 ** (4 - dprec)


P = ParamSpec("P")


def reduce_cover(
    rect: Rectangle | Parallelogram,
    cover_list: list[str | tuple[str, Circle]],
    gens: dict | None = None,
    max_n: Integer_t = 10,
    fix_circles: list[str] | None = None,
    verbose: int = 0,
    return_words: bool = False,
) -> list[tuple[str, Circle] | str]:
    r"""
    Given a list of group elements that cover a rectangle find a smaller covering subset.

    INPUT:

    - ``rect`` -- a Rectangle or Parallelogram to be covered
    - ``cover_list`` -- list of words (strings) or tuples ``(word, Circle)``
    - ``gens`` -- (optional) dictionary of generators, required if cover_list contains strings
    - ``max_n`` -- (default: 10) maximum subdivision level for coverage checks
    - ``fix_circles`` -- (optional) list of words that should not be removed
    - ``verbose`` -- (default: 0) verbosity level
    - ``return_words`` -- (default: False) whether to return words instead of circles
    OUTPUT:

    A reduced list of tuples ``(word, Circle)`` that still covers the rectangle.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import reduce_cover
        sage: from maass_forms_klein.hyperbolic_space.types import Circle, Parallelogram
        sage: p = Parallelogram(base=vector((0,0)), v1=vector((0.1,0)), v2=vector((0,0.1)))
        sage: c1 = ('a', Circle(center=vector((0.05, 0.05)), radius=0.2))
        sage: reduce_cover(p, [c1])
        [('a', Circle(center=(0.0500000000000000, 0.0500000000000000), radius=0.200000000000000))]

    """
    # First remove circles that are duplicates of other circles
    cover_list = remove_duplicate_circles(cover_list, gens)
    cover_list.sort(key=lambda x: x[1].radius, reverse=True)
    # Then remove circles that are contained in other circles
    cover_list = remove_contained_circles(cover_list, gens=gens)
    if not fix_circles:
        fix_circles = []
    new_list = deepcopy(cover_list)
    # Check that the original list is covering the parallelogram
    if not parallelogram_covered_by_circles(
        rect, [w[1] for w in new_list], scaling_factor=0.99, n_max=max_n
    ):
        raise ArithmeticError("The original list is not covering the parallelogram.")

    # Now we test and see which circles we can remove and still cover the parallelogram
    # Starting with smaller circles
    cover_list.sort(key=lambda x: x[1].radius)
    for x, c in cover_list:
        if len(new_list) == 1:
            break
        if x in fix_circles:
            continue
        new_list.remove((x, c))
        try:
            test = parallelogram_covered_by_circles(
                rect, [w[1] for w in new_list], scaling_factor=0.99, n_max=max_n
            )
        except ArithmeticError:
            test = False
        if test:
            continue
        else:
            new_list.append((x, c))
    new_list.sort(key=lambda x: word_list_sort_key(x[0]))
    if return_words:
        return [w[0] for w in new_list]
    return new_list


def _ensure_tuples(
    cover_list: list[str | tuple[str, Circle]], gens: dict | None = None
) -> list[tuple[str, Circle]]:
    r"""
    Ensure that the cover list is a list of tuples ``(word, Circle)``.

    If ``cover_list`` contains plain strings, convert them to tuples using ``gens``
    to compute the corresponding invariant circles.

    INPUT:

    - ``cover_list`` -- list of strings or tuples ``(word, Circle)``
    - ``gens`` -- (optional) dictionary of generators, required if cover_list contains strings

    OUTPUT:

    A list of tuples ``(word, Circle)``.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import _ensure_tuples
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c = Circle(center=vector((0, 0)), radius=1)
        sage: _ensure_tuples([('a', c)])
        [('a', Circle(center=(0, 0), radius=1))]

    """
    if not all(isinstance(x, tuple) for x in cover_list):
        if not gens:
            raise ValueError("If cover_list is not a list of tuples, gens must be provided.")
        cover_list = [(w, word_to_circle(w, gens)) for w in cover_list]
    return cover_list


def remove_duplicate_circles(
    cover_list: list[tuple[str, Circle] | str],
    gens: dict | None = None,
    include_parabolic: bool = False,
) -> list[tuple[str, Circle]]:
    """
    Given a list of words remove duplicates corresponding to the same
    invariant circles.

    INPUT:

    - ``cover_list`` -- (list) list of words
    - ``gens`` -- (dictionary) generators of the group
    - ``include_parabolic`` -- (bool) if True include parabolic circles in the list

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import remove_duplicate_circles
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: gens = {'a': matrix([[0, 1], [1, 0]]), 'b': matrix([[0, 1], [1, 0]])}
        sage: c1 = 'a', Circle(center=vector((1, 1)), radius=1)
        sage: c2 = 'a', Circle(center=vector((1, 1)), radius=1)
        sage: c3 = 'a', Circle(center=vector((2, 1)), radius=1)
        sage: c4 = 'a', Circle(center=vector((2, 1)), radius=1.0000000001)  # Nearly identical to c3
        sage: c5 = 'a', Circle(center=vector((3, 1)), radius=0.5)
        sage: c6 = 'a', Circle(center=vector((3, 1)), radius=0.5)
        sage: l = [c1, c2, c3, c4, c5, c6]
        sage: result = remove_duplicate_circles(l)
        sage: len(result)
        3
        sage: [x[1].center for x in result]
        [(1, 1), (2, 1), (3, 1)]

    Test with empty list::

        sage: remove_duplicate_circles([])
        []
    """
    cover_list = _ensure_tuples(cover_list, gens)
    if not cover_list:
        return []
    # Get 'infinity' cut-off
    elt = cover_list[0][1].radius
    if hasattr(elt, "base_ring"):
        eps = elt.base_ring().epsilon()
        if eps == 0:  # Handle rings like Integer Ring with epsilon = 0
            eps = 1e-10
    else:
        eps = 1e-10
    max_word_len = max([len(c[0]) for c in cover_list]) + 2
    eps = eps * 2**max_word_len
    infinity_value = 1 / eps

    cover_list.sort(key=lambda x: (x[1].center[0], x[1].center[1], x[1].radius))
    new_list = []
    n = 0
    while n < len(cover_list):
        if not include_parabolic and cover_list[n][1].radius >= infinity_value:
            n += 1
            continue
        cn = cover_list[n]
        new_list.append(cn)
        while n < len(cover_list) and cn[1].within_epsilon(cover_list[n][1], eps):
            n += 1
    return new_list


def remove_contained_circles(
    circle_list: list[tuple[str, Circle] | str],
    gens: dict | None = None,
    include_parabolic: bool = False,
) -> list[tuple[str, Circle]]:
    """
    Given a list of words remove duplicates corresponding to the same
    invariant circles.

    INPUT:

    - ``circle_list`` -- (list) list of words
    - ``gens`` -- (dictionary) generators of the group
    - ``include_parabolic`` -- (bool) if True include parabolic circles in the list

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import remove_contained_circles
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c1 = ('a', Circle(center=vector((0, 0)), radius=1))
        sage: c2 = ('b', Circle(center=vector((0, 0)), radius=2))
        sage: c3 = ('c', Circle(center=vector((5, 5)), radius=0.5))
        sage: result = remove_contained_circles([c1, c2, c3])
        sage: len(result)
        2

    """
    circle_list = _ensure_tuples(circle_list, gens)
    if not circle_list:
        return []
    radii = [c[1].radius for c in circle_list]
    eps = max(get_epsilon(r) for r in radii)
    if eps == 0:
        eps = 2**-53
    max_word_len = max([len(c[0]) for c in circle_list]) + 2
    eps = eps * 2**max_word_len
    infinity_value = 1 / eps
    circle_list.sort(key=lambda x: x[1].radius, reverse=True)
    new_list = []
    for w, c in circle_list:
        if not include_parabolic and c.radius >= infinity_value:
            continue
        if not circle_in_circle_in_list(c, [x[1] for x in new_list]):
            new_list.append((w, c))
    return new_list


def remove_non_intersecting_circles(
    circle_list: list[tuple[str, Circle]],
    para: Parallelogram,
    gens: dict | None = None,
) -> list[tuple[str, Circle]]:
    """
    Given a list of words remove duplicates corresponding to the same
    invariant circles.

    INPUT:

    - ``circle_list`` -- (list) list of words
    - ``para`` -- (Parallelogram)
    - ``gens`` -- (dictionary) generators of the group

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import (
        ....:     remove_non_intersecting_circles)

        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: gens = {'a': matrix([[0, 1], [1, 0]]), 'b': matrix([[0, 1], [1, 0]])}
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
        sage: para = Parallelogram(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: circle_list = [(0, Circle(center=vector((1, 1)), radius=1)),
        ....:                (1, Circle(center=vector((2, 1)), radius=1)),
        ....:                (2, Circle(center=vector((3, 1)), radius=1))]
        sage: result = remove_non_intersecting_circles(circle_list, para, gens)
        sage: len(result)
        2

    """
    circle_list = _ensure_tuples(circle_list, gens)
    circle_list = [x for x in circle_list if x[1].radius != Infinity]
    circle_list.sort(key=lambda x: x[1].radius, reverse=True)
    new_list = []
    for w, c in circle_list:
        if parallelogram_intersect_circle(para, c):
            new_list.append((w, c))
    return new_list


def circle_is_contained(c1: Circle, c2: Circle) -> bool:
    """
    Return True if the circle with centre c1 and radius r1 is contained in the circle
    with centre c2 and radius r3

    EXAMPLES:
        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import circle_is_contained
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c1 = Circle(center=vector((0, 0)), radius=1)
        sage: c2 = Circle(center=vector((0, 0)), radius=2)
        sage: circle_is_contained(c1, c1)
        True
        sage: circle_is_contained(c1, c2)
        True
        sage: c3 = Circle(center=vector((0, 0)), radius=0.5)
        sage: circle_is_contained(c1, c3)
        False
        sage: c4 = Circle(center=vector((0, 0)), radius=1.5)
        sage: c5 = Circle(center=vector((0, 0.5)), radius=1)
        sage: circle_is_contained(c5, c4)
        True
    """
    if not isinstance(c1, Circle) or not isinstance(c2, Circle):
        raise ValueError("Input should be circles")
    vx = c1.center[0] - c2.center[0]
    vy = c1.center[1] - c2.center[1]
    if vx == 0 and vy == 0:
        return c1.radius <= c2.radius
    absv = sqrt(vx**2 + vy**2)
    vnx_c1radius = vx / absv * c1.radius
    vny_c1radius = vy / absv * c1.radius
    abs_plus = sqrt((vx + vnx_c1radius) ** 2 + (vy + vny_c1radius) ** 2)
    abs_minus = sqrt((vx - vnx_c1radius) ** 2 + (vy - vny_c1radius) ** 2)
    return abs_plus <= c2.radius and abs_minus <= c2.radius


def circle_intersects_circle(c1: Circle, c2: Circle) -> bool:
    """
    Check if circle c1 intersects circle c2.

    INPUT:

    - ``c1`` -- a Circle
    - ``c2`` -- a Circle

    OUTPUT:

    ``True`` if the two circles intersect (including touching), ``False`` otherwise.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import circle_intersects_circle
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c1 = Circle(center=vector((0, 0)), radius=1)
        sage: c2 = Circle(center=vector((1, 0)), radius=1)
        sage: circle_intersects_circle(c1, c2)
        True
        sage: c3 = Circle(center=vector((3, 0)), radius=1)
        sage: circle_intersects_circle(c1, c3)
        False

    """
    vx = c1.center[0] - c2.center[0]
    vy = c1.center[1] - c2.center[1]
    return sqrt(vx**2 + vy**2) <= c1.radius + c2.radius


def circle_approx_equal(c1: Circle, c2: Circle, prec: Real_t = 1e-10) -> bool:
    """
    Return True if the circle c1 is approximately equal to the circle c2

    EXAMPLES:
        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import circle_approx_equal
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c1 = Circle(center=vector((0, 0)), radius=1)
        sage: c2 = Circle(center=vector((0, 0)), radius=2)
        sage: circle_approx_equal(c1, c1)
        True
        sage: circle_approx_equal(c1, c2)
        False
        sage: c3 = Circle(center=vector((0, 0)), radius=0.5)
        sage: circle_approx_equal(c1, c3)
        False
        sage: c4 = Circle(center=vector((0, 0)), radius=1.5)
        sage: c5 = Circle(center=vector((0, 0.5)), radius=1)
        sage: circle_approx_equal(c5, c4)
        False
        sage: c6 = Circle(center=vector((0, 0)), radius=1.5+1e-15)
        sage: circle_approx_equal(c6, c4)
        True
    """
    if not isinstance(c1, Circle) or not isinstance(c2, Circle):
        raise ValueError("Input should be circles")
    v = vector(c1.center) - vector(c2.center)
    return abs(v) <= prec and abs(c1.radius - c2.radius) <= prec


def invariant_sphere_is_contained(A1, A2):
    r"""
    Check if the invariant sphere of matrix ``A1`` is contained in the invariant
    sphere of matrix ``A2``.

    INPUT:

    - ``A1`` -- a 2x2 matrix in PSL(2,C)
    - ``A2`` -- a 2x2 matrix in PSL(2,C)

    OUTPUT:

    ``True`` if the invariant circle of ``A1`` is contained in that of ``A2``.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import (
        ....:     invariant_sphere_is_contained)
        sage: A1 = matrix(CC, [[0, -1], [1, 0]])
        sage: A2 = matrix(CC, [[0, -2], [0.5, 0]])
        sage: invariant_sphere_is_contained(A1, A1)
        True
        sage: invariant_sphere_is_contained(A1, A2)
        True

    """
    c = -A1[1, 1] / A1[1, 0]
    r1 = 1 / abs(A1[1, 0])
    c1 = Circle(center=vector((real(c), imag(c))), radius=r1)
    c = -A2[1, 1] / A2[1, 0]
    r2 = 1 / abs(A2[1, 0])
    c2 = Circle(center=vector((real(c), imag(c))), radius=r2)
    return circle_is_contained(c1, c2)


def rectangle_in_circle(rect: Rectangle, circle: Circle) -> bool:
    """
    Check if rectangle given by center and side lengths.
    is contained in circle given by center and radius.

    EXAMPLE:
        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import rectangle_in_circle
        sage: from maass_forms_klein.hyperbolic_space.types import Rectangle, Circle
        sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))

        sage: rectangle_in_circle(rect, Circle(center=vector((0, 0)), radius=1))
        False
        sage: rectangle_in_circle(rect, Circle(center=vector((0, 0)), radius=2))
        True
    """
    for sgn_x in [-1, 1]:
        corner_x = rect.center[0] + sgn_x * rect.sides[0] / 2
        for sgn_y in [-1, 1]:
            corner_y = rect.center[1] + sgn_y * rect.sides[1] / 2
            if abs(vector((corner_x, corner_y)) - circle.center) > circle.radius:
                return False
    return True


def is_point_in_rectangle(rect, point):
    """
    Check if point is inside a rectangle.

    INPUT:

    - ``rect`` -- a Rectangle
    - ``point`` -- a 2-dimensional vector or tuple

    OUTPUT:

    ``True`` if the point lies inside the rectangle (within machine epsilon).

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import is_point_in_rectangle
        sage: from maass_forms_klein.hyperbolic_space.types import Rectangle
        sage: rect = Rectangle(base=vector(RR, (0, 0)),
        ....:     v1=vector(RR, (1, 0)), v2=vector(RR, (0, 1)))
        sage: is_point_in_rectangle(rect, vector(RR, (0.3, 0.3)))
        True
        sage: is_point_in_rectangle(rect, vector(RR, (2.0, 0.3)))
        False

    """
    x, y = point
    cRx, cRy = rect.center
    lx, ly = rect.sides
    # allow for ploss of precision in matrix elements and rectangle
    eps = (2**4) * x.base_ring().epsilon()
    rectangle_left = cRx - lx / 2 - eps
    rectangle_right = cRx + lx / 2 + eps
    rectangle_bottom = cRy - ly / 2 - eps
    rectangle_top = cRy + ly / 2 + eps
    return rectangle_left <= x <= rectangle_right and rectangle_bottom <= y <= rectangle_top


def split_rectangle(rectangle: Rectangle, n1: Integer_t = 2, n2: Integer_t = 2) -> list[Rectangle]:
    """
    Split a rectangle into a set of smaller rectangles of the same size.

    INPUT:

    - ``rectangle_center`` -- Center of rectangle

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import split_rectangle
        sage: from maass_forms_klein.hyperbolic_space.types import Rectangle
        sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))

        sage: rects = list(split_rectangle(rect, 2, 2))
        sage: len(rects)
        4
        sage: rects[0].base
        (0, 0)
        sage: rects[0].v1
        (1/2, 0)
        sage: rects[0].v2
        (0, 1/2)
    """
    rectangles = split_parallelogram(rectangle, n1, n2)
    return [Rectangle(base=rect.base, v1=rect.v1, v2=rect.v2) for rect in rectangles]


@cached_function
def matrix_to_circle(mat: Matrix) -> Circle:
    r"""
    Find the invariant circle of a matrix.

    Note: If lower-left entry is 0 the invariant circle is a vertical line parallel
    with the imaginary axis.

    INPUT:

    - ``mat`` -- a 2x2 matrix

    OUTPUT:

    A Circle representing the invariant circle of the matrix.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import matrix_to_circle
        sage: A = matrix(CC, [[0, -1], [1, 0]])
        sage: c = matrix_to_circle(A)
        sage: c.radius
        1.00000000000000
        sage: c.center # abstol 1e-100
        (0.000000000000000, 0.000000000000000)
        sage: matrix_to_circle(identity_matrix(CC, 2)).radius
        +Infinity

    """
    if not isinstance(mat, Matrix):
        raise ValueError(f"Need matrix input. Got: {mat}.")
    if mat[1, 0] == 0:
        return Circle(center=vector((0, 0)), radius=Infinity)
    if hasattr(mat.base_ring(), "prec"):
        prec = mat.base_ring().prec()
    else:
        prec = 53
    if hasattr(mat[1, 0], "n"):
        c, d = mat[1, 0].n(prec), mat[1, 1].n(prec)
    else:
        CF = ComplexField(prec)
        c, d = CF(mat[1, 0]), CF(mat[1, 1])
    radius = 1 / abs(c)
    center = -d / c
    return Circle(center=vector((real(center), imag(center))), radius=radius)


def rectangle_covered_by_matrices(
    rect: Rectangle, matrices, scaling_factor=1.0, n_max=10000
) -> bool:
    """
    Return True if rectangle is covered by the invariant circles given by the list of matrices
    scaled with `scaling_factor`.

    INPUT:

    - ``rect`` -- (tuple of real numbers) center of the rectangle
    - ``sides`` -- (tuple of real numbers) sides of the rectangle
    - ``matrices`` -- (list of matrices) list of matrices

    EXAMPLES:
        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import (
        ....:     rectangle_covered_by_matrices)

        sage: from maass_forms_klein.hyperbolic_space.types import Rectangle
        sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: rectangle_covered_by_matrices(rect, [identity_matrix(2)])
        Traceback (most recent call last):
        ...
        ValueError: Matrix has infinite invariant circle

    """

    if not isinstance(rect, Rectangle):
        raise ValueError("rect must be given as a Rectangle")

    if not isinstance(matrices, list):
        matrices = [matrices]
    if not all(isinstance(x, Matrix) for x in matrices):
        raise ValueError(f"Input should be a list of matrices. Got: {matrices}")
    max_radius = max([abs(matrix_to_circle(g).radius) for g in matrices])
    if max_radius == Infinity:
        raise ValueError("Matrix has infinite invariant circle")
    max_side_length = max_radius / RR(2.0).sqrt()
    n = 1
    while n < n_max:
        all_are_covered = True
        for rect_small in split_rectangle(rect, n):
            if not one_rectangle_covered_by_matrices(
                rect_small, matrices, scaling_factor=scaling_factor
            ):
                all_are_covered = False
                break
        if all_are_covered:
            return True
        if rect.v1.norm() / n < max_side_length and rect.v2.norm() / n < max_side_length:
            return False
        n += 1
    raise ArithmeticError("Could not determine if rectangle is covered in given number of steps.")


def one_rectangle_covered_by_matrices(
    rect: Rectangle, matrices: list[Matrix], scaling_factor: Real_t = 1.0
):
    """
    Return True if rectangle is covered by the invariant circles given by the list of matrices
    scaled by the scaling factor.

    INPUT:

    - ``center`` -- (tuple of real numbers) center of the rectangle
    - ``sides`` -- (tuple of real numbers) sides of the rectangle
    - ``matrices`` -- (list of matrices) list of matrices

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import (
        ....:     one_rectangle_covered_by_matrices)
        sage: from maass_forms_klein.hyperbolic_space.types import Rectangle
        sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1/2, 0)), v2=vector((0, 1/2)))
        sage: one_rectangle_covered_by_matrices(rect, [matrix([[0, -1],[1, 0]])])
        True
        sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: one_rectangle_covered_by_matrices(rect, [matrix([[0, -1],[1, 0]])])
        False

    """
    if not isinstance(matrices, list):
        matrices = [matrices]
    for g in matrices:
        cr = matrix_to_circle(g)
        cr = cr.scale(scaling_factor)
        if rectangle_in_circle(rect, cr):
            return True
    return False


def display_rectangle_and_matrices(rect, matrices, **kwargs):
    r"""
    Create a plot of a rectangle and the invariant circles of the given matrices.

    INPUT:

    - ``rect`` -- a Rectangle
    - ``matrices`` -- list of 2x2 matrices
    - ``**kwargs`` -- optional keyword arguments passed to the circle plot
      (``alpha``, ``thickness``, ``color``)

    OUTPUT:

    A SageMath graphics object.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import (
        ....:     display_rectangle_and_matrices)
        sage: from maass_forms_klein.hyperbolic_space.types import Rectangle
        sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: p = display_rectangle_and_matrices(rect, [matrix(CC, [[0, -1],[1, 0]])])

    """
    center = rect.center
    sides = rect.sides
    plot = polygon(
        [
            [center[0] - sides[0] / 2, center[1] - sides[1] / 2],
            [center[0] - sides[0] / 2, center[1] + sides[1] / 2],
            [center[0] + sides[0] / 2, center[1] + sides[1] / 2],
            [center[0] + sides[0] / 2, center[1] - sides[1] / 2],
        ],
        alpha=0.5,
        thickness=0.2,
    )
    for g in matrices:
        cr = matrix_to_circle(g)
        plot += circle(
            tuple(cr.center),
            cr.radius,
            alpha=kwargs.get("alpha", 0.5),
            thickness=kwargs.get("thickness", 0.2),
            color=kwargs.get("color", "red"),
        )
    return plot


def circle_in_list(new_matrix_or_circle, existing_matrices_or_circles):
    """
    Find if a matrix has an invariant circle within a given list up to
    machine precision.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import circle_in_list
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c1 = Circle(center=vector((0, 0)), radius=1)
        sage: circle_in_list(c1, [c1])
        True
        sage: c2 = Circle(center=vector((0, 0)), radius=1.0)
        sage: c3 = Circle(center=vector((1, 0)), radius=1)
        sage: circle_in_list(c2, [c2, c3])
        True
        sage: circle_in_list(c2, [c3])
        False

    """
    if isinstance(new_matrix_or_circle, Matrix):
        new_matrix_or_circle = matrix_to_circle(new_matrix_or_circle)
    if isinstance(new_matrix_or_circle.radius, float):
        eps = 2 ** (4 - 53)
    else:
        eps = 2**4 * new_matrix_or_circle.radius.base_ring().epsilon()
    if not isinstance(existing_matrices_or_circles[0], Circle):
        convert = True
    else:
        convert = False
    for c in existing_matrices_or_circles:
        if convert:
            c = matrix_to_circle(c)
        if c.within_epsilon(new_matrix_or_circle, eps):
            return True
    return False


def circle_in_circle_in_list(new_matrix_or_circle, existing_matrices_or_circles):
    """
    Check if a circle (or invariant circle of a matrix) is contained in any
    circle in the given list.

    INPUT:

    - ``new_matrix_or_circle`` -- a Circle or a 2x2 matrix
    - ``existing_matrices_or_circles`` -- list of Circles or matrices

    OUTPUT:

    ``True`` if the circle is contained in one of the circles in the list.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import circle_in_circle_in_list
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c1 = Circle(center=vector((0, 0)), radius=0.5)
        sage: c2 = Circle(center=vector((0, 0)), radius=2)
        sage: circle_in_circle_in_list(c1, [c2])
        True
        sage: circle_in_circle_in_list(c2, [c1])
        False

    """
    if isinstance(new_matrix_or_circle, Matrix):
        new_matrix_or_circle = matrix_to_circle(new_matrix_or_circle)
    for c in existing_matrices_or_circles:
        if isinstance(c, Matrix):
            c = matrix_to_circle(c)
        if circle_is_contained(new_matrix_or_circle, c):
            return True
    return False


def rectangle_intersect_circle(rect: Rectangle, circ: Circle) -> bool:
    """
    Check if a rectangle (parallel with the coordinate axes) and a circle intersect.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import (
        ....:     rectangle_intersect_circle)
        sage: from maass_forms_klein.hyperbolic_space.types import Rectangle, Circle
        sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: rectangle_intersect_circle(rect, Circle(center=vector((0, 0)), radius=1))
        True
        sage: rectangle_intersect_circle(rect, Circle(center=vector((0, 0)), radius=2))
        True
        sage: rectangle_intersect_circle(rect, Circle(center=vector((0, 0)), radius=0.1))
        True
        sage: rectangle_intersect_circle(rect, Circle(center=vector((0, 0)), radius=5))
        True

    """
    cx = circ.center[0]
    cy = circ.center[1]
    r = circ.radius
    cRx, cRy = rect.center
    lx, ly = rect.sides
    rectangle_left = cRx - lx / 2
    rectangle_right = cRx + lx / 2
    rectangle_bottom = cRy - ly / 2
    rectangle_top = cRy + ly / 2
    closest_x = max(rectangle_left, min(cx, rectangle_right))
    closest_y = max(rectangle_bottom, min(cy, rectangle_top))
    return (cx - closest_x) ** 2 + (cy - closest_y) ** 2 <= r**2


def is_circle_covered_by_circles(
    c: Circle, circles: list[Circle], scaling_factor: Real_t = 1.0, n_max: Integer_t = 100
) -> bool:
    r"""
    Return True if circle ``c`` is covered by the circles in the list.

    INPUT:

    - ``c`` -- a Circle
    - ``circles`` -- list of Circles
    - ``scaling_factor`` -- (default: 1.0) scaling factor for the covering circles
    - ``n_max`` -- (default: 100) maximum subdivision level

    OUTPUT:

    ``True`` if ``c`` is covered by the union of the circles in the list.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import (
        ....:     is_circle_covered_by_circles)
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c = Circle(center=vector((0, 0)), radius=0.5)
        sage: c2 = Circle(center=vector((0, 0)), radius=1)
        sage: is_circle_covered_by_circles(c, [c2])
        True

    """
    # First check if c is inside one of them
    for cother in circles:
        if circle_is_contained(c, cother):
            return True
    # Otherwise we do a subdivision algorithm over a parallelogram cover of the circle.
    para = Parallelogram(
        base=c.center - vector((c.radius, c.radius)),
        v1=vector((0, 2 * c.radius)),
        v2=vector((2 * c.radius, 0)),
    )
    # Ignore circles that do not intersect the rectangle
    circles = [cother for cother in circles if parallelogram_intersect_circle(para, cother)]
    # Then we look at smaller parallelograms until we either
    # 1. find a small parallelogram that is in c but none of the other circles, or
    # 2. find that all small parallelograms that intersect c also are
    #    inside one of the other circles.
    n = 1
    while n < n_max:
        all_are_covered = True
        split = split_parallelogram(para, n1=n, n2=n)
        for para_small in split:
            if not parallelogram_intersect_circle(para_small, c):
                continue
            # If the small parallelogram does not intersect any of the circles
            # then it is clearly not covered
            if not any(parallelogram_intersect_circle(para_small, c) for c in circles):
                return False
            # If it is not covered by any circle then we need to divide further
            if not any(
                parallelogram_in_circle(c, para_small, scaling_factor=scaling_factor)
                for c in circles
            ):
                all_are_covered = False
                break
        if all_are_covered:
            return True
        n += 1
    raise ArithmeticError("Could not determine if circle is covered in given number of steps.")


def point_in_circle(p: Vector, c: Circle, scaling_factor: Real_t = 1.0) -> bool:
    """
    Return True if the point ``p`` is inside the circle ``c`` (optionally scaled).

    INPUT:

    - ``p`` -- a 2-dimensional vector
    - ``c`` -- a Circle
    - ``scaling_factor`` -- (default: 1.0) scaling factor applied to the radius

    OUTPUT:

    ``True`` if the point lies inside the (scaled) circle.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import point_in_circle
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c = Circle(center=vector((0, 0)), radius=1)
        sage: point_in_circle(vector((0.5, 0)), c)
        True
        sage: point_in_circle(vector((2, 0)), c)
        False

    """
    return sqrt((p[0] - c.center[0]) ** 2 + (p[1] - c.center[1]) ** 2) <= c.radius * scaling_factor


def point_covered_by_circles(
    p: Vector, circles: list[Circle], scaling_factor: Real_t = 1.0
) -> bool:
    r"""
    Return True if the point ``p`` is inside at least one of the given circles.

    INPUT:

    - ``p`` -- a 2-dimensional vector
    - ``circles`` -- list of Circles
    - ``scaling_factor`` -- (default: 1.0) scaling factor applied to each circle radius

    OUTPUT:

    ``True`` if the point lies inside at least one circle.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import point_covered_by_circles
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: c1 = Circle(center=vector((0, 0)), radius=1)
        sage: c2 = Circle(center=vector((3, 0)), radius=1)
        sage: point_covered_by_circles(vector((0.5, 0)), [c1, c2])
        True
        sage: point_covered_by_circles(vector((1.5, 0)), [c1, c2])
        False

    """
    return any(point_in_circle(p, c, scaling_factor) for c in circles)


def display_parallelogram_and_circles(
    par: Parallelogram, circles: list[Circle | Matrix | str] | None = None, **kwargs: P.kwargs
):
    r"""
    Create a plot of a parallelogram and a list of circles (or invariant circles of matrices).

    INPUT:

    - ``par`` -- a Parallelogram
    - ``circles`` -- (optional) list of Circles, matrices, or words
    - ``**kwargs`` -- optional keyword arguments: ``alpha``, ``thickness``, ``color``,
      ``pcolor``, ``gens``

    OUTPUT:

    A SageMath graphics object.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.geometry_utils import (
        ....:     display_parallelogram_and_circles)
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram, Circle
        sage: par = Parallelogram(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: c = Circle(center=vector((0.5, 0.5)), radius=0.3)
        sage: p = display_parallelogram_and_circles(par, [c])

    """
    alpha = kwargs.get("alpha", 0.5)
    thickness = kwargs.get("thickness", 0.2)
    color = kwargs.get("color", "red")
    pcolor = kwargs.get("pcolor", "blue")
    plot = polygon(par.vertices, alpha=alpha, thickness=thickness, color=pcolor)
    if not circles:
        circles = []
    for cr in circles:
        if isinstance(cr, Matrix):
            cr = matrix_to_circle(cr)
        if isinstance(cr, str):
            cr = word_to_circle(cr, kwargs.get("gens", []))
        if cr.radius == Infinity:
            continue
        plot += circle(cr.center, cr.radius, alpha=alpha, thickness=thickness, color=color)
    return plot
