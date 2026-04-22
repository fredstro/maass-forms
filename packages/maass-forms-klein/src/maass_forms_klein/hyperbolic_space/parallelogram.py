"""
Functions for working with translation domains that are parallelograms.
We represent a parallelogram as a triple of vectors (v1,v2,b)
where b is the base and v1 and v2 span the parallelogram sides.
"""
from maass_forms_klein.hyperbolic_space.word_utils import translation_tuple_to_word
from sage.all import RR
from sage.categories.sets_cat import cartesian_product
from sage.functions.other import floor
from sage.functions.trig import cos, sin
from sage.misc.functional import sqrt
from sage.modules.free_module_element import vector
from sage.rings.infinity import Infinity
from sage.structure.element import Vector, Matrix
# Define types locally to avoid import chain issues
from sage.rings.real_mpfr import RealNumber
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from typing import Union
Real_t = Union[RealNumber, Integer, Rational, int, float]
Integer_t = Union[Integer, int]

from maass_forms_klein.hyperbolic_space.types import Parallelogram, Circle, Line, Rectangle



def reduce_in_parallelogram(p: Parallelogram, points: list[Vector],
                            return_translations: bool = False) -> tuple | list[Vector]:
    """
    Reduce a list of points by vectors spanning the parallelogram into a
     fundamental domain corresponding to the given base of the parallelogram.

    INPUT:
    - ``p`` - Parallelogram
    - ``points`` - list of points

    EXAMPLES:

        sage: from maass_forms_klein.hyperbolic_space.parallelogram import reduce_in_parallelogram
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
        sage: from maass_forms_klein.hyperbolic_space.types import Point
        sage: p = Parallelogram(base=(0,0), v1=(0,1), v2=(1,0))
        sage: pts = [Point(0,0), Point(0,1), Point(1,0), Point(1,1)]
        sage: reduce_in_parallelogram(p, pts)
        [(0, 0), (0, 0), (0, 0), (0, 0)]
        sage: pts = [Point(0,0), Point(0,-1), Point(-1,0), Point(-1,-1)]
        sage: reduce_in_parallelogram(p, pts)
        [(0, 0), (0, 0), (0, 0), (0, 0)]
        sage: pts = [Point(2,0), Point(0,2), Point(0.5,0.5), Point(2,2)]
        sage: reduce_in_parallelogram(p, pts)
        [(0, 0), (0, 0), (0.500000000000000, 0.500000000000000), (0, 0)]
        sage: pts = [Point(-2,0), Point(0,-2), Point(-0.5,-0.5), Point(-2,-2)]
        sage: reduce_in_parallelogram(p, pts)
        [(0, 0), (0, 0), (0.500000000000000, 0.500000000000000), (0, 0)]
        sage: v1 = (3.46410161513775, 0.000000000000000)
        sage: v2 = (0.000000000000000, 1.00000000000000)
        sage: par = Parallelogram(base=(-1.73205080756888, -0.500000000000000), v1=v1, v2=v2)
        sage: reduce_in_parallelogram(par, [vector((3,0))]) # abs tol 8e-15
        [(-0.464101615137754, 0.000000000000000)]
        sage: reduce_in_parallelogram(par, [(3,0)]) # abs tol 8e-15
        [(-0.464101615137754, 0.000000000000000)]

    """
    V = (p.v1.base_ring() ** 2).span_of_basis((p.v1, p.v2))
    base_coords = V.coordinates(p.base)
    reduced_points = []
    translations = []
    if not isinstance(points, (list, tuple)):
        points = [points]
    for point in points:
        if not isinstance(point, Vector):
            point = vector(point)
        point_coords = V.coordinates(point)
        t0 = floor(point_coords[0] - base_coords[0])
        t1 = floor(point_coords[1] - base_coords[1])
        reduced_points.append(point - t0 * p.v1 - t1 * p.v2)
        if return_translations:
            translations.append((t0, t1))
    if return_translations:
        return reduced_points, translations
    return reduced_points


def reduce_list_to_parallelogram(p: Parallelogram, word_circle_list: list[tuple[str, Circle]]):
    """
    Reduce a list of circles by vectors spanning the parallelogram into a
     fundamental domain corresponding to the given base of the parallelogram.

    INPUT:
    - ``p`` - Parallelogram
    - ``word_circle_list`` - list of circles

    EXAMPLES:

        sage: from maass_forms_klein.hyperbolic_space.parallelogram import reduce_list_to_parallelogram
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: from maass_forms_klein.hyperbolic_space.types import Point
        sage: p = Parallelogram(base=(0,0), v1=(0,1), v2=(1,0))
        sage: circle_list = [(0, Circle(center=Point(0,0), radius=1))]
        sage: reduce_list_to_parallelogram(p, circle_list)
        [('0', Circle(center=(0, 0), radius=1))]
    """
    centers, translations = reduce_in_parallelogram(p, [x[1].center for x in word_circle_list],
                                                                 True)
    return [
        (f"{x[0]}{translation_tuple_to_word(translations[n], {0: 'L', 1: 'M'})}",
         Circle(center=centers[n], radius=x[1].radius)) for n, x in enumerate(word_circle_list)]


def point_on_boundary_of_parallelogram(point: Vector, p: Parallelogram, verbose: bool = False)\
        -> bool:
    """
    Check if point is on the boundary of the parallelogram.

    :param point:
    :param p:
    :return:

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.utils import point_on_boundary_of_parallelogram
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
        sage: v1 = (3.46410161513775, 0.000000000000000)
        sage: v2 = (0.000000000000000, 1.00000000000000)
        sage: par = Parallelogram(base=vector((0, 0)), v1=v1, v2=v2)
        sage: point_on_boundary_of_parallelogram(vector((5, 1.0)), rect)
        False
        sage: point_on_boundary_of_parallelogram(vector((2, 1.5)), rect)
        False
        sage: point_on_boundary_of_parallelogram(vector((2, 1.0)), rect)
        True

    """
    V = (p.v1.base_ring() ** 2).span_of_basis((p.v1, p.v2))
    base_coords = V.coordinates(p.base)
    if not isinstance(point, Vector):
        point = vector(point)
    point_coords = V.coordinates(point)
    x0 = point_coords[0] - base_coords[0]
    x1 = point_coords[1] - base_coords[1]
    prec = p.v1.base_ring().epsilon() * 8
    if min(abs(x0), abs(x0 - 1)) <= prec and -prec <= x1 <= 1 + prec:
        return True
    if min(abs(x1), abs(x1 - 1)) <= prec and -prec <= x0 <= 1 + prec:
        return True
    return False


def parallelogram_intersect_circle(p: Parallelogram, cr: Circle,
                                   scaling_factor: Real_t = 1.0) -> bool:
    """
    Check if a given parallelogram intersects a given circle (possibly scaled to account for
    numerical errors).

    INPUT:
    - ``cr`` - Circle
    - ``p`` - Parallellogram
    - ``scaling_factor`` - scaling the circle this amount.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
        sage: from maass_forms_klein.hyperbolic_space.utils import parallelogram_intersect_circle
        sage: p = Parallelogram(base=vector((0, 0)), v1=(3.46410161513775, 0.000000000000000), v2=(0.000000000000000, 1.00000000000000))
        sage: cr = Circle(center=vector((0, 0)), radius=1)
        sage: parallelogram_intersect_circle(p, cr)
        True
        sage: cr = Circle(center=vector((0, 0)), radius=100)
        sage: parallelogram_intersect_circle(p, cr)
        True
        sage: cr = Circle(center=vector((0.1, 0.1)), radius=0.01)
        sage: parallelogram_intersect_circle(p, cr)
        True
        sage: base = vector((0.0895715621873992, -0.147916720404305))
        sage: v1 = vector((-0.168820510688814, 0.281139756031150))
        sage: v2 = vector((-0.0103226136859846, 0.0146936847774607))
        sage: p = Parallelogram(base=base, v1=v1, v2=v2)
        sage: r = 0.00432823543907501
        sage: cr = Circle(center=(0.0700167159440167, -0.110113017437223), radius=r)
        sage: parallelogram_intersect_circle(p, cr)
        True
        sage: p = Parallelogram(base=vector((0, 0)), v1=(2.0, 0.0), v2=(0.0, 0.5))
        sage: cr = Circle(center=vector((1.0, 0.5)), radius=0.51)
        sage: parallelogram_intersect_circle(p, cr)
        True

    """
    # Check if the c + r * n is in the parallelogram
    # where n is in the direction of the orthogonal projection onto the spanning vectors
    n1 = perpendicular(p.v1)
    n2 = perpendicular(p.v2)
    r = cr.radius * scaling_factor
    circle_intersects_parallelogram = any(
        any(point_in_parallelogram(cr.center + ni * ri, p)
            for ri in [r, -r]) for ni in [n1, n2])
    if circle_intersects_parallelogram:
        return True
    vertices = p.base, p.base + p.v1, p.base + p.v2, p.base + p.v1 + p.v2
    from maass_forms_klein.hyperbolic_space.geometry_utils import point_in_circle
    parallelogram_in_circle = any(point_in_circle(pt, cr, scaling_factor=scaling_factor)
                                  for pt in vertices)
    if parallelogram_in_circle:
        return True
    # Quick check for exclusion
    return not circle_outside_parallelogram(cr, p, scaling_factor=scaling_factor)
    # Else do a brute-force check
    twopi = RR.pi() * 2
    pts = [
        cr.center + vector((cos(twopi * n / 50),
                           sin(twopi * n / 50))) * r for n in range(50)]
    for pt in pts:
        if point_in_parallelogram(pt, p):
            return True
    return False


def point_in_parallelogram(c: Vector, p: Parallelogram) -> bool:
    """
    Check if the point `c` is inside the parallelogram spanned by `v1` and `v2`

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
        sage: from maass_forms_klein.hyperbolic_space.parallelogram import point_in_parallelogram
        sage: cr = Circle(center=vector((0,0)), radius=1)
        sage: p = Parallelogram(base=vector((-0.5,-0.5)), v1=vector((1,0)), v2=vector((0,1)))
        sage: point_in_parallelogram(cr.center, p)
        True
        sage: point_in_parallelogram(cr.center + vector((0.5,0.5)), p)
        True
        sage: point_in_parallelogram(cr.center + vector((1,1)), p)
        False
        sage: point_in_parallelogram(cr.center + vector((1,-1)), p)
        False
        sage: point_in_parallelogram(cr.center + vector((-1,1)), p)
        False

    """
    cx, cy = p.coordinates(c)
    return -1/2 <= cx <= 1/2 and -1/2 <= cy <= 1/2

    L1, L2, L3, L4 = parallelogram_to_lines(p)
    a1 = orthogonal_projection_on_line(c, L1)
    a2 = orthogonal_projection_on_line(c, L2)
    b1 = orthogonal_projection_on_line(c, L3)
    b2 = orthogonal_projection_on_line(c, L4)
    in1 = points_between(a1, c, a2) or points_between(a2, c, a1)
    in2 = points_between(b1, c, b2) or points_between(b2, c, b1)
    return in1 and in2

def circle_outside_parallelogram(c: Circle, p: Parallelogram, scaling_factor: Real_t = 1) -> bool:
    """
    Return True if circle is outside of parallelogram and False if it can not be determined
    (in particular ``False`` does not necessarily imply that the circle intersects the parallelogram).

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.types import Circle
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
        sage: from maass_forms_klein.hyperbolic_space.parallelogram import circle_outside_parallelogram
        sage: cr = Circle(center=vector((0,0)), radius=1)
        sage: p = Parallelogram(base=vector((-0.5,-0.5)), v1=vector((1,0)), v2=vector((0,1)))
        sage: circle_outside_parallelogram(cr, p)
        False
        sage: circle_outside_parallelogram(cr, p)
        False
        sage: c = Circle(center=(3.8888888888888893, -5.0), radius=0.200000000000000)
        sage: v1 = vector((-7.28367233850816, 8.31121964982996))
        sage: v2 = vector((3.46820643039263, -3.44937511512084))
        sage: p = Parallelogram(base=(3.8888888888888893, -5.0), v1=v1, v2=v2)
        sage: circle_outside_parallelogram(c, p)
        False
        """
    # Quickest check first: see if the circle is disjoint to the circumscribed circle of the parallelogram
    from maass_forms_klein.hyperbolic_space.geometry_utils import circle_intersects_circle
    if not circle_intersects_circle(c, p.circumscribed_circle):
        return True

    n1 = perpendicular(p.v1)
    n2 = perpendicular(p.v2)
    r = c.radius * scaling_factor
    # Check if all circle projections on the directions perpendicular to the
    # parallelogram basis are on the same side of the outside of the parallelogram
    coord11 = p.coordinates(c.center + n1 * r)
    coord12 = p.coordinates(c.center - n1 * r)
    for i in range(2):
        if coord11[i] < -1 / 2 and coord12[i] < -1 / 2:
            return True
        if coord11[i] > 1 / 2 and coord12[i] > 1 / 2:
            return True
    coord21 = p.coordinates(c.center + n2 * r)
    coord22 = p.coordinates(c.center - n2 * r)
    for i in range(2):
        if coord21[i] < -1 / 2 and coord22[i] < -1 / 2:
            return True
        if coord21[i] > 1 / 2 and coord22[i] > 1 / 2:
            return True
    # If one coordinate is in [-1/2, 1/2] and the other is not completely outside then we do intersect
    if abs(coord11[0]) <= 1 / 2 and abs(coord12[0]) <= 1 / 2  and \
            (coord11[1]+1/2) * (coord12[1]+1/2) <= 0:
        return False
    if abs(coord21[0]) <= 1 / 2 and abs(coord22[0]) <= 1 / 2  and \
            (coord21[1]+1/2) * (coord22[1]+1/2) <= 0:
        return False

    #if abs(coord11[0]) <= 1 / 2 or abs(coord11[1]) <= 1 / 2:
    # If the circle is not entirely in a hyperplane then it can only intersect the parallelogram
    # if it contains one of the vertices.
    # from maass_forms_klein.hyperbolic_space.geometry_utils import point_in_circle
    # if not any(point_in_circle(v, c, scaling_factor) for v in p.vertices):
    #     return True
    return False

def perpendicular(v: Vector) -> Vector:
    """
    Calculate a normalised vector perpendicular to v.

    INPUT:

    - ``v`` - Vector

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.types import Vector
        sage: from maass_forms_klein.hyperbolic_space.parallelogram import perpendicular
        sage: v = vector((1, 0))
        sage: perpendicular(v)
        (0, 1)
        sage: v = vector((1, 1))
        sage: perpendicular(v)
        (1/2*sqrt(2), -1/2*sqrt(2))
        sage: v = vector((1.0, 1.0))
        sage: perpendicular(v)
        (0.707106781186547, -0.707106781186547)
    """
    if v[1] == 0 and v[0] == 0:
        raise ValueError("Zero vector")
    elif v[1] == 0:
        return vector((v[1], 1))
    n1x = 1
    n1y = -v[0] / v[1]
    norm = sqrt(1 + n1y * n1y)
    return vector((1 / norm, n1y / norm))

    # if v[1] == 0:
    #     return vector((0, 1))
    # n1 = vector((1, -v[0] / v[1]))
    # n1 = n1 / n1.norm()
    # return n1


def orthogonal_projection_on_line(a: Vector, L: Line) -> Vector:
    """
    Compute the orthogonal projection of a point on a line.
    """
    v = L.direction
    n = perpendicular(v)
    return intersection_point_two_lines(Line(base=a, direction=n), L)


def parallelogram_to_lines(p: Parallelogram) -> tuple[Line, Line, Line, Line]:
    """
    Find the lines given by the sides of a parallelogram.

    EXAMPLES:
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
        sage: from maass_forms_klein.hyperbolic_space.parallelogram import parallelogram_to_lines
        sage: p = Parallelogram(vector((0, 0)), vector((1, 0)), vector((0, 1)))
        sage: l1, l2, l3, l4 = parallelogram_to_lines(p)
        sage: l1
        Line(base=(0, 0), direction=(1, 0))
        sage: l2
        Line(base=(0, 1), direction=(1, 0))
        sage: l3
        Line(base=(0, 0), direction=(0, 1))
        sage: l4
        Line(base=(1, 0), direction=(0, 1))
    """
    v1 = p.v1
    v2 = p.v2
    l1 = Line(base=p.base, direction=v1)
    l2 = Line(base=p.base + v2, direction=v1)
    l3 = Line(base=p.base, direction=v2)
    l4 = Line(base=p.base + v1, direction=v2)
    return l1, l2, l3, l4


def intersection_point_two_lines(l1: Line, l2: Line) -> Vector:
    a = l1.base
    v = l1.direction
    c = l2.base
    n = l2.direction
    if n[0] != 0:
        t = (n[0] * (c[1] - a[1]) + n[1] * (a[0] - c[0])) / (n[0] * v[1] - n[1] * v[0])
    else:
        t = (c[0] - a[0]) / v[0]
    return a + t * v


def points_between(A: Vector, B: Vector, C: Vector) -> bool:
    """
    Check if B is on a straight line between A and C.

    EXAMPLES:

        sage: from maass_forms_klein.hyperbolic_space.types import vector
        sage: from maass_forms_klein.hyperbolic_space.parallelogram import points_between
        sage: A = vector((-1.00000000000000, -0.0773502691896260))
        sage: B = vector((3.17542648054294, -0.0773502691896260))
        sage: C = vector((0.000000000000000, -0.0773502691896260))
        sage: points_between(A, B, C)
        False
    """
    return ((A[0] <= B[0] <= C[0] or C[0] <= B[0] <= A[0]) and
            (A[1] <= B[1] <= C[1] or C[1] <= B[1] <= A[1]))


def parallelogram_in_circle(cr: Circle, p: Parallelogram, scaling_factor: Real_t = 1) -> bool:
    """
    Check if the parallelogram spanned by v1 and v2 is contained in a given circle


    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.types import Circle, Parallelogram
        sage: from maass_forms_klein.hyperbolic_space.parallelogram import parallelogram_in_circle
        sage: cr = Circle(center=vector((0,0)), radius=1)
        sage: p = Parallelogram(base=vector((0,0)), v1=vector((1,0)), v2=vector((0,1)))
        sage: parallelogram_in_circle(cr, p)
        False
        sage: parallelogram_in_circle(cr, p, scaling_factor=1.42)
        True
        sage: cr = Circle(center=vector((0.5,0.5)), radius=1)
        sage: parallelogram_in_circle(cr, p)
        True

    """

    v1 = p.v1
    v2 = p.v2
    return all(sqrt((corner[0] - cr.center[0]) ** 2 + (corner[1] - cr.center[1]) ** 2)
               <= cr.radius * scaling_factor
               for corner in [p.base, p.base+v1, p.base + v2, p.base + v1 + v2])


def split_parallelogram(p: Parallelogram | Rectangle, n1: Integer_t = 2, n2: Integer_t = 2):
    """
    Split a parallelogram into a list of sub-parallelograms of size n1 x n2

    INPUT:

    - ``p`` -- parallelogram
    - ``n1`` -- number of sub-parallelograms in the first (x-) dimension
    - ``n2`` -- number of sub-parallelograms in the second (y-) dimension

    EXAMPLE:
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
        sage: from maass_forms_klein.hyperbolic_space.parallelogram import split_parallelogram
        sage: p = Parallelogram(base=vector((0,0)), v1=vector((1,0)), v2=vector((0,1)))
        sage: split_parallelogram(p, n1=1, n2=2)
        [Parallelogram(base=(0, 0), v1=(1, 0), v2=(0, 1/2)),
         Parallelogram(base=(0, 1/2), v1=(1, 0), v2=(0, 1/2))]
    """
    v1_split = p.v1 / n1
    v2_split = p.v2 / n2
    return [Parallelogram(v1=v1_split, v2=v2_split, base=h1 * v1_split + h2 * v2_split + p.base)
            for h1, h2 in cartesian_product([range(n1), range(n2)])]


def parallelogram_covered_by_matrices(p: Parallelogram, matrices: list[Matrix],
                                      scaling_factor=1.0, n_max=1000):
    if not isinstance(matrices, list) or not all(isinstance(x, Matrix) for x in matrices):
        raise ValueError(
            f"Input should be a list of matrices. Got: {matrices} type={type(matrices[0])}")
    from maass_forms_klein.hyperbolic_space.geometry_utils import matrix_to_circle
    circles = [matrix_to_circle(m) for m in matrices]
    return parallelogram_covered_by_circles(p, circles, scaling_factor=scaling_factor,
                                            n_max=n_max)


def parallelogram_covered_by_circles(p: Parallelogram, circles: list[Circle],
                                     scaling_factor: Real_t = 1.0,
                                     verbose: int = 0,
                                     n_min: Integer_t = 2,
                                     n_max: Integer_t = 1000):
    """
    Return True if parallelogram is covered by the circles scaled with `scaling_factor`.

    INPUT:

    - ``p`` -- parallelogram
    - ``circles`` -- list of circles

    EXAMPLES:
        sage: from maass_forms_klein.hyperbolic_space.types import Circle, Parallelogram
        sage: from maass_forms_klein.hyperbolic_space.parallelogram import parallelogram_covered_by_circles
        sage: cr = Circle(center=vector((0,0)), radius=1)
        sage: p = Parallelogram(base=vector((0,0)), v1=vector((1,0)), v2=vector((0,1)))
        sage: parallelogram_covered_by_circles(p, [cr])
        False
        sage: parallelogram_covered_by_circles(p, [cr], scaling_factor=2)
        True

    TODO: KMake it a bit smarter and only bisect the rectangles that are not covered
    """
    if not isinstance(circles, list) or not all(isinstance(x, Circle) for x in circles):
        raise ValueError(f"Input should be a list of circles. Got: {circles} type={type(circles[0])}")
    if not circles:
        return True  # Empty list of circles covers nothing (vacuously true for any condition)
    max_radius = max([c.radius for c in circles])
    if max_radius == Infinity:
        raise ValueError("Matrix has infinite invariant circle")
    # if the parallelogram does not intersect any circle then it is not covered
    if not any(parallelogram_intersect_circle(p, c, scaling_factor=scaling_factor)
               for c in circles):
        return False
    # If the parallelogram is covered by one of the circles then we are done.
    if any(parallelogram_in_circle(c, p, scaling_factor=scaling_factor) for c in circles):
        return True
    n = 1
    # If the corners of the parallelogram are not contained in the circles scaled by 0.999
    # then this is a quick out as it is either not covered at all or it is just n the boundary.
    # Note: for numerical stability we need a slightly larger cover.
    from maass_forms_klein.hyperbolic_space.geometry_utils import point_covered_by_circles
    if not all(point_covered_by_circles(v, circles,
                                        scaling_factor=scaling_factor * 0.9999)
               for v in p.vertices):
        return False
    ratio_x = max(1, int(p.sides[0] / p.sides[1]))
    ratio_y = max(1, int(p.sides[1] / p.sides[0]))
    not_covered = []
    n = n_min
    # Check if we have too small rectangle:
    mind = 0.001
    for c1 in circles:
        for c2 in circles:
            if c2 == c1:
                continue
            d = abs(abs(c1.center - c2.center) - (c1.radius + c2.radius))
            if 1e-10 < d < mind:
                mind = d
    if verbose > 0:
        print("-"* 30)
        print("Check parallelogram cover for:", p)
        print("mind=", mind, "ratio", ratio_x, ratio_y)
        print("n_min=", n, "n_max=", n_max)
    while n <= n_max:
        if p.sides[0] < n * ratio_x * mind or p.sides[1] < n * ratio_y * mind:
            print("too small:", p.sides[0], n * ratio_x * mind, p.sides[1], n * ratio_y * mind)
            break
        all_are_covered = True
        if verbose:
            print("-="*10, n * ratio_x, n * ratio_y)
        split = split_parallelogram(p, n1=n * ratio_x, n2=n * ratio_y)
        for para_small in split:

            # If the small parallelogram does not intersect any of the circles
            # then it is clearly not covered
            if not any(parallelogram_intersect_circle(para_small, c, scaling_factor=scaling_factor)
                       for c in circles):
                if verbose:
                    print("not any")
                return False
            # If it is not covered by any circle then we need to divide further
            if not any(parallelogram_in_circle(c, para_small, scaling_factor=scaling_factor)
                       for c in circles):
                if verbose:
                    print(para_small)
                    print("all are not covered")
                    print("circles=", circles)
                    checks = [parallelogram_in_circle(c, para_small, scaling_factor=scaling_factor) for c in circles]
                    print("in:", checks)

                all_are_covered = False
                not_covered.append(para_small)
                #break
        if all_are_covered:
            if verbose:
                print("all are covered")
            return True
        all_are_covered = True
        if verbose:
            print("not_covered=", len(not_covered))
        for para_small in not_covered:
            if verbose:
                print("Check Not covered:", para_small)
            test = parallelogram_covered_by_circles(para_small, circles,
                                                    scaling_factor=scaling_factor,
                                                    verbose=verbose,
                                                    n_max=2)
            if not test:
                all_are_covered = False
                break
        if all_are_covered:
            if verbose:
                print("all are covered")
            return True
        n += 1
    raise ArithmeticError("Could not determine if rectangle is covered in given number of steps.")
