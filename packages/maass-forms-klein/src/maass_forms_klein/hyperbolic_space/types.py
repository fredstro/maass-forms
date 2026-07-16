from dataclasses import dataclass, field
from typing import NoReturn, Union

# Define Real_t locally to avoid import chain issues
from sage.rings.real_mpfr import RealNumber
from sage.rings.integer import Integer
from sage.rings.rational import Rational

from sage.matrix.constructor import matrix
from sage.misc.functional import sqrt
from sage.modules.free_module import FreeModule_ambient
from sage.modules.free_module_element import vector
from sage.structure.element import Vector

Real_t = Union[RealNumber, Integer, Rational, int, float]


@dataclass(frozen=True)
class Point:
    x: Real_t
    y: Real_t

    def __iter__(self):
        """
        Iterate over the coordinates of the point.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Point
            sage: p = Point(1, 2)
            sage: list(p)
            [1, 2]
        """
        yield self.x
        yield self.y


@dataclass(eq=True, unsafe_hash=True)
class Parallelogram:
    base: Vector
    v1: Vector
    v2: Vector
    vector_space: FreeModule_ambient = field(init=False)

    def __post_init__(self) -> NoReturn:
        """Initialize internal fields from base, v1, and v2.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
            sage: P = Parallelogram(base=(0, 0), v1=(1, 0), v2=(0, 1))
            sage: P.base
            (0, 0)
        """
        self.base = vector(self.base)
        self.v1 = vector(self.v1)
        self.v2 = vector(self.v2)
        self.base.set_immutable()
        self.v1.set_immutable()
        self.v2.set_immutable()
        self.vector_space = (self.v1.base_ring() ** 2).span_of_basis((self.v1, self.v2))
        base_coords = self.vector_space.coordinates(self.base)
        # Set coordinates to go between -1/2 and 1/2
        self.base_coords = [base_coords[0] + 1 / 2, base_coords[1] + 1 / 2]
        radius = max((self.v1 + self.v2).norm(), (self.v1 - self.v2).norm()) / 2
        self.circumscribed_circle = Circle(self.center, radius)

    @property
    def center(self):
        """
        Center of a parallelogram.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
            sage: from sage.modules.free_module_element import vector
            sage: P = Parallelogram(base=vector([0, 0]), v1=vector([2, 0]), v2=vector([0, 2]))
            sage: P.center
            (1, 1)
        """
        return self.base + self.v1 / 2 + self.v2 / 2

    @property
    def sides(self):
        """
        Side lengths of a parallelogram.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
            sage: from sage.modules.free_module_element import vector
            sage: P = Parallelogram(base=vector([0, 0]), v1=vector([3, 0]), v2=vector([0, 4]))
            sage: P.sides
            [3, 4]
        """
        return [self.v1.norm(), self.v2.norm()]

    @property
    def vertices(self):
        """
        Return the four vertices of the parallelogram.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
            sage: from sage.modules.free_module_element import vector
            sage: P = Parallelogram(base=vector([0, 0]), v1=vector([1, 0]), v2=vector([0, 1]))
            sage: P.vertices
            [(0, 0), (1, 0), (1, 1), (0, 1)]
        """
        return [self.base, self.base + self.v1, self.base + self.v1 + self.v2, self.base + self.v2]

    @property
    def area(self):
        """
        Return the signed area of the parallelogram.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
            sage: from sage.modules.free_module_element import vector
            sage: P = Parallelogram(base=vector([0, 0]), v1=vector([3, 0]), v2=vector([0, 4]))
            sage: P.area
            12
        """
        return matrix([self.v1, self.v2]).determinant()

    def point_on_boundary(self, point: Vector) -> bool:
        """
        Check if a point is on the boundary of the parallelogram.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
            sage: from sage.modules.free_module_element import vector
            sage: P = Parallelogram(base=vector([0, 0]), v1=vector([1, 0]), v2=vector([0, 1]))
            sage: P.point_on_boundary(vector([0, 0]))
            True
        """
        x, _y = point
        return (
            x == self.base[0]
            or x == self.base + self.v1
            or x == self.base + self.v1 + self.v2
            or x == self.base + self.v2
        )

    def to_json(self):
        """
        Return a JSON-serializable dictionary representation.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
            sage: from sage.modules.free_module_element import vector
            sage: P = Parallelogram(base=vector([0, 0]), v1=vector([1, 0]), v2=vector([0, 1]))
            sage: P.to_json()
            {'base': [0, 0], 'v1': [1, 0], 'v2': [0, 1]}
        """
        return {
            "base": list(self.base),
            "v1": list(self.v1),
            "v2": list(self.v2),
        }

    def coordinates(self, x: Vector | Point) -> list:
        """
        Return the coordinates of ``x`` relative to the parallelogram.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
            sage: from sage.modules.free_module_element import vector
            sage: P = Parallelogram(base=vector([0, 0]), v1=vector([2, 0]), v2=vector([0, 2]))
            sage: P.coordinates(vector([1, 1]))  # doctest: +SKIP
            [0.0, 0.0]
        """
        coords = self.vector_space.coordinates(x)
        return [coords[0] - self.base_coords[0], coords[1] - self.base_coords[1]]

    def __repr__(self) -> str:
        """
        Custom repr that excludes vector_space field.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram
            sage: from sage.modules.free_module_element import vector
            sage: P = Parallelogram(base=vector([0, 0]), v1=vector([1, 0]), v2=vector([0, 1]))
            sage: P
            Parallelogram(base=(0, 0), v1=(1, 0), v2=(0, 1))
        """
        return f"Parallelogram(base={self.base}, v1={self.v1}, v2={self.v2})"


@dataclass(eq=True, unsafe_hash=True, order=True)
class Circle:
    center: Vector = field(init=True)
    radius: Real_t

    def __post_init__(self) -> NoReturn:
        """Initialize and make center vector immutable.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Circle
            sage: c = Circle(center=(0, 0), radius=1)
            sage: c.center
            (0, 0)
        """
        self.center = vector(self.center)
        self.center.set_immutable()

    def scale(self, scale_factor: Real_t) -> "Circle":
        """
        Return a new circle scaled by the given factor.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Circle
            sage: from sage.modules.free_module_element import vector
            sage: c = Circle(center=vector([0, 0]), radius=1)
            sage: c.scale(2).radius
            2
        """
        return Circle(center=self.center, radius=self.radius * scale_factor)

    def within_epsilon(self, other: "Circle", epsilon: Real_t) -> bool:
        """
        Check if self is within epsilon of other circle.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Circle
            sage: from sage.modules.free_module_element import vector
            sage: c1 = Circle(center=vector([0, 0]), radius=1)
            sage: c2 = Circle(center=vector([0, 0]), radius=1)
            sage: c1.within_epsilon(c2, 0.01)
            True
        """
        diff_radius = abs(self.radius - other.radius)
        if diff_radius > epsilon:
            return False
        diff_center = sqrt(
            (other.center[0] - self.center[0]) ** 2 + (other.center[1] - self.center[1]) ** 2
        )
        return diff_center <= epsilon
        # return abs(self.center - other.center) < epsilon


@dataclass
class Face:
    r"""
    A face of the Ford domain: a visible isometric hemisphere together with its
    combinatorial data.

    A face is carried by the isometric hemisphere of a group element ``g``; the
    hemisphere projects to the disc ``circle`` in the boundary plane. ``word`` is
    a word representing ``g`` in the generators, ``inverse_word`` the word of the
    paired face (carried by the hemisphere of ``g^{-1}``, i.e. the class matched
    modulo the translation lattice), ``witness`` a boundary point where the face
    is visible on the Ford floor, and ``margin`` the strict-visibility margin
    (positive for a genuine face).

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.types import Face, Circle
        sage: from sage.modules.free_module_element import vector
        sage: f = Face(word='b', circle=Circle(center=vector((0, 0)), radius=1),
        ....:          inverse_word='B')
        sage: f.word
        'b'
        sage: f.inverse_word
        'B'
        sage: f.circle.radius
        1
    """

    word: str
    circle: Circle
    inverse_word: str | None = None
    witness: Vector | None = None
    margin: Real_t | None = None

    def __repr__(self) -> str:
        r"""
        Concise representation of a face.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Face, Circle
            sage: from sage.modules.free_module_element import vector
            sage: Face(word='b', circle=Circle(center=vector((0, 0)), radius=1),
            ....:      inverse_word='B')
            Face(word='b', inverse_word='B', radius=1)
        """
        return (
            f"Face(word={self.word!r}, inverse_word={self.inverse_word!r}, "
            f"radius={self.circle.radius})"
        )


@dataclass(frozen=True, eq=True)
class Line:
    base: Vector
    direction: Vector

    def __post_init__(self) -> NoReturn:
        """Make base and direction vectors immutable.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Line
            sage: L = Line(base=vector((0, 0)), direction=vector((1, 0)))
            sage: L.base
            (0, 0)
        """
        self.base.set_immutable()
        self.direction.set_immutable()


@dataclass
class Rectangle(Parallelogram):
    base: Vector
    v1: Vector
    v2: Vector

    def __post_init__(self) -> NoReturn:
        """Initialize rectangle via parent Parallelogram.

        EXAMPLES::

            sage: from maass_forms_klein.hyperbolic_space.types import Rectangle
            sage: R = Rectangle(base=(0, 0), v1=(1, 0), v2=(0, 1))
            sage: R.base
            (0, 0)
        """
        super().__post_init__()
