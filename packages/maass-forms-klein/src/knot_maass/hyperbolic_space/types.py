from dataclasses import dataclass, field, astuple
from typing import NoReturn

from knot_maass.modform.utils import Real_t
from sage.matrix.constructor import matrix
from sage.misc.functional import sqrt
from sage.modules.free_module_element import vector
from sage.structure.element import Vector

@dataclass(frozen=True)
class Point:
    x: Real_t
    y: Real_t

    def __iter__(self):
        yield self.x
        yield self.y

@dataclass(eq=True, unsafe_hash=True)
class Parallelogram:
    base: Vector
    v1: Vector
    v2: Vector

    def __post_init__(self) -> NoReturn:
        self.base = vector(self.base)
        self.v1 = vector(self.v1)
        self.v2 = vector(self.v2)
        self.base.set_immutable()
        self.v1.set_immutable()
        self.v2.set_immutable()
        self.vector_space = (self.v1.base_ring() ** 2).span_of_basis((self.v1, self.v2))
        base_coords = self.vector_space.coordinates(self.base)
        # Set coordinates to go between -1/2 and 1/2
        self.base_coords = [base_coords[0] + 1/2, base_coords[1] + 1/2]
        radius = max((self.v1 + self.v2).norm(), (self.v1 - self.v2).norm()) / 2
        self.circumscribed_circle = Circle(self.center, radius)
    @property
    def center(self):
        """
        Center of a parallelogram.

        """
        return self.base + self.v1 / 2 + self.v2 / 2

    @property
    def sides(self):
        """
        Side lengths of a parallelogram

        """
        return [self.v1.norm(), self.v2.norm()]

    @property
    def vertices(self):
        return [self.base, self.base + self.v1, self.base + self.v1 + self.v2, self.base + self.v2]

    @property
    def area(self):
        return matrix([self.v1, self.v2]).determinant()

    def point_on_boundary(self, point: Vector) -> bool:
        """
        Check if a point is on the boundary of the parallelogram.
        """
        x, y = point
        return x == self.base[0] or x == self.base + self.v1 or x == self.base + self.v1 + self.v2 or x == self.base + self.v2

    def to_json(self):
        return {
            "base": list(self.base),
            "v1": list(self.v1),
            "v2": list(self.v2),
        }

    def coordinates(self, x: Vector | Point) -> list:
        coords = self.vector_space.coordinates(x)
        return [coords[0] - self.base_coords[0], coords[1] - self.base_coords[1]]


@dataclass(eq=True, unsafe_hash=True, order=True)
class Circle:
    center: Vector = field(init=True)
    radius: Real_t

    def __post_init__(self) -> NoReturn:
        self.center = vector(self.center)
        self.center.set_immutable()

    def scale(self, scale_factor: Real_t) -> "Circle":
        return Circle(center=self.center, radius=self.radius * scale_factor)

    def within_epsilon(self, other: "Circle", epsilon: Real_t) -> bool:
        """
        Check if self is within epsilon of other circle
        """
        diff_radius = abs(self.radius - other.radius)
        if diff_radius > epsilon:
            return False
        diff_center = sqrt((other.center[0] - self.center[0]) ** 2 +
                           (other.center[1] - self.center[1]) ** 2)
        return diff_center <= epsilon
        # return abs(self.center - other.center) < epsilon


@dataclass(frozen=True, eq=True)
class Line:
    base: Vector
    direction: Vector

    def __post_init__(self) -> NoReturn:
        self.base.set_immutable()
        self.direction.set_immutable()


@dataclass
class Rectangle(Parallelogram):
    base: Vector
    v1: Vector
    v2: Vector

    def __post_init__(self) -> NoReturn:
        self.base.set_immutable()
        self.v1.set_immutable()
        self.v2.set_immutable()

