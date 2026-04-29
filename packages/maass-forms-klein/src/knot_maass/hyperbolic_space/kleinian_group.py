r"""
Class representing a Kleinian group, i.e. a discrete subgroup of PSL(2,C)
acting discretely on the upper half-space.

"""
import json
import string
from typing import ParamSpec

from knot_maass.hyperbolic_space.utils import find_covering_generators
from knot_maass.hyperbolic_space.upper_half_space import UpperHalfSpaceElement__class
from math import ceil, floor

from knot_maass.hyperbolic_space.utils import Integer_t, Real_t, get_lattice_values
from knot_maass.hyperbolic_space.geometry_utils import split_rectangle, rectangle_in_circle, \
    matrix_to_circle
from knot_maass.hyperbolic_space.word_utils import find_inverse_word, normalize_word, expand_word, \
    word_to_element
from knot_maass.utils.json_converters import matrix_to_json, matrix_from_json
from sage.categories.groups import Groups
from sage.functions.other import imag, real
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.cachefunc import cached_method
from sage.modular.cusps_nf import NFCusp
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexField, ComplexNumber
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.qqbar import QQbar
from sage.rings.real_mpfr import RealField
from sage.structure.element import Matrix, Element, Vector
from sage.groups.matrix_gps.linear import LinearMatrixGroup_generic
from sage.rings.number_field.number_field_base import NumberField

from knot_maass.hyperbolic_space.types import Rectangle, Parallelogram
from snappy import Manifold

from knot_maass.hyperbolic_space.word_utils import change_letters_in_word

from knot_maass.hyperbolic_space.utils import find_side_pairing_gens

P = ParamSpec('P')


# noinspection PyArgumentList
class KleinianGroup_class(LinearMatrixGroup_generic):
    r"""
    Class representing a Kleinian group, i.e. a discrete subgroup of PSL(2,C)
    acting discretely on the upper half-space.

    TODO: write examples.
    TODO: interface from snappy
    """
    def __init__(self, x: Integer_t | tuple | list = None, **kwargs: P.kwargs) -> None:
        r"""

        INPUT:
        - ``x`` -- integer, order in imaginary quadratic number field, list of generators

        NOTES: Generators are represented as complex matrices.
        TODO: Do we need to separate parabolic and hyperbolic?
        EXAMPLE:

            sage: from knot_maass.hyperbolic_space.kleinian_group import KleinianGroup_class
            sage: K = KleinianGroup_class(-4)
            sage: K
            Bianchi Group: Q(sqrt(-4))
        """
        data = kwargs.get('data', {})
        if isinstance(data, str):
            data = json.loads(data)
        self._gens = [matrix_from_json(g) for g in data.get('_gens', [])]
        self._named_gens = {k: matrix_from_json(g) for k, g in
                               data.get('_named_gens', {}).items()}
        self._covering_generators = {k: matrix_from_json(g) for k, g in
                                        data.get('_covering_generators', {}).items()}
        self._named_gens = {k: matrix_from_json(g) for k, g in
                            data.get('_named_gens', {}).items()}
        self._side_pairing_gens = data.get('_side_pairing_gens', [])
        self._covering_generators_words = data.get('_covering_generators_words', [])
        self._prec = data.get('_prec', None)
        self._init_prec = None
        basis_matrix = data.get('_translation_lattice', None)
        if basis_matrix:
            basis_matrix = matrix_from_json(basis_matrix)
            prec = basis_matrix.parent().base_ring().prec()
            self._translation_lattice = (RealField(prec)**2).span_of_basis(basis_matrix)
        else:
            self._translation_lattice = None
        self._manifold = kwargs.get('manifold', data.get('_manifold', ''))
        self._name_string = data.get('_name_string', '')
        self._latex_string = data.get('_latex_string', '')
        if self._manifold and not self._name_string:
            self._name_string = f'KleinianGroup from Knot: {self._manifold}'
        if self._manifold and not self._latex_string:
            self._latex_string = f'Kleinian Group from Knot: ${self._manifold}$'
        if self._gens:
            base_ring = self._gens[0].base_ring()
        elif isinstance(x, (int, Integer)) and x < 0 and ZZ(x).is_fundamental_discriminant():
            base_ring = QuadraticField(x)
            if base_ring.class_number() > 1:
                raise NotImplementedError("Only class number 1 is supported.")
            prec = kwargs.get('prec', 53)
            MS = MatrixSpace(ComplexField(prec), 2, 2)
            b1, b2 = base_ring.ring_of_integers().basis()
            L = MS([[1, b1.complex_embedding(prec)], [0, 1]])
            M = MS([[1, b2.complex_embedding(prec)], [0, 1]])
            S = MS([[0, -1], [1, 0]])
            self._init_prec = prec
            self._gens = [L, M, S, L ** -1, M ** -1]
            self._gens.append(MS([[0, -1], [1, 0]]))
            A = S
            B = S
            self._named_gens = {
                'L': L, 'M': M,
                'l': L ** -1, 'm': M ** -1,
                'A': A, 'a': A**-1, 'B': B, 'b': B**-1
            }
            if not self._name_string:
                self._name_string = f'Bianchi Group: Q(sqrt({x}))'
            if not self._latex_string:
                self._latex_string = rf'Bianchi Group: $\mathbb{{Q}}(\sqrt{{{x}}})$'
        elif isinstance(x, (list, tuple)):
            x = dict(x)
            if not all(isinstance(g, Matrix) for g in x.values()):
                raise NotImplementedError(f"Can not make a Kleinian group from {x}")
            self._gens = list(x.values())
            self._named_gens = x
            # In this case it is something like ComplexField
            base_ring = self._gens[0].base_ring()
        else:
            raise NotImplementedError
        super().__init__(degree=Integer(2), base_ring=base_ring,
                         special=True,
                         sage_name=self._name_string,
                         latex_string=self._latex_string,
                         category=Groups().Infinite(),
                         invariant_form=None)

    def _cache_key(self):
        return json.dumps(self.to_json())

    def name(self):
        return self._name_string

    def to_json(self, **kwargs: P.kwargs) -> dict:
        return {
            '_gens': [matrix_to_json(g) for g in self._gens],
            '_covering_generators_words': self._covering_generators_words,
            '_named_gens': {k: matrix_to_json(g) for k, g in self._named_gens.items()},
            '_covering_generators': {k: matrix_to_json(g) for k, g in
                                     self._covering_generators.items()},
            '_name_string': self._name_string,
            '_latex_string': self._latex_string,
            '_manifold': self._manifold,
            '_translation_lattice': matrix_to_json(self.translation_lattice().basis_matrix()),
            'type': 'KleinianGroup',
        }

    def __eq__(self, other: 'KleinianGroup_class') -> bool:
        return self.to_json() == other.to_json()
    @classmethod
    def from_json(cls, data: dict | str) -> 'KleinianGroup':
        r"""
        Create a KleinianGroup object from json

        INPUT:

        - ``data`` - json data

        EXAMPLES::

            sage: from knot_maass.hyperbolic_space.kleinian_group import KleinianGroup_class
            sage: import json
            sage: K = KleinianGroup_class(-4)
            sage: K.from_json(json.dumps(K.to_json())) == K
            True
        """
        if isinstance(data, dict):
            data = json.dumps(data)
        return KleinianGroup_class(data=data)

    def __repr__(self):
        return self._name_string or f'Kleinian Group ({self.base_ring()})'

    @cached_method
    def generators(self, include_parabolic: bool = True, prec: int = 0):
        """
        Return a (in general not minimal) list of generators.
        :return:
        """
        if prec:
            CF = ComplexField(prec=prec)
        return [g.change_ring(CF) if prec else g for g in self._gens
                if include_parabolic or g.trace() ** 2 != 4]

    @cached_method
    def named_generators(self, include_parabolic: bool = True,
                         only_parabolic: bool = False,
                         prec: int = 0):
        """
        Return a (in general not minimal) dictionary of names and generators.
        Note: This is useful for word problems.

        INPUT:

        - ``include_parabolic`` -- include parabolic generators
        - ``only_parabolic`` -- only parabolic generators
        - ``prec`` -- precision

        EXAMPLES::


        """
        if only_parabolic and not include_parabolic:
            return {}
        if prec:
            CF = ComplexField(prec=prec)
        return {name: g.change_ring(CF) if prec else g for name, g in self._named_gens.items()
                if (include_parabolic or g.trace() ** 2 != 4) and
                   (not only_parabolic or g.trace() ** 2 == 4)
                }

    def covering_generators(self, include_parabolic=False, prec=53):
        """
        A set of generators such that the corresponding fundamental domain has oly one vertex.


        """
        if self._covering_generators:
            return self._covering_generators
        gens = self.named_generators(prec=prec)
        parallelogram = self.translation_fundamental_domain()
        if not self._covering_generators_words:
            word_list = []
            # Try to get an initial list of words from database
            try:
                from knot_maass.database.models import Word
                if self._manifold and Word.objects(label=self._manifold):
                    w = (Word.objects(label=self._manifold).
                         order_by('-reduced_cover,-reduced_to_fd').first())
                    word_list = w.words
            except ImportError:
                pass
            self._covering_generators_words = find_covering_generators(parallelogram, gens,
                                                                       covering_list=word_list)
            self._covering_generators_words = [expand_word(x[0]) for x in self._covering_generators_words]

        gens = {x: self.word_to_element(x, prec=prec) for x in self._covering_generators_words}
        if include_parabolic:
            gens.update(self.named_generators(only_parabolic=True, include_parabolic=True,
                                              prec=prec))
        self._covering_generators = gens
        return gens

    def generators_parabolic(self):
        """
        Return a dictionary of
        :return:
        """
        return self.generators(include_parabolic=True)

    def side_pairing_gens(self):
        if not self.manifold():
            raise NotImplementedError("Only implemented for Kleinian groups defined by manifolds")
        if not self._side_pairing_gens:
            self._side_pairing_gens = find_side_pairing_gens(self.manifold())
        return self._side_pairing_gens

    @cached_method
    def word_to_element(self, word, prec=53):
        return word_to_element(word, self.named_generators(prec=prec))

    def translation_lattice(self, prec=53):
        r"""
        Pull back a point in the upper half-space to an element of the fundamental
        domain.
        """
        if self._translation_lattice and prec == self._translation_lattice.base_ring().prec():
            return self._translation_lattice
        if isinstance(self.base_ring(), NumberField) and self.base_ring().discriminant() == -4:
            basis = [[1, 0], [0, 1]]
        else:
            L = self.named_generators()['L']
            M = self.named_generators()['M']
            tL = L[0][0] * L[0][1]
            tM = M[0][0] * M[0][1]
            basis = [[real(tL.n(prec)), imag(tL.n(prec))],
                     [real(tM.n(prec)), imag(tM.n(prec))]]
        self._translation_lattice = (RealField(prec)**2).span_of_basis(basis)
        if self._translation_lattice.rank() != 2:
            raise ArithmeticError(f"Translation lattice has rank"
                                  f" {self._translation_lattice.rank()}. "
                                  f"Expected 2.")
        return self._translation_lattice

    def translation_fundamental_domain(self, base: Vector = None):
        r"""
        A choice of fundamental domain for the translations.

        INPUT:

        - ``base`` -- base point of the fundamental domain (default = None)
                          default position is centered at (0,0)

        EXAMPLES::

            sage: from knot_maass.hyperbolic_space.kleinian_group import KleinianGroup
            sage: from knot_maass.hyperbolic_space.types import Rectangle
            sage: K = KleinianGroup(-4)
            sage: K.translation_fundamental_domain()
            Parallelogram(base=(-0.500000000000000, -0.500000000000000), v1=(1.00000000000000, ...
            sage: K.translation_fundamental_domain(base=vector([0, 0]))
            Parallelogram(base=(0, 0), v1=(1.00000000000000, 0.000000000000000),...
        """
        v1, v2 = self.translation_lattice().basis()
        if base is None:
            base = - v1 / 2 - v2 / 2
        return Parallelogram(base=base, v1=v1, v2=v2)

    def _init_from_manifold(self, M: 'Manifold'):
        r"""
        TODO: Implement this

        INPUT:
        - ``M`` -- Snappy manifold
        """
        raise NotImplementedError

    def pullback(self, z: UpperHalfSpaceElement__class, check=True, use_exact=False):
        r"""
        Pull back a point in the upper half-space to an element of the fundamental
        domain.
        """
        if self.base_ring().discriminant() == -4:
            w, map = self._pullback_gaussian_integers(z)
        else:
            w, map, word = self._pullback_general(z)
        # Use exact matrices (slower)
        if use_exact:
            map = word_to_element(word, self._named_gens)
            z_mapped = z.action(map, check=check)
        else:
            z_mapped = w
        # if check and (z_mapped - w).norm() > 16 * w[0].base_ring().epsilon():
        #     raise ArithmeticError(f"Point `{z}` is not pulled back correctly: "
        #                           f"{z.action(map)} != {w}"
        #                           f" norm={(z.action(map) - w).norm()} > "
        #                           f"{16 * w[0].base_ring().epsilon()}")
        return z_mapped, map

    @cached_method
    def T(self, a=1):
        """
        Return the element T^a = ( 1 & a // 0 & 1 ) of self.

        INPUT:

        - ``a`` -- integer in number field (default=1)

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: H.T()
            [1 1]
            [0 1]
            sage: u0,u1=H.base_ring().number_field().unit_group().gens()
            sage: H.T(u0)
            [ 1 -1]
            [ 0  1]
            sage: H.T(u0*u1)
            [          1 1/2*a - 1/2]
            [          0           1]
        """
        return self([1, a, 0, 1])

    def _is_reduced_translation(self, z: UpperHalfSpaceElement__class, eps: Real_t = 0):
        lattice_coordinates = self.translation_lattice().coordinates(list(z.z()))
        return all(-1 / 2 - eps <= v <= 1 / 2 + eps for v in lattice_coordinates)

    def reduce_by_translations(self, z: UpperHalfSpaceElement__class):
        r"""
        Reduce z with respect to the lattice generated by the cusp translations.

        INPUT:

        - ``z`` -- UpperHalfSpaceElement

        OUTPUT:

        - ``(z, (t1, t2))`` - a reduced point and a matrix

        EXAMPLES::

            sage: from knot_maass.all import KleinianGroup, UpperHalfSpaceElement
            sage: K=KleinianGroup(-4)
            sage: z = UpperHalfSpaceElement([1, 2, 3])
            sage: w, A, s = K.reduce_by_translations(z)
            sage: w
            0.000000000000000 + 0.000000000000000i + 3.00000000000000j
            sage: A
                [ 1.00000000000000 -1.00000000000000 - 2.00000000000000*I]
                [ 0.000000000000000                       1.00000000000000]
            sage: s
            'lmm'
        """
        lattice_coordinates = self.translation_lattice().coordinates(list(z.z()))
        t1 = - floor(lattice_coordinates[0] + 1 / 2)
        t2 = - floor(lattice_coordinates[1] + 1 / 2)
        basis_matrix = self.translation_lattice().basis_matrix()
        translation = vector([t1, t2]) * basis_matrix
        translation = ComplexField(basis_matrix.base_ring().prec())(list(translation))
        word = ''
        word += 'L' if t1 > 0 else 'l' * abs(t1)
        word += 'M' if t2 > 0 else 'm' * abs(t2)
        return z.translate(translation), matrix([[1, translation], [0, 1]]), word

    def _is_reduced_reflections(self, z: UpperHalfSpaceElement__class, eps: Real_t = 1e-15):
        for g in self.covering_generators(include_parabolic=False, prec=53).values():
            circle = matrix_to_circle(g)
            dist = (z - circle.center).norm()
            if dist < circle.radius - eps:
                return False
        return True

    def reduce_by_reflections(self, z: UpperHalfSpaceElement__class, eps: Real_t = 0):
        eps = eps or 2**(8 - z[0].prec())
        max_height = z.y()
        max_w = z
        max_g = None
        max_name = None
        for name, g in self.covering_generators(include_parabolic=False, prec=z[0].prec()).items():
            c = matrix_to_circle(g)
            if (z - c.center).norm() < c.radius - eps:
                w = z.action(g, check=False)
                if w.y() > max_height:
                    max_height = w.y()
                    max_w = w
                    max_g = g
                    max_name = name
        if max_g:
            return max_w, max_g, max_name
        return z, matrix([[1, 0], [0, 1]]), ""

    def _is_reduced(self, z: UpperHalfSpaceElement__class, eps: Real_t = 0):
        return self._is_reduced_translation(z, eps=eps) and self._is_reduced_reflections(z, eps=eps)


    def _pullback_gaussian_integers(self, z: UpperHalfSpaceElement__class) -> (
            tuple)[UpperHalfSpaceElement__class, matrix]:
        r"""
        Special case of Gaussian integers

        TODO: Implement this and write examples.
        NOTE: Could be made faster by e.g. calling a Cython function.

        EXAMPLES:

            sage: from knot_maass.all import KleinianGroup, UpperHalfSpaceElement
            sage: K=KleinianGroup(-4)
            sage: K._covering_generators_words = ['A', 'B']
            sage: z = UpperHalfSpaceElement([1.0, 2, 3.5])
            sage: w, A = K._pullback_gaussian_integers(z)
            sage: w
            0.000000000000000 + 0.000000000000000i + 3.50000000000000j
            sage: A
                [ 1.00000000000000 -1.00000000000000 - 2.00000000000000*I]
                [ 0.000000000000000                       1.00000000000000]
            sage: z.action(A) == w
            True
            sage: z = UpperHalfSpaceElement([1.0, 2.0, 0.5])
            sage: w, A = K._pullback_gaussian_integers(z)
            sage: w
            0.000000000000000 + 0.000000000000000i + 2.00000000000000j
            sage: A
                [ 0.000000000000000 -1.00000000000000]
                [ 1.00000000000000 -1.00000000000000 - 2.00000000000000*I]
            sage: z.action(A) == w
            True


        """
        if isinstance(z[0], Integer_t):
            prec = 53
        elif hasattr(z[0], 'parent'):
            prec = z[0].parent().prec()
        else:
            raise ValueError(f"Can not find precision for z={z}")
        CF = ComplexField(prec=prec)
        mat = matrix(CF, [[1, 0], [0, 1]])
        inversion = matrix(CF, [[0, -1], [1, 0]])
        eps = 8 * CF.epsilon()
        while not self._is_reduced(z, eps=eps):
            if not self._is_reduced_translation(z, eps=eps):
                z, A, _ = self.reduce_by_translations(z)
                mat = A * mat
            if z.norm() < 1 - eps:
                z = z.action(inversion, check=False)
                mat = inversion * mat
        if not self._is_reduced(z, eps=eps):
            raise ArithmeticError(f"Could not reduce point: z={z}")
        return z, mat

    def _pullback_general(self, z: UpperHalfSpaceElement__class, max_iterations: int = 1000) -> (
            tuple)[UpperHalfSpaceElement__class, matrix]:
        r"""
        Special case of Gaussian integers

        TODO: Implement this and write examples.
        NOTE: Could be made faster by e.g. calling a Cython function.

        EXAMPLES:

            sage: from knot_maass.all import KleinianGroup, UpperHalfSpaceElement
            sage: K=KleinianGroup(-4)
            sage: z = UpperHalfSpaceElement([1.0, 2, 3.5])
            sage: w, A = K._pullback_gaussian_integers(z)
            sage: w
            0.000000000000000 + 0.000000000000000i + 3.50000000000000j
            sage: A
                [ 1.00000000000000 -1.00000000000000 - 2.00000000000000*I]
                [ 0.000000000000000                       1.00000000000000]
            sage: z.action(A) == w
            True
            sage: z = UpperHalfSpaceElement([1.0, 2.0, 0.5])
            sage: w, A = K._pullback_gaussian_integers(z)
            sage: w
            0.000000000000000 + 0.000000000000000i + 2.00000000000000j
            sage: A
                [ 0.000000000000000 -1.00000000000000]
                [ 1.00000000000000 -1.00000000000000 - 2.00000000000000*I]
            sage: z.action(A) == w
            True


        """
        if isinstance(z[0], Integer_t):
            prec = 53
        elif hasattr(z[0], 'parent'):
            prec = z[0].parent().prec()
        else:
            raise ValueError(f"Can not find precision for z={z}")
        CF = ComplexField(prec=prec)
        mat = matrix(CF, [[1, 0], [0, 1]])
        inversion = matrix(CF, [[0, -1], [1, 0]])
        eps = 8 * CF.epsilon()
        n = 0
        word = ""
        while not self._is_reduced(z, eps=eps):
            z, A, s = self.reduce_by_translations(z)
            mat = A * mat
            word = s + word
            z, g, s = self.reduce_by_reflections(z, eps=eps)
            mat = g * mat
            word = s + word
            if n > max_iterations:
                break
            n += 1
        if n >= max_iterations:
            raise ArithmeticError(f"Could not reduce point: z={z}")
        return z, mat, word

    def cusp_normaliser(self, cusp: NFCusp | ComplexNumber):
        if isinstance(cusp, NFCusp) and cusp.number_field() != self.base_ring():
            raise ValueError("Cusp must be defined over the base ring of the Kleinian group.")
        if isinstance(cusp, NFCusp):
            return cusp.ABmatrix()
        return matrix(self.base_ring(), [[1, 0], [cusp ** -1, 1]])

    @cached_method
    def dual_translation_lattice_vectors(self, M: Integer_t,
                                         return_indices: bool = False) -> list | tuple[list, list]:
        B = self.translation_lattice().basis_matrix()
        dual_lattice_basis = tuple((B ** -1).transpose())
        return get_lattice_values(dual_lattice_basis, - M, M + 1, return_indices=return_indices)

    def check_translation_cell_coverage(self) -> bool:
        b1, b2 = self.translation_lattice().basis()
        if b1.imag() != 0 or b2.real() != 0:
            raise ArithmeticError("Translation basis should be parallel with standard basis."
                                  f"Basis: {b1}, {b2}")
        center = (0.0, 0.0)
        sides = (b1.real(), b2.imag())
        gens = self.generators(include_parabolic=False)
        for n in range(1, 10, 2):
            rectangles = split_rectangle(center, sides, n)
            for rect in rectangles:
                for g in gens:
                    c, r = matrix_to_circle(g)
                    if not rectangle_in_circle(rect[0], rect[1], c, r):
                        print("Rectangle", rect, "not in circle", c, r)

    def manifold(self) -> Manifold:
        if not self._manifold:
            raise ValueError("This group is not initialized from a manifold.")
        return Manifold(self._manifold)

class KleinianGroupElement__class(Element):

    _is_parabolic = None
    _is_hyperbolic = None
    _is_elliptic = None
    _is_loxodromic = None

    def __init__(self, x, parent, *args: P.args, **kwargs: P.kwargs):
        if isinstance(x, Matrix):
            self._matrix = x
        self._parent = parent
        # Account for loss of precision in matrix operations.
        self._epsilon = 2 * x.base_ring().epsilon()

    def _matrix_(self):
        return self._matrix

    def trace(self):
        return self._matrix.trace()

    def _trace_is_real(self):
        return abs(self.trace().imag()) < self._epsilon

    def is_parabolic(self):
        if self._is_parabolic is None:
            self._is_parabolic = abs(self.trace()**2 - 4) < self._epsilon
        return self._is_parabolic

    def is_hyperbolic(self):
        if self._is_hyperbolic is None:
            self._is_hyperbolic = self._trace_is_real() and \
                                   (self.trace()**2).real() - 4 > self._epsilon
        return self._is_hyperbolic

    def is_elliptic(self):
        if self._is_hyperbolic is None:
            self._is_hyperbolic = self._trace_is_real() and \
                                  (self.trace()**2).real() - 4 < -self._epsilon
        return self._is_hyperbolic

    def is_loxodromic(self):
        if self._is_loxodromic is None:
            self._is_loxodromic = not self._trace_is_real()
        return self._is_loxodromic

    def fixed_points(self):
        if self.is_parabolic():
            a, b, c, d = list(self._matrix)

            if c == 0 or abs(c) < self._epsilon:
                if isinstance(self.base_ring(), NumberField):
                    return NFCusp(self.base_ring(), [1, 0])
                else:
                    return Infinity
            if isinstance(self.base_ring(), NumberField):
                return NFCusp(self.base_ring(), [(a-d), 2*c])
            else:
                return (a - d) / (c * 2)
        raise NotImplementedError("Fixed point only implemented for parabolic elements")


def KleinianGroup(*args: P.args, **kwargs: P.kwargs):
    """
    Construct a Kleinian group from given input.

    INPUT:
    - ``args`` -- arguments for the constructor
    - ``kwargs`` -- keyword arguments for the constructor

    EXAMPLES::

        sage: from knot_maass.hyperbolic_space.kleinian_group import KleinianGroup
        sage: G = KleinianGroup(4)
    """
    if len(args) == 0:
        return KleinianGroup_class(**kwargs)
    if isinstance(args[0], Integer):
        return KleinianGroup_class(*args, **kwargs)
    if isinstance(args[0], (dict, str)):
        return KleinianGroup_class.from_json(args[0])
    try:
        from snappy import Manifold
    except ImportError:
        Manifold = None
    if isinstance(args[0], Manifold):
        return KleinianGroup__from_manifold(args[0], **kwargs)
    return KleinianGroup_class(*args, **kwargs)


def KleinianGroup__from_manifold(manifold: 'Manifold', **kwargs: P.kwargs):
    """
    Create a Kleinian group from a manifold.

    INPUT:
    - ``manifold`` -- the manifold

    EXAMPLES::

        sage: from knot_maass.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from snappy import Manifold
        sage: M = Manifold('4_1')
        sage: G = KleinianGroup(M)
    """
    # TF = manifold.trace_field_gens()
    # K, g, _ = TF.find_field(prec=kwargs.get('prec', 300), degree=kwargs.get('degree', 20),
    #                         optimize=True)
    # Get translations
    #T = manifold.cusp_neighborhood().all_translations()
    HME = manifold.holonomy_matrix_entries()
    try:
        K, g, entries = HME.find_field(prec=300, degree=kwargs.get('degree', 20), optimize=True)
    except TypeError:
        entries = [x(kwargs.get('prec', 100)) for x in HME.list()]
        K = ComplexField(kwargs.get('prec', 100))
        print("eps=", K.epsilon())
        pass
    G = manifold.fundamental_group()
    nf_gens = {
        'a': matrix(K, [[entries[0], entries[1]], [entries[2], entries[3]]]),
        'b': matrix(K, [[entries[4], entries[5]], [entries[6], entries[7]]])
    }
    nf_gens['A'] = nf_gens['a'].inverse()
    nf_gens['B'] = nf_gens['b'].inverse()
    if len(HME.list()) == 12:
        nf_gens['c'] = matrix(K, [[entries[8], entries[9]], [entries[10], entries[11]]])
        nf_gens['C'] = nf_gens['c'].inverse()
    # Parabolics
    L = word_to_element(G.longitude(), nf_gens)
    M = word_to_element(G.meridian(), nf_gens)
    if abs(L.trace()**2 - 4) > 1e-10 or abs(M.trace()**2 - 4) > 1e-10:
        raise ArithmeticError("Longitude and meridian are not parabolic!")
    if isinstance(K, NumberField):
        cuspL = NFCusp(K, [L[0][0] - L[1][1], 2 * L[1][0]])
        cuspM = NFCusp(K, [M[0][0] - M[1][1], 2 * M[1][0]])
        if cuspL != cuspM:
            raise ArithmeticError("Longitude and meridian do not have the same fixed point!")
        normaliser = matrix(2, 2, cuspL.ABmatrix())
    else:
        cuspL = [L[0][0] - L[1][1], 2 * L[1][0]]
        cuspM = [M[0][0] - M[1][1], 2 * M[1][0]]
        if abs(cuspL[0] * cuspM[1] - cuspL[1] * cuspM[0]) > 1e-16:
            raise ArithmeticError("Longitude and meridian do not have the same fixed point!")
        normaliser = matrix(2, 2, [[cuspL[0], 0], [cuspL[1], 1/cuspL[0]]])
    normalised_gens_nf = {k: normaliser.inverse() * g * normaliser for (k, g) in nf_gens.items()}
    normalised_gens_nf['L'] = normaliser.inverse() * L * normaliser
    normalised_gens_nf['l'] = normalised_gens_nf['L'].inverse()
    normalised_gens_nf['M'] = normaliser.inverse() * M * normaliser
    normalised_gens_nf['m'] = normalised_gens_nf['M'].inverse()
    # Try to find side-pairing generators
    side_pairing_gens = find_side_pairing_gens(manifold)
    return KleinianGroup_class(tuple(normalised_gens_nf.items()),
                               manifold=manifold.name(),
                               side_pairing_gens=side_pairing_gens)


# def KleinianGroup__from_manifold2(manifold: 'Manifold', **kwargs: P.kwargs):
#     """
#     Create a Kleininan group from a manifold.
#     :param args:
#     :param kwargs:
#     :return:
#     """
#     # First extract the raw side-pairing transformations
#     eps = 2 ** (8 - 53)
#     F = FundamentalPolyhedronEngine.from_manifold_and_shapes(manifold,
#                                                              manifold.tetrahedra_shapes('rect'),
#                                                              normalize_matrices=True)
#     parabolic_elements = []
#     other_elements = []
#     vertices = []
#
#     for k, A in F.mcomplex.GeneratorMatrices.items():
#         if k == 0:
#             continue
#         if abs(A.trace()**2 - 4) < eps:
#             parabolic_elements.append(A)
#             a, b, c, d = A.list()
#             vertices.append((a - d) / (2 * c))
#         else:
#             other_elements.append(A)
#     # other_elements_inverse = [B**-1 for B in other_elements]
#     # parabolic_elements_inverse = [B**-1 for B in parabolic_elements]
#     if len(parabolic_elements) < 2:
#         raise ValueError("Manifold must have a parabolic of rank 2!")
#     # Now we need to sort out the cusps and map them to the first cusp
#     identity = SL(2, CC).one()
#     mappings = [identity]
#     cusp1 = UpperHalfSpaceElement([vertices[0][0], vertices[0][1], 0])
#     for c in vertices[1:]:
#         cusp = UpperHalfSpaceElement([c[0], c[1], 0])
#         for B in  [identity] + other_elements: #+ other_elements_inverse:
#             w = cusp.action(B)
#             if (w - cusp1).norm() < eps:
#                 print("mapping", cusp, "->",cusp1, B)
#                 mappings.append(B)
#                 break
#     mapped_parabolics = [m * parabolic_elements[i] * m.inverse() for i, m in enumerate(mappings)]
#
#     print("parabolics:")
#     for p in parabolic_elements:
#         print(p)
#         print("fixed=", (p[0, 0] - p[1, 1]) / (2 * p[1, 0]))
#         print("-")
#     print("mappings=")
#     for m in mappings:
#         print(m)
#     print("mapped_parabolics=")
#     for p in mapped_parabolics:
#         print(p)
#         print("fixed=", (p[0,0]-p[1,1])/(2*p[1,0]))
#         print("-")
#     if cusp1.z() == 0:
#         normaliser = matrix(CC, [[0, -1], [1, 0]])
#     else:
#         normaliser = matrix(CC, [[1, 0], [1/cusp1.z(), 1]])
#     T1 = normaliser.inverse() * parabolic_elements[0] * normaliser
#     print("T1=", T1)
#     L = T1[0, 1]
#     # scaling = matrix([[L.sqrt(), 0], [0, 1/L.sqrt()]])
#     # print("scaling=", scaling)
#     # print("ST1=", scaling.inverse() * T1 * scaling)
#     # print("ST1=", scaling * T1 * scaling.inverse())
#     # normaliser = normaliser * scaling
#     generators = [normaliser.inverse() * g * normaliser for g in mapped_parabolics]
#     generators += [normaliser.inverse() * g * normaliser for g in other_elements]
#     return generators
#
#     # TF = manifold.trace_field_gens()
#     # K, g, _ = TF.find_field(prec=kwargs.get('prec', 100), degree=kwargs.get('degree', 10),
#     #                         optimize=True)
#     # # Get translations
#     # #T = manifold.cusp_neighborhood().all_translations()
#     # HME = manifold.holonomy_matrix_entries()
#     # K, g, entries = HME.find_field(prec=100, degree=10, optimize=True)
#     # A = matrix(K, [[entries[0], entries[1]], [entries[2], entries[3]]])
#     # B = matrix(K, [[entries[4], entries[5]], [entries[6], entries[7]]])
#     # P = A**-1 * B**-1
#     # has_parabolic = None
#     # for p_test in [A, B, P]:
#     #     if abs(p_test.trace()**2 - 4) < 1e-10:
#     #         has_parabolic = p_test
#     # if not has_parabolic:
#     #     raise ValueError("Manifold must have a parabolic holonomy matrix")
#     # print(has_parabolic, type(has_parabolic))
#     # a, b, c, d = has_parabolic.list()
#     # cusp = NFCusp(K, [a-d, 2*c])
#     # print("hsa_parabolic=", has_parabolic.n(53))
#     # print("cusp=",ComplexField(53)((a-d)/(2*c)))
#     # a1, b1, c1, d1 = cusp.ABmatrix()
#     # normaliser = matrix(K, [[a1, b1], [c1, d1]])
#     # print("normaliser=", normaliser.n(53))
#     # gens_original = [has_parabolic]
#     # if has_parabolic == A:
#     #     gens_original.append(B)
#     # else:
#     #     gens_original.append(A)
#     # generators = [normaliser.inverse() * g * normaliser for g in gens_original]
#     # return generators, normaliser

def KleinianGroupElement(*args: P.args, **kwargs: P.kwargs):
    return KleinianGroupElement__class(*args, **kwargs)

