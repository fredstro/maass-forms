import logging
import string
from typing import TYPE_CHECKING, ParamSpec, Union

from sage.combinat.words.words import FiniteWords
from sage.groups.free_group import FreeGroup
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.real_mpfr import RealNumber

from maass_forms_klein.hyperbolic_space.geometry_utils import (
    circle_approx_equal,
    circle_is_contained,
    point_in_circle,
    reduce_cover,
    remove_contained_circles,
    remove_duplicate_circles,
    remove_non_intersecting_circles,
)
from maass_forms_klein.hyperbolic_space.parallelogram import (
    parallelogram_covered_by_circles,
    parallelogram_intersect_circle,
    point_on_boundary_of_parallelogram,
    reduce_in_parallelogram,
)
from maass_forms_klein.hyperbolic_space.types import Circle, Parallelogram, Rectangle
from maass_forms_klein.hyperbolic_space.word_utils import (
    change_letters_in_word,
    is_equivalent_mod_parabolic,
    reduce_word,
    translation_tuple_to_word,
    word_list_sort_key,
    word_to_circle,
)

if TYPE_CHECKING:
    from snappy import Manifold

Real_t = Union[RealNumber, Integer, Rational, int, float]
Integer_t = Union[Integer, int]

P = ParamSpec("P")
log = logging.getLogger(__name__)


def get_lattice_values(
    lattice_basis: tuple[tuple[Real_t, Real_t], tuple[Real_t, Real_t]],
    Q0: Integer_t,
    Q1: Integer_t = None,
    shift: bool = False,
    return_indices: bool = False,
) -> list | tuple[list, list]:
    r"""
    Returns a list of lattice values v = [i * v0 + j * v1 for i in [-Q0, Q1-1], j in [-Q0, Q1-1]]

    INPUT:
        ``lattice_basis`` -- (list of real numbers) basis of the scaled lattice
        ``Q0`` -- integer
        ``Q1`` -- integer (default = None => Q1=Q0)

    EXAMPLES:
        sage: from maass_forms_klein.hyperbolic_space.utils import get_lattice_values
        sage: get_lattice_values(((1,0), (0,1)), 0, 2)
        [(0, 0), (0, 1), (1, 0), (1, 1)]
        sage: get_lattice_values(((1,0), (0,1)), 1)
        [(0, 0), (-1, 0), (0, -1), (-1, -1)]
    """
    if Q1 is None:
        Q1 = abs(Q0)
        Q0 = -abs(Q0)
    if Q0 > Q1:
        raise ValueError(f"Q0=`{Q0}` must be smaller than Q1=`{Q1}`")
    if shift:
        shiftv = 1 / 2
    else:
        shiftv = 0
    shiftv = 0
    cartesian_product = [(v1, v2) for v1 in range(Q0, Q1) for v2 in range(Q0, Q1)]
    cartesian_product.sort(key=lambda x: abs(complex(x[0], x[1])))
    lattice_values = [
        vector(lattice_basis[0]) * (v[0] - shiftv) + vector(lattice_basis[1]) * (v[1] - shiftv)
        for v in cartesian_product
    ]
    if return_indices:
        return lattice_values, cartesian_product
    else:
        return lattice_values


def filter_list_mod_parabolics(list_of_words, gens):
    """
    Filter a list of words representing group elements modulo parabolic elements (on the left),
    i.e. a ~ b if a = t * b for a parabolic t.

    INPUT:

    - ``list_of_words`` -

    INPUT:

        sage: from maass_forms_klein.hyperbolic_space.utils import (  # doctest: +ELLIPSIS
        ....:     filter_list_mod_parabolics)

        ...
        sage: filter_list_mod_parabolics(['a', 'b'], {'a': matrix([[0, 1], [1, 0]]),
        ....:                                            'b': matrix([[0, 1], [1, 0]])})
        ['a']

    """
    new_list = []
    for x in list_of_words:
        equivalent_element_exists = False
        for y in new_list:
            if is_equivalent_mod_parabolic(x, y, gens):
                equivalent_element_exists = True
                break
        if not equivalent_element_exists:
            new_list.append(x)
    return new_list


def find_covering_generators(
    rect: Rectangle | Parallelogram,
    named_gens: dict,
    covering_list: list | None = None,
    check: bool = True,
    itermax: int = 4,
    **kwargs: P.kwargs,
) -> list[tuple[str, Circle]] | list[str]:
    """
    Find a set of generators that cover the rectangle.

    Note 1: We assume that the generators are named in the form 'a', 'A', 'b', 'B', 'ab', 'ba',
        'AB', 'BA', 'bA', 'aB', 'Ba', 'Ab'.
    Note: WE do this by translating up/down left/right with parabolic elements.

    EXAMPLES:
        sage: from maass_forms_klein.hyperbolic_space.utils import find_covering_generators
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram, Rectangle
        sage: s = matrix([[0, -1], [1, 0]])
        sage: t = matrix([[1, 1], [0, 1]])
        sage: l = matrix([[1, I], [0, 1]])
        sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: find_covering_generators(rect, {'a': s,'b':s,'B':s,'A':s})
        Traceback (most recent call last):
        ...
        ArithmeticError: Could not find a covering list.
        sage: gens = {'a': s,'b':s,'B':s,'A':s, 'm':t, 'M': t**-1, 'l':l, 'L':l**-1}
        sage: find_covering_generators(rect, gens)
        [('a',
          Circle(center=(-0.000000000000000, 0.000000000000000), radius=1.00000000000000)),
         ('aL',
          Circle(center=(0.000000000000000, 1.00000000000000), radius=1.00000000000000)),
         ('aM',
          Circle(center=(1.00000000000000, 0.000000000000000), radius=1.00000000000000)),
         ('aLM',
          Circle(center=(1.00000000000000, 1.00000000000000), radius=1.00000000000000))]
        sage: rect = Rectangle(base=vector((-0.5, -0.5)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: find_covering_generators(rect, gens)
         [('a',
          Circle(center=(-0.000000000000000, 0.000000000000000), radius=1.00000000000000))]
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import (
        ....:     KleinianGroup__from_manifold)

        sage: from snappy import Manifold
        sage: M = Manifold("4_1")
        sage: G = KleinianGroup__from_manifold(M)
        sage: gens = G.named_generators(prec=53)
        sage: v1 = vector((3.46410161513775, 0))
        sage: rect = Parallelogram(base=vector((0, 0)), v1=v1, v2=vector((0, 1)))
        sage: find_covering_generators(rect, gens)
         [('AB',
              Circle(center=(1.11022302462516e-16, 1.92296268638356e-16), radius=1.00000000000000)),
             ('Bl',
              Circle(center=(1.73205080756888, 3.33066907387547e-16), radius=1.00000000000000)),
             ('bl',
              Circle(center=(3.46410161513775, 2.22044604925031e-16), radius=1.00000000000000)),
             ('abM',
              Circle(center=(0.866025403784438, 0.500000000000001), radius=1.00000000000000)),
             ('BAb',
              Circle(center=(-9.99200722162641e-16, 1.00000000000000), radius=1.00000000000000)),
             ('BAl',
              Circle(center=(2.59807621135332, 0.499999999999999), radius=1.00000000000000)),
             ('bal',
              Circle(center=(1.73205080756888, 1.00000000000000), radius=1.00000000000000)),
             ('bml',
              Circle(center=(3.46410161513775, 1.00000000000000), radius=1.00000000000000))]
        sage: rect =  G.translation_fundamental_domain()
        sage: find_covering_generators(rect, gens, return_words=True)
        ['ABab', 'BAB', 'BAbM', 'BAba', 'bMBl', 'baMl', 'bmaa']
    """
    verbose = kwargs.get("verbose", False)
    if not isinstance(named_gens, dict):
        raise ValueError("named_gens must be a dictionary.")
    if not isinstance(rect, (Rectangle, Parallelogram)):
        raise ValueError(f"rect must be a Rectangle or Parallelogram: {type(rect)}")
    if not covering_list:
        # Find some starting words
        covering_list = [(x, word_to_circle(x, named_gens)) for x in named_gens.keys()]
    test = False
    n = 0

    # Check if the reduced list works:
    while not test and n < itermax:
        # First do a preliminary reduction for circles that are covered by others.from
        if verbose > 0:
            log.debug("Check list of length: ", len(covering_list))
        if verbose > 0 and len(covering_list) < 1500:
            log.debug(
                "new check list= %s",
                [x[0] if isinstance(x, tuple) else x for x in covering_list],
            )
        elif verbose > 0:
            log.debug("trunc check list=", list(covering_list)[0:10])
        covering_list = remove_duplicate_circles(covering_list, named_gens)
        if verbose > 0 and len(covering_list) < 1500:
            log.debug("after remove dup, length=", len(covering_list))
            log.debug("new dup list=", [x[0] for x in covering_list])
        elif verbose > 0:
            log.debug("trunc dup list=", list(covering_list)[0:10])
        if verbose > 0:
            log.debug("Removed duplicated length: ", len(covering_list))
        try:
            covering_list_contained = remove_contained_circles(covering_list, named_gens)
            if verbose > 0:
                log.debug("Removed contained circles: ", len(covering_list))
            if verbose > 0:
                log.debug("list=", [x[0] for x in covering_list])
                log.debug("rec=", rect)

            covering_list_intersect = remove_non_intersecting_circles(covering_list_contained, rect)
            test = parallelogram_covered_by_circles(
                rect,
                [x[1] for x in covering_list_intersect],
                scaling_factor=0.999,
                n_max=kwargs.get("n_max", 20),
                verbose=verbose - 1 if verbose > 0 else 0,
            )
            if verbose > 0:
                log.debug("Test: ", test)
            if verbose > 0 and len(covering_list_intersect) < 200:
                log.debug("new list=", [x[0] for x in covering_list_intersect])
            if test:
                covering_list = covering_list_intersect
        except ArithmeticError as e:
            log.info(str(e))
            pass
        # Update the list
        if test:
            break
        logging.debug("Increase word length.")
        new_gens = []
        for g in named_gens:
            ginv = g.swapcase()
            new_gens += [x[0] + g for x in covering_list if x[0][-1] != ginv]
        new_gens = [(x, word_to_circle(x, named_gens)) for x in new_gens]
        covering_list += new_gens
        if verbose > 0:
            log.debug("len new gens=", len(new_gens))
            log.debug("len new list=", len(covering_list))
        n += 1
    if not test:
        raise ArithmeticError("Could not find a covering list.")
    return_words = kwargs.pop("return_words", False)
    res = reduce_cover(rect, covering_list, named_gens, **kwargs)
    if return_words:
        return sorted([x[0] for x in res])
    return res


def find_covering_generators2(
    rect: Rectangle | Parallelogram, named_gens: dict, check: bool = True, **kwargs: P.kwargs
) -> list:
    """
    Find a set of generators that cover the rectangle.

    Note 1: We assume that the generators are named in the form 'a', 'A', 'b', 'B', 'ab', 'ba',
        'AB', 'BA', 'bA', 'aB', 'Ba', 'Ab'.
    Note: WE do this by translating up/down left/right with parabolic elements.

    EXAMPLES:

        sage: from maass_forms_klein.hyperbolic_space.utils import find_covering_generators2
        sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram, Rectangle
        sage: s = matrix([[0, -1], [1, 0]])
        sage: t = matrix([[1, 1], [0, 1]])
        sage: l = matrix([[1, I], [0, 1]])
        sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: find_covering_generators2(rect, {'a': s,'b':s,'B':s,'A':s})
        Traceback (most recent call last):
        ...
        ValueError: Incomplete generators...
        sage: gens = {'a': s,'b':s,'B':s,'A':s, 'm':t, 'M': t**-1, 'l':l, 'L':l**-1}
        sage: find_covering_generators2(rect, gens, return_words=True)
        ['A', 'AL', 'AM', 'ALM']
        sage: rect = Rectangle(base=vector((-0.5, -0.5)), v1=vector((1, 0)), v2=vector((0, 1)))
        sage: find_covering_generators2(rect, gens, return_words=True)
        ['A']
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import (
        ....:     KleinianGroup__from_manifold)

        sage: from snappy import Manifold
        sage: M = Manifold("4_1")
        sage: G = KleinianGroup__from_manifold(M)
        sage: gens = G.named_generators(prec=53)
        sage: v1 = vector((3.46410161513775, 0))
        sage: rect = Parallelogram(base=vector((0, 0)), v1=v1, v2=vector((0, 1)))
        sage: find_covering_generators2(rect, gens, return_words=True) # doctest: +SKIP
        ['b', 'bm', 'BAl', 'aabMMM', 'baabMM', 'bABBlm', 'BABBlmm', 'bABBlmm']
        sage: rect =  G.translation_fundamental_domain()
        sage: find_covering_generators2(rect, gens, return_words=True)
        ['BAB', 'bab', 'BABm', 'babM', 'bABBm', 'aabMMM', 'aabMMML']
    """
    verbose = kwargs.get("verbose", False)
    if not isinstance(named_gens, dict):
        raise ValueError("named_gens must be a dictionary.")
    if not isinstance(rect, (Rectangle, Parallelogram)):
        raise ValueError(f"rect must be a Rectangle or Parallelogram: {type(rect)}")
    all_gens = ["a", "A", "b", "B", "m", "M", "l", "L"]
    if set(named_gens) != set(all_gens):
        raise ValueError(f"Incomplete generators. Need all of {all_gens}")
    word_gens = "aAbB"
    orig_list = []
    for word_length in range(1, kwargs.get("max_word_length", 5)):
        for word in FiniteWords(word_gens).iterate_by_length(word_length):
            # Reduce words first
            orig_list.append(reduce_word(str(word)))
    orig_list = list(set(orig_list))
    # Make a first filter
    max_radius = 10 * max(1, max(v.norm() for v in rect.vertices))
    if verbose > 0:
        log.debug("orig_list=", len(orig_list), orig_list)
    circles = {x: word_to_circle(x, named_gens) for x in orig_list}
    # Detect parabolic elements numerically
    orig_list = [reduce_word(x) for x in orig_list if abs(circles[x].radius) < max_radius]
    orig_list = list(set(orig_list))
    orig_list.sort(key=word_list_sort_key)
    method = kwargs.get("method", "default")
    if verbose > 0:
        log.debug("orig_list=", len(orig_list), orig_list)
    new_list = []
    for x in orig_list:
        cx = circles[x]
        add_x = True
        # Check if circle is different from new circles
        for y in new_list:
            cy = circles[y]
            if cy == cx or circle_approx_equal(cy, cx):
                add_x = False
                break
        if not add_x:
            continue
        # Check if circle is inside other new circles
        for y in new_list:
            cy = circles[y]
            if circle_is_contained(cx, cy):
                add_x = False
                break
        if add_x:
            new_list.append(x)
    new_list.sort(key=word_list_sort_key)
    if verbose > 0:
        log.debug("new_list=", len(new_list), new_list)
    # Now add reductions of these circles to the fundamental domain.
    new_reduced_list = []
    new_list_centers = [circles[x].center for x in new_list]
    if method == "default":
        reduced_list_centers, translations = reduce_in_parallelogram(rect, new_list_centers, True)
    else:
        # Make more reductions
        numpts = kwargs.get("numpts", 50)
        newer_list = []
        reduced_list_centers = []
        translations = []
        new_list_centers = []
        for x in new_list:
            c = circles[x]
            newer_list += [x] * (2 * numpts + 1) * 2
            pts = [
                c.center + i * c.radius / numpts * vector((1, 0))
                for i in range(-numpts, numpts + 1)
            ]
            pts += [
                c.center + i * c.radius / numpts * vector((0, 1))
                for i in range(-numpts, numpts + 1)
            ]
            reduced_list_centers_x, translations_x = reduce_in_parallelogram(rect, pts, True)
            if x == "a":
                for _n, xx in enumerate(pts):
                    log.debug(x, c, xx, point_in_circle(xx, c, scaling_factor=1.01))
            reduced_list_centers += reduced_list_centers_x
            translations += translations_x
            new_list_centers += pts

        new_list = newer_list
    # new_list_centers = [circles[x].center for x in new_list]
    for n, x in enumerate(new_list):
        if verbose > 0:
            log.debug(
                "n= %s x= %s %s %s %s",
                n,
                x,
                new_list_centers[n],
                translations[n],
                reduced_list_centers[n],
            )
        # Check
        if not point_in_circle(new_list_centers[n], circles[x], scaling_factor=1.01):
            if verbose > 0:
                log.debug("Not in circle")
            raise ArithmeticError("Point is not in circle.")
        w = translation_tuple_to_word(translations[n], {0: "L", 1: "M"})
        new_word = f"{x}{w}"
        if verbose > 0:
            log.debug("new_words=", new_word, word_to_circle(new_word, named_gens).center)
        new_circle = word_to_circle(new_word, named_gens)
        if not parallelogram_intersect_circle(rect, new_circle, scaling_factor=1):
            if verbose > 0:
                log.debug("Not intersecting:")
                log.debug("new word=", new_word, new_circle)
            continue
        new_reduced_list.append(new_word)
        # If the center is a vertex of the parallelogram, we add all translations of it
        if x == "ab":
            if verbose > 0:
                log.debug(
                    "on bd= %s",
                    point_on_boundary_of_parallelogram(reduced_list_centers[n], rect),
                )
                log.debug(
                    "on bd= %s",
                    point_on_boundary_of_parallelogram(new_list_centers[n], rect),
                )
        if verbose:
            log.debug("r=", reduced_list_centers[n], rect)
            log.debug("on bd=", point_on_boundary_of_parallelogram(reduced_list_centers[n], rect))
        if point_on_boundary_of_parallelogram(reduced_list_centers[n], rect):
            if verbose > 0:
                log.debug("pt on boundary=", reduced_list_centers[n])
            for w1 in ["L", "M", "m", "l", "LM", "lm"]:
                if w1 == w:
                    continue
                other_translated_word = f"{new_word}{w1}"
                circle = word_to_circle(other_translated_word, named_gens)
                if x == "ab":
                    if verbose > 0:
                        log.debug("Add word", other_translated_word, circle.center)
                        log.debug("on bd=", point_on_boundary_of_parallelogram(circle.center, rect))
                if point_on_boundary_of_parallelogram(circle.center, rect):
                    if verbose > 0:
                        log.debug("Add word", other_translated_word, circle.center)
                    new_reduced_list.append(other_translated_word)
    new_reduced_list = [reduce_word(x) for x in new_reduced_list]
    new_reduced_list = list(set(new_reduced_list))
    if verbose > 0:
        log.debug("new_reduced_list=", len(new_reduced_list), new_reduced_list)
    # min_radius = min([abs(word_to_circle(m, named_gens).radius) for m in origlist])
    try:
        test = parallelogram_covered_by_circles(
            rect,
            [word_to_circle(x, named_gens) for x in new_reduced_list],
            scaling_factor=0.999,
            n_max=kwargs.get("n_max", 20),
            verbose=verbose - 1 if verbose > 0 else 0,
        )
    except ArithmeticError as e:
        log.info(str(e))
        test = False
    if test:
        covering_generators = reduce_cover(rect, new_reduced_list, named_gens)
    else:
        raise ArithmeticError(
            f"The parallelogram {rect} is not covered by the circles:{new_reduced_list}"
        )
    # Remove duplicates
    covering_generators = [reduce_word(w) for w in covering_generators]
    covering_generators = list(set(covering_generators))
    # Reduce the cover and see if still covering
    if verbose > 0:
        log.debug("covering_generators1=", covering_generators)
    covering_generators = reduce_cover(rect, covering_generators, named_gens)
    if verbose > 0:
        log.debug("covering_generators2=", covering_generators)
    if covering_generators:
        if check and not check_generators(covering_generators, named_gens=named_gens):
            # First attempt we add back one of the generators
            covering_generators += ["a"]
            if verbose:
                log.debug("Check failed first. Adding back 'a'")
            if not check_generators(covering_generators, named_gens=named_gens):
                raise ArithmeticError(
                    f"Covering generators {covering_generators} do not generatethe same group."
                )

    covering_generators.sort(key=word_list_sort_key)
    return covering_generators


# def find_covering_generators2(rect: Rectangle | Parallelogram, named_gens: dict,
#                              check: bool = True,
#                              **kwargs: P.kwargs) -> list:
#     """
#     Find a set of generators that cover the rectangle.
#
#     Note 1: We assume that the generators are named in the form 'a', 'A', 'b', 'B', 'ab', 'ba',
#         'AB', 'BA', 'bA', 'aB', 'Ba', 'Ab'.
#     Note: WE do this by translating up/down left/right with parabolic elements.
#
#     EXAMPLES:
#         sage: from maass_forms_klein.hyperbolic_space.utils import find_covering_generators2
#         sage: from maass_forms_klein.hyperbolic_space.types import Parallelogram, Rectangle
#         sage: s = matrix([[0, -1], [1, 0]])
#         sage: t = matrix([[1, 1], [0, 1]])
#         sage: l = matrix([[1, I], [0, 1]])
#         sage: rect = Rectangle(base=vector((0, 0)), v1=vector((1, 0)), v2=vector((0, 1)))
#         sage: find_covering_generators2(rect, {'a': s,'b':s,'B':s,'A':s})
#         Traceback (most recent call last):
#         ...
#         ValueError: Incomplete generators. Need all of ['a', 'A', 'b', 'B', 'm', 'M', 'l', 'L']
#         sage: gens = {'a': s,'b':s,'B':s,'A':s, 'm':t, 'M': t**-1, 'l':l, 'L':l**-1}
#         sage: find_covering_generators2(rect, gens)
#         ['a', 'laL', 'maM', 'lmaML']
#         sage: from maass_forms_klein.hyperbolic_space.kleinian_group import (
#         ....:     KleinianGroup__from_manifold)
#         sage: from snappy import Manifold
#         sage: M = Manifold("4_1")
#         sage: G = KleinianGroup__from_manifold(M)
#         sage: gens = G.named_generators(prec=53)
#         sage: v1 = vector((3.46410161513775, 0))
#         sage: rect = Parallelogram(base=vector((0, 0)), v1=v1, v2=vector((0, 1)))
#         sage: find_covering_generators2(rect, gens)
#         ['b', 'LBl', 'Lbl', 'Mbm', 'LBAl', 'Lbal', 'mabM', 'LMbml']
#         sage: rect =  G.translation_fundamental_domain()
#         sage: find_covering_generators2(rect, gens)
#         ['B', 'b', 'BA', 'LBl', 'mabM', 'mBAM', 'mmabMM']
#
#     """
#     verbose = kwargs.get('verbose', False)
#     if not isinstance(named_gens, dict):
#         raise ValueError("named_gens must be a dictionary.")
#     if not isinstance(rect, (Rectangle, Parallelogram)):
#         raise ValueError(f"rect must be a Rectangle or Parallelogram: {type(rect)}")
#     all_gens = ['a', 'A', 'b', 'B', 'm', 'M', 'l', 'L']
#     if set(named_gens) != set(all_gens):
#         raise ValueError(f"Incomplete generators. Need all of {all_gens}")
#     origlist = ['a', 'A', 'b', 'B']
#     origlist += [normalize_word(x)
#                            for x in FiniteWords("".join(origlist)).iterate_by_length(2)]
#     origlist = [x for x in origlist
#                 if word_to_circle(x, named_gens).radius != Infinity]
#     # origlist = ['a', 'A', 'b', 'B', 'ab', 'ba', 'AB', 'BA', 'bA', 'aB', 'Ba', 'Ab']
#     min_radius = min([abs(word_to_circle(m, named_gens).radius) for m in origlist])
#
#     current_matrices = {x: word_to_element(x, named_gens) for x in origlist}
#     # filter out circles that have radius smaller than 1 or are parabolic
#     current_matrices = {x: m for x, m in current_matrices.items() if m[1][0] != 0 and
#                         matrix_to_circle(m).radius >= kwargs.get('matrix_min', 0)}
#     origlist = list(current_matrices)
#     if verbose:
#         print("ORIG LIST=", origlist)
#     origlist = filter_list_mod_parabolics(origlist, gens=named_gens)
#     origlist.sort(key=word_list_sort_key)
#     if verbose:
#         print("reduced 1 list=", origlist)
#     F = FiniteWords('lmLM')
#     current_words_list = deepcopy(origlist)
#     current_matrices_list = [word_to_element(x, named_gens) for x in origlist]
#     covering_generators = None
#     # What is the max word length needed?
#     # we need to be able to translate all circles along the fundamental domain
#     max_length = ceil( max(rect.sides) / min_radius / 2)
#     print("max length=", max_length)
#     current_words_list_len = 0
#     for n in range(1, kwargs.get('max_length', max_length)+1):
#         if verbose:
#             print("word length=", n)
#         new_list = []
#         for w in F.iterate_by_length(n):
#             w = normalize_word(w)
#             new_list += [str(w) + x + find_inverse_word(str(w)) for x in origlist]
#         # Remove duplicates
#         reduced_list = list(set([reduce_word(w) for w in new_list]))
#         if verbose:
#             print("reduced_list1=", reduced_list)
#         reduced_list = filter_list_mod_parabolics(reduced_list, gens=named_gens)
#         if verbose:
#             print("reduced_list2=", reduced_list)
#
#         reduced_list.sort(key=word_list_sort_key)
#         new_list_matrix = {x: word_to_element(x, named_gens) for x in reduced_list}
#         # Filter out generators that correspond to existing matrices
#         new_list_matrix = {x: m for x, m in new_list_matrix.items() if
#                            not circle_in_list(m, current_matrices_list)}
#         # Filter out circles that are not intersecting the rectangle
#         # if verbose:
#         #     print("new_list1=", new_list_matrix)
#         if isinstance(rect, Rectangle):
#             new_list_matrix = {x: m for x, m in new_list_matrix.items() if
#                                m != identity_matrix(2) and
#                                m != -identity_matrix(2) and
#                                rectangle_intersect_circle(rect, matrix_to_circle(m))
#                                }
#         else:
#             new_list_matrix = {x: m for x, m in new_list_matrix.items() if
#                                m != identity_matrix(2) and
#                                m != -identity_matrix(2) and
#                                parallelogram_intersect_circle(rect, matrix_to_circle(m))
#                                }
#         # if verbose:
#         #     print("new_list2=", new_list_matrix)
#         current_words_list += new_list_matrix.keys()
#         current_words_list = filter_list_mod_parabolics(current_words_list, gens=named_gens)
#         current_words_list.sort()
#         if verbose:
#             print("xcurrent_words_list=", current_words_list)
#         # if len(current_words_list) == current_words_list_len:
#         #     raise ArithmeticError("No new words!")
#         current_words_list_len = len(current_words_list)
#         current_matrices = {x: word_to_element(x, named_gens) for x in current_words_list}
#         # To avoid vertices in C we check for coverage by slightly shrunken circles.
#         try:
#             test = parallelogram_covered_by_circles(rect, [word_to_circle(x, named_gens)
#                                                            for x in current_words_list],
#                                                     scaling_factor=0.999,
#                                                     n_max=kwargs.get('n_max', 20),
#                                                     verbose=verbose-1 if verbose > 0 else 0)
#             if verbose:
#                 print("test=", test)
#         except ArithmeticError as e:
#             print(e)
#             test = False
#         if test:
#             covering_generators = reduce_cover(rect, list(current_matrices.keys()), named_gens)
#             break
#     if covering_generators:
#         if check and not check_generators(covering_generators, named_gens=named_gens):
#             # First attempt we add back one of the generators
#             covering_generators += ['a']
#             if verbose:
#                 print("Check failed first. Adding back 'a'")
#             if not check_generators(covering_generators, named_gens=named_gens):
#                 raise ArithmeticError(f"Covering generators {covering_generators} do not generate"
#                                   f"the same group.")
#
#         return covering_generators
#
#     raise ArithmeticError("Could not find covering set of generators using translations.")


def check_generators(sub_gens: list[str], named_gens: dict | None = None) -> bool:
    """
    Check that the list of new generators still generates the same group, i.e. that
    it contains the a and b generators.

    INPUT:
    - sub_gens: list of strings in 'a', 'b', 'A', 'B', 'm', 'M', 'l', 'L'
      corresponding to generators

    Note:

    EXAMPLES::

    sage: from maass_forms_klein.hyperbolic_space.utils import check_generators
    sage: from sage.matrix.constructor import matrix
    sage: from sage.rings.complex_mpfr import ComplexField
    sage: # Test check_generators with mock named generators
    sage: CC = ComplexField(53)
    sage: mock_gens = {
    ....:     'a': matrix(CC, [[1, 1], [0, 1]]),
    ....:     'b': matrix(CC, [[1, 0], [1, 1]])
    ....: }
    sage: new_gens = ['b', 'ab', 'ba']
    sage: check_generators(new_gens, named_gens=mock_gens)
    True

    """
    if "a" in named_gens and "b" in named_gens and named_gens["a"] != named_gens["b"]:
        FG = FreeGroup(["a", "b", "m", "l"])
        a, b, m, l = FG.gens()  # noqa: E741
    elif "a" in named_gens:
        FG = FreeGroup(["a", "m", "l"])
        a, m, l = FG.gens()  # noqa: E741
        b = FG.one()
    else:
        raise ValueError("Need at least 'a' in named gens.")
    gens_dict = {"a": a, "b": b, "A": a**-1, "B": b**-1, "m": m, "M": m**-1, "l": l, "L": l**-1}
    new_gens = [m, l]
    for word in sub_gens:
        g = FG.one()
        for letter in word:
            g *= gens_dict.get(letter)
        if g != FG.one():
            new_gens.append(g)
    SUBG = FG.subgroup(new_gens)
    return a in SUBG and b in SUBG


def find_side_pairing_gens(manifold: "Manifold") -> list[str]:
    """
    Find the side pairing generators for a manifold.

    INPUT:
    - ``manifold`` -- the manifold

    EXAMPLES::

        sage: # Skipping examples that depend on snappy and Kleinian groups  # doctest: +SKIP
        sage: from maass_forms_klein.hyperbolic_space.kleinian_group import KleinianGroup
        sage: from maass_forms_klein.hyperbolic_space.utils import find_side_pairing_gens
        sage: from snappy import Manifold
        sage: M = Manifold('4_1')
        sage: find_side_pairing_gens(M)
        ['AB', 'ABB', 'B', 'BAB', 'b', 'bABB', 'ba', 'bab', 'babABB', 'bba', 'bbaB', 'bbaBAB']
    """
    # Try to find side-pairing generators
    D = manifold.dirichlet_domain(include_words=True)
    side_pairing_dirichlet = D.pairing_words()
    G = manifold.fundamental_group()
    # These are given in terms of the original generators, so
    # we need to translate to a and b etc.
    az = string.ascii_lowercase
    original_gens = {az[n]: x for n, x in enumerate(G.original_generators())}
    return sorted(
        [reduce_word(change_letters_in_word(x, original_gens)) for x in side_pairing_dirichlet]
    )
