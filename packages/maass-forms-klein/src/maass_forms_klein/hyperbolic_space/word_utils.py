import re

from maass_forms_klein.hyperbolic_space.types import Circle
from sage.matrix.constructor import matrix
from sage.structure.element import Matrix


def find_inverse_word(word: str) -> str:
    """
    Invert word by interchanging upper and lower case and revert.

    INPUT:

    - ``word`` -- (string) word

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.word_utils import find_inverse_word
        sage: find_inverse_word('aB')
        'bA'
        sage: find_inverse_word('AA')
        'aa'
    """
    return word.swapcase()[::-1]


def is_equivalent_mod_parabolic(w1, w2, gens):
    """
    Check if two words are equivalent modulo parabolic elements (on the left),
    i.e. a ~ b if a = t * b for a parabolic t.

    INPUT:

    - ``w1`` -- (string) first word
    - ``w2`` -- (string) second word
    - ``gens`` -- (dictionary) generators of the group

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.word_utils import is_equivalent_mod_parabolic
        sage: is_equivalent_mod_parabolic('a', 'b', {'a': matrix([[0, 1], [1, 0]]),
        ....:                                            'b': matrix([[0, 1], [1, 0]])})
        True


    """
    mat = word_to_element(w1, gens) * word_to_element(w2, gens) ** -1
    try:
        eps = mat.base_ring().epsilon()
    except NotImplementedError:
        eps = 0
    if abs(mat[1][0]) <= 2 * eps:
        return True
    return False


def word_list_sort_key(word):
    """
    Sort key for list of words so that we first sort by word length, then alphanum (case-insensitive)
    and finally by case.

    """
    return len(word), word.lower(), word

def word_to_str(word):
    """
    Convert word to string.

    """
    return str(word)

def normalize_word(word):
    """
    Normalize a word.

    """
    lw = list(word)
    lw.sort(key=lambda x: x.lower())
    return "".join(lw)


def reduce_word(word, n: int = 0):
    """
    Reduce a word with respect to inverses only assuming that inverse elements
    are represented by upper/lower versions of the same letter.

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.word_utils import reduce_word
        sage: reduce_word('abAmlB')
        'abAmlB'
        sage: reduce_word('aABb')
        ''
        sage: reduce_word('abAaB')
        'a'
    """
    len_word_in = len(word)
    if len_word_in <= 1:
        return word
    if n and 2 * n > len_word_in:
        return word
    # hard-coded list is quicker ...
    for pair in ["aA", "bB", "lL", "mM", "Aa", "Bb", "Ll", "Mm"]:
        word = word.replace(pair, "")
    # If no change in word we just return it
    if len_word_in == len(word):
        return word
    # Otherwise, we might need to do another reduction
    return reduce_word(word, n + 1)

def expand_parentheses(word):
    """
    Expand parentheses in a word.

    INPUT:

    - ``word`` -- (string) word

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.word_utils import expand_parentheses
        sage: expand_parentheses('(ab)^2a^-1')
        'ababa^-1'
        sage: expand_parentheses('(a(ab)^2)^2')
        'aababaabab'
    """
    if ("(" in word and ")" not in word) or (")" in word and "(" not in word):
        raise ValueError(f"Unmatched parentheses in {word}")
    if "(" not in word:
        return word
    for w, n in re.findall(r"\(([^()]*)\)\^(-?\d+)", word):
        wnew = expand_parentheses(w)
        if int(n) < 0:
            wnew = wnew.swapcase()
        word = word.replace(f"({w})^{n}", wnew * abs(int(n)))
    return expand_parentheses(word)

def expand_word(word: str) -> str:
    """
    Convert a word of the form (ab)^2*a^-1 to abababA etc.

    INPUT:

    - ``word`` -- (string) word

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.word_utils import expand_word
        sage: expand_word('abAmlB')
        'abAmlB'
        sage: expand_word('a^3*b^-1*a*b^-10')
        'aaaBaBBBBBBBBBB'
    """
    word = word.replace("*", "")
    word = expand_parentheses(word)
    # Replace any a^2 with aa etc.
    replacements = re.findall(r"(\w)\^(\d+)", word)
    # Need to make sure that we match b^10 before b^1
    replacements.sort(key=lambda x: len(x[1]), reverse=True)
    for w, n in replacements:
        word = word.replace(f"{w}^{n}", w * int(n))
    # Replace any a^-2 with AA etc.
    replacements = re.findall(r"(\w)\^(-\d+)", word)
    replacements.sort(key=lambda x: len(x[1]), reverse=True)
    for w, n in replacements:
        winv = w.swapcase()
        word = word.replace(f"{w}^{n}", winv * abs(int(n)))
    return word
def word_to_element(word: str, gens: dict) -> Matrix:
    """
    Convert a word to a matrix.

    INPUT:

    - ``word`` -- (string) word
    - ``gens`` -- (dictionary) generators of the group

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.word_utils import word_to_element
        sage: word_to_element('a', {'a': matrix([[0, 1], [1, 0]])})
        [0 1]
        [1 0]
        sage: from maass_forms_klein.hyperbolic_space.word_utils import word_to_element
        sage: word_to_element('a^3', {'a': matrix([[0, 1], [1, 0]])})
        [0 1]
        [1 0]

    """
    if not gens:
        raise ValueError("Need non-empty generators")
    g = matrix([[1, 0], [0, 1]])
    word = expand_word(word)
    for w in word:
        try:
            g = g * gens[w]
        except KeyError:
            raise ValueError(f"Generator dict has no key `{w}`")
    return g


def word_to_circle(word: str, gens: dict) -> Circle:
    matrix = word_to_element(word, gens)
    from maass_forms_klein.hyperbolic_space.geometry_utils import matrix_to_circle
    return matrix_to_circle(matrix)


def translation_tuple_to_word(t: tuple, gens: dict) -> str:
    r"""
    Map a tuple (t_1, ..., t_n) to a string g1..g1g2...g2...gn...gn)
    where each gi is repeated |t_i| times and if t_i is negative then the case of gi is swapped to indicate an inverse.

    EXAMPLE:

    sage: from maass_forms_klein.hyperbolic_space.word_utils import translation_tuple_to_word
    sage: translation_tuple_to_word( (1,-2),{0:'L',1:'M'})
    'Lmm'
    sage: translation_tuple_to_word( (1,2),{0:'L',1:'M'})
    'LMM'
    """
    return "".join(
        [(gens[i] if t[i] > 0 else gens[i].swapcase()) * abs(t[i]) for i in range(len(t))])


def change_letters_in_word(word: str, new_names: dict) -> str:
    """
    Change the letters in a word using a dictionary and use the convention
    that the upper-case letters are inverses to lower case letters.

    INPUT:

    - ``word`` -- (string) word
    - ``new_names`` -- (dictionary) dictionary of new names

    EXAMPLES::

        sage: from maass_forms_klein.hyperbolic_space.word_utils import change_letters_in_word
        sage: change_letters_in_word('aC', {'a': 'A', 'c': 'B'})
        'Ab'
    """

    new_word = []
    for letter in word:
        if letter.isupper():
            new_letter = new_names[letter.lower()].swapcase()
        else:
            new_letter = new_names[letter]
        new_word.append(new_letter)
    return "".join(new_word)