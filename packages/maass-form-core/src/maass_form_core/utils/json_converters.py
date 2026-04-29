"""
JSON serialization utilities for SageMath mathematical objects.

Provides converters for rings, ring elements, matrices, number fields,
and complex numbers to and from JSON-compatible dictionaries.

EXAMPLES::

    sage: from maass_form_core.utils.json_converters import ring_to_json, ring_from_json
    sage: ring_to_json(RealField(53))
    {'__type__': 'ring', 'name': 'RealField', 'prec': 53}
    sage: ring_from_json({'name': 'RealField', 'prec': 53})
    Real Field with 53 bits of precision
"""

import json
from json import JSONDecodeError
from typing import Any

from sage.all import ZZ
from sage.matrix.constructor import matrix
from sage.rings.complex_mpfr import ComplexField, ComplexField_class, ComplexNumber
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing, IntegerRing_class
from sage.rings.number_field.number_field import NumberField, NumberField_generic
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.rational import Rational
from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RealField, RealField_class, RealNumber
from sage.structure.element import Matrix, Vector


def get_prec(x: Any) -> int:
    r"""
    Return the precision of a numerical value.

    INPUT:

    - ``x`` -- a numerical value

    OUTPUT:

    - int; the precision in bits, or 0 for exact types, or 53 as default

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import get_prec
        sage: get_prec(RR(1.0))
        53
        sage: get_prec(RealField(100)(1.0))
        100
        sage: get_prec(ZZ(1))
        0
        sage: get_prec(QQ(1/2))
        0
        sage: get_prec(1)
        0
        sage: get_prec(1.0)
        53
    """
    if isinstance(x, (RealNumber, ComplexNumber)):
        return int(x.prec())
    if isinstance(x, (Integer, int, Rational)):
        return 0
    return 53


class SageJSONEncoder(json.JSONEncoder):
    r"""
    Custom JSON encoder for SageMath types.

    Handles Integer, Real, Complex, Rational, NumberField elements,
    matrices, vectors, and dicts with non-string keys.

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import SageJSONEncoder
        sage: import json
        sage: json.loads(json.dumps(ZZ(42), cls=SageJSONEncoder))['value']
        '42'
    """

    def default(self, o):
        if isinstance(o, str):
            return o
        if isinstance(
            o,
            (
                IntegerRing_class,
                RealField_class,
                ComplexField_class,
                RationalField,
                NumberField_generic,
            ),
        ):
            return ring_to_json(o)
        if isinstance(o, (float, int, complex)):
            return o
        if isinstance(o, (Integer, RealNumber, ComplexNumber, Rational, NumberFieldElement)):
            return ring_element_to_json(o)
        if isinstance(o, Matrix):
            return matrix_to_json(o)
        if isinstance(o, (Vector, tuple)):
            return list(o)
        if isinstance(o, dict):
            return dict_to_json(o)
        return super().default(o)


def decode_function(obj: dict | str) -> Any:
    r"""
    Decode a typed JSON object back to a SageMath object.

    INPUT:

    - ``obj`` -- dict or str; a JSON object potentially containing '__type__' key

    OUTPUT:

    - The decoded SageMath object, or the original object if not a typed JSON

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import decode_function
        sage: decode_function({'__type__': 'matrix',
        ....:                   'base_ring': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
        ....:                    'entries': [['1', '2'], ['3', '4']]})
        [1 2]
        [3 4]
        sage: decode_function('hello')
        'hello'
        sage: decode_function({"__type__": "element",
        ....:                   "parent": {"__type__": "ring", "name": "IntegerRing", "prec": 0},
        ....:                   "value": "1"})
        1
    """
    if isinstance(obj, dict) and "__type__" in obj:
        if obj["__type__"] == "ring":
            return ring_from_json(obj)
        elif obj["__type__"] == "matrix":
            return matrix_from_json(obj)
        elif obj["__type__"] == "element":
            return ring_element_from_json(obj)
    if isinstance(obj, dict):
        new_dict = {}
        for k, v in obj.items():
            try:
                k_new = json.loads(k, object_hook=decode_function)
            except (JSONDecodeError, TypeError):
                k_new = k
            try:
                v_new = json.loads(v, object_hook=decode_function)
            except TypeError:
                v_new = v
            except JSONDecodeError:
                v_new = v
            if isinstance(k_new, list):
                k_new = tuple(k_new)
            new_dict[k_new] = v_new
    if isinstance(obj, list):
        return tuple([json.loads(x, object_hook=decode_function) for x in obj])
    return obj


# --- Ring serialization ---


def ring_to_json(F) -> dict:
    r"""
    JSON representation of a base ring.

    INPUT:

    - ``F`` -- a SageMath ring (IntegerRing, RationalField, RealField, ComplexField,
      or NumberField)

    OUTPUT:

    - dict; JSON-serializable dictionary with '__type__', 'name', and ring-specific fields

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import ring_to_json
        sage: ring_to_json(RationalField())
        {'__type__': 'ring', 'name': 'RationalField', 'prec': 0}
        sage: ring_to_json(RealField(53))
        {'__type__': 'ring', 'name': 'RealField', 'prec': 53}
        sage: ring_to_json(ComplexField(53))
        {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53}
        sage: ring_to_json(ZZ)
        {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0}
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: ring_to_json(F)
        {'__type__': 'ring',
         'field': {'__type__': 'ring',
          'embedding': None,
          'names': ['z'],
          'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2'},
         'name': 'NumberField'}
    """
    if isinstance(F, RealField_class):
        return {"name": "RealField", "prec": int(F.prec()), "__type__": "ring"}
    if isinstance(F, ComplexField_class):
        return {"name": "ComplexField", "prec": int(F.prec()), "__type__": "ring"}
    if isinstance(F, RationalField):
        return {"name": "RationalField", "prec": 0, "__type__": "ring"}
    if isinstance(F, IntegerRing_class):
        return {"name": "IntegerRing", "prec": 0, "__type__": "ring"}
    if isinstance(F, NumberField_generic):
        return {"name": "NumberField", "field": number_field_to_json(F), "__type__": "ring"}
    raise ValueError(f"Unsupported base ring {F}")


def ring_from_json(data: dict | str):
    r"""
    Construct ring from JSON data.

    INPUT:

    - ``data`` -- dict or str; JSON representation of a ring

    OUTPUT:

    - A SageMath ring (RealField, ComplexField, RationalField, IntegerRing, or NumberField)

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import ring_from_json, number_field_to_json
        sage: ring_from_json({'name': 'RealField', 'prec': 53})
        Real Field with 53 bits of precision
        sage: ring_from_json({'name': 'ComplexField', 'prec': 53})
        Complex Field with 53 bits of precision
        sage: ring_from_json({'name': 'RationalField'})
        Rational Field
        sage: ring_from_json({'name': 'IntegerRing'})
        Integer Ring
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: ring_from_json({'name': 'NumberField', 'field': number_field_to_json(F)})
        Number Field in z with defining polynomial x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2
    """
    if isinstance(data, (RealField_class, ComplexField_class, RationalField, IntegerRing_class)):
        return data
    if isinstance(data, str):
        data = json.loads(data)
    if data["name"] == "RealField":
        return RealField(data["prec"])
    if data["name"] == "ComplexField":
        return ComplexField(data["prec"])
    if data["name"] == "RationalField":
        return RationalField()
    if data["name"] == "IntegerRing":
        return IntegerRing()
    if data["name"] == "NumberField":
        return number_field_from_json(data["field"])
    raise ValueError(f"Unsupported base ring {data}")


# --- Ring element serialization ---


def ring_element_to_json(data: Any) -> dict:
    r"""
    JSON representation of a ring element.

    INPUT:

    - ``data`` -- a ring element (Integer, Real, Complex, Rational, or NumberField element)

    OUTPUT:

    - dict; JSON-serializable dictionary with '__type__', 'parent', and 'value' keys

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import ring_element_to_json
        sage: ring_element_to_json(ZZ(1))
        {'__type__': 'element',
         'parent': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
         'value': '1'}
        sage: ring_element_to_json(RR(3.14))['parent']['name']
        'RealField'
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: ring_element_to_json(F.gen())['value']
        'z'
    """
    base_ring = None
    if isinstance(data, (Integer, int)):
        base_ring = {"name": "IntegerRing", "prec": get_prec(data), "__type__": "ring"}
    if isinstance(data, (RealNumber, float)):
        base_ring = {"name": "RealField", "prec": get_prec(data), "__type__": "ring"}
    if isinstance(data, (ComplexNumber, complex)):
        base_ring = {"name": "ComplexField", "prec": get_prec(data), "__type__": "ring"}
    if isinstance(data, Rational):
        base_ring = {"name": "RationalField", "prec": get_prec(data), "__type__": "ring"}
    if isinstance(data, NumberFieldElement):
        base_ring = number_field_to_json(data.parent())
    return {"parent": base_ring, "value": str(data), "__type__": "element"}


def ring_element_from_json(data: dict | str) -> Any:
    r"""
    Convert JSON representation of a ring element to a SageMath object.

    INPUT:

    - ``data`` -- dict or str; JSON representation of a ring element

    OUTPUT:

    - A SageMath ring element

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import ring_element_from_json
        sage: ring_element_from_json({'parent': {'name': 'IntegerRing', 'prec': 0}, 'value': '1'})
        1
        sage: ring_element_from_json({'parent': {'name': 'RealField', 'prec': 53}, 'value': '3.14'})
        3.14000000000000
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: z = F.gen()
        sage: ring_element_from_json({'parent': {'field': {'names': ['z'],
        ....: 'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2'},
        ....: 'name': 'NumberField'},
        ....: 'value': 'z^2'})
        z^2
    """
    if not isinstance(data, dict):
        data = json.loads(data)
    parent = ring_from_json(data["parent"])
    return parent(data["value"])


# --- Matrix serialization ---


def matrix_to_json(m) -> dict:
    r"""
    JSON representation of a matrix.

    INPUT:

    - ``m`` -- a SageMath matrix

    OUTPUT:

    - dict; JSON-serializable dictionary with '__type__', 'base_ring', and 'entries' keys

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import matrix_to_json
        sage: m = matrix(ZZ, [[1, 2], [3, 4]])
        sage: matrix_to_json(m)
        {'__type__': 'matrix',
         'base_ring': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
         'entries': [['1', '2'], ['3', '4']]}
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: z = F.gen()
        sage: d = matrix_to_json(matrix(F, [[z^2, z^3], [z^4, z^5]]))
        sage: d['base_ring']['name']
        'NumberField'
    """
    return {
        "base_ring": ring_to_json(m.base_ring()),
        "entries": [[str(x) for x in row] for row in list(m)],
        "__type__": "matrix",
    }


def matrix_from_json(data: dict | str):
    r"""
    Convert JSON representation of a matrix to a SageMath matrix.

    INPUT:

    - ``data`` -- dict or str; JSON representation of a matrix

    OUTPUT:

    - A SageMath matrix

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import matrix_from_json
        sage: matrix_from_json({'base_ring': {'name': 'IntegerRing', 'prec': 0},
        ....:                    'entries': [['1', '2'], ['3', '4']]})
        [1 2]
        [3 4]
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: matrix_from_json({'base_ring': {'field': {'names': ['z'],
        ....: 'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2'},
        ....: 'name': 'NumberField'},
        ....: 'entries': [['z^2', 'z^3'], ['z^4', 'z^5']]})
        [z^2 z^3]
        [z^4 z^5]
    """
    if isinstance(data, str):
        data = json.loads(data)
    return matrix(ring_from_json(data["base_ring"]), data["entries"])


# --- Number field serialization ---


def number_field_to_json(nf: NumberField_generic) -> dict:
    r"""
    JSON representation of a number field.

    INPUT:

    - ``nf`` -- a NumberField

    OUTPUT:

    - dict; JSON-serializable dictionary with polynomial, names, and optional embedding

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import number_field_to_json
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: number_field_to_json(F)
        {'__type__': 'ring',
         'embedding': None,
         'names': ['z'],
         'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2'}
        sage: gen = CC('-0.3122516446211537 + 1.026735437750787*I')
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z', embedding=gen)
        sage: d = number_field_to_json(F)
        sage: d['embedding'] is not None
        True
    """
    return {
        "polynomial": str(nf.polynomial()),
        "names": list(nf._names),
        "__type__": "ring",
        "embedding": str(nf.gen_embedding().n(53)) if nf.gen_embedding() else None,
    }


def number_field_from_json(data: dict | str) -> NumberField_generic:
    r"""
    Create a number field from JSON data.

    INPUT:

    - ``data`` -- dict or str; JSON representation of a number field

    OUTPUT:

    - A NumberField

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import number_field_from_json
        sage: number_field_from_json({'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2',
        ....:                          'names': ['a']})
        Number Field in a with defining polynomial x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2
        sage: from maass_form_core.utils.json_converters import number_field_to_json
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: number_field_from_json(number_field_to_json(F)) == F
        True
        sage: gen = CC('-0.3122516446211537 + 1.026735437750787*I')
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z', embedding=gen)
        sage: number_field_from_json(number_field_to_json(F)) == F
        True
    """
    if isinstance(data, str):
        data = json.loads(data)
    emb = data.get("embedding", None)
    return NumberField(
        ZZ["x"](data["polynomial"]),
        names=tuple(data["names"]),
        embedding=ComplexField(53)(emb) if emb else None,
    )


# --- Complex number serialization ---


def is_json_number(data: Any) -> bool:
    r"""
    Check if the input is a JSON representation of a complex number.

    INPUT:

    - ``data`` -- any; the data to check

    OUTPUT:

    - Boolean; True if the data is a dict with exactly 'prec' and 'val' keys

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import is_json_number
        sage: is_json_number({'prec': 53, 'val': '1.0'})
        True
        sage: is_json_number({'foo': 1})
        False
        sage: is_json_number('not a dict')
        False
    """
    return isinstance(data, dict) and "prec" in data and "val" in data and len(data.keys()) == 2


def complex_number_to_json(s) -> dict:
    r"""
    Convert a ComplexNumber to a JSON-serializable dictionary.

    INPUT:

    - ``s`` -- ComplexNumber or coercible value

    OUTPUT:

    - dict; JSON representation with 'prec' and 'val' keys

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import complex_number_to_json
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CF = ComplexField(53)
        sage: d = complex_number_to_json(CF(1 + 2*I))
        sage: d['prec']
        53
        sage: 'val' in d
        True
    """
    if is_json_number(s):
        return s
    CC = ComplexField(53)
    if not isinstance(s, ComplexNumber):
        s = CC(s)
    return {"prec": s.parent().prec(), "val": str(s)}


def complex_number_from_json(json_complex: dict | str) -> ComplexNumber:
    r"""
    Create a ComplexNumber from JSON data.

    INPUT:

    - ``json_complex`` -- dict or str; a dictionary with 'prec' and 'val' keys,
      or a JSON string representing such a dictionary

    OUTPUT:

    - ComplexNumber; the reconstructed complex number

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import complex_number_from_json
        sage: complex_number_from_json({'prec': 53, 'val': '1.0 + 2.0*I'})
        1.00000000000000 + 2.00000000000000*I
        sage: complex_number_from_json('{"prec": 53, "val": "1.0"}')
        1.00000000000000
        sage: complex_number_from_json({'prec': 100, 'val': '1.0'})
        1.0000000000000000000000000000

    TESTS::

        sage: complex_number_from_json('not a json')
        Traceback (most recent call last):
        ...
        JSONDecodeError: ...
    """
    if isinstance(json_complex, str):
        json_complex = json.loads(json_complex)
    try:
        return ComplexField(json_complex["prec"])(json_complex["val"])
    except TypeError:
        return ComplexField(json_complex)


def complex_tuple_to_json(complex_tuple) -> list[dict]:
    r"""
    Convert a tuple of complex numbers to a JSON-serializable list.

    INPUT:

    - ``complex_tuple`` -- tuple or list of ComplexNumber

    OUTPUT:

    - list of dicts; JSON representation of the complex numbers

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import complex_tuple_to_json
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CF = ComplexField(53)
        sage: t = (CF(1), CF(2 + I))
        sage: j = complex_tuple_to_json(t)
        sage: len(j)
        2
        sage: j[0]['prec']
        53
    """
    return [complex_number_to_json(s) for s in complex_tuple]


def complex_tuple_from_json(json_list: list[dict] | str) -> tuple:
    r"""
    Convert a JSON list to a tuple of complex numbers.

    INPUT:

    - ``json_list`` -- list of dicts or a string; JSON representation of complex numbers

    OUTPUT:

    - Tuple of ComplexNumber objects

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import complex_tuple_to_json, complex_tuple_from_json
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CF = ComplexField(53)
        sage: t = (CF(1), CF(2))
        sage: j = complex_tuple_to_json(t)
        sage: complex_tuple_from_json(j)
        (1.00000000000000, 2.00000000000000)
        sage: complex_tuple_from_json('[{"prec": 53, "val": "1.0"}]')
        (1.00000000000000,)
    """
    if isinstance(json_list, str):
        json_list = json.loads(json_list)
    return tuple(complex_number_from_json(s) for s in json_list)


# --- Dict serialization (non-string keys) ---


def to_json_dict_key(key) -> str:
    r"""
    Convert a dict key to a JSON-compatible string key.

    INPUT:

    - ``key`` -- a dict key (string, tuple, or SageMath element)

    OUTPUT:

    - str; JSON-safe string representation of the key

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import to_json_dict_key
        sage: to_json_dict_key('hello')
        'hello'
        sage: to_json_dict_key((int(1), int(2)))
        '[1, 2]'
    """
    if isinstance(key, str):
        return key
    return json.dumps(key, cls=SageJSONEncoder)


def dict_to_json(data: dict) -> dict:
    r"""
    Transform a Python dict to a JSON-safe dict with string keys.

    INPUT:

    - ``data`` -- dict; a Python dictionary (may have non-string keys)

    OUTPUT:

    - dict; a dictionary with string keys and JSON-safe values

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import dict_to_json, dict_from_json
        sage: d = {(int(1), int(2)): int(1), int(2): int(2)}
        sage: dict_from_json(dict_to_json(d)) == d
        True
        sage: d = {(1, 2): 3, 4: 5}
        sage: dict_from_json(dict_to_json(d)) == d
        True
    """
    new_dict = {}
    for key, value in data.items():
        new_key = to_json_dict_key(key)
        if isinstance(value, dict):
            value = dict_to_json(value)
        else:
            value = SageJSONEncoder().default(value)
        new_dict[new_key] = value
    return new_dict


def dict_from_json(data: str | dict) -> dict:
    r"""
    Transform a JSON dict back to a Python dict, restoring non-string keys.

    INPUT:

    - ``data`` -- str or dict; JSON data

    OUTPUT:

    - dict; Python dictionary with restored keys

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import dict_from_json
        sage: dict_from_json('{"a": 1, "b": 2}')
        {'a': 1, 'b': 2}
        sage: dict_from_json('{"1": 1, "2": 2}')
        {1: 1, 2: 2}
        sage: dict_from_json('{"[1,2]": 1, "2": 2}')
        {(1, 2): 1, 2: 2}
    """
    if isinstance(data, str):
        data = json.loads(data)
    if "__type__" in data:
        return decode_function(data)
    new_dict = {}
    for key, value in data.items():
        try:
            new_key = json.loads(key, object_hook=decode_function)
            if isinstance(new_key, list):
                new_key = tuple(
                    [
                        json.loads(x, object_hook=decode_function) if isinstance(x, str) else x
                        for x in new_key
                    ]
                )
        except json.decoder.JSONDecodeError:
            new_key = key
        new_dict[new_key] = decode_function(value)
    return new_dict


# --- Coefficient dict serialization ---


def coefficient_dict_to_json(coeff_dict: dict) -> dict:
    r"""
    Convert a coefficient dict (tuple keys, complex values) to JSON format.

    INPUT:

    - ``coeff_dict`` -- dict; keys are integer tuples, values are complex numbers

    OUTPUT:

    - dict; JSON-serializable with stringified keys and complex number dicts as values

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import coefficient_dict_to_json
        sage: data = { (0, 0): 1, (1, 1): 1, (2, 2): 1 }
        sage: coefficient_dict_to_json(data)
        {'[0, 0]': {'prec': 53, 'val': '1.00000000000000'},
         '[1, 1]': {'prec': 53, 'val': '1.00000000000000'},
         '[2, 2]': {'prec': 53, 'val': '1.00000000000000'}}
        sage: from sage.rings.complex_mpfr import ComplexField
        sage: CF = ComplexField(103)
        sage: data = { (0, 0): 1, (1, 1): CF(1), (2, 2): CF(0,1) }
        sage: d = coefficient_dict_to_json(data)
        sage: d['[2, 2]']['prec']
        103
    """
    if not coeff_dict:
        return {}
    if all(isinstance(x, str) and is_json_number(x) for x in coeff_dict.items()):
        return coeff_dict
    return {
        json.dumps([int(ki) for ki in k]): complex_number_to_json(v) for k, v in coeff_dict.items()
    }


def coefficient_dict_from_json(data: dict | str) -> dict:
    r"""
    Convert a coefficient dict from JSON format.

    INPUT:

    - ``data`` -- dict or str; JSON data

    OUTPUT:

    - dict; coefficient dictionary with tuple keys and ComplexNumber values

    EXAMPLES::

        sage: from maass_form_core.utils.json_converters import (
        ....:     coefficient_dict_from_json, coefficient_dict_to_json)
        sage: data = { (0, 0): 1, (1, 1): 1, (2, 2): 1 }
        sage: json_data = coefficient_dict_to_json(data)
        sage: coefficient_dict_from_json(json_data)
        {(0, 0): 1.00000000000000, (1, 1): 1.00000000000000, (2, 2): 1.00000000000000}
    """
    if isinstance(data, str):
        data = json.loads(data)
    return {tuple(json.loads(k)): complex_number_from_json(v) for k, v in data.items()}
