import json
from json import JSONDecodeError
from typing import Any

from sage.all import ZZ
from sage.matrix.constructor import matrix
from sage.rings.complex_mpfr import ComplexField, ComplexField_class
from sage.rings.number_field.number_field import NumberField, NumberField_generic
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.rational import Rational
from sage.rings.rational_field import RationalField
from sage.rings.integer_ring import IntegerRing_class, IntegerRing
from sage.rings.real_mpfr import RealField, RealField_class

from maass_forms_klein.modform.utils import Integer_t, Real_t, get_prec, Complex_t
from sage.structure.element import Matrix, Vector


class SageJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, str):
            return o
        if isinstance(o, (IntegerRing_class, RealField_class, ComplexField_class,
                          RationalField, NumberField_generic)):
            return ring_to_json(o)
        if isinstance(o, (float, int, complex)):
            return o
        if isinstance(o, (Integer_t, Real_t, Complex_t, Rational)):
            return ring_element_to_json(o)
        if isinstance(o, Matrix):
            return matrix_to_json(o)
        if isinstance(o, (Vector, tuple)):
            return list(o)
        if isinstance(o, dict):
            return dict_to_json(o)
        return super(self, SageJSONEncoder).default(o)


def decode_function(obj: dict | str) -> Any:
    r"""
    Decode a function from JSON or dict.

    EXAMPLES::

        sage: from maass_forms_klein.utils.json_converters import decode_function
        sage: decode_function({'__type__': 'matrix',
        ....:                   'base_ring': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
        ....:                    'entries': [['1', '2'], ['3', '4']]})
        [1 2]
        [3 4]
        sage: decode_function('[[1, 2], [3, 4]]')
        '[[1, 2], [3, 4]]'
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
            print(k, v)
            try:
                k_new = json.loads(k, object_hook=decode_function)
            except JSONDecodeError:
                k_new = k
            try:
                v_new = json.loads(v, object_hook=decode_function)
            except TypeError:
                v_new = v
            if isinstance(k_new, list):
                k_new = tuple(k_new)
            new_dict[k_new] = v_new
    if isinstance(obj, list):
        return tuple([json.loads(x, object_hook=decode_function) for x in obj])
    return obj


def matrix_to_json(m):
    """
    JSON representation of a matrix.

    EXAMPLES::
        sage: from maass_forms_klein.utils.json_converters import matrix_to_json
        sage: m = matrix(ZZ, [[1, 2], [3, 4]])
        sage: matrix_to_json(m)
        {'__type__': 'matrix',
         'base_ring': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
         'entries': [['1', '2'], ['3', '4']]}
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: z = F.gen()
        sage: m = matrix(F, [[z^2, z^3], [z^4, z^5]])
        sage: matrix_to_json(m)
        {'__type__': 'matrix',
         'base_ring': {'__type__': 'ring',
          'field': {'__type__': 'ring',
           'embedding': None,
           'names': ['z'],
           'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2'},
          'name': 'NumberField'},
         'entries': [['z^2', 'z^3'], ['z^4', 'z^5']]}
    """
    return {
        "base_ring": ring_to_json(m.base_ring()),
        "entries": [[str(x) for x in row] for row in list(m)],
        "__type__": "matrix"
    }

def ring_element_to_json(data: Any) ->dict:
    """
    JSON representation of a ring element.

    INPUT:

    - ``data`` -- a ring element

    EXAMPLES::

        sage: from maass_forms_klein.utils.json_converters import ring_element_to_json
        sage: ring_element_to_json(1)
        {'__type__': 'element',
         'parent': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
         'value': '1'}
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: ring_element_to_json(F.gen())
        {'__type__': 'element',
         'parent': {'__type__': 'ring',
          'embedding': None,
          'names': ['z'],
          'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2'},
         'value': 'z'}
    """
    base_ring = None
    # print("data", data, type(data))
    if isinstance(data, Integer_t):
        base_ring = {"name": "IntegerRing", "prec": get_prec(data), "__type__": "ring"}
    if isinstance(data, Real_t):
        base_ring = {"name": "RealField", "prec": get_prec(data), "__type__": "ring"}
    if isinstance(data, Complex_t):
        base_ring = {"name": "ComplexField", "prec": get_prec(data), "__type__": "ring"}
    if isinstance(data, Rational):
        base_ring = {"name": "RationalField", "prec": get_prec(data), "__type__": "ring"}
    if isinstance(data, NumberFieldElement):
        base_ring = number_field_to_json(data.parent())
    return {
        "parent": base_ring,
        "value": str(data),
        "__type__": "element"
    }


def ring_element_from_json(data: dict) -> Any:
    """
    Convert JSON representation of a ring element to a sage object.

    INPUT:

    - ``data`` -- a JSON object
    EXAMPLES::

        sage: from maass_forms_klein.utils.json_converters import ring_element_from_json
        sage: ring_element_from_json({'parent': {'name': 'IntegerRing', 'prec': 0}, 'value': '1'})
        1
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


def matrix_from_json(data):
    """
    Convert JSON representation of a matrix to a sage matrix.

    EXAMPLES:
        sage: from maass_forms_klein.utils.json_converters import matrix_from_json
        sage: matrix_from_json({'base_ring': {'name': 'IntegerRing', 'prec': 0}, 'entries': [['1', '2'], ['3', '4']]})
        [1 2]
        [3 4]
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: z = F.gen()
        sage: matrix_from_json({'base_ring': {'field': {'names': ['z'],
        ....: 'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2'},
        ....: 'name': 'NumberField'},
        ....: 'entries': [['z^2', 'z^3'], ['z^4', 'z^5']]})
        [z^2 z^3]
        [z^4 z^5]
    """
    return matrix(ring_from_json(data["base_ring"]), data["entries"])


def ring_to_json(F):
    """
    JSON representation of a base ring.

    EXAMPLES:
        sage: from maass_forms_klein.utils.json_converters import ring_to_json
        sage: ring_to_json(RationalField())
        {'__type__': 'ring', 'name': 'RationalField', 'prec': 0}
        sage: ring_to_json(matrix(QQ, [[1, 2], [3, 4]]).base_ring())
        {'__type__': 'ring', 'name': 'RationalField', 'prec': 0}
        sage: ring_to_json(RealField(53))
        {'__type__': 'ring', 'name': 'RealField', 'prec': 53}
        sage: ring_to_json(ComplexField(53))
        {'__type__': 'ring', 'name': 'ComplexField', 'prec': 53}
        sage: ring_to_json(ZZ)
        {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0}
        sage: ring_to_json(matrix(ZZ, [[1, 2], [3, 4]]).base_ring())
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
        return {"name": "RationalField", "prec": (0), "__type__": "ring"}
    if isinstance(F, IntegerRing_class):
        return {"name": "IntegerRing", "prec": (0), "__type__": "ring"}
    if isinstance(F, NumberField_generic):
        return {"name": "NumberField", "field": number_field_to_json(F), "__type__": "ring"}
    raise ValueError(f"Unsupported base ring {F}")


def ring_from_json(data: dict | str) -> (RealField_class | ComplexField_class | RationalField |
                                         NumberField_generic | IntegerRing_class):
    """
    Construct ring from json data.

    INPUT:

    - ``data``: JSON data as string or dict.

    EXAMPLES:
        sage: from maass_forms_klein.utils.json_converters import ring_from_json, number_field_to_json
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


def number_field_to_json(nf: NumberField_generic) -> dict:
    """
    JSON representation of number field.

    NOTE: Any information about embeddings is ignored.

    EXAMPLES:

        sage: from maass_forms_klein.utils.json_converters import number_field_to_json
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z')
        sage: number_field_to_json(F)
        {'__type__': 'ring',
         'embedding': None,
         'names': ['z'],
         'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2'}
        sage: gen = CC('-0.3122516446211537 + 1.026735437750787*I')
        sage: F = NumberField(x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2, 'z', embedding=gen)
        sage: number_field_to_json(F)
        {'__type__': 'ring',
         'embedding': '-0.312251644621154 + 1.02673543775079*I',
         'names': ['z'],
         'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2'}

    """
    return {"polynomial": str(nf.polynomial()), "names": list(nf._names), "__type__": "ring",
            "embedding": str(nf.gen_embedding().n(53)) if nf.gen_embedding() else None}


def number_field_from_json(data: dict | str) -> NumberField_generic:
    """
    Create number field from json data.

    INPUT:

    - ``data``: JSON data as string or dict.

    EXAMPLES:
        sage: from maass_forms_klein.utils.json_converters import number_field_from_json
        sage: number_field_from_json({'polynomial': 'x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2', 'names': ['a']})
        Number Field in a with defining polynomial x^8 + 2*x^6 + 3*x^4 + 3*x^2 + 2
        sage: from maass_forms_klein.utils.json_converters import number_field_to_json
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

    return NumberField(ZZ["x"](data["polynomial"]), names=tuple(data["names"]),
                       embedding=ComplexField(53)(emb) if emb else None)


def dict_from_json(data: str | dict) -> dict:
    """
    Transform a json dict to a python dict, which can have non-str keys.

    INPUT:

    - ``data``: JSON data as string or dict.

    EXAMPLES::

        sage: from maass_forms_klein.utils.json_converters import dict_from_json
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
                new_key = tuple([json.loads(x, object_hook=decode_function)
                                 if isinstance(x, str) else x
                                 for x in new_key])
        except json.decoder.JSONDecodeError:
            new_key = key
        # if isinstance(value, dict):
        #     value = dict_from_json(value)
        new_dict[new_key] = decode_function(value)
    return new_dict

def to_json_dict_key(key):
    if isinstance(key, str):
        return key
    return json.dumps(key, cls=SageJSONEncoder)

def dict_to_json(data: dict) -> str:
    """
    Transform a python dict to a json dict, which can only have str keys.

    INPUT:

    - ``data``: Python dict

    EXAMPLES::

        sage: from maass_forms_klein.utils.json_converters import dict_to_json
        sage: dict_to_json({'a': 1, 'b': 2})
        {'a': {'__type__': 'element',
          'parent': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
          'value': '1'},
         'b': {'__type__': 'element',
          'parent': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
          'value': '2'}}
        sage: dict_to_json({int(1): int(1), 2: int(2)})
        {'1': 1,
         '{"parent": {"name": "IntegerRing", "prec": 0,  "__type__": "ring"}, "value": "2", ...
        sage: dict_to_json({(int(1), int(2)): int(1), int(2): 'a'})
        {'2': 'a', '[1, 2]': 1}
        sage: dict_to_json({(1, 2): 1, 2: 2})
        {'[{"parent": {"name": "IntegerRing", "prec": 0, "__type__": "ring"}, "value": "1", ...
              'parent': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
              'value': '1'},
             '{"parent": {"name": "IntegerRing", "prec": 0, "__type__": "ring"}, "value": "2", ...
              'parent': {'__type__': 'ring', 'name': 'IntegerRing', 'prec': 0},
              'value': '2'}}
        sage: from maass_forms_klein.utils.json_converters import dict_from_json
        sage: d = {(int(1), int(2)): int(1), int(2): int(2)}
        sage: dict_from_json(dict_to_json(d)) == d
        True
        sage: d = {(1, 2): 3, 4: 5}
        sage: dict_from_json(dict_to_json(d))
        {(1, 2): 3, 4: 5}
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
        # elif isinstance(value, Integer_t):
        #     value = int(value)
        # elif isinstance(value, Real_t):
        #     value = float(value)
        # print("new key=", new_key, "value=", value)
        # try:
        new_dict[new_key] = value
        # except TypeError as e:
        #     print(e)
        #     new_dict[str(new_key)] = value
        # print("new_dict", new_dict)
    return new_dict
