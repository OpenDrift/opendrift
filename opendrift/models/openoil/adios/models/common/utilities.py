"""
Tools for helping make our data models.

So far: making dataclasses read/writable as JSON
"""


def something(val):
    '''
        much like python's "Truthy" and Falsey", but we want some values
        to not be false like zero, for instance
    '''
    return ((val == 0) or (val is not None) and val)


def dataclass_to_json(cls):
    """
    class decorator that adds the ability to save a dataclass as JSON

    All fields must be either JSON-able Python types or
    have be a type with a _to_json method
    """

    @classmethod
    def from_py_json(cls, py_json, allow_none=False):
        """
        classmethod to create a dataclass from json compatible python data
        structure.
        """
        # if there is a pre-processor - run it
        if hasattr(cls, "_pre_from_py_json"):
            py_json = cls._pre_from_py_json(py_json)

        arg_dict = {}

        if py_json is None and allow_none is True:
            # the parent object defined an attribute with a default of None
            # We could actually allow other default types, but this one is
            # common
            return py_json

        for fieldname, fieldobj in cls.__dataclass_fields__.items():
            if fieldname in py_json:
                allow_none = True if fieldobj.default is None else False

                try:  # see if it's "one of ours"
                    arg_dict[fieldname] = (fieldobj.type
                                           .from_py_json(py_json[fieldname],
                                                         allow_none=allow_none)
                                           )
                except AttributeError:
                    # it's not, so we just use the value
                    arg_dict[fieldname] = py_json[fieldname]
                except TypeError as err:
                    raise TypeError(f'TypeError in '
                                    f'{cls.__name__}._from_py_json(): '
                                    f'field: {fieldname}\n'
                                    f'{err.args}')

        obj = cls(**arg_dict)
        return obj

    def py_json(self, sparse=True):
        """
        function to convert a dataclass to json compatible python

        :param sparse=True: If sparse is True, only non-empty fields will be
                            written. If False, then all fields will be
                            included.
        """
        json_obj = {}
        for fieldname in self.__dataclass_fields__.keys():
            val = getattr(self, fieldname)
            try:  # convert to json
                val = val.py_json(sparse=sparse)
            except AttributeError:
                pass

            if not sparse:
                json_obj[fieldname] = val
            elif something(val):
                json_obj[fieldname] = val

        return json_obj

    def validate(self):
        """
        Function to validate a dataclass with fields that have validate
        methods.  The validate methods are expected to return a list of
        validation messages.

        The top-level validator extends the existing list
        """
        # I have NO idea how this happens, but it does ?!?!?
        if self is None:
            return []

        # first see if there is a "private" one:
        if hasattr(self, '_validate'):
            messages = self._validate()
        else:
            messages = []

        for fieldname, fieldobj in self.__dataclass_fields__.items():
            value = getattr(self, fieldname)
            if hasattr(fieldobj, 'type') and hasattr(fieldobj.type, 'validate'):
                messages.extend(fieldobj.type.validate(value))
            # try:
            #     # validate with the type's validate method
            #     messages.extend(fieldobj.type.validate(value))
            # except AttributeError as err:  # This one doesn't have a validate method.
            #     print("\nAttributeError:", err)
            #     pass

        return messages

    def __setattr__(self, name, val):
        try:
            _fieldobj = self.__dataclass_fields__[name]
        except KeyError:
            raise AttributeError(f"You can only set existing attributes: "
                                 f"{self.__class__.__name__}.{name} "
                                 "does not exist")
        self.__dict__[name] = val

    def __repr__(self):
        atts = ((att, getattr(self, att)) for att in self.__dataclass_fields__.keys())
        atts = (f'{att}={val!r}' for att, val in atts if val is not None)
        return f'{self.__class__.__name__}({", ".join(atts)})'

    cls.py_json = py_json

    cls.from_py_json = from_py_json

    if hasattr(cls, "validate"):
        cls._validate = cls.validate
    cls.validate = validate

    cls.__setattr__ = __setattr__
    cls.__repr__ = __repr__

    return cls


class JSON_List(list):
    """
    just like a list, but with the ability to turn it into JSON

    A regular list can only be converted to JSON if it has
    JSON-able objects in it.

    Note: must be subclassed, and the item_type attribute set
    """
    item_type = None

    def py_json(self, sparse=True):
        json_obj = []
        for item in self:
            try:
                json_obj.append(item.py_json(sparse))
            except AttributeError:
                json_obj.append(item)

        return json_obj

    @classmethod
    def from_py_json(cls, py_json, allow_none=False):
        """
        create a JSON_List from json array of objects
        that may be json-able.
        """
        if cls.item_type is None:
            raise TypeError("You can not reconstruct a list of unknown type")

        jl = cls()  # an empty JSON_List

        # loop through contents
        for item in py_json:
            jl.append(cls.item_type.from_py_json(item))

        return jl

    def __repr__(self):
        return f"{self.__class__.__name__}({list.__repr__(self)})"
