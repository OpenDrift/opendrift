"""
class to handle the version number
"""
from dataclasses import dataclass, field
import functools

from ..common.utilities import dataclass_to_json


class VersionError(ValueError):
    pass


# Note: not using @dataclass_to_json -- defining py_json and from_pyjson instead
@functools.total_ordering
# @dataclass
class Version:
    """
    class to represent a version with major, minor, and patch components
    """
    def __init__(self, major, minor=0, patch=0):
        """
        initialize and instance object

        """
        if isinstance(major, str) and minor == 0 and patch == 0:
            self.__init__(*self._parts_from_string(major))
        else:
            self.major = int(major)
            self.minor = int(minor)
            self.patch = int(patch)

    def __repr__(self):
        return f"Version({self.major}, {self.minor}, {self.patch})"

    def __str__(self):
        return f"{self.major}.{self.minor}.{self.patch}"

    def py_json(self, sparse=True):
        # return f"{self.major}.{self.minor}.{self.patch}"
        return str(self)

    @staticmethod
    def _parts_from_string(string):
        try:
            return tuple(int(part) for part in string.split("."))
        except Exception:  # anything goes wrong, we'll raise a new error
            raise ValueError(f"{string} not a valid value for Version string")

    @classmethod
    def from_py_json(cls, py_json, allow_none=False):
        return cls(*cls._parts_from_string(py_json))

    def __eq__(self, other):
        return (self.major == other.major
                and self.minor == other.minor
                and self.patch == other.patch)

    def __gt__(self, other):
        return ((self.major, self.minor, self.patch)
                > (other.major, other.minor, other.patch))
