"""
cleanup modules

all cleanup classes are registered here

so far:
 - compute density from API

possible things to cleanup:

* normalizing units for various things

* putting distillation cuts in order

"""

# you need to import a cleanup class here to get it registered

from .cleanup import Cleanup

from .density import FixAPI

ALL_CLEANUPS = []
CLEANUP_MAPPING = {}  # mapping to get cleanup object from its ID

for name, obj in dict(locals()).items():
    # print("in init, obj:")
    # print(obj)
    if (not name.startswith("_")
        and isinstance(obj, type)
        and issubclass(obj, Cleanup)
        and (obj is not Cleanup)):

        ALL_CLEANUPS.append(obj)
        CLEANUP_MAPPING[obj.ID] = obj






