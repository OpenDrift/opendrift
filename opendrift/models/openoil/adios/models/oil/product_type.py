
"""
Class that represents the product type

With validation

Also maintains the products types and labels mapping
"""

from pathlib import Path
import csv

from ...util.many_many import ManyMany
from ..common.validators import EnumValidator
from .validation.warnings import WARNINGS


class TypeLabelsMap(ManyMany):
    """
    class to maintain a many to many relationship between product types and labels

    The ``.product_types`` attribute is a mapping with the labels as keys,
    and product types as values.

    The ``.labels`` attribute is a mapping with the product type as keys,
    and the associated labels as values.
    """

    product_types = ManyMany.right
    labels = ManyMany.left


def load_from_csv_file(filepath=None):
    """
    Loads the product types and labels mapping from a CSV file

    :param filepath=None: The name of the file to load from. If not
                          provided, it will look for:
                          "product_types_and_labels.csv" next to this module

    """

    if filepath is None:
        filepath = Path(__file__).parent / "product_types_and_labels.csv"

    with open(filepath, newline='', encoding="utf-8") as csvfile:
        # skip the header:
        while True:
            line = csvfile.readline()
            if line.strip().startswith("DATA TABLE:"):
                break

        ptypes_labels = {}
        reader = csv.reader(csvfile, dialect='excel')
        labels = next(reader)[2:]
        for row in reader:
            pt = row[0]
            for i, val in enumerate(row[2:]):
                if val.strip():
                    ptypes_labels.setdefault(pt, set()).add(labels[i])
            # add the pt to the labels
            # not doing this anymore.
            # ptypes_labels.setdefault(pt, set()).add(pt)
        return ptypes_labels


PRODUCT_TYPE_LABEL_MAPPING = load_from_csv_file()
PRODUCT_TYPES = tuple(PRODUCT_TYPE_LABEL_MAPPING)

types_to_labels = TypeLabelsMap(PRODUCT_TYPE_LABEL_MAPPING)

# I'd much rather not hard-code this, but how else to do it?
# this is used by the validation code
DOESNT_NEED_API = set(('Refined Product NOS',
                       'Refinery Intermediate',
                       'Solvent',
                       'Natural Plant Oil',
                       'Other'))


class ProductType(str):
    _valid_types = PRODUCT_TYPES
    _validator = EnumValidator(_valid_types,
                               WARNINGS["W003"],
                               case_insensitive=True)

    @classmethod
    def validate(cls, value):
        if not value:
            return [WARNINGS["W002"]]
        return cls._validator(value)







