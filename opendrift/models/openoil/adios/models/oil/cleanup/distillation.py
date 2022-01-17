"""
cleanups that work with distillation data

There could be more, but for now, it sets the unit_type for the cuts.

"""
import math
from operator import itemgetter

import unit_conversion as uc

from .cleanup import Cleanup

from ....computation.physical_properties import Density


class FixCutUnitType(Cleanup):

    """
    sets the unit type for the cut fractions, based on the distillation type
    """

    ID = "002"  # maybe a better way to keep track of these?

    mapping = {None: 'concentration',
               'mass fraction': 'massfraction',
               'volume fraction': 'volumefraction',
               }


    def check(self):
        """
        checks to see if there is something to fix

        returns: flag, msg

        if nothing is needed, flag is None
        if something can be cleaned up, flag is True
        if something is wrong, but can not be cleaned up, flag is False

        fixme: -- maybe cleanup and validation should be better integrated?
        """
        needs_fixing = False
        for ss in self.oil.sub_samples:
            dist_data = ss.distillation_data
            for cut in dist_data.cuts:
                if cut.fraction.unit_type is not self.mapping[dist_data.type]:
                    needs_fixing = True
                    break
        if needs_fixing:
            return (True, f"Cleanup: {self.ID}: "
                    "distillation cuts have wrong unit_type -- can be corrected")
        else:
            return None, "distillation cuts have correct unit_type"


    def cleanup(self):
        """
        run this particular cleanup option
        """
        for ss in self.oil.sub_samples:
            dist_data = ss.distillation_data
            for cut in dist_data.cuts:
                cut.fraction.unit_type = self.mapping[dist_data.type]


