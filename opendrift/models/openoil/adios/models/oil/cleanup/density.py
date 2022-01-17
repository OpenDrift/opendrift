"""
cleanups that work with density
"""
import math
from operator import itemgetter

import unit_conversion as uc

from .cleanup import Cleanup

from ....computation.physical_properties import Density

class FixAPI(Cleanup):
    """
    adds (or replaces) the API value, from the density measurements

    NOTE: this could be extended to interpolate, but it that actually needed?
          There is code in the computation.physical_properties package to help, if needed.
    """
    ID = "001"

    def check(self):
        """
        checks to see if there is something to fix

        returns: flag, msg

        if nothing is needed, flag is None
        if something can be cleaned up, flag is True
        if something is wrong, but can not be cleaned up, flag is False

        fixme: -- maybe cleanup and validation should be better integrated?
        """
        API = self.oil.metadata.API

        # densities = oil.sub_samples[0].physical_properties.densities
        if API is None:

            density = self.find_density_near_15C()
            if density:
                return (True, f"Cleanup: {self.ID}: No API value provided for {self.oil.oil_id}"
                               " -- can be computed from density")
            else:
                return (False, f"Cleanup: {self.ID}: No API value provided for {self.oil.oil_id}"
                                " -- can NOT be computed from density")

        return None, "API is fine"

    def cleanup(self):
        """
        run this particular cleanup option

        :param oil: an Oil object to act on

        :param do_it=False: flag to tell the cleanup to do its thing. If False,
                            the method returns a message. If True, the action is
                            taken, and the Oil object is altered.

        :returns: a message of what could be done, or what was done.
        """
        density_at_15 = self.find_density_near_15C()

        if density_at_15:
            API = uc.convert("density", "kg/m^3", "API", density_at_15)
            self.oil.metadata.API = round(API, 2)
            return f"Cleanup: {self.ID}: Set API for {self.oil.oil_id} to {API}."

    def check_for_valid_api(self):
        """
        check is the API value is already valid
        """
        API = self.oil.metadata.API

        density_at_15 = self.find_density_near_15C()

        if uc.convert("density", "kg/m^3", "API", density_at_15) == API:
            return True
        else:
            return False

    def find_density_near_15C(self):
        """
        Returns the density (in kg/m3)

        It will interpolate and extrapolate as needed
        """
        try:
            return Density(self.oil).at_temp(uc.convert("C", "K", 15))
        except ValueError:
            return None
