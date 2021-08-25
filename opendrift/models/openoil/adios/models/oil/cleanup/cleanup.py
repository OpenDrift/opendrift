#!/usr/bin/env python

"""
classes for options to clean up Records
"""

class Cleanup:

    # the ID is so that we can know which cleanup method this is
    # subclasses should define this!
    ID = None

    def __init__(self, oil):
        """
        initialize a Cleanup object with an oil object

        :param oil: the oil object you want to clean up
        :type oil: adios_db.models.oil.oil.OIl
        """
        self.oil = oil

    def cleanup(self, oil, do_it=False):
        """
        run this particular cleanup option

        :param oil: an Oil object to act on

        :param do_it=False: flag to tell the cleanup to do its thing. If False,
                            the method returns a message. If True, the action is
                            taken, and the Oil object is altered.

        :returns: a message of what could be done, or what was done.
        """
        raise NotImplementedError
