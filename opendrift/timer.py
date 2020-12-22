from collections import OrderedDict
from datetime import datetime, timedelta

class Timeable:
    """
    Utility class for measuring total time spent in various steps in a class
    throughout program execution.
    """
    __timers__ = None
    __timing__ = None

    @property
    def timers(self):
        if self.__timers__ is None:
            self.__timers__ = OrderedDict()

        return self.__timers__

    @property
    def timing(self):
        if self.__timing__ is None:
            self.__timing__ = OrderedDict()

        return self.__timing__

    def timer_start(self, category):
        if category not in self.timing:
            self.timing[category] = timedelta(0)
        self.timers[category] = datetime.now()

    def timer_end(self, category):
        if self.timers[category] is not None:
            self.timing[category] += datetime.now() - self.timers[category]
        self.timers[category] = None

