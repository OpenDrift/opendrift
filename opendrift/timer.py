from collections import OrderedDict
from datetime import datetime, timedelta

class Timeable:
    """
    Utility class for measuring total time spent in various steps in a class
    throughout program execution.
    """
    timers = OrderedDict()
    timing = OrderedDict()

    def timer_start(self, category):
        if category not in self.timing:
            self.timing[category] = timedelta(0)
        self.timers[category] = datetime.now()

    def timer_end(self, category):
        if self.timers[category] is not None:
            self.timing[category] += datetime.now() - self.timers[category]
        self.timers[category] = None

