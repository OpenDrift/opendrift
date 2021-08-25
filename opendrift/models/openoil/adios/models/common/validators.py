

class EnumValidator:
    """
    validator for Enum values: a value that can only be one of a set

    """
    def __init__(self, valid_items, err_msg, case_insensitive=False):
        """
        :param valid_items: list of valid items -- can be anything that an `in`
                            test will work for.
        :param err_msg: The error message that should be used on failure.
                        Should be a format string that takes two parameters:
                        item and valid_items

        :param case_insensitive=False: whether you want the test to be
                                       case-insensitive.
                                       Only works for string values, of course.
        """
        if case_insensitive:
            valid_items = [item.lower() for item in valid_items]

        self.valid_items = valid_items
        self.err_msg = err_msg
        self.case_insensitive = case_insensitive

    def __call__(self, item):
        if self.case_insensitive:
            try:
                item = item.lower()
            except AttributeError:
                pass  # so non-strings will then fail the in test, but not crash.

        if item not in self.valid_items:
            return [self.err_msg.format(item, self.valid_items)]
        else:
            return []


class FloatRangeValidator:
    """
    Validator for float values that can only be a given range

    range is inclusive (<= and >=)
    """
    def __init__(self, min_value, max_value, err_msg=None):
        """
        :param min: minimum value allowed

        :param max: maximum value allowed

        :param err_msg: The error message that should be used on failure.
                        Should be a format string that takes three parameters:
                        default is:
                            "ValidationError: {} is not between {} and {}"
        """
        self.min = min_value
        self.max = max_value

        if err_msg is None:
            self.err_msg = "ValidationError: {} is not between {} and {}"
        else:
            self.err_msg = err_msg

    def __call__(self, value):
        try:
            value = float(value)
        except ValueError:
            return [self.err_msg.format(value, self.min, self.max)]

        if not self.min <= value <= self.max:
            return [self.err_msg.format(value, self.min, self.max)]
        else:
            return []

class YearValidator:
    """
    Validator for float values that can only be a given range

    range is inclusive (<= and >=)
    """
    def __init__(self, min_year, max_year, err_msg=None):
        """
        :param min: minimum year allowed

        :param max: maximum year allowed

        :param err_msg: The error message that should be used on failure.
                        Should be a format string that takes three parameters:
                        default is:
                            "ValidationError: {} is not between {} and {}"
        """
        self.min = min_year
        self.max = max_year

        if err_msg is None:
            self.err_msg = "ValidationError: {} is not a year between {} and {}"
        else:
            self.err_msg = err_msg

    def __call__(self, value):
        try:
            value = int(value)
        except ValueError:
            return [self.err_msg.format(value, self.min, self.max)]

        if not self.min <= value <= self.max:
            return [self.err_msg.format(value, self.min, self.max)]
        else:
            return []




