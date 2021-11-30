# coding=utf-8

# Guards exceptions

class AbstractConstructionException(Exception):
    def __init__(self, *args):
        super().__init__(*args)


class ConstraintNotFoundException(Exception):
    pass  # TODO: Fill the exception


class GuardNotFoundException(Exception):
    """
    Raised when a guards set is empty
    """
    pass


class ClockNotFoundException(Exception):
    """
    Raised when the clocks number is under 1.
    """
    pass


class IndexOutOfRange(Exception):
    """
    Raised when a clock or reset index is greater than the number of clocks
    of the timed Automaton.
    """

    def __init__(self, guard, type_object):
        if type_object == "clock":
            super().__init__("Clock index in guard: " + str(guard) +
                             "out of range")
        elif type_object == "reset":
            super().__init__("Reset index: " + str(guard) +
                             "out of range")


class LocationNotFoundException(Exception):
    """
    Raised when a location is not in the locations set.
    """

    def __init__(self, loc1, loc2):
        if loc1 is None:
            super().__init__("Final location: " + str(loc2) + " Not found")
        elif loc2 is None:
            super().__init__("Starting location: " + str(loc1) + " Not found")
        else:
            super().__init__("Location: " + str(loc1) + " or location: " +
                             str(loc2) + " Not found")


class NegativeIntervalException(Exception):
    pass


class NegativeClockIndexException(Exception):
    pass


class FutureLocationNotFound(Exception):
    pass


class UnknownGuardType(Exception):
    pass


class ConstraintTypeFail(Exception):
    pass


class WrongType(Exception):
    def __init__(self, right_type, element):
        super().__init__(element + " should be " + right_type)


class IllegalClosedArgument(Exception):
    def __init__(self):
        super().__init__("Illegal closed argument")


class NotMergeableOrDisjointIntervals(Exception):
    def __init__(self, int_1, int_2):
        super().__init__("The two intervals " + int_1 + " and " + int_2 + " are not mergeables or disjoints")


class NotIncludedInterval(Exception):
    def __init__(self, int_1, int_2):
        super().__init__("Interval " + int_1 + " is not included in interval" + int_2)


class EmptyInterval(Exception):
    pass


class BoundException(Exception):
    def __init__(self, lower_bound, upper_bound):
        super().__init__("left side " + lower_bound + " of interval must be <= right side" + upper_bound)
