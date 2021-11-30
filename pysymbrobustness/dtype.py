import math
from fractions import Fraction
from typing import Union, List

from pysymbrobustness.ta import exceptions as exceptions

Delay = Union[int, float, Fraction]  # float to allow math.inf
Valuation = List[Delay]
Reset = List[int]
Action = str
Sampling_Step = Union[int, Fraction]
Location = Union[int, str]


def check_delay_type(delay: Delay) -> None:
    """
    check if delay is compatible with the Interval class: i.e it should be integer, fraction, or math.inf
    :param delay:
    :return:
    """
    if type(delay) is not Fraction and type(delay) is not int and not math.isinf(delay):
        raise exceptions.WrongType(element="delay=" + str(delay),
                                   right_type="integer or Fraction, or math.inf")