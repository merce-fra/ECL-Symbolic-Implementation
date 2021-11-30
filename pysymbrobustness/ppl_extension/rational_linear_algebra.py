# coding=utf-8
"""
==================================================
rational_linear_algebra module
==================================================
This module provide  an extension from the Linear_Expression class from pplpy, such that linear expression with
rational coefficients can be represented

Classes:
------
ExtendedFraction
SymetricPolynomial: Union[pplLinearExpression, math.inf, -math.inf) (float that is math.isinf)
LinearFunction: already made by ppl.LinearExpression
1dPiecewieseLinearFunction: a couple of interval and a ppl.LinearExpression or float
PiecewieseLinearFunction

Methods:
------
"""

from __future__ import annotations  # For forward reference typing

import math
import numpy as np
from typing import List, Union, Tuple
import ppl as ppl
from fractions import Fraction
from gmpy2 import mpz
import pysymbrobustness.ppl_extension.linear_algebra as linear_alg


def fraction_as_lcm_tuple_list(
        fraction_list: List[Union[int, Fraction, mpz, ...]]) -> Tuple[
    List[int], int]:
    for elt in fraction_list:
        if type(elt) != int and type(elt) != Fraction:
            raise TypeError(str(elt) + " should be Rational instance")
    dens = [int(Fraction(elt).denominator) for elt in fraction_list]
    lcm = np.lcm.reduce(dens)
    return [int(elt * lcm) for elt in fraction_list], lcm


class AbstractExpression(object):
    pass


class RationalLinearExpression(AbstractExpression):
    __slots__ = ["linear_expression", "lcm"]

    def __init__(self,
                 homogeneous_terms: Tuple[Union[int, Fraction], ...] = None,
                 inhomogeneous_term: Union[int, Fraction] = 0, div: int = 1):

        if not homogeneous_terms:
            # if type(inhomogeneous_term) != int and type(inhomogeneous_term) != Fraction:
            #     raise Exception("WrongType")
            constant = Fraction(inhomogeneous_term, div)
            self.linear_expression = ppl.Linear_Expression(constant.numerator)
            self.lcm = constant.denominator
        else:
            terms_list = [Fraction(elt, div) for elt in homogeneous_terms]
            terms_list.append(Fraction(inhomogeneous_term, div))
            # test if inhomogeneous_term = None ou 0
            int_list, self.lcm = fraction_as_lcm_tuple_list(terms_list)
            constant = int_list[len(int_list) - 1]  # verifier si il n'existe pas une pythonerie (l[-1] fonctionne)
            int_list.pop()

            self.linear_expression = ppl.Linear_Expression(int_list, constant)

    def is_equal_to(self, other: RationalLinearExpression):
        if type(other) != RationalLinearExpression:
            return False
        return self.lcm == other.lcm and linear_alg.equality(self.linear_expression, other.linear_expression)

    @staticmethod
    def linear_expression_as_rational_linear_expression(
            le: ppl.Linear_Expression):
        return RationalLinearExpression(
            homogeneous_terms=le.coefficients(),
            inhomogeneous_term=le.inhomogeneous_term())

    @staticmethod
    def constant_as_rational_linear_expression(
            constant: Union[Fraction, int, mpz]):
        return RationalLinearExpression(inhomogeneous_term=constant)

    def integer_coefficients(self):
        return self.linear_expression.coefficients()

    def coefficients(self) -> Tuple[Fraction]:
        return tuple(
            Fraction(elt, self.lcm) for elt in self.integer_coefficients()
        )

    def integer_coefficient(self, v: ppl.Variable) -> int:
        return self.linear_expression.coefficient(v)

    def coefficient(self, v: ppl.Variable) -> Fraction:
        return Fraction(self.linear_expression.coefficient(v), self.lcm)

    def inhomogeneous_term(self) -> Fraction:
        return Fraction(self.linear_expression.inhomogeneous_term(), self.lcm)

    def integer_inhomogeneous_term(self):
        return self.linear_expression.inhomogeneous_term()

    def __add__(self,
                other: Union[RationalLinearExpression, ppl.Linear_Expression,
                             int, Fraction, float]) \
            -> Union[RationalLinearExpression, float]:

        other_ty = type(other)
        if other_ty == float and not math.isinf(other):
            raise NotImplementedError
        elif (other_ty == float and math.isinf(other)) \
                or other_ty == InfiniteExpression:
            return other
        elif other_ty == int or other_ty == Fraction:
            return self.__add__(RationalLinearExpression(inhomogeneous_term=other))
        elif other_ty == ppl.Linear_Expression:
            return self.__add__(RationalLinearExpression(
                homogeneous_terms=other.coefficients(),
                inhomogeneous_term=other.inhomogeneous_term(),
                div=1
            ))
        elif other_ty == RationalLinearExpression:
            self_int_coeff = self.integer_coefficients()
            self_lcm = self.lcm
            other_int_coeff = other.integer_coefficients()
            other_lcm = other.lcm

            self_len = len(self_int_coeff)
            other_len = len(other_int_coeff)

            max_len = max(self_len, other_len)
            if self_len < max_len:
                self_int_coeff += tuple(0 for i in range(self_len, other_len))
            elif other_len < max_len:
                other_int_coeff += tuple(0 for i in range(other_len, self_len))

            new_homogeneous_terms = tuple(
                self_int_coeff[i] * other_lcm + other_int_coeff[i] * self_lcm
                for i in range(max_len)
            )

            new_div = self_lcm * other_lcm
            new_inhomogeneous_term = self.integer_inhomogeneous_term() * other_lcm + \
                                     other.integer_inhomogeneous_term() * self_lcm

            return RationalLinearExpression(
                homogeneous_terms=new_homogeneous_terms,
                inhomogeneous_term=new_inhomogeneous_term,
                div=new_div
            )

        else:
            raise NotImplementedError

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return RationalLinearExpression(
            homogeneous_terms=self.integer_coefficients(),
            inhomogeneous_term=self.integer_inhomogeneous_term(),
            div=-self.lcm
        )

    def __sub__(self, other: Union[RationalLinearExpression,
                                   ppl.Linear_Expression, int, Fraction, float]) \
            -> Union[RationalLinearExpression, float]:
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other: Union[int, Fraction, float]) -> Union[RationalLinearExpression, float]:
        if type(other) == int or type(other) == Fraction:
            new_homogeneous_terms = tuple(int(homo_term) * other for homo_term in self.integer_coefficients())
            return RationalLinearExpression(
                homogeneous_terms=new_homogeneous_terms,
                inhomogeneous_term=int(self.integer_inhomogeneous_term()) * other,
                div=self.lcm
            )
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        return self.__mul__(other)

    def __gt__(self, other):
        f = self - other
        f_linear = type(f) == RationalLinearExpression
        if f_linear:
            return f.linear_expression > 0 * ppl.Variable(0)
        elif type(f) == InfiniteExpression:
            return f > 0
        # return self.__sub__(other).linear_expression > 0 * ppl.Variable(0)

    def __ge__(self, other):
        f = self - other
        f_linear =  type(f) == RationalLinearExpression
        if f_linear:
            return f.linear_expression >= 0 * ppl.Variable(0)
        elif type(f) == InfiniteExpression:
            return f >= 0

        # return self.__sub__(other).linear_expression >= 0 * ppl.Variable(0)

    def __lt__(self, other):
        f = self.__sub__(other)
        f_linear = type(f) == RationalLinearExpression
        if f_linear:
            return self.__sub__(other).linear_expression < 0 * ppl.Variable(0)
        elif type(f) == InfiniteExpression:
            return f < 0

    def __le__(self, other):

        f = self.__sub__(other)
        f_linear = type(f) == RationalLinearExpression
        if f_linear:
            return self.__sub__(other).linear_expression <= 0 * ppl.Variable(0)
        elif type(f) == InfiniteExpression:
            return f <= 0

    def __eq__(self, other):
        f = self.__sub__(other)
        f_linear = type(f) == RationalLinearExpression
        if f_linear:
            return f.linear_expression == 0
        elif type(f) == InfiniteExpression:
            return f == 0


    # __ne__ is not implemented, check if an error NotImplementedError is return in unit test

    def __repr__(self) -> str:
        return self.__str__()

    def __str__(self):
        """
        return a string representation of the rational linear expression.
        Output:
        A string
        Examples:

        """
        res = "linear expression: {le} lcm = {lcm}".format(
            le=str(self.linear_expression),
            lcm=self.lcm
        )
        return res

    def all_homogeneous_terms_are_zero(self) -> bool:
        return self.linear_expression.all_homogeneous_terms_are_zero()

    def is_zero(self) -> bool:
        return self.integer_inhomogeneous_term() == 0 and self.all_homogeneous_terms_are_zero()

    def space_dimension(self) -> int:
        return self.linear_expression.space_dimension()

    def permute_space_dimensions(self, cycle: List[int]):
        self.linear_expression.permute_space_dimensions(cycle)

    def remove_space_dimensions(self, variable_set: ppl.Variables_Set):
        self.linear_expression.remove_space_dimensions(variable_set)

    def reduce(self):
        """
        reduce the number of useless zeros of the coefficients of the RLE (the last one)
        """
        coeffs = list(self.integer_coefficients())
        while len(coeffs) > 0 and coeffs[-1] == 0:

            coeffs.pop(-1)
        self.linear_expression = ppl.Linear_Expression(coeffs,
                                                       self.integer_inhomogeneous_term())

    def partial_evaluation(self, variable: ppl.Variable, value: Union[Fraction, int, RationalLinearExpression]) -> \
            RationalLinearExpression:

        var_rle = self.linear_expression_as_rational_linear_expression(variable * self.integer_coefficient(variable))
        var_rle.lcm = self.lcm
        coeff = self.coefficient(variable)

        return self - var_rle + coeff * value


class InfiniteExpression(AbstractExpression):
    __slots__ = ["is_positive"]

    def __init__(self, sign: bool = True):
        if type(sign) == bool:
            self.is_positive = sign
        else:
            raise TypeError(str(sign) + " should be a boolean")

    def __repr__(self):
        return str(self.inhomogeneous_term())

    def inhomogeneous_term(self) -> float:
        return math.inf if self.is_positive else -math.inf

    @staticmethod
    def is_zero() -> bool:
        return False

    def is_equal_to(self, other: InfiniteExpression):
        if type(other) != InfiniteExpression:
            return False
        else:
            return self.is_positive == other.is_positive

    @staticmethod
    def all_homogeneous_terms_are_zero() -> bool:
        return False

    # Arithmetic operations

    def __add__(self,
                other: Union[InfiniteExpression, RationalLinearExpression, ppl.Linear_Expression,
                             int, Fraction]
                ) -> InfiniteExpression:

        other_ty = type(other)
        if other_ty in [RationalLinearExpression, ppl.Linear_Expression, int, Fraction]:
            return self
        elif other_ty == InfiniteExpression and self.is_positive == other.is_positive:
            return self
        else:
            raise NotImplementedError

    def __neg__(self):
        return InfiniteExpression(not self.is_positive)

    def __sub__(self, other):
        return self.__add__(-other)

    def __mul__(self, other: Union[int, Fraction, RationalLinearExpression, InfiniteExpression]
                ) -> InfiniteExpression:
        other_ty = type(other)
        if other_ty in [Fraction, int]:
            if other == 0:
                raise NotImplementedError
            return self if other > 0 else - self

        elif other_ty in [ppl.Linear_Expression, RationalLinearExpression]:
            # if self.all_homogeneous_terms_are_zero:
            #     return self.__mul__(other.inhomogeneous_term())
            # else:
            raise NotImplementedError
        elif other_ty == InfiniteExpression:
            return InfiniteExpression(sign=self.is_positive == other.is_positive)

    # Reverse operations

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return -self.__add__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    # Creation of constraints

    def __gt__(self, other):
        if type(other) != InfiniteExpression and type(other) != RationalLinearExpression:
            return NotImplementedError
        if self.is_positive:
            if type(other) == InfiniteExpression and other.is_positive:
                return 1 < 0 * ppl.Variable(0)
            else:
                return 1 > 0 * ppl.Variable(0)
        else:
            return 1 < 0 * ppl.Variable(0)

    def __ge__(self, other):
        if self.is_positive:
            return 1 >= 0 * ppl.Variable(0)
        else:
            return 1 <= 0 * ppl.Variable(0)

    def __lt__(self, other):
        if type(other) != InfiniteExpression and type(other) != RationalLinearExpression:
            return NotImplementedError
        if not self.is_positive and type(other) == InfiniteExpression and not other.is_positive:
            return 1 < 0 * ppl.Variable(0)
        elif self.is_positive:
            return 1 < 0 * ppl.Variable(0)
        else:
            return 1 > 0 * ppl.Variable(0)

    def __le__(self, other):
        if self.is_positive:
            return 1 <= 0 * ppl.Variable(0)
        else:
            return 1 >= 0 * ppl.Variable(0)

    def __eq__(self, other):
        if type(other) == InfiniteExpression:
            if self.is_positive == other.is_positive:
                return 1 > 0 * ppl.Variable(0)
            else:
                return 1 < 0 * ppl.Variable(0)
