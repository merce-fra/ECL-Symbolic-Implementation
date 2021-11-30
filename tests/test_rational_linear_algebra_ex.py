# coding=utf-8
"""
==================================================
test_rational_linear_algebra_ex
==================================================
This module provide examples for the unit test of the module rational_linear_algebra

Type of examples: rational_linear_expression and tuple of fractions
------
Empty list, tuple or arguments (*)
Tuples with more than two elements, includings:
- only integers (*) (positive, negative, and zeros)
- only fractions (*) (same)
- only zeros (*)
- mixing fractions and integers (*)
- mixing fractions, integers and zero elements (at differences index: begin and end included, multiple zero) (*)
- mixing fractions, integers, zero elements and positive and negative elements (*)
- having div =0 (error expected) (*)
"""
import math

import ppl as ppl
import pytest
import pysymbrobustness.ppl_extension.rational_linear_algebra as rla
from fractions import Fraction

# Examples of empty RLE


@pytest.fixture(autouse=True)
def linear_expression_empty():
    return rla.RationalLinearExpression()


# Standard examples with two elements


@pytest.fixture(autouse=True)  # Standard example with an integer and a fraction
def rla_with_an_integer():
    return rla.fraction_as_lcm_tuple_list([Fraction(1, 5), 2])
    # Expected: ([1,10], div=5)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_an_integer():
    return rla.RationalLinearExpression(
        homogeneous_terms=tuple([Fraction(1, 5)]),
        inhomogeneous_term=2)
    # Expected: x+10, div=5


# Zeros examples

@pytest.fixture(autouse=True)  # Standard example with only (several) zeros
def rla_with_only_zeros():
    return rla.fraction_as_lcm_tuple_list([0, 0, 0, 0])
    # Expected: ([0,0,0,0], div=1)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_only_zeros():
    return rla.RationalLinearExpression(
        homogeneous_terms=tuple([0, 0, 0, 0]),
        inhomogeneous_term=0)

    # Expected: (0,0,0,0), div=1


@pytest.fixture(autouse=True)  # Standard example with only one zero
def rla_with_only_a_zero():
    return rla.fraction_as_lcm_tuple_list([0])
    # Expected: ([0], div=1)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_only_a_zero():
    return rla.RationalLinearExpression(
        homogeneous_terms=tuple([0]))
    # Expected: 0, div=1, should be equivalent tu linear_expression_empty


# Examples with only integers

@pytest.fixture(autouse=True)  # Standard example with positive integers
def rla_with_only_integers():
    return rla.fraction_as_lcm_tuple_list([4, 6, 23, 8])
    # Expected: ([4, 6, 23, 8], div=1)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_only_integers():
    return rla.RationalLinearExpression(
        homogeneous_terms=(4, 6, 23),
        inhomogeneous_term=8)
    # Expected: (4x0 + 6x1 + 23x2 + 8, div=1)


@pytest.fixture(autouse=True)  # Standard example with positive and negative integers
def rla_with_neg_and_pos_integers():
    return rla.fraction_as_lcm_tuple_list([-4, 6, 23, -8])
    # Expected: ([-4, 6, 23, -8], div=1)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_pos_and_neg_integers():
    return rla.RationalLinearExpression(
        homogeneous_terms=(-4, 6, 23),
        inhomogeneous_term=-8)
    # Expected: (-4x0 + 6x1 + 23x2 - 8, div=1)


@pytest.fixture(autouse=True)  # Standard example with positive and negative integers and zeros
def rla_with_pos_neg_integers_and_zeros():
    return rla.fraction_as_lcm_tuple_list([-4, 0, 6, 0, 0, -23, 8])
    # Expected: ([-4, 0, 6, 0, 0, -23, 8], div=1)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_pos_neg_integers_and_zeros():
    return rla.RationalLinearExpression(
        homogeneous_terms=(-4, 0, 6, 0, 0, -23),
        inhomogeneous_term=8)
    # Expected: (-4x0 + 6x2 - 23x5 + 8, div=1)


# Examples with only fractions


@pytest.fixture(autouse=True)  # Standard example with only fractions
def rla_with_only_fractions():
    return rla.fraction_as_lcm_tuple_list([Fraction(2, 7), Fraction(11, 3), Fraction(2, 33)])
    # Expected ([2*33, 11*11*7, 2*7], div = 33*7)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_only_fractions():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(2, 7), Fraction(11, 3)),
        inhomogeneous_term=Fraction(2, 33))
    # Expected: (2*33x0 + 11*11*7x1 + 2*7, div = 33*7)


@pytest.fixture(autouse=True)  # Example with permutations
def rla_with_only_fractions_permutation():
    return rla.fraction_as_lcm_tuple_list([Fraction(11, 3), Fraction(2, 7), Fraction(2, 33)])
    # Expected ([11*11*7, 2*33, 2*7], div = 33*7)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_only_fractions_permutation():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(11, 3), Fraction(2, 7)),
        inhomogeneous_term=Fraction(2, 33))
    # Expected: (2*33x1 + 11*11*7x0 + 2*7, div = 33*7)


@pytest.fixture(autouse=True)  # Associated RLE with non default div, goal: check that it does not changes
def linear_expression_with_only_fractions_div_33():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(2*33, 7), Fraction(11*11, 1)),
        inhomogeneous_term=Fraction(2, 1),
        div=33)
    # Expected (11*11*7x0 + 2*33*x1 + 2*7, div = 33*7)


@pytest.fixture(autouse=True)  # Associated RLE with non default div that should be reduced
def linear_expression_with_only_fractions_div_1650():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(2*1650, 7), Fraction(11*1650, 3)),
        inhomogeneous_term=Fraction(2*1650, 33),
        div=1650)
    # Expected (11*11*7x0 + 2*33*x1 + 2*7, div = 33*7)


@pytest.fixture(autouse=True)  # Associated RLE with zero inhomogeneous term
def linear_expression_with_only_fractions_zero_inhomogeneous_term():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(2, 7), Fraction(11, 3), Fraction(2, 33)),
        inhomogeneous_term=0)
    # Expected: (2*33x0 + 11*11*7x1 + 2*7x2, div = 33*7)


@pytest.fixture(autouse=True)  # Adding negative element (only fractions)
def rla_with_negative_fractions():
    return rla.fraction_as_lcm_tuple_list([Fraction(2, 7), Fraction(-11, 3), Fraction(-2, 33)])
    # Expected ([2*33, -11*11*7, -2*7], div = 33*7)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_negative_fractions():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(2, 7), Fraction(-11, 3)),
        inhomogeneous_term=Fraction(-2, 33))
    # Expected: (2*33x0 + 11*11*7x1 - 2*7, div = 33*7)


# Mixed example with integers and fractions


@pytest.fixture(autouse=True)  # Total mix: zero, fractions, integers, negative, positive
def rla_with_mixed_int_and_fractions_and_zeros():
    return rla.fraction_as_lcm_tuple_list(
        [Fraction(1, 5), Fraction(-4, 5), -Fraction(1, 4), 0, 5, -10, 0, 4, Fraction(3, 6)]
    )
    # Expected ([1*4, -4*4, -5, 0, 5*20, -10*20, 0, 4*20, 10 ], div = 5*4 = 20)


@pytest.fixture(autouse=True)  # Associated RLE
def linear_expression_with_mixed_int_and_fractions_and_zeros():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(1, 5), Fraction(-4, 5), -Fraction(1, 4), 0, 5, -10, 0, 4),
        inhomogeneous_term=Fraction(3, 6))
    # Expected (1*4x0 -4*4x1 -5x2 + 5*20x4 -10*20x5 + 4*20x7 + 10 , div = 20)


@pytest.fixture(autouse=True)  # RLE with a constant only
def linear_expression_with_only_constant():
    return rla.RationalLinearExpression(
        inhomogeneous_term=Fraction(3, 6))
    # Expected  (1 , div = 2)


@pytest.fixture(autouse=True) # RLE with a negative constant only
def linear_expression_with_only_neg_constant():
    return rla.RationalLinearExpression(
        inhomogeneous_term=Fraction(-15, 34))
    # Expected (-15 , div = 34)


# Tests for add, radd, sub, rsub, mul and rmul


@pytest.fixture(autouse=True)
def rational_linear_expression_0():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7)))


@pytest.fixture(autouse=True)
def rational_linear_expression_1():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(1, 3), 0,  Fraction(-2, 5)),
        inhomogeneous_term=1)


@pytest.fixture(autouse=True)
def rational_linear_expression_2():
    return rla.RationalLinearExpression(
        homogeneous_terms=(0, 1, 0),
        inhomogeneous_term=Fraction(1, 5))


@pytest.fixture(autouse=True)
def l0_plus_l1():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(4, 3), Fraction(1, 2), Fraction(-4, 35)),
        inhomogeneous_term=1)


@pytest.fixture(autouse=True)
def l0_plus_l2():
    return rla.RationalLinearExpression(
        homogeneous_terms=(1, Fraction(3, 2), Fraction(2, 7)),
        inhomogeneous_term=Fraction(1, 5))


@pytest.fixture(autouse=True)
def l1_plus_l2():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(1, 3), 1, Fraction(-2, 5)),
        inhomogeneous_term=Fraction(6, 5))


@pytest.fixture(autouse=True)
def l0_minus_l1():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(2, 3), Fraction(1, 2), Fraction(24, 35)),
        inhomogeneous_term=-1)


@pytest.fixture(autouse=True)
def l0_minus_l2():
    return rla.RationalLinearExpression(
        homogeneous_terms=(1, Fraction(-1, 2), Fraction(2, 7)),
        inhomogeneous_term=Fraction(-1, 5))


@pytest.fixture(autouse=True)
def l1_minus_l2():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(1, 3), -1, Fraction(-2, 5)),
        inhomogeneous_term=Fraction(4, 5))


# Rational Linear expression to reduce

@pytest.fixture(autouse=True)
def linear_expression_with_three_zeros():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(1, 3), -1, 0, Fraction(-2, 5), 0, 0 ,0 ),
        inhomogeneous_term=Fraction(4, 5))


@pytest.fixture(autouse=True)
def linear_expression_reduced():
    return rla.RationalLinearExpression(
        homogeneous_terms=(Fraction(1, 3), -1, 0, Fraction(-2, 5)),
        inhomogeneous_term=Fraction(4, 5))


@pytest.fixture(autouse=True)
def plus_inf_fct():
    return rla.InfiniteExpression(sign=True)


@pytest.fixture(autouse=True)
def minus_inf_fct():
    return rla.InfiniteExpression(sign=False)





