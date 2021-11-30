# coding=utf-8
import pytest
from gmpy2 import mpz
from numpy.polynomial import Polynomial

from tests.test_rational_linear_algebra_ex import *
from pysymbrobustness.ppl_extension.rational_linear_algebra import RationalLinearExpression as rle
from pysymbrobustness.ppl_extension.rational_linear_algebra import InfiniteExpression as inf_fct
from fractions import Fraction


class TestMethods:
    def test_fraction_as_lcm_tuple_list_exception(self):
        # Empty list, tuple or arguments

        with pytest.raises(TypeError):
            rla.fraction_as_lcm_tuple_list()

        with pytest.raises(TypeError):
            rla.fraction_as_lcm_tuple_list(())

        with pytest.raises(TypeError):
            rla.fraction_as_lcm_tuple_list([])

        # List of elements containing floats

        with pytest.raises(TypeError):
            rla.fraction_as_lcm_tuple_list([Fraction(1, 5), 2, 0.33])

    def test_fraction_as_lcm_tuple_list(self, rla_with_an_integer, rla_with_only_zeros, rla_with_only_a_zero,
                                        rla_with_only_integers, rla_with_neg_and_pos_integers,
                                        rla_with_pos_neg_integers_and_zeros, rla_with_mixed_int_and_fractions_and_zeros,
                                        rla_with_only_fractions, rla_with_only_fractions_permutation,
                                        rla_with_negative_fractions):
        # Only integers
        assert rla_with_only_integers == ([4, 6, 23, 8], 1)  # Testing with only integers
        assert rla_with_neg_and_pos_integers == ([-4, 6, 23, -8], 1)  # Adding neg integers
        assert rla_with_pos_neg_integers_and_zeros == ([-4, 0, 6, 0, 0, -23, 8], 1)  # Adding zeros

        # Only zeros
        assert rla_with_only_zeros == ([0, 0, 0, 0], 1)  # Testing with only zeros
        assert rla_with_only_a_zero == ([0], 1)  # Testing with one zeros

        # Only Fractions
        assert rla_with_only_fractions == ([2 * 33, 11 * 11 * 7, 2 * 7], 33 * 7)
        assert rla_with_negative_fractions == ([2 * 33, -11 * 11 * 7, -2 * 7], 33 * 7)

        # Mixing Fraction & integers
        assert rla_with_an_integer == ([1, 10], 5)  # Testing a mix of an integer and a fraction
        assert rla_with_mixed_int_and_fractions_and_zeros == ([4, -16, -5, 0, 100, -200, 0, 80, 10], 20)
        assert rla_with_only_fractions_permutation == ([11 * 11 * 7, 2 * 33, 2 * 7], 33 * 7)  # Testing permutations


class TestRationalLinearExpression:

    def test_init_exception(self):
        # Rational linear expression with floats

        with pytest.raises(TypeError):
            rla.RationalLinearExpression(
                homogeneous_terms=(Fraction(1, 5), 2, 0.33), inhomogeneous_term=Fraction(2, 33)
            )

        with pytest.raises(TypeError):
            rla.fraction_as_lcm_tuple_list([0.44, 0.3445, 0.345345])

        with pytest.raises(TypeError):
            rla.RationalLinearExpression(
                homogeneous_terms=(0.44, 0.3445, 0.345345), binhomogeneous_term=Fraction(2, 33)
            )

        with pytest.raises(TypeError):
            rla.RationalLinearExpression(inhomogeneous_term=0.345345)

        with pytest.raises(TypeError):
            rla.fraction_as_lcm_tuple_list([Fraction(1, 5), 2, math.inf])

        with pytest.raises(TypeError):
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 5), 2, math.inf))

        with pytest.raises(TypeError):
            rla.RationalLinearExpression(
                homogeneous_terms=(Fraction(1, 5), 2), inhomogeneous_term=math.inf)

        # Changing div for 0

        with pytest.raises(ZeroDivisionError):
            rla.RationalLinearExpression(
                homogeneous_terms=(Fraction(1, 5), Fraction(-4, 5), -Fraction(1, 4), 0, 5, -10, 0, 4),
                inhomogeneous_term=Fraction(3, 6),
                div=0)

    def test_init(self, linear_expression_empty, linear_expression_with_an_integer, linear_expression_with_only_zeros,
                  linear_expression_with_only_a_zero, linear_expression_with_only_integers,
                  linear_expression_with_pos_and_neg_integers, linear_expression_with_pos_neg_integers_and_zeros,
                  linear_expression_with_only_fractions, linear_expression_with_only_fractions_permutation,
                  linear_expression_with_only_fractions_div_33, linear_expression_with_only_fractions_div_1650,
                  linear_expression_with_only_fractions_zero_inhomogeneous_term,
                  linear_expression_with_negative_fractions,
                  linear_expression_with_mixed_int_and_fractions_and_zeros, linear_expression_with_only_constant,
                  linear_expression_with_only_neg_constant):
        # Only integers:
        # 1) Testing empty linear expressions
        assert linear_expression_empty.linear_expression.is_equal_to(ppl.Linear_Expression())
        assert linear_expression_empty.lcm == 1

        # 2) Only zeros
        assert linear_expression_with_only_zeros.linear_expression.is_equal_to(ppl.Linear_Expression([0, 0, 0, 0], 0))
        assert linear_expression_with_only_zeros.lcm == 1

        assert linear_expression_with_only_a_zero.linear_expression.is_equal_to(ppl.Linear_Expression([0], 0))
        assert linear_expression_with_only_a_zero.lcm == 1
        # 3) pos integers
        assert linear_expression_with_only_integers.linear_expression.is_equal_to(ppl.Linear_Expression([4, 6, 23], 8))
        assert linear_expression_with_only_integers.lcm == 1
        # 4) Adding neg int
        assert linear_expression_with_pos_and_neg_integers.linear_expression.is_equal_to(
            ppl.Linear_Expression([-4, 6, 23], -8)
        )
        assert linear_expression_with_pos_and_neg_integers.lcm == 1
        # # Adding zeros
        assert linear_expression_with_pos_neg_integers_and_zeros.linear_expression.is_equal_to(
            ppl.Linear_Expression([-4, 0, 6, 0, 0, -23], 8)
        )
        assert linear_expression_with_pos_neg_integers_and_zeros.lcm == 1
        # Only Fractions
        assert linear_expression_with_only_fractions.linear_expression.is_equal_to(
            ppl.Linear_Expression([2 * 33, 11 * 11 * 7], 2 * 7)
        )
        assert linear_expression_with_only_fractions.lcm == 33 * 7

        assert linear_expression_with_only_fractions_permutation.linear_expression.is_equal_to(
            ppl.Linear_Expression([11 * 11 * 7, 2 * 33], 2 * 7)
        )
        assert linear_expression_with_only_fractions_permutation.lcm == 33 * 7

        assert linear_expression_with_only_fractions_div_33.linear_expression.is_equal_to(
            ppl.Linear_Expression([2 * 33, 11 * 11 * 7], 2 * 7)
        )
        assert linear_expression_with_only_fractions_div_33.lcm == 33 * 7

        assert linear_expression_with_only_fractions_div_1650.linear_expression.is_equal_to(
            ppl.Linear_Expression([2 * 33, 11 * 11 * 7], 2 * 7)
        )
        assert linear_expression_with_only_fractions_div_1650.lcm == 33 * 7

        assert linear_expression_with_only_fractions_zero_inhomogeneous_term.linear_expression.is_equal_to(
            ppl.Linear_Expression([2 * 33, 11 * 11 * 7, 2 * 7], 0))
        assert linear_expression_with_only_fractions_zero_inhomogeneous_term.lcm == 33 * 7

        assert linear_expression_with_negative_fractions.linear_expression.is_equal_to(
            ppl.Linear_Expression([2 * 33, -11 * 11 * 7], -2 * 7)
        )
        assert linear_expression_with_negative_fractions.lcm == 33 * 7
        # Mixing int and Fractions
        assert linear_expression_with_an_integer.linear_expression.is_equal_to(ppl.Linear_Expression([1], 10))
        assert linear_expression_with_an_integer.lcm == 5

        assert linear_expression_with_mixed_int_and_fractions_and_zeros.linear_expression.is_equal_to(
            ppl.Linear_Expression([4, -16, -5, 0, 100, -200, 0, 80], 10))
        assert linear_expression_with_mixed_int_and_fractions_and_zeros.lcm == 20
        # Constant linear expressions
        assert linear_expression_with_only_constant.linear_expression.is_equal_to(ppl.Linear_Expression([], 1))
        assert linear_expression_with_only_constant.lcm == 2

        assert linear_expression_with_only_neg_constant.linear_expression.is_equal_to(ppl.Linear_Expression([], -15))
        assert linear_expression_with_only_neg_constant.lcm == 34

    def test_linear_expression_as_rational_linear_expression(self):
        # Only integers:
        # 1) Testing empty linear expressions
        LE = rle.linear_expression_as_rational_linear_expression(ppl.Linear_Expression())
        assert LE.linear_expression.is_equal_to(ppl.Linear_Expression())
        assert LE.lcm == 1

        # 2) Only zeros
        LE = rle.linear_expression_as_rational_linear_expression(ppl.Linear_Expression([0, 0, 0, 0], 0))
        assert LE.linear_expression.is_equal_to(ppl.Linear_Expression([0, 0, 0, 0], 0))
        assert LE.lcm == 1

        LE = rle.linear_expression_as_rational_linear_expression(ppl.Linear_Expression([0], 0))
        assert LE.linear_expression.is_equal_to(ppl.Linear_Expression([0], 0))
        assert LE.lcm == 1
        # 3) pos integers
        LE = rle.linear_expression_as_rational_linear_expression(ppl.Linear_Expression([4, 6, 23], 8))
        assert LE.linear_expression.is_equal_to(ppl.Linear_Expression([4, 6, 23], 8))
        assert LE.lcm == 1
        # 4) Adding neg int
        LE = rle.linear_expression_as_rational_linear_expression(ppl.Linear_Expression([-4, 6, 23], -8))
        assert LE.linear_expression.is_equal_to(ppl.Linear_Expression([-4, 6, 23], -8))
        assert LE.lcm == 1
        # # Adding zeros
        LE = rle.linear_expression_as_rational_linear_expression(ppl.Linear_Expression([-4, 0, 6, 0, 0, -23], 8))
        assert LE.linear_expression.is_equal_to(ppl.Linear_Expression([-4, 0, 6, 0, 0, -23], 8))
        assert LE.lcm == 1

    def test_constant_as_rational_linear_expression(self):
        cste_RLE = rle.constant_as_rational_linear_expression(0)
        assert cste_RLE.linear_expression.is_equal_to(ppl.Linear_Expression([], 0))
        assert cste_RLE.lcm == 1

        cste_RLE = rle.constant_as_rational_linear_expression(134)
        assert cste_RLE.linear_expression.is_equal_to(ppl.Linear_Expression([], 134))
        assert cste_RLE.lcm == 1

        cste_RLE = rle.constant_as_rational_linear_expression(-1513)
        assert cste_RLE.linear_expression.is_equal_to(ppl.Linear_Expression([], -1513))
        assert cste_RLE.lcm == 1

        cste_RLE = rle.constant_as_rational_linear_expression(Fraction(-4, 6))
        assert cste_RLE.linear_expression.is_equal_to(ppl.Linear_Expression([], -2))
        assert cste_RLE.lcm == 3

        cste_RLE = rle.constant_as_rational_linear_expression(Fraction(1))
        assert cste_RLE.linear_expression.is_equal_to(ppl.Linear_Expression([], 1))
        assert cste_RLE.lcm == 1

        cste_RLE = rle.constant_as_rational_linear_expression(Fraction(mpz(34), 3))
        assert cste_RLE.linear_expression.is_equal_to(ppl.Linear_Expression([], 34))
        assert cste_RLE.lcm == 3

        with pytest.raises(TypeError):
            rle.constant_as_rational_linear_expression(0.235245)

        with pytest.raises(TypeError):
            rle.constant_as_rational_linear_expression(math.inf)

        with pytest.raises(TypeError):
            rle.constant_as_rational_linear_expression(0.235245)

    def test_integer_coefficients(self, linear_expression_empty,
                                  linear_expression_with_mixed_int_and_fractions_and_zeros):
        assert linear_expression_with_mixed_int_and_fractions_and_zeros.integer_coefficients() == \
               (1 * 4, -4 * 4, -5, 0, 5 * 20, -10 * 20, 0, 4 * 20)

        assert linear_expression_empty.integer_coefficients() == ()

    def test_rational_coefficients(self, linear_expression_empty, linear_expression_with_only_zeros,
                                   linear_expression_with_mixed_int_and_fractions_and_zeros):
        assert linear_expression_empty.coefficients() == ()
        assert linear_expression_with_mixed_int_and_fractions_and_zeros.coefficients() == \
               (Fraction(1, 5), Fraction(-4, 5), -Fraction(1, 4), 0, 5, -10, 0, 4)
        assert linear_expression_with_only_zeros.coefficients() == (0, 0, 0, 0)

    def test_rational_inhomogeneous_terms(self, linear_expression_empty, linear_expression_with_only_zeros,
                                          linear_expression_with_mixed_int_and_fractions_and_zeros):
        assert linear_expression_empty.inhomogeneous_term() == 0
        assert linear_expression_with_mixed_int_and_fractions_and_zeros.inhomogeneous_term() == \
               Fraction(3, 6)
        assert linear_expression_with_only_zeros.inhomogeneous_term() == 0

    def test_integer_inhomogeneous_term(self, linear_expression_empty, linear_expression_with_only_zeros,
                                        linear_expression_with_mixed_int_and_fractions_and_zeros):
        assert linear_expression_empty.integer_inhomogeneous_term() == 0
        assert linear_expression_with_mixed_int_and_fractions_and_zeros.integer_inhomogeneous_term() == \
               10
        assert linear_expression_with_only_zeros.integer_inhomogeneous_term() == 0

    def test_all_homogeneous_terms_are_zero(self, linear_expression_empty, linear_expression_with_only_zeros,
                                            linear_expression_with_only_a_zero,
                                            linear_expression_with_mixed_int_and_fractions_and_zeros):
        assert linear_expression_empty.all_homogeneous_terms_are_zero()
        assert not linear_expression_with_mixed_int_and_fractions_and_zeros.all_homogeneous_terms_are_zero()
        assert linear_expression_with_only_zeros.all_homogeneous_terms_are_zero()
        assert linear_expression_with_only_a_zero.all_homogeneous_terms_are_zero()

    def test_is_zero(self, linear_expression_empty, linear_expression_with_only_zeros,
                     linear_expression_with_only_a_zero,
                     linear_expression_with_mixed_int_and_fractions_and_zeros):
        assert linear_expression_empty.is_zero()
        assert not linear_expression_with_mixed_int_and_fractions_and_zeros.is_zero()
        assert linear_expression_with_only_zeros.is_zero()
        assert linear_expression_with_only_a_zero.is_zero()

    def test_is_equal_to(self, linear_expression_empty, linear_expression_with_only_zeros,
                         linear_expression_with_only_a_zero,
                         linear_expression_with_only_fractions,
                         linear_expression_with_only_fractions_div_33,
                         linear_expression_with_only_fractions_div_1650):
        assert not linear_expression_empty.is_equal_to(linear_expression_with_only_a_zero)
        assert linear_expression_empty.is_equal_to(rla.RationalLinearExpression(None))
        clone_zero = rla.RationalLinearExpression(homogeneous_terms=tuple([0, 0, 0, 0]), inhomogeneous_term=0)
        assert linear_expression_with_only_zeros.is_equal_to(clone_zero)
        assert linear_expression_with_only_fractions.is_equal_to(linear_expression_with_only_fractions_div_33)
        assert linear_expression_with_only_fractions_div_33.is_equal_to(linear_expression_with_only_fractions_div_1650)

    def test_space_dimension(self, linear_expression_empty, linear_expression_with_only_zeros,
                             linear_expression_with_only_a_zero,
                             linear_expression_with_mixed_int_and_fractions_and_zeros):
        assert linear_expression_empty.space_dimension() == 0
        assert linear_expression_with_mixed_int_and_fractions_and_zeros.space_dimension() == 8
        assert linear_expression_with_only_zeros.space_dimension() == 4
        assert linear_expression_with_only_a_zero.space_dimension() == 1

    def test_permute_space_dimensions(self, linear_expression_empty, linear_expression_with_only_zeros,
                                      linear_expression_with_only_a_zero,
                                      linear_expression_with_only_fractions,
                                      linear_expression_with_only_fractions_permutation):
        clone_zero = linear_expression_with_only_zeros
        clone_zero.permute_space_dimensions([0, 1, 2, 3])
        assert linear_expression_with_only_zeros.is_equal_to(clone_zero)
        new_linear_expression_with_only_fractions = linear_expression_with_only_fractions
        new_linear_expression_with_only_fractions.permute_space_dimensions([0, 1])
        assert new_linear_expression_with_only_fractions.is_equal_to(
            linear_expression_with_only_fractions_permutation)

    def test_remove_space_dimension(self, linear_expression_empty, linear_expression_with_only_zeros,
                                    linear_expression_with_only_a_zero,
                                    linear_expression_with_mixed_int_and_fractions_and_zeros):
        new_linear_expression_with_only_a_zero = linear_expression_with_only_zeros
        new_linear_expression_with_only_a_zero.remove_space_dimensions(ppl.Variables_Set(1, 3))
        assert new_linear_expression_with_only_a_zero.is_equal_to(linear_expression_with_only_a_zero)

    def test_addition(self, linear_expression_empty, linear_expression_with_only_zeros,
                      linear_expression_with_an_integer, linear_expression_with_only_a_zero,
                      linear_expression_with_mixed_int_and_fractions_and_zeros,
                      rational_linear_expression_0, rational_linear_expression_1, rational_linear_expression_2,
                      l0_minus_l1, l0_minus_l2, l0_plus_l1, l0_plus_l2, l1_plus_l2, l1_minus_l2):
        assert (linear_expression_empty + linear_expression_with_mixed_int_and_fractions_and_zeros).is_equal_to(
            linear_expression_with_mixed_int_and_fractions_and_zeros)
        assert (linear_expression_with_only_a_zero + linear_expression_with_mixed_int_and_fractions_and_zeros). \
            is_equal_to(linear_expression_with_mixed_int_and_fractions_and_zeros)
        assert (linear_expression_with_only_zeros + linear_expression_with_mixed_int_and_fractions_and_zeros). \
            is_equal_to(linear_expression_with_mixed_int_and_fractions_and_zeros)
        # Verify that the number of zeros is counted
        assert not (linear_expression_with_only_zeros + linear_expression_with_an_integer).is_equal_to(
            linear_expression_with_an_integer)
        assert (linear_expression_with_only_a_zero + linear_expression_with_an_integer).is_equal_to(
            linear_expression_with_an_integer)

        assert (rational_linear_expression_0 + rational_linear_expression_1).is_equal_to(l0_plus_l1)

        assert (rational_linear_expression_0 + rational_linear_expression_2).is_equal_to(l0_plus_l2)

        assert (rational_linear_expression_1 + rational_linear_expression_2).is_equal_to(l1_plus_l2)

        assert (rational_linear_expression_1 + rational_linear_expression_0).is_equal_to(l0_plus_l1)

        assert (rational_linear_expression_2 + rational_linear_expression_0).is_equal_to(l0_plus_l2)

        assert (rational_linear_expression_2 + rational_linear_expression_1).is_equal_to(l1_plus_l2)

        assert rational_linear_expression_0 + math.inf == math.inf

        assert rational_linear_expression_0 - math.inf == - math.inf

        assert (rational_linear_expression_0 + 10).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7)),
                                         inhomogeneous_term=10))

        assert (rational_linear_expression_0 + ppl.Linear_Expression([], 10)).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7)),
                                         inhomogeneous_term=10))

        assert (rational_linear_expression_0 + Fraction(-11, 45)).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7)),
                                         inhomogeneous_term=Fraction(-11, 45)))

        with pytest.raises(NotImplementedError):
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7)),
                                         inhomogeneous_term=Fraction(-11, 45)).__add__(0.254543534)

        with pytest.raises(NotImplementedError):
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7)),
                                         inhomogeneous_term=Fraction(-11, 45)).__add__(Polynomial([1, 2, 3]))

    def test_reverse_addition(self, linear_expression_empty, linear_expression_with_only_zeros,
                              linear_expression_with_an_integer, linear_expression_with_only_a_zero,
                              linear_expression_with_mixed_int_and_fractions_and_zeros,
                              rational_linear_expression_0):
        # Testing __radd__ function
        assert (linear_expression_with_mixed_int_and_fractions_and_zeros + linear_expression_empty).is_equal_to(
            linear_expression_with_mixed_int_and_fractions_and_zeros)
        assert (
                linear_expression_with_mixed_int_and_fractions_and_zeros + linear_expression_with_only_a_zero).is_equal_to(
            linear_expression_with_mixed_int_and_fractions_and_zeros)
        assert (
                linear_expression_with_mixed_int_and_fractions_and_zeros + linear_expression_with_only_zeros).is_equal_to(
            linear_expression_with_mixed_int_and_fractions_and_zeros)
        # Verify that the number of zeros is counted
        assert not (linear_expression_with_an_integer + linear_expression_with_only_zeros).is_equal_to(
            linear_expression_with_an_integer)
        assert (linear_expression_with_an_integer + linear_expression_with_only_a_zero).is_equal_to(
            linear_expression_with_an_integer)

        assert math.inf + rational_linear_expression_0 == math.inf

        assert (10 + rational_linear_expression_0).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7))
                                         , inhomogeneous_term=10))

        assert (Fraction(-11, 45) + rational_linear_expression_0).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7)),
                                         inhomogeneous_term=Fraction(-11, 45)))

        with pytest.raises(NotImplementedError):
            (0.254543534 + rla.RationalLinearExpression(homogeneous_terms=(2, 1, Fraction(4, 7))))

    def test_substraction(self, rational_linear_expression_1, rational_linear_expression_2,
                          rational_linear_expression_0, l0_minus_l2, l1_minus_l2, l0_minus_l1):
        assert (rational_linear_expression_0 - rational_linear_expression_1).is_equal_to(l0_minus_l1)

        assert (rational_linear_expression_0 - rational_linear_expression_2).is_equal_to(l0_minus_l2)

        assert (rational_linear_expression_1 - rational_linear_expression_2).is_equal_to(l1_minus_l2)

        assert (rational_linear_expression_1 - rational_linear_expression_0).is_equal_to(-l0_minus_l1)

        assert (rational_linear_expression_2 - rational_linear_expression_0).is_equal_to(-l0_minus_l2)

        assert (rational_linear_expression_2 - rational_linear_expression_1).is_equal_to(-l1_minus_l2)

        assert rational_linear_expression_0 - (+math.inf) == -math.inf

        assert rational_linear_expression_0 - (-math.inf) == math.inf

        assert (rational_linear_expression_0 - 10).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7)),
                                         inhomogeneous_term=-10))

        assert (rational_linear_expression_0 - Fraction(-11, 45)).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(1, 1), Fraction(1, 2), Fraction(2, 7)),
                                         inhomogeneous_term=Fraction(11, 45)))

        with pytest.raises(NotImplementedError):
            rla.RationalLinearExpression(homogeneous_terms=(2, 1, Fraction(4, 7))) - 0.2454564

    def test_reverse_substraction(self, rational_linear_expression_1, rational_linear_expression_0,
                                  rational_linear_expression_2):
        assert (+math.inf) - rational_linear_expression_0 == math.inf

        assert (-math.inf) - rational_linear_expression_0 == -math.inf

        assert (10 - rational_linear_expression_0).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(-1, 1), Fraction(-1, 2), Fraction(-2, 7)),
                                         inhomogeneous_term=+10))

        assert (Fraction(11, 45) - rational_linear_expression_0).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(Fraction(-1, 1), Fraction(-1, 2), Fraction(-2, 7)),
                                         inhomogeneous_term=Fraction(11, 45)))

        with pytest.raises(NotImplementedError):
            0.245456 - rla.RationalLinearExpression(homogeneous_terms=(2, 1, Fraction(4, 7)))

    def test_multiplication(self, rational_linear_expression_0):
        assert (2 * rational_linear_expression_0).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(2, 1, Fraction(4, 7))))

        with pytest.raises(NotImplementedError):
            rational_linear_expression_0 * math.inf

        assert (rational_linear_expression_0 * 2).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(2, 1, Fraction(4, 7))))

        assert (rational_linear_expression_0 * 0).is_equal_to(
            rla.RationalLinearExpression(homogeneous_terms=(0, 0, 0)))

        with pytest.raises(NotImplementedError):
            0.2454564 * rla.RationalLinearExpression(homogeneous_terms=(0, 0, 0))

        with pytest.raises(NotImplementedError):
            rla.RationalLinearExpression(homogeneous_terms=(0, 0, 0)) * 0.45346

    def test_constraints(self, rational_linear_expression_0, linear_expression_with_only_zeros):
        assert (rational_linear_expression_0 <= 0).is_equivalent_to(
            14 * ppl.Variable(0) + 7 * ppl.Variable(1) + 4 * ppl.Variable(2) <= 0)

        assert (rational_linear_expression_0 > 0).is_equivalent_to(
            14 * ppl.Variable(0) + 7 * ppl.Variable(1) + 4 * ppl.Variable(2) > 0)

        assert (rational_linear_expression_0 < 0).is_equivalent_to(
            14 * ppl.Variable(0) + 7 * ppl.Variable(1) + 4 * ppl.Variable(2) < 0)

        assert not (linear_expression_with_only_zeros <= 0).is_equivalent_to(0 * ppl.Variable(0) <= 1)
        # Not the same space dimension...

        assert (linear_expression_with_only_zeros <= 1).is_equivalent_to(0 * ppl.Variable(0) + 0 * ppl.Variable(1) +
                                                                         0 * ppl.Variable(2) + 0 * ppl.Variable(3) <= 2)

        assert (linear_expression_with_only_zeros >= 0).is_tautological()

        assert (linear_expression_with_only_zeros == 0).is_tautological()

        assert (linear_expression_with_only_zeros >= -1).is_equivalent_to(0 * ppl.Variable(0) + 0 * ppl.Variable(1) +
                                                                          0 * ppl.Variable(2) + 0 * ppl.Variable(
            3) <= 2)

        assert (linear_expression_with_only_zeros == 1).is_equivalent_to(0 * ppl.Variable(0) + 0 * ppl.Variable(1) +
                                                                         0 * ppl.Variable(2) + 0 * ppl.Variable(3) == 1)

    def test_reduce(self, linear_expression_with_only_zeros, linear_expression_empty,
                    linear_expression_with_three_zeros,
                    linear_expression_reduced):
        linear_expression_with_three_zeros.reduce()
        assert linear_expression_with_three_zeros.is_equal_to(linear_expression_reduced)

        linear_expression_with_only_zeros.reduce()
        assert linear_expression_with_only_zeros.is_equal_to(linear_expression_empty)

    def test_partial_evaluation(self):
        rle = rla.RationalLinearExpression((1, 0))
        evaluation = rla.RationalLinearExpression((0, -1), 1)
        res = rle.partial_evaluation(variable=ppl.Variable(0),
                                     value=evaluation)
        assert res == evaluation


class TestInfiniteExpression:
    def test_init(self):
        assert inf_fct(sign=True).is_positive
        assert not inf_fct(sign=False).is_positive
        with pytest.raises(TypeError):
            inf_fct(5)
        with pytest.raises(TypeError):
            inf_fct([4, 6, 7], 4)

    def test_inhomogeneous_term(self, plus_inf_fct, minus_inf_fct):
        assert plus_inf_fct.inhomogeneous_term() == math.inf
        assert minus_inf_fct.inhomogeneous_term() == -math.inf

    def test_is_zero(self, plus_inf_fct, minus_inf_fct):
        assert not plus_inf_fct.is_zero()
        assert not minus_inf_fct.is_zero()

    def test_all_homogeneous_terms_are_zero(self):
        assert not inf_fct(sign=True).all_homogeneous_terms_are_zero()
        assert not inf_fct(sign=False).all_homogeneous_terms_are_zero()

    def test_is_equal_to(self, plus_inf_fct, minus_inf_fct, rational_linear_expression_1):
        assert plus_inf_fct.is_equal_to(plus_inf_fct)
        assert not plus_inf_fct.is_equal_to(minus_inf_fct)
        assert not minus_inf_fct.is_equal_to(plus_inf_fct)
        assert minus_inf_fct.is_equal_to(minus_inf_fct)
        assert not plus_inf_fct.is_equal_to(rational_linear_expression_1)
        assert not minus_inf_fct.is_equal_to(rational_linear_expression_1)

    def test_add(self, plus_inf_fct, minus_inf_fct, rational_linear_expression_1):
        assert (plus_inf_fct + rational_linear_expression_1).is_equal_to(plus_inf_fct)
        assert (minus_inf_fct + rational_linear_expression_1).is_equal_to(minus_inf_fct)
        assert (rational_linear_expression_1 + minus_inf_fct).is_equal_to(minus_inf_fct)
        assert (minus_inf_fct + minus_inf_fct).is_equal_to(minus_inf_fct)
        assert (plus_inf_fct + plus_inf_fct).is_equal_to(plus_inf_fct)
        assert (plus_inf_fct + 1).is_equal_to(plus_inf_fct)
        assert (1 + plus_inf_fct).is_equal_to(plus_inf_fct)

        with pytest.raises(NotImplementedError):
            plus_inf_fct + minus_inf_fct

        with pytest.raises(NotImplementedError):
            minus_inf_fct + plus_inf_fct

        with pytest.raises(NotImplementedError):
            minus_inf_fct + math.inf

    def test_neg(self, plus_inf_fct, minus_inf_fct):
        assert (- plus_inf_fct).is_equal_to(minus_inf_fct)
        assert (- minus_inf_fct).is_equal_to(plus_inf_fct)

    def test_sub(self, plus_inf_fct, minus_inf_fct, rational_linear_expression_1):
        assert (plus_inf_fct - rational_linear_expression_1).is_equal_to(plus_inf_fct)
        assert (minus_inf_fct - rational_linear_expression_1).is_equal_to(minus_inf_fct)
        assert (rational_linear_expression_1 - minus_inf_fct).is_equal_to(plus_inf_fct)
        assert (minus_inf_fct - plus_inf_fct).is_equal_to(minus_inf_fct)
        assert (plus_inf_fct - minus_inf_fct).is_equal_to(plus_inf_fct)
        assert (plus_inf_fct - 1).is_equal_to(plus_inf_fct)
        assert (1 - plus_inf_fct).is_equal_to(minus_inf_fct)

    def test_mult(self, plus_inf_fct, minus_inf_fct, rational_linear_expression_1):
        assert (2 * plus_inf_fct).is_equal_to(plus_inf_fct)
        assert (2 * minus_inf_fct).is_equal_to(minus_inf_fct)
        assert (-2 * plus_inf_fct).is_equal_to(minus_inf_fct)
        assert (-3 * minus_inf_fct).is_equal_to(plus_inf_fct)
        assert (Fraction(3, 4) * minus_inf_fct).is_equal_to(minus_inf_fct)

        with pytest.raises(NotImplementedError):
            0 * minus_inf_fct
        with pytest.raises(NotImplementedError):
            0 * plus_inf_fct

        with pytest.raises(NotImplementedError):
            minus_inf_fct * 0
        with pytest.raises(NotImplementedError):
            plus_inf_fct * 0

        with pytest.raises(NotImplementedError):
            Fraction(0, 4) * minus_inf_fct

        with pytest.raises(NotImplementedError):
            rla.RationalLinearExpression([0], 4) * plus_inf_fct

        assert (minus_inf_fct * plus_inf_fct).is_equal_to(minus_inf_fct)
        assert (plus_inf_fct * plus_inf_fct).is_equal_to(plus_inf_fct)
        assert (minus_inf_fct * minus_inf_fct).is_equal_to(plus_inf_fct)

    def test_inequalities(self, minus_inf_fct, plus_inf_fct, rational_linear_expression_1):
        assert (plus_inf_fct >= rational_linear_expression_1).is_tautological()
        assert (plus_inf_fct <= rational_linear_expression_1).is_inconsistent()
        assert (minus_inf_fct <= rational_linear_expression_1).is_tautological()
        assert (minus_inf_fct >= rational_linear_expression_1).is_inconsistent()
        # assert (minus_inf_fct - 1 <= 0).is_tautological()

        assert (plus_inf_fct <= rational_linear_expression_1).is_inconsistent()
        assert (plus_inf_fct >= rational_linear_expression_1).is_tautological()
        assert (minus_inf_fct < rational_linear_expression_1).is_tautological()
        assert (minus_inf_fct > rational_linear_expression_1).is_inconsistent()
        # assert (minus_inf_fct - 1 < 0).is_tautological()

        assert (minus_inf_fct == minus_inf_fct).is_tautological()





    #
    # def test_eq(self, infinite, minus_infinite):
    #     assert (infinite.__eq__(ppl.Linear_Expression([4, 5], 3))).is_inconsistent()
    #     assert (minus_infinite.__eq__(
    #         ppl.Linear_Expression([4, 5], 3))).is_inconsistent()
    #     assert (infinite.__eq__(infinite)).is_tautological()
    #     assert (minus_infinite.__eq__(minus_infinite)).is_tautological()
    #
    # def test_gt(self, infinite, minus_infinite):
    #     assert (infinite.__gt__(ppl.Linear_Expression([4, 5], 3))).is_tautological()
    #     assert (minus_infinite.__gt__(
    #         ppl.Linear_Expression([4, 5], 3))).is_inconsistent()
    #     assert (infinite.__gt__(infinite)).is_tautological()
    #     assert (infinite.__gt__(minus_infinite)).is_tautological()
    #     assert (minus_infinite.__gt__(minus_infinite)).is_inconsistent()
    #
    # def test_lt(self, infinite, minus_infinite):
    #     assert (infinite.__lt__(ppl.Linear_Expression([4, 5], 3))).is_inconsistent()
    #     assert (minus_infinite.__lt__(
    #         ppl.Linear_Expression([4, 5], 3))).is_tautological()
    #     assert (infinite.__lt__(infinite)).is_inconsistent()
    #     assert (infinite.__lt__(minus_infinite)).is_inconsistent()
    #     assert (minus_infinite.__lt__(minus_infinite)).is_tautological()
