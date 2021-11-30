# coding=utf-8
from fractions import Fraction

import ppl
import pytest
from tests.test_technical_lemma_ex import *
import pysymbrobustness.ppl_extension.technical_lemma as tech
import pysymbrobustness.ppl_extension.rational_linear_algebra as rla
import ppl as ppl


class TestSubFunctionMethods:

    def test_get_coefficient(self):
        with pytest.raises(TypeError):
            tech.get_coefficient(rla.InfiniteExpression(sign=True), ppl.Variable(0))

        with pytest.raises(TypeError):
            tech.get_coefficient(ppl.Linear_Expression([5, 2]), ppl.Variable(1))

        fct = rla.RationalLinearExpression((3, Fraction(-4, 5), -9, 0, 1))
        v_2_coeff, b_2 = tech.get_coefficient(fct=fct, v=ppl.Variable(2))
        v_0_coeff, b_0 = tech.get_coefficient(fct=fct, v=ppl.Variable(0))
        v_end_coeff, b_end = tech.get_coefficient(fct=fct, v=ppl.Variable(4))
        assert v_0_coeff == Fraction(3)
        assert b_0 == rla.RationalLinearExpression((0, Fraction(-4, 5), -9, 0, 1))
        assert v_2_coeff == Fraction(-9, 1)
        assert b_2 == rla.RationalLinearExpression((3, Fraction(-4, 5), 0, 0, 1))
        assert v_end_coeff == Fraction(1)
        assert b_end == rla.RationalLinearExpression((3, Fraction(-4, 5), -9, 0, 0))

    def test_compute_permissiveness(self):
        pass

    def test_opti_f_inf_g_inf(self):
        pass

    def test_name_to_find(self):
        pass

    def test_a_c_same_sign(self):
        pass

    def test_opti_f_inf_g_rle(self):
        pass

    def test_opti_f_rle_g_inf(self):
        pass

    def test_opti_f_rle_g_rle(self):
        pass


class TestPlfOptimization():

    def test_plf_optimization(self):
        pass
