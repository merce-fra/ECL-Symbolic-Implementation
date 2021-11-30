# coding=utf-8
import pytest
import ppl as ppl
import pysymbrobustness.ppl_extension.linear_algebra as la
from tests.test_linear_algebra_ex import *
class TestConstraintsMethods:
    # Works in console, see later...
    # def test_cs_from_list(self,
    #                       empty_cs_list, cs_list_0, cs_list_1,
    #                       empty_constraint_system, constraint_system_0, constraint_system_1):
    #     assert la.cs_from_list(empty_cs_list) == empty_constraint_system
    #     assert la.cs_from_list(cs_list_0) == constraint_system_0
    #     assert la.cs_from_list(cs_list_1) == constraint_system_1
    #
    # def test_cs_linear_expressions(self,
    #                                empty_constraint_system, constraint_system_0, constraint_system_1,
    #                                empty_linear_expression, linear_expressions_0, linear_expressions_1):
    #
    #     list_0 = la.cs_linear_expressions(constraint_system_0)
    #     assert la.equality(list_0[0], linear_expressions_0[0])
    #     # assert la.equality(list_0[1], linear_expressions_0[1])
    #     # assert la.equality(list_0[2], linear_expressions_0[2])
    #     #
    #     # list_1 = la.cs_linear_expressions(constraint_system_1)
    #     # assert la.equality(list_1[0], linear_expressions_1[0])

    # def test_absolute_value_equality(self,
    #                                  linear_expression_0, linear_expression_0_mirror,
    #                                  linear_expression_1, linear_expression_1_mirror):
    #     assert la.absolute_value_equality(linear_expression_0, linear_expression_0_mirror) == (True, True)
    #     assert la.absolute_value_equality(linear_expression_0, linear_expression_0) == (False, True)
    #     assert la.absolute_value_equality(linear_expression_1, linear_expression_1) == (False, True)
    #     assert la.absolute_value_equality(linear_expression_1, linear_expression_1_mirror) == (True, True)
    #     assert la.absolute_value_equality(ppl.Linear_Expression([], 0), ppl.Linear_Expression([], -0)) == (False, True)
    #     assert la.absolute_value_equality(ppl.Linear_Expression([], 7), ppl.Linear_Expression([], 7)) == (False, True)
    #     assert la.absolute_value_equality(ppl.Linear_Expression([], 7), ppl.Linear_Expression([], -7)) == (True, True)
    #     assert not la.absolute_value_equality(linear_expression_0, linear_expression_1) == (False, True)
    #     assert not la.absolute_value_equality(ppl.Linear_Expression([], 7), linear_expression_1) == (False, True)

    def test_equality(self, linear_expression_0, linear_expression_0_mirror,
                      linear_expression_1, linear_expression_1_mirror):
        assert not la.equality(linear_expression_0, linear_expression_0_mirror)
        assert la.equality(linear_expression_0, linear_expression_0)
        assert la.equality(linear_expression_1, linear_expression_1)
        assert not la.equality(linear_expression_1, linear_expression_1_mirror)
        assert la.equality(ppl.Linear_Expression([], 0),ppl.Linear_Expression([], -0))
        assert la.equality(ppl.Linear_Expression([], 7),ppl.Linear_Expression([], 7))
        assert not la.equality(ppl.Linear_Expression([], 7), ppl.Linear_Expression([], -7))
        assert not la.equality(linear_expression_0, linear_expression_1)
        assert not la.equality(ppl.Linear_Expression([], 7), linear_expression_1)


