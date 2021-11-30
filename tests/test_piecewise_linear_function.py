# coding=utf-8
import pytest
import ppl as ppl
from tests.test_piecewise_linear_function_ex import *
import ppl as ppl
x = ppl.Variable(0)
y = ppl.Variable(1)


class TestSubSplines:
    def test_init(self):
        pass

    def test_minimum(self):
        pass

    def test_maximum(self):
        pass

    def test_min_max_op_spline(self):
        pass


class TestSplines:
    def test_init(self):
        pass

    def test_add_sub_spline(self, spline_0, sub_spline_with_strict_equality,
                            sub_spline_with_empty_poly_0,
                            sub_spline_with_empty_poly_1, sub_spline_0):
        old_spline_len = len(spline_0.sub_splines)
        spline_0.add_sub_spline(sub_spline_with_strict_equality)
        new_spline_len = len(spline_0.sub_splines)

        assert old_spline_len == new_spline_len

        old_spline_len = len(spline_0.sub_splines)
        spline_0.add_sub_spline(sub_spline_with_strict_equality)
        new_spline_len = len(spline_0.sub_splines)

        assert old_spline_len == new_spline_len

        old_spline_len = len(spline_0.sub_splines)
        spline_0.add_sub_spline(sub_spline_with_strict_equality)
        new_spline_len = len(spline_0.sub_splines)

        assert old_spline_len == new_spline_len

        old_spline_len = len(spline_0.sub_splines)
        spline_0.add_sub_spline(sub_spline_0)
        new_spline_len = len(spline_0.sub_splines)

        assert old_spline_len + 1 == new_spline_len

    def test_functions(self):
        # TODO: not urgent
        pass

    def test_polyhedrons(self):
        # TODO: not urgent
        pass

    def test_min(self):
        pass

    def test_max(self):
        pass

    def test_get_constraints(self):
        pass

    def test_perfect_partition_on_a_sub_spline(self):
        pass

    def test_perfect_partition(self):
        pass
