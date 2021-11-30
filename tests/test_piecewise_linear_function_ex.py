# coding=utf-8
from fractions import Fraction
import pysymbrobustness.ppl_extension.piecewise_linear_function as plf

import pytest
import ppl as ppl
import pysymbrobustness.ppl_extension.rational_linear_algebra as rla
import pysymbrobustness.ppl_extension.linear_algebra as la

x = ppl.Variable(0)
y = ppl.Variable(1)
z = ppl.Variable(2)


@pytest.fixture(autouse=True)
def formats_polyhedrons_0():
    """ formats partition examples for the computation of permissiveness for
    the timed automata with two clocks x,y
    and two linear transition 0 <= x <= 1, 0 <= y <= 1"""
    c_0 = [x - y >= 0, x <= 1, y <= 1]
    c_1 = [x - y >= 0, x <= 1, y <= 1]

    return ppl.C_Polyhedron(la.cs_from_list(c_0)), \
           ppl.C_Polyhedron(la.cs_from_list(c_1))


@pytest.fixture(autouse=True)
def formats_functions_0():
    """ formats functions examples for the computation of permissiveness for
    the timed automata with two clocks x,y and two linear transition 0 <= x
    <= 1, 0 <= y <= 1"""
    f_0 = rla.RationalLinearExpression((0, Fraction(-1, 4)), Fraction(1, 4))
    f_1 = rla.RationalLinearExpression((Fraction(-1, 4), 0), Fraction(1, 4))
    return f_0, f_1


@pytest.fixture(autouse=True)
def spline_0(formats_functions_0, formats_polyhedrons_0):
    c_0, c_1 = formats_polyhedrons_0
    f_0, f_1 = formats_functions_0
    return plf.Spline([plf.SubSpline(polyhedron=c_0, function=f_0),
                       plf.SubSpline(polyhedron=c_1, function=f_1)])


@pytest.fixture(autouse=True)
def global_sub_splines_0(formats_functions_0):
    c = [x <= 1, y <= 1]
    P = ppl.C_Polyhedron(la.cs_from_list(c))
    f_0, f_1 = formats_functions_0

    return [plf.SubSpline(polyhedron=P, function=f_0),
            plf.SubSpline(polyhedron=P, function=f_1)]


@pytest.fixture(autouse=True)
def sub_spline_with_strict_equality():
    c = [x <= 1, y == 1, x >= 0]
    P = ppl.C_Polyhedron(la.cs_from_list(c))
    f = rla.RationalLinearExpression((0, Fraction(-1, 4)), Fraction(1, 4))

    return plf.SubSpline(polyhedron=P, function=f)


@pytest.fixture(autouse=True)
def sub_spline_with_empty_poly_0():
    c = [x <= 1, x >= 2, y <= 1]
    P = ppl.C_Polyhedron(la.cs_from_list(c))
    f = rla.RationalLinearExpression((0, Fraction(-1, 4)), Fraction(1, 4))

    return plf.SubSpline(polyhedron=P, function=f)


@pytest.fixture(autouse=True)
def sub_spline_with_empty_poly_1():
    P = ppl.C_Polyhedron(3, 'empty')
    f = rla.RationalLinearExpression((0, Fraction(-1, 4)), Fraction(1, 4))

    return plf.SubSpline(polyhedron=P, function=f)


@pytest.fixture(autouse=True)
def sub_spline_0():
    P = ppl.C_Polyhedron(4, 'universe')
    f = rla.RationalLinearExpression((0, Fraction(-1, 4)), Fraction(1, 4))
    return plf.SubSpline(polyhedron=P, function=f)
