# coding=utf-8
import pytest
import ppl as ppl
import pysymbrobustness.ppl_extension.polyhedra as polyhedra
from tests.test_polyhedra_ex import *


class TestPolyhedraMethods:
    def test_c_polyhedron_constructor(self):
        P_0 = polyhedra.c_polyhedron_constructor(4, [])
        assert P_0 == ppl.C_Polyhedron(4, 'universe')

        x = ppl.Variable(0)
        y = ppl.Variable(1)
        cs_1 = [x - y >= 2, x <= 4]
        P_0 = polyhedra.c_polyhedron_constructor(2, cs_1)

        cs_1_bis = ppl.Constraint_System()
        cs_1_bis.insert(x - y >= 2)
        cs_1_bis.insert(x <= 4)
        assert P_0 == ppl.C_Polyhedron(cs_1_bis)

    def test_interior_polyhedron(self, polyhedron_open_0_associated, polyhedron_closed_0,
                                 polyhedron_with_equalities_0, polyhedron_with_equalities_1,
                                 polyhedron_with_equalities_2, polyhedron_with_strict_inequalities,
                                 polyhedron_with_strict_inequalities_open):
        assert polyhedra.polyhedron_interior(polyhedron_closed_0) == polyhedron_open_0_associated
        assert polyhedra.polyhedron_interior(polyhedron_with_equalities_0) == ppl.NNC_Polyhedron(3, 'empty')
        assert polyhedra.polyhedron_interior(polyhedron_with_equalities_2) == ppl.NNC_Polyhedron(5, 'empty')
        assert polyhedra.polyhedron_interior(polyhedron_with_strict_inequalities) == \
               polyhedron_with_strict_inequalities_open

    def test_contains_pure_equalities(self, polyhedron_with_equalities_0, polyhedron_with_equalities_1,
                                      polyhedron_with_equalities_2, polyhedron_closed_0):
        assert polyhedra.contains_pure_equalities(polyhedron_with_equalities_0)
        assert polyhedra.contains_pure_equalities(polyhedron_with_equalities_1)
        assert polyhedra.contains_pure_equalities(polyhedron_with_equalities_2)
        assert not polyhedra.contains_pure_equalities(polyhedron_closed_0)
