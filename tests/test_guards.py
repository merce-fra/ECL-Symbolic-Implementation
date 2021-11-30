# coding=utf-8
import pytest
import ppl as ppl
from pysymbrobustness.ta.guards import number_of_constraints
from tests.test_guards_ex import *


class TestGuardsMethods:
    def test_number_of_constraints(self, empty_cs_polyhedra, len_1_cs_polyhedra, len_1_cs_again_polyhedra,
                                   len_2_cs_polyhedra, len_3_cs_polyhedra):
        assert number_of_constraints(empty_cs_polyhedra) == 0
        assert number_of_constraints(len_1_cs_polyhedra) == 1
        assert number_of_constraints(len_1_cs_again_polyhedra) == 1
        assert number_of_constraints(len_2_cs_polyhedra) == 2
        assert number_of_constraints(len_3_cs_polyhedra) == 3

    def test_eq(self, len_1_cs_polyhedra, len_1_cs_again_polyhedra):
        assert len_1_cs_polyhedra == len_1_cs_again_polyhedra
