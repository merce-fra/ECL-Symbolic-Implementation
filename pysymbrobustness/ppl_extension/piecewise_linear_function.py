# coding=utf-8
"""
==================================================
piecewise linear function module
==================================================


Classes:
------


Methods:
------
"""

from __future__ import annotations  # For forward reference typing

import itertools
from copy import deepcopy
from typing import List, Union, Tuple, Iterable
import ppl as ppl
import pysymbrobustness.ppl_extension.rational_linear_algebra as rla
import pysymbrobustness.ppl_extension.linear_algebra as linear_alg
from pysymbrobustness.ppl_extension import polyhedra

LinearFunction = Union[rla.InfiniteExpression,
                       rla.RationalLinearExpression,
                       rla.AbstractExpression]


class SubSpline(object):
    __slots__ = ["function", "polyhedron"]

    def __init__(self, function: LinearFunction,
                 polyhedron: ppl.C_Polyhedron):
        if type(function) != rla.RationalLinearExpression \
                and type(function) != rla.InfiniteExpression \
                and type(function) != rla.AbstractExpression:
            raise TypeError(str(function) + " should be either a " +
                            "RationalLinearExpression or an InfiniteExpression")
        self.function = function
        self.polyhedron = polyhedron

    def _minmax(self, other: SubSpline, isMax: bool) -> Spline:
        """

        @param other:
        @type other:
        @param isMax:
        @type isMax:
        @return:
        @rtype:
        """
        final_spline = Spline()
        f0 = self.function
        f1 = other.function

        p0 = self.polyhedron
        p1 = other.polyhedron

        n = self.polyhedron.space_dimension()

        if polyhedra.disjoint_by_interior(P=p0, Q=p1):
            return Spline(SubSpline(rla.RationalLinearExpression(),
                                     ppl.C_Polyhedron(n, 'empty')))

        cs_f0 = ppl.Constraint_System()
        cs_f1 = ppl.Constraint_System()

        for constraint in itertools.chain(p0.constraints(), p1.constraints()):
            cs_f0.insert(constraint)
            cs_f1.insert(constraint)

        cs_f0.insert(f1 <= f0 if isMax else f0 <= f1)
        cs_f1.insert(f0 <= f1 if isMax else f1 <= f0)

        final_spline.add_sub_spline(SubSpline(f0, polyhedra.c_polyhedron_constructor(dimension=n, cs_list=list(cs_f0))))
        final_spline.add_sub_spline(SubSpline(f1, polyhedra.c_polyhedron_constructor(dimension=n, cs_list=list(cs_f1))))

        return final_spline

    def minimum(self: SubSpline, other: SubSpline) -> Spline:
        """
        Input:
        :param other: SubSpline
        Output:
        a Spline, that represents the minimum between the two SubSplines
        self and piece1 over the intersection of their two polyhedron entry
        definition sets
        -----
        Example:
        self = (x>=0 ,x - y)
        piece1 (x>=3, 1 + x)
        output: the Spline [( y<= 1, 1+x), (y>= 1, x-y)]
        """
        return self._minmax(other, isMax=False)

    def maximum(self: SubSpline, other: SubSpline) -> Spline:

        """
        Input:
        :param other: SubSpline
        Output:
        a Spline, that represents the maximum between the two SubSplines
        self and piece1 over the intersection of their two polyhedron entry
        definition sets
        -----
        Example:
        self = (x>=0 ,x - y)
        piece1 (x>=3, 1 + x)
        Output: the Spline [( y>= 1, 1+x), (y<= 1, x-y)]
        """
        return self._minmax(other, isMax=True)

    def is_equal_to(self, other: SubSpline):
        if type(self.function) == type(other.function):
            return self.polyhedron == other.polyhedron and self.function.is_equal_to(other.function)

        return False


class Spline(object):
    __slots__ = ["sub_splines"]

    def __init__(self, sub_spline_list: Union[Iterable[SubSpline], SubSpline] = None):
        """

        """

        self.sub_splines = []
        if sub_spline_list is None:
            return

        if type(sub_spline_list) == SubSpline:
            self.add_sub_spline(sub_spline_list)
            return

        for sub in sub_spline_list:
            self.add_sub_spline(sub)

    def add_sub_spline(self, sub_spline: SubSpline):
        """

        @param sub_spline:
        @type sub_spline:
        @return:
        @rtype:
        """

        if type(sub_spline) != SubSpline:
            raise TypeError
        elif polyhedra.contains_pure_equalities(sub_spline.polyhedron):
            pass  # It is just a frontier or a point
        elif any(sub_spline.is_equal_to(self_sub_spline) for self_sub_spline in self.sub_splines):
            pass  # This sub-spline already exists
        elif sub_spline.polyhedron.is_empty():
            pass  # Empty polyhedron, empty sub-spline
        elif type(sub_spline.function) == rla.InfiniteExpression and not sub_spline.function.is_positive:
            pass  # - infinite function, useless to add.
        else:
            self.sub_splines.append(sub_spline)

    def pop_sub_spline(self, index: int):
        self.sub_splines.pop(index)

    def is_equal_to(self, other: Spline):
        """
        Verifies if two spline have the same partition system, with the same function labeled
        """
        for sub_self in self.sub_splines:
            if not any(sub_self.is_equal_to(other=s) for s in other.sub_splines):
                return False

        for sub_other in other.sub_splines:
            if not any(sub_other.is_equal_to(other=s) for s in self.sub_splines):
                return False

        return True

    def fusion(self, other: Spline):
        for sub in other.sub_splines:
            self.add_sub_spline(sub)

    @property
    def functions(self) -> List[LinearFunction]:
        return [sub.function for sub in self.sub_splines]

    @property
    def polyhedrons(self) -> List[ppl.C_Polyhedron]:
        return [sub.polyhedron for sub in self.sub_splines]

    def partition(self, other: Spline) -> Spline:
        """
        Return the same Spline as self, but over-partitioned with other's subspline's polyhedron constraints.
        Warning: this function only works on Spline that are perfect partition with closed polyhedrons.
        """
        self_new_spline = Spline()
        other_constraints = Spline.list_const_and_sym_constr(spline=other)

        for sub_spline in self.sub_splines:
            new_spline = Spline.sub_partition(
                constraints_list=other_constraints,
                polyhedron_list=[sub_spline.polyhedron],
                label=sub_spline.function)

            for new_sub_spline in new_spline.sub_splines:
                self_new_spline.add_sub_spline(sub_spline=new_sub_spline)

        return self_new_spline

    @staticmethod
    def list_const_and_sym_constr(spline: Spline) -> List[Tuple[ppl.Constraint, ppl.Constraint]]:
        """
        return all the constraints of a spline,
        if the spline is a perfect partition of closed polyhedron.
        These constraint are sorted with grouped constraints of the form (le>=0 , le<=0)
        """
        constraint_list = []

        for sub in spline.sub_splines:
            constraints = sub.polyhedron.constraints()

            for constraint in constraints:
                le = linear_alg.constraint_as_le(constraint)
                c_0 = (le >= 0)
                c_1 = (le <= 0)
                # not_add_c_0 = Exists c in constraint_list tq c0 <=> c
                not_add_c_0 = any(
                    c_0.is_equivalent_to(c[0]) or c_0.is_equivalent_to(c[1])
                    for c in constraint_list
                )
                if not not_add_c_0:
                    constraint_list.append((c_0, c_1))

        return constraint_list

    @staticmethod
    def sub_partition(
            constraints_list: List[Tuple[ppl.Constraint, ppl.Constraint]],
            polyhedron_list: List[ppl.C_Polyhedron],
            label: LinearFunction):
        for constraint_0, constraint_1 in constraints_list:
            to_replace_polyhedron = []
            for poly in polyhedron_list:
                if poly.relation_with(constraint_0).implies(
                        ppl.Poly_Con_Relation.strictly_intersects()):
                    new_poly_0 = deepcopy(poly)
                    new_poly_1 = deepcopy(poly)
                    new_poly_0.add_constraint(constraint_0)
                    new_poly_1.add_constraint(constraint_1)
                    to_replace_polyhedron.append(new_poly_0)
                    to_replace_polyhedron.append(new_poly_1)
                else:
                    to_replace_polyhedron.append(poly)
            polyhedron_list = to_replace_polyhedron

        sub = [SubSpline(polyhedron=poly, function=label) for poly in polyhedron_list]
        return Spline(sub)

    @staticmethod
    def infinite(spline_0: Spline, spline_1: Spline, final_spline: Spline, op_is_min: bool):
        to_remove_index = []
        for index, sub_spline in enumerate(spline_0.sub_splines):
            if all(polyhedra.disjoint_by_interior(P=sub_spline.polyhedron, Q=poly)
                   for poly in spline_1.polyhedrons):
                if not op_is_min:
                    final_spline.add_sub_spline(sub_spline=sub_spline)

                to_remove_index.append(index)

        to_keep = []
        next_remove = 0
        for i, s in enumerate(spline_0.sub_splines):
            # Warning: was <=, check...
            if next_remove < len(to_remove_index) and i == to_remove_index[next_remove]:
                next_remove += 1
                continue
            to_keep.append(s)
        spline_0.sub_splines = to_keep

    def operator(self, other: Spline, op_is_min: bool):
        final_spline = Spline()
        self_partitioned_spline = self.partition(other=other)
        other_partitioned_spline = other.partition(other=self)

        # Dealing with disjoint sub_spline (both spline)

        Spline.infinite(spline_0=self_partitioned_spline, spline_1=other_partitioned_spline,
                        final_spline=final_spline, op_is_min=op_is_min)
        Spline.infinite(spline_0=other_partitioned_spline, spline_1=self_partitioned_spline,
                        final_spline=final_spline, op_is_min=op_is_min)

        # Dealing with the rest of sub_spline
        # Now the union of the polyhedrons of the sub_splines of self_partitioned_spline and other_partitioned_spline is
        # the same

        for sub_spline_0, sub_spline_1 in itertools.product(self_partitioned_spline.sub_splines,
                                                            other_partitioned_spline.sub_splines):
            res = sub_spline_0.minimum(sub_spline_1) if op_is_min else sub_spline_0.maximum(sub_spline_1)

            for sub_spline in res.sub_splines:
                final_spline.add_sub_spline(sub_spline)

        return final_spline

    # TODO: test if the two spline are - math.inf value...
    def maximum(self, other: Spline):
        return self.operator(other=other, op_is_min=False)

    def maximum_list(self, others: Iterable[Spline]):
        current = deepcopy(self)
        for other in others:
            current = current.maximum(other=other)
        return current

    @staticmethod
    def maximum_of_splines(spline_list: List[Spline]):
        if len(spline_list) == 0:
            return Spline()
        initial_spline = spline_list.pop()
        return initial_spline.maximum_list(others=spline_list)

    def minimum(self, other: Spline):
        return self.operator(other=other, op_is_min=True)

    def minimum_list(self, others: Iterable[Spline]):
        current = deepcopy(self)
        for other in others:
            current = current.minimum(other=other)
        return current

    @staticmethod
    def minimum_of_splines(spline_list: List[Spline]):
        initial_spline = spline_list.pop()
        return initial_spline.minimum_list(others=spline_list)