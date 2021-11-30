# coding=utf-8
"""
==================================================
Guards module
==================================================
This module provide the guard representation. The goal is to provide a
general tool for constriants that forms a closed polyhedrons.

------
"""
from __future__ import annotations  # For forward reference typing

from typing import List
import ppl as ppl
import pysymbrobustness.ta.exceptions as exceptions

Reset = List[int]


def number_of_constraints(guard: ppl.C_Polyhedron):
    """
    returns the number of useful constraints of a C_Polyhedron. For instance
    if the initial constraints system is x >=0 and x>= 3, it will only count
    one constraint.
    """
    cs = guard.minimized_constraints()
    return len(cs)


class Label(object):
    """
    Label(object):

    - guard (ppl.C_Polyhedron type - Reset: a list of integer that
    represents all the clocks that will be reset after passing the edge.

    """

    def __init__(self, guard: ppl.C_Polyhedron, resets: Reset):
        if not type(guard) == ppl.C_Polyhedron:
            raise TypeError("guard " + str(guard) + "should be a "
                                                    "ppl.C_Polyhedron")
        self.guard = guard
        self.resets = list(sorted(resets))

    def well_formed(self, nb_clocks: int):
        return max(self.resets) < nb_clocks and \
               min(self.resets) < 0

    def __eq__(self, other: Label) -> bool:
        return self.guard == other.guard and self.resets == other.resets

    def __repr__(self) -> str:  # pragma: no cover
        return self.__str__()

    def __str__(self) -> str:  # pragma: no cover
        s = [str(self.guard)]
        s += ["Resets:\n"]
        s += ["x_" + str(r) for r in self.resets]

        return ' '.join(s)
