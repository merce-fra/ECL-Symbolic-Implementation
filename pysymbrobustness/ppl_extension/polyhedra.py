# coding=utf-8
"""
==================================================
polyhedra module
==================================================
Classes:
------
Methods:
------
"""
from __future__ import annotations  # For forward reference typing

from typing import List, Tuple, Optional, Union
import ppl as ppl
import numpy as np
import pysymbrobustness.ppl_extension.linear_algebra as la
import pyny3d.geoms as pyny


def c_polyhedron_constructor(dimension: int, cs_list: Optional[List[ppl.Constraint]] = None) \
        -> ppl.C_Polyhedron:
    """
    Input:
    :param dimension: integer
    :param cs_list: a list of pplpy constraint
    Output:
    return a pplpy closed polyhedra (ppl.C_Polyhedron) in QQ^n where n =
    dimension and that contains exactly all the constraints of cs_list.
    -----
    Examples:
    TODO
    """
    if cs_list is None:
        cs_list = []

    P = ppl.C_Polyhedron(dimension, 'universe')
    for i in range(dimension):
        P.add_constraint(ppl.Variable(i) >= 0)
    for constraint in cs_list:
        P.add_constraint(constraint)
    return P


def polyhedron_interior(poly: Union[ppl.C_Polyhedron, ppl.NNC_Polyhedron]) -> \
        ppl.NNC_Polyhedron:
    """
    Returns the interior polyhedron of a non necessary (or a closed) polyhedron.
    Input:
    :param poly: A closed or a non-necessary closed polyhedron from pplpy
    Output:
    A non-necessary closed polyhedron from pplpy.
    If poly does not contains any equality constraints, the output polyhedron
    constains all the strict inequalities from poly's
    constraints systems and the non-strict inequalities, but as strict
    inequalities.
    Examples:
    P.constraints() = {x0>=0, -x0+1>=0, x1>0, -x1+1==0, -x2+1>=0, x2>=0}
    polyhedron_interior(P) = NNC_Polyhedron(P.space_dimension(), 'empty')

    Q.constraints() = {x0>=0, -x0+1>=0, x1>0, -x1+1>=0, -x2+1>=0, x2>=0}
    polyhedron_interior(P).constraints() = {x0>0, -x0+1>0, x1>0, -x1+1>0,
    -x2+1>0, x2>0}

    """
    initial_cs = poly.minimized_constraints()
    interior_cs = ppl.Constraint_System()
    for constraint in initial_cs:
        if not constraint.is_tautological() or not constraint.is_inconsistent():
            if constraint.is_equality():
                return ppl.NNC_Polyhedron(poly.space_dimension(), 'empty')
            # Indeed, le == 0 is equivalent to le >= 0 and le <=0 , the constraint that should then be add
            # is le >0 and le <0, which gives directly an empty polyhedron
            elif constraint.is_nonstrict_inequality():
                interior_cs.insert(la.constraint_as_le(constraint) > 0)
            elif constraint.is_strict_inequality():
                interior_cs.insert(constraint)

    interior_poly = ppl.NNC_Polyhedron(poly.space_dimension(), 'universe')
    interior_poly.add_constraints(interior_cs)
    return interior_poly


def disjoint_by_interior(P: Union[ppl.NNC_Polyhedron, ppl.C_Polyhedron],
                         Q: Union[ppl.NNC_Polyhedron, ppl.C_Polyhedron]) -> \
        bool:
    """
    Returns true if the interior of the two input polyhedrons are disjoints.
    Input:
    param: P: a ppl.Polyhedron
    param: Q a ppl.Polyhedron
    Output:
    param: a boolean, True if P and Q are disjoints, or if their interiors are. False otherwise.
    """
    return P.is_disjoint_from(Q) or \
           polyhedron_interior(P).is_disjoint_from(polyhedron_interior(Q))


def contains_pure_equalities(P: Union[ppl.NNC_Polyhedron, ppl.C_Polyhedron]) \
        -> bool:

    return any(constraint.is_equality() for constraint in P.constraints())


def rank(P: Union[ppl.NNC_Polyhedron, ppl.C_Polyhedron]) -> int:
    """
    returns the number of useful constraints of a C_Polyhedron. For instance
    if the initial constraints system is x >=0 and x>= 3, it will only count
    one constraint.
    """
    cs = P.minimized_constraints()
    return len(cs)


def show(P: Union[ppl.NNC_Polyhedron, ppl.C_Polyhedron]):
    if (P.space_dimension() == 2 or P.space_dimension() == 3) and all(g.is_point() for g in P.generators()):
        gens = [g.coefficients() for g in P.generators()]
        polyhedron = pyny.Polygon(np.array(gens))
        polyhedron.plot('b')
    else:
        raise NotImplementedError("show only implemented for 2 or 3-dimension polyhedron")
