from functools import reduce
from typing import List, Optional, Tuple

import ppl


def cs_from_list(constraint_list: Optional[List[ppl.Constraint]]) -> \
        ppl.Constraint_System:
    """
    Input: a list of pplpy constraints
    Output: a pplpy constraint system, that contains all the constraints of
    the list of constraints
    -----
    Example:
        If the list is constraint_list = [ x - y <=0 , y >= 45, z <= 4],
        the output will be the constraint system that contains exactly
        x - y <=0 , y >= 45, z <= 4
    """
    cs = ppl.Constraint_System()
    for constraint in constraint_list:
        cs.insert(constraint)
    return cs


def cs_linear_expressions(constraint_cs: ppl.Constraint_System) -> \
        List[ppl.Linear_Expression]:
    """
    Input: a pplpy constraint system
    Output: a list of the pplpy linear expression associated (in terms of
    coefficients)
    -----
    Example on one constraint:
    c = [ x - y >=0 ], the associated linear_expression will be [x-y]
    c = [ x - y >= 3], the associated linear expression will be [x-y-3]
    c = [x - y <= -4], the associated linear expression will be [-x+y-4]
    c = [x - z >= 0], the associated linear expression will be [x-z]
    """
    L = list(map(constraint_as_le, constraint_cs))
    return L


def constraint_as_le(c: ppl.Constraint) -> ppl.Linear_Expression:
    """
        Input: a pplpy constraint
        Output: a pplpy linear expression associated (in terms of
        coefficients)
        -----
        Example on one constraint:
        c = x - y >=0 , the associated linear_expression will be [x-y]
        c = x - y >= 3, the associated linear expression will be [x-y-3]
        c = x - y <= -4, the associated linear expression will be [-x+y-4]
        c = x - z >= 0, the associated linear expression will be [x-z]
        """
    return reduce(lambda x, xi: x + c.coefficient(xi) * xi,
                  map(lambda i: ppl.Variable(i), range(c.space_dimension())),
                  c.inhomogeneous_term())


def non_strict_inverse(constraint: ppl.Constraint):
    le = constraint_as_le(constraint)
    return le <= 0


def equality(le_0: ppl.Linear_Expression, le_1:ppl.Linear_Expression) -> bool:
    """
    returns True if the linear_expression linear_expression_0 =
    linear_expression_1. False otherwise.
    """
    if le_0.space_dimension() != le_1.space_dimension():
        return False

    if le_0.inhomogeneous_term() != le_1.inhomogeneous_term():
        return False

    for i, coefficient in enumerate(le_1.coefficients()):
        v = ppl.Variable(i)
        if coefficient != le_0.coefficient(v):
            return False

    return True
