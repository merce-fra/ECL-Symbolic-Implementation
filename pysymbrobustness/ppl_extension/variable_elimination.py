# coding=utf-8
from __future__ import annotations  # For forward reference typing

from copy import deepcopy
from itertools import product, chain
from typing import Tuple, List

import numpy as np
import ppl as ppl
import pysymbrobustness.ppl_extension.linear_algebra as lin_alg
import pysymbrobustness.ppl_extension.piecewise_linear_function as plf
import pysymbrobustness.ppl_extension.rational_linear_algebra as rla
from pysymbrobustness.ppl_extension import polyhedra


def sorted_index(cs: ppl.Constraint_System, resets: List[int]) -> Tuple[
    List[int], List[int], List[int], List[int]]:
    negative_p = []
    zero_p = []
    positive_p = []

    sum_of_coefficients = []
    for i, cs_i in enumerate(cs):
        cs_coeff = list(cs_i.coefficients())

        for reset_index in resets:
            cs_coeff[reset_index] = 0

        res = sum(cs_coeff)
        sum_of_coefficients.append(res)
        if res > 0:
            positive_p.append(i)
        elif res < 0:
            negative_p.append(i)
        else:
            zero_p.append(i)
    return negative_p, zero_p, positive_p, sum_of_coefficients


def sorted_ind(le: List[ppl.Linear_Expression]) -> Tuple[List[int], ...]:
    negative = []
    positive = []
    zero = []
    sum_of_coeff = []

    for i, le_i in enumerate(le):

        le_coeff = list(le_i.coefficients())

        res = sum(le_coeff)
        sum_of_coeff.append(res)

        if res > 0:
            positive.append(i)
        elif res < 0:
            negative.append(i)
        else:
            zero.append(i)
    return negative, zero, positive, sum_of_coeff


def reset_le(le: ppl.Linear_Expression, resets: List[int]) \
        -> ppl.Linear_Expression:
    res = deepcopy(le)
    for index in resets:
        temp = le.coefficient(ppl.Variable(index)) * ppl.Variable(index)
        res = res - temp
    return res


def truncated_le_from_cs(poly: ppl.C_Polyhedron, resets: List[int]) -> \
        List[ppl.Linear_Expression]:
    return [reset_le(le=lin_alg.constraint_as_le(c), resets=resets)
            for c in poly.constraints()]


def entry_set(poly_p: ppl.C_Polyhedron, poly_q: ppl.C_Polyhedron, resets:
List[int], dim: int) -> ppl.C_Polyhedron:
    if poly_p.space_dimension() != dim or poly_q.space_dimension() != dim:
        raise Exception

    # Get the constraint system of each polyhedron and it's linear expressions

    le_P = truncated_le_from_cs(poly=poly_p, resets=resets)
    le_Q = truncated_le_from_cs(poly=poly_q, resets=resets)

    # Get Z, N, and P for poly_p and poly_q:

    negative_p, zero_p, positive_p, sum_p = sorted_ind(le=le_P)
    negative_q, zero_q, positive_q, sum_q = sorted_ind(le=le_Q)

    # Forms phi inequalities:

    phi_p = [le_P[i] >= 0 for i in chain(negative_p, zero_p)]
    phi_q = [le_Q[i] >= 0 for i in chain(negative_q, zero_q)]

    phi = phi_p + phi_q

    # Forms the lists of inequalities that use to depend on alpha and beta

    # Squeezing alpha
    list_p = []
    for i, j in product(negative_p, positive_p):
        list_p.append(le_P[i]*sum_q[j] - le_Q[j]*sum_p[i] >= 0)

    # Squeezing beta
    list_q = []
    for i, j in product(negative_q, positive_q):
        list_q.append(le_Q[i]*sum_q[j] - le_Q[j]*sum_q[i] >= 0)

    # Inequalities providing that the upper bound of beta is greater then
    # the lower bound of alpha

    list_q_p = []
    for i, j in product(negative_q, positive_p):
        list_q_p.append(le_Q[i] * sum_q[j] - le_Q[j] * sum_q[i] >= 0)


    # Constraint system

    cs = ppl.Constraint_System()

    for constraint in chain(phi, list_p, list_q, list_q_p):
        cs.insert(constraint)

    res_polyhedra = polyhedra.c_polyhedron_constructor(dimension=dim,
                                                       cs_list=list(cs))

    return res_polyhedra

#
# def entry_set(poly_p: ppl.C_Polyhedron, poly_q: ppl.C_Polyhedron,
#               resets: List[int], dim: int) -> ppl.C_Polyhedron:
#     # Check if P and Q share the same dimension set:
#     if poly_p.space_dimension() != poly_q.space_dimension():
#         raise Exception
#
#     # Get the constraint system of each polyhedron
#
#     cs_P = poly_p.constraints()
#     cs_Q = poly_q.constraints()
#
#     # Get N_p, Z_p, P_p
#
#     negative_p, zero_p, positive_p, sum_p = sorted_index(cs=cs_P,
#                                                          resets=resets)
#     negative_q, zero_q, positive_q, sum_q = sorted_index(cs=cs_Q,
#                                                          resets=resets)
#
#     # Forms the inequalities:
#
#     l_0 = [cs_P[i] for i in
#            negative_p]  # All the constraints of poly_p s.t. the reseted_sum of their coefficients <=0
#     l_1 = [cs_Q[i] for i in
#            negative_q]  # all the constraint  of poly_p s.t. the reseted_sum of their coefficients <=0
#     l_2 = []
#     l_3 = []
#     l_4 = []
#     phi = [cs_P[i] for i in zero_p] + [cs_Q[i] for i in zero_q] + \
#           [cs_P[i] for i in negative_q] + [cs_Q[i] for i in negative_q]
#     # Adding constraint that does not depends on alpha or beta
#
#     # For every index i in N_p and j in P_p, we know that -cs_P[i]/sum_p[i] >= -cs_P[j]/sum_p[j]
#     # and -cs_P[i]/sum_p[i] >=0
#     # Therefore, we had the constraint cs_P[i] >=0 and +cs_P[i]*sum_p[j] - cs_P[j]*sum_p[i] >= 0 (multiplying by
#     # sum_p[i] change the sign of the inequality)
#     for neg_p_index, pos_p_index in product(negative_p, positive_p):
#         le_phi_neg_p = lin_alg.constraint_as_le(cs_P[neg_p_index])
#         sum_neg_p = sum_p[neg_p_index]
#
#         le_phi_pos_p = lin_alg.constraint_as_le(cs_P[pos_p_index])
#         sum_pos_p = sum_p[pos_p_index]
#
#         l_2.append(le_phi_neg_p * sum_pos_p - le_phi_pos_p * sum_neg_p >= 0)
#
#     # For every index i in N_q, j in P_q and k in P_p, we know that:
#     # 1) -cs_Q[i]/sum_q[i] >= -cs_Q[j]/sum_q[j] (which gives +cs_Q[i]*sum_q[j] - cs_Q[j]*sum_q[i] >= 0)
#     # 2) -cs_Q[i]/sum_q[i] >= 0 (which gives cs_Q[i] >= 0)
#     # 3)
#     for neg_q_index, pos_q_index, pos_p_index in product(negative_q,
#                                                          positive_q,
#                                                          positive_p):
#         le_phi_neg_q = lin_alg.constraint_as_le(cs_Q[neg_q_index])
#         sum_neg_q = sum_q[neg_q_index]
#
#         le_phi_pos_p = lin_alg.constraint_as_le(cs_P[pos_p_index])
#         sum_pos_p = sum_p[pos_p_index]
#
#         le_phi_pos_q = lin_alg.constraint_as_le(cs_Q[pos_q_index])
#         sum_pos_q = sum_q[pos_q_index]
#
#         l_3.append(le_phi_neg_q * sum_pos_q - le_phi_pos_q * sum_neg_q >= 0)
#         l_4.append(
#             le_phi_pos_q * sum_pos_p - le_phi_pos_p * sum_pos_q <= 0)  # The lower bound of alpha should be lower
#         # Than the lower bound of beta
#
#     cs = ppl.Constraint_System()
#
#     for constraint in chain(l_0, l_1, l_2, l_3, l_4, phi):
#         cs.insert(constraint)
#
#     res_polyhedra = polyhedra.c_polyhedron_constructor(dimension=dim,
#                                                        cs_list=list(cs))
#     # res_polyhedra = ppl.C_Polyhedron(dim, 'universe')
#     # res_polyhedra.add_constraints(cs)
#
#     return res_polyhedra


def guard_poly(guard: ppl.C_Polyhedron, dim: int) -> ppl.C_Polyhedron:
    """

    @param guard:
    @type guard:
    @param dim:
    @type dim:
    @return:
    @rtype:
    """
    cs_g = guard.constraints()
    negative, zero, positive, sum = sorted_index(cs=cs_g, resets=[])

    cs_g_le = [lin_alg.constraint_as_le(c) for c in cs_g]
    phi = [cs_g_le[i] >= 0 for i in zero]

    l = [-cs_g_le[p] * sum[n] + cs_g_le[n] * sum[p] >= 0 for p, n in
         product(positive, negative)]
    up_guard = [cs_g_le[n] >= 0 for n in negative]
    cs = ppl.Constraint_System()
    for constraint in chain(phi, l, up_guard):
        cs.insert(constraint)

    res_polyhedra = polyhedra.c_polyhedron_constructor(dimension=dim,
                                                       cs_list=list(cs))
    return res_polyhedra


def alpha_beta_squeezing(poly_p: ppl.C_Polyhedron, poly_q: ppl.C_Polyhedron,
                         resets: List[int]) -> \
        Tuple[List[rla.RationalLinearExpression], List[
            rla.RationalLinearExpression],
              List[rla.RationalLinearExpression], List[
                  rla.RationalLinearExpression]]:
    """
    @param poly_p: a target polyhedron where the point v+alpha[r] is.
    @type poly_p: a ppl.C_polyhedron.
    @param poly_q: a target polyhedron where the point v+beta[r] is.
    @type poly_q: a ppl.C_polyhedron.
    @param resets: the list of resets.
    @type resets: a list of integer.
    @return: Four list of rational linear expressions, low_alpha, up_alpha, low_beta, up_beta such that for a fixed v,
    the alpha such that v+alpha[r] is in poly_p is the alpha such that low_alpha <= alpha <= up_alpha and the beta
    such that v+beta[r] is in poly_q is the beta such that low_beta <= beta <= up_beta.
    @rtype: a tuple of four list of rla.RationalLinearExpression.
    """
    # Check if P and Q share the same dimension set:
    dim = poly_p.space_dimension()

    if poly_p.space_dimension() != poly_q.space_dimension():
        raise Exception

    le_P = [le >= 0 for le in truncated_le_from_cs(poly=poly_p, resets=resets)]
    le_Q = [le >= 0 for le in truncated_le_from_cs(poly=poly_q, resets=resets)]

    reset_poly_P = polyhedra.c_polyhedron_constructor(dimension=dim, cs_list=le_P)
    reset_poly_Q = polyhedra.c_polyhedron_constructor(dimension=dim, cs_list=le_Q)

    # Get the constraint system of each polyhedron

    cs_P = reset_poly_P.constraints()
    cs_Q = reset_poly_Q.constraints()

    # Get N_p, Z_p, P_p

    negative_p, zero_p, positive_p, sum_p = sorted_index(cs=cs_P,
                                                         resets=resets)
    negative_q, zero_q, positive_q, sum_q = sorted_index(cs=cs_Q,
                                                         resets=resets)

    def inner_function(constraint: ppl.Constraint,
                       sum: int) -> rla.RationalLinearExpression:
        """
        @param constraint: a constraint in pplpy, associated with the linear expression le
        @type constraint: a ppl.Constraint
        @param sum: The sum of all the coefficients of the linear expression which coefficients values the same as the
        linear expression of constraint if the variable is not in the reset list, and 0 otherwise.
        @type sum: an int
        @return: the rational linear expression that values -le/sum
        @rtype: A rational linear expression
        """


        coefficients = tuple(
            - c for c in lin_alg.constraint_as_le(constraint).coefficients())
        inhomogeneous_term = - lin_alg.constraint_as_le(
            constraint).inhomogeneous_term()
        rle = rla.RationalLinearExpression(homogeneous_terms=coefficients,
                                           inhomogeneous_term=inhomogeneous_term,
                                           div=sum)
        return rle

    zero_function = rla.RationalLinearExpression(homogeneous_terms=0,
                                                 inhomogeneous_term=0,
                                                 div=1) # To add that alpha and beta should be superior than 0

    low_alpha = [inner_function(constraint=cs_P[i], sum=sum_p[i]) for i in
                 positive_p]
    low_alpha.append(zero_function)
    low_beta = [inner_function(constraint=cs_Q[i], sum=sum_q[i]) for i in
                positive_q]
    low_beta.append(zero_function)
    up_alpha = [inner_function(constraint=cs_P[i], sum=sum_p[i]) for i in
                negative_p]
    up_beta = [inner_function(constraint=cs_Q[i], sum=sum_q[i]) for i in
               negative_q]

    return low_alpha, up_alpha, low_beta, up_beta


def maximal_intervals(poly_p: ppl.C_Polyhedron, poly_q: ppl.C_Polyhedron,
                      resets: List[int],
                      guard: ppl.C_Polyhedron, dim: int) -> Tuple[
    plf.Spline, plf.Spline, plf.Spline, plf.Spline]:
    """
    @param dim:
    @type dim:
    @param poly_p: A target polyhedron where the point v+alpha[r] is.
    @type poly_p: a pplpy closed polyhedron ppl.C_Polyhedron
    @param poly_q: A target polyhedron where the point v+beta[r] is.
    @type poly_q: a pplpy closed polyhedron ppl.C_Polyhedron
    @param resets: the list the reset of the variables.
    @type resets: a list of integers, possibly empty.
    @param guard: the polyhedron that represents every valuation that fits the guard.
    @type guard: a pplpy closed polyhedron ppl.C_Polyhedron
    @return: The two set of bounds min_alpha, max_alpha, min_beta, max_beta such that I_alpha = [min_alpha, max_alpha]
    and I_beta = [min_beta, max_beta].
    I_alpha represents, for a fixed configuration (l,v), the intervals of possible alpha such that v+alpha[r] is in
    poly_p.
    Same for I_beta but replace alpha by beta and poly_p by poly_q.
    @rtype: a tuple of four Spline.
    """
    # Entry polyhedron set

    entry_poly = entry_set(poly_p=poly_p, poly_q=poly_q, resets=resets,
                           dim=dim)
    entry_poly_base = guard_poly(guard=guard, dim=dim)
    entry_poly.intersection_assign(entry_poly_base)

    # From the target polyhedrons

    low_alpha, up_alpha, low_beta, up_beta = alpha_beta_squeezing(
        poly_p=poly_p, poly_q=poly_q, resets=resets)

    # From the guards

    low_alpha_g, up_alpha_g, low_beta_g, up_beta_g = alpha_beta_squeezing(
        poly_p=guard, poly_q=guard, resets=[])

    def inner_function(eq_guard: List[rla.RationalLinearExpression],
                       eq: List[rla.RationalLinearExpression],
                       is_min: bool) -> plf.Spline:
        """
        Apply the minimum or maximum to a list of linear expression, with the entry_poly as a definition entry set.
        @param eq_guard: a list of linear expression with rational coefficients
        @type eq_guard: a list of relation linear expressions
        @param eq: a list of linear expression with rational coefficients
        @type eq: a list of relation linear expressions
        @param is_min: a boolean that is true is the minimum must be applied and False if the maximum must be applied
        @type is_min: a bool
        @return: The minimum (or maximum if is_min is False) between the two list of spline with the same entry set
        entry_poly and with the linear expressions eq_guard, eq, and + infinity (or 0 if is_min is False)
        @rtype: a spline.
        """
        list_splines = list(
            map(lambda f: plf.Spline(
                plf.SubSpline(polyhedron=entry_poly, function=f)),
                chain(eq, eq_guard)))

        fct_res = rla.RationalLinearExpression(
            inhomogeneous_term=0) if not is_min else rla.InfiniteExpression(
            sign=True)
        res = plf.Spline(plf.SubSpline(polyhedron=entry_poly,
                                       function=fct_res))
        if is_min:
            res = res.minimum_list(others=list_splines)
        else:
            res = res.maximum_list(others=list_splines)
        return res

    min_alpha = inner_function(eq_guard=low_alpha_g, eq=low_alpha,
                               is_min=False)
    min_beta = inner_function(eq_guard=low_beta_g, eq=low_beta, is_min=False)
    max_alpha = inner_function(eq_guard=up_alpha_g, eq=up_alpha, is_min=True)

    max_beta = inner_function(eq_guard=up_beta_g, eq=up_beta, is_min=True)

    return min_alpha, max_alpha, min_beta, max_beta
