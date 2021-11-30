# coding=utf-8
from __future__ import annotations  # For forward reference typing

from copy import copy, deepcopy
from itertools import product
from typing import Tuple, List, Union

import numpy as np
import ppl as ppl
import pysymbrobustness.ppl_extension.linear_algebra as constraints
import pysymbrobustness.ppl_extension.rational_linear_algebra as rla
import pysymbrobustness.ppl_extension.piecewise_linear_function as plf
import pysymbrobustness.ppl_extension.linear_algebra as cstr
from fractions import Fraction


def _get_coefficient(
        f: Union[rla.RationalLinearExpression, rla.InfiniteExpression],
        alpha: ppl.Variable) -> Tuple[Fraction, rla.RationalLinearExpression, rla.RationalLinearExpression]:
    if type(f) == rla.RationalLinearExpression:
        # Getting a

        a = Fraction(f.coefficient(alpha))  # type=fraction

        # Getting the RLA a*alpha (alpha_term)

        alpha_id = alpha.id()

        alpha_list = [0 for i in range(alpha_id)]
        alpha_list.append(a)

        alpha_term = rla.RationalLinearExpression(tuple(alpha_list))

        # Getting b (type = RLA)
        b = f - alpha_term

        return a, alpha_term, b


def g_inf_min(a: Fraction, b: rla.RationalLinearExpression,
              alpha_opt: rla.RationalLinearExpression, beta_opt: rla.RationalLinearExpression,
              entry_polyhedron: ppl.C_Polyhedron) -> plf.Spline:
    h = plf.Spline(plf.SubSpline(polyhedron=entry_polyhedron,
                                 function=beta_opt - alpha_opt))
    f = plf.Spline(plf.SubSpline(polyhedron=entry_polyhedron, function=a * alpha_opt + b))
    return f.minimum(h)
    # return plf.minimum_spline([f, h])


def f_inf_min(c: Fraction, d: rla.RationalLinearExpression,
              alpha_opt: rla.RationalLinearExpression, beta_opt: rla.RationalLinearExpression,
              entry_polyhedron: ppl.C_Polyhedron) -> plf.Spline:
    h = plf.Spline(plf.SubSpline(polyhedron=entry_polyhedron,
                                 function=beta_opt - alpha_opt))
    g = plf.Spline(plf.SubSpline(polyhedron=entry_polyhedron, function=c * beta_opt + d))
    return g.minimum(h)
    # return plf.minimum_spline([g, h])


def f_g_min(a: Fraction, b: rla.RationalLinearExpression, c: Fraction, d: rla.RationalLinearExpression,
            alpha_opt: rla.RationalLinearExpression, beta_opt: rla.RationalLinearExpression,
            entry_polyhedron: ppl.C_Polyhedron) -> plf.Spline:
    h = plf.Spline(plf.SubSpline(polyhedron=entry_polyhedron,
                                 function=beta_opt - alpha_opt))
    f = plf.Spline(plf.SubSpline(polyhedron=entry_polyhedron, function=a * beta_opt + b))
    g = plf.Spline(plf.SubSpline(polyhedron=entry_polyhedron, function=c * beta_opt + d))
    return f.minimum_list([g, h])


def plf_optimization(
        f: Union[rla.RationalLinearExpression, rla.InfiniteExpression],
        g: Union[rla.RationalLinearExpression, rla.InfiniteExpression],
        entry_polyhedron: ppl.C_Polyhedron,
        alpha: ppl.Variable,
        beta: ppl.Variable,
        m_inf_alpha: rla.RationalLinearExpression,
        m_sup_alpha: Union[
            rla.RationalLinearExpression, rla.InfiniteExpression],
        m_inf_beta: rla.RationalLinearExpression,
        m_sup_beta: Union[
            rla.RationalLinearExpression, rla.InfiniteExpression]):
    """
    D = [m_a, M_a] x [m_b, M_b], all RLE or infinite
    f = a*alpha + b where b is a RLE or f = + inf, or - inf (Infinite function)
    g = c*beta + d where d is a RLE or g = + inf or - inf (Infinite function)

    """
    f_is_inf, g_is_inf = False, False
    beta_opt = 0
    alpha_opt = 0
    perm_value = 0

    # Getting the coefficient or the infinite boolean of f

    if type(f) == rla.RationalLinearExpression:
        a, alpha_term, b = _get_coefficient(f, alpha)
    elif type(f) == rla.InfiniteExpression:
        f_is_inf = True
    else:
        raise TypeError(str(
            f) + "should be either an InfiniteExpression or a RationalLinearExpression")

    # Getting the coefficient or the infinite boolean of f

    if type(g) == rla.RationalLinearExpression:
        c, beta_term, d = _get_coefficient(f=g, alpha=beta)
    elif type(g) == rla.InfiniteExpression:
        g_is_inf = True
    else:
        raise TypeError(str(
            g) + "should be either an InfiniteExpression or a RationalLinearExpression")

    # Dealing with infinite intervals
    # TODO: deal with the evaluation of the permissiveness when beta = + inf
    if not f_is_inf and not g_is_inf and type(m_sup_beta) == rla.InfiniteExpression and type(
            m_sup_alpha) == rla.InfiniteExpression:
        return f_g_min(a=a, b=b, c=c, d=d,
                       alpha_opt=m_inf_alpha,
                       beta_opt=m_sup_beta,
                       entry_polyhedron=entry_polyhedron)
    elif not f_is_inf and not g_is_inf and type(m_sup_beta) == rla.InfiniteExpression and type(
            m_sup_alpha) == rla.RationalLinearExpression:
        if a >= 0:
            return f_g_min(a=a, b=b, c=c, d=d,
                           alpha_opt=m_sup_alpha,
                           beta_opt=m_sup_beta,
                           entry_polyhedron=entry_polyhedron)
        elif a <= 0:
            return f_g_min(a=a, b=b, c=c, d=d,
                           alpha_opt=m_inf_alpha,
                           beta_opt=m_sup_beta,
                           entry_polyhedron=entry_polyhedron)
    elif not f_is_inf and not g_is_inf and type(m_sup_beta) == rla.RationalLinearExpression and type(
            m_sup_alpha) == rla.InfiniteExpression:
        if c >= 0:
            return f_g_min(a=a, b=b, c=c, d=d,
                           alpha_opt=m_inf_alpha,
                           beta_opt=m_sup_beta,
                           entry_polyhedron=entry_polyhedron)
        elif c <= 0:
            pass
            # Copy another case...

    # When the intervals are not infinite: dealing with infinite permissiveness
    spline = plf.Spline()
    if f_is_inf and g_is_inf:  # alpha = m_a, Beta = M_b, perm_value = M_b - m_a
        alpha_opt = m_inf_alpha
        beta_opt = m_sup_beta
        sub_spline = plf.SubSpline(polyhedron=entry_polyhedron,
                                   function=beta_opt - alpha_opt)
        spline.add_sub_spline(sub_spline)
    elif f_is_inf and not g_is_inf:  # alpha = m_a

        if c >= 0:

            # beta = M_b, perm_value = min(M_b - m_a, cM_b +d)
            return f_inf_min(c=c, d=d,
                             alpha_opt=m_inf_alpha,
                             beta_opt=m_sup_beta,
                             entry_polyhedron=entry_polyhedron)

        elif c <= 0:

            middle_rle = (d + m_inf_alpha) * (Fraction(1, 1 - c))

            # If (d+m_a)/(1-c) <= max(m_a,m_b) and m_a <= m_b ->
            # beta = max(m_a,m_b) = m_b

            constraints_0 = [
                middle_rle <= m_inf_alpha,
                middle_rle <= m_inf_beta,
                m_inf_alpha <= m_inf_beta
            ]
            polyhedra_0 = deepcopy(entry_polyhedron)
            polyhedra_0.add_constraints(cstr.cs_from_list(constraints_0))
            sub_spline_0 = f_inf_min(c=c, d=d,
                                     alpha_opt=m_inf_alpha,
                                     beta_opt=m_inf_beta,
                                     entry_polyhedron=polyhedra_0)

            # If (d+m_a)/(1-c) <= max(m_a,m_b) and m_a >= m_b ->
            # beta = max(m_a,m_b) = m_a

            constraints_1 = [
                middle_rle <= m_inf_alpha,
                middle_rle <= m_inf_beta,
                m_inf_alpha >= m_inf_beta
            ]

            polyhedra_1 = deepcopy(entry_polyhedron)
            polyhedra_1.add_constraints(cstr.cs_from_list(constraints_1))

            sub_spline_1 = f_inf_min(c=c, d=d,
                                     alpha_opt=m_inf_alpha,
                                     beta_opt=m_inf_alpha,
                                     entry_polyhedron=polyhedra_1)

            # If (d+m_a)/(1-c) >= M_b -> beta = M_b
            constraints_2 = [middle_rle >= m_sup_beta]

            polyhedra_2 = deepcopy(entry_polyhedron)
            polyhedra_2.add_constraints(cstr.cs_from_list(constraints_2))
            sub_spline_2 = f_inf_min(c=c, d=d,
                                     alpha_opt=m_inf_alpha,
                                     beta_opt=m_sup_beta,
                                     entry_polyhedron=polyhedra_2)

            # If max(m_a,m_b) <= (d+m_a)/(1-c) <= M_b -> beta = (d+m_a)/(1-c)

            constraints_3 = [
                middle_rle <= m_sup_beta,
                middle_rle >= m_inf_beta,
                middle_rle >= m_inf_alpha
            ]

            polyhedra_3 = deepcopy(entry_polyhedron)
            polyhedra_3.add_constraints(cstr.cs_from_list(constraints_3))
            sub_spline_3 = f_inf_min(c=c, d=d,
                                     alpha_opt=m_inf_alpha,
                                     beta_opt=middle_rle,
                                     entry_polyhedron=polyhedra_3)

            return plf.Spline(
                [
                    sub_spline_0,
                    sub_spline_1,
                    sub_spline_2,
                    sub_spline_3
                ]
            )

        # Total: perm_value min(beta-alpha, c*beta+d)
    elif not f_is_inf and g_is_inf:  # beta = M_b

        if a <= 0:
            # alpha = m_a
            return g_inf_min(a=a, b=b,
                             alpha_opt=m_inf_alpha,
                             beta_opt=m_sup_beta,
                             entry_polyhedron=entry_polyhedron)

        elif a >= 0:

            # If (M_b - b)/(a+1) <= m_a : alpha = m_a
            middle_rle = (m_sup_beta - b) * Fraction(1, a + 1)

            constraints_0 = [
                middle_rle <= m_inf_alpha
            ]
            polyhedra_0 = deepcopy(entry_polyhedron)
            polyhedra_0.add_constraints(cstr.cs_from_list(constraints_0))
            sub_spline_0 = g_inf_min(a=a, b=b,
                                     alpha_opt=m_inf_alpha,
                                     beta_opt=m_sup_beta,
                                     entry_polyhedron=polyhedra_0)

            # If m_a <= (M_b - b)/(a+1) <= M_a : alpha = (M_b -b)/(a+1)

            constraints_1 = [
                middle_rle <= m_sup_alpha,
                middle_rle >= m_inf_alpha,
                middle_rle <= m_sup_beta
            ]
            polyhedra_1 = deepcopy(entry_polyhedron)
            polyhedra_1.add_constraints(cstr.cs_from_list(constraints_1))
            sub_spline_1 = g_inf_min(a=a, b=b,
                                     alpha_opt=middle_rle,
                                     beta_opt=m_sup_beta,
                                     entry_polyhedron=polyhedra_1)

            # If (M_b-b)/(a+1) >= M_a: alpha = min(M_a, M_b)

            constraints_2 = [
                middle_rle >= m_sup_beta,
                m_sup_alpha >= m_sup_beta
            ]
            polyhedra_2 = deepcopy(entry_polyhedron)
            polyhedra_2.add_constraints(cstr.cs_from_list(constraints_2))
            sub_spline_2 = g_inf_min(a=a, b=b,
                                     alpha_opt=m_sup_alpha,
                                     beta_opt=m_sup_beta,
                                     entry_polyhedron=polyhedra_2)

            constraints_3 = [
                middle_rle >= m_sup_beta,
                m_sup_alpha <= m_sup_beta
            ]
            polyhedra_3 = deepcopy(entry_polyhedron)
            polyhedra_3.add_constraints(cstr.cs_from_list(constraints_3))
            sub_spline_3 = g_inf_min(a=a, b=b,
                                     alpha_opt=m_sup_beta,
                                     beta_opt=m_sup_beta,
                                     entry_polyhedron=polyhedra_3)
            return plf.Spline(
                [
                    sub_spline_0,
                    sub_spline_1,
                    sub_spline_2,
                    sub_spline_3
                ]
            )

        # Total: perm_value = min( M_b - m_a, am_a +b)
    elif not f_is_inf and not g_is_inf:

        if a <= 0 and c >= 0:
            # alpha = m_a
            # beta = M_b
            return f_g_min(a=a, b=b, c=c, d=d,
                           alpha_opt=m_inf_alpha,
                           beta_opt=m_sup_beta,
                           entry_polyhedron=entry_polyhedron)

        elif a >= 0 and c >= 0:
            # Beta = M_b

            # If (M_b - b)/(a+1) <= m_a : alpha = m_a
            middle_rle = (m_sup_beta - b) * Fraction(1, a + 1)

            constraints_0 = [
                middle_rle <= m_inf_alpha
            ]
            polyhedra_0 = deepcopy(entry_polyhedron)
            polyhedra_0.add_constraints(cstr.cs_from_list(constraints_0))
            sub_spline_0 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_inf_alpha,
                                   beta_opt=m_sup_beta,
                                   entry_polyhedron=polyhedra_0)

            # If m_a <= (M_b - b)/(a+1) <= M_a : alpha = (M_b -b)/(a+1)

            constraints_1 = [
                middle_rle <= m_sup_alpha,
                middle_rle >= m_inf_alpha,
                middle_rle <= m_sup_beta
            ]
            polyhedra_1 = deepcopy(entry_polyhedron)
            polyhedra_1.add_constraints(cstr.cs_from_list(constraints_1))
            sub_spline_1 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=middle_rle,
                                   beta_opt=m_sup_beta,
                                   entry_polyhedron=polyhedra_1)

            # If (M_b-b)/(a+1) >= M_a: alpha = min(M_a, M_b)

            constraints_2 = [
                middle_rle >= m_sup_beta,
                m_sup_alpha >= m_sup_beta
            ]
            polyhedra_2 = deepcopy(entry_polyhedron)
            polyhedra_2.add_constraints(cstr.cs_from_list(constraints_2))
            sub_spline_2 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_sup_alpha,
                                   beta_opt=m_sup_beta,
                                   entry_polyhedron=polyhedra_2)

            constraints_3 = [
                middle_rle >= m_sup_beta,
                m_sup_alpha <= m_sup_beta
            ]
            polyhedra_3 = deepcopy(entry_polyhedron)
            polyhedra_3.add_constraints(cstr.cs_from_list(constraints_3))
            sub_spline_3 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_sup_beta,
                                   beta_opt=m_sup_beta,
                                   entry_polyhedron=polyhedra_3)
            return plf.Spline(
                [
                    sub_spline_0,
                    sub_spline_1,
                    sub_spline_2,
                    sub_spline_3
                ]
            )

        elif a <= 0 and c <= 0:
            # alpha = m_a

            middle_rle = (d + m_inf_alpha) * (Fraction(1, 1 - c))

            # If (d+m_a)/(1-c) <= max(m_a,m_b) and m_a <= m_b ->
            # beta = max(m_a,m_b) = m_b

            constraints_0 = [
                middle_rle <= m_inf_alpha,
                middle_rle <= m_inf_beta,
                m_inf_alpha <= m_inf_beta
            ]
            polyhedra_0 = deepcopy(entry_polyhedron)
            polyhedra_0.add_constraints(cstr.cs_from_list(constraints_0))
            sub_spline_0 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_inf_alpha,
                                   beta_opt=m_inf_beta,
                                   entry_polyhedron=polyhedra_0)

            # If (d+m_a)/(1-c) <= max(m_a,m_b) and m_a >= m_b ->
            # beta = max(m_a,m_b) = m_a

            constraints_1 = [
                middle_rle <= m_inf_alpha,
                middle_rle <= m_inf_beta,
                m_inf_alpha >= m_inf_beta
            ]

            polyhedra_1 = deepcopy(entry_polyhedron)
            polyhedra_1.add_constraints(cstr.cs_from_list(constraints_1))

            sub_spline_1 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_inf_alpha,
                                   beta_opt=m_inf_alpha,
                                   entry_polyhedron=polyhedra_1)

            # If (d+m_a)/(1-c) >= M_b -> beta = M_b
            constraints_2 = [middle_rle >= m_sup_beta]

            polyhedra_2 = deepcopy(entry_polyhedron)
            polyhedra_2.add_constraints(cstr.cs_from_list(constraints_2))
            sub_spline_2 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_inf_alpha,
                                   beta_opt=m_sup_beta,
                                   entry_polyhedron=polyhedra_2)

            # If max(m_a,m_b) <= (d+m_a)/(1-c) <= M_b -> beta = (d+m_a)/(1-c)

            constraints_3 = [
                middle_rle <= m_sup_beta,
                middle_rle >= m_inf_beta,
                middle_rle >= m_inf_alpha
            ]

            polyhedra_3 = deepcopy(entry_polyhedron)
            polyhedra_3.add_constraints(cstr.cs_from_list(constraints_3))
            sub_spline_3 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_inf_alpha,
                                   beta_opt=middle_rle,
                                   entry_polyhedron=polyhedra_3)

            return plf.Spline(
                [
                    sub_spline_0,
                    sub_spline_1,
                    sub_spline_2,
                    sub_spline_3
                ]
            )
        elif a >= 0 and c <= 0:

            # Set the partial evaluations

            f_m_inf_alpha = f.partial_evaluation(alpha, m_inf_alpha)
            f_m_inf_beta = f.partial_evaluation(alpha, m_inf_beta)
            f_m_sup_alpha = f.partial_evaluation(alpha, m_sup_alpha)
            f_m_sup_beta = f.partial_evaluation(alpha, m_sup_beta)

            g_m_inf_alpha = g.partial_evaluation(beta, m_inf_alpha)
            g_m_inf_beta = g.partial_evaluation(beta, m_inf_beta)
            g_m_sup_alpha = g.partial_evaluation(beta, m_sup_alpha)
            g_m_sup_beta = g.partial_evaluation(beta, m_sup_beta)

            # If f <= g and f<= h AT (PARTIAL EVALUATION) [min(M_a, M_b), M_b]: alpha = min(M_a, M_b), beta = M_b
            # Case M_a <= M_b
            constraints_0 = [
                m_sup_alpha <= m_sup_beta,
                f_m_sup_alpha <= g_m_sup_beta,
                f_m_sup_alpha <= m_sup_beta - m_sup_alpha
            ]
            polyhedra_0 = deepcopy(entry_polyhedron)
            polyhedra_0.add_constraints(cstr.cs_from_list(constraints_0))
            sub_spline_0 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_sup_alpha,
                                   beta_opt=m_sup_beta,
                                   entry_polyhedron=polyhedra_0)

            # Case M_a <= M_b
            constraints_1 = [
                m_sup_alpha >= m_sup_beta,
                f_m_sup_beta <= g_m_sup_beta,
                f_m_sup_beta <= 0
            ]
            polyhedra_1 = deepcopy(entry_polyhedron)
            polyhedra_1.add_constraints(cstr.cs_from_list(constraints_1))
            sub_spline_1 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_sup_beta,
                                   beta_opt=m_sup_beta,
                                   entry_polyhedron=polyhedra_1)
            # If g <= h, g<= f at [m_a, max(m_a,m_b)]: alpha=m_a, beta = max(m_a, m_b)
            # Case m_a <= m_b
            # m_a <= m_b, f(m_a, max(m_a, m_b)) <= g(m_a, max(m_a, m_b)) and h(m_a, m_b)
            constraints_2 = [
                m_inf_alpha <= m_inf_beta,
                f_m_inf_alpha <= g_m_inf_beta,
                g_m_inf_beta <= m_inf_beta - m_inf_alpha
            ]
            polyhedra_2 = deepcopy(entry_polyhedron)
            polyhedra_2.add_constraints(cstr.cs_from_list(constraints_2))
            sub_spline_2 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_inf_alpha,
                                   beta_opt=m_inf_beta,
                                   entry_polyhedron=polyhedra_2)

            # Case m_a >= m_b
            constraints_3 = [
                m_inf_alpha >= m_inf_beta,
                f_m_inf_alpha <= g_m_inf_alpha,
                g_m_inf_beta <= 0
            ]
            polyhedra_3 = deepcopy(entry_polyhedron)
            polyhedra_3.add_constraints(cstr.cs_from_list(constraints_3))
            sub_spline_3 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_inf_alpha,
                                   beta_opt=m_inf_alpha,
                                   entry_polyhedron=polyhedra_3)

            # If h <= f , g at (m_a, M_b): alpha = m_a, beta = M_b
            constraints_4 = [
                m_sup_beta - m_inf_alpha <= f_m_inf_alpha,
                m_sup_beta - m_inf_alpha <= g_m_sup_beta
            ]
            polyhedra_4 = deepcopy(entry_polyhedron)
            polyhedra_4.add_constraints(cstr.cs_from_list(constraints_4))
            sub_spline_4 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_inf_alpha,
                                   beta_opt=m_sup_beta,
                                   entry_polyhedron=polyhedra_4)

            # T_a = [d-b(1-c)]/[(a+1)(1-c)-1], T_b [d(a+1)-b]/[(a+1)(1-c) -1]

            intersection_f_g_h = [(d - b * (1 - c)) * Fraction(1, (a + 1) * (1 - c) - 1),
                                  (d * (a + 1) - b) * Fraction(1, (a + 1) * (1 - c) - 1)]

            # If T_b >= M_b: alpha = (M_b - b)/(a+1), beta = M_b
            constraints_5 = [
                intersection_f_g_h[1] >= m_sup_beta
            ]
            polyhedra_5 = deepcopy(entry_polyhedron)
            polyhedra_5.add_constraints(cstr.cs_from_list(constraints_5))
            sub_spline_5 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=(m_sup_beta - b) * Fraction(1, a + 1),
                                   beta_opt=m_sup_beta,
                                   entry_polyhedron=polyhedra_5)

            # If T_a <= m_a: alpha = m_a, beta = (m_a+d)/(1-c)
            constraints_6 = [
                intersection_f_g_h[0] <= m_inf_alpha
            ]
            polyhedra_6 = deepcopy(entry_polyhedron)
            polyhedra_6.add_constraints(cstr.cs_from_list(constraints_6))
            sub_spline_6 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=m_inf_alpha,
                                   beta_opt=(m_inf_alpha + d) * Fraction(1, 1 - c),
                                   entry_polyhedron=polyhedra_6)

            # If ad <= bc, g<= f,h at [min(m_b, M_a), m_b]: alpha = (cm_b + d -b)/a, beta = m_b
            # Case M_a <= m_b NB: the case M_a >= m_b is useless
            constraints_7 = [
                m_sup_alpha <= m_inf_beta,
                a * d - b * c <= 0,
                g_m_inf_beta <= f_m_sup_alpha,
                g_m_inf_beta <= m_inf_beta - m_sup_alpha
            ]
            polyhedra_7 = deepcopy(entry_polyhedron)
            polyhedra_7.add_constraints(cstr.cs_from_list(constraints_7))
            sub_spline_7 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=(c * m_inf_beta + d - b) * Fraction(1, a),
                                   beta_opt=m_inf_beta,
                                   entry_polyhedron=polyhedra_7)

            # If ad <= bc :
            # # g<= f,h at [M_a, max(m_b, M_a)]: alpha = (d-b)/(a-c), beta =(d-b)/(a-c)
            constraints_8 = [
                m_sup_alpha <= m_inf_beta,
                a * d - b * c <= 0,
                g_m_inf_beta <= f_m_sup_alpha,
                g_m_inf_beta <= m_inf_beta - m_sup_alpha
            ]
            polyhedra_8 = deepcopy(entry_polyhedron)
            polyhedra_8.add_constraints(cstr.cs_from_list(constraints_8))
            sub_spline_8 = f_g_min(a=a, b=b, c=c, d=d,
                                   alpha_opt=(c * m_inf_beta + d - b) * Fraction(1, a),
                                   beta_opt=m_inf_beta,
                                   entry_polyhedron=polyhedra_8)

            # If ad <= bc: else
            # TODO: WARNING: check formats and understand the strange equations...
            # If ad >= bc :
            # # If T_b <= m_b: alpha = (1-c)m_b - d, beta = m_b

            constraints_1_0 = [
                a * d - b * c >= 0,
                intersection_f_g_h[1] <= m_sup_beta
            ]
            polyhedra_1_0 = deepcopy(entry_polyhedron)
            polyhedra_1_0.add_constraints(cstr.cs_from_list(constraints_1_0))
            sub_spline_1_0 = f_g_min(a=a, b=b, c=c, d=d,
                                     alpha_opt=(1 - c) * m_inf_beta - d,
                                     beta_opt=m_inf_beta,
                                     entry_polyhedron=polyhedra_1_0)

            # # If T_a >= M_a: alpha = M_a, beta = (a+1)M_a + b
            constraints_1_1 = [
                a * d - b * c >= 0,
                intersection_f_g_h[0] >= m_sup_alpha
            ]
            polyhedra_1_1 = deepcopy(entry_polyhedron)
            polyhedra_1_1.add_constraints(cstr.cs_from_list(constraints_1_1))
            sub_spline_1_1 = f_g_min(a=a, b=b, c=c, d=d,
                                     alpha_opt=m_sup_alpha,
                                     beta_opt=(a + 1) * m_sup_alpha + b,
                                     entry_polyhedron=polyhedra_1_1)

            # # Else: alpha = T_a, beta = T_b
            constraints_1_2 = [
                a * d - b * c >= 0,
                intersection_f_g_h[1] >= m_sup_beta,
                intersection_f_g_h[0] >= m_sup_alpha
            ]
            polyhedra_1_2 = deepcopy(entry_polyhedron)
            polyhedra_1_2.add_constraints(cstr.cs_from_list(constraints_1_2))
            sub_spline_1_2 = f_g_min(a=a, b=b, c=c, d=d,
                                     alpha_opt=intersection_f_g_h[0],
                                     beta_opt=intersection_f_g_h[1],
                                     entry_polyhedron=polyhedra_1_2)

            return plf.Spline(
                [
                    sub_spline_0,
                    sub_spline_1,
                    sub_spline_2,
                    sub_spline_3,
                    sub_spline_4,
                    sub_spline_5,
                    sub_spline_6,
                    sub_spline_7,
                    sub_spline_8,
                    # sub_spline_9,
                    sub_spline_1_0,
                    sub_spline_1_1,
                    sub_spline_1_2
                ]
            )

        # If ad <= bc :
        # # g<= f,h at [min(m_b, M_a), m_b]: alpha = (cm_b + d -b)/a, beta = m_b
        # # g<= f,h at [M_a, max(m_b, M_a)]: alpha = (d-b)/(a-c), beta =(d-b)/(a-c)
        # # else: alpha = M_a, beta = (aM_a + b - d)/c
        # If ad <= bc :
        # # If T_b <= m_b: alpha = (1-c)m_b - d, beta = m_b
        # # If T_a >= M_a: alpha = M_a, beta = (a+1)M_a + b
        # # Else: alpha = T_a, beta = T_b
