# coding=utf-8
from __future__ import annotations  # For forward reference typing

from copy import deepcopy
from typing import Tuple, List, NamedTuple

import ppl
import pysymbrobustness.ppl_extension.linear_algebra as lin_alg
import pysymbrobustness.ppl_extension.rational_linear_algebra as rla
import pysymbrobustness.ppl_extension.piecewise_linear_function as plf
from fractions import Fraction


class CaseData(NamedTuple):
    c_list: List[ppl.Constraint]
    alpha_opt: plf.LinearFunction
    beta_opt: plf.LinearFunction


def get_coefficient(fct: rla.RationalLinearExpression) -> Tuple[
    Fraction, rla.RationalLinearExpression]:
    """
    From a function fct = b, returns the tuple a, b where a is the sum of the coefficients of b.
    @param fct:
    @type fct:
    @return:
    @rtype:
    """
    if type(fct) != rla.RationalLinearExpression:
        raise TypeError

    b = fct
    a = sum(fct.coefficients())

    return a, b


#
# def get_coefficient(fct: rla.RationalLinearExpression, v: ppl.Variable) -> \
#         Tuple[Fraction, rla.RationalLinearExpression]:
#     """
#     @param fct: a function of the form v_coeff*v + b, where v_coeff is a fraction, v a ppl.Variable and b a
#     rla.RationalLinearExpression that does not depend on v.
#     @type fct:
#     @param v:
#     @type v:
#     @return:
#     @rtype:
#     """
#
#     if type(fct) != rla.RationalLinearExpression:
#         raise TypeError
#     index = v.id()
#     v_coefficient = fct.coefficient(v)
#     coeffs = [0 for i in range(v.id()+1)]
#     coeffs[index] = fct.coefficient(v)
#     v_term = rla.RationalLinearExpression(tuple(coeffs))
#     b = fct - v_term
#
#     return v_coefficient, b


def compute_permissiveness(poly: ppl.C_Polyhedron, f: plf.LinearFunction,
                           g: plf.LinearFunction,
                           alpha_opt: plf.LinearFunction,
                           beta_opt: plf.LinearFunction):
    if (type(f) != rla.InfiniteExpression and type(
            f) != rla.RationalLinearExpression) or \
            (type(g) != rla.InfiniteExpression and type(
                g) != rla.RationalLinearExpression):
        raise TypeError

    f_linear = type(f) == rla.RationalLinearExpression
    g_linear = type(g) == rla.RationalLinearExpression

    h_opt = plf.Spline(plf.SubSpline(polyhedron=poly,
                                     function=beta_opt - alpha_opt))

    if f_linear and g_linear:
        a, b = get_coefficient(f)
        c, d = get_coefficient(g)
        f_opt = plf.Spline(
            plf.SubSpline(polyhedron=poly, function=a * alpha_opt + b))
        g_opt = plf.Spline(
            plf.SubSpline(polyhedron=poly, function=c * beta_opt + d))
        return f_opt.minimum_list([g_opt, h_opt])

    elif not f_linear and not g_linear:
        return h_opt
    else:
        a, b = get_coefficient(
            f if f_linear else g)
        # alpha if f_linear else beta)
        opt = alpha_opt if f_linear else beta_opt
        spline = plf.Spline(
            plf.SubSpline(polyhedron=poly, function=a * opt + b))
        return h_opt.minimum(spline)


def opti_f_inf_g_inf(f: plf.LinearFunction, g: plf.LinearFunction,
                     entry_poly: ppl.C_Polyhedron,
                     entry_interval: Tuple[
                         plf.LinearFunction, ...]) -> plf.Spline:
    min_alpha, max_alpha, min_beta, max_beta = entry_interval

    return compute_permissiveness(poly=entry_poly, f=f, g=g,
                                  alpha_opt=min_alpha, beta_opt=max_beta)


def name_to_find(cs_list: List[ppl.Constraint], poly: ppl.C_Polyhedron,
                 f: plf.LinearFunction, g: plf.LinearFunction,
                 alpha_opt: plf.LinearFunction,
                 beta_opt: plf.LinearFunction) -> plf.Spline:
    cs = lin_alg.cs_from_list(cs_list)
    p_0 = deepcopy(poly)
    p_0.add_constraints(cs)
    return compute_permissiveness(poly=p_0, f=f, g=g, alpha_opt=alpha_opt,
                                  beta_opt=beta_opt)


def a_c_same_sign(f: plf.LinearFunction, g: plf.LinearFunction,
                  entry_poly: ppl.C_Polyhedron,
                  entry_interval: Tuple[plf.LinearFunction, ...],
                  negative: bool) -> plf.Spline:
    min_alpha, max_alpha, min_beta, max_beta = entry_interval
    ac, bd = get_coefficient(g) if negative else get_coefficient(f)

    middle_rle = (min_alpha + bd) if negative else (max_beta - bd)
    middle_factor = Fraction(1, 1 - ac) if negative else Fraction(1, 1 + ac)
    rle_in_the_middle = middle_rle * middle_factor

    permissiveness = plf.Spline()
    min_f, max_f = (min_beta, max_alpha) if negative else (min_alpha, max_beta)

    constraint_list = [
        ([max_f <= rle_in_the_middle], max_f),
        ([min_f <= rle_in_the_middle, rle_in_the_middle <= max_f], rle_in_the_middle),  #NB: Generates the 1/2 function for the (0,0), (1,1), (1,0) triangle
        ([rle_in_the_middle <= min_f], min_f)
    ]

    for constraints, opt in constraint_list:
        alpha_opt = min_alpha if negative else opt
        beta_opt = opt if negative else max_beta

        spline = name_to_find(constraints, poly=entry_poly, f=f, g=g,
                              alpha_opt=alpha_opt, beta_opt=beta_opt)
        permissiveness.fusion(other=spline)

    return permissiveness


def opti_f_inf_g_rle(f: plf.LinearFunction, g: plf.LinearFunction,
                     entry_poly: ppl.C_Polyhedron,
                     entry_interval: Tuple[
                         plf.LinearFunction, ...]) -> plf.Spline:
    min_alpha, max_alpha, min_beta, max_beta = entry_interval
    c, d = get_coefficient(g)
    if c >= 0:
        return compute_permissiveness(poly=entry_poly, f=f, g=g,
                                      alpha_opt=min_alpha, beta_opt=max_beta)
    elif c <= 0:
        # Same case a f and g rational linear expression, and c <= 0 and a <= 0
        return a_c_same_sign(f, g, entry_poly, entry_interval, negative=True)


def opti_f_rle_g_inf(f: plf.LinearFunction, g: plf.LinearFunction,
                     entry_poly: ppl.C_Polyhedron,
                     # alpha: ppl.Variable, beta: ppl.Variable,
                     entry_interval: Tuple[
                         plf.LinearFunction, ...]) -> plf.Spline:
    a, b = get_coefficient(f)
    min_alpha, max_alpha, min_beta, max_beta = entry_interval
    if a <= 0:
        return compute_permissiveness(poly=entry_poly, f=f, g=g,
                                      alpha_opt=min_alpha, beta_opt=max_beta)
    elif a >= 0:
        # Same case a f and g rational linear expression, and c >= 0 and a >= 0
        return a_c_same_sign(f, g, entry_poly, entry_interval, negative=False)


def a_pos_c_neg(f: plf.LinearFunction, g: plf.LinearFunction,
                entry_poly: ppl.C_Polyhedron,
                entry_interval: Tuple[plf.LinearFunction, ...]):
    a, b = get_coefficient(f)
    c, d = get_coefficient(g)
    min_alpha, max_alpha, min_beta, max_beta = entry_interval

    # Representation of the entry set trapeze A-B-E-C-D

    A = (min_alpha, max_beta)
    B = (max_alpha, max_beta)
    C = (min_alpha, min_beta)

    def d_pt(alpha_sup_beta: bool):
        return (max_alpha, min_beta) if alpha_sup_beta else (min_beta,
                                                             min_beta)

    def e_pt(alpha_sup_beta: bool):
        return (max_alpha, max_alpha) if alpha_sup_beta else (max_alpha,
                                                              min_beta)

    D_0 = d_pt(False)
    D_1 = d_pt(True)
    E_0 = e_pt(False)
    E_1 = e_pt(True)

    T_den: Fraction = Fraction(1, (a + 1) * (1 - c) - 1)
    T_alpha, T_beta = (d - b * (1 - c)) * T_den, (d * (a + 1) - b) * T_den

    def _fct_eval(is_f: bool, is_g: bool, point: Tuple[plf.LinearFunction,
                                                       plf.LinearFunction]):
        """
        evaluation f(X) (if is_f) or g(Y) (if is_g) or h(X,Y) = Y-X (
        otherwise) in point=(X,Y)
        """
        if is_f:
            return a * point[0] + b  # f(X) = aX+b
        elif is_g:
            return c * point[1] + d  # g(Y) = cY+D
        else:
            return point[1] - point[0]  # h(X,Y) = Y-X

    def g_ev(P: Tuple[plf.LinearFunction, plf.LinearFunction]):
        return _fct_eval(False, True, P)

    def f_ev(P: Tuple[plf.LinearFunction, plf.LinearFunction]):
        return _fct_eval(True, False, P)

    def h_ev(P: Tuple[plf.LinearFunction, plf.LinearFunction]):
        return _fct_eval(False, False, P)

    case_0 = CaseData(
        c_list=[
            f_ev(B) <= g_ev(B),
            f_ev(B) <= h_ev(B)  # f minimal on B
        ],
        alpha_opt=max_alpha,
        beta_opt=max_beta
    )

    case_1 = CaseData(
        c_list=[
            g_ev(C) <= f_ev(C),
            g_ev(C) <= h_ev(C)  # g minimal on C
        ],
        alpha_opt=min_alpha,
        beta_opt=min_beta
    )

    case_2 = CaseData(
        c_list=[
            h_ev(A) <= f_ev(A),
            h_ev(A) <= g_ev(A)  # h minimal on A
        ],
        alpha_opt=min_alpha,
        beta_opt=max_beta
    )

    case_3 = CaseData(
        c_list=[T_beta >= max_beta,  # no g zone
                f_ev(A) <= h_ev(A),
                f_ev(A) <= g_ev(A),  # f minimal on A
                h_ev(B) <= f_ev(B),
                h_ev(B) <= g_ev(B)  # h minimal on B
                ],
        alpha_opt=(max_beta - b) * Fraction(1, a + 1),
        beta_opt=max_beta
    )

    case_4 = CaseData(
        c_list=[T_alpha <= min_alpha,  # no f zone
                g_ev(A) <= h_ev(A),
                g_ev(A) <= f_ev(A),  # g minimal on A
                h_ev(C) <= g_ev(C),
                h_ev(C) <= f_ev(C)  # h minimal on C
                ],
        alpha_opt=min_alpha,
        beta_opt=(min_alpha + d) * Fraction(1, 1 - c)
    )

    # D = D_0 <=> alpha_sup_beta = False
    case_5_0 = CaseData(
        c_list=[
                f_ev(C) <= g_ev(C),
                f_ev(C) <= h_ev(C),  # f minimal on C
                min_beta >= max_alpha,  # alpha_sup_beta = False
                g_ev(D_0) <= f_ev(D_0),
                g_ev(D_0) <= h_ev(D_0)  # g minimal on D_0
                ],
        alpha_opt=(1 - c) * min_beta - d,
        beta_opt= min_beta
    )

    # D = D_1

    case_5_1 = CaseData(
        c_list=[
                f_ev(C) <= h_ev(C),
                f_ev(C) <= g_ev(C),  # f minimal on C
                min_beta <= max_alpha,  # alpha_sup_beta = True
                g_ev(D_1) <= f_ev(D_1),
                g_ev(D_1) <= h_ev(D_1)  # g minimal on D_1
                ],
        alpha_opt=(1 - c) * min_beta - d,
        beta_opt=min_beta
    )

    # E = E_0 <=> alpha_sup_beta = False
    case_6_0 = CaseData(
        c_list=[g_ev(B) <= h_ev(B),
                g_ev(B) <= f_ev(B),  # g minimal on B
                min_beta >= max_alpha,  # alpha_sup_beta = False
                f_ev(E_0) <= g_ev(E_0),
                f_ev(E_0) <= h_ev(E_0)  # f minimal on E_0
                ],
        alpha_opt=max_alpha,
        beta_opt=min_beta
    )

    # E = E_1
    case_6_1 = CaseData(
        c_list=[g_ev(B) <= h_ev(B),
                g_ev(B) <= f_ev(B),  # g minimal on B
                min_beta <= max_alpha,  # alpha_sup_beta = True
                f_ev(E_1) <= g_ev(E_1),
                f_ev(E_1) <= h_ev(E_1)  # f minimal on E_1
                ],
        alpha_opt=max_alpha,
        beta_opt=max_alpha
    )

    # E = E_0, D = D_0 <=> alpha_sup_beta = False
    case_7_0 = CaseData(
        c_list=[f_ev(D_0) <= h_ev(D_0),
                f_ev(D_0) <= g_ev(D_0),  # f minimal on D
                min_beta >= max_alpha,  # alpha_sup_beta = False
                g_ev(E_0) <= f_ev(E_0),
                g_ev(E_0) <= h_ev(E_0)  # g minimal on E_0
                ],
        alpha_opt=(b - d) * Fraction(1, c - a),
        beta_opt=(b - d) * Fraction(1, c - a)
    )

    # E = E_1; D = D_1
    case_7_1 = CaseData(
        c_list=[f_ev(D_1) <= h_ev(D_1),
                f_ev(D_1) <= g_ev(D_1),  # f minimal on D
                min_beta <= max_alpha,  # alpha_sup_beta = True
                g_ev(E_1) <= f_ev(E_1),
                g_ev(E_1) <= h_ev(E_1)  # g minimal on E_1
                ],
        alpha_opt=(b - d) * Fraction(1, c - a),
        beta_opt=(b - d) * Fraction(1, c - a)
    )

    # E = E_0 <=> alpha_sup_beta = False
    case_8_0 = CaseData(
        c_list=[T_alpha >= max_alpha,
                g_ev(B) <= h_ev(B),
                g_ev(B) <= f_ev(B),  # g minimal on B
                min_beta >= max_alpha,  # alpha_sup_beta = False
                h_ev(E_0) <= f_ev(E_0),
                h_ev(E_0) <= g_ev(E_0)  # h minimal on E_0
                ],
        alpha_opt=max_alpha,
        beta_opt=(a+1) * max_alpha + b
    )

    # E = E_1
    case_8_1 = CaseData(
        c_list=[T_alpha >= max_alpha,
                g_ev(B) <= h_ev(B),
                g_ev(B) <= f_ev(B),  # g minimal on B
                min_beta <= max_alpha,  # alpha_sup_beta = True
                h_ev(E_1) <= f_ev(E_1),
                h_ev(E_1) <= g_ev(E_1)  # h minimal on E_1
                ],
        alpha_opt=max_alpha,
        beta_opt=(a+1) * max_alpha + b
    )

    # D = D_0 <=> alpha_sup_beta = False
    case_9_0 = CaseData(
        c_list=[T_beta <= min_beta,
                f_ev(C) <= g_ev(C),
                f_ev(C) <= h_ev(C),  # f minimal on C
                min_beta >= max_alpha,  # alpha_sup_beta = False
                h_ev(D_0) <= f_ev(D_0),
                h_ev(D_0) <= g_ev(D_0)  # h minimal on D_0
                ],
        alpha_opt=c * min_beta * Fraction(1, a) + ( d - b) * Fraction(1, a),
        beta_opt=min_beta
    )

    # D = D_1

    case_9_1 = CaseData(
        c_list=[T_beta <= min_beta,
                f_ev(C) <= g_ev(C),
                f_ev(C) <= h_ev(C),  # g minimal on C
                min_beta <= max_alpha,  # alpha_sup_beta = True
                h_ev(D_1) <= f_ev(D_1),
                h_ev(D_1) <= g_ev(D_1)  # g minimal on D_1
                ],
        alpha_opt=c * min_beta * Fraction(1, a) + (d - b) * Fraction(1, a),
        beta_opt=min_beta
    )


    """
    case_5 = CaseData(
        c_list=[T_beta <= T_alpha,
                min_beta >= max_alpha,  # alpha_sup_beta = False
                g_ev(D_0) <= f_ev(D_0), g_ev(D_0) <= h_ev(D_0),  # g minimal on D
                f_ev(C) <= g_ev(C), f_ev(C) <= h_ev(C)  # f minimal on C
                ],
        alpha_opt=(c * min_beta + d - b) * Fraction(1, a),
        beta_opt=min_beta
    )

    case_6 = CaseData(
        c_list=[T_beta <= T_alpha,
                min_beta <= max_alpha,  # alpha_sup_beta = True
                g_ev(D_1) <= f_ev(D_1), g_ev(D_1) <= h_ev(D_1),  # g minimal on D
                f_ev(C) <= g_ev(C), f_ev(C) <= h_ev(C)  # f minimal on C
                ],
        alpha_opt=(c * min_beta + d - b) * Fraction(1, a),
        beta_opt=min_beta
    )

    case_7 = CaseData(
        c_list=[T_beta <= T_alpha,
                g_ev(B) <= f_ev(B), g_ev(B) <= h_ev(B),  # g minimal on B
                min_beta <= max_alpha,  # alpha_sup_beta = True
                f_ev(E_1) <= g_ev(E_1), f_ev(E_1) <= h_ev(E_1)],  # f minimal on E
        alpha_opt=max_alpha,
        beta_opt=(a * max_alpha + b - d) * Fraction(1, c)
    )

    case_8 = CaseData(
        c_list=[T_beta <= T_alpha,
                g_ev(B) <= f_ev(B), g_ev(B) <= h_ev(B),  # g minimal on B
                min_beta >= max_alpha,  # alpha_sup_beta = False
                f_ev(E_0) <= g_ev(E_0), f_ev(E_0) <= h_ev(E_0)],  # f minimal on E
        alpha_opt=max_alpha,
        beta_opt=(a * max_alpha + b - d) * Fraction(1, c)
    )

    case_9 = CaseData(
        c_list=[T_alpha <= T_beta,
                T_beta <= min_beta,
                min_beta >= max_alpha,
                g_ev(C) <= f_ev(C), g_ev(C) <= h_ev(C),  # g minimal on C
                h_ev(D_0) <= g_ev(D_0), h_ev(D_0) <= f_ev(D_0)
                ],
        alpha_opt=(1 - c) * min_beta - d,
        beta_opt=min_beta
    )

    case_1_0 = CaseData(
        c_list=[T_alpha <= T_beta,
                T_beta <= min_beta,
                min_beta <= max_alpha,
                g_ev(C) <= f_ev(C), g_ev(C) <= h_ev(C),  # g minimal on C
                h_ev(D_1) <= g_ev(D_1), h_ev(D_1) <= f_ev(D_1)
                ],
        alpha_opt=(1 - c) * min_beta - d,
        beta_opt=min_beta
    )

    case_1_1 = CaseData(
        c_list=[
            T_alpha <= T_beta,
            T_alpha >= max_alpha,
            min_beta >= max_alpha,
            f_ev(B) <= g_ev(B), f_ev(B) <= h_ev(B),  # f minimal on B
            h_ev(E_0) <= f_ev(E_0), h_ev(E_0) <= g_ev(E_0)  # h minimal on E

        ],
        alpha_opt=max_alpha,
        beta_opt=(a + 1) * max_alpha + b
    )

    case_1_2 = CaseData(
        c_list=[
            T_alpha <= T_beta,
            T_alpha >= max_alpha,
            min_beta <= max_alpha,  # alpha_sup_beta = True
            f_ev(B) <= g_ev(B), f_ev(B) <= h_ev(B),  # f minimal on B
            h_ev(E_1) <= f_ev(E_1), h_ev(E_1) <= g_ev(E_1)  # h minimal on E

        ],
        alpha_opt=max_alpha,
        beta_opt=(a + 1) * max_alpha + b
    )
    """
    case_1_0 = CaseData(
        c_list=[
            max_beta >= T_beta, T_beta >= min_beta,
            max_alpha >= T_alpha, T_alpha >= min_alpha,
            T_beta >= T_alpha
        ],
        alpha_opt=T_alpha,
        beta_opt=T_beta
    )
    cases = [case_0, case_1, case_2, case_3, case_4, case_5_0, case_5_1, case_6_0, case_6_1, case_7_0, case_7_1,
             case_8_0, case_8_1, case_9_0, case_9_1, case_1_0]

    permissiveness = plf.Spline()
    for constraints, alpha_opt, beta_opt in cases:
        spline = name_to_find(constraints, poly=entry_poly, f=f, g=g,
                              alpha_opt=alpha_opt, beta_opt=beta_opt)
        permissiveness.fusion(other=spline)

    return permissiveness


def opti_f_rle_g_rle(f: plf.LinearFunction, g: plf.LinearFunction,
                     entry_poly: ppl.C_Polyhedron,
                     entry_interval: Tuple[
                         plf.LinearFunction, ...]) -> plf.Spline:
    a, b = get_coefficient(f)
    c, d = get_coefficient(g)
    min_alpha, max_alpha, min_beta, max_beta = entry_interval

    if a <= 0 and c >= 0:
        return compute_permissiveness(poly=entry_poly,
                                      f=f, g=g,
                                      alpha_opt=min_alpha, beta_opt=max_beta)

    elif (a >= 0 and c >= 0) or (a <= 0 and c <= 0):
        return a_c_same_sign(f, g, entry_poly, entry_interval, negative=a <= 0)
    elif a > 0 and c < 0:
        return a_pos_c_neg(f, g, entry_poly, entry_interval)


def plf_optimization(
        f: plf.LinearFunction, g: plf.LinearFunction,
        entry_poly: ppl.C_Polyhedron,
        entry_set: Tuple[
            plf.SubSpline, plf.SubSpline, plf.SubSpline, plf.SubSpline]) -> plf.Spline:
    """
    @param f:
    @type f:
    @param g:
    @type g:
    @param entry_poly:
    @type entry_poly:
    @param entry_set:
    @type entry_set:
    @return:
    @rtype:
    """
    entry_interval: Tuple[plf.LinearFunction, ...] = tuple(
        map(lambda x: x.function, entry_set))

    copy_entry_poly = deepcopy(entry_poly)
    for sub in entry_set:
        copy_entry_poly.intersection_assign(sub.polyhedron)

    min_alpha, max_alpha, min_beta, max_beta = entry_interval

    copy_entry_poly.add_constraints(lin_alg.cs_from_list([
        min_alpha <= min_beta, max_alpha <= max_beta]))
    # entry_poly = entry_poly.intersection_assign()

    if type(f) == rla.InfiniteExpression and type(g) == rla.InfiniteExpression:
        return opti_f_inf_g_inf(f, g, copy_entry_poly, entry_interval)
    elif type(f) == rla.InfiniteExpression and type(
            g) == rla.RationalLinearExpression:
        return opti_f_inf_g_rle(f, g, copy_entry_poly, entry_interval)
    elif type(f) == rla.RationalLinearExpression and type(
            g) == rla.InfiniteExpression:
        return opti_f_rle_g_inf(f, g, copy_entry_poly, entry_interval)
    elif type(f) == rla.RationalLinearExpression and type(
            g) == rla.RationalLinearExpression:
        return opti_f_rle_g_rle(f, g, copy_entry_poly, entry_interval)
    else:
        raise TypeError
