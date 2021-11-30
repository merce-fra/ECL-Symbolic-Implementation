import ppl as ppl
import pysymbrobustness.ppl_extension.linear_algebra as la
import pysymbrobustness.ppl_extension.rational_linear_algebra as rla
import pysymbrobustness.ppl_extension.piecewise_linear_function as plf
from fractions import Fraction


def main():
    x = ppl.Variable(0)
    y = ppl.Variable(1)
    z = ppl.Variable(2)

    # c = la.cs_from_list([x <= 1, y <= 1])
    # f = rla.RationalLinearExpression((0, Fraction(-1, 4)), Fraction(1, 4))
    # g = rla.RationalLinearExpression((Fraction(-1, 4), 0), Fraction(1, 4))
    # P = ppl.C_Polyhedron(c)
    #
    # sub_spline_0 = plf.SubSpline(polyhedron=P, function=f)
    # sub_spline_1 = plf.SubSpline(polyhedron=P, function=g)
    #
    # min_spline = sub_spline_0.minimum(sub_spline_1)
    # max_spline = sub_spline_0.maximum(sub_spline_1)
    # print("minimum=")
    # print(min_spline.sub_splines[0].function)
    # print(min_spline.sub_splines[0].polyhedron.constraints())
    # print(min_spline.sub_splines[1].function)
    # print(min_spline.sub_splines[1].polyhedron.constraints())
    #
    # print("maximum=")
    # print(max_spline.sub_splines[0].function)
    # print(max_spline.sub_splines[0].polyhedron.constraints())
    # print(max_spline.sub_splines[1].function)
    # print(max_spline.sub_splines[1].polyhedron.constraints())
    #
    # print("test of maximum between splines")
    #
    # new_spline = min_spline.min_max_op(max_spline, True)
    # print(new_spline)
    # print(new_spline.sub_splines[0].function)
    # print(new_spline.sub_splines[0].polyhedron.constraints())
    # print(new_spline.sub_splines[1].function)
    # print(new_spline.sub_splines[1].polyhedron.constraints())

    es_f = ppl.C_Polyhedron(la.cs_from_list([x >= 0, y >= 0, x <= 1, y <= 1, 1 + x <= 2 * y]))
    es_g = ppl.C_Polyhedron(la.cs_from_list([x >= 0, y >= 0, x <= 1, y <= 1, 1 + x >= 2 * y]))
    f = rla.RationalLinearExpression((0, -1), 1)
    g = rla.RationalLinearExpression((Fraction(-1, 2), 0), Fraction(1, 2))

    es_h = ppl.C_Polyhedron(la.cs_from_list([x >= 0, y >= 0, x <= 1, y <= 1, x >= y]))
    es_i = ppl.C_Polyhedron(la.cs_from_list([x >= 1, x <= 2, y >= 0, y <= 1, -1 + x >= y]))
    es_j = ppl.C_Polyhedron(la.cs_from_list([x >= 0, y >= 0, x >= 1, x <= 2, y >= 0, y <= 1, -1 + x <= y]))
    h = rla.RationalLinearExpression((1, -1), 0)
    i = rla.RationalLinearExpression((-1, 0), 2)
    j = rla.RationalLinearExpression((0, -1), 1)

    ex_1_sub_spline_0 = plf.Spline([plf.SubSpline(polyhedron=es_f, function=f),
                                    plf.SubSpline(polyhedron=es_g, function=g)])
    ex_1_sub_spline_1 = plf.Spline([plf.SubSpline(polyhedron=es_h, function=h),
                                    plf.SubSpline(polyhedron=es_i, function=i),
                                    plf.SubSpline(polyhedron=es_j, function=j)])

    final = ex_1_sub_spline_0.operator(ex_1_sub_spline_1, False)
    print(final.functions)
    print(final)
    for i in range(7):
        print(final.sub_splines[i].function)
        print(final.sub_splines[i].polyhedron.lin_alg())


if __name__ == '__main__':
    main()
