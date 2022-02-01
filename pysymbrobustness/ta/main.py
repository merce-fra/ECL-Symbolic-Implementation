# coding=utf-8
# import math
import dataclasses
import random
import time
from fractions import Fraction
from pathlib import Path
from typing import List

import pysymbrobustness.ppl_extension.rational_linear_algebra as rla
import pysymbrobustness.ppl_extension.linear_algebra as lin_alg
import pysymbrobustness.ppl_extension.piecewise_linear_function as plf
import pysymbrobustness.ppl_extension.polyhedra as polyhedra
import ppl as ppl
import pysymbrobustness.ta.timedauto as timedauto
import pysymbrobustness.ta.guards as guards
import pysymbrobustness.ta.creators as creators
import pysymbrobustness.permissiveness_computation.explorer as explorer
from pysymbrobustness.ta.misc_debug import visual, latex_output

Variable = ppl.Variable
Constraint_System = ppl.Constraint_System
C_Polyhedron = ppl.C_Polyhedron
Constraint = ppl.Constraint

x = ppl.Variable(0)
y = ppl.Variable(1)


def polyhedron_example(constraint_list) -> C_Polyhedron:
    cs = Constraint_System()
    for constr in constraint_list:
        cs.insert(constr)
    return C_Polyhedron(cs)


def formats_polyhedrons_0():
    """ formats partition examples for the computation of permissiveness for
    the timed automata with two clocks x,y
    and two linear transition 0 <= x <= 1, 0 <= y <= 1"""
    c_0 = [y <= x, x <= 1, y <= 1, x >= 0, y >= 0]
    c_1 = [y + 1 >= x, x <= 2, y <= 1, x >= 1, y >= 0]
    c_2 = [y + 1 <= x, x <= 2, y <= 1, x >= 1, y >= 0]

    return ppl.C_Polyhedron(lin_alg.cs_from_list(c_0)), \
           ppl.C_Polyhedron(lin_alg.cs_from_list(c_1)), \
           ppl.C_Polyhedron(lin_alg.cs_from_list(c_2))


def formats_functions_0():
    """ formats functions examples for the computation of permissiveness for
    the timed automata with two clocks x,y and two linear transition 0 <= x
    <= 1, 0 <= y <= 1"""
    f_0 = rla.RationalLinearExpression((1, -1), 0)
    f_1 = rla.RationalLinearExpression((0, -1), 1)
    f_2 = rla.RationalLinearExpression((-1, 0), 2)
    return f_0, f_1, f_2


def spline_0():
    c_0, c_1, c_2 = formats_polyhedrons_0()
    f_0, f_1, f_2 = formats_functions_0()
    return plf.Spline([plf.SubSpline(polyhedron=c_0, function=f_0),
                       plf.SubSpline(polyhedron=c_1, function=f_1),
                       plf.SubSpline(polyhedron=c_2, function=f_2)])


def formats_polyhedrons_1():
    """ formats partition examples for the computation of permissiveness for
    the timed automata with two clocks x,y
    and two linear transition 0 <= x <= 1, 0 <= y <= 1"""
    c_0 = [x >= 0, y >= 0, x <= 1, 2 * y >= 1 + x, y <= 1]
    c_1 = [x >= 0, y >= 0, x <= 1, 2 * y <= 1 + x, y <= 1]

    return ppl.C_Polyhedron(lin_alg.cs_from_list(c_0)), \
           ppl.C_Polyhedron(lin_alg.cs_from_list(c_1))


def formats_functions_1():
    """ formats functions examples for the computation of permissiveness for
    the timed automata with two clocks x,y and two linear transition 0 <= x
    <= 1, 0 <= y <= 1"""
    f_0 = rla.RationalLinearExpression((0, -1), 1)
    f_1 = rla.RationalLinearExpression((Fraction(-1, 2), 0), Fraction(1, 2))
    return f_0, f_1


def spline_1():
    c_0, c_1 = formats_polyhedrons_1()
    f_0, f_1 = formats_functions_1()
    return plf.Spline([plf.SubSpline(polyhedron=c_0, function=f_0),
                       plf.SubSpline(polyhedron=c_1, function=f_1)])


"""
timed_automaton_creator({
        "transitions": [
            {
                "start_location": 0,
                "end_location": 1,
                "data": [{
                    "action": "a",
                    "guard": {
                        "type": "linear",
                        "constraints": [
                            {
                                "type": "linear",
                                "data": {
                                    P: ppl.C_Polyhedron()
                                },
                            },
                            {
                                "type": "linear",
                                "data": {
                                    Q : ppl.C_polyhedron()
                                },
                            }
                        ]
                    },
                    "resets": [1],
                }]
            },
        ],
        "init_location": 0,
        "goal_location": 1,
        "number_clocks": 2,
    })
    """


def two_transition_timed_auto(guard_0: ppl.C_Polyhedron, reset_0: guards.Reset,
                              guard_1: ppl.C_Polyhedron,
                              reset_1: guards.Reset) -> timedauto.TimedAutomaton:
    return creators.timed_automaton_creator(
        {
            "transitions": [{
                "start_location": 0,
                "end_location": 1,
                "data": [
                    {
                        "action": "a_0",
                        "guard": guard_0,
                        "resets": reset_0
                    }]
            },
                {
                    "start_location": 1,
                    "end_location": 2,
                    "data": [
                        {
                            "action": "a_1",
                            "guard": guard_1,
                            "resets": reset_1
                        }
                    ]
                }],
            "init_location": 0,
            "goal_location": 2,
            "number_clocks": 2,

        },

    )


def three_transitions_timed_auto(
        guard_0: ppl.C_Polyhedron,
        reset_0: guards.Reset,
        guard_1: ppl.C_Polyhedron,
        reset_1: guards.Reset,
        guard_2: ppl.C_Polyhedron,
        reset_2: guards.Reset
) -> timedauto.TimedAutomaton:
    return creators.timed_automaton_creator(
        {
            "transitions": [{
                "start_location": 0,
                "end_location": 1,
                "data": [
                    {
                        "action": "a_0",
                        "guard": guard_0,
                        "resets": reset_0
                    }]
            },
                {
                    "start_location": 1,
                    "end_location": 2,
                    "data": [
                        {
                            "action": "a_1",
                            "guard": guard_1,
                            "resets": reset_1
                        }
                    ]
                },
                {
                    "start_location": 2,
                    "end_location": 3,
                    "data": [
                        {
                            "action": "a_2",
                            "guard": guard_2,
                            "resets": reset_2
                        }
                    ]
                }
            ],
            "init_location": 0,
            "goal_location": 3,
            "number_clocks": 2,

        },

    )


def formats_timed_automaton_0() -> timedauto.TimedAutomaton:
    guard = ppl.C_Polyhedron(
        lin_alg.cs_from_list([x >= 0, x <= 1, y >= 0, y <= 1]))

    reset = []
    return two_transition_timed_auto(guard_0=guard, reset_0=reset,
                                     guard_1=guard, reset_1=reset)


def formats_timed_automaton_1() -> timedauto.TimedAutomaton:
    guard_0 = ppl.C_Polyhedron(
        lin_alg.cs_from_list([x >= 0, x <= 1, y >= 0, y <= 1]))
    reset_0 = [1]
    guard_1 = ppl.C_Polyhedron(
        lin_alg.cs_from_list([x >= 1, y >= 1, x <= 2, y <= 2]))
    reset_1 = []
    return two_transition_timed_auto(guard_0=guard_0, reset_0=reset_0,
                                     guard_1=guard_1, reset_1=reset_1)


def formats_timed_automaton_2() -> timedauto.TimedAutomaton:
    guard_0 = ppl.C_Polyhedron(
        lin_alg.cs_from_list([x >= 0, x <= 1, y >= 0, y <= 1]))
    reset_0 = [1]
    guard_1 = ppl.C_Polyhedron(
        lin_alg.cs_from_list([x >= 1, y >= 0, x <= 2, y <= 1]))
    reset_1 = []
    return two_transition_timed_auto(guard_0=guard_0, reset_0=reset_0,
                                     guard_1=guard_1, reset_1=reset_1)


def formats_timed_automaton_3() -> timedauto.TimedAutomaton:
    guard_0 = ppl.C_Polyhedron(lin_alg.cs_from_list([y >= 0, y <= 1]))
    reset_0 = [1]
    guard_1 = ppl.C_Polyhedron(
        lin_alg.cs_from_list([x >= 1, y >= 0, x <= 2, y <= 1]))
    reset_1 = []
    return two_transition_timed_auto(guard_0=guard_0, reset_0=reset_0,
                                     guard_1=guard_1, reset_1=reset_1)


def formats_timed_automaton_4() -> timedauto.TimedAutomaton:
    guard_0 = ppl.C_Polyhedron(lin_alg.cs_from_list([x >= 0, x <= 2, y >= 0,
                                                     y <= 2]))
    reset_0 = [1]
    guard_1 = ppl.C_Polyhedron(
        lin_alg.cs_from_list([x >= 0, y >= 0, x <= 4, y <= 2]))
    reset_1 = []
    guard_2 = ppl.C_Polyhedron(
        lin_alg.cs_from_list([y >= 0, y <= 1]))
    reset_2 = []
    return three_transitions_timed_auto(
        guard_0=guard_0, reset_0=reset_0,
        guard_1=guard_1, reset_1=reset_1,
        guard_2=guard_2, reset_2=reset_2
    )

def formats_timed_automaton_5() -> timedauto.TimedAutomaton:
    guard_0 = ppl.C_Polyhedron(lin_alg.cs_from_list([
        x >= 0, x <= 1,
    ]))
    reset_0 = []
    guard_1 = ppl.C_Polyhedron(lin_alg.cs_from_list([
        y >= 0, y <= 1
    ]))
    reset_1 = [1]
    guard_2 = ppl.C_Polyhedron(lin_alg.cs_from_list([
        x >= 1, x <= 2,
        y >= 0, y <= 1,
    ]))
    reset_2 = []
    return three_transitions_timed_auto(
        guard_0=guard_0, reset_0=reset_0,
        guard_1=guard_1, reset_1=reset_1,
        guard_2=guard_2, reset_2=reset_2
    )


def formats_timed_automata_0_n_transitions(n: int) -> timedauto.TimedAutomaton:
    guard = ppl.C_Polyhedron(
        lin_alg.cs_from_list([x >= 0, x <= 1, y >= 0, y <= 1]))
    reset = []

    t = [{
        "start_location": i,
        "end_location": i + 1,
        "data": [
            {
                "action": "a_0",
                "guard": guard,
                "resets": reset
            }]
    } for i in range(n)]
    return creators.timed_automaton_creator(
        {
            "transitions": t,
            "init_location": 0,
            "goal_location": n,
            "number_clocks": 2,

        },

    )


@dataclasses.dataclass
class AssociationFunctionPolyhedra:
    function: plf.LinearFunction
    poly_list: List[ppl.C_Polyhedron]


def timed_automata_n_transitions(nb_c: int, guards: List[ppl.C_Polyhedron],
                                 resets: List[List[int]]) -> timedauto.TimedAutomaton:
    n = len(guards)
    t = [{
        "start_location": i,
        "end_location": i + 1,
        "data": [
            {
                "action": "a_0",
                "guard": guards[i],
                "resets": resets[i]
            }]
    } for i in range(n)]
    return creators.timed_automaton_creator(
        {
            "transitions": t,
            "init_location": 0,
            "goal_location": n,
            "number_clocks": nb_c,

        },

    )


def experiment_clocks(nb_transition: int, nb_clock: int):
    variables = [Variable(i) for i in range(nb_clock)]

    guard = polyhedra.c_polyhedron_constructor(
        dimension=nb_clock,
        cs_list=
        [variables[i] >= 0 for i in range(nb_clock)]
        + [variables[i] <= 1 for i in range(nb_clock)]

    )

    guards = [guard for i in range(nb_transition)]
    resets = [[] for i in range(nb_transition)]
    # resets = [[0], [], [1], []]
    nb_reset = 0

    while nb_reset <= nb_clock // 2:
        t = random.randint(0, nb_transition - 1)
        c = random.randint(0, nb_clock - 1)

        if len(resets[t]) == 0:
            resets[t].append(c)
            nb_reset += 1

    ta = timed_automata_n_transitions(nb_c=nb_clock,
                                      guards=guards,
                                      resets=resets)

    print(f"Resets : {resets}")

    ex = explorer.Explorer(ta=ta)

    t0 = time.time()
    ex.explorer()
    t1 = time.time()

    filepath = Path('./../../experiment_random_resets.txt')

    with open(filepath, 'w') as file:
        file.write(f"resets : {resets}\n"
                   f"runtime : {t1 - t0}\n")

        file.write("\n\n----- PERMISSIVENESS FUNCTIONS -----\n\n")
        for node in ta.nodes:
            spline = ta.nodes[node]['Permissiveness_function']
            file.write(f"Node : {node}\n")
            visual_spline = visual(spline)

            file.write("[\n")
            for sub_spline in visual_spline:
                file.write(f"{sub_spline},\n")
            file.write("]\n\n")

    pass


def experiment_tikz():
    tas = [
        # formats_timed_automaton_5(),
        formats_timed_automaton_0(),
        formats_timed_automaton_2(),
        formats_timed_automaton_3()
    ]

    for i, ta in enumerate(tas):
        print(i)
        ex = explorer.Explorer(ta=ta)
        perm = ex.all_nodes_permissiveness()[0]

        for j, perm in ex.all_nodes_permissiveness().items():
            if j != 2:
                latex_output(perm, Path(f'./test_{i}_l{j}.tex').absolute())

def experiment_long_computation():
    guard = polyhedra.c_polyhedron_constructor(
        dimension=3,
        cs_list=
        [ppl.Variable(i) >= 0 for i in range(3)]
        + [ppl.Variable(i) <= 1 for i in range(3)]

    )
    guards = [guard for i in range(4)]

    # resets_0 = [[], [1], [], [2]]
    resets_1 = [[2], [1], [], []]
    ta = timed_automata_n_transitions(nb_c=3,
                                      guards=guards,
                                      resets=resets_1)

    print(f"Resets : {resets_1}")

    ex = explorer.Explorer(ta=ta)

    t0 = time.time()
    ex.explorer()
    t1 = time.time()
    perm = ex.permissiveness()
    filepath = Path('./../../experiment_random_resets.txt')

    with open(filepath, 'w') as file:
        file.write(f"resets : {resets_1}\n"
                   f"runtime : {t1 - t0}\n")

    return perm





def main():
    x = Variable(0)
    y = Variable(1)
    t = Variable(2)

    # print(plf.Spline.list_const_and_sym_constr(spline=spline_0()))
    # new_spline = spline_1().partition(spline_0())
    # print(len(new_spline.sub_splines))
    #
    # for sspline in new_spline.sub_splines:
    #     print("new:")
    #     print(sspline.polyhedron.constraints())
    #     print(sspline.function)
    #
    # new_spline_1 = spline_0().partition(spline_1())
    # print(len(new_spline_1.sub_splines))
    #
    # for sspline in new_spline_1.sub_splines:
    #     print("new:")
    #     print(sspline.polyhedron.constraints())
    #     print(sspline.function)
    #
    # poly_new = formats_polyhedrons_0()
    # for poly in poly_new:
    #     print(poly.minimized_constraints())
    #
    # new_spline = spline_0().maximum(other=spline_1())
    # print(new_spline)
    # print(len(new_spline.sub_splines))
    # print("Test maximum")
    # for sspline in new_spline.sub_splines:
    #     print(sspline.polyhedron.constraints())
    #     print(sspline.function)
    #
    # gen = new_spline.sub_splines[0].polyhedron.generators()
    # print(gen[0].coefficients())
    #
    # min_spline = spline_0().minimum(spline_1())
    # print("Test minimum")
    # for sspline in min_spline.sub_splines:
    #     print(sspline.polyhedron.constraints())
    #     print(sspline.function)

    # TEST OF TECHNICAL LEMMA

    # f=
    # g=

    # polyhedra.show(new_spline.sub_splines[0].polyhedron)

    # TEST OF THE EXPLORER

    ta = formats_timed_automata_0_n_transitions(2), \
         formats_timed_automaton_1(), \
         formats_timed_automaton_2(), \
         formats_timed_automaton_3(), \
         formats_timed_automaton_4()
    # TIMING RUNTIME TEST
    """
    timings = []

    for i in [1, 2, 3, 4]:
        ex = explorer.Explorer(ta=ta[i])

        t0 = time.time()
        ex.explorer()
        t1 = time.time()

        timings.append(t1 - t0)

    pprint.pprint({i: t for i,t in enumerate(timings)})
    """
    # PRINT THE SPLINE ON TIKZ INTERFACE
    ex = explorer.Explorer(ta=ta[3])
    ex, perm = ex, ex.all_nodes_permissiveness()[1]

    latex_output(perm, Path('./test.tex').absolute())
    """

    # CONVEX HULL TEST

    ta = formats_timed_automaton_4()
    ex = explorer.Explorer(ta=ta)
    ex, perm = ex, ex.all_nodes_permissiveness()[0]

    # L: List[AssociationFunctionPolyhedra] = []
    # for sub in perm.sub_splines:
    #     P = sub.polyhedron
    #     f = sub.function
    #     found_assoc = False
    #
    #     for assoc in L:
    #         if f.is_equal_to(assoc.function):
    #             assoc.poly_list.append(P)
    #             found_assoc = True
    #             break
    #
    #     if not found_assoc:
    #         assoc = AssociationFunctionPolyhedra(function=f, poly_list=[P])
    #         L.append(assoc)

    # first_assoc = L[3]
    # print(first_assoc.poly_list)
    # P_init = ppl.C_Polyhedron(2, 'empty')
    # for poly in first_assoc.poly_list:
    #     P_init.poly_hull_assign(poly)
    # print(P_init.constraints())

    # ex = explorer.Explorer(ta=ta[4])
    #
    # index = 0
    #
    # ex, perm = ex, ex.all_nodes_permissiveness()[0]
    # print(perm)
    #
    # latex_output(perm, Path('./test.tex').absolute())
    #
    # print(len(perm.sub_splines))
    # for sub in perm.sub_splines:
    #     print(sub.polyhedron.generators())
    #     print(sub.function)

    # #
    print("Permissiveness of all nodes")

    # Time complexity evaluation
    """
    """
    max_transition = 1_000_000
    step = 1_000
    timing_tab = {}
    for i in range(step, max_transition+step, step):
        # start of operation
        timing_tab[i] = timeit.timeit(
            "explorer.Explorer(formats_timed_automata_0_n_transitions( " + str(i) + "))",
            number=1,
            globals=globals())
        # print(bt)
        # pprint.pprint(explorer.Trace.compute_trace_permissiveness(bt))
        # end of operation

    pprint.pprint( timing_tab )
    """
    """
    # perm = explorer.Explorer(formats_timed_automata_0_n_transitions(7)).permissiveness()
    # for sub in perm.sub_splines:
    #     print(sub.polyhedron.constraints())
    #     print(sub.function)

    #
    # perms = ex.all_nodes_permissiveness()
    # for node in ta[index].nodes:
    #     perm = perms[node]
    #     print("------------------------")
    #     print("permissiveness of node " + str(node))
    #     for i, sub in enumerate(perm.sub_splines):
    #         print("--------------")
    #         print("Piece number " + str(i))
    #         print(sub.polyhedron.constraints())
    #         print(sub.function)

    # entry_poly = polyhedra.c_polyhedron_constructor(2, [x >=0, y>=0])
    # guard = polyhedra.c_polyhedron_constructor(2, [ x >=0 , x<=1, y>=0, y<=1])
    #
    # max_int = var_elim.maximal_intervals(poly_p=entry_poly, poly_q=entry_poly, guard=guard, resets=[], dim=2)

    # g = ppl.C_Polyhedron(lin_alg.cs_from_list([x >= 1, x <= 2, y >=0, y <= 1]))
    # guard_poly = var_elim.guard_poly(guard=g, dim=2)
    # print(guard_poly.constraints())

    # P = polyhedra.c_polyhedron_constructor(dimension=2, cs_list=[-x + 2 >= 0,
    #                                                              x - y >= 0,
    #                                                              y - 1 >= 0])
    # R = polyhedra.c_polyhedron_constructor(dimension=2, cs_list=[x >= 0,
    #                                                              x >= y,
    #                                                              y <= 1])
    # S = polyhedra.c_polyhedron_constructor(dimension=2, cs_list=[- x + 2 >=
    #                                                              0,
    #                                                              y - x >= 0,
    #                                                              y - 1 >= 0])
    # print("P = " + str(P.generators()))
    # n = 2
    # r = [1]
    #
    # Q = var_elim.entry_set(poly_p=R, poly_q=R, resets=r, dim=n)
    # print("Q non minimized constraints = "+str(Q.constraints()))
    # print("Q minimized constraints = "+str(Q.minimized_constraints()))
    # print("Q generators = "+str(Q.generators()))
    """

if __name__ == "__main__":
    # execute only if run as a script
    # main()
    # experiment_tikz()
    # TEST OF THREE CLOCKS TA (RUNTIME HIGH !!!) Random TA.
    # experiment_clocks(4, 3)
    perm = experiment_long_computation()