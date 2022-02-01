# coding=utf-8
"""
==================================================
explorer
==================================================
This module provides a bakward algorithm, presentend in FORMSTS 2020 by Emily Clement, Thierry Jeron, Nicolas Markey and
David Mentre, to compute the permissiveness of a timed automata

Classes:
------

Methods:
------
"""

import itertools
from functools import reduce
from itertools import product
import math
from pathlib import Path
from typing import Optional, List, NamedTuple, Union, Dict, Tuple

import ppl as ppl

from pysymbrobustness.ppl_extension import polyhedra
from pysymbrobustness.ta.misc_debug import LatexSplineVisualisation
from pysymbrobustness.ta.timedauto import TimedAutomaton, Configuration
import networkx as nx
import pysymbrobustness.ppl_extension.piecewise_linear_function as plf
import pysymbrobustness.ppl_extension.variable_elimination as var_elim
import pysymbrobustness.ppl_extension.technical_lemma as tech
import pysymbrobustness.ppl_extension.rational_linear_algebra as rla

Location = Union[int, str]
Permissiveness = plf.Spline
PERMISSIVENESS_FUNCTION = "Permissiveness_function"


class Explorer(object):
    """
    TODO
    """

    def __init__(self, ta: TimedAutomaton,
                 to_print: bool = False):
        self.ta = ta
        self.to_print = to_print
        self.executed = False

    def detect_cycle(self) -> List[List[Optional[Location]]]:
        """
        This function take a TA and gives all the basic cycle of this TA
        Tips: use networkx
        """
        return nx.cycle_basis(self.ta)

    def compute_permissiveness(self, location: Location) -> Permissiveness:

        successors = self.ta.successors(location)
        dict_permissiveness = {suc: self.compute_permissiveness_suc(location, suc) for suc in successors}

        return plf.Spline.maximum_of_splines(spline_list=list(dict_permissiveness.values()))

    def compute_permissiveness_suc(self, location: Location, successor: Location) -> Permissiveness:

        to_plot = False

        if to_plot:
            full_perm_visual = LatexSplineVisualisation(Path('./full_permissiveness.tex'))
        # Get the future permissiveness of successor. Type: Permissiveness (=plf.Spline)
        future_perm: plf.Spline = self.ta.nodes[successor][PERMISSIVENESS_FUNCTION]

        # Get the transition between location and successors and the label/reset/guard and number of clocks
        edge = self.ta[location][successor]
        label = list(edge.values())[0]
        reset = label.resets
        guard = label.guard
        n = self.ta.number_clocks  # ppl.Variable(0)....ppl.Variable(n-1) correspond to the clock of the automata

        # Name two variable alpha and beta whose index are just after the clock's indexes
        # alpha = ppl.Variable(n)
        # beta = ppl.Variable(n + 1)

        # a list containing the permissiveness of location, which maximizes  the permissiveness of this list.
        to_max_permissiveness = []

        # Let's denote suc(d) the configuration of the successor of location when applying the delay "d".
        # future_perm is partitioned into multiple polyhedron and the successors end up in one of them
        # We iterate on the cartesian product of this partition. For each couple of polyhedron, we consider
        # suc(alpha) and suc(beta) in this couple of polyhedron and compute the minimum of their permissiveness.
        # We compute this list to_max_permissiveness in the following loop:
        # Typ: List of Permissiveness (=List[plf.Spline])
        # Let's chose two couple of (h_alpha, f) and (h_beta, g)

        # Applying Fourier-Motzkin algorithm to determine the valuation v s.t there exists an alpha and a beta s.t
        # alpha < beta and v+alpha and v+beta fit the guard's constraints.
        entry_poly_base = var_elim.guard_poly(guard=guard, dim=n)

        for n_sub_debug, (sub_alpha, sub_beta) in enumerate(product(future_perm.sub_splines, repeat=2)):

            if to_plot:
                sub_visual = LatexSplineVisualisation(Path(f'./{n_sub_debug}_permissiveness.tex'))
            # Precising the type of sub_alpha and sub_beta for static analysis.
            sub_alpha: plf.SubSpline
            sub_beta: plf.SubSpline

            # Applying Fourier-Motzkin algorithm to determine the valuation v s.t. there exist an alpha and a beta s.t.
            # alpha < beta and v+alpha[r] is in the polyhedron of sub_alpha (resp v+beta[r] / sub_beta).
            entry_poly: ppl.C_Polyhedron = var_elim.entry_set(poly_p=sub_alpha.polyhedron, poly_q=sub_beta.polyhedron,
                                                              resets=reset, dim=n)
            entry_poly.intersection_assign(entry_poly_base)
            # entry_poly.intersection_assign(guard)

            # Applying partially Fourier-Motzkin algorithm to determine all the possible alpha (resp beta)
            #  for a fixed arbitrary v in entry_poly.
            # Type of max_interval: a tuple of four Splines. The first (resp last) two splines are the bounds of the
            # interval of possible alpha (resp beta).
            # For each n-dimension polyhedron, the value of the bound is the associated plf.LinearFunction
            max_interval = var_elim.maximal_intervals(poly_p=sub_alpha.polyhedron, poly_q=sub_beta.polyhedron,
                                                      resets=reset, guard=guard,
                                                      dim=n)  # NB: Check the max_interval, suspicious on the x1 variable (infinite while blocked by guard)

            # Get the value of f = alpha*a + b, g = beta*c + d which represent the value of the permissiveness of
            # the suc(alpha) and sub(beta)

            f = self.next_valuation(fct=sub_alpha.function, reset=reset)
            g = self.next_valuation(fct=sub_beta.function, reset=reset)

            # Initiate the permissiveness list for each polyhedron where the bound are not piecewise-affine function,
            # but affine functions.
            to_max_max_permissiveness = []
            # A tuple of four elements, each element corresponds to the SubSplines of the associated Spline.
            sub_splines = tuple(map(lambda x: x.sub_splines, max_interval))

            # *subsplines transforms the tuple into four arguments (each are a list of SubSplines)
            for n_debug, entry_set in enumerate(product(*sub_splines)):
                # entry_set = sub_alpha_min, sub_alpha_max, sub_beta_min, sub_beta_max
                entry_set: Tuple[plf.SubSpline, plf.SubSpline, plf.SubSpline, plf.SubSpline]
                permissiveness: plf.Spline = tech.plf_optimization(f=f, g=g, entry_poly=entry_poly, entry_set=entry_set)
                # TODO: if not branch-free, add something...
                if len(permissiveness.sub_splines) > 0:
                    to_max_max_permissiveness.append(permissiveness)

                    if to_plot:
                        sub_visual.add_spline(permissiveness, caption=f"{n_debug}th entry set")

            # Compute to_max_permissiveness, the permissiveness that maximizes the permissiveness among all
            # the possible [alpha,beta]
            maximum_spline = plf.Spline.maximum_of_splines(to_max_max_permissiveness)

            if to_plot:
                sub_visual.add_spline(maximum_spline, caption=f"fusion spline")
                sub_visual.output()
            if len(maximum_spline.sub_splines) > 0:
                to_max_permissiveness.append(maximum_spline)

                if to_plot:
                    full_perm_visual.add_spline(maximum_spline, caption=f"{n_sub_debug}th spline")

        # Compute the permissiveness that is maximized among the choice of the zone where (l,v+alpha[r]) (resp beta)
        # will end up in.
        permissiveness: plf.Spline = plf.Spline.maximum_of_splines(spline_list=to_max_permissiveness)

        if to_plot:
            full_perm_visual.add_spline(permissiveness, caption=f"Global permissiveness")
            full_perm_visual.output()
        # Return the final permissiveness, type: Permissiveness (plf.Spline)
        return permissiveness

    # def fct_plus_delay(self, delay: ppl.Variable, fct: plf.LinearFunction) -> plf.LinearFunction:
    #     fct_linear = type(fct) == rla.RationalLinearExpression
    #     n = self.ta.number_clocks
    #     if delay.id() < n:
    #         raise IndexError(str(delay) + "'s index should be superior to the clock indexes")
    #     if type(fct) != rla.InfiniteExpression and not fct_linear:
    #         raise TypeError
    #     if not fct_linear:
    #         return fct
    #
    #     sum_coeff = sum(fct.coefficients())
    #     coeffs = [0 for i in range(delay.id())] + [sum_coeff]
    #     delay_term = rla.RationalLinearExpression(tuple(coeffs))
    #     return fct + delay_term

    def next_valuation(self, fct: plf.LinearFunction, reset: List[int]) -> plf.LinearFunction:
        fct_linear = type(fct) == rla.RationalLinearExpression
        if type(fct) != rla.InfiniteExpression and type(fct) != rla.RationalLinearExpression:
            raise TypeError
        if not fct_linear:
            return fct

        reseted_fct = rla.RationalLinearExpression(homogeneous_terms=tuple(0 if index in reset
                                                                           else fct.coefficient(ppl.Variable(index))
                                                                           for index in range(fct.space_dimension())),
                                                   inhomogeneous_term=fct.inhomogeneous_term())

        return reseted_fct

    def start_value(self, node: Location) -> Permissiveness:
        poly = polyhedra.c_polyhedron_constructor(dimension=self.ta.number_clocks)
        return plf.Spline(plf.SubSpline(polyhedron=poly,
                                        function=rla.InfiniteExpression(sign=node == self.ta.goal_location)))

    def explorer(self):
        """
        Basic explorer for general graph.
        TODO: an optimisation with topological sort and recursive treatment of cycles.
        """

        if self.executed:
            return

        if self.ta.is_acyclic():
            self.acyclic_explorer()
            return

        # Initialisation the permissiveness sequence
        for node in self.ta.nodes:
            node[PERMISSIVENESS_FUNCTION] = self.start_value(node)

        # Computation of fixed point of the permissiveness sequence for each nodes
        nodes = [self.ta.goal_location]
        while len(nodes) != 0:
            current_node = nodes.pop()
            Perm = self.compute_permissiveness(current_node)
            if not self.ta.nodes[current_node][PERMISSIVENESS_FUNCTION].is_equal_to(Perm):  # Search for the fixed point
                self.ta.nodes[current_node][PERMISSIVENESS_FUNCTION] = Perm
                for p in self.ta.predecessors(current_node):
                    nodes.append(p)
        self.executed = True

    def compute_final_permissiveness(self, location: Location):
        self.explorer()
        return self.ta.nodes[location][PERMISSIVENESS_FUNCTION]

    def all_nodes_permissiveness(self):
        self.explorer()
        dict_permissiveness = {node: self.compute_final_permissiveness(location=node) for node in self.ta.nodes}
        return dict_permissiveness

    def permissiveness(self):
        return self.compute_final_permissiveness(location=self.ta.init_location)

    def acyclic_explorer(self):
        """
        Basic explorer for acyclic graph with a (reverse) topological sort
        """

        # Do not execute the explorer if it has already been executed
        if self.executed:
            return

        # Initialisation of the permissiveness value for each nodes

        for node in self.ta.nodes:
            self.ta.nodes[node][PERMISSIVENESS_FUNCTION] = self.start_value(node)

        # Reverse the graph to apply a classic topological sort
        simple_ta = nx.DiGraph(self.ta.edges())
        new_ta = simple_ta.reverse(copy=False)
        sorted_nodes = nx.topological_sort(new_ta)

        # Computation of the final permissiveness

        for current_node in sorted_nodes:
            if current_node != self.ta.goal_location:
                import time
                t0 = time.time()
                Perm = self.compute_permissiveness(current_node)
                t1 = time.time()
                print(f"{current_node} : {t1 - t0:<10.6f} | {len(Perm.functions)}")
                self.ta.nodes[current_node][PERMISSIVENESS_FUNCTION] = Perm

        self.executed = True
