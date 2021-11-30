# coding=utf-8
"""
==================================================
Timed Automaton module
==================================================
This module provide a number of objects and methods useful to manipulate timed
automaton.

Classes:
------
Transitions
TimedAutomaton

Methods:
TODO: enumerate all methods.
------
"""

from __future__ import annotations  # For forward reference typing

import math
from collections import namedtuple
from copy import copy
from operator import eq
import networkx as nx

from pysymbrobustness.dtype import Location, Action
from typing import List, Set, Optional, Dict, Tuple
import matplotlib.pyplot as plt
import pysymbrobustness.ta.guards as guards

plt.switch_backend("Qt5Agg")

Edge = namedtuple("Edge", ["start_location", "end_location", "data"])
Configuration = namedtuple("Configuration", ["location", "valuation"])
Label = guards.Label
Data = Dict[Action, Label]
DiGraphEdge = Tuple[Location, Location, Data]


def edge_form(edge: Edge) -> DiGraphEdge:
    """Convert an edge into a diGraph-friendly form

    :param
    ----------
    edge: an edge as defined in the nametuple Edge

    :return:
    ----------
    e: an edge of a graph.
    """
    return edge.start_location, edge.end_location, edge.data


def edges_form(edges_list: List[Edge]) -> List[DiGraphEdge]:
    """Convert all the edges into a DiGraph-friendly form

    :param
    ----------
    edges_list: a list of Edges.

    :return:
    ----------
    e: a list of edges.
    """
    return list(map(edge_form, edges_list))


class TimedAutomaton(nx.DiGraph):
    """
    Timed Automaton class
    """

    class _TAFlags(object):
        """
        Immutable key flags dictionary
        Flags:
        well_formed: To Be documented, check basic property needed to have a
        well formed TA.
        single_action: True if For a label 'action' and a location l,
        the timed automaton has at one edge, starting from l, labeled 'action'
        branch_free: True if the timed automaton has at most one edge for each
        location.
        acyclic: True if the timed automaton does not contain any cycle (
        loop transition from a location to itself)

        """
        _flags = ["well_formed", "single_action",
                  "branch_free", "acyclic"]
        __slots__ = ["well_formed", "single_action",
                     "branch_free", "acyclic", "overwrite"]

        def __init__(self, overwrite: Optional[Set[str]]):
            self.overwrite = overwrite
            self.reset()

        def reset(self):
            for flag in self._flags:
                if flag in self.overwrite:
                    self.__setattr__(flag, True)
                else:
                    self.__setattr__(flag, False)

    def __init__(self, transitions: List[Edge], init_location: Location,
                 goal_location: Location,
                 number_clocks: int, overwrite: Optional[Set[str]] = None):
        super().__init__()
        if overwrite is None:
            overwrite = []
        self.add_edges_from(edges_form(transitions))
        self.init_location = init_location
        self.goal_location = goal_location
        self.number_clocks = number_clocks
        self._flags = self._TAFlags(overwrite)

        # self.is_well_formed()

    # Methods:

    def _invalidate_flags_on_change(self):
        self._flags.reset()

    def __eq__(self, other: TimedAutomaton) -> bool:
        """Compare two ta"""
        if type(other) != TimedAutomaton:
            return False

        return (nx.is_isomorphic(self, other, node_match=eq,
                                 edge_match=eq) and
                self.init_location == other.init_location and
                self.goal_location == other.goal_location and
                self.number_clocks == other.number_clocks)
        # and self.is_well_formed() == other.is_well_formed()

    def __str__(self) -> str:  # pragma: no cover
        return "The timed automaton has " + str(self.number_clocks) + \
               " clocks, " + str(len(self.edges)) + \
               " transitions and " + str(len(self.nodes)) + \
               " locations, (l_i)_i and ends at l_" + \
               str(len(self.nodes)) + ". Its locations are :" + str(
            self.nodes) + " and its transitions are: " + str(
            self.edges)

    def __repr__(self) -> str:  # pragma: no cover
        """ Print the information about the transition """
        return self.__str__()

    def graphical_print(self):  # pragma: no cover
        """
        Print with graphic representation of all locations in the timed
        automaton, and the non-label transitions (transitions where guards
        are not printed.)
        """
        nx.draw(self, with_labels=True)
        plt.show()

    def add_transition(self, start: Location, end: Location, label: Label,
                       action: Action):
        """ Add a transition to the directed graph of the timed automaton.

        :param
        ----------
        transition: a Transition. """
        self.add_edge(start, end, **{action: label})
        self._invalidate_flags_on_change()

    def add_transitions(self, edges: List[Edge]):
        """ Add a transition to the directed graph of the timed automaton.

        :param
        ----------
        edges: a Edge list. """
        self.add_edges_from(edges_form(edges))
        self._invalidate_flags_on_change()
        # Warning to put if a transition between these two location already
        # exists: it erase the existing one.
        # Clean le precedent label

    def change_start_location(self, location: Location):
        """change the starting location of the timed automaton.

        :param
        ----------
        location: an integer, that represents the starting location.
        """
        self.init_location = location
        self._invalidate_flags_on_change()

    def change_end_location(self, location: Location):
        """change the end location of the timed automaton.

        :param
        ----------
        location: an integer, that represents the end location.
        """
        self.goal_location = location
        self._invalidate_flags_on_change()

    def change_clock_index(self, index: int):
        """ TODO: document and change this fct
        """
        self.number_clocks = index
        self._invalidate_flags_on_change()

    # def is_well_formed(self) -> bool:
    #     """check if the timed automaton is well-formed and if so switch the
    #     well-formed flag.
    #     """
    #     if self._flags.well_formed:
    #         return True
    #
    #     if self.init_location not in self.nodes:
    #         raise exceptions.LocationNotFoundException(self.init_location, None)
    #     if self.goal_location not in self.nodes:
    #         raise exceptions.LocationNotFoundException(None, self.goal_location)
    #     if self.number_clocks < 1:
    #         raise exceptions.ClockNotFoundException
    #
    #     for (s_location, f_location) in self.edges:
    #         # NB: like this : transition.check_guards(self.number_clocks)
    #         for label in self[s_location][f_location].values():
    #             label.well_formed(self.number_clocks)
    #
    #     self._flags.well_formed = True
    #     return True

    def is_single_action(self) -> bool:
        """
        This function verifies for each nodes, that every edges had
        different action label.
        :return:
        """
        for source in self.nodes:
            actions = set()
            for value in self[source].values():
                for action in value.keys():
                    if action in actions:
                        self._flags.single_action = False
                        return False
                    else:
                        actions.add(action)
        self._flags.single_action = True
        return True

    def is_acyclic(self) -> bool:
        if self._flags.acyclic:
            return True
        boolean = len(list(nx.simple_cycles(nx.DiGraph(self.edges())))) != 0
        if boolean:
            self._flags.acyclic = False
            return False
        else:
            self._flags.acyclic = True
            return True

    def is_branch_free(self) -> bool:
        for source in self.nodes:
            if len(list(self.successors(source))) > 1:
                self._flags.branch_free = False
                return False
        self._flags.branch_free = True
        self._flags.acyclic = True
        self._flags.single_action = True
        return True
