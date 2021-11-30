# coding=utf-8
from __future__ import annotations  # For forward reference typing
from typing import List

import ppl as ppl
import pysymbrobustness.ta.timedauto as timed_auto
import pysymbrobustness.ta.guards as guards
import pysymbrobustness.ta.exceptions as exceptions

"""
Module to generate the object of the ta library from JSON-like format
(linearConstraint, LinearGuard,Label, Edge, and TimedAutomaton)
"""


def label_creator(guard_data: ppl.C_Polyhedron, resets: List[int]):
    """
    Generate from a JSON-like syntax a Label.
    :param: guard_data: a list of dict (see linear_guard_creator for the syntax)
    :param: resets: a list of integers
    :return: a Label
    :type: "guard_data" should be a guard creator type, and resets should be a list of integer.
    ------
    Example:
        label_creator(
                        guard_data={
                            "type": "linear",
                            "constraints": [{
                                "type": "linear",
                                "data": {
                                    'clock_index': 0,
                                    'lower_bound': 0,
                                    'upper_bound': 1,
                                }
                            }]
                        },
                        resets=[]
                    )
    """
    guard = guard_data

    if guard is None:
        raise exceptions.UnknownGuardType(guard_data["type"]) # NB Verify that

    return guards.Label(
        guard=guard,
        resets=resets
    )


def edge_creator(serialized):
    """
    Generate from a JSON-like syntax a linearGuard.
    :param: serialized: a dict
    :return: an edge
    :type: "start_location" and "end_location" should be integers or str,
            "data" should be a label creator type
            and resets should be a list of integer.
    ------
    Example:
        edge_creator(
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
                                    "lower_bound": 1,
                                    "upper_bound": 2,
                                    "clock_index": 0,
                                },
                            },
                            {
                                "type": "linear",
                                "data": {
                                    "lower_bound": 0,
                                    "upper_bound": 1,
                                    "clock_index": 1,
                                },
                            }
                        ]
                    },
                    "resets": [1],
                }]
            }
                    )
    """
    return timed_auto.Edge(
        start_location=serialized["start_location"],
        end_location=serialized["end_location"],
        data={
            data["action"]: label_creator(
                guard_data=data["guard"],
                resets=data["resets"]) for data in serialized["data"]
        },
    )


def timed_automaton_creator(serialized):
    """
    Generate from a JSON-like syntax a linearGuard.
    :param serialized: a dict
    :return: a TimedAutomaton
    :type: : "Transitions" should contain edge_creator input type
            "init_location" and "goal_location" should contain integer or str
            and "number_clocks" should contains integer.
    -----
    Example:

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
    overwrite = serialized.get("overwrite", None)

    return timed_auto.TimedAutomaton(
        transitions=[edge_creator(t) for t in serialized["transitions"]],
        init_location=serialized["init_location"],
        goal_location=serialized["goal_location"],
        number_clocks=serialized["number_clocks"],
        overwrite=overwrite)
