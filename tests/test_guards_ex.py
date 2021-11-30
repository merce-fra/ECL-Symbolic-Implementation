import pytest
import ppl as ppl
import pysymbrobustness.ppl_extension.linear_algebra as constraints


@pytest.fixture(autouse=True)
def empty_cs_polyhedra():
    empty_cs = ppl.Constraint_System()
    return ppl.C_Polyhedron(empty_cs)


@pytest.fixture(autouse=True)
def len_1_cs_polyhedra():
    x = ppl.Variable(0)
    y = ppl.Variable(1)
    cs = ppl.Constraint_System()
    cs.insert(x - y >= 0)
    cs.insert(x - y >= 3)
    return ppl.C_Polyhedron(cs)


@pytest.fixture(autouse=True)
def len_1_cs_again_polyhedra():
    x = ppl.Variable(0)
    y = ppl.Variable(1)
    cs_len_1_2 = [
        x - y >= 3
    ]
    return ppl.C_Polyhedron(constraints.cs_from_list(cs_len_1_2))


@pytest.fixture(autouse=True)
def len_2_cs_polyhedra():
    x = ppl.Variable(0)
    y = ppl.Variable(1)
    cs_2 = [
        x - y >= 8,
        x >= 0
    ]
    return ppl.C_Polyhedron(constraints.cs_from_list(cs_2))


@pytest.fixture(autouse=True)
def len_3_cs_polyhedra():
    x = ppl.Variable(0)
    y = ppl.Variable(1)
    z = ppl.Variable(2)
    cs_3 = [
        x >= 0,
        y >= 9,
        z >= 5
    ]
    return ppl.C_Polyhedron(constraints.cs_from_list(cs_3))
