# coding=utf-8
import ppl as ppl
import pytest
import pysymbrobustness.ppl_extension.polyhedra as poly

x = ppl.Variable(0)
y = ppl.Variable(1)
z = ppl.Variable(2)


@pytest.fixture(autouse=True)
def polyhedron_closed_0():
    cs = ppl.Constraint_System()
    cs.insert(x <= 1)
    cs.insert(x >= 0)
    cs.insert(y <= 2)
    cs.insert(y >= 0)
    cs.insert(z <= x - y)
    cs.insert(z >= 0)

    return ppl.C_Polyhedron(cs)


@pytest.fixture(autouse=True)
def polyhedron_open_0_associated():
    cs = ppl.Constraint_System()
    cs.insert(x < 1)
    cs.insert(x > 0)
    cs.insert(y < 2)
    cs.insert(y > 0)
    cs.insert(z < x - y)
    cs.insert(z > 0)

    return ppl.NNC_Polyhedron(cs)


@pytest.fixture(autouse=True)
def polyhedron_with_equalities_0():
    cs = ppl.Constraint_System()
    cs.insert(x <= 1)
    cs.insert(y <= 2)
    cs.insert(y >= 0)
    cs.insert(z == x - y)
    cs.insert(z >= 0)

    return ppl.C_Polyhedron(cs)


@pytest.fixture(autouse=True)
def polyhedron_with_equalities_1():
    cs = ppl.Constraint_System()
    cs.insert(x == 1)
    cs.insert(y <= 2)
    cs.insert(y >= 0)
    cs.insert(z == x - y)
    cs.insert(z >= 0)

    return ppl.C_Polyhedron(cs)


@pytest.fixture(autouse=True)
def polyhedron_with_equalities_2():
    P = ppl.C_Polyhedron(5, 'universe')
    cs = ppl.Constraint_System()
    cs.insert(x == 1)
    cs.insert(y <= 2)
    cs.insert(y >= 0)
    cs.insert(z == x - y)
    cs.insert(z >= 0)
    P.add_constraints(cs)
    return P

@pytest.fixture(autouse=True)
def polyhedron_with_strict_inequalities():
    cs = ppl.Constraint_System()
    cs.insert(x > 1)
    cs.insert(y <= 2)
    cs.insert(y >= 0)
    cs.insert(z < x - y)
    cs.insert(z >= 0)
    return ppl.NNC_Polyhedron(cs)

@pytest.fixture(autouse=True)
def polyhedron_with_strict_inequalities_open():
    cs = ppl.Constraint_System()
    cs.insert(x > 1)
    cs.insert(y < 2)
    cs.insert(y > 0)
    cs.insert(z < x - y)
    cs.insert(z > 0)
    return ppl.NNC_Polyhedron(cs)

