import pytest
import ppl as ppl

x = ppl.Variable(0)
y = ppl.Variable(1)
z = ppl.Variable(2)


# Examples of list of constraints
@pytest.fixture(autouse=True)
def empty_cs_list():
    return []


@pytest.fixture(autouse=True)
def empty_linear_expression():
    return ppl.Linear_Expression([], 0)


@pytest.fixture(autouse=True)
def empty_constraint_system():
    return ppl.Constraint_System()


@pytest.fixture(autouse=True)
def cs_list_0():
    return [x - y >= 0, 3 + y >= 1, 4 + x + z <= 2]


@pytest.fixture(autouse=True)
def linear_expressions_0():
    return [x - y, 2 + y, -2 -x - z]


@pytest.fixture(autouse=True)
def constraint_system_0():
    cs = ppl.Constraint_System()
    cs.insert(x - y >= 0)
    cs.insert(3 + y >= 1)
    cs.insert(4 + x + z <= 2)
    return cs


@pytest.fixture(autouse=True)
def cs_list_1():
    return [x >= 0]


@pytest.fixture(autouse=True)
def linear_expressions_1():
    return [x]


@pytest.fixture(autouse=True)
def constraint_system_1():
    cs = ppl.Constraint_System()
    cs.insert(x >= 0)
    return cs


# List of linear expressions

@pytest.fixture(autouse=True)
def linear_expression_0():
    return x - y


@pytest.fixture(autouse=True)
def linear_expression_0_mirror():
    return -x + y


@pytest.fixture(autouse=True)
def linear_expression_1():
    return 3 * x - y + 5


@pytest.fixture(autouse=True)
def linear_expression_1_mirror():
    return -3 * x + y - 5
