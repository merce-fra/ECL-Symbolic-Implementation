import math
from functools import reduce
from pathlib import Path
from typing import List, Tuple, Dict
from math import atan2, gcd
from dataclasses import dataclass
import re

import ppl

from pysymbrobustness.ppl_extension.piecewise_linear_function import Spline, \
    SubSpline
from pysymbrobustness.ppl_extension.rational_linear_algebra import \
    RationalLinearExpression


@dataclass
class Point(object):
    x: int
    y: int
    lcm: int

    @property
    def rx(self):
        return self.x / self.lcm

    @property
    def ry(self):
        return self.y / self.lcm

    def __add__(self, other):
        new_gcd = gcd(self.lcm, other.lcm)
        new_lcm = self.lcm * other.lcm // new_gcd

        return Point(
            self.x * (other.lcm // new_gcd) + other.x * (self.lcm // new_gcd),
            self.y * (other.lcm // new_gcd) + other.y * (self.lcm // new_gcd),
            new_lcm)

    def __neg__(self):
        return Point(-self.x, -self.y, self.lcm)

    def __sub__(self, other):
        return self + (-other)

    def __floordiv__(self, other):
        return Point(self.x, self.y, self.lcm * other)

    def __str__(self):
        return f"({2 * self.x} / {self.lcm}, {2 * self.y} / {self.lcm})"

    def __repr__(self):
        return f"Point({self.x}, {self.y}, {self.lcm})"


def visual(perm: Spline) -> List[
    Tuple[RationalLinearExpression, ppl.constraint.Constraint_System]]:
    return list(map(lambda p: (p.function, p.polyhedron.generators()),
                    perm.sub_splines))


COLORS = ['blue!50', 'green!50',
          'yellow!50!black', 'violet',
          'blue!75', 'green!75',
          'teal', 'lime!50!black', 'olive',
          'cyan', 'white', 'gray']
COLORS += [f'red!{j}' for j in range(1, 101, 1)]
COLORS += [f'blue!{j}' for j in range(1, 101, 1) if j != 50]
COLORS += [f'green!{j}' for j in range(1, 101, 1) if j != 50 and j != 75]


def transform_points_order(coords: List[Point]) -> List[Point]:
    p1 = coords[0]
    p2 = coords[1]
    middle_coord = (p1 + p2) // 2

    new_coords = list(map(lambda c: c - middle_coord, coords))
    angles = list(enumerate(map(lambda c: atan2(c.rx, c.ry), new_coords)))
    sorted_angles = sorted(angles, key=lambda t: t[1])

    return [coords[i] for i, _ in sorted_angles]


def generator_to_points(sub: SubSpline) -> List[Point]:
    point_list = []

    for point in sub.polyhedron.generators():
        lcm = int(point.divisor().digits())
        x = int(point.coefficients()[0].digits())
        y = int(point.coefficients()[1].digits())

        point_list.append(Point(x, y, lcm))

    return point_list


@dataclass
class SubSplineTikzData(object):
    n: int
    color: str
    sub_spline: SubSpline

    @staticmethod
    def func_repr(sub_spline: SubSpline) -> Tuple[int, int, int, int]:
        f = sub_spline.function
        return (f.integer_coefficient(ppl.Variable(0)),
                f.integer_coefficient(ppl.Variable(1)),
                f.integer_inhomogeneous_term(),
                f.lcm)


class LatexSplineVisualisation(object):
    def __init__(self, output_file: Path):
        self.output_file: Path = output_file.absolute()
        self.latex: str = ""
        self.colors: Dict[Tuple[int, int, int, int], str] = {}
        self.index_color: int = 0

        self.latex_intro()

    def latex_intro(self):
        self.latex += "\\documentclass{article}\n" \
                      "\\usepackage{amsmath}\n" \
                      "\\usepackage{tikz}\n" \
                      "\\usepackage{array}\n" \
                      "\\usepackage{caption}\n" \
                      "\n" \
                      "\\begin{document}\n"

    def latex_new_tikz(self, xsize: int = 2, ysize: int = 2):
        xbar = "0/0"
        for x in range(1, xsize + 1):
            xbar += f",{2 * x}/{x}"

        ybar = "0/0"
        for y in range(1, ysize + 1):
            ybar += f",{2 * y}/{y}"
        self.latex += "\\begin{figure}\n" \
                      "\\begin{tikzpicture}\n" \
                      "\\tikzstyle{reg}=[minimum width=2cm, minimum height=2cm, rectangle]\n" \
                      "\n" \
                      "% Draw the scale\n" \
                      f"\\draw[black,->] (-0.2,0) -- ({2 * xsize + 0.25},0) node[right] {{$ x $}};\n" \
                      f"\\draw[black,->] (0,-0.25) -- (0,{2 * ysize + 0.25}) node[above] {{$ y $}};\n" \
                      f"\\draw[loosely dotted] (0,0) grid ({2 * xsize},{2 * ysize});\n" \
                      f"\\foreach \\x/\\xtext in {{ {xbar} }}\n" \
                      f"\\draw[shift={{(\\x,0)}}] (0pt,4pt) -- (0pt,-4pt) " \
                      f"node[below] {{$\\xtext$}};\n" \
                      f"\\foreach \\y/\\ytext in {{ {ybar} }}\n" \
                      f"\\draw[shift={{(0,\\y)}}] (4pt,0pt) -- (-4pt," \
                      f"0pt) node[left] {{$\\ytext$}};\n"

    def latex_tikz_end(self, caption: str = ""):
        self.latex += "\\end{tikzpicture}\n" + \
                      (
                          f"\\caption{{ {caption} }}\n" if caption != "" else "\n") + \
                      "\\end{figure}\n"

    def latex_outro(self):
        self.latex += "\\end{document}\n"

    def legend_intro(self):
        self.latex += "$$\n" \
                      "\\begin{array}{|l||l|}\n" \
                      "\\hline\n" \
                      "\\text{Permissiveness function} & \\text{" \
                      "Associated cells} \\\\\n" \
                      "\\hline\n"

    def legend_add_element(self, data: List[SubSplineTikzData]):
        x = ppl.Variable(0)
        y = ppl.Variable
        f_repr = data[0].sub_spline.function
        # f_repr = SubSplineTikzData.func_repr(data[0].sub_spline)
        le_repr = str(f_repr)
        le_repr = le_repr.replace('linear expression:', ' ')
        le_repr = le_repr.replace('x0', 'x')
        le_repr = le_repr.replace('x1', 'y')
        le_repr = re.sub(' lcm = [0-9]*', '', le_repr)
        """
        if f_repr[0] == 1:
            le_repr += "x +"
        elif f_repr[0] == -1:
            le_repr += "-x +"
        elif f_repr[0] != 0:
            le_repr += str(f_repr[0]) + "\cdot x +"

        if f_repr[1] == 1:
            le_repr += "y +"
        elif f_repr[1] == -1:
            le_repr += "-y +"
        elif f_repr[1] !=0 :
            le_repr += str(f_repr[1]) + "\cdot y + "

        if str(f_repr[2]) != 0:
            le_repr += str(f_repr[2])
        if le_repr == str():
            le_repr += str("0")
        """
        function_repr = "\\left( {le} \\right) / {lcm}".format(
            le=le_repr,
            lcm=f_repr.lcm  # f_repr[3]
        )
        list_num = str(data[0].n)
        for d in data[1:]:
            list_num += f", {d.n}"

        self.latex += f"{function_repr} & {list_num} \\\\\n" \
                      "\\hline\n"

    def legend_outro(self):
        self.latex += "\\end{array}\n" \
                      "$$\n"

    def add_spline(self, spline: Spline, caption: str = ""):
        tikz = ""
        max_x = 1  # The maximum x coordonate
        max_y = 1  # The maximum y coordonate
        function_map: Dict[
            Tuple[int, int, int, int], List[SubSplineTikzData]] = {}
        tikz_data: List[SubSplineTikzData] = []

        # First extract useful information (colors, existing functions, max_x and max_y)

        for sub in spline.sub_splines:
            f_repr = SubSplineTikzData.func_repr(sub)
            color = None

            if f_repr in self.colors:
                color = self.colors[f_repr]
            else:
                color = COLORS[self.index_color]
                self.index_color += 1
                self.colors[f_repr] = color

            data = SubSplineTikzData(
                sub_spline=sub,
                color=color,
                n=0
            )

            tikz_data.append(data)

            if f_repr not in function_map:
                function_map[f_repr] = []

            function_map[f_repr].append(data)

            for point in generator_to_points(sub):
                max_x = max(max_x, point.rx)
                max_y = max(max_y, point.ry)

        max_x = int(math.ceil(max_x))
        max_y = int(math.ceil(max_y))

        tikz_data.sort(key=lambda d: SubSplineTikzData.func_repr(d.sub_spline))

        self.latex_new_tikz(max_x, max_y)

        # Generate the tikz polygone and the index for each subspline
        for n, sub_data in enumerate(tikz_data):
            sub_data.n = n
            color = sub_data.color

            tikz += f"\\draw[fill={color}, opacity=0.5] "
            points = transform_points_order(
                generator_to_points(sub_data.sub_spline))
            for p in points:
                tikz += f"{p} -- "

            tikz += "cycle;\n"

            # Isobarycenter computation
            isobar = sum(points, start=Point(0, 0, 1)) // len(points)
            tikz += f"\\node at {isobar} {{  {n}  }};\n"

        self.latex += tikz
        self.latex_tikz_end(caption)

        # Add the legend
        self.legend_intro()

        for datas in function_map.values():
            self.legend_add_element(datas)

        self.legend_outro()

    def output(self):
        self.latex_outro()
        with open(self.output_file, 'w') as f:
            f.write(self.latex)


def latex_output(spline, file):
    visualisation = LatexSplineVisualisation(file)
    visualisation.add_spline(spline, "test")
    visualisation.output()


if __name__ == "__main__":
    pass
