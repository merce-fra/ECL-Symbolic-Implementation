from pathlib import Path
from typing import List, Tuple
from math import atan2, gcd
from dataclasses import dataclass

import ppl
from pysymbrobustness.ppl_extension.piecewise_linear_function import Spline, SubSpline
from pysymbrobustness.ppl_extension.rational_linear_algebra import RationalLinearExpression


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

        return Point(self.x * (other.lcm // new_gcd) + other.x * (self.lcm // new_gcd),
                     self.y * (other.lcm // new_gcd) + other.y * (self.lcm // new_gcd),
                     new_lcm)

    def __neg__(self):
        return Point(-self.x, -self.y, self.lcm)

    def __sub__(self, other):
        return self + (-other)

    def __floordiv__(self, other):
        return Point(self.x, self.y, self.lcm * other)

    def __str__(self):
        return f"({self.x} / {self.lcm}, {self.y} / {self.lcm})"

    def __repr__(self):
        return f"Point({self.x}, {self.y}, {self.lcm})"


def visual(perm: Spline) -> List[Tuple[RationalLinearExpression, ppl.constraint.Constraint_System]]:
    return list(map(lambda p: (p.function, p.polyhedron.generators()), perm.sub_splines))


COLORS = ['blue!50', 'green!50',
          'red!50', 'yellow!50!black', 'violet',
          'blue!75', 'green!75', 'red!75',
          'teal', 'lime!50!black', 'olive',
          'cyan']


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


class LatexSplineVisualisation(object):
    def __init__(self, output_file: Path):
        self.output_file = output_file.absolute()
        self.latex = ""
        self.colors = []
        self.index_color = 0

        self.latex_intro()

    def latex_intro(self):
        self.latex += "\\documentclass{article}\n" \
                      "\\usepackage{tikz}\n" \
                      "\\usepackage{caption}\n" \
                      "\n" \
                      "\\begin{document}\n"

    def latex_new_tikz(self):
        self.latex += "\\begin{figure}\n" \
                      "\\begin{tikzpicture}[scale=4]\n" \
                      "\\draw[black,->] (0,0) -- (2,0);\n" \
                      "\\draw[black,->] (0,0) -- (0,2);\n"

    def latex_tikz_end(self, caption=""):
        self.latex += "\\end{tikzpicture}\n" \
                      f"\\caption{{ {caption} }}\n" \
                      "\\end{figure}\n"

    def latex_outro(self):
        self.latex += "\\end{document}\n"

    def add_spline(self, spline: Spline, caption: str = ""):
        used_colors = {}
        tikz = ""
        color_map = {}

        self.latex_new_tikz()

        # Assign a color to each sub_spline
        # Adding a new color if there is none already existing

        for sub_n, sub in enumerate(spline.sub_splines):
            color_of_sub = None
            for color, sub_color in self.colors:
                if sub_color.function.is_equal_to(sub.function):
                    color_of_sub = color
                    break

            if color_of_sub is None:
                # No color have been found
                color_of_sub = COLORS[self.index_color]
                self.colors.append((color_of_sub, sub))
                self.index_color += 1

            color_map[sub_n] = color_of_sub
            if color_of_sub not in used_colors:
                used_colors[color_of_sub] = sub

        # Generate the tikz polygone and the index for each subspline
        for sub_n, sub in enumerate(spline.sub_splines):
            color = color_map[sub_n]

            tikz += f"\\draw[fill={color}, opacity=0.25] "
            points = transform_points_order(generator_to_points(sub))
            for p in points:
                x = f"{p.x} / {p.lcm}"
                y = f"{p.y} / {p.lcm}"
                tikz += f"({x}, {y}) -- "

            tikz += "cycle;\n"

            # Isobarycenter computation
            isobar = sum(points, start=Point(0, 0, 1)) // len(points)
            tikz += f"\\node at {isobar} {{ \\tiny{{ {sub_n} }} }};\n"

        # Add the legend
        for i, color in enumerate(used_colors):
            sub = used_colors[color]
            function = sub.function
            function_repr = "{le} / {lcm}".format(
                le=str(function.linear_expression),
                lcm=function.lcm
            )
            tikz += f"\\node[{color}] at (2.5,{i} / 2) {{ {function_repr} }};\n"

        self.latex += tikz
        self.latex_tikz_end(caption)

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
