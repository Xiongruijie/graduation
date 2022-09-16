import math
import random

from matplotlib import pyplot as plt
from shapely.geometry import Polygon, Point


def distance(a, b):
    return Point(a).distance(Point(b))


class Square(object):
    def __init__(self):
        self.x0, self.y0 = random.random(), random.random()
        theta = random.randint(0, 90) * math.pi / 180  # Angle of rotation
        self.x1 = self.x0 + (0.1 * math.cos(theta))
        self.x2 = self.x1 + (0.1 * math.cos((90 * math.pi/180) + theta))
        self.x3 = self.x2 + (0.1 * math.cos((180 * math.pi/180) + theta))
        self.y1 = self.y0 + (0.1 * math.sin(theta))
        self.y2 = self.y1 + (0.1 * math.sin((90 * math.pi/180) + theta))
        self.y3 = self.y2 + (0.1 * math.sin((180 * math.pi/180) + theta))
        self.corners = ((self.x0, self.y0), (self.x1, self.y1),
                        (self.x2, self.y2), (self.x3, self.y3))

    @property
    def center(self):
        """(x, y) of the center of the polygon."""
        return Polygon(self.corners).centroid.coords[0]

    @property
    def half_diag(self):
        """The distance of 1/2 the shape's diagonal (center-to-corner)."""
        p0, p1, p2, p3 = self.corners
        return 0.5 * distance(p0, p1) * math.sqrt(2)


def test_overlap(square1, square2):
    """Do two shapes overlap?

    Note this is a 'conservative' test.  May return True if they do not
    (false positive), but will never return False if they do (false negative).
    """

    # Distance between two centers
    ctc = distance(square1.center, square2.center)
    # Sum of half-diagonals
    halfdiags = square1.half_diag + square2.half_diag
    res = ctc < halfdiags
    return res


class Squares(object):
    def __init__(self):
        self.squares = []

    def add_square(self):
        new_square = Square()
        if not self.squares:
            # Initial empty list/container - just add without any tests
            self.squares.append(new_square)
        else:
            while True:
            # Test that new_square overlaps with existing
                res = [test_overlap(square, new_square) for square in self.squares]
                if any(res):
                    # We have at least 1 case of overlap (1 True)
                    new_square = Square()
                else:
                    # Safe to add
                    self.squares.append(new_square)
                    break

    def plot_squares(self):
        for square in self.squares:
            (x0, y0), (x1, y1), (x2, y2), (x3, y3) = square.corners
            plt.plot([x0, x1, x2, x3, x0], [y0, y1, y2, y3, y0])