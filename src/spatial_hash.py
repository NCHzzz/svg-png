"""
Spatial hash grid for fast nearest-point and range queries on 2D coordinates.

Used throughout the pipeline: dedupe (two-track removal), fragment pruning,
T-connection detection, and junction endpoint registration.
"""

import math
from collections import defaultdict

_CELL = 8.0


class PointGrid:
    """Spatial hash for nearest-point and predicate queries."""

    def __init__(self):
        self.grid = defaultdict(list)

    def add(self, x, y, tag=None):
        self.grid[(int(x // _CELL), int(y // _CELL))].append((x, y, tag))

    def near(self, x, y, radius, exclude_tag=None):
        """Return (dist, x, y, tag) of nearest point within radius, or None."""
        gx, gy = int(x // _CELL), int(y // _CELL)
        reach = int(radius // _CELL) + 1
        best = None
        for dgy in range(-reach, reach + 1):
            for dgx in range(-reach, reach + 1):
                for (ox, oy, tag) in self.grid.get((gx + dgx, gy + dgy), ()):
                    if exclude_tag is not None and tag == exclude_tag:
                        continue
                    d = math.hypot(ox - x, oy - y)
                    if d <= radius and (best is None or d < best[0]):
                        best = (d, ox, oy, tag)
        return best

    def any_match(self, x, y, radius, pred):
        """True if any point within radius satisfies pred(tag)."""
        gx, gy = int(x // _CELL), int(y // _CELL)
        reach = int(radius // _CELL) + 1
        rr = radius * radius
        for dgy in range(-reach, reach + 1):
            for dgx in range(-reach, reach + 1):
                for (ox, oy, tag) in self.grid.get((gx + dgx, gy + dgy), ()):
                    if (ox - x) ** 2 + (oy - y) ** 2 <= rr and pred(tag):
                        return True
        return False
