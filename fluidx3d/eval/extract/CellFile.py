# -*- coding: utf-8 -*-
"""
File to house the  class.

Created on Thu Jan 30 13:21:13 2025

@author: Richard Kellnberger
"""

from pathlib import Path
from typing import Self

import numpy as np
import pyvista as pv


class CellFile:

    def __init__(self, path: Path):
        self.path = path
        self.cached = False

    def cache(self) -> Self:
        if self.cached:
            return self
        data = pv.read(self.path)
        self.points = data.points

        self.com = np.average(self.points, 0)
        self.p1 = self.points[1] - self.com
        self.p5 = self.points[5] - self.com

        self.cached = True
        return self

    def getSpecialPoints(self):
        if not self.cached:
            self.cache()
        return self.p1, self.p5

    @staticmethod
    def getSpecialParallel(path: Path):
        data = pv.read(path)
        points = data.points

        com = np.average(points, 0)
        p1 = points[1] - com
        p5 = points[5] - com
        return p1, p5
