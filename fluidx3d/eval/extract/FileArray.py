# -*- coding: utf-8 -*-
"""
File to house the  class.

Created on Thu Jan 30 13:26:48 2025

@author: Richard Kellnberger
"""

import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pyvista as pv
from natsort import natsorted

from fluidx3d.eval.extract.CellFile import CellFile

symbols = {"cell": "cells"}


class FileArray:

    def __init__(self, path: Path):
        self.path = path
        self.symbol = symbols[self.path.name]
        self.classType = CellFile if self.symbol == "cells" else None
        files = natsorted(os.listdir(path))
        self.vtkIntervall = int(
            files[1].split("_")[1].split(".")[0]
        )  # TODO this could be determined further up, but is fast anyway
        self.numberFiles = len(files)
        self.fileCache = [None] * self.numberFiles
        self.fullyCached = False

    def cacheAll(self):
        paths = [
            self.path / f"{self.symbol}_{n*self.vtkIntervall}.vtk" for n in range(self.numberFiles)
        ]

        with Pool() as p:
            self.fileCache = p.map(self.classType, paths)
        self.fullyCached = True

    def getFile(self, i: int):
        if not self.fileCache[i]:
            self.fileCache[i] = self.classType(
                self.path / f"{self.symbol}_{i*self.vtkIntervall}.vtk"
            )  # type: ignore
            self.fullyCached = None not in self.fileCache
        return self.fileCache[i]

    def getSpecialParallelWrapper(self, n):
        return CellFile.getSpecialParallel(self.path / f"{self.symbol}_{n*self.vtkIntervall}.vtk")

    def getSpecialPoints(self):  # TODO make sensible
        R = self.getFile(1).p1[0]  # TODO make sensible

        with Pool() as p:
            p1s, p5s = np.swapaxes(
                p.map(self.getSpecialParallelWrapper, range(self.numberFiles)), 0, 1
            )

        return p1s / R, p5s / R
