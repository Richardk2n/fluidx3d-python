# -*- coding: utf-8 -*-
"""
File to house the  class.

Created on Thu Jan 30 13:26:48 2025

@author: Richard Kellnberger
"""

import os
from multiprocessing import Pool
from pathlib import Path

from natsort import natsorted

from fluidx3d.eval.extract.CellFile import CellFile
from fluidx3d.eval.extract.FieldFile import FieldFile


class FileArray:

    def __init__(
        self,
        path: Path,
        symbol: str,
        conversionFactor: int,
        numberComponents: int,
    ):
        self.path = path
        self.symbol = symbol
        self.conversionFactor = conversionFactor
        self.numberComponents = numberComponents

    def scout(self) -> None:
        files = natsorted(os.listdir(self.path))
        self.numberFiles = len(files)

        if self.numberFiles > 1:
            self.vtkIntervall = int(
                files[1].split("_")[1].split(".")[0]
            )  # this could be determined further up, but is fast anyway
        else:
            self.vtkIntervall = 0

        self.files: list[CellFile | FieldFile] = []
        for i in range(self.numberFiles):
            if self.symbol == "cells":
                self.files.append(CellFile(self.path / f"{self.symbol}_{i*self.vtkIntervall}.vtk"))
            else:
                self.files.append(
                    FieldFile(
                        self.path / f"{self.symbol}_{i*self.vtkIntervall}.vtk",
                        self.path.name,
                        self.numberComponents,
                        self.conversionFactor,
                    )
                )

    def cache(self) -> None:
        with Pool() as p:
            if self.symbol == "cells":
                self.files = p.map(CellFile.cache, self.files)
            else:
                self.files = p.map(FieldFile.cache, self.files)

    def mask(self, mask) -> None:
        if self.path.name == "flags" or self.symbol == "cells":
            return
        with Pool() as p:
            self.files = p.starmap(FieldFile.mask, zip(self.files, [mask] * self.numberFiles))

    def __getitem__(self, i) -> CellFile | FieldFile:
        return self.files[i]


"""
        with Pool() as p:
            self.fileCache = p.map(self.classType, paths)

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
"""
