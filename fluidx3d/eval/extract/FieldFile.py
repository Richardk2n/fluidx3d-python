# -*- coding: utf-8 -*-
"""
File to house the   class.

Created on Mon Feb 17 16:32:48 2025

@author: Richard Kellnberger
"""
# from multiprocessing import Pipe, Process
# from multiprocessing.connection import Connection
from pathlib import Path
from typing import Self

import numpy as np
import numpy.ma as ma
import pyvista as pv

"""
class FieldFile(Process):

    def __init__(
        self,
        path: Path,
        name: str,
        numberComponents: int,
        conversionFactor: int,
        connection: Connection,
    ):
        super().__init__()
        self.path = path
        self.name = name
        self.numberComponents = numberComponents
        self.conversionFactor = conversionFactor
        self.connection = connection
        self.running = True

    def run(self) -> None:
        data = pv.read(self.path)
        self.dims = data.dimensions
        Lx, Ly, Lz = self.dims
        field = data.get_array(self.name) * self.conversionFactor
        self.field = np.reshape(field, (Lx, Ly, Lz, self.numberComponents), "F")
        while self.running:
            command: str = self.connection.recv()
            match command:
                case "getDimensions":
                    self.connection.send(self.dims)
                case _:
                    print(f"Unexpected command: {command}")

    def stop(self) -> None:
        self.running = False
"""


class FieldFile:

    def __init__(
        self,
        path: Path,
        name: str,
        numberComponents: int,
        conversionFactor: int,
    ):
        super().__init__()
        self.path = path
        self.name = name
        self.numberComponents = numberComponents
        self.conversionFactor = conversionFactor
        self.cached = False

    def cache(self) -> Self:
        if self.cached:
            return self
        data = pv.read(self.path)
        self.dims = data.dimensions
        Lx, Ly, Lz = self.dims
        field = data.get_array(self.name) * self.conversionFactor
        self.field = np.reshape(field, (Lx, Ly, Lz, self.numberComponents), "F")
        self.cached = True
        return self

    def getField(self):
        if not self.cached:
            self.cache()
        return self.field

    def mask(self, mask) -> Self:
        if not self.cached:
            self.cache()
        self.field = ma.masked_array(self.field, mask=mask)
        return self
