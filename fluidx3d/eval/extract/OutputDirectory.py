# -*- coding: utf-8 -*-
"""
File to house the   class.

Created on Mon Feb 17 17:53:07 2025

@author: Richard Kellnberger
"""

import json
import os
from pathlib import Path

import numpy as np
from natsort import natsorted

from fluidx3d.eval.extract.FileArray import FileArray

TYPE_W = 0b00000001  # Simple bounce-back wall
TYPE_VW = 0b00000010  # Velocity boundary wall
TYPE_PW = 0b00000100  # Pressure boundary wall
TYPE_MW = 0b00001000  # Massflow boundary wall (velocity boundary with fixed mass flow)
TYPE_OW = 0b00010000  # Outflow boundary wall (velocity boundary with arbitrary mass flow)
TYPE_SW = 0b00100000  # Free slip boundary wall

TYPE_F = 0b10000000  # fluid // TODO remove // ATTENTION also recycled for force evaluation

symbols = {
    "cell": "cells",
    "inout": "inout",
    "velocity": "u",
    "density": "ρ",
    "flags": "flags",
    "strainRateTensor": "D",
    "polymerConformationTensor": "C",
}
conversion = {
    "cell": [0, 0, 0],
    "inout": [0, 0, 0],
    "velocity": [0, 1, -1],
    "density": [1, -3, 0],
    "flags": [0, 0, 0],
    "strainRateTensor": [0, 0, -1],
    "polymerConformationTensor": [0, 0, 0],
}
numberComponents = {
    "cell": 0,
    "inout": 1,
    "velocity": 3,
    "density": 1,
    "flags": 1,
    "strainRateTensor": 6,
    "polymerConformationTensor": 6,
}


class OutputDirectory:

    def __init__(self, path: Path):
        self.path = path
        self.scout()
        self.masks: dict[int, np.ndarray] = {}

    def scout(self) -> None:
        jsonFile = self.path / "parameters.json"
        f = open(jsonFile)
        parameters = json.load(f)
        f.close()

        if "internal" in parameters:  # TODO proper versions
            internal = parameters["internal"]
            if "SI" in internal:
                self.SI = internal["SI"]
            else:
                print("Missing SI in parameters.json")
        else:
            print("Output directory too old: Missing internal in parameters.json")
            if "info" in parameters and "conversions" in parameters["info"]:
                print("Falling back to conversions")
                L0 = parameters["info"]["conversions"]["L0"]
                T0 = parameters["info"]["conversions"]["T0"]
                rho0 = parameters["info"]["conversions"]["rho0"]

                self.SI = {"kg": rho0 * L0**3, "m": L0, "s": T0}

        vtkPath = self.path / "vtkfiles"
        arrays = [f for f in natsorted(os.listdir(vtkPath)) if (vtkPath / f).is_dir()]

        self.numberArrays = len(arrays)

        self.arrays: list[FileArray] = []
        self.ids: dict[str, int] = {}

        for i in range(self.numberArrays):
            name = arrays[i]
            self.ids[name] = i

            kg, m, s = conversion[name]
            fac = self.SI["kg"] ** kg * self.SI["m"] ** m * self.SI["s"] ** s
            if name == "density":
                fac /= self.SI["ReScale"]
            self.arrays.append(
                FileArray(vtkPath / name, symbols[name], fac, numberComponents[name])
            )
            self.arrays[-1].scout()

    def cache(self):
        for array in self.arrays:
            array.cache()

    def mask(self):
        for name in self.ids:
            nc = numberComponents[name]
            mask = self.getMask(nc)
            self.arrays[self.ids[name]].mask(mask)

    def getMask(self, numberComponents: int):
        if numberComponents in self.masks:
            return self.masks[numberComponents]
        elif "flags" in self.ids:
            if 1 not in self.masks:
                self.masks[1] = self.flags[-1].getField().astype(int) & TYPE_W
                self.Lx, self.Ly, self.Lz, _ = self.masks[1].shape
            if numberComponents == 1:
                return self.masks[1]

            self.masks[numberComponents] = (
                np.ones((self.Lx, self.Ly, self.Lz, numberComponents)) * self.masks[1]
            )
            return self.masks[numberComponents]
        else:
            print("Cannot mask without flags")

    def getDims(self):
        try:
            if self.Lx:
                return self.Lx, self.Ly, self.Lz
        except Exception:
            pass

        for name in self.ids:
            if name != "cell":
                self.Lx, self.Ly, self.Lz, _ = self.arrays[self.ids[name]][0].getField().shape
                return self.Lx, self.Ly, self.Lz
        print("Cannot get dimensions of non existent fields")

    def getExtent(self, fac: float):
        Lx, Ly, Lz = self.getDims()
        m = self.SI["m"]
        Ex = (-0.5 * m * fac, (Lx - 1 + 0.5) * m * fac)
        Ey = (-Ly / 2 * m * fac, Ly / 2 * m * fac)
        Ez = (-Lz / 2 * m * fac, Lz / 2 * m * fac)
        return Ex, Ey, Ez

    def __getattr__(self, name: str):
        if name in symbols:
            return self.arrays[self.ids[name]]
        else:
            raise Exception(f"VTK Array {name} does not exist!")
