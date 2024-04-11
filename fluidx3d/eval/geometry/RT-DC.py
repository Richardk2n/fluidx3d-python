# -*- coding: utf-8 -*-
"""
File to house the RT-DC class.

Created on Thu Apr 11 10:20:20 2024

@author: Richard Kellnberger
"""

from typing import Dict, List, Tuple

import numpy as np


class RT_DC:

    def __init__(self, widthSI: float, lengthSI: float, paddingSI: float):
        """
        Init the RT-DC class.

        Parameters
        ----------
        width : int
            The width of the small part of the channel.
        length : int
            The length of the small part of the channel.
        padding : int
            The length of one side of the straight large channel.

        Returns
        -------
        None.

        """
        self.masks: Dict[int, np.ndarray] = {}
        self.widthSI: float = widthSI
        self.lengthSI: float = lengthSI
        self.paddingSI: float = paddingSI

    def generateMasks(self, numberComponents: List[int], dims: Tuple[int, int, int], L0: float = 1):
        Lx, Ly, Lz = dims
        width: int = np.round(self.widthSI / L0)
        length: int = np.round(self.lengthSI / L0)
        padding: int = np.round(self.paddingSI / L0)
        mask = np.zeros(dims)

        for x in range(Lx):
            for y in range(Ly):
                for z in range(Lz):
                    if x < padding:
                        if y == 0 or y == Ly - 1 or z == 0 or z == Lz - 1:
                            mask[x, y, z] = 1
                    elif x < padding + width:
                        if (
                            y < x - (padding - 2)
                            or y >= Ly - x + (padding - 2)
                            or z == 0
                            or z == Lz - 1
                        ):
                            mask[x, y, z] = 1
                    elif x < padding + width + length:
                        if y < width + 1 or y >= 2 * width + 1 or z == 0 or z == Lz - 1:
                            mask[x, y, z] = 1
                    elif x < padding + 2 * width + length:
                        if (
                            y < (padding + width + length) - x + (width + 1)
                            or y >= x - (padding + width + length) + (2 * width + 1)
                            or z == 0
                            or z == Lz - 1
                        ):
                            mask[x, y, z] = 1
                    else:
                        if y == 0 or y == Ly - 1 or z == 0 or z == Lz - 1:
                            mask[x, y, z] = 1

        expandedMask = np.expand_dims(mask, -1)
        for c in numberComponents:
            h = np.ones((Lx, Ly, Lz, c))
            self.masks[c] = expandedMask * h

    def getMask(self, numberComponents: int) -> np.ndarray:
        return self.masks[numberComponents]


guck = RT_DC(20e-6, 300e-6, 400e-6)
steffen = RT_DC(40e-6, 300e-6, 400e-6)
