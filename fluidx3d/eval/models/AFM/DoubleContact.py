# -*- coding: utf-8 -*-
"""
File to house the  class.

Created on Tue Sep 24 13:52:03 2024

@author: Richard Kellnberger
"""


import numpy as np

from fluidx3d.eval.models.AFM.Hertz import Hertz


class DoubleContact:

    def __init__(self, R1, R2, R3=np.inf, poissonRatio=0.48):
        R12 = 1 / (1 / R1 + 1 / R2)
        R32 = 1 / (1 / R3 + 1 / R2)
        self.k = R32 ** (1 / 3) / (R32 ** (1 / 3) + R12 ** (1 / 3))
        self.hertz = Hertz(R1, R2, poissonRatio)

    def deformation(self, force, youngsModulus, offset):
        return self.hertz.deformation(force, youngsModulus, offset) / self.k

    def force(self, deformation, youngsModulus, offset):
        return self.hertz.force(deformation, youngsModulus, offset) * self.k ** (3 / 2)
