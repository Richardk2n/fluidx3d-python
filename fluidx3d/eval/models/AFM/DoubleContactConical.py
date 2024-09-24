# -*- coding: utf-8 -*-
"""
File to house the  class.

Created on Tue Sep 24 14:21:58 2024

@author: Richard Kellnberger
"""


import numpy as np

from fluidx3d.eval.models.AFM.Hertz import Hertz
from fluidx3d.eval.models.AFM.HertzConical import HertzConical


class DoubleContactConical:

    def __init__(self, coneHalfAngle, R2, R3=np.inf, poissonRatio=0.48):
        self.hertzConical = HertzConical(coneHalfAngle, poissonRatio)
        self.hertz = Hertz(R3, R2, poissonRatio)

    def deformation(self, force, youngsModulus, offset):
        return self.hertzConical.deformation(force, youngsModulus, offset) + self.hertz.deformation(
            force, youngsModulus, offset
        )
