# -*- coding: utf-8 -*-
"""
File to house the  class.

Created on Tue Sep 24 13:06:25 2024

@author: Richard Kellnberger
"""


import numpy as np


class HertzConical:

    def __init__(self, coneHalfAngle, poissonRatio=0.48):
        self.poissonRatio = poissonRatio
        self.prefactor = (np.pi / 2 * (1 - self.poissonRatio**2) / np.tan(coneHalfAngle)) ** (1 / 2)

    def deformation(self, force, youngsModulus, offset):
        return self.prefactor * (force / youngsModulus) ** (1 / 2) + offset

    def force(self, deformation, youngsModulus, offset):
        return youngsModulus * ((deformation - offset) / self.prefactor) ** 2
