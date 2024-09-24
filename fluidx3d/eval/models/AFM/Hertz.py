# -*- coding: utf-8 -*-
"""
File to house the  class.

Created on Tue Sep 24 12:43:24 2024

@author: Richard Kellnberger
"""


import numpy as np


class Hertz:

    def __init__(self, R1, R2, poissonRatio=0.48):
        R12 = 1 / (1 / R1 + 1 / R2)
        self.poissonRatio = poissonRatio
        self.prefactor = (3 / 4 * (1 - self.poissonRatio**2) / np.sqrt(R12)) ** (2 / 3)

    def deformation(self, force, youngsModulus, offset):
        return self.prefactor * (force / youngsModulus) ** (2 / 3) + offset

    def force(self, deformation, youngsModulus, offset):
        return youngsModulus * ((deformation - offset) / self.prefactor) ** (3 / 2)
