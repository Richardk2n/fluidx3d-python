# -*- coding: utf-8 -*-
"""
File to house the  class.

Created on Mon Apr 22 15:01:49 2024

@author: Richard Kellnberger
"""

import numpy as np


class PoiseuilleFlow:
    numberCoords = 0

    def __init__(self, dynamicViscosity: float):
        self.dynamicViscosity = dynamicViscosity

    def proportionalityFactorVelocityProfile(self, pos):
        pass

    def proportionalityFactorVolumeFlow(self):
        pass

    # Here is where the planned quantity system will help clean up functions a lot
    def getCenterVelocityFromPressureGradiant(self, pressureGradient: float):
        return (
            self.proportionalityFactorVelocityProfile(np.zeros(self.numberCoords))
            * pressureGradient
        )

    def getCenterVelocityFromVolumeFlow(self, volumeFlow: float):
        return self.getCenterVelocityFromPressureGradiant(
            self.getPressureGradientFromVolumeFlow(volumeFlow)
        )

    def getPressureGradientFromCenterVelocity(self, centerVelocity: float):
        return centerVelocity / self.proportionalityFactorVelocityProfile(
            np.zeros(self.numberCoords)
        )

    def getPressureGradientFromVolumeFlow(self, volumeFlow: float):
        return volumeFlow / self.proportionalityFactorVolumeFlow()

    def getVolumeFlowFromPressureGradiant(self, pressureGradient: float):
        return self.proportionalityFactorVolumeFlow() * pressureGradient

    def getVolumeFlowFromCenterVelocity(self, centerVelocity: float):
        return self.getVolumeFlowFromPressureGradiant(
            self.getPressureGradientFromCenterVelocity(centerVelocity)
        )

    def getVelocityProfileFromPressureGradient(
        self, pressureGradient: float, pos
    ):  # pos is the radius for Planar and Circular and y, z (from center) for Rectangular, Square and Oval
        return self.proportionalityFactorVelocityProfile(pos) * pressureGradient

    def getVelocityProfileFromVolumeFlow(
        self, volumeFlow: float, pos
    ):  # pos is the radius for Planar and Circular and y, z (from center) for Rectangular, Square and Oval
        return self.getVelocityProfileFromPressureGradient(
            self.getPressureGradientFromVolumeFlow(volumeFlow), pos
        )

    def getVelocityProfileFromCenterVelocity(
        self, centerVelocity: float, pos
    ):  # pos is the radius for Planar and Circular and y, z (from center) for Rectangular, Square and Oval
        return self.getVelocityProfileFromPressureGradient(
            self.getPressureGradientFromCenterVelocity(centerVelocity), pos
        )


class PlanarPoiseuilleFlow(PoiseuilleFlow):
    numberCoords = 1

    def proportionalityFactorVelocityProfile(self, pos):
        return ((self.size / 2) ** 2 - pos[0] ** 2) / (2 * self.dynamicViscosity)

    def proportionalityFactorVolumeFlow(self):
        return self.size**3 / (12 * self.dynamicViscosity)

    def __init__(self, dynamicViscosity: float, size: float):
        super().__init__(dynamicViscosity)
        self.size = size


class RectangularPoiseuilleFlow(PoiseuilleFlow):
    numberCoords = 2

    def proportionalityFactorVelocityProfile(self, pos):
        sum_ = 0
        """
        We have to sum to infinity here. However this converges so fast, that 3 would already be
        practically infinit. As this does not get executed ofted, we might as well sum to 100.
        """
        for n in range(1, 101):
            beta = (2 * n - 1) * np.pi / self.sideA
            sum_ -= (
                (2 * n - 1) ** -3
                * np.cosh(beta * pos[1])
                / np.cosh(beta * self.sideB / 2)
                * (-1) ** (n % 2)
                * np.cos(beta * pos[0])
            )

        return ((self.sideA / 2) ** 2 - pos[0] ** 2 - 8 * self.sideA**2 * sum_ / np.pi**3) / (
            2 * self.dynamicViscosity
        )

    def proportionalityFactorVolumeFlow(self):
        sum_ = 0
        """
        We have to sum to infinity here. However this converges so fast, that 3 would already be
        practically infinit. As this does not get executed ofted, we might as well sum to 100.
        """
        for n in range(1, 101):
            beta = (2 * n - 1) * np.pi * self.sideB / self.sideA
            sum_ += (2 * n - 1) ** -5 * np.tanh(beta / 2)

        return (
            self.sideA**3
            * (self.sideB / 12 - 16 * self.sideA * sum_ / np.pi**5)
            / self.dynamicViscosity
        )

    def __init__(self, dynamicViscosity: float, sideA: float, sideB: float):
        super().__init__(dynamicViscosity)
        self.sideA = sideA
        self.sideB = sideB


class SquarePoiseuilleFlow(PoiseuilleFlow):
    numberCoords = 2

    def proportionalityFactorVelocityProfile(self, pos):
        sum_ = 0
        """
        We have to sum to infinity here. However this converges so fast, that 3 would already be
        practically infinit. As this does not get executed ofted, we might as well sum to 100.
        """
        for n in range(1, 101):
            beta = (2 * n - 1) * np.pi / self.side
            sum_ -= (
                (2 * n - 1) ** -3
                * np.cosh(beta * pos[1])
                / np.cosh(beta * self.side / 2)
                * (-1) ** (n % 2)
                * np.cos(beta * pos[0])
            )  # Can be transformed into a form clearly indifferent to a<->b switch, but is longer

        return ((self.side / 2) ** 2 - pos[0] ** 2 - 8 * self.side**2 * sum_ / np.pi**3) / (
            2 * self.dynamicViscosity
        )

    def proportionalityFactorVolumeFlow(self):
        sum_ = 0
        """
        We have to sum to infinity here. However this converges so fast, that 3 would already be
        practically infinit. As this does not get executed ofted, we might as well sum to 100.
        """
        for n in range(1, 101):
            beta = (2 * n - 1) * np.pi
            sum_ += (2 * n - 1) ** -5 * np.tanh(beta / 2)

        return self.side**4 * (1 / 12 - 16 * sum_ / np.pi**5) / self.dynamicViscosity

    def __init__(self, dynamicViscosity: float, side: float):
        super().__init__(dynamicViscosity)
        self.side = side


class OvalPoiseuilleFlow(PoiseuilleFlow):
    numberCoords = 2

    def proportionalityFactorVelocityProfile(self, pos):
        return (
            (self.semiAxisA * self.semiAxisB) ** 2
            / (2 * self.dynamicViscosity * (self.semiAxisA**2 + self.semiAxisB**2))
            * (1 - (pos[0] / self.semiAxisA) ** 2 - (pos[1] / self.semiAxisB) ** 2)
        )

    def proportionalityFactorVolumeFlow(self):
        return (
            np.pi
            * (self.semiAxisA * self.semiAxisB) ** 3
            / (4 * self.dynamicViscosity * (self.semiAxisA**2 + self.semiAxisB**2))
        )

    def __init__(self, dynamicViscosity: float, semiAxisA: float, semiAxisB: float):
        super().__init__(dynamicViscosity)
        self.semiAxisA = semiAxisA
        self.semiAxisB = semiAxisB


class CircularPoiseuilleFlow(PoiseuilleFlow):
    numberCoords = 1

    def proportionalityFactorVelocityProfile(self, pos):
        return (self.radius**2 - pos[0] ** 2) / (4 * self.dynamicViscosity)

    def proportionalityFactorVolumeFlow(self):
        return np.pi * self.radius**4 / (8 * self.dynamicViscosity)

    def __init__(self, dynamicViscosity: float, radius: float):
        super().__init__(dynamicViscosity)
        self.radius = radius


# Only sub classes are for public consumption
__all__ = [
    "PlanarPoiseuilleFlow",
    "RectangularPoiseuilleFlow",
    "SquarePoiseuilleFlow",
    "OvalPoiseuilleFlow",
    "CircularPoiseuilleFlow",
]
