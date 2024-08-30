# -*- coding: utf-8 -*-
"""
File to house the   class.

Created on Thu Aug 29 15:42:48 2024

@author: Richard Kellnberger
"""


import matplotlib.pyplot as plt
import numpy as np

from fluidx3d.eval.style import cm

channel = 0
pipe = 1


class PowerLaw:
    def __init__(self, eta_p, lambda_p, eta_s, power):
        self.eta_p = eta_p
        self.gd0 = 1/lambda_p
        self.eta_s = eta_s
        self.power = power

    def prepareVelocityProfile(self, R, G, j):  # TODO check j
        if j == 0:
            raise Exception("Not implemented")

        def gd(r):
            e = 1 / (1 - self.power)

            return -((G * r / (2 * self.eta_p * self.gd0**self.power)) ** e)

        self.gd = gd

        def u(r):
            e = (2 - self.power) / (1 - self.power)

            return (
                1
                / e
                * (G / (2 * self.eta_p * self.gd0**self.power)) ** (1 / (1 - self.power))
                * (R**e - r**e)
            )

        self.u = u

        def Q():
            e = (4 - 3 * self.power) / (1 - self.power)

            return (
                np.pi
                * (G / (2 * self.eta_p * self.gd0**self.power)) ** (1 / (1 - self.power))
                / e
                * R**e
            )

        self.Q = Q()

    def eta(self, gd):
        gd = np.abs(gd)

        return np.clip(
            self.eta_p / (gd / self.gd0) ** self.power + self.eta_s, 0, self.eta_p + self.eta_s
        )

    def u(self, _):
        raise Exception("Forgot to call prepareVelocityProfile before u")

    def gd(self, _):
        raise Exception("Forgot to call prepareVelocityProfile before gd")

    def plot(self, r):
        u = self.u(r)

        plt.figure(figsize=(15.5 * cm, 15.5 / 2 * cm))
        plt.plot(r * 1e6, u, label=r"$v_x$")
        plt.xlabel(r"$r/\unit{\micro\meter}$")
        plt.ylabel(r"$v_x/\unit{\meter\per\second}$")
        plt.legend()
        plt.show()

        gd = self.gd(r)

        plt.figure(figsize=(15.5 * cm, 15.5 / 2 * cm))
        plt.plot(r * 1e6, gd, label=r"$\dot{\gamma}$")
        plt.xlabel(r"$r/\unit{\micro\meter}$")
        plt.ylabel(r"$\dot{\gamma}/\unit{\per\second}$")
        plt.legend()
        plt.show()

        eta = self.eta(gd)
        plt.figure(figsize=(15.5 * cm, 15.5 / 2 * cm))
        plt.plot(r * 1e6, eta, label=r"$\eta$")
        plt.xlabel(r"$r/\unit{\micro\meter}$")
        plt.ylabel(r"$\eta/\unit{\pascal\second}$")
        plt.legend()
        plt.show()

        rho = 1e3
        Re = rho * u * np.max(r) / eta

        plt.figure(figsize=(15.5 * cm, 15.5 / 2 * cm))
        plt.plot(r * 1e6, Re, label=r"$Re$")
        plt.xlabel(r"$r/\unit{\micro\meter}$")
        plt.ylabel(r"$Re$")
        plt.legend()
        plt.show()

        plt.figure(figsize=(15.5 * cm, 15.5 / 2 * cm))
        plt.plot(r * 1e6, eta * gd, label=r"$\sigma$")
        plt.xlabel(r"$r/\unit{\micro\meter}$")
        plt.ylabel(r"$\sigma/\unit{\pascal}$")
        plt.legend()
        plt.show()
