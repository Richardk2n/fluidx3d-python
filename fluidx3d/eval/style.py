# -*- coding: utf-8 -*-
"""
File to house the style defenition.

Created on Fri Apr 14 14:49:37 2023

@author: Richard Kellnberger
"""

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 12})
plt.rcParams.update({"text.usetex": True})
plt.rcParams.update({"figure.autolayout": True})
plt.rcParams.update(
    {
        "text.latex.preamble": "".join(
            [
                r"\usepackage[locale=US]{siunitx}"
                r"\sisetup{per-mode=fraction}"
                r"\sisetup{separate-uncertainty=true}"
                r"\usepackage{amsmath}"
                r"\newcommand{\abs}[1]{\lvert#1\rvert}"
                r"\newcommand{\ma}[1]{\underline{#1}}"
                r"\DeclareMathOperator{\Tr}{Tr}"
            ]
        )
    }
)
plt.rcParams["image.origin"] = "lower"


def rgb(r, g, b):
    return r / 255, g / 255, b / 255


activeStyle = "light"
cm = 1 / 2.54
tabColors = list(mcolors.TABLEAU_COLORS.keys())


def color(index):
    index %= len(tabColors)
    return tabColors[index]


def style(style: str, r, g, b):
    global activeStyle
    activeStyle = style
    plt.rcParams["figure.facecolor"] = rgb(r, g, b)
    plt.rcParams["axes.facecolor"] = rgb(r, g, b)
    plt.rcParams["savefig.facecolor"] = rgb(r, g, b)


def dark_():
    style("dark", 51, 51, 51)


def dark():
    plt.style.use("dark_background")
    dark_()


def light_():
    style("light", 255, 255, 255)


def light():
    plt.style.use("default")
    light_()


def savefig(name: str, **kwargs):
    global activeStyle
    if activeStyle == "light":
        plt.savefig(f"{name}#light.pdf", **kwargs)
    elif activeStyle == "dark":
        plt.savefig(f"{name}#dark.pdf", **kwargs)
    else:
        plt.savefig(f"{name}.pdf", **kwargs)


def styled(plotFunction):
    with plt.style.context("default"):
        light_()
        plotFunction()
        plt.close()
    with plt.style.context("dark_background"):
        dark_()
        plotFunction()
        plt.show()


def styledS(plotFunction):
    with plt.style.context("default"):
        light_()
        plotFunction()
        plt.close()
    with plt.style.context("dark_background"):
        dark_()
        plotFunction()
        plt.close()


def figure(**kwargs):
    return plt.figure(figsize=(15.5 * cm, 15.5 / 2 * cm), **kwargs)


dark()
