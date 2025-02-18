# -*- coding: utf-8 -*-
"""
File to house the   class.

Created on Tue Feb 18 16:25:29 2025

@author: Richard Kellnberger
"""

import os
from pathlib import Path


def separate():
    path = Path(os.getcwd())
    back = Path("")
    while path.name != "python":
        back = path.name / back
        path = path.parent
    return path.parent, back


def python():
    p, _ = separate()
    return p / "python"


def data():
    p, _ = separate()
    return p / "data"


def plots():
    p, _ = separate()
    return p / "plots"


def cPython():
    p, b = separate()
    return p / "python" / b


def cData():
    p, b = separate()
    return p / "data" / b


def cPlots():
    p, b = separate()
    return p / "plots" / b
