# -*- coding: utf-8 -*-
"""
File to house the   class.

Created on Tue Apr  8 14:50:35 2025

@author: Richard Kellnberger
"""

import numpy as np


def parseError(ps, cov):
    errs = np.sqrt(np.diag(cov))
    pNew = []
    errNew = []
    for p, err in zip(ps, errs):
        power = -np.floor(np.log10(err))
        pNew.append(np.round(p * 10**power) / 10**power)
        errNew.append(np.round(err * 10**power) / 10**power)
    return np.asarray(pNew), np.asarray(errNew)
