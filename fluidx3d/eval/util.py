# -*- coding: utf-8 -*-
"""
File to house the   class.

Created on Wed Mar 26 15:00:30 2025

@author: Richard Kellnberger
"""

import numpy as np


def promote(A):
    """
    Help broadcasting by inserting unneccessary .

    Assume you have a Tensor field with shape (Lx, Ly, Lz, d) and intend to multiply it
    by a scalar field with shape (Lx, Ly, Lz).
    The scalar field needs to be promoted to (Lx, Ly, Lz, 1) first.
    Otherwise the operation fails.

    Parameters
    ----------
    A : np.ndarray
        Array to be promoted.

    Returns
    -------
    np.ndarray
        Promoted array.

    """
    return np.expand_dims(A, -1)


def strainRateTensorToShearRate(D):
    """
    Convert the strain-rate tensor to the shear-rate.

    According to https://doi.org/10.1103/PhysRevApplied.19.064061
    equation (2)

    Parameters
    ----------
    D : np.ndarray
        The strain-rate tensor.

    Returns
    -------
    gd : np.ndarray
        The shear-rate.

    """
    diag = D[..., :3]
    offDiag = D[..., 3:]

    gd = np.sqrt(2 * (np.sum(diag**2, -1) + 2 * np.sum(offDiag**2, -1)))
    return gd


def stressTensorToVonMisesStress(S):
    """
    Convert the Cauchy stress tensor into the von Mises stress.

    According to https://doi.org/10.1103/PhysRevApplied.19.064061
    equation (11)

    Parameters
    ----------
    S : np.ndarray
        The Cauchy stress tensor.

    Returns
    -------
    sigmaVM : np.ndarray
        The von Mises stress.

    """
    diag = S[..., :3]
    deviatoric = S - promote(np.sum(diag, -1) / 3)

    deviatoricDiag = deviatoric[..., :3]
    deviatoricOffDiag = deviatoric[..., 3:]
    sigmaVM = np.sqrt(
        3 / 2 * (np.sum(deviatoricDiag**2, -1) + 2 * np.sum(deviatoricOffDiag**2, -1))
    )
    return sigmaVM
