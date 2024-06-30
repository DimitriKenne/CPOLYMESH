# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 17:49:32 2022

Updated on May 30, 2024

@author: Dimitri Jordan Kenne
-------

Authors:
-------
    - D. J. Kenne (Code Implementation) <dimitri.kenne@doctoral.uj.edu.pl>
    - Alvise Sommariva <alvise@math.unipd.it>
    - Marco Vianello  <marcov@math.unipd.it>

Thank you to:
------------ 
    Prof. Dr. hab. Leokadia Biales-Ciez for her help in developing this work.

License:
-------
This project is licensed under the MIT License - 
see the LICENSE file for details.

"""


# ------------------------------------------------------------------------------
# ------------------------------import packages---------------------------------
# pylint: disable=invalid-name
import numpy as np
from scipy.linalg import qr
from stabilized_vandermonde_matrix_constructor.cvand import Cvand
from discrete_orthogonal_polynomials_constructor_and_evaluator.cdop import Cdop

# ------------------------------------------------------------------------------


def AFP(deg_interp, adm_mesh):
    """
    Generate the deg_interp-th set of Approximate Fekete Points (AFP) using QR
    factorisation.

    Parameters
    ----------
    deg_interp : int
        The degree of interpolation.
    adm_mesh : array
        It is an admissible mesh of degree d, from which the AFP are extracted.


    Returns
    -------
    fekete_pts : array
        The d-th set of AFP.

    """
    # Compute the rectangular Vandermonde matrix and
    # Orthogonalize it (2 iterations)
    Q = Cdop(adm_mesh, deg_interp)[0]

    # Select the indices of the approximate Fekete points
    P = qr(Q.T, pivoting=True)[2]
    ind = P[:deg_interp+1]

    # Extraction of AFP
    fekete_pts = np.array(adm_mesh)[ind]
    return fekete_pts
