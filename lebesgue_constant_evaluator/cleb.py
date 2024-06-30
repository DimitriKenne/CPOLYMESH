# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 23:21:10 2024

Updated on June 1, 2024

@author: Dimitri Jordan Kenne


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
from discrete_orthogonal_polynomials_constructor_and_evaluator.cdop import Cdopeval

# ------------------------------------------------------------------------------


def Cleb(n, X, Z):
    '''
    Computes on a control set Z the maximum of the Lebesgue function of 
    interpolation on a set X with card(X) = n + 1 or least-squares with 
    card(X) > n + 1.

    Parameters
    ----------
    n : int
        The degree of the polynomial approximant
    X : array
        The interpolation node set
    Z : array
        The control set which is a discretization of the domain.

    Returns
    -------
    float.

    '''

    W, Q = Cdopeval(Z, X, n)
    return np.linalg.norm(W@Q.conj().T, ord=np.inf)
