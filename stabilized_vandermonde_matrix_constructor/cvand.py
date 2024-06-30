# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 17:49:32 2022

Updated on May 30, 2024

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

# ------------------------------------------------------------------------------


def Cvand(pts, center=0, radius=1, number_columns=None):
    """
    Compute the vandermonde-like matrix at the points in pts using
    a shifted monomial basis i.e. monomial basis of degree "deg", centered in 
    zB and with radius "delta", that is the (k+1) column is given by 
    "(pts-zB)^k/delta".

    Parameters
    ----------
    pts : array
        1D array containing complex numbers at which the Vandermonde
        matrix is computed.
    grid: array, optional
        It is a discretisation of a domain that contains the points in pts.
        The default is None.
    number_colum : int, optional
        The number of columns for the matrix. Each column corresponds to
        a polynomial in the polynomial basis. If it is not provided then
        number_columns = len(pts). The default is None.

    Returns
    -------
    V : array
        (len(pts), number_colum) array.

    """
    if number_columns is None:
        number_columns = len(pts)
    if center is None or radius is None:
        center = np.average(pts)
        radius = max(abs(pts-center))

    V = np.array([shifted_monomial(d, pts, center, radius)
                  for d in np.arange(number_columns)]).T
    return V


def shifted_monomial(deg, x, center=0, radius=1):
    '''  Evaluate at x, a shifted monomial of degree 'deg', that is 
    ((x-center)/radius)**deg   
    '''
    x = np.array(x)
    return ((x-center)/radius)**deg
