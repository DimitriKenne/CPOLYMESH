# -*- coding: utf-8 -*-
"""
Created on Thu May 30 06:30:43 2024

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
from stabilized_vandermonde_matrix_constructor.cvand import Cvand

# ------------------------------------------------------------------------------


def Cdop(pts, deg, domain=None):
    '''
    Computes a discrete orthogonal polynomial basis on the set of points pts 
    and return the corresponding unitary matrix Q and two matrices R_1 and 
    R_2 which define the transformation matrix T = inv(R_1)*inv(R_2) 

    Parameters
    ----------
    pts : list or array
        The set of nodes for polynomial approximation.
    deg : int
        The degree of approximation.
    domain : PolyCurve, optional
        See domaine-structure.polynomial_curves for description

    Returns
    -------
    Q : (len(pts), deg) array
        The unitary matrix.
    R_1 : (deg, deg) array

    R_2 : (len(pts), deg) array

    Z_c : complex
        The center of a disk that encloses domain.
    r : float
        The radius of a disk that encloses domain.

    '''
    if domain is None:
        Z_c = np.average(pts)
        r = max(abs(pts-Z_c))
    else:
        Z_c, r = domain.disk_enclosing()

    V = Cvand(pts, center=Z_c, radius=r, number_columns=deg+1)
    Q_1, R_1 = np.linalg.qr(V)
    Q, R_2 = np.linalg.qr(Q_1)
    return (Q, R_1, R_2, Z_c, r)


def Cdopeval(evaluation_set, pts, deg, domain=None):
    '''
    Evaluates the orthogonal basis, computed via Cdop, on the target complex 
    set evaluation_set

    Parameters
    ----------
    evaluation_set : list or array
        The set of evaluation.
    pts : list or array
        The set of nodes for polynomial approximation.
    deg : int
        The degree of approximation.
    domain : PolyCurve, optional
        See domaine-structure.polynomial_curves for description

    Returns
    -------
    W : (deg+1, deg+1) array
        The evaluation of the orthogonal basis on evaluation_set.
    Q : (len(pts), deg+1) array
        This is the unitary matrix

    '''
    Q, R_1, R_2, Z_c, r = Cdop(pts, deg, domain)
    V = Cvand(evaluation_set, center=Z_c, radius=r, number_columns=deg+1)
    W_1 = right_division(V, R_1)
    W = right_division(W_1, R_2)
    return (W, Q)


def right_division(B, A):
    '''
    Solve XA = B using lstsq on the transposed system

    Parameters
    ----------
    B : array

    A : array


    Returns
    -------
    X : array

    '''
    # Since B/A in MATLAB means X*A = B, we transpose both A and B to
    # use lstsq
    # A.T * X.T = B.T => X.T = np.linalg.lstsq(A.T, B.T, rcond=None)
    # this works if A.T is a square matrice or is (m,n) matrix with m>n

    X_T = np.linalg.lstsq(A.T, B.T, rcond=None)[0]
    X = X_T.T  # Transpose back to the original shape

    # We could also use following line instead
    # X = B@ np.linalg.pinv(A)

    return X
