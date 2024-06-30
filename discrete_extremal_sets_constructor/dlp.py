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
# pylint: disable=unused-variable

import numpy as np
from scipy.linalg import lu
from stabilized_vandermonde_matrix_constructor.cvand import Cvand
from discrete_orthogonal_polynomials_constructor_and_evaluator.cdop import Cdop
# ------------------------------------------------------------------------------


def DLP(deg_interp, adm_mesh, lu_factorization=True,
        first_point=None):
    """

    Generate the deg_interp+1 first points of a discrete Leja sequence from an
    Admissible mesh (AM)

    Parameters
    ----------
    deg_interp : int
        The degree of interpolation.
    adm_mesh : array
        The AM. It elements are complex numbers.
    lu_factorisation : bool, optional
        If True then The LU factorization with standard row pivoting
        is use for extracting the points. Otherwise, the points a extract
        by maximizing the product of distances between them and the previous
        points.
        The default is True.
    first_point : complex, optional
        The first point of the sequence. It is used when lu_factorization is
        set to False. If it is not specified the first point will be a point 
        in adm_mesh which has the highest imaginary part. The default is None.

    Returns
    -------
    discrete_leja_pts : array

    """

    adm_mesh_size = len(adm_mesh)
    assert adm_mesh_size > deg_interp, "len(adm_mesh) must be >= deg_interp"
    if lu_factorization is True:
        # Create the rectangular Vandermonde matrix and
        # Orthogonalize it (2 iterations)
        V = Cdop(adm_mesh, deg_interp)[0]

        # LU factorisation of the rectangular Vandermonde matrix and get the
        # permutation vector
        P = lu(V)[0]
        perm_vector = P.T@np.arange(0, adm_mesh_size)

        # Extract the indices from perm_vector and extract the Leja points
        indices = perm_vector[:deg_interp+1].astype(int)
        discrete_leja_pts = np.array(adm_mesh)[indices]
    else:
        print('Discrete Leja points computed without using LU factoriszation')
        # Extraction without using the LU factorisation
        if first_point is None:
            i_0 = np.argmax([a.imag for a in adm_mesh])
            first_point = adm_mesh[i_0]
        discrete_leja_pts = np.array([first_point])
        indices = np.array([0])
        for j in range(1, deg_interp+1):
            Z = [np.abs(np.prod(x-discrete_leja_pts)) for x in adm_mesh]
            argmax = np.argmax(Z)
            indices = np.append(indices, argmax)
            discrete_leja_pts = np.append(discrete_leja_pts, adm_mesh[argmax])
    return discrete_leja_pts
