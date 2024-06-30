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
# ------------------------------------------------------------------------------


def PLP(deg_interp, adm_mesh, first_point=None):
    """
    Generate the d+1 first points of a pseudo Leja sequence from an Admissible
    mesh (AM)

    Parameters
    ----------
    deg_interp : int
        The degree of interpolation.
    adm_mesh : list
        Contains all the AMs of degree 1,2,...,d. adm_mesh[i] is an AM
        of degree i+1.
        An AM here is an 1D array containing complex number. The first point
        is x0 if it is provided as parameter. Otherwise,
        x0 is the point in adm_mesh[0] with highest imaginary part. The i-th
        point is selected from adm_mesh[i-2]
    first_point : complex, optional
        The first point of the sequence. If it is not provided the first point
        will be a point in adm_mesh which has the highest imaginary part.
        The default is None.

    Returns
    -------
    pseudo_leja_pts : array

    """
    assert len(adm_mesh) >= deg_interp, "len(adm_mesh) must be >= deg_interp"

    if first_point is None:
        i_0 = np.argmax([a.imag for a in adm_mesh[0]])
        first_point = adm_mesh[0][i_0]
    pseudo_leja_pts = np.array([first_point])
    for k in range(0, deg_interp):
        adm_mesh_of_degree_k = adm_mesh[k]
        Z = [np.abs(np.prod(x-pseudo_leja_pts)) for x in adm_mesh_of_degree_k]

        # We choose the index of the maximum
        argmax = np.argmax(Z)
        next_pseudo_leja_pt = adm_mesh_of_degree_k[argmax]
        pseudo_leja_pts = np.append(pseudo_leja_pts, next_pseudo_leja_pt)
    return pseudo_leja_pts
