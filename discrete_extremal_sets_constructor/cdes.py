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
# pylint: disable=dangerous-default-value
# pylint: disable=consider-using-enumerate

from discrete_extremal_sets_constructor.plp import PLP
from discrete_extremal_sets_constructor.dlp import DLP
from discrete_extremal_sets_constructor.afp import AFP
from discrete_extremal_sets_constructor.leja_pts_unit_circle import leja_pts_unit_circle
from polynomial_mesh_constructor.cpom import Cpom

# ------------------------------------------------------------------------------


def Cdes(deg, domain, adm_mesh_param=2, pts_type=['DLP', 'PLP', 'AFP'],
         gauss_cheb_mesh=True):
    """
    Compute three set of interpolation nodes. We have three options:
    pseudo Leja points ('PLP'), discrete Leja points ('DLP')
    and approximate Fekete points ('AFP')

    Parameters
    ----------
    deg : int
        The degree of interpolation.
    domain : PolyParamCurve or UnionPolyParamCurve
        It is a parametrisation of the domain.
        See the definitions in define_polynomial_curves.py
    adm_mesh_param : float, optional
        It is the AM factor. The default is 2.
    pts_type : str or list of str, optional
        It is the type of the nodes.Choose between 'PLP', DLP', 'AFP'. 
        The default is ['DLP', 'PLP', 'AFP']. 
    gauss_cheb_mesh : bool, optional
        If True then the Gauss-Chebyshev points will be used. Otherwise it is
        the Chebyshev_lobatto points that are used. The default is True.

    Returns
    -------
    interp_pts: list of arrays
            The i-th array contains the interpolation nodes (complex number) 
            corresponding to pts_type[i].
    adm_mesh: list or array of complex numbers
        This is the admissible mesh of degree deg where the interpolation points 
        were eventually extracted.

    """

    adm_mesh = Cpom(deg, domain, adm_mesh_param, gauss_cheb_mesh,
                    only_deg=False)[0]
    if isinstance(pts_type, str):
        pts_type = [pts_type]
    interp_pts = []
    for i in range(len(pts_type)):
        if pts_type[i] == 'Leja_pts_unit_circle':
            leja_pts = leja_pts_unit_circle(deg)
            interp_pts.append(leja_pts)
        if pts_type[i] == 'PLP':
            pseudo_leja_pts = PLP(deg, adm_mesh)
            interp_pts.append(pseudo_leja_pts)
        if pts_type[i] == 'DLP':
            discrete_leja_pts = DLP(deg, adm_mesh[-1])
            interp_pts.append(discrete_leja_pts)
        if pts_type[i] == 'AFP':
            approx_fekete_pts = AFP(deg, adm_mesh[-1])
            interp_pts.append(approx_fekete_pts)
        if pts_type[i] == 'lsqp':
            interp_pts.append(adm_mesh[-1])
    return interp_pts, adm_mesh[-1]
