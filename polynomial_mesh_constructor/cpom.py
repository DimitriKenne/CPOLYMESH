# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 13:25:41 2023

Updated on May 30, 2024

@author: Dimitri Jordan Kenne

This module focuses on the computation of admissible meshes for domains
within the complex plane, bounded by curves or arcs defined by algebraic or
trigonometric polynomial parametrizations. Chebyshev points are utilized
for this task.


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


# -----------------------------------------------------------------------------
# ----------------------- Import the necessary packages------------------------


# pylint: disable=invalid-name
from domains_structure.polynomial_curves import UnionPolyCurves
import numpy as np
import copy


# -----------------------------------------------------------------------------
# ----------------------------AM for polynomial curves--------------------------


def Cpom(deg, domain, adm_mesh_param=2, gauss_cheb_mesh=True,
         only_deg=True):
    """

    Compute an admissible mesh (AM) for domain, where domain must be an 
    (algebraic or trigonometric) polynomial curve
    or a union of polynomial curves (all algebraic or all trigonometric or 
    mixte of algebraic and trigonometric).
    ------------------------------------------

    Parameters
    ----------
    deg : int
        The degree of the mesh.

    domain: PolyCurve
        It is a the domain. 
        See domain_structure.py for description.

    adm_mesh_param : float, optional
        The AM factor refers to the constant in the AM inequality. 
        The greater this constant, the larger the mesh becomes. 
        The default is 2.

    gauss_cheb_mesh : bool, optional
        If set to True, Gauss-Chebyshev points are utilized. If not, 
        Chebyshev-Lobatto points are employed. The default setting is True.


    only_deg : bool, optional
        If set to True, a single set of points representing the AM of degree 
        deg is returned. If False, a list containing deg sets of points is 
        returned, where the i-th element corresponds to the AM of degree i. 
        The default setting is True.

    Returns
    -------
    tuple (2 elements)

        mesh: array or list
            It is an array if only_deg = True and a list of arrays otherwise. 
            mesh[0] is the set of points that form the AM. 
        adm_mesh_const: float
            It is the constant associated with the AM
"""
    if isinstance(domain.degree,  list):  # UnionPolyCurses case
        curves = domain.curves
        mesh = Cpom(deg, curves[0], adm_mesh_param, gauss_cheb_mesh, only_deg)
        for i in range(1, len(curves)):
            mesh_new = Cpom(deg, curves[i], adm_mesh_param,
                            gauss_cheb_mesh, only_deg)
            mesh = union_adm_meshes(mesh, mesh_new, only_deg)
    else:
        if domain.curve_type == 'alg':
            mesh = adm_mesh_alg_poly_curve(deg, domain, adm_mesh_param,
                                           gauss_cheb_mesh, only_deg)
        if domain.curve_type == 'trig':
            mesh = adm_mesh_trig_poly_curve(deg, domain, adm_mesh_param,
                                            gauss_cheb_mesh, only_deg)
    return mesh


# -----------------------------------------------------------------------------
# -----------------------AM for algebraic polynomial curves--------------------


def adm_mesh_alg_poly_curve(deg, domain, adm_mesh_param=2,
                            gauss_cheb_mesh=True, only_deg=True):
    """


    Compute an AM for domain, where domain must be an algebraic polynomial 
    parametric curve
    ------------------------------------------

    Parameters
    ----------
    deg : int
        The degree.

    domain: PolyCurve
        It is a the domain. 
        See domain_structure.py for description.

    adm_mesh_param : float, optional
        The AM factor refers to the constant in the AM inequality. 
        The greater this constant, the larger the mesh becomes. 
        The default is 2.

    gauss_cheb_mesh : bool, optional
        If set to True, Gauss-Chebyshev points are utilized. If not, 
        Chebyshev-Lobatto points are employed. The default setting is True.


    only_deg : bool, optional
        If set to True, a single set of points representing the AM of degree 
        deg is returned. If False, a list containing deg sets of points is 
        returned, where the i-th element corresponds to the AM of degree i. 
        The default setting is True.

    Returns
    -------
    tuple (2 elements)

        mesh: array or list
            It is an array if only_deg = True and a list of arrays otherwise. 
            mesh[0] is the set of points that form the AM. 
        adm_mesh_const: float
            It is the constant associated with the AM


    """
    assert domain.curve_type == "alg", \
        "The domain in question is not characterized by an algebraic \
        polynomial curve, but rather by a trigonometric one."

    parametrization = domain.polynomial
    k = domain.degree
    a, b = domain.interval
    adm_mesh_const = 1 / np.cos(np.pi/(2*adm_mesh_param))
    if only_deg:
        if gauss_cheb_mesh:
            cheb_pts = gauss_cheb_pts(adm_mesh_param*deg*k, a=a, b=b)
        else:
            cheb_pts = cheb_lobato_pts(adm_mesh_param*deg*k, a=a, b=b)
        # Compute the image of Chebyshev points using the given parametrization.
        mesh = np.array([parametrization(t) for t in cheb_pts])
        return (mesh, adm_mesh_const)

    # In case only_deg is False
    meshes = []
    for n in range(1, deg+1):
        mesh = adm_mesh_alg_poly_curve(n, domain, adm_mesh_param,
                                       gauss_cheb_mesh)
        meshes.append(mesh[0])
    return (meshes, adm_mesh_const)


# -----------------------------------------------------------------------------
# -------------------AM for trigonometric polynomial curves--------------------


def adm_mesh_trig_poly_curve(deg, domain, adm_mesh_param=2,
                             gauss_cheb_mesh=True, only_deg=True):
    """

    Compute an AM for the curve domain
    ------------------------------------------

    Parameters
    ----------
    deg : int
        The degree.

    domain: PolyCurve
        It is a the domain. 
        See domain_structure.py for description.

    adm_mesh_param : float, optional
        The AM factor refers to the constant in the AM inequality. 
        The greater this constant, the larger the mesh becomes. 
        The default is 2.

    gauss_cheb_mesh : bool, optional
        If set to True, Gauss-Chebyshev points are utilized. If not, 
        Chebyshev-Lobatto points are employed. The default setting is True.


    only_deg : bool, optional
        If set to True, a single set of points representing the AM of degree 
        deg is returned. If False, a list containing deg sets of points is 
        returned, where the i-th element corresponds to the AM of degree i. 
        The default setting is True.

    Returns
    -------
   tuple (2 elements)

        mesh: array or list
            It is an array if only_deg = True and a list of arrays otherwise. 
            mesh[0] is the set of points that form the AM. 
        adm_mesh_const: float
            It is the constant associated with the AM

    """

    assert domain.curve_type == "trig", \
        "The domain in question is not characterized by an trigonometric\
            polynomial curve, but rather by a algebraic one."

    parametrization = domain.polynomial
    k = domain.degree
    a, b = domain.interval
    adm_mesh_const = 1 / np.cos(np.pi/(2*adm_mesh_param))
    if only_deg:
        if gauss_cheb_mesh:
            cheb_pts = gauss_cheb_pts(2*adm_mesh_param*deg*k)
        else:
            cheb_pts = cheb_lobato_pts(2*adm_mesh_param*deg*k)
        pre_mesh = sigma(cheb_pts, a, b)
        mesh = np.array([parametrization(t) for t in pre_mesh])
        return (mesh, adm_mesh_const)
    # In case only_deg is False
    meshes = []
    for n in range(1, deg+1):
        mesh = adm_mesh_trig_poly_curve(n, domain, adm_mesh_param,
                                        gauss_cheb_mesh)
        meshes.append(mesh[0])
    return (meshes, adm_mesh_const)


# ----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# --------------------------Gauss-Chebyshev points-----------------------------

def gauss_cheb_pts(N, a=-1, b=1):
    """ 
    Generate the N-th Gauss-Chebyshev points (N points) inside the segment 
    [a,b] 
    ---------------------------------------------------------------------------

    Parameter: 
        --------------
        N: int
           N is the number of points that will be generated
        a, b: complex  
           There are the vertices of the segment to be considered 

    return:
        ----------
        X: array
          The number in X are the Chebyshev points

    """

    k = np.arange(1, N+1)
    X = 0.5*(a+b) + 0.5*(b-a)*np.cos((2*k-1)*np.pi/(2*N))
    X.sort()
    return X


# -----------------------------------------------------------------------------
# -----------------------Chebyshev-Lobatto points------------------------------


def cheb_lobato_pts(N, a=-1, b=1):
    """ 
    Generate the N-th Chebyshev-Lobatto points (N+1 points) inside the segment 
    [a,b] 
    ---------------------------------------------------------------------------

    Parameters: 
        --------------
        N: int
           N+1 is the number of points that will be generated
        a, b: complex numbers 
           There are the vertices of the segment to be considered 

    returns:
        ----------
        X: array
          the number in X are the Chebyshev points

    """

    k = np.arange(0, N+1)
    X = 0.5*(a+b) + 0.5*(b-a)*np.cos(k*np.pi/N)
    X.sort()
    return X

# -----------------------------------------------------------------------------


def sigma(u, a, b):
    """
    Will be useful for the computation of an AM in the case of trigonometric 
    polynomial curve

    """

    R = 2*np.arcsin(u*np.sin(0.25*(b-a))) + 0.5*(a+b)
    return R

# ------------------------------------------------------------------------------


def union_list(lst1, lst2):
    """
    Union of to lists without repetition of an element

    Parameters
    ----------
    lst1 : list

    lst2 : list

    Returns
    -------
    final_list : list

    """

    final_list = list(set(lst1) | set(lst2))
    return final_list


def union_adm_meshes(A1, A2, only_deg=True):
    """
    Unify two admissible meshes.

    Parameters
    ----------
    A1 : array or list
        It contains the first AM with its constant C.

    A2 : array or list
        It contains the first AM with its constant C.
        A2 must have the same properties with A1.

    only_deg : bool, optional
        If A1 and A2 are arrays containing complex number write True.
        But if A1 and A2 are list of array then write False. 
        The default is True.

    Returns
    -------
    A : list (2 elements)

        A[0]: array or list
            It is an array if only_deg = True and a list of arrays otherwise.

        A[1]: float
            It is the constant associated with the AM

    """
    C = max(A1[1], A2[1])
    if only_deg:
        B = union_list(A1[0], A2[0])
    else:
        B1 = A1[0]
        B2 = A2[0]
        assert len(B1) == len(B2), "B1 and B2 must have the same lenght"

        B = []
        for j in range(len(B1)):
            D = union_list(B1[j], B2[j])
            B.append(np.array(D))
    return (B, C)


# -----------------------------------------------------------------------------
