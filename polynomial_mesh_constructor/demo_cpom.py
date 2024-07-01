# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 13:25:41 2023

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

from domains_structure.examples_of_domains import define_domain
from polynomial_mesh_constructor.cpom import Cpom
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------


def demo_cpom(deg, domain, adm_mesh_param=4, only_deg=True):
    '''
    This demo shows how to compute an admissible mesh 


    Parameters
    ----------
    deg : int
        The degree of the admissible mesh
    domain : PolyCurve
        Choose a number between 1 and 2. The default is 1.
    adm_mesh_param : int, optional
        The admissible mesh factor. The default is 4.
    only_deg : bool, optional
        If set to True, a single set of points representing the AM of degree 
        deg is returned. If False, a list containing deg sets of points is 
        returned, where the i-th element corresponds to the AM of degree i. 
        The default setting is True.

    Returns
    -------
    None.

    '''

    # Compute the admissible mesh
    adm_mesh = Cpom(deg, domain, adm_mesh_param, only_deg=only_deg)
    print('\t The admissible mesh is\n', adm_mesh[0])

    # Plot the admissible mesh
    if only_deg == True:
        plot_adm_mesh(adm_mesh[0], deg=deg)
    else:
        for d, mesh in enumerate(adm_mesh[0]):
            plot_adm_mesh(mesh, deg=d+1)

    # Display the result
    
    if only_deg == True:
        card_adm_mesh = len(adm_mesh[0])
    else:
       card_adm_mesh = len(adm_mesh[0][-1])
    statistics = {
        'domain': domain.name,
        'degree': deg,
        'AM factor': adm_mesh_param,
        'AM constant': adm_mesh[1],
        'AM size': card_adm_mesh
        }
    print('\n \t Demo of computing admissible meshes')
    print("="*60)
    for key, value in statistics.items():
        print(f"{key:<20} {value:<40}")
    print("="*60)


def plot_adm_mesh(pts, deg=None, display=True):
    '''
    Plot the complex points inside pts

    Parameters
    ----------
    pts : array
        The set of points to plot.
    display : bool, optional
        Specify whether to plt.show() or not. The default is True.

    Returns
    -------
    None.

    '''
    X = [x.real for x in pts]
    Y = [x.imag for x in pts]
    plt.plot(X, Y, 'bo', label='AM', markersize=1)
    if display == True:
        if deg is not None:
            plt.title(f'AM of degree {deg}')
        plt.legend()
        plt.axis('equal')
        plt.show()

# -----------------------------------------------------------------------------


# Set the degree of approximation
deg = 10

# To select a prebuilt domain, choose a number from 0 to 30 and input it into
# the define_domain function.
domain = define_domain(30)

# The admissible mesh factor can be modified directly in the call of the demo
# function provided below.
demo_cpom(deg, domain, adm_mesh_param=4, only_deg=True)
