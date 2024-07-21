# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 13:25:41 2023

Updated on July 21, 2024

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
import os

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
    
    # Create a directory for figures if it doesn't exist
    if not os.path.exists('figures'):
        os.makedirs('figures')
    
    # Plot the admissible mesh
    if only_deg:
        plot_adm_mesh(adm_mesh[0], deg=deg, save=True, 
                      filename='figures/adm_mesh_deg_{}.png'.format(deg))
    else:
        for d, mesh in enumerate(adm_mesh[0]):
            plot_adm_mesh(mesh, deg=d+1, save=True, 
                          filename='figures/adm_mesh_deg_{}.png'.format(d+1),
                          display_block=False)
        plt.show()
    
    # Display the result
    card_adm_mesh = len(adm_mesh[0]) if only_deg else len(adm_mesh[0][-1])
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

def plot_adm_mesh(pts, deg=None, save=False, filename=None, display=True, display_block=True):
    '''
    Plot the complex points inside pts

    Parameters
    ----------
    pts : array
        The set of points to plot.
    save : bool, optional
        Specify whether to save the plot. The default is False.
    filename : str, optional
        The filename to save the plot as. Required if save is True.

    Returns
    -------
    None.
    '''
    X = [x.real for x in pts]
    Y = [x.imag for x in pts]
    plt.figure()
    plt.plot(X, Y, 'bo', label='AM', markersize=1)
    if deg is not None:
        plt.title(f'AM of degree {deg}')
    plt.legend()
    plt.axis('equal')
    
    if save:
        if filename is None:
            raise ValueError("Filename must be specified if save is True.")
        plt.savefig(filename)
    
    if display:    
        plt.show(block=display_block)

# -----------------------------------------------------------------------------
# Set the degree of approximation
deg = 3

# To select a prebuilt domain, choose a number from 0 to 31 and input it into
# the define_domain function.
domain = define_domain(30)

# The admissible mesh factor can be modified directly in the call of the demo
# function provided below.
demo_cpom(deg, domain, adm_mesh_param=4, only_deg=True)
