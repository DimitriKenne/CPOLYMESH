# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 17:20:13 2024

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
from domains_structure.find_best_grid import find_best_grid
from discrete_extremal_sets_constructor.cdes import Cdes
from lebesgue_constant_evaluator.cleb import Cleb
from polynomial_mesh_constructor.cpom import Cpom
import matplotlib.pyplot as plt
import numpy as np
import os

# ------------------------------------------------------------------------------


def demo(deg, domain, pts_type=['DLP', 'PLP', 'AFP', 'lsqp'],
              adm_mesh_param=4, adm_mesh_param_ep=2):
    '''
    This demo shows the computation of Lebesgue constants for polynomial 
    interpolation and discrete least squares projectors for 
    degrees= 1,..., deg. It also plot the considered interpolation points.

    Parameters
    ----------
    deg : int
        The degree of approximation.
    domain : PolyCurve
        See domain_structure.polynomial_curves for description.
    pts_type : str or list, optional
        The type of interpolation nodes. 
        The default is ['DLP', 'PLP', 'AFP', 'lsqp'].
    adm_mesh_param : TYPE, optional
       The admissible factor involved in estimating the Lebesgue constant is 
       significant. The larger it is, the better the estimation of the Lebesgue 
       constants will be. The default value is 4.
    adm_mesh_param_ep : int, optional
        The admissible mesh factor for extracting the interpolation nodes
        The default is 2.

    Returns
    -------
    None.

    '''
    n = len(pts_type)
    leb = [[] for _ in range(n)]  # will contain the Lebesgue constants

    

    for d in range(1, deg+1):
        # Compute an Am for evaluating the Lebesgue constant
        mesh = Cpom(d, domain, adm_mesh_param)

        # Compute the interpolation sets of nodes
        interp_pts, A = Cdes(d, domain, adm_mesh_param_ep, pts_type)

        # Compute Lebesgue constants
        print('\n Demo on computing Lebesgue constant\n')
        for i in range(n):
            leb_const = Cleb(d, interp_pts[i], mesh[0])
            leb[i].append(leb_const)

            # display the statistics
            rel_error_est = mesh[1]-1
            abs_error_est = (mesh[1]-1)*leb_const
            
            if pts_type[i] == 'PLP':
                interp_nodes =  'Pseudo Leja points'
            if pts_type[i] == 'DLP':
                interp_nodes = 'discrete Leja points'
            if pts_type[i] == 'AFP':
                interp_nodes = 'approximate Fekete points'
            if pts_type[i] == 'lsqp':
                interp_nodes =  'Discrete least square approximation'
            
            statistics = {
                'domain': domain.name,
                'Interp. nodes': interp_nodes,
                'degree': d,
                'Lebesgue Constant AM': leb_const,
                'Relative err. est.': rel_error_est,
                'Absolute err. est.': abs_error_est,
                'Mesh Leb. const. degree': d,
                'Mesh Leb. const. factor': adm_mesh_param,
                'Mesh Leb. const. constant': mesh[1],
                'Card. mesh Leb. const.': len(mesh[0]),
                'number of interpolation points': len(interp_pts[i]),
                'Mesh interp. nodes degree': d,
                'Mesh interp. nodes factor': adm_mesh_param_ep,
                'Card mesh for extracting the nodes': len(A),
                }
            
            print("="*80)
            for key, value in statistics.items():
               print(f"{key:<40} {value:<40}")
            print("="*80)

    
    # Define the folder where you want to save the plots
    output_folder = "figures"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Plot the extremal points for the degree deg
    Y_1 = [y.real for y in A]
    Y_2 = [y.imag for y in A]

    ## Create a figure with subplots
    m = len([pt for pt in pts_type if pt != "lsqp"]) + 1 
    rows, cols = find_best_grid(m)
    fig, axs = plt.subplots( rows, cols, figsize=(12, 8))
    axs = axs.ravel()  # Flatten the array of axes to make indexing easier
    k = 0 # To track the plots already added to our figure

    for i, X in enumerate(interp_pts):
        if pts_type[i] != 'lsqp':
            X_1 = [x.real for x in X]
            X_2 = [x.imag for x in X]
            
            # Plot the extremal points of type pts_type[i]
            axs[k].plot(X_1, X_2, 'ro', label=pts_type[i], markersize=5)

            # Plot the admissible mesh on the same figure
            axs[k].plot(Y_1, Y_2, 'b.', label='AM', markersize=1)
            axs[k].legend()
            axs[k].axis('equal')
            # axs[k].set_title('Extremal Points')
            k +=1
            
    # Plot Lebesgue constants
    degrees = np.arange(1, deg + 1)
    for i in range(len(pts_type)):
        axs[k].plot(degrees, leb[i], fmt[2], label = pts_type[i], markersize=5)
    axs[k].set_xlabel('Degree')
    axs[k].legend()
    axs[k].axis('equal')
    #axs[k].set_ylim(top=deg)
    axs[k].set_xlim([0,deg+1])
    axs[k].set_title('Lebesgue constant')
    axs[k].grid(True)
    
    # Adjust layout
    fig.tight_layout()

    # Save the figure
    plot_path = os.path.join(output_folder, 'demo_plots.png')
    fig.savefig(plot_path)

    # Optionally display the figure
    plt.show()



fmt = ["-^", '-*', '-o', '-s']
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Set the degree of approximation
deg = 10

# To select a domain, choose a number from 0 to 31 and input it into the
# define_domain function.
domain = define_domain(31)

# The admissible mesh factor can be modified directly in the call of the demo
# function provided below.

demo(deg, domain, pts_type=['DLP', 'PLP', 'AFP', 'lsqp'],
          adm_mesh_param=4)
