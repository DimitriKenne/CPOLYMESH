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
import matplotlib.pyplot as plt
from polynomial_mesh_constructor.cpom import Cpom
from discrete_extremal_sets_constructor.cdes import Cdes
from lebesgue_constant_evaluator.cleb import Cleb
from domains_structure.examples_of_domains import define_domain

# ------------------------------------------------------------------------------



def demo_cleb(deg, domain, pts_type=['DLP', 'PLP', 'AFP', 'lsqp'],
              adm_mesh_param=4, adm_mesh_param_ep=2):
    '''
    This demo shows the computation of Lebesgue constants for polynomial 
    interpolation and discrete least squares projectors for degree deg

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
    leb = []  # will contain the Lebesgue constants
    
    
    # Compute an Am for evaluating the Lebesgue constant
    mesh, c = Cpom(deg, domain, adm_mesh_param)

    # Compute the interpolation sets of nodes and the admissible mesh of degree
    # deg from which they (DLP or AFP or the last PLP) were extracted
    interp_pts, A = Cdes(deg, domain, adm_mesh_param_ep, pts_type)
    
    print('\n Demo on computing Lebesgue constant\n')
    print(f'\n Domain \t {domain.name}')
    for i in range(n):
        leb_const = Cleb(deg, interp_pts[i], mesh)
        leb.append(leb_const)

        # display the statistics
        rel_error_est = c-1
        abs_error_est = (c-1)*leb_const
        
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
            'degree': deg,
            'Lebesgue Constant AM': leb_const,
            'Relative err. est.': rel_error_est,
            'Absolute err. est.': abs_error_est,
            'Card. admiss. mesh Lesb. const.': len(mesh),
            'number of interpolation points': len(interp_pts[i]),
            'Card mesh for extracting the nodes': len(A)
            }
        
        print("="*80)
        for key, value in statistics.items():
            print(f"{key:<40} {value:<40}")
        print("="*80)


    

fmt = ["-^", '-*', '-o', '-s']
# -----------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Set the degree of approximation
deg = 10

# To select a domain, choose a number from 0 to 30 and input it into the
# define_domain function.
domain = define_domain(20)

# The admissible mesh factor can be modified directly in the call of the demo
# function provided below.

demo_cleb(deg, domain, pts_type=['DLP', 'PLP', 'AFP', 'lsqp'],
          adm_mesh_param=4)