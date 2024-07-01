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

from cfit import Cfit
from discrete_extremal_sets_constructor.cdes import Cdes
from domains_structure.examples_of_domains import define_domain
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------


def demo_cfit(deg, f, domain, pts_type=['DLP', 'PLP', 'AFP', 'lsqp'],
              adm_mesh_param=2,
              grid_size=100):
    '''
    This demo compute the error of approximation for each degree d=1,..., deg,
    for polynomial interpolation at DLP, PLP, AFP extremal points as well as for
    the discrete least squares polynomial of approximation (lsqp). The lsqp is 
    performed on the whole admissible mesh of a fixed degree.
    
    Parameters 
    ----------
    deg : int
        The degree of polynomial approximation
    domain : PolyCurve
        See domain_structure.polynomial_curves for description.
    pts_type : str or list, optional
        The type of interpolation nodes. 
        The default is ['DLP', 'PLP', 'AFP', 'lsqp'].
        The admissible mesh factor. The default is 2.
    adm_mesh_param : TYPE, optional
       The admissible mesh factor. The larger it is, the larger will be the 
       admissible used for extracting extremal sets. The default value is 2.
    grid_size : int, optional
        This is the size of grid that will be used for evaluating the error of 
        approximation.

    Returns
    -------
    None.

    '''
    if isinstance(pts_type, str):
        pts_type = [pts_type]
    n = len(pts_type)
    domain_grid = domain.grid(grid_size)
    fdomain_grid = np.array([f(x) for x in domain_grid])

    errors = [[] for _ in range(n)]  # will contain  log_10(errors)
    # will contain log_10(errors_on_interp_set_pts)
    errors_0 = [[] for _ in range(n)]
    for d in range(1, deg+1):
        # Compute the interpolation sets of nodes
        interp_pts, A = Cdes(d, domain, adm_mesh_param, pts_type=pts_type)
        
        # Plot the extremal points for the degree d=deg
        if d == deg:
            Y_1 = [y.real for y in A]
            Y_2 = [y.imag for y in A]
            for i, X in enumerate(interp_pts):
                if pts_type[i] != 'lsqp':
                    X_1 = [x.real for x in X]
                    X_2 = [x.imag for x in X]
                    plt.plot(X_1, X_2, 'ro', label=pts_type[i], markersize=5)
                    plt.plot(Y_1, Y_2, 'b.', label='AM', markersize=1)
                    plt.legend()
                    plt.axis('equal')
                    plt.show()
        for i, X in enumerate(interp_pts):
            fX = [f(x) for x in X]
            L = Cfit(d, X, fX, domain_grid)
            L_0 = Cfit(d, X, fX, X)
            error = np.max(abs(L-fdomain_grid))
            error_0 = np.max(abs(L_0-fX))
            errors[i].append(np.log10(error))
            errors_0[i].append(np.log10(error_0))

            if pts_type[i] == 'lsqp':
                print('\n Demo discrete least square polynomial')
            else:
                print('\n Demo polynomial interpolation')
                nodes_type= pts_type[i]
            statistics = {
                'interp_points' : nodes_type,
                'degree' : d,
                'Approximation error': error,
                'Approx. error on interp. pts set': error_0
                }
            print("="*80)
            for key, value in statistics.items():
               print(f"{key:<40} {value:<40}")
            print("="*80)
           
    # Plot the domain
        
    # Plot approximation errors
    degrees = np.arange(1, deg+1)
    for i in range(len(pts_type)):
        plt.plot(degrees, errors[i], fmt[i], label=pts_type[i], markersize=5)
    plt.xlabel('degree')
    plt.ylabel('log_10(error)')
    plt.legend()
    plt.axis('equal')
    plt.show()


fmt = ["^", '*', 'o', 's']
# ------------------------------------------------------------------------------

# Define a function for approximation


def f(x): return np.exp(x-2)


# Set the degree of approximation
deg = 10

# To select a prebuilt domain, choose a number from 1 to 19 and input it into
# the define_domain function.
domain = define_domain(31)

# The admissible mesh factor can be modified directly in the call of the demo
# function provided below.


demo_cfit(deg, f, domain, pts_type=['DLP', 'PLP', 'AFP', 'lsqp'],
          adm_mesh_param=2, grid_size=1000)
