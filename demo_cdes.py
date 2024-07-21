# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 17:49:32 2022

Updated on July 21, 2024

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

from domains_structure.examples_of_domains import define_domain
from domains_structure.find_best_grid import find_best_grid
from discrete_extremal_sets_constructor.cdes import Cdes
import matplotlib.pyplot as plt
import numpy as np
import os

# ------------------------------------------------------------------------------


def demo_cdes(deg, domain, pts_type=np.array(['DLP', 'PLP', 'AFP']),
                    adm_mesh_param=2):
    interp_pts, adm_mesh = Cdes(deg, domain, adm_mesh_param=adm_mesh_param,
                                pts_type=pts_type)
        
    # Ensure the 'figures' directory exists
    if not os.path.exists('figures'):
        os.makedirs('figures')

    # Create subplots
    m = len(pts_type)
    rows, cols = find_best_grid(m)
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    # Plot the adm. mesh 
    Y_1 = [y.real for y in adm_mesh]
    Y_2 = [y.imag for y in adm_mesh]

    # Plot the interp_pts points on each subplot
    for i, X in enumerate(interp_pts):
        X_1 = [x.real for x in X]
        X_2 = [x.imag for x in X]
        ax = axs[i // 2, i % 2]
        ax.plot(X_1, X_2, 'ro', label=pts_type[i], markersize=5)
        ax.plot(Y_1, Y_2, 'b.', label='AM', markersize=1)
        ax.legend()
        ax.axis('equal')

    # Hide the unused subplot (bottom right) if there are only 3 plots
    if len(interp_pts) < 4:
        fig.delaxes(axs[1, 1])

    # Adjust layout
    plt.tight_layout()

    # Save the figure
    plt.savefig('figures/cdes_plots.png')

    # Show the figure
    plt.show()


# Set the degree of approximation
deg = 10

# To select a prebuilt domain, choose a number from 0 to 31 and input it into
# the define_domain function.
domain = define_domain(5)


demo_cdes(deg, domain)
