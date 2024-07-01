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

from discrete_extremal_sets_constructor.cdes import Cdes
from domains_structure.examples_of_domains import define_domain
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------------------


def demo_cdes(deg, domain, pts_type=np.array(['DLP', 'PLP', 'AFP']),
                    adm_mesh_param=2):
    interp_pts, adm_mesh = Cdes(deg, domain, adm_mesh_param=adm_mesh_param,
                                pts_type=pts_type)
    Y_1 = [y.real for y in adm_mesh]
    Y_2 = [y.imag for y in adm_mesh]
    for i, X in enumerate(interp_pts):
        X_1 = [x.real for x in X]
        X_2 = [x.imag for x in X]
        plt.plot(X_1, X_2, 'ro', label=pts_type[i], markersize=5)
        plt.plot(Y_1, Y_2, 'b.', label='AM', markersize=1)
        plt.legend()
        plt.axis('equal')
        plt.show()


# Set the degree of approximation
deg = 50

# To select a prebuilt domain, choose a number from 0 to 30 and input it into
# the define_domain function.
domain = define_domain(5)


demo_cdes(deg, domain)
