# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 02:40:17 2023

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
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring
# pylint: disable=too-many-arguments
#
import numpy as np
import matplotlib.pyplot as plt


class PolyCurve:
    '''
    This class represents a polynomial curve or arc, either algebraic or 
    trigonometric.

    Attributes:
    -----------
    polynomial: function
        Defines the curve.
    interval: tuple or list (2 elements)
        The definition interval of the polynomial, where 
        a=interval[0] < b=interval[1].
    degree: int
        The degree of the polynomial.
    curve_type: str
        Indicates if the polynomial is 'alg' (algebraic) or 'trig' 
        (trigonometric).
    name: str, optional
        The curve's name.

    Methods:
    --------
    grid(N):
        Creates a grid of N points along the curve, mapping the uniform grid 
        on [a,b].
    disk_enclosing():
        Determines the center and radius of a disk that encloses the curve.
    plot_curve():
        Plots the curve.

    '''

    def __init__(self, polynomial, interval, degree, curve_type, name=None):
        '''
        Initializes a polynomial curve given a parametrization (polynomial),
        its interval of definition, its degree and its polynomial type. 
        The name of the curve is optional 

        '''
        self.name = name
        self.polynomial = polynomial
        self.degree = degree
        self.interval = interval
        self.curve_type = curve_type

    def grid(self, N):
        '''
        Generate a grid for self as image of the uniform grid of  self.interval

        Parameters
        ----------
        N : int
            The size of the uniform grid used.

        Returns
        -------
        grid : list of complex numbers

        '''
        a, b = self.interval
        T = np.linspace(a, b, N)
        Z = [self.polynomial(t) for t in T]
        return Z

    def disk_enclosing(self, N=50):
        '''Determine the center and the radius of the cercle enclosing self'''
        grid = self.grid(N)
        Zc = np.average(grid)
        Zr = max(abs(grid-Zc))
        if Zr == 0:
            Zr = 1
        return (Zc, Zr)

    def plot_curve(self, N=50, color='m', display=True):
        '''Plot the curve defined by self'''
        grid = self.grid(N)
        X = [x.real for x in grid]
        Y = [x.imag for x in grid]

        plt.plot(X, Y, color)

        if display:
            print('I am here', self.name)
            plt.title(self.name)
            plt.axis('equal')
            plt.show()


class UnionPolyCurves(PolyCurve):
    '''
    This subclass of PolyCurve represents a union of polynomial curves or arcs, 
    either algebraic or trigonometric.

    Attributes:
    -----------
    list_of_curves: list
        This is the list of polynomial curves (instances of Polycurve) which
        form the domain
    name: str, optional
        The curve's name.

    Methods:
    --------
    Inherit all methods of the PolyCurve class.

    '''

    def __init__(self, list_of_polycurves, name=None):
        '''
        Initializes a union of polynomial curves given a parametrization 
        (polynomial), its interval of definition, its degree and its polynomial 
        type. The name of the curve is optional 

        '''
        polynomials = [curve.polynomial for curve in list_of_polycurves]
        intervals = [curve.interval for curve in list_of_polycurves]
        degrees = [curve.degree for curve in list_of_polycurves]
        curve_types = [curve.curve_type for curve in list_of_polycurves]
        super().__init__(polynomials, intervals, degrees, curve_types, name)
        self.curves = list_of_polycurves

    def grid(self, N):
        '''
        Generate a grid for self as image of the uniform grid of  self.interval

        Parameters
        ----------
        N : int
            The size of the uniform grid used.

        Returns
        -------
        grid : list of complex numbers

        '''
        grid = []
        for curve in self.curves:
            grid = list(set(grid) | set(curve.grid(N)))
        return grid

    def plot_curve(self, N=100, color='m', display=True):
        '''Plot the curve'''
        for curve in self.curves:
            curve.plot_curve(N, display=False)

        if display:
            plt.title(self.name)
            plt.axis('equal')
            plt.show()

# -----------------------------------------------------------------------------
