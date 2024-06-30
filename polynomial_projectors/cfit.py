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

from discrete_orthogonal_polynomials_constructor_and_evaluator.cdop import Cdopeval

# ------------------------------------------------------------------------------


def Cfit(deg, X, f, Y):
    '''
    Given a sample column array f = f (X) of a function at a finite complex set 
    X with card(X) ≥ n + 1, computes the polynomial projector coefficients in 
    an orthogonal polynomial basis at X and evaluates the projector Ln f at a 
    target complex set Y. n is the degree of the polynomial projector.

    Parameters
    ----------
    X : array
        The interpolation points set.
    f : array
        The image of the interpolation points set through a function.
    Y : array
        Evaluation set for the approximation polynomial.

    Returns
    -------
    L: (len(Y), ) array
        The evaluation of the polynomial projector on Y

    '''
    assert len(X) == len(f), "X and f must have the same size"
    assert deg <= len(X), "len(X) must be greater or equal to deg"

    #
    W, Q = Cdopeval(Y, X, deg)
    L = W@Q.conj().T@f

    return L
