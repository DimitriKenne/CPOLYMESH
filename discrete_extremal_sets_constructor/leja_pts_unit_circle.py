# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 17:49:32 2022

Updated on May 30, 2024

@author: Dimitri Jordan Kenne
-------

"""


# ------------------------------------------------------------------------------
# ------------------------------import packages---------------------------------
# pylint: disable=invalid-name

import numpy as np
# ------------------------------------------------------------------------------


def leja_pts_unit_circle(deg):
    '''
    Compute the Leja points on the disk starting from 1

    Parameters
    ----------
    deg : int
        The degree of interpolation.
    Returns
    -------
    E : TYPE
        DESCRIPTION.

    '''
    E = [1]
    for n in range(1, deg+1):
        v = binary_form(n)
        s = len(v)
        # create a list of 2^{-k} k = s-1, ..., 0
        u = 1/2**(np.array(range(s-1, -1, -1)))
        # compute e_n
        e = np.exp(1j*np.pi*np.dot(u, v))
        # print(u, v, np.dot(u, v))
        E.append(e)
    return E


def binary_form(n, length=None):
    '''
    Convert a number to its binary decomposition

    Parameters
    ----------
    n : int
    length : int, optional
        Tells how many digit of the decomposition to return.
        The default is None.

    Returns
    -------
    list
        The elements of the list form the binary representation of n.

    '''
    # Convert the number to binary and remove the '0b' prefix
    binary_str = bin(n)[2:]
    # If a length is specified, pad the binary string with leading zeros
    if length:
        binary_str = binary_str.zfill(length)
    # Convert the binary string into a list of characters
    binary_list = list(binary_str)
    # Convert each character in the list to an integer
    return [int(digit) for digit in binary_list]
