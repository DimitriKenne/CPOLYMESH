# -*- coding: utf-8 -*-
"""
Created on Sat July  20 2024

Updated on Jul7 20, 2024

@author: Dimitri Jordan Kenne
"""

#------------------------------------------------------------------------------
import math

def find_best_grid(m):
    """
    Finds the best grid dimensions for m plots.
    
    Parameters:
    m (int): Number of plots
    
    Returns:
    tuple: (rows, cols)
    """
    # Start with the square root of m
    side = math.ceil(math.sqrt(m))
    
    # Start with assuming square grid
    rows = cols = side
    
    # Adjust rows and columns to accommodate m plots
    while rows * cols >= m:
        cols -= 1
    cols += 1
    
    return rows, cols



