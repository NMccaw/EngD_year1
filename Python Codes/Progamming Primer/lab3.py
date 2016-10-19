# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 16:24:23 2016

@author: nm2e16
"""
import math
g = 9.81


def swing_time(L):
    """
    Returns the swing time
    """
    T = (2*math.pi)*math.sqrt(L/g)
    return T


def range_squared(n):
    """
    Returns the range squared
    """
    r = []
    for i in range(n):
        r.append(i**2)
    return r


def count(element, seq):
    """
    returns the number of times an element displays in a sequence
    """
    c = 0
    for i in range(len(seq)):
        if seq[i] == element:
            c = c+1
    return c
