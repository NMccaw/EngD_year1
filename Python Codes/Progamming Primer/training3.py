# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 15:48:31 2016

@author: nm2e16
"""
import math


def degree(x):
    """
    Returns the angle in radians
    """

    rad = x*(360/(2*math.pi))

    return rad


def min_max(x):
    """
    Returns the minimum and maximum of a list as a tuple
    """
    mn = min(x)
    mx = max(x)

    return mn, mx


def geometric_mean(xs):
    """
    Returns the geometric mean of a list
    """
    mu = 1
    for i in xs:
        mu = mu*i
    gm = (mu)**(1/len(xs))
    return gm
