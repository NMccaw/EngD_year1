# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 12:12:29 2016

@author: nm2e16
"""


def box_volume(a, b, c):
    """
    Returns the volume of a box
    """
    bv = a * b * c
    return bv


def fall_time(h):
    """
    Returns the time for an object to fall a defined distance
    """
    g = 9.81
    time = (2*(h/g))**0.5
    return time


def interval_point(a, b, x):
    """
    Returns how far to go from points a and b
    """
    distogo = a + (abs(a-b)*x)
    return distogo


def impact_velocity(h):
    """
    Returns the impact velocity
    """
    g = 9.81
    time = (2*(h/g))**0.5
    v = g*time
    return v


def signum(x):
    """
    Returns a value of a 1,0 or -1 if given conditions are met.
    """
    if x > 0:
        ans = 1
    elif x == 0:
        ans = 0
    else:
        ans = -1

    return ans