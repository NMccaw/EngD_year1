# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 13:16:23 2016

@author: nm2e16
"""
import math


def seconds2days(n):
    """returns the number of days given the number of seconds"""
    days = n / 60 / 60 / 24
    return days


def box_surface(a, b, c):
    """Returns the surface area of a cuboid"""
    A1 = a*b
    A2 = c*a
    A3 = c*b
    s_area = (A1*2)+(A2*2)+(A3*2)
    return s_area


def triangle_area(a, b, c):
    """Returns the area of a triangle"""
    s = (a+b+c)/2
    A = math.sqrt(s*(s-a)*(s-b)*(s-c))
    return A
