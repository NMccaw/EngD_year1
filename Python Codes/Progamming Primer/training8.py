# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 09:32:46 2016

@author: nm2e16
"""

import math


def f1(x):
    """
    returns the output of a function
    """
    return math.cos(2*math.pi*x)*math.exp(-x**2)


def f2(x):
    """
    returns the output of a function
    """
    return math.log(x+2.1)


def positive_places(f, xs):
    """
    returns the numbers of a list which are positive after a function
    """
    num_list = [res for res in xs if f(res) > 0]
    return num_list
