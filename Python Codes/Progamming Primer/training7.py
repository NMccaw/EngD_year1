# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 15:54:54 2016

@author: nm2e16
"""


def count_chars(s):
    """
    Returns the number of times characters occurs in a string in dictionary
    form
    """
    d = {}
    for char in s:
        if char not in d.keys():
            d[char] = s.count(char)
    return d


def derivative(f, x):
    """
    returns the central fininte difference of a function
    """
    eps = 10**-6
    der = (f(x+(eps/2))-f(x-(eps/2)))/eps
    return der


def derivative(f, x, eps=10**-6):
    """
    returns the central fininte difference of a function with variable accuracy
    """
    der = (f(x+(eps/2))-f(x-(eps/2)))/eps
    return der
