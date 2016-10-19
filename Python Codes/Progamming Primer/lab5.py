# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 17:24:51 2016

@author: nm2e16
"""


def vector_product3(a, b):
    """
    Returns the vector product
    """
    (a1, a2, a3) = a
    (b1, b2, b3) = b
    print(a1)

    return [a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1]


def seq_mult_scalar(a, s):
    """
    Returns a multiplied list
    """
    value_list = []
    for an in a:

        value = an*s
        print(value)
        value_list.append(value)

    return value_list


def powers(n, k):
    """
    Returns the power of a list
    """
    value_list = []
    for i in range(k+1):
        value = n**i
        print(value)
        value_list.append(value)
    return value_list


def traffic_light(load):
    """
    returns a colour given a set of insrtructions
    """
    if load < 0.7:
        return 'green'
    elif load >= 0.7 and load < 0.9:
        return 'amber'
    elif load >= 0.9:
        return 'red'
