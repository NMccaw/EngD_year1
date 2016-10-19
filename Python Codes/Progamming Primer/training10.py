# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 11:33:48 2016

@author: Nicholas
"""
import numpy as np


def model(t, Ti, Ta, c):
    """
    Returns a model equation
    """
    T = (Ti-Ta)*np.exp(-t/c)+Ta
    return T
