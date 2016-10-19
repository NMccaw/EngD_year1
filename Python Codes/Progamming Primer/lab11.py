# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 15:50:33 2016

@author: Nicholas
"""
import math
from scipy.optimize import fmin
import scipy.integrate
import scipy.interpolate


def f(x):
    """
    returns the result of the function
    """
    return math.exp(-x**2)/(1+x**2)+(2*math.cos(x)**2)/(1+(x-4)**2)


def make_multiplier(factor):
    """
    returns factor multiplied by x
    """
    return lambda x: factor * x


def integrate_f_from0(b):
    """
    returns the integral of a function from 0 to assigned value b
    """
    integral, err = scipy.integrate.quad(f, 0, b)
    return integral


def g(x):
    return -f(x)


def find_max_f():
    """
    returns the maximum value of a function
    """
    fmax = fmin(g, 2)
    return fmax[0]


def f2(x):
    return (math.exp(-x**2)/(1+x**2)+(2*math.cos(x)**2)/(1+(x-4)**2))-1


def find_f_equals_1():
    """
    returns the value of x where f(x)=1
    """
    f1 = scipy.optimize.brentq(f2, 0, -10)
    return f1


def lin_int(xs, ys):
    """
    returns a callable function which returns an array using linear
    interpolation
    """
    return scipy.interpolate.interp1d(xs, ys)


def make_oscillator(frequency):
    """
    returns the sin of time *frequency
    """
    return lambda t: math.sin(t*frequency)
