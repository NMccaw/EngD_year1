# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 11:43:06 2016

@author: nm2e16
"""
import math
import pylab
import scipy.optimize


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


def f3(x):
    return (math.cos(2*math.pi*x)*math.exp(-x**2)-(math.log(x+2.1)))


def create_plot_data(f, xmin, xmax, n):
    xs = [(xmin+(i*((xmax-xmin)/(n-1)))) for i in range(n)]
    ys = [f(x) for x in xs]
    return (xs, ys)


def myplot():
    (x1, y1) = create_plot_data(f1, -2, 2, 1001)
    (x2, y2) = create_plot_data(f2, -2, 2, 1001)
    line1 = pylab.plot(x1, y1, label='Function 1')
    line2 = pylab.plot(x2, y2, label='Function 2')
    pylab.xlabel('x')
    pylab.legend()
    pylab.show
    pylab.savefig('plot', format='png')
    pylab.savefig('plot', format='pdf')


def find_cross():
    cross = scipy.optimize.brentq(f3, 0, 1)
    return cross


def reverse_dic(d):
    new_values = [key for key in d.keys()]
    new_key = [val for val in d.values()]
    print(new_key)
    print(new_values)
    r = {}
    for i in range(len(new_key)):
        r[new_key[i]] = new_values[i]
    return r
