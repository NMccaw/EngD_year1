# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 11:43:15 2016

@author: Nicholas
"""

import numpy as np
import scipy.optimize
import pylab


def model(t, Ti, Ta, c):
    """
    Returns a model equation
    """
    return (Ti-Ta)*np.exp(-t/c)+Ta


def read2coldata(filename):
    """
    reads data from a file. Returns data in array form
    """
    f = open(filename, 'rt')
    data = f.readlines()
    array1 = []
    array2 = []
    for dat in data:
        val1, val2 = dat.split()
        array1.append(float(val1))
        array2.append(float(val2))

    return np.asarray(array1), np.asarray(array2)


def extract_parameters(ts, Ts):
    """
    Returns the initial temp, average temp andcooling constant c
    """
    p, pcov = scipy.optimize.curve_fit(model, np.array(ts), np.array(Ts),
                                       p0=[5, 5, 5])
    return Ti, Ta, c


def plot(ts, Ts, Ti, Ta, c):
    """Input Parameters:

      ts : Data for time (ts)
                (numpy array)
      Ts : data for temperature (Ts)
                (numpy arrays)
      Ti : model parameter Ti for Initial Temperature
                (number)
      Ta : model parameter Ta for Ambient Temperature
                (number)
      c  : model parameter c for the time constant
                (number)

    This function will create plot that shows the model fit together
    with the data.

    Function returns None.
    """

    pylab.plot(ts, Ts, 'o', label='data')
    fTs = model(ts, Ti, Ta, c)
    pylab.plot(ts, fTs, label='fitted')
    pylab.legend()
    pylab.show()


def f3(t, Ti, Ta, c):
    return (Ti-Ta)*np.exp(-t/c)+Ta-60


def sixty_degree_time(Ti, Ta, c):
    """
    Returns the time taken to reach 60C
    """
    t60 = scipy.optimize.newton(f3, x0=300, args=(Ti, Ta, c))
    return t60

"""
ts, Ts = read2coldata('time_temp.txt')
Ti, Ta, c = extract_parameters(ts, Ts)
print("Model parameters are Ti=%f C, Tf=%fC," % (Ti, Ta))
print("                     time constant=%fs" % c)
waittime = sixty_degree_time(Ti, Ta, c)
print("The drink reaches 60 degrees after %s seconds = %f minutes"
      % (waittime, waittime / 60.))
"""
