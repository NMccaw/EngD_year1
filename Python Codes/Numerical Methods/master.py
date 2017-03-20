"""
master.py
"""

import numpy as np
import pylab as py
from scipy.optimize import newton
from scipy.integrate import odeint


def L(y, ydot, t, alpha=7/4, beta=5):
    """
    Parameters
    ----------
    y: (N+2,) numpy array of floats
        The output of the machine
        ydot: (N+2,) array of numpy foats
            The time
    t: (N+2,) numpy array of floats
        The time the lagrangian is required for
    dt: float
        The timestep length
    alpha: float (optional)
        Leading constant of proportionality
    beta: float (optional)
        Trailing constant of proportionality
    """
    Lag = alpha * ydot ** 2 + beta * (t ** 2 - 1)  * ydot ** 3 - y
    return Lag

def dLdy(q, t, h=1e-3):
    return (L(q[0] + h, q[1], t) - L(q[0] - h, q[1], t)) / (2 * h)

def dLdydot(y, ydot, t, h=1e-3):
    return (L(y, ydot + h, t) - L(y, ydot - h, t)) / (2 * h)

def ddLdydotdydot(q, t, h=1e-3):
    func = dLdydot
    return (func(q[0], q[1] + h, t) - func(q[0], q[1] - h, t)) / (2 * h)

def ddLdtdydot(q, t, dt, h=1e-3):
    func = dLdydot
    return (func(q[0], q[1], t + h) - func(q[0], q[1], t - h)) / (2 * dt)

def ddLdydydot(q, t, h=1e-3):
    func = dLdydot
    return (func(q[0] + h, q[1], t) - func(q[0] - h, q[1], t)) / (2 * h)

def f(q, t, dt):
    dqdx = np.zeros_like(q)
    dqdx[0] = q[1]
    p1 = dLdy(q, t)
    p2 = ddLdtdydot(q, t, dt)
    p3 = ddLdydydot(q, t)
    p4 = ddLdydotdydot(q, t)
    dqdx[1] = (p1 - p2 - q[1] * p3) / p4
    return dqdx

def phi(z, t, dt):
    soln = odeint(f, [1, z], t, args=(dt,))
    y_boundary = 0.9
    return -(y_boundary - soln[-1, 0])

def shooting(t, dt):
    z_guess = 1
    z_proper = newton(phi, z_guess, args=(t, dt))
    soln = odeint(f, [1, z_proper], t, args=(dt,))
    return soln[:, 0]



if __name__ == "__main__":

    N = 50
    t, dt = np.linspace(0, 1, N + 2, retstep=True)
    ys = shooting(t, dt)

    py.plot(t, ys)
    py.xlabel(r"$t$")
    py.ylabel(r"$y$")
    py.title("Optimal turn-off")

