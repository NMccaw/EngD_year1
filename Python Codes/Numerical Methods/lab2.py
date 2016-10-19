# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 09:31:44 2016

@author: Nicholas
"""
import numpy as np
import scipy
x=np.linspace(0,10)

def f(s):
    return np.sin(s)**2

C_sol,err= scipy.integrate.quad(f,0,10)
C=1+C_sol


def f_prime(y,x):
    return -C*y



sol=scipy.integrate.odeint(f_prime,1,x)

t= np.linspace(0,500)



sol_x=scipy.integrate.odeint(x_prime,1,t)
sol_y=scipy.integrate.odeint(y_prime,0,t)



r=np.sqrt(sol_x**2+sol_y**2)
plt.plot(t,r)

def x_prime(y,t):
    return -y

def y_prime(x,t):
    return x
