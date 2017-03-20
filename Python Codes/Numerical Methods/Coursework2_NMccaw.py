#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 17:17:45 2016

@author: Nicholas

Justification for shooting method.

Here the Shooting method and Relaxation method will be described breifly then
justification for using the shooting method will be made described.

Firstly the shooting method involves creating an initial value problem (IVP) 
from the boundary value problem (BVP). This is done by refomulating the IVP 
as y" = f(t,y',y") with initial conditions y(0) = 1 and the y'(0) = z ,where z
is an unknown value. The value of z is found by finding the value of y(1) = y(1)
of the BVP. In this code this is done by minimizing the function phi(z) = y(1) - 0.9
using newton minimization. Once the value of z is found the system is integrated
using odeint to to find the value of y.

For relaxation the a grid evenly spaced over the interval is placed. All derivatives
are approximated using finite differences with the end points defined by the boundary
conditions. A linear system is then constructed from the finite difference 
approximation equations. The linear system is of the form Ty=F which, given the 
Dirchlet boundary conditions can be solved. However the problem given is non-linear.
A linear system cannot therefore be constructed. Instead a guess for the y points
is produced then the function f is calculated, using the y vector, is then calculated
using a linear system. The process is repeated until the residuals of the system 
are minimized. To minimize the residuals the Newtons method is required. However
this requires the construction of the jacobian matrix. The construction of the 
jacobian matrix requires the partial differential of the residual w.r.t y. This 
will add additional computational cost and errors due to finite differences.

Therefore for the system in the problem the shooting method was choosen due to
the simplicity of implementation. If, however, the shooting method did not work 
the relaxation method would be implemented. The shooting method would fail if:
the root finding process didn't converge or if the IVP was unstable.

For the alpha = 7/4, beta = 5 case it was found that the value of y was over 1.
This suggests that for the stated condition for the optimum output drop down
the machine will need to be over 100%. This is not practical. Therefore the 
mathematical problem has not been posed correctly to solve the original problem.
"""

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import newton
from matplotlib import pyplot as plt
import pytest

def shooting(z,t,alpha,beta,h,f,y_0_boundary):
    """
    Shooting functionl. Uses odeint to compute the y function. Takes the initial
    conditions for the Initial Value problem.
    
    Parameters
    ----------
    
    z: float
        value of y'(0)
    t: array
        time series
    alpha: float
            constant to define penalty function
    beta: float
        constant to define penalty function
    h: float
        step size used for integration
    y_0_boundary: float
                value of y(0)
    
    
    Returns
    -------
    array of float
    the y values of the function. Computed using odeint
    
    """ 
    assert type(z) == np.float64 or np.float, ' z is not of type float'
    assert type(t) == np.array or np.ndarray, 't is not of type np.array'
    assert type(alpha) == np.float64 or np.float, ' alpha is not of type float'
    assert type(beta) == np.float64 or np.float, ' beta is not of type float'
    assert type(h) == np.float64 or np.float, ' h is not of type float'

    return odeint(f,[y_0_boundary,z],t,args = (alpha,beta,h))

def phi(z,alpha,beta,h,f,y_0_boundary,y_boundary = 0.9):
    """
    Phi function. Function output to be minimized by newton minimization. Takes
    in a guess for y'(0) and subtracts from known boundary at y(1). This value 
    should be minimized such that for the initial y'(0) the final value of y
    is equal to the boundary condition at y(1)
    
    Parameters
    ----------
    
    z: float
        value of y'(0)
    alpha: float
            constant to define penalty function
    beta: float
        constant to define penalty function
    h: float
        step size used for integration
    
    
    Returns
    -------
    float
    the last value of y computed by shooting - the boundary at y(1). 
    """
    assert type(z) == np.float64 or np.float, ' z is not of type float'
    assert type(alpha) == np.float64 or np.float, ' alpha is not of type float'
    assert type(beta) == np.float64 or np.float, ' beta is not of type float'
    assert type(h) == np.float64 or np.float, ' h is not of type float'
    # compute y using shooting method
    y = shooting(z,[0,1],alpha,beta,h,f,y_0_boundary)
    # subtract y(1) of shooting method and the value the boundary value problem states
    return y[-1,0] - y_boundary
    

    
def L_fun(t,q,alpha,beta,h,y):
    """
    Function used to define the action
    
    Parameters
    ----------
    t: float
        time
    q: float
        value of y'
    alpha: float
            constant to define penalty function
    beta: float
        constant to define penalty function
    h: float
        step size used for integration
    y: float
        value of y
    
    
    Returns
    -------
    L: float
        the value of the action equation for given y' , y and t
    """
    assert type(q) == np.float64 or np.float, ' q is not of type float'
    assert type(t) == np.float64 or np.float, ' t is not of type float'
    assert type(alpha) == np.float64 or np.float, ' alpha is not of type float'
    assert type(beta) == np.float64 or np.float, ' beta is not of type float'
    assert type(h) == np.float64 or np.float, ' h is not of type float'
    assert type(y) == np.float64 or np.float, ' y is not of type float'
    # define L
    L = (alpha * q**2) +(beta*(t**2 -1) * q**3) - y
    
    return L
    

def dLdy_fun(q,t,alpha,beta,h,y):
    """
    Function used to define the derivative of the action wrt y
    
    Parameters
    ----------
    q: array
        value of y and y'
    t: float
        time
    alpha: float
            constant to define penalty function
    beta: float
        constant to define penalty function
    h: float
        step size used for integration
    y: float
        value of y
    
    
    Returns
    -------
    dLdy: float
        the value of the derivative of the action equation wrt y for given y' , y and t
    """
    assert type(q) == np.float64 or np.float, ' q is not of type float'
    assert type(t) == np.float64 or np.float, ' t is not of type float'
    assert type(alpha) == np.float64 or np.float, ' alpha is not of type float'
    assert type(beta) == np.float64 or np.float, ' beta is not of type float'
    assert type(h) == np.float64 or np.float, ' h is not of type float'
    assert type(y) == np.float64 or np.float, ' y is not of type float'
    # central differencing of L wrt y
    dLdy = (L_fun(t,q[1],alpha,beta,h,q[0]+h) - L_fun(t,q[1],alpha,beta,h,q[0]-h)) / (2*h)
    
    return dLdy
   
    
def d2Ldydq_fun(q,t,alpha,beta,h,y):
    """
    Function used to define the derivative of the action wrt y and y'
    
    Parameters
    ----------
    q: array
        value of y and y'
    t: float
        time
    alpha: float
            constant to define penalty function
    beta: float
        constant to define penalty function
    h: float
        step size used for integration
    y: float
        value of y
    
    
    Returns
    -------
    d2Ldydq: float
        the value of the derivative of the action equation wrt y and y' for given y' , y and t
    """
    assert type(q) == np.array or np.ndarray, ' q is not of type array'
    assert type(t) == np.float64 or np.float, ' t is not of type float'
    assert type(alpha) == np.float64 or np.float, ' alpha is not of type float'
    assert type(beta) == np.float64 or np.float, ' beta is not of type float'
    assert type(h) == np.float64 or np.float, ' h is not of type float'
    assert type(y) == np.float64 or np.float, ' y is not of type float'
    # central differencing of L wrt y and ydot
    d2Ldydq = (dLdq_fun(t,q[1],alpha,beta,h,q[0]+h) - dLdq_fun(t,q[1],alpha,beta,h,q[0]-h)) / (2*h)

    return d2Ldydq

def d2Ldtdq_fun(q,t,alpha,beta,h,y):  
    """
    Function used to define the derivative of the action wrt t and y'
    
    Parameters
    ----------
    q: array
        value of y and y'
    t: float
        time
    alpha: float
            constant to define penalty function
    beta: float
        constant to define penalty function
    h: float
        step size used for integration
    y: float
        value of y
    
    
    Returns
    -------
    : float
        the value of the derivative of the action equation wrt t and y' for given y' , y and t
    """
    assert type(q) == np.array or np.ndarray, ' q is not of type array'
    assert type(t) == np.float64 or np.float, ' t is not of type float'
    assert type(alpha) == np.float64 or np.float, ' alpha is not of type float'
    assert type(beta) == np.float64 or np.float, ' beta is not of type float'
    assert type(h) == np.float64 or np.float, ' h is not of type float'
    assert type(y) == np.float64 or np.float, ' y is not of type float'
    # central differencing of L wrt t and ydot
    return (dLdq_fun(t+h,q[1],alpha,beta,h,q[0])-dLdq_fun(t-h,q[1],alpha,beta,h,q[0]))/(2*h)
    
def dLdq_fun(t,q,alpha,beta,h,y):
    """
    Function used to define the derivative of the action wrt y'
    
    Parameters
    ----------
    q: float
        value of y'
    t: float
        time
    alpha: float
            constant to define penalty function
    beta: float
        constant to define penalty function
    h: float
        step size used for integration
    y: float
        value of y
    
    
    Returns
    -------
    : float
        the value of the derivative of the action equation wrt y' for given y' , y and t
    """
    assert type(q) == np.float64 or np.float, ' q is not of type float'
    assert type(t) == np.float64 or np.float, ' t is not of type float'
    assert type(alpha) == np.float64 or np.float, ' alpha is not of type float'
    assert type(beta) == np.float64 or np.float, ' beta is not of type float'
    assert type(h) == np.float64 or np.float, ' h is not of type float'
    assert type(y) == np.float64 or np.float, ' y is not of type float'
    # central differencing of L wrt  ydot
    return (L_fun(t,(q+h),alpha,beta,h,y) - L_fun(t,(q-h),alpha,beta,h,y)) / (2*h)
    
def d2Ldqdq_fun(q,t,alpha,beta,h,y):
    """
    Function used to define the 2nd derivative of the action wrt y'
    
    Parameters
    ----------
    q: array
        value of y and y'
    t: float
        time
    alpha: float
            constant to define penalty function
    beta: float
        constant to define penalty function
    h: float
        step size used for integration
    y: float
        value of y
    
    
    Returns
    -------
    : float
        the value of the 2nd derivative of the action equation wrt y' for given y' , y and t
    """
    assert type(q) == np.array or np.ndarray, ' q is not of type array'
    assert type(t) == np.float64 or np.float, ' t is not of type float'
    assert type(alpha) == np.float64 or np.float, ' alpha is not of type float'
    assert type(beta) == np.float64 or np.float, ' beta is not of type float'
    assert type(h) == np.float64 or np.float, ' h is not of type float'
    assert type(y) == np.float64 or np.float, ' y is not of type float'
    # central differencing of L wrt ydot and ydot
    return (dLdq_fun(t,q[1]+h,alpha,beta,h,q[0]) - dLdq_fun(t,q[1]-h,alpha,beta,h,q[0]) )/(2*h)
    
    
    
def f(q,t,alpha,beta,h):
    """
    Function used to define the function which y'' is equal to.
    
    Parameters
    ----------
    q: array
        value of y and y'
    t: float
        time
    alpha: float
            constant to define penalty function
    beta: float
        constant to define penalty function
    h: float
        step size used for integration
    y: float
        value of y
    
    
    Returns
    -------
    dqdt: array
        the value of y' and y''
    """
    assert type(q) == np.array or np.ndarray, ' q is not of type array'
    assert type(t) == np.float64 or np.float, ' t is not of type float'
    assert type(alpha) == np.float64 or np.float, ' alpha is not of type float'
    assert type(beta) == np.float64 or np.float, ' beta is not of type float'
    assert type(h) == np.float64 or np.float, ' h is not of type float'

    #q = [y,ydot]
    dqdt = np.zeros_like(q)
    #dqdt = [ydot,ydotdot]
    dqdt[0] = q[1]
    y = q[0]
    dLdy1  = dLdy_fun(q,t,alpha,beta,h,y)
    d2Ldydq = d2Ldydq_fun(q,t,alpha,beta,h,y)
    d2Ldtdq = d2Ldtdq_fun(q,t,alpha,beta,h,y)
    d2Ldqdq = d2Ldqdq_fun(q,t,alpha,beta,h,y)
    # define ydotdot = f(t,y,ydot)
    dqdt[1] = (dLdy1 - (q[1]* d2Ldydq) - d2Ldtdq) / d2Ldqdq
    
    return dqdt
 
    
    

if __name__ == '__main__':
    pytest.main('-v Coursework2_NMccaw.py')
    # set values for constants
    t, dt = np.linspace(0,1,retstep = True)
    alpha = 5
    beta = 5
    h = 1e-3
    y_0_boundary = 1
    y_boundary = 0.9
    # find ydot(0) using newton minimization
    z_true = newton(phi,1.0,args = (alpha,beta,h,f,y_0_boundary,y_boundary))
    # compute y using shooting method
    y =shooting(z_true,t,alpha,beta,h,f,y_0_boundary)
    #plots
    plt.figure(1)
    plt.plot(t,y[:,0])
    plt.ylim(0.88,1.02)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$y$")
    plt.title(r'Optimal machine turn-down $\alpha = {}$ $\beta = {}$'.format(alpha,beta))
    # set values for constants
    alpha = 7/4
    beta = 5
    h = 1e-3
    # find ydot(0) using newton minimization
    z_true = newton(phi,1.0,args = (alpha,beta,h,f,y_0_boundary,y_boundary))
    # compute y using shooting method
    y =shooting(z_true,t,alpha,beta,h,f,y_0_boundary)  
    # plots
    plt.figure(2)
    plt.plot(t,y[:,0])
    plt.ylim(0.88,1.02)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$y$")
    plt.title(r'Optimal machine turn-down $\alpha = {}$ $\beta = {}$'.format(alpha,beta))
     
    y_convergence = np.zeros((8,50))
    h_convergence = np.zeros(8)
    error = np.zeros(8)
    for i in range(8):
        # set values of constants
        t = np.linspace(0,1,50)
        alpha = 7/4
        beta = 5        
        h = 0.05/(2**i)
        y_0_boundary = 1
        y_boundary = 0.9
        # find ydot(0) using newton minimization
        z_true = newton(phi,1.0,args = (alpha,beta,h,f,y_0_boundary,y_boundary))
        # compute y using shooting method
        y =shooting(z_true,t,alpha,beta,h,f,y_0_boundary)
        y_convergence[i,:] = y[:,0] 
        h_convergence[i] = h
        # compute error
        error[i] = np.linalg.norm(y_convergence[i,:]-y_convergence[i-1,:],axis= 0) 
    #plots
    plt.figure(5)
    m ,p= np.polyfit(np.log(h_convergence[1:]),np.log(error[1:]),1)
    line_fit_y = np.exp(p)*h_convergence[1:]**m
    plt.loglog(h_convergence[1:],line_fit_y,label = 'Line gradient = {0:.2f}'.format(m))  
    plt.loglog(h_convergence[1:],error[1:],'kx',label = 'Numerical points')
    plt.xlabel('h')
    plt.ylabel('Error')
    plt.title('Error of solution relative to high accuracy solution \n with changing step size')
    plt.legend()

    

##############Tests##############
   
def test_L_fun():
    alpha = 5
    q = 1 
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    L = L_fun(t,q,alpha,beta,h,y)
    assert np.allclose(L,0.3),'L function not defined correctly'
    
def test_dLdy_fun():
    alpha = 5
    q = [0.95,1]
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    dLdy = dLdy_fun(q,t,alpha,beta,h,y)
    assert np.allclose(dLdy,-1),'dLdy function not defined correctly'
    
def test_d2Ldydq_fun():
    alpha = 5
    q = [0.95,1]
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    d2Ldydq = d2Ldydq_fun(q,t,alpha,beta,h,y)
    assert np.allclose(d2Ldydq,0),'d2Ldydq function not defined correctly'
    
def test_d2Ldtdq_fun():
    alpha = 5
    q = [0.95,1]
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    d2Ldtdq = d2Ldtdq_fun(q,t,alpha,beta,h,y)
    
    assert np.allclose(d2Ldtdq,15,rtol=1e-2),'d2Ldtdq function not defined correctly'\

def test_dLdq_fun():
    alpha = 5
    q = 1 
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    dLdq = dLdq_fun(t,q,alpha,beta,h,y)

    assert np.allclose(dLdq,-1.25,rtol=1e-2),'dLdq function not defined correctly'
    
def test_d2Ldqdq_fun():
    alpha = 5
    q = [0.95,1] 
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    d2Ldqdq = d2Ldqdq_fun(q,t,alpha,beta,h,y)
    assert np.allclose(d2Ldqdq,-12.5,rtol=1e-2),'d2Ldqdq function not defined correctly'

def test_L_fun_2():
    alpha = 0.5
    q = 2 
    beta = 2
    t = 1
    h = 1e-3
    y = 0.9
    L = L_fun(t,q,alpha,beta,h,y)
    assert np.allclose(L,1.1),'L function not defined correctly'
    
def test_dLdy_fun_2():
    alpha = 0.5
    q = [0.9,2] 
    beta = 2
    t = 1
    h = 1e-3
    y = 0.9
    dLdy = dLdy_fun(q,t,alpha,beta,h,y)
    assert np.allclose(dLdy,-1),'dLdy function not defined correctly'
    
def test_d2Ldydq_fun_2():
    alpha = 0.5
    q = [0.9,2] 
    beta = 2
    t = 1
    h = 1e-3
    y = 0.9
    d2Ldydq = d2Ldydq_fun(q,t,alpha,beta,h,y)
    assert np.allclose(d2Ldydq,0),'d2Ldydq function not defined correctly'
    
def test_d2Ldtdq_fun_2():
    alpha = 0.5
    q = [0.9,2] 
    beta = 2
    t = 1
    h = 1e-3
    y = 0.9
    d2Ldtdq = d2Ldtdq_fun(q,t,alpha,beta,h,y)
    
    assert np.allclose(d2Ldtdq,48,rtol=1e-2),'d2Ldtdq function not defined correctly'\

def test_dLdq_fun_2():
    alpha = 0.5
    q = 2 
    beta = 2
    t = 1
    h = 1e-3
    y = 0.9
    dLdq = dLdq_fun(t,q,alpha,beta,h,y)

    assert np.allclose(dLdq,2,rtol=1e-2),'dLdq function not defined correctly'
    
def test_d2Ldqdq_fun_2():
    alpha = 0.5
    q = [0.9,2] 
    beta = 2
    t = 1
    h = 1e-3
    y = 0.9
    d2Ldqdq = d2Ldqdq_fun(q,t,alpha,beta,h,y)
    assert np.allclose(d2Ldqdq,1,rtol=1e-2),'d2Ldqdq function not defined correctly'
    
    
def test_shooting_known():
    z = 0
    alpha =1
    h = 1e-3
    beta = 5
    y_0_boundary = 0
    t_series = np.linspace(0,1)
    def f_known(q,t,alpha,beta,h):
        dqdt = np.zeros_like(q)
        dqdt[0] = q[1]
        dqdt[1] = 2
        return dqdt
    y = shooting(z,t_series,alpha,beta,h,f_known,y_0_boundary)
    assert np.allclose(y[:,0],t_series**2)

def test_phi_known():
    z = 0
    alpha =1
    beta= 1
    h = 1e-3
    y_0_boundary = 0
    y_boundary = 1
    def f_known(q,t,alpha,beta,h):
        dqdt = np.zeros_like(q)
        dqdt[0] = q[1]
        dqdt[1] = 2
        return dqdt
    zero= phi(z,alpha,beta,h,f_known,y_0_boundary,y_boundary)

    assert np.allclose(zero,0)

def test_raise_error_L_fun_str():
    alpha = 'hi'
    q = 1 
    beta = 5
    t = 0.5
    h = 1e-3
    y = 'hi'
    with pytest.raises(TypeError):
        L_fun(t,q,alpha,beta,h,y)
    
def test_raise_error_dLdy_fun_str():
    alpha = 'hi'
    q = [1,0.9] 
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    with pytest.raises(TypeError):
        dLdy_fun(q,t,alpha,beta,h,y)
        
def test_raise_error_d2Ldydq_fun_str():
    alpha = 'hi'
    q = [1,0.9] 
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    with pytest.raises(TypeError):
       d2Ldydq_fun(q,t,alpha,beta,h,y)
    
def test_raise_error_d2Ldtdq_fun_str():
    alpha = 'hi'
    q = [1,0.9] 
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    with pytest.raises(TypeError):
       d2Ldtdq_fun(q,t,alpha,beta,h,y)

def test_raise_error_dLdq_fun_str():
    alpha = 'hi'
    q = 1 
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    with pytest.raises(TypeError):
       dLdq_fun(t,q,alpha,beta,h,y)
       
def test_raise_error_d2Ldqdq_fun_str():
    alpha = 'hi'
    q = [1,3] 
    beta = 5
    t = 0.5
    h = 1e-3
    y = 0.95
    with pytest.raises(TypeError):
       d2Ldqdq_fun(q,t,alpha,beta,h,y)
 
def test_raise_error_f_str():
    alpha = 5
    q = [1,3] 
    beta = 'hi'
    t = 0.5
    h = 1e-3
    
    with pytest.raises(TypeError):
       f(q,t,alpha,beta,h)
      

       
       
      




