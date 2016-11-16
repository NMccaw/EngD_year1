# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 11:30:28 2016

@author: Nicholas

Solutions to Coursework 1. Using 3rd order Runge-Kutta method and Gauss-Radau 
Runge-Kutta method to find integral of the Modified Prothero-Robinson Problem.
"""
import numpy as np
import scipy.optimize
from matplotlib import pyplot
from numba import jit
import warnings
#options =[gamma,eps,omega]


@jit
def f(t,q,options):
    """
    
    Define the Modified Prothero-Robinson Problem
    
    Parameters
    ----------
    
    t : float 
        time at timestep
    q : 2 X 1 array
        x and y values of previous timestep      
    options: 3X1 array
        gives either stiff or non stiff values for gamma , epsilon and omega
    
        
    Returns
    -------
    qout: 2 X 1 array
        x and y values at timestep
    """
    assert type(t) == np.float or int or np.float64, 't is not type float'
    assert type(q) == np.ndarray, 'q is not type array'
    assert len(q) == 2, 'q is not of length 2'
    assert type(options) == np.ndarray or list, 'options is not of type nd.array'
    assert len(options) == 3,'3 values not given for options'
    
    # unpack q to give x and y values
    x,y = q
    # unpack options to set gamma, epsilon and omega values
    gamma, epsilon,omega = options
    # compute scaling matrix
    scaling = np.matrix([[gamma, epsilon], [epsilon, -1]])
    # compute 2nd matrix of equation
    equ1 = np.matrix([[(-1 + x**2 - np.cos(t)) / (2 * x)],
                          [(-2 + y**2 - np.cos(omega * t)) / (2 * y)]])
    # compute 3rd matrix of equation
    equ2 = np.matrix([[(np.sin(t) / (2 * x))],
                           [(omega * np.sin(omega * t)) / (2 * y)]])
    # bring all matrices together
    df = scaling * equ1 - equ2
    # return the correct portions of matrices
    qout = np.array([df[0, 0], df[1,0]])
    
    
    return qout 

def RK3_step(f,t,q,dt,options):
    """
    
    Runge-Kutta 3 algorithm. Takes in value of q uses 3 steps of k and weights
    using simpsons rule. Returns q at time step t+dt
    
    Parameters
    ----------
    f: function
        function the method is solving
    t : integer 
        time at timestep
    q : 2 X 1 array
        x and y values of previous timestep 
    dt : integer 
        timestep
    options: 3X1 array
        gives either stiff or non stiff values for gamma , epsilon and omega
    
        
    Returns
    -------
    q: 2 X 1 array
        x and y values at  t + dt
    """
    assert type(t) == np.float or  np.float64, 't is not type float'
    
    assert type(options) == np.ndarray or list, 'options is not of type nd.array'
    assert len(options) == 3,'3 values not given for options'
    assert type(dt) == np.float or  np.float64, 'dt is not type float'
    assert hasattr(f,'__call__'),'f must be a callable function'
    # compute k1. 1st step
    k1=f(t,q,options)
    # compute k2. 2nd step
    k2=f(t+(0.5*dt),q+(dt/2)*k1,options)
    # compute k3. 3rd step
    k3=f(t+dt,q+dt*(-k1+(2*k2)),options)
    # compute value of q given weighting of k1 , k2 and k3 
    q=q+(dt/6)*(k1+(4*k2)+k3)
    
    return q

    
def MyGRRK3_step(f,t,q,dt,options):
    """
    
    Gauss-Radau Runge-Kutta method. Take in value of q and iterates to q at 
    timestep t+dt. Method is implicit so use fsolve to find values of k1 and k2
    simutaneously.
    
    Parameters
    ----------
    f: function
        function the method is solving
    t : integer 
        time at timestep
    q : 2 X 1 array
        x and y values of previous timestep 
    dt : integer
        timestep
    options: 3X1 array
        gives either stiff or non stiff values for gamma , epsilon and omega
    
        
    Returns
    -------
    q: 2 X 1 array
        x and y values at  t + dt
    """
    assert type(t) == np.float or  np.float64, 't is not type float'
    assert type(q) == np.ndarray  
    assert type(options) == np.ndarray or list, 'options is not of type nd.array'
    assert len(options) == 3,'3 values not given for options'
    assert type(dt) == np.float or  np.float64, 'dt is not type float'
    assert hasattr(f,'__call__'),'f must be a callable function'
    # initial guess for K1
    K1 = f(t+(dt/3),q,options)
    # initiial guess for k2
    K2 = f(t+dt,q,options)
    # bring together
    K = np.array([K1, K2])
    
    assert len(K) != 0,'K has no values'
    # use if statements to reshape for different size problems
    if len(q) >= 2:
        m,n=K.shape
        # reshape to 4,0 array
        K=np.reshape(K,newshape=(m*n,))
    

    def F(K):
        """
        
        Function for use in fsolve
        
        Parameters
        ----------
        K : 2X2 array
            array of gueses for k1 and k2
        
            
        Returns
        -------
        F2: 4 X 0 array
            values of k1 and k2
        """
        assert type(K)==np.ndarray,'K is not an array'
        len_K = len(K)
        assert len(K) != 0,'K has no values'
        
        # set array to find k1 and k2
        F1 = np.array([[f(t+(dt/3),q+((dt/12)*((5*K[0:len(K)//2])-K[len(K)//2:])),options)],[f(t+dt,q+((dt/4)*((3*K[0:len(K)//2])+K[len(K)//2:])),options)]])
        
        if len(K) >= 4:
            m,l,n = F1.shape
            # define function to find the zero 
            F2 = K - np.reshape(F1,newshape=(m*n,))
        if len(K) < 4:
            # define function to find the zero 
            F2 = K - np.reshape(F1,newshape=(len(K)))
       
        return F2
    
    len_K = len(K)
    # use fsolve to find the values of k such that F = 0  
    kout= scipy.optimize.fsolve(F,K,full_output = False)
    #supress warnings
    warnings.filterwarnings('ignore', 'The iteration is not making good progress')
    # set k1 to first 2 values of kout
    k1=kout[0:len(K)//2]
    
    k2=kout[len(K)//2:]
    
    # calculate q n+1 with weighting 
    q=q+((dt/4)*((3*k1)+k2))
    
 
    return q


def task3():
    """
    
    Task 3 function. Returns plots of exact solution and computed solution 
    from RK3 and GRRK3 methods for non-stiff options
    
    Parameters
    ----------
    none
    
        
    Returns
    -------
    subplot 1: exact solution, GRRK and RK of non-stiff case. x-component
    subplot 2: exact solution, GRRK and RK of non-stiff case. y-component
    """
    t_end=1
    dt=0.05
    # set number of steps
    steps=round(t_end/dt)   
    # initialise data
    q_series_exp= np.zeros((steps+1,2))
    
    q_series_imp= np.zeros((steps+1,2))
    
    t_series=np.zeros((steps+1))
    # set gamma, epsilon, omega values
    options=[-2,0.05,5] # not-stiff
    
    # initial data
    q_series_exp[0,0]= np.sqrt(2)
    q_series_exp[0,1]= np.sqrt(3)
    q_series_imp[0,0]= np.sqrt(2)
    q_series_imp[0,1]= np.sqrt(3)
    t_series[0]=0
    t_start = 0
    #set time steps
    xs = np.linspace(t_start, t_end, steps+1)
    # calculate exact values
    exact_x = np.sqrt(1+np.cos(xs))
    exact_y = np.sqrt(2+np.cos(options[2]*xs))
    
    
    for step in range(1,len(xs)):
        t = xs[step]
        # set values of q
        q_exp=q_series_exp[step-1,:]
        q_imp=q_series_imp[step-1,:]
        # calculate values of q from methods
        q_series_exp[step,:]= RK3_step(f,t-dt,q_exp,dt,options)

        q_series_imp[step,:]= MyGRRK3_step(f,t-dt,q_imp,dt,options)
        # set time series for plotting
        t_series[step] = t    
        
    # plots
    fig = pyplot.figure(1)
    ax = fig.add_subplot(121)
    ax.plot(t_series,q_series_exp[:,0], '.', label='RK')
    ax.plot(t_series,q_series_imp[:,0], '.', label='GRRK')
    ax.plot(t_series,exact_x, label='Exact')
    pyplot.title('x-solution with RK and GRRK \n method approximations')
    pyplot.xlabel('t')
    pyplot.ylabel('x')
    
    
    pyplot.legend()
    
    ax2 = fig.add_subplot(122)
    ax2.plot(t_series,q_series_exp[:,1], '.', label='RK')
    ax2.plot(t_series,q_series_imp[:,1], '.', label='GRRK')
    ax2.plot(t_series,exact_y, label='Exact')
    pyplot.xlabel('t')
    pyplot.ylabel('y')
    pyplot.title('y-solution with RK and GRRK \n methods approximations')
    fig.subplots_adjust(wspace=.5)
    pyplot.legend()

    return 
    
def task4():
    """
    
    Task 4 function. Returns plot of error with timestep. Shows computed points 
    of error with given timestep for each method using non- stiff options with
    polyfit to show convergence rate of 3
    
    Parameters
    ----------
    none
    
        
    Returns
    -------
    plot: error against timestep of RK and GRRK method with lines of best fit
    for each
    """
    # initialise data for error
    err_RK= np.zeros(8)
    err_GRRK= np.zeros(8)
    t_end = 1
    dt_series = np.zeros((8))
    
    for j in range(8):
        # set for different values of dt
        dt= 0.1/(2**j)
        steps=round(t_end/dt)
        # set time values
        xs = np.linspace(0, 1, steps+1)
        # initialise data
        q_series_exp= np.zeros((steps+1,2))        
        q_series_imp= np.zeros((steps+1,2))
        t_series=np.zeros((steps+1))
        #options=[-2e5,0.5,20] # stiff
        options=[-2,0.05,5] # not-stiff
        
        # set intial data
        q_series_exp[0,0]= np.sqrt(2)
        q_series_exp[0,1]= np.sqrt(3)
        q_series_imp[0,0]= np.sqrt(2)
        q_series_imp[0,1]= np.sqrt(3)
        t_series[0]=0
        
        t = 0
        # set exact values
        exact_y = np.sqrt(2+np.cos(options[2]*xs))
        dt_series[j] = dt

        for step in range(1,len(xs)):
            t = xs[step]
            
            # set q to previous timestep
            q_exp=q_series_exp[step-1,:]
            q_imp=q_series_imp[step-1,:]

            # calculate q using methods
            q_series_exp[step,:]= RK3_step(f,t-dt,q_exp,dt,options)   
            q_series_imp[step,:]= MyGRRK3_step(f,t-dt,q_imp,dt,options)
      
            t_series[step] = t
            
            # calculate error and store for different dt
            err_RK[j] = dt*np.sum(np.abs(q_series_exp[:,1]- exact_y))
            assert type(err_RK[j])==np.float or np.float64,'error is not of type np.float'
    
            err_GRRK[j] = dt*np.sum(np.abs(q_series_imp[:,1]- exact_y))
            assert type(err_GRRK[j])==np.float or np.float64,'error is not of type np.float'
    assert len(err_RK)==8, 'There are not enough error values'
    assert len(err_GRRK)==8, 'There are not enough error values'
    
    # use polyfit to calculate gradient and y intercept of error convergence
    m_exp ,p_exp= np.polyfit(np.log(dt_series),np.log(err_RK),1)
    m_imp ,p_imp= np.polyfit(np.log(dt_series),np.log(err_GRRK),1)
    
    # plots
    fig = pyplot.figure(2)
    pyplot.loglog(dt_series,err_RK,'r+',label = 'Error of RK method')
    pyplot.loglog(dt_series,err_GRRK,'k*',label = 'Error of GRRK method')
    line_fit_x = dt_series
    line_fit_y_exp = np.exp(p_exp)*dt_series**m_exp
    line_fit_y_imp = np.exp(p_imp)*dt_series**m_imp
    pyplot.loglog(line_fit_x,line_fit_y_exp,label = 'RK: Gradient = %.2f' %m_exp)
    pyplot.loglog(line_fit_x,line_fit_y_imp,label = 'GRRK:Gradient = %.2f' %m_imp)
    pyplot.xlabel(r'$\Delta t$')
    pyplot.ylabel('Error')
    pyplot.title('Error convergence of RK3 and GRRK3 for non-stiff case')
    pyplot.legend()
    
    
def task5():
    """
    
    Task 5 function. Returns plots of exact solution and computed solution 
    from RK3 method for stiff options
    
    Parameters
    ----------
    none
    
        
    Returns
    -------
    subplot 1: exact solution and RK of stiff case. x-component. Unstable for RK
    subplot 2: exact solution and RK of stiff case. y-component. Unstable for RK
    """
    t_end=1
    dt=0.001
    
    steps=round(t_end/dt)
    
    # initialise arrays
    q_series_exp = np.zeros((steps+1,2))
    
    t_series=np.zeros((steps+1))
    # set options to stiff
    options=[-2e5,0.5,20] # stiff
    
    # initial data
    q_series_exp[0,0]= np.sqrt(2)
    q_series_exp[0,1]= np.sqrt(3)
    t_series[0]=0
    t = 0
    
    xs = np.linspace(0, 1, steps+1)
    #exact solutions
    exact_x = np.sqrt(1+np.cos(xs))
    exact_y = np.sqrt(2+np.cos(options[2]*xs))
    
    
 
    for step in range(1,len(xs)):
        t = xs[step]
        
        # set q to previous timestep
        q_exp=q_series_exp[step-1,:]
        # calculate q at next time step via RK3 method        
        q_series_exp[step,:]= RK3_step(f,t-dt,q_exp,dt,options=[-2e5,0.5,20])
        
          
        t_series[step] = t
             
    #plots   
    fig = pyplot.figure(3)
    ax = fig.add_subplot(121)
    ax.plot(t_series,q_series_exp[:,0],  label='RK')
    ax.plot(t_series,exact_x, label='exact')
    pyplot.xlabel('t')
    pyplot.ylabel('x')
    ax.set_ylim([1.24,1.5])
    fig.subplots_adjust(wspace=.5)
    pyplot.title('x-solution with RK method \n approximations')
    pyplot.legend()
    
    ax2 = fig.add_subplot(122)
    ax2.plot(t_series,q_series_exp[:,1],  label='RK')
    ax2.plot(t_series,exact_y, label='exact')
    ax2.set_ylim([1,1.8])
    pyplot.xlabel('t')
    pyplot.ylabel('y')
    pyplot.title('y-solution with RK method \n approximations')
    pyplot.legend()

    
def task6():
    """
    
    Task 6 function. Returns plots of exact solution and computed solution 
    from GRRK3 method for stiff options
    
    Parameters
    ----------
    none
    
        
    Returns
    -------
    subplot 1: exact solution and GRRK of stiff case. x-component. stable for GRRK
    subplot 2: exact solution and GRRK of stiff case. y-component. stable for GRRK
    """
    t_end=1
    dt=0.005
    
    steps=round(t_end/dt)
    
    
    q_series_imp= np.zeros((steps+1,2))
    
    t_series=np.zeros((steps+1))
    # set options to stiff values
    options=[-2e5,0.5,20] # stiff
    #initial data
    q_series_imp[0,0]= np.sqrt(2)
    q_series_imp[0,1]= np.sqrt(3)
    t_series[0]=0
    t = 0
    
    xs = np.linspace(0, 1, steps+1)
    #exact solutions
    exact_x = np.sqrt(1+np.cos(xs))
    exact_y = np.sqrt(2+np.cos(options[2]*xs))
    
    
 
    for step in range(1,len(xs)):
        t = xs[step]
        
        # set q to value at previous timestep
        q_imp=q_series_imp[step-1,:]
        # calculate q from GRRK3 method
        q_series_imp[step,:]= MyGRRK3_step(f,t-dt,q_imp,dt,options)
              
        t_series[step] = t
        
    
    # plots
    fig = pyplot.figure(4)
    ax = fig.add_subplot(121)
    ax.plot(t_series,q_series_imp[:,0], '.', label='GRRK')
    ax.plot(t_series,exact_x, label='Exact')
    pyplot.legend()
    pyplot.xlabel('t')
    pyplot.ylabel('x')
    pyplot.title('x-solution with GRRK method \n approximations')
    
    
    ax2 = fig.add_subplot(122)
    ax2.plot(t_series,q_series_imp[:,1], '.', label='GRRK')
    ax2.plot(t_series,exact_y, label='Exact')
    pyplot.xlabel('t')
    pyplot.ylabel('y')
    fig.subplots_adjust(wspace=.5)
    pyplot.title('y-solution with GRRK method \n approximations')
    pyplot.legend()
            
            
            
def task7():
    """
    
    Task 7 function. Returns plot of error with timestep. Shows computed points 
    of error with given timestep for GRRK3 with stiff options method with 
    polyfit to show convergence rate of 3
    
    Parameters
    ----------
    none
    
        
    Returns
    -------
    plot: error against timestep of  GRRK method with lines of best fit
    
    """
    # intitial error data
    err= np.zeros(8)
    t_end = 1
    dt_series = np.zeros((8))
    for j in range(8):
        # set dt to various values
        dt= 0.05/(2**j)
        steps=round(t_end/dt)
        xs = np.linspace(0, 1, steps+1)
 
        # intialise arrays
        q_series_exp= np.zeros((steps+1,2))        
        q_series_imp= np.zeros((steps+1,2))
        t_series=np.zeros((steps+1))
        # set options to stiff values
        options=[-2e5,0.5,20] # stiff
        #options=[-2,0.05,5] # not-stiff
        
        # initial data
        q_series_exp[0,0]= np.sqrt(2)
        q_series_exp[0,1]= np.sqrt(3)
        q_series_imp[0,0]= np.sqrt(2)
        q_series_imp[0,1]= np.sqrt(3)
        t_series[0]=0
        
        t = 0
        exact_y = np.sqrt(2+np.cos(options[2]*xs))
        dt_series[j] = dt
        
        
     
        for step in range(1,len(xs)):
            t = xs[step]
            
    
            # set q to value at previous timestep
            q_imp=q_series_imp[step-1,:]
                        
            # calculate q from GRRK3 method
            q_series_imp[step,:]= MyGRRK3_step(f,t-dt,q_imp,dt,options)
      
            t_series[step] = t
                     
            
            # calculate errors
            err[j] = dt*np.sum(np.abs(q_series_imp[:,1]- exact_y))
            assert type(err[j])==np.float or np.float64,'error is not of type np.float'
    assert len(err)==8, 'There are not enough error values'
    
    #use polyfit to obtain gradient and y intercept of error convergence
    m ,p= np.polyfit(np.log(dt_series),np.log(err),1)
    
    # plots
    fig = pyplot.figure(5)
    pyplot.loglog(dt_series,err,'r+',label = 'Computed Value')
    line_fit_x = dt_series
    line_fit_y = np.exp(p)*dt_series**m
    pyplot.loglog(line_fit_x,line_fit_y,label = 'Gradient = %.2f' %m)
    pyplot.xlabel(r'$\Delta t$')
    pyplot.ylabel('Error')
    pyplot.title('Error convergence of GRRK3 for stiff case')
    pyplot.legend()
            


### Tests ###

import pytest
# test if the function gets close to the exact value
def test_RK_function_close():
    q_test = RK3_step(f,0.01,np.array([np.sqrt(2),np.sqrt(3)]),0.01,options=[-2,0.05,5])
    tol = 1e-4  
    assert np.abs(q_test[0]) <= np.abs(np.sqrt(1+np.cos(0.01)) + tol)
    assert np.abs(q_test[1]) <= np.abs(np.sqrt(2+np.cos(5*0.01)) + tol)
# test if the algorithm gives correct values for known function
def test_known_function():
    f_known = lambda x,q,options:np.array([x**2,x])
    q_known_RK = RK3_step(f_known,0,np.array([0,0]),1,[0,0,0]) 
    q_known_GRRK = MyGRRK3_step(f_known,0,np.array([0,0]),1,[0,0,0])
    tol = 1e-4  
    assert q_known_RK[0] <= 1/3 + tol 
    assert q_known_RK[1] <= 1/2 + tol
    assert q_known_GRRK[0] <= 1/3 + tol 
    assert q_known_GRRK[1] <= 1/2 + tol
    
# test if algorithm can handle a 1 value problem
def test_1_value_function():
    f_known = lambda x,q,options:x**2
    q_known_RK = RK3_step(f_known,0,np.array([0]),1,[0,0,0])
    q_known_GRRK = MyGRRK3_step(f_known,0,np.array([0]),1,[0,0,0])
    tol = 1e-4
    assert q_known_RK <= 1/3 + tol     
    assert q_known_GRRK <= 1/3 + tol 
# test if algorithm can handle a 3 value problem  
def test_3_value_function():
    f_known = lambda x,q,options:np.array([x**2,x,x**3])
    q_known_RK = RK3_step(f_known,0,np.array([0,0,0]),1,[0,0,0]) 
    q_known_GRRK = MyGRRK3_step(f_known,0,np.array([0,0,0]),1,[0,0,0])
    tol = 1e-1
    

    assert q_known_RK[0] <= 1/3 + tol 
    assert q_known_RK[1] <= 1/2 + tol
    assert q_known_RK[2] <= 1/4 + tol
    assert q_known_GRRK[0] <= 1/3 + tol 
    assert q_known_GRRK[1] <= 1/2 + tol
    assert q_known_GRRK[2] <= 1/4 + tol
# test if algorithm can handle a 4 value problem
def test_4_value_function():
    f_known = lambda x,q,options:np.array([x**2,x,x**3,x**4])
    q_known_RK = RK3_step(f_known,0,np.array([0,0,0,0]),1,[0,0,0]) 
    q_known_GRRK = MyGRRK3_step(f_known,0,np.array([0,0,0,0]),1,[0,0,0])   
    tol = 1e-1
    assert q_known_RK[0] <= 1/3 + tol 
    assert q_known_RK[1] <= 1/2 + tol
    assert q_known_RK[2] <= 1/4 + tol
    assert q_known_RK[2] <= 1/5 + tol
    assert q_known_GRRK[0] <= 1/3 + tol 
    assert q_known_GRRK[1] <= 1/2 + tol
    assert q_known_GRRK[2] <= 1/4 + tol
    assert q_known_GRRK[3] <= 1/5 + tol
    
    

if __name__ == '__main__':    
    test_RK_function_close()
    test_known_function()
    test_1_value_function()
    test_3_value_function()
    test_4_value_function()
    task3()
    task4()
    task5()
    task6()
    task7()
    
    