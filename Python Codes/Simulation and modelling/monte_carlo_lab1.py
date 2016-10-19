# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:25:14 2016

@author: Nicholas
"""

"""
import numpy as np

upp_bound=10
lower_bound=0

xn=upp_boudx*np.random.random - lower_bound

print(xn)

"""
import numpy as np
import pylab
from scipy import special
from matplotlib import pyplot
unit=2
def f2(x):
    return np.sqrt((unit**2)-(x**2))



def mc_integrate(f2, domain_x, domain_y, N = 10000):
    """
    Monte Carlo integration function: to be completed.

    computes the integral by scattering random points across the domain.
    Returns the integral by multipying the area of the domain by the points
    under curve divided by total points

    Result, for the given f, should be around 1.46.
    """
    import numpy.random
    ran_num= numpy.random.random((N,2))


    ran_x = ((domain_x[1]-domain_x[0])* ran_num[:,0]) + domain_x[0]

    ran_y = ((domain_y[1]-domain_y[0])* ran_num[:,1]) + domain_y[0]


    integrand= f2(ran_x)

    p=sum(ran_y < integrand)/N


    area=(domain_x[1]-domain_x[0])*(domain_y[1]-domain_y[0])

    I=p*area


    return I

def err_N():
    err=[]
    Ni=[]
    for i in range(17):
        N=2**i
        Ni.append(N)

        err_v=abs(mc_integrate(f2, [0,2], [0,2],N)-np.pi)
        err.append(err_v)

    p=pylab.polyfit(np.log(Ni),np.log(err),1)

    line=np.log(Ni)*p[0]+(p[1])
    print(p)
    pylab.plot((Ni),np.exp(line),'b')
    pylab.loglog()
    pylab.plot((Ni),(err),'ro')
    pylab.loglog()
    pylab.show


def mv_integrate(f, domain_x, N = 10000):
    """
    Mean value Monte Carlo integration: to be completed

    """
    import numpy.random
    ran_num= numpy.random.random((N,1))
    sum_fx=0

    ran_x = ((domain_x[1]-domain_x[0])* ran_num[:,0]) + domain_x[0]
    for i in range(N):
        sum_fx=sum_fx+f(ran_x[i])
    I=((domain_x[1]-domain_x[0])/N)*sum_fx
    return I


def volume_hypersphere(ndim=3):
    return np.pi**(float(ndim)/2.0) / special.gamma(float(ndim)/2.0 + 1.0)


def f3(x):
    return np.sqrt(sum(x**2))

def mc_integrate_multid(f, domain, N = 10000):
    """
    Monte Carlo integration in arbitrary dimensions
    (read from the size of the domain): computes the integreal of a multi-
    dimensional object by monte carlo integration.
    """
    import numpy.random
    ran_num= numpy.random.random((N,len(domain)+1))
    ran_x=[]
    for i,d in enumerate(domain):



        ran_x.append(((d[1]-d[0])* ran_num[:,i]) + d[0])


    ran_y = ((1+1)* ran_num[:,-1]) + -1
    print(ran_x[0])
    print(len(d))

    integrand=f3(np.array(ran_x))
    print(integrand)
    p=sum(integrand <= 1)/N
    print(p)
    vol=(d[1]-d[0])**3
    print(vol)
    print(d[0])
    I=vol*p

    return I

def f4(x):
    return 1/(np.exp(x)+1)

def f5(x):
    return (np.power(x,-0.5))/(np.exp(x)+1)

def imp_samp(f,domain_x, N = 10000):
    """
    Importance sampling:
    an integral has a power distibution function. this is used as a weighting
    function such that random points are generated with a greater focus within
    this dustrubution.
    """


    import numpy.random
    ran_num= numpy.random.power(0.5,N)
    print(ran_num)
    p=f(ran_num)
    print(p)
    I=(2/N)*sum(p)
    return I

def q(x):
    return 1/(3*x**(2/3))


def p(x):
    return 1/(2*np.sqrt(x))

def Kq(N=10000):
    """
    Returns the values of x which are kept after the rejection sampling
    random numbers are distributed along x from the known function Kq(x).
    Random numbers are then generated along a uniform distribution in y. If the
    uniform distribution is less than the Kq then the random uniform
    distribution values are kept.
    """
    x=np.linspace(0,1,N)

    x_keep=[]
    y_keep=[]
    import numpy.random
    k=1.6
    ran_num=numpy.random.power(1/3,N)
    p_x=p(ran_num)

    integ=q(ran_num)
    ran_num_U=numpy.random.random(N)
    ran_num_u_y=1*ran_num_U


    prob= sum((p_x/(k*integ))< ran_num_u_y)
    for i in range(N):
        #if (p_x[i]/(k*integ[i])) > ran_num_u_y[i]:
        if (p_x[i]/(k*integ[i]))  > ran_num_u_y[i]:
            x_keep.append(ran_num[i])
            y_keep.append(ran_num_u_y[i])

    print(len(x_keep))

    #pyplot.semilogy(x_keep,y_keep,'ro')

    pyplot.semilogy(x,p(x))
    pyplot.hist(x_keep,bins=50,normed=True)
    proba=sum(f5(x_keep)< y_keep )/(len(x_keep))
    I=proba

    return I

def imp_samp2(f,domain_x, N = 10000):


    import numpy.random
    ran_num= Kq()

    p=f(ran_num)

    I=(1/N)*sum(p)
    return I