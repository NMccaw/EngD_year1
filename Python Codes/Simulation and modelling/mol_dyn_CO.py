# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 17:47:00 2016

@author: Nicholas
"""


import sympy
import numpy as np
from matplotlib import pyplot as plt


r_eq=1.1283
D_0=258.9
alpha=2.302

#r,r_eq,D_0,alpha=sympy.symbols('r,r_eq,D_0,alpha')
r=sympy.symbols('r')

v=D_0*(1-sympy.exp(-alpha*(r-r_eq)))**2
v


dv=sympy.diff(v)
dv

(D_0,r)=sympy.symbols('D_0,r')

dv_dt=sympy.lambdify(r,dv,'numpy')








    #x=np.array([[0,0,0],[2*(1/6),0,0]])
    #rij=np.zeros_like(n)



def accel(x):
    N=len(x)
    a=np.zeros_like(x)
    for i in range(N):

        rij=x[i,:]-x[:,:]
        #rij[np.abs(rij)>L/2]-=np.sign(rij[np.abs(rij)>L/2])*L




        for j in range(i+1,N):

            mag=np.sqrt(np.dot(rij[j,:],rij[j,:]))


            d_phi_dr=dv_dt(mag)

            a[i,:]-=rij[j,:]*d_phi_dr/1
            a[j,:]+=rij[j,:]*d_phi_dr/1.332029




    return a

def velo_vel():
    x=np.array([[1.0],[2.11]])

    del_t=0.0001
    steps=10000#int(1/del_t)
    a=accel(x)
    v=np.zeros_like(a)
    x_plot1=np.zeros((steps))
    x_plot2=np.zeros((steps))
    for step in range(steps):



        x=x + del_t*v +0.5*del_t**2*a


        #x=peri_bound(x,L)


        v_star=v+0.5*del_t*a

        a=accel(x)


        v=v_star+0.5*del_t*a
        x_plot1[step]=x[0,:]
        x_plot2[step]=x[1,:]


    plt.plot(x_plot1)
    plt.plot(x_plot2)


    return x[0]



velo_vel()