# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 14:55:26 2016

@author: Nicholas
"""
import scipy.optimize
import numpy as np
from matplotlib import pyplot as plt
import sympy
sympy.init_printing()

r_eq=1.1283
Do=258.9
alpha=2.302
k=2743
def harmonic(r):
    V=0.5*k*(r-r_eq)**2
    return V

def kratzer(r):
    V=Do*((r-r_eq)/r)**2
    return V

def Morse(r):
    V=Do*(1-np.exp(-alpha*(r-r_eq)))**2
    return V



def mol_dyn_l2(f,r_guess):

    min_V=scipy.optimize.minimize(f,r_guess)
    return min_V

r=np.linspace(1,1.5)

plt.figure(1)
plt.axis([1, 1.5, 0, 200])
plt.plot(r,harmonic(r),'b')


plt.plot(r,kratzer(r),'g')



plt.plot(r,Morse(r),'r')

plt.legend(('Harmonic','Kratzer','Morse'))


def symp_Morse():


    del_t=0.001
    steps=1/0.001
    r_eq=1.1283
    D_0=258.9
    alpha=2.302




    return Do*(1-sympy.exp(-alpha*(r-r_eq)))**2

def diff_symp_Morse():
    v=symp_Morse()
    return sympy.diff(v)

r=sympy.symbols('r')















def mol_dyn2():


    #Rc=L/2
    #a=np.zeros_like(n)
    x=np.array([[0,0,0],[1.12,0]])
    #rij=np.zeros_like(n)
    N=len(x)




    for i in range(N):

        rij=np.sqrt(np.dot(x[0],x[1]))
        #rij[np.abs(rij)>L/2]-=np.sign(rij[np.abs(rij)>L/2])*L




        for j in range(i+1,N):
            #mag=np.dot(rij[j,:],rij[j,:])


            d_phi_dr=sympy.lambdify(r,symp_Morse(),'numpy')
            d_phi = d_phi_dr(rij)

            a[i,:]+=rij[j,:]*d_phi
            a[j,:]-=rij[j,:]*d_phi

    return a


def velo_vel():
    x=np.array([[0,0],[1.12,0]])
    del_t=0.001
    steps=int(1/del_t)
    a=np.zeros_like(x)
    v=np.zeros_like(x)
    for step in range(steps):

        x=x + del_t*v +0.5*del_t**2*a


        #x=peri_bound(x,L)


        v_star=v+0.5*del_t*a

        a=mol_dyn2(x,N,L)


        v=v_star+0.5*del_t*a
        x_plot[step,:,:]=x
        T=temp(v,N,L)
        T_plot.append(T)
        time[step]=step*del_t


