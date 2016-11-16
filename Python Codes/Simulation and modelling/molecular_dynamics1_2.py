# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:42:45 2016

@author: Nicholas
"""


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def d_phi(r,N,Rc):


    d_phi=24*(((2*((1/r)**7))-((1/r)**4)))


    return d_phi

def r_ij(n):
    rij=particle_i-n
    for i in range(N):
        rij[i][abs(rij[i])>L/2] -= np.sign(rij[i][abs(rij[i])>L/2])*L
    mag=np.linalg.norm(r,axis=0)
    print(mag)
    return mag,rij



def mol_dyn():
    L=5
    N=2
    Rc=L/2
    #n= np.random.random(size=(N,3))
    n=np.array([[0,0,0],[2**(1/6),0,0]])
    I=np.random.randint(N)
    particle_i=n[I,:]

    a=np.zeros((N,3))

    #rij=np.zeros[N,3]

    for i in range(N):

        rij=n[i,:]-n[:,:]

        rij[np.abs(rij)>L/2]-=np.sign(rij[np.abs(rij)>L/2])*L





        for j in range(i+1,N):
            mag=np.dot(rij[j,:],rij[j,:])


            d_phi_dr=d_phi(mag,N,Rc)

            a[i,:]+=rij[j,:]*d_phi_dr
            a[j,:]-=rij[j,:]*d_phi_dr

    return a




def mol_dyn2(n,N,L):

    N=len(n)
    Rc=L/2
    a=np.zeros_like(n)
    #x=np.array([[0,0,0],[2*(1/6),0,0]])
    #rij=np.zeros_like(n)




    for i in range(N):

        rij=n[i,:]-n[:,:]
        rij[np.abs(rij)>L/2]-=np.sign(rij[np.abs(rij)>L/2])*L




        for j in range(i+1,N):
            mag=np.dot(rij[j,:],rij[j,:])


            d_phi_dr=d_phi(mag,N,Rc)

            a[i,:]+=rij[j,:]*d_phi_dr
            a[j,:]-=rij[j,:]*d_phi_dr

    return a

def peri_bound(x,L):
    x[x<0]+=L
    x[x>L]-=L
    return x


def time_evo():
    del_t=0.005
    steps=100
    L=6.1984
    Rc=0.5*L

    #x=np.zeros((t/del_t,2,3))


    #x=np.array([[0,0,0],[2**(1/3),0,0]])
    x=np.loadtxt('input.dat',dtype=float)





    #x=np.array([[0,0,0],[2**(1/6),0,0]])
    N=len(x)
    x_plot=np.zeros((steps,N,3))

    mag=np.zeros(N)
    a=np.zeros((N,3))
    v=np.zeros((N,3))
    time=np.zeros(steps)
    T_plot=[]

    for step in range(steps):
        x=x + del_t*v +0.5*del_t**2*a


        x=peri_bound(x,L)


        v_star=v+0.5*del_t*a

        a=mol_dyn2(x,N,L)


        v=v_star+0.5*del_t*a
        x_plot[step,:,:]=x
        T=temp(v,N,L)
        T_plot.append(T)
        time[step]=step*del_t


    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs=x_plot[:,100,0],ys=x_plot[:,100,1],zs=x_plot[:,100,2])

    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d')
    for i in range(N):
        ax.plot(xs=x_plot[:,i,0],ys=x_plot[:,i,1],zs=x_plot[:,i,2])
    print(T.shape)
    fig=plt.figure(3)

    plt.plot(time, T_plot)



    return

def temp(v,N,L):
    for i in range(N):
        E_kin=0.5*(L**2)*np.dot(v[i,:],v[i,:])
        T=2/(3*N)*np.sum(E_kin)
    return T






