# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:42:45 2016

@author: Nicholas
"""


import numpy as np




def d_phi(r):
    for i in range(N):
        if r < Rc:
            d_phi=24*(2*(((1/r)**14)-((1/r)**8)))
        else:
            d_phi=0
    return d_phi

def r_ij(n):
    rij=particle_i-n
    for i in range(N):
        rij[i][abs(rij[i])>L/2] -= np.sign(rij[i][abs(rij[i])>L/2])*L
    mag=np.linalg.norm(r,axis=0)
    print(mag)
    return mag,rij



def mol_dyn():

    N=2
    Rc=2.5
    #n= np.random.random(size=(N,3))
    n=np.array([[0,0,0],[2*(1/6),0,0]])
    I=np.random.randint(N)
    particle_i=n[I,:]
    L=1
    a=np.zeros((N,3))
    mag=np.zeros(N)
    #rij=np.zeros[N,3]

    for i in range(N):

        rij=n[i,:]-n[:,:]
        rij[np.abs(rij)>L/2]-=np.sign(rij[np.abs(rij)>L/2])*L




        for j in range(i+1,N):
            mag=np.dot(rij[j,:],rij[j,:])


            d_phi_dr=d_phi(mag)

            a[i,:]+=rij[j,:]*d_phi_dr
            a[j,:]-=rij[j,:]*d_phi_dr

    return a




def mol_dyn2(n):

    N=len(n)
    Rc=2.5
    I=np.random.randint(N)
    #particle_i=n[I,:]
    L=1
    rij=np.zeros_like(n)



    for i in range(N):

        rij=n[i,:]-n[:,:]
        rij[np.abs(rij)>L/2]-=np.sign(rij[np.abs(rij)>L/2])*L




        for j in range(i+1,N):
            mag=np.dot(rij[j,:],rij[j,:])


            d_phi_dr=d_phi(mag)

            a[i,:]+=rij[j,:]*d_phi_dr
            a[j,:]-=rij[j,:]*d_phi_dr

    return a
def peri_bound(x,L):
    x[x<0]+=L
    x[x>L]-=L
    return x


def time_evo():
    del_t=0.01
    t=1
    x=np.zeros((t/del_t,2,3))
    print(x[1,:,0])

    x[0,0,:]=np.array([0,0,0])
    x[0,1,:]=np.array([2*(1/6),0,0])
    print(x[:,1,0])

    #x=np.array([[0,0,0],[2*(1/6),0,0],[0,0,0])
    N=len(x)

    mag=np.zeros(N)
    a=np.zeros((N,3,t/del_t))
    v=np.zeros((N,3,t/del_t))

    for n in range(t):
        x[n+1,:,:]=x[n,:,:] +del_t*v[n,:,:] +0.5*del_t**2*a[n,:,:]

        v_star=v[n,:,:]+0.5*del_t*a[n,:,:]

        a[n+1,:,:]=mol_dyn2(x[:,:,n+1])
        v[n+1,:,:]=v_star+0.5*del_t*a[n+1,:,:]

    return a,v,x




