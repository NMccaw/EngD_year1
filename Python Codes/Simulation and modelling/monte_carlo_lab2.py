# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 13:33:53 2016

@author: Nicholas
"""

def U_cal(r,L):
    for i in range(N):
        for j in range(i):
            if r[i,j] < 0.5*L:
                u[i,j]=4*((1/(r[i,j]**12))-(1/(r[i,j]**6)))
            else:
                u[i,j]=0

    return u


def r_mag(n,part_state):


    for i in range(N):
        for j in range(i):
            r[i,j]= n[i]-n[j]



    #r=(part_state-n)




            r[i,j][abs(r[i,j])>L/2] -= np.sign(r[i,j][abs(r[i,j])>L/2])*L

    mag=np.linalg.norm(r,axis=2)



        #r_ij[i,j]=mag[i,j][mag[i,j]>0]



    return mag





N=100
steps=5000
a=1
rho=a/10
L=10
    #(rho/N)**(-1/3)
T=2
#n=np.random.random(size=(N,3))*L
#u=np.zeros(N)
beta=0.1
#states=np.random.rand(N,3)*L
#I=np.random.randint(N)
#part_state=n[I,:]

#mag=np.zeros((N,N))

E_keep=np.zeros(steps)
I=np.random.randint(N)
n=np.random.random(size=(N,3))*L
part_state=n[I,:]
#r_ij_w=np.zeros(N-1)

u=np.zeros((N,N))
r=np.zeros((N,N,3))
r_ij=np.zeros((N,N,3))

for step in range(steps):

    r_ij=r_mag(n,part_state)

    u=U_cal(r_ij,L)

    E=np.sum(u)


    part_state_new=np.random.rand(3)*L



    #for i in range(0,I):
    #*r_ij**2)
    #for i in range(I+1,N):
    #    delE_2 = sum(u*r_ij[i]**2)
    n_new=n.copy()

    n_new[I,:]=np.random.random(size=(3))*L

    r_ij_new=r_mag(n_new,part_state_new)


    u_new=U_cal(r_ij_new,L)

    E_new=np.sum(u_new)
    delE=E_new-E
    #delE=delE_1+delE_2

    if delE <0 or np.exp(-beta*(delE))>np.random.random():
            #part_state += change
        r_ij=r_mag(n_new,part_state_new)
        u_new=U_cal(r_ij,L)
        u=u_new
            #for i in range(0,I):
        E=np.sum(u)
            #for i in range(I+1,N):
            #    delE_2 = sum(u*r_ij[i]**2)
            #delE=delE_1+delE_2
    E_keep[step]=E



print(r_ij)
print(u)
print(E_keep[-1])

import matplotlib.pyplot as plt
plt.plot(E_keep)