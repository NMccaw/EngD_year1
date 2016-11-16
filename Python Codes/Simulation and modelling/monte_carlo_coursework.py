# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:30:48 2016

@author: Nicholas
"""

import numpy as np
import matplotlib.pyplot as plt





def Ising():
    N = 40
    steps =int(2e6)
    S = np.random.choice([-1,1],size=(N,N))
    beta = 0.44
    E=np.zeros((N,N))
    E_plot = np.zeros(steps+1)
    E_tosum=np.zeros((N,N))


    for i in range(N):
        for j in range(N):




            E_tosum[i,j] = (S[bc_i(i,N),bc_j(j,N)]*(S[bc_i(i-1,N),bc_j(j,N)]+\
            S[bc_i(i+1,N),bc_j(j,N)]+S[bc_i(i,N),bc_j(j-1,N)]+\
            S[bc_i(i,N),bc_j(j+1,N)]))

    E_plot[0]=-1*np.sum(E_tosum)





    for step in range(1,steps+1):

        i = np.random.randint(N-1)
        j = np.random.randint(N-1)



        S_new=S.copy()

        #S_new[i,j]=S[i,j] * -1


        del_E = del_E_cal(i,j,S_new)


        if del_E < 0 or np.random.rand() < np.exp(-beta*(del_E))  :

            S_new[i,j]=S[i,j] * -1

            S[i,j]=S_new[i,j]
            #del_E=del_E_cal(i,j,S)


            E_plot[step]=E_plot[step-1] + del_E



        else:

            E_plot[step]=E_plot[step-1]


    plt.plot(E_plot)
    M=(1/N**2)*np.sum(S)

    return E_plot,S

def del_E_cal(i,j,S):

    del_E = 2 * S[i,j] * (S[i-1,j] + S[i+1,j] + S[i,j-1] + S[i,j+1])
    return del_E

def bc_i(i,N):
    if i>=N:
        i=0
    return i

def bc_j(j,N):
    if j>=N:
        j=0
    return j



def bc2(i,j,S,N):

    if i==0:
        S[i-1,j]=S[N-1,j]
    if i==N:
        S[N+1,j]=S[0,j]
    if j==0:
        S[i,j-1]=S[i,N-1]
    if j==N:
        S[i,N+1]=S[0,j]


    return S


