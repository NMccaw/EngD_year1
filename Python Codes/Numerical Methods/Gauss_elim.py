# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 15:54:27 2016

@author: Nicholas
"""

import numpy as np

def gauss_elim():

    A=np.array([[3,0,1],[6,2,4],[9,2,6]])
    b=np.array([4,10,15])
    N=len(A)
    print(A)
    M= np.hstack((A,b.reshape(N,1)))
    for j in range(0,N):
        for i in range(j+1,N):
            M[i,:]=M[i,:]-((M[i,j]/M[j,j]*M[j,:]))

    return M