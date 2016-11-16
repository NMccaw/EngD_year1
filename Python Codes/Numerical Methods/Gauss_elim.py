# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 15:54:27 2016

@author: Nicholas
"""

import numpy as np

def gauss_elim():

    A = np.eye(2)
    b = np.array([1.0, 2.0])
    N=len(A)
    print(A)
    M= np.hstack((A,b.reshape(N,1)))
    print(M)
    for j in range(0,N):
        for i in range(j+1,N):
            M[i,:]=M[i,:]-((M[i,j]/M[j,j]*M[j,:]))

    return M