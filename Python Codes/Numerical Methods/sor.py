# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:33:32 2016

@author: Nicholas
"""

import numpy as np

def sor():
    A=np.array([[3,1],[1,3]])
    n=len(A)
    b=np.array([5,7])
    A_aug=np.zeros((n,n))
    b_aug=np.zeros(n)
    steps=15
    x=np.zeros((n))
    c=np.zeros((n))
    omega=1.035

    for i in range(n):
        A_aug[:,i]=A[:,i]/A[i,i]
        b_aug[i]=b[i]/A[i,i]


    Au=np.triu(np.identity(n)-A_aug)
    Al=np.tril(np.identity(n)-A_aug)

    for step in range(1,steps):
        c=np.dot(Al,x)+np.dot(Au-np.identity(n),x)+b_aug

        x=x+(omega*c)

    print(x)