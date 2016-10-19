# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 11:37:17 2016

@author: Nicholas
"""
import numpy as np
def jacobi():
    A=np.array([[3,1],[1,3]])
    b=np.array([5,7])
    n=len(A)

    x=np.zeros(n)

    A_aug=np.zeros((n,n))
    b_aug=np.zeros(n)
    for i in range(n):
        A_aug[:,i]=A[:,i]/A[i,i]
        b_aug[i]=b[i]/A[i,i]
    Al=np.tril(-A_aug)

    Au=np.triu(-A_aug)
    P=(np.identity(n)-A_aug)


    for k in range(100):

        x=(np.dot(P,x))+b_aug


    return x
