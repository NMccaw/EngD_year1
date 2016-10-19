# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 13:55:57 2016

@author: Nicholas
"""
import numpy as np


def gauss_seidal():
    A=np.array([[3,1],[1,3]])
    n=len(A)
    b=np.array([5,7])
    A_aug=np.zeros((n,n))
    b_aug=np.zeros(n)
    steps=20
    x=np.zeros((n))

    for i in range(n):
        A_aug[:,i]=A[:,i]/A[i,i]
        b_aug[i]=b[i]/A[i,i]


    Au=np.triu(np.identity(n)-A_aug)
    Al=np.tril(np.identity(n)-A_aug)

    print(Au)


    for step in range(1,steps):
        x=np.dot(Al,x)+np.dot(Au,x)+b_aug
    print(x)

