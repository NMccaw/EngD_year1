# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 15:11:28 2016

@author: Nicholas
"""
import numpy as np


def Decomposition_meth():
    """
    Note we want to calculate the rows of u and the columns of l because they
    are upper and lower diagonal matricies
    """
    A=np.array([[2,1,-1],[4,1,0],[-2,-3,8]]) #Here we are setting up the
    n=len(A)    #                               Matricies
    l=np.identity(n) #  set the diagonals of l to 1
    u=np.zeros((n,n))#
    for k in range(n):
        u[k,k]=A[k,k]-np.dot(l[k,:],u[:,k]) #first calculate the diagonals of u
        for j in range(k+1,n):
            u[k,j]=(A[k,j]-np.dot(l[k,:],u[:,j]))/l[k,k] #calculate the rows of u
        for i in range(k+1,n):
            l[i,k]=(A[i,k]-np.dot(l[i,:],u[:,k]))/u[k,k] #calculate the columns
    return l, u