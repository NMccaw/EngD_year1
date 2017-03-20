#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 16:15:28 2017

@author: Nicholas
"""
import numpy as np
from matplotlib import pyplot as plt
def shape_func(xi):
    N = np.zeros(2)
    N[0] = 0.5*(1+xi)
    N[1]= 0.5*(1-xi)
    return N
    
def shape_func_dev(xi):
    dNdxi = np.zeros(2)
    dNdxi[0] = -0.5
    dNdxi[1] = 0.5
    return dNdxi
    
def finite_element_1d(N_elements,q,L):
    E = 5e6
    nu = 0.3
    G =  E/(2*(1+nu))
    kappa = 5/6
    
    h = 2
    I = (h**3)/12
    EI = E*I
    
    
    # Even grid
    nodes, dx = np.linspace(0, L, N_elements+1, retstep=True)
    
    # Location matrix
    LM = np.zeros((2, N_elements), dtype=np.int)
    for e in range(N_elements):
        LM[0, e] = e
        LM[1, e] = e+1
    
    # Global stiffness matrix and force vector
    K = np.zeros((2*(N_elements+1), 2*(N_elements+1)))
    
    F = np.zeros((2*(N_elements+1),))
    dxdxi = dx/2
    integration_points = [-np.sqrt(1/3),np.sqrt(1/3)]
    # Loop over elements
    for e in range(N_elements):
        for xi in integration_points:
            N = shape_func(xi)
            dNdxi = shape_func_dev(xi)
            dNdx = dNdxi * (1/dxdxi)
            B = np.zeros((2,4))
            B[0,2:] = dNdx
            A,C = np.meshgrid([LM[:,e], LM[:,e] + N_elements+1],[LM[:,e], LM[:,e] + N_elements+1])
            
            K[A,C] += np.dot(np.transpose(B),B) * EI * dxdxi
            
            
            F[LM[:,e]] += N * q[e] * dxdxi
    
            
            B = np.zeros((2,4))
            B[1,:2] = dNdx
            B[1,2:] = -N
            K[A,C] += np.dot(np.transpose(B),B) * kappa * h * G * dxdxi
    
    #Boundary Conditions
    fixedNodeW = [0]
    fixedNodeTheta = [0+N_elements+1]
    
    free_dof = np.delete(np.arange(0,2*(N_elements+1)),[fixedNodeW,fixedNodeTheta])
    
    
    active_nodes_A,active_nodes_B = np.meshgrid(free_dof,free_dof)
    
    U = np.linalg.solve(K[active_nodes_A,active_nodes_B],F[free_dof])
    
    
    
    return U,free_dof
  
N_elements = 100
# q is force per unit length 
q = np.zeros(N_elements)
L = 10
elements_per_unit_length = N_elements//L
q[:] = -1000


#U,free_dof = finite_element_1d(N_elements,q,L)
#plt.plot(free_dof[:N_elements],U[:N_elements])

'''
def BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_BEMT,Cl):
    q = np.zeros(N_elements_FEA)
    element_scale_factor = N_elements_FEA/N_elements_BEMT
    residual_factor = element_scale_factor - int(element_scale_factor)
    
#    if residual_factor > 0:
#        for i in range(N_elements_BEMT):
#            if i == N_elements_BEMT - 1:
#                q[i*element_scale_factor:element_scale_factor+(i*element_scale_factor)] = Cl[i]
#                break
#            q[i*element_scale_factor:element_scale_factor+(i*element_scale_factor)] = Cl[i]
#            q[(element_scale_factor+(i*N_elements_BEMT)) + 1] = residual_factor * Cl[i] +(1-residual_factor)\
#             * Cl[i+1]
#            if i == N_elements_BEMT - 1:
#                q[i*element_scale_factor:element_scale_factor+(i*element_scale_factor)] = Cl[i]
    print(residual_factor)
    for i in range(N_elements_BEMT):
        if residual_factor > 0:
            q[i*element_scale_factor:element_scale_factor+(i*element_scale_factor)] = Cl[i]
            q[(element_scale_factor+(i*N_elements_BEMT)) + 1] = residual_factor * Cl[i] +(1-residual_factor) * Cl[i+1]
            temp = 1-residual_factor
            residual_factor = residual_factor - (1-temp)
            print(residual_factor)

        else:
            q[i*element_scale_factor:element_scale_factor+(i*element_scale_factor)] = Cl[i]       
    if residual_factor == 0.0:
        for i in range(N_elements_BEMT):
            q[i*element_scale_factor:element_scale_factor+(i*element_scale_factor)] = Cl[i]
    print(Cl)
    print(q)
        
def test_BEMT_FEA_mesh_interp():
    N_elements_FEA = 100
    N_elements_BEMT = 6
    Cl = np.linspace(0,1000,N_elements_BEMT)
    BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_BEMT,Cl)
    
test_BEMT_FEA_mesh_interp()
'''

def BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_BEMT,Cl,L):
    LM_BEMT = np.zeros((2,N_elements_BEMT))
    LM_FEA = np.zeros((2,N_elements_FEA))
    for e in range(N_elements_BEMT):
        LM_BEMT[0,e] = e*L/N_elements_BEMT
        LM_BEMT[1,e] = (e+1)*L/N_elements_BEMT

    for e in range(N_elements_FEA):
        LM_FEA[0,e] = e*L/N_elements_FEA
        LM_FEA[1,e] = (e+1)*L/N_elements_FEA

    q = np.zeros(N_elements_FEA)
    residual = 0
    
    
    for j in range(N_elements_FEA):
        
        k = j 
        for i in range(N_elements_BEMT):
                          
            if LM_FEA[0,k] >= LM_BEMT[0,i] and LM_FEA[1,k] <= LM_BEMT[1,i]:
                q[k] = Cl[i]
                k+=1
            elif LM_FEA[0,k] < LM_BEMT[1,i] and LM_FEA[1,k] > LM_BEMT[1,i]:
                residual = abs(LM_FEA[1,j] - LM_BEMT[1,i])
                q[k] = residual*Cl[i] + (1-residual)*Cl[i+1]
                k+=1
    return q
                
                
    
def test_BEMT_FEA_mesh_interp():
    N_BEMT = 4
    N_FEA = 10
    Cl = np.linspace(1,10,N_BEMT)
    
    L= 10
    q = BEMT_FEA_mesh_interp(N_FEA,N_BEMT,Cl,L)
    assert np.all(q == np.array([1,1,2.5,4,4,7,7,8.5,10,10]))
test_BEMT_FEA_mesh_interp()
    
    
    
    
def test_timoshenko_uniform_load():
    N_elements = 1000
    L = 10
    # q is pressure per element area not force
    q = np.zeros(N_elements)
    q[:] = -1000
    U,free_dof = finite_element_1d(N_elements,q,L)
    
    assert np.allclose(U[N_elements-1],-0.3906)
test_timoshenko_uniform_load()   