# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 09:02:38 2016

@author: Nicholas
"""
import sympy
import scipy.optimize
import numpy as np
from matplotlib import pyplot as plt


D_0=101.9188
alpha=2.567
k_theta=328.645606
k_r_theta=-211.4672
k_rr=111.70765
r_OH_eq=1
r_HH_eq=1.633

O=np.array([0,0,0])
H1=np.array([0.8,0.6,0])
H2=np.array([-0.8,0.6,0])





(x_H1,y_H1,z_H1,x_H2,y_H2,z_H2,x_O1,y_O1,z_O1)=sympy.symbols('x_H1,y_H1,z_H1,x_H2,y_H2,z_H2,x_O1,y_O1,z_O1')

r_OH1=sympy.sqrt((x_H1-x_O1)**2+(y_H1-y_O1)**2+(z_H1-z_O1)**2)
r_OH2=sympy.sqrt((x_H2-x_O1)**2+(y_H2-y_O1)**2+(z_H2-z_O1)**2)
r_HH1=sympy.sqrt((x_H1-x_H2)**2+(y_H1-y_H2)**2+(z_H1-z_H2)**2)

delta_r_HH = (r_HH1-r_HH_eq)
delta_r_OH1 = (r_OH1-r_OH_eq)
delta_r_OH2 = (r_OH2-r_OH_eq)

#(r_OH1,r_OH2 ,r_HH1,r_HH2)=sympy.symbols('r_OH1,r_OH2,r_HH1,r_HH2')

V=(D_0*(1-sympy.exp(alpha*(r_OH1-r_OH_eq)))**2+D_0*(1-sympy.exp(alpha*(r_OH2-r_OH_eq)))**2)+(0.5*k_theta*(r_HH1-r_HH_eq)**2)+(k_r_theta*(r_HH1-r_HH_eq)*((r_OH1-r_OH_eq)+(r_OH2-r_OH_eq))**2)+(k_rr*(r_OH1-r_OH_eq)*(r_OH2-r_OH_eq))


pot=sympy.lambdify((x_H1,y_H1,z_H1,x_H2,y_H2,z_H2 ,x_O1,y_O1,z_O1),V,'numpy')

def pot_cal(H1):
    return np.array(pot(H1[0], H1[1], 0, H1[2], H1[3], 0,0,0,0 ))


opt=scipy.optimize.minimize(pot_cal, [0.8,0.6,-0.8,0.6] ,tol=1e-4)


fO=(-V.diff(x_O1),-V.diff(y_O1), -V.diff(z_O1))
fH1=(-V.diff(x_H1), -V.diff(y_H1),-V.diff(z_H1))
fH2=(-V.diff(x_H2), -V.diff(y_H2),-V.diff(z_H2))

dv_dO=sympy.lambdify((x_H1,y_H1,z_H1,x_H2,y_H2,z_H2 ,x_O1,y_O1,z_O1),fO,'numpy')
dv_dH1=sympy.lambdify((x_H1,y_H1,z_H1,x_H2,y_H2,z_H2 ,x_O1,y_O1,z_O1),fH1,'numpy')
dv_dH2=sympy.lambdify((x_H1,y_H1,z_H1,x_H2,y_H2,z_H2 ,x_O1,y_O1,z_O1),fH2,'numpy')


def fO(H1,H2,O):
    return np.array(dv_dO(H1[0], H1[1], H1[2],H2[0], H2[1], H2[2], O[0], O[1], O[2] ))

def fH1(H1,H2,O):
    return np.array(dv_dH1(H1[0], H1[1], H1[2],H2[0], H2[1], H2[2], O[0], O[1], O[2] ))

def fH2(H1,H2,O):
    return np.array(dv_dH2(H1[0], H1[1], H1[2],H2[0], H2[1], H2[2], O[0], O[1], O[2] ))





def accel(x,mass):
    N=len(x)
    a=np.zeros_like(x)







    a[0,:]=fH1(x[0,:],x[1,:],x[2,:])/mass[1]


    a[1,:]=fH2(x[0,:],x[1,:],x[2,:])/mass[1]


    a[2,:]=fO(x[0,:],x[1,:],x[2,:])/mass[0]




    return a


def velo_vel():
    x=np.array([[opt.x[0],opt.x[1],0],[opt.x[2],opt.x[3],0],[0,0,0]])
    print(x[2,:])
    mass=[15.999,1.008]
    del_t=0.001
    steps=1000#int(1/del_t)
    a=np.zeros_like(x)
    a=accel(x,mass)
    v=np.zeros_like(a)
    positions = np.zeros((steps, x.shape[0], x.shape[1]))
    positions[0,:,:] = x.copy()

    for step in range(steps):



        x=x + del_t*v +0.5*del_t**2*a


        #x=peri_bound(x,L)


        v_star=v+0.5*del_t*a

        a=accel(x,mass)


        v=v_star+0.5*del_t*a
        positions[step,:,:] = x.copy()







    return positions



positions= velo_vel()





def plot():
    fig = plt.figure()

    ax1 = fig.add_subplot(131, projection='3d')
    ax1.plot3D(positions[:,0,0], positions[:,1,0], positions[:,2,0], 'b^-')
    ax1.set_title('Oxygen');
    ax1.set_xlabel(r'$x$');ax1.set_ylabel(r'$y$');ax1.set_zlabel(r'$z$')


    ax2 = fig.add_subplot(132, projection='3d')
    ax2.plot3D(positions[:,1,0]-positions[0,1,0], positions[:,1,1]-positions[0,1,1],positions[:,1,2]-positions[0,1,2], 'g*-');
    ax2.set_title('Hydrogen 1');
    ax2.set_xlabel(r'$x$');ax2.set_ylabel(r'$y$');ax2.set_zlabel(r'$z$')
    #ax2.view_init(90)

    ax3 = fig.add_subplot(133, projection='3d')
    ax3.plot3D(positions[:,2,0]-positions[0,2,0], positions[:,2,1]-positions[0,2,1], positions[:,2,2]-positions[0,2,2], 'g*-');
    ax3.set_title('Hydrogen 2');
    ax3.set_xlabel(r'$x$');ax3.set_ylabel(r'$y$');ax3.set_zlabel(r'$z$')
    #ax3.view_init(x)
