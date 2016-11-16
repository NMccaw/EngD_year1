# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 17:05:46 2016

@author: Nicholas
"""
import numpy as np
from matplotlib import pyplot as plt
import sympy
from mpl_toolkits.mplot3d import Axes3D


def d_phi_external(r,N,Rc,ei,ej,eps,sigma):
    
    """
    Compute the external potential derivative
    
    Parameters
    ----------
    
    r : Nx3 array
        distance array from the atom to all other atoms
    N : integer
        Number of atoms       
    Rc: integer
        cut-off radius
    ei,ej: float
            charge of atom
    sigma: float
            coupling between atoms in different molecules
        
    Returns
    -------
    
    d_phi: array
            potential derivative
    """
    #external potential derivative equation
    if r <= Rc:
        d_phi=-((ei*ej)/(4*np.pi*r**2))-(24*sigma*eps*(((2*((sigma/r)**14))-((sigma/r)**8))))
    else:
        d_phi=0


    return d_phi


def mol_dyn2(x,N,L,e,eps,sigma,mass):

    """
    Compute the molecular dynamics acceleration
    
    Parameters
    ----------
    
    x : Nx3 array
        positions of atoms
    N : integer
        Number of atoms       
    L : integer
        size of box
    e: array
        charge of all atoms
    sigma: float
            coupling between atoms in different molecules
    mass: array
        mass of oxygen and two hydrogen atoms
        
    Returns
    -------
    
    a: array
        acceleration matrix
            
    """
    
    Rc=5
    # set accelerations to zero such that accelerations are only calculated 
    # for the distance between molecules.
    a=np.zeros_like(x)
    rij=np.zeros_like(x)

    for i in range(N):


        rij=x[i,:]-x[:,:]
        # include the closest atom including the boundary condiitons 
        rij[np.abs(rij)>L/2] -= np.sign(rij[np.abs(rij)>L/2])*L


        for j in range(i+1,N):
            mag=np.sqrt(np.dot(rij[j,:],rij[j,:]))



            # set correct e
            ei=e[i]
            ej=e[j]

            # use modulus function to choose correct sigma and epsilon 
            if i%3 == j%3:
                if i%3 == 0:
                    sigma_ij=3.166
                    eps_ij=0.1554
                else:
                    sigma_ij=0
                    eps_ij=0

            else:
                sigma_ij=0
                eps_ij=0






            d_phi_dr=d_phi_external(mag,N,Rc,ei,ej,eps_ij,sigma_ij)

            a[i,:]-=rij[j,:]*d_phi_dr/mass[i%3]

            a[j,:]+=rij[j,:]*d_phi_dr/mass[j%3]

    return a

def peri_bound(x,L):
    """
      Enforces periodic boundary conditions
    
    Parameters
    ----------
    
    x : Nx3 array
        positions of atoms      
    L : integer
        size of box
        
    Returns
    -------
    
    x :  array
        position matrix  
    """
    x[x<0]+=L
    x[x>L]-=L
    return x




def time_evo():
    """
      Sets up molecules in box and applies time evolution of molecules
    
    Parameters
    ----------
    
    none
        
    Returns
    -------
    
    x_plot : array
            position matrix with time
        
    a : array
        accelration of atoms
        
    T_plot : array
            temperature of atoms with time
    time: array
            all time steps in array
    """
    # set time step
    del_t=0.001
    t=5
    steps=int(t/del_t)
    L=7
    Rc=5


    x=np.zeros((24,3))
    # set up box of atoms
    # Oxygens
    x[0,:] = [1.75,1.75,1.75]
    x[3,:] = [5.25,1.75,1.75]
    x[6,:] = [1.75,5.25,1.75]
    x[9,:] = [1.75,1.75,5.25]
    x[12,:] = [5.25,5.25,1.75]
    x[15,:] = [5.25,1.75,5.25]
    x[18,:] = [1.75,5.25,5.25]
    x[21,:] = [5.25,5.25,5.25]
    # Hydrogens
    x[1,:] = x[0,:] + [0.8,0.6,0]
    x[4,:] = x[3,:] + [0.8,0.6,0]
    x[7,:] = x[6,:] + [0.8,0.6,0]
    x[10,:] = x[9,:] + [0.8,0.6,0]
    x[13,:] = x[12,:] + [0.8,0.6,0]
    x[16,:] = x[15,:] + [0.8,0.6,0]
    x[19,:] = x[18,:] + [0.8,0.6,0]
    x[22,:] = x[21,:] + [0.8,0.6,0]
    # Hydrogens
    x[2,:] = x[0,:] + [-0.8,0.6,0]
    x[5,:] = x[3,:] + [-0.8,0.6,0]
    x[8,:] = x[6,:] + [-0.8,0.6,0]
    x[11,:] = x[9,:] + [-0.8,0.6,0]
    x[14,:] = x[12,:] + [-0.8,0.6,0]
    x[17,:] = x[15,:] + [-0.8,0.6,0]
    x[20,:] = x[18,:] + [-0.8,0.6,0]
    x[23,:] = x[21,:] + [-0.8,0.6,0]

    # set appropriate e
    e=np.array([-0.82,0.41,0.41,-0.82,0.41,0.41,-0.82,0.41,0.41,-0.82,0.41,0.41
            ,-0.82,0.41,0.41,-0.82,0.41,0.41,-0.82,0.41,0.41,-0.82,0.41,0.41])




    eps=np.array([0.1554,0,0])


    sigma=np.array([3.166,0,0])

    mass=np.array([15.999,1.008,1.008])






  
    N=len(x)
    x_plot=np.zeros((steps,N,3))

    # set up acceleration arrays
    a1=np.zeros_like(x)
    a2=np.zeros_like(x)

    a=np.zeros_like(x)
    v=np.zeros_like(x)
    time=np.zeros(steps)
    T_plot=[]

    for step in range(steps):
        
        # velocity verlet section
        x=x + del_t*v +0.5*del_t**2*a


        x=peri_bound(x,L)


        v_star=v+0.5*del_t*a

        a1=mol_dyn2(x,N,L,e,eps,sigma,mass)
        a2= accel_int(x,mass)

        a=a1+a2

        v=v_star+0.5*del_t*a
        x_plot[step,:,:]=x
        T=temp(v,N,L)
        T_plot.append(T)
        time[step]=step*del_t


    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs=x_plot[:,0,0],ys=x_plot[:,0,1],zs=x_plot[:,0,2])

    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d')
    for i in range(N):
        ax.plot(xs=x_plot[:,i,0],ys=x_plot[:,i,1],zs=x_plot[:,i,2])

    fig = plt.figure(10)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs=x_plot[:,7,0],ys=x_plot[:,7,1],zs=x_plot[:,7,2])


    fig=plt.figure(3)
    plt.plot(time, T_plot)

    fig=plt.figure(4)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs=x_plot[:,0,0],ys=x_plot[:,0,1],zs=x_plot[:,0,2])
    ax.plot(xs=x_plot[:,1,0],ys=x_plot[:,1,1],zs=x_plot[:,0,2])
    ax.plot(xs=x_plot[:,2,0],ys=x_plot[:,2,1],zs=x_plot[:,0,2])

    fig=plt.figure(5)
    plt.plot(time, x_plot[:,0,1])
    plt.plot(time, x_plot[:,1,1])
    plt.plot(time, x_plot[:,2,1])


    return x_plot,a,T_plot,time

def temp(v,N,L):
    """
      Calculates temperature of system
    
    Parameters
    ----------
    
    v: array
        velocity of each atom array 
    N : integer
        number of atoms in system
    L : integer
        size of box
        
    Returns
    -------
    
    T : integer
        temperature of system
    """


    E_kin=np.abs(v)**2

    T=(L**2)/(3*N)*np.sum(E_kin)
    return T

#######Internal accelerations ##########




#set sympy symbols
(x_H1,y_H1,z_H1,x_H2,y_H2,z_H2,x_O1,y_O1,z_O1)=sympy.symbols('x_H1,y_H1,z_H1,x_H2,y_H2,z_H2,x_O1,y_O1,z_O1')

# internal constants
D_0=101.9188
alpha=2.567
k_theta=328.645606
k_r_theta = -211.4672
k_rr=111.70765
r_OH_eq=1
r_HH_eq=1.633

# distance between internal atoms
r_OH1=sympy.sqrt((x_H1-x_O1)**2+(y_H1-y_O1)**2+(z_H1-z_O1)**2)
r_OH2=sympy.sqrt((x_H2-x_O1)**2+(y_H2-y_O1)**2+(z_H2-z_O1)**2)
r_HH1=sympy.sqrt((x_H1-x_H2)**2+(y_H1-y_H2)**2+(z_H1-z_H2)**2)


 
# potential equation insympy form
V=(D_0*(1-sympy.exp(alpha*(r_OH1-r_OH_eq)))**2+D_0*(1-sympy.exp(alpha*(r_OH2-r_OH_eq)))**2)+(0.5*k_theta*(r_HH1-r_HH_eq)**2)+(k_r_theta*(r_HH1-r_HH_eq)*((r_OH1-r_OH_eq)+(r_OH2-r_OH_eq))**2)+(k_rr*(r_OH1-r_OH_eq)*(r_OH2-r_OH_eq))

# turn the potential into a function
pot=sympy.lambdify((x_H1,y_H1,z_H1,x_H2,y_H2,z_H2 ,x_O1,y_O1,z_O1),V,'numpy')

#turn pot into a function to unpack array
def pot_cal(x):
    return np.array(pot(x[1,0],x[1,1],x[1,2], x[2,0],x[2,1],x[2,2],0,0,0))


opt=scipy.optimize.minimize(pot_cal, [0.8,0.6,-0.8,0.6] ,tol=1e-4)


# differentiate potential equation w.r.t appropriate variables
fO=(-V.diff(x_O1),-V.diff(y_O1), -V.diff(z_O1))
fH1=(-V.diff(x_H1), -V.diff(y_H1),-V.diff(z_H1))
fH2=(-V.diff(x_H2), -V.diff(y_H2),-V.diff(z_H2))
# turn the differentiated potential into functions
dv_dO=sympy.lambdify((x_H1,y_H1,z_H1,x_H2,y_H2,z_H2 ,x_O1,y_O1,z_O1),fO,'numpy')
dv_dH1=sympy.lambdify((x_H1,y_H1,z_H1,x_H2,y_H2,z_H2 ,x_O1,y_O1,z_O1),fH1,'numpy')
dv_dH2=sympy.lambdify((x_H1,y_H1,z_H1,x_H2,y_H2,z_H2 ,x_O1,y_O1,z_O1),fH2,'numpy')

# unpack arrays
def fO(H2,H1,O):
    return np.array(dv_dO(H1[0],H1[1],H1[2], H2[0],H2[1],H2[2],O[0],O[1],O[2]))

def fH1(H2,H1,O):
    return np.array(dv_dH1(H1[0],H1[1],H1[2], H2[0],H2[1],H2[2],O[0],O[1],O[2]))

def fH2(H2,H1,O):
    return np.array(dv_dH2(H1[0],H1[1],H1[2], H2[0],H2[1],H2[2],O[0],O[1],O[2]))


def accel_int(x,mass):
    
    """
     Computes the internal acceleration 
    
    Parameters
    ----------
    
    x : array 
        positions of all atoms
    mass: array 
            mass of oxygen, hydrogen and hydrogen
        
    Returns
    -------
        
    a : array
        accelration of atoms by internal forces   
    
    """
    N=len(x)
    a=np.zeros_like(x)




    for mol in range(0,N,3):

        
        # calculate acceleration of Oxygen
        a[mol,:]=fO(x[mol+2,:],x[mol+1,:],x[mol,:])/mass[0]

        # calculate acceleration of Hydrogens
        a[mol+1,:]=fH1(x[mol+2,:],x[mol+1,:],x[mol,:])/mass[1]


        a[mol+2,:]=fH2(x[mol+2,:],x[mol+1,:],x[mol,:])/mass[1]




    return a



x_plot,a,T_plot, time = time_evo()



plt.plot( x_plot[:,6,0]-x_plot[0,6,0])
plt.plot( x_plot[:,6,1]-x_plot[0,6,1])
plt.plot( x_plot[:,6,2]-x_plot[0,6,2])

plt.plot( x_plot[:,7,0]-x_plot[0,7,0])
plt.plot( x_plot[:,7,1]-x_plot[0,7,1])
plt.plot( x_plot[:,7,2]-x_plot[0,7,2])

plt.plot( x_plot[:,8,0]-x_plot[0,8,0])
plt.plot( x_plot[:,8,1]-x_plot[0,8,1])
plt.plot( x_plot[:,8,2]-x_plot[0,8,2])

plt.plot( x_plot[:,0,2]-x_plot[0,0,2])
plt.plot( x_plot[:,1,2]-x_plot[0,1,2])
plt.plot( x_plot[:,2,2]-x_plot[0,2,2])



