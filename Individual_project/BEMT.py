#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 11:07:03 2016

Blade Element Momentum Theory code. Written by N.McCaw 2016. Based on Fortran 
code by Molland and Turnock.

This code follows the algorithm set out in figure 15.7 on P351 of 'Ship Resistance
and Propuslsion' by A.Molland and S.Turnock

Estimated to be 42x faster than Matlab version

@author: Nicholas
"""


import numpy as np
from matplotlib import pyplot as plt
import pytest

def BEMT(J,P_D,BA_ratio,N_blades,chord,relative_pitch_distribution,
                        lift_curve_slope,root_thickness):
    '''
    Main Blade Element Momentum Theory Function. Uses algorithm form referenced
    book.
    
    Parameters
    ----------
    
    J: float
        advance ratio
    P_D: float
        pitch diameter ratio at 0.7 radius
    BA_ratio: float
            blade area ratio
    N_blades: int
            number of blades
    chord: array of floats
        chord of blade at different radial positions
    relative_pitch_distribution: array of float
        pitch distribution at radial positions relative to pitch distibution at
        x_R = 0.7
    lift_curve_slope: float
                lift curve slope of blade
    root_thickness: float
                thickness at blade root
    
    Returns
    -------
    Total_Thrust: float
                total thrust coefficient along blade
    Total_Torque: float
                total torque coefficient along blade
    Open_eta: float
            open water efficiency
    '''
    assert J <= 1.1 and J >= 0.3,' Advance ratio must be between 0.3 and 1.1 \
    but has value {}'.format(J)
    assert P_D >= J+0.03 and P_D >= 0.6 and P_D <= 1.2, 'Pitch/Diameter Ratio at \
    0.7R is outside limits for this J value'
    assert BA_ratio >= 0.4 and BA_ratio <= 0.8,' Blade area ratio is out limits.\
    Must be between 0.4 and 0.8 but has value {}'.format(BA_ratio)
    assert N_blades >= 3 and N_blades <=5,'Number of blades must be between 3 \
and 5 but is of value {}'.format(N_blades)
    assert type(J) == np.float or np.float64,'J is not of type float but is of type {}'.format(type(J))
    assert type(P_D) == np.float or np.float64,'P_D is not of type float but is of type {}'.format(type(P_D))
    assert type(BA_ratio) == np.float or np.float64,'BA_ratio is not of type float but is of type {}'.format(type(BA_ratio))
    assert type(N_blades) == int,'N_blades is not of type int but is of type {}'.format(type(N_blades))
    assert type(chord) == np.array or np.ndarray,'chord is not of type array but is of type {}'.format(type(chord))
    assert type(relative_pitch_distribution) == np.array or np.ndarray,'relative_pitch_distribution is not of type array but is of type {}'.format(type(relative_pitch_distribution))
    assert type(lift_curve_slope) == np.float or np.float64,'lift_curve_slope is not of type float but is of type {}'.format(type(lift_curve_slope))
    assert type(root_thickness) == np.float or np.float64,'root_thickness is not of type float but is of type {}'.format(type(root_thickness))
    # create arrays
    x_R = np.zeros_like(chord)
    phi_plus_alhpa = np.zeros_like(chord)
    tan_psi = np.zeros_like(chord)
    local_P_D = np.zeros_like(chord)
    K_series = np.zeros_like(chord)
    chord_diameter = np.zeros_like(chord)
    thickness_chord = np.zeros_like(chord)
    KT_dx = np.zeros_like(chord)
    Cl = np.zeros_like(chord)
    Cd = np.zeros_like(chord)
    KQ_dx = np.zeros_like(chord)
    phi = np.zeros_like(chord)
    alpha = np.zeros_like(chord)
    K2 = np.zeros_like(chord)
    MC = np.zeros_like(chord)
    MT =  np.zeros_like(chord)
    CC= np.zeros_like(chord)
    
    LI = 0.8*12.67
    LM = 0.8*11.3
    VV = 0
    for i in range(1,9):
        # calculate local pitch
        local_P_D[i] = P_D * relative_pitch_distribution[i]
        # set value of x_R for each element
        x_R[i] = (i+1)/10
        # compute chord/diameter from blade geometry
        chord_diameter[i]=chord[i]*BA_ratio*4.0/(N_blades*0.5)
        # compute Nc/D for later use
        Nc_D = N_blades*chord_diameter[i]
        # thickness distribution from rootthickness
        thickness_distribution=root_thickness-(root_thickness*0.935)*x_R[i]
        # thickness/chord
        thickness_chord[i]=thickness_distribution/chord_diameter[i]
        # set initial angle of attack to 0
        alpha[i] = 0.0
        # set number of alpha iterations to 0
        alpha_iterations = 0
        # set convergence flag to 0
        alpha_converge = 0
        # start of alpha convergenc loop
        while alpha_iterations < 2000 and alpha_converge == 0: 
            alpha_iterations += 1
            # calculate tan of inflow angle
            tan_psi[i] = J / (np.pi*x_R[i])
            # induced flow angle plus angle of attack
            phi_plus_alhpa[i] = np.arctan2(local_P_D[i],(np.pi*x_R[i]))
            # inflow angle
            phi[i] = phi_plus_alhpa[i] - alpha[i]
            # ideal efficieny
            eta_ideal = tan_psi[i] / np.tan(phi[i])
            # initially set efficiency to ideal  
            eta  = 0.9 * eta_ideal
            gamma = 0.0 # zero for ideal efficiency
            KF = np.tan(phi[i])*x_R[i]
            # Goldstein Correction
            SF=N_blades/(2.0*np.tan(phi[i])*x_R[i])-0.5
            F1=np.cosh(SF)
            F2=np.cosh(SF*x_R[i])
            F3=F2/F1
            F4=np.arccos(F3)
            K=2.0*F4/np.pi
            # store factor in array
            K_series[i] = K
            # set number of efficiency iterations to 0
            eta_iterations = 0
            eta_converge = 0
            while eta_iterations < 2000 and eta_converge == 0:
                # axial inflow factor
                a = (1 - eta_ideal)/(eta_ideal + (tan_psi[i]**2)/eta)
                # local thrust coefficient
                KT_dx[i] = np.pi*(J**2) *x_R[i]*K*a*(1+a)
                # circumferential inflow factor (a')
                a_p =1-eta_ideal*(1+a)
                # lift coefficient
                Cl[i] = KT_dx[i]/(((np.pi**2)/4)*(N_blades*chord_diameter[i])*(x_R[i]**2)*((1-a_p)**2) \
                       * (1/np.cos(phi[i]))*(1-np.tan(phi[i])*np.tan(gamma)))
                
                # compute drag coefficient from angle of attack and geometry
                F6=0.0107+(alpha[i]+1.0)*(-0.0015+alpha[i]*(0.0015+0.000965*(alpha[i]-1.0)));
                F7=-0.03833+(alpha[i]+1.0)*(0.0133+alpha[i]*(-0.015-0.01166*(alpha[i]-1.0)));
                F8=0.8193+(alpha[i]+1.0)*(-0.0138+alpha[i]*(0.0903+0.079*(alpha[i]-1.0)));
                F9=-3.076+(alpha[i]+1.0)*(-0.0728+alpha[i]*(-0.3162-0.2437*(alpha[i]-1.0)));
                Cd[i]=F6+thickness_chord[i]*(F7+(thickness_chord[i]-0.06)*(F8+F9*(thickness_chord[i]-0.12)));
                
                # update gamma
                gamma = np.arctan2(Cd[i],Cl[i])
                # compute efficieny
                eta_calculated = tan_psi[i]/(np.tan((phi[i]+gamma)))
                eta_iterations += 1
                if np.abs(eta_calculated - eta) < 0.001:
                    # if old efficieny = new efficiency
                    # set convergence flag to 1 to stop iteration
                    eta_converge = 1
                    # camber correction
                    # geometric 
                    K2[1]=1.0+2.857*(BA_ratio-0.4)*(BA_ratio-0.6)*(BA_ratio-0.9)
                    K2[2]=1.19+(BA_ratio-0.4)*(BA_ratio-0.6)*(0.267+(BA_ratio-0.9)*0.1665)
                    K2[3]=1.3+(BA_ratio-0.4)*(0.1+(BA_ratio-0.6)*(1+(BA_ratio-0.9)*0.43))
                    K2[4]=1.54+(BA_ratio-0.4)*(0.15+(BA_ratio-0.6)*(1.1+0.286*(BA_ratio-0.9)))
                    K2[5]=1.67+(BA_ratio-0.4)*(0.55+(BA_ratio-0.6)*(0.1667+(BA_ratio-0.9)*2.9465))
                    K2[6]=1.8+(BA_ratio-0.4)*(0.75+(BA_ratio-0.6)*(0.3+(BA_ratio-0.9)*2.835))
                    K2[7]=1.8+(BA_ratio-0.4)*(1.0+(BA_ratio-0.6)*(1.333+(BA_ratio-0.9)*1.905))
                    K2[8]=1.75+(BA_ratio-0.4)*(1.25+(BA_ratio-0.6)*(1.5+(BA_ratio-0.9)*3.55))
                    # from camber correction curves
                    # KF is lambda from goldstein correction
                    U1 = -0.65*KF*KF+1.1*KF+0.664
                    U2 = (0.85+(KF-0.3)*(-4.0+(KF-0.4)*(15.42-47.95*(KF-0.5))))
                    U2 = -0.09+(KF-0.2)*U2
                    U3 = (1.375+(KF-0.3)*(-3.75+(KF-0.4)*(20.85-75.7875*(KF-0.5))))
                    U3 = -0.2+(KF-0.2)*U3;
                    K1 = U1+(BA_ratio-0.4)*(U2+U3*(BA_ratio-0.8))
                    # camber correction
                    CC[i] = K1*K2[i]
                    
                    # compute angle of attack using camber correction
                    # MC = camber/chord
                    # MT = camber/thickness
                    # multiply by pi/180 as lift curve slope is degrees
                    if VV==1:
                        MC[i]=MT[i]*thickness_chord[i];
                        alpha_calculated=(CC[i]*Cl[i]/LI-MC[i])*LM/0.1097+(MC[i]*LI/lift_curve_slope)
                        
                    elif VV==5:
                        MC[i]=0.5*thickness_chord[i]
                        alpha_calculated =(CC[i]*Cl[i]/LI-MC[i])*LM/0.1097+(Cl[i]/lift_curve_slope)
                        alpha_calculated = alpha_calculated*np.pi/180
                        MT[i]=MC[i]/thickness_chord[i]
                        
                    else:
                        alpha_calculated = Cl[i]/lift_curve_slope
                        alpha_calculated = alpha_calculated*np.pi/180
                        MC[i]=Cl[i]/LI
                        MC[i]=MC[i]*CC[i]
                        
                        
                        if MC[i]<(0.5*thickness_chord[i]):
                            MT[i]=MC[i]/thickness_chord[i]
                            
                        else:
                            VV=5
                            MC[i]=0.5*thickness_chord[i]
                            alpha_calculated = (CC[i]*Cl[i]/LI-MC[i])*LM/0.1097+(Cl[i]/lift_curve_slope)
                            alpha_calculated = alpha_calculated*np.pi/180
                            MT[i]=MC[i]/thickness_chord[i]                
                else: 
                    # set the new efficiency to calculated efficiency
                    eta = eta_calculated        
                
            if np.abs((alpha_calculated - alpha[i])/alpha_calculated)<0.1:
                # compute local torque coefficient
                KQ_dx[i]=4.935*J*(x_R[i]**3)*K*a_p*(1+a)
                # set convergence flag to 1 to stop iteration
                alpha_converge = 1
            else:
                alpha[i] = (alpha_calculated + alpha[i]) / 2.0
                
        
    # set tip coefficients to 0
    x_R[9]=1.0;
    KT_dx[9]=0.0
    KQ_dx[9]=0.0
    K_series[9]=0.0
    local_P_D[9]=local_P_D[8]
    alpha[9]=alpha[8]
    
    print(KT_dx)
    # Integrate over blade to compute total coefficients. Using simpsons rule
    Total_Thrust=(KT_dx[1]+2.0*(KT_dx[3]+KT_dx[5]+KT_dx[7])+4.0*(KT_dx[2]+KT_dx[4]+KT_dx[6]+KT_dx[8]))   
    Total_Thrust=0.1*Total_Thrust/3.0  
    Total_Torque=(KQ_dx[1] + 2.0*(KQ_dx[3]+KQ_dx[5]+KQ_dx[7])+4.0*(KQ_dx[2]+KQ_dx[4]+KQ_dx[6]+KQ_dx[8]))  
    Total_Torque=(0.1*Total_Torque)/3.0  
    # open water efficiency
    Open_eta=(J*Total_Thrust)/(2.0*3.1416*Total_Torque)
    
    return Total_Thrust , Total_Torque, Open_eta
if __name__ == '__main__':
    pytest.main('-v BEMT.py')
    
    J=0.8      #Advance Ratio
    P_D=0.95      # Pitch/Diamter Ratio
    BA_ratio=0.4   #  Blade Area Ratio
    N_blades=4    # number of blades
    x_R = []
    local_P_D = []
    chord = []
    chord_diameter = []
    thickness_distibution = []
    f = open('propgeom.txt', 'r',encoding= 'iso-8859-1')
    
 
    lines = 0
    data = f.readlines()
    
    
    for i in range(len(data)):
        if 'Number of Blade:' in data[i]:
            N_blades = data[i].split(' ')
            N_blades = int(N_blades[-1])
        if 'Propeller Diameter' in data[i]:
            D = data[i].split(' ')
            D = float(D[-2])
        if 'r/R' in data[i]:
            columns = data[i:]
            for j in range(len(data) - i-1): 
                values = list(filter(None,columns[j+1].split(' ')))
                x_R.append(float(values[0]))
                local_P_D.append(float(values[1]))
                chord_diameter.append(float(values[4]))
                thickness_distibution.append(float(values[6]))
                
    chord = [i*D for i in chord_diameter]      
    print(chord)  
    
    #Set chord at spanwise positions
    chord = np.zeros(10)
    chord[1] = 0.208    # chord at x_R = 0.2
    chord[2] = 0.241    # chord at x_R = 0.3
    chord[3] = 0.263    #       .
    chord[4] = 0.276    #       .
    chord[5] = 0.279    #       .
    chord[6] = 0.269    #       .
    chord[7] = 0.241    #       .
    chord[8] = 0.184    # chord at x_R = 0.9
    # Pitch distribution along blade length relative to pitch at R = 0.7
    relative_pitch_distribution = np.zeros_like(chord)
    relative_pitch_distribution[1] = 0.888  # relative pitch at x_R = 0.2
    relative_pitch_distribution[2] = 1.008  # relative pitch at x_R = 0.3
    relative_pitch_distribution[3] = 1.055  #           .
    relative_pitch_distribution[4] = 1.060  #           .
    relative_pitch_distribution[5] = 1.039  #           .
    relative_pitch_distribution[6] = 1.000  #           .
    relative_pitch_distribution[7] = 0.948  #           .
    relative_pitch_distribution[8] = 0.888  # relative pitch at x_R = 0.9

    # blade parameters
    lift_curve_slope = 0.97 # lift curve slope 
    root_thickness = 0.045 #thickness_distibution[0]*D  # blade thickness at root
    
    TT,TQ,TE = BEMT(J,P_D,BA_ratio,N_blades,chord,relative_pitch_distribution,
                        lift_curve_slope,root_thickness)
    
    print('Thrust coefficient = {}'.format(TT))
    print('Torque coefficient = {}'.format(TQ))
    print('Open water efficiency = {}'.format(TE))
    
    
    ################# Tests ####################
  
def test_fortran1():
    J=0.8       #Advance Ratio
    P_D=0.95      # Pitch/Diamter Ratio
    BA_ratio=0.4   #  Blade Area Ratio
    N_blades=4    # number of blades
    
    #Set chord at spanwise positions
    chord = np.zeros(10)
    chord[1] = 0.208    # chord at x_R = 0.2
    chord[2] = 0.241    # chord at x_R = 0.3
    chord[3] = 0.263    #       .
    chord[4] = 0.276    #       .
    chord[5] = 0.279    #       .
    chord[6] = 0.269    #       .
    chord[7] = 0.241    #       .
    chord[8] = 0.184    # chord at x_R = 0.9
    # Pitch distribution along blade length relative to pitch at R = 0.7
    relative_pitch_distribution = np.zeros_like(chord)
    relative_pitch_distribution[1] = 0.888  # relative pitch at x_R = 0.2
    relative_pitch_distribution[2] = 1.008  # relative pitch at x_R = 0.3
    relative_pitch_distribution[3] = 1.055  #           .
    relative_pitch_distribution[4] = 1.060  #           .
    relative_pitch_distribution[5] = 1.039  #           .
    relative_pitch_distribution[6] = 1.000  #           .
    relative_pitch_distribution[7] = 0.948  #           .
    relative_pitch_distribution[8] = 0.888  # relative pitch at x_R = 0.9
    # blade parameters
    lift_curve_slope = 0.97 # lift curve slope 
    root_thickness = 0.045  # blade thickness at root
    
    TT,TQ,TE = BEMT(J,P_D,BA_ratio,N_blades,chord,relative_pitch_distribution,
                        lift_curve_slope,root_thickness)
    assert np.allclose(TT,0.113,atol = 1e-3)
    assert np.allclose(TQ,0.019,atol = 1e-3)
    assert np.allclose(TE,0.746,atol = 1e-2)
    
def test_fortran2():
    J=0.7       #Advance Ratio
    P_D=0.8      # Pitch/Diamter Ratio
    BA_ratio=0.8   #  Blade Area Ratio
    N_blades=3    # number of blades
    
    #Set chord at spanwise positions
    chord = np.zeros(10)
    chord[1] = 0.208    # chord at x_R = 0.2
    chord[2] = 0.241    # chord at x_R = 0.3
    chord[3] = 0.263    #       .
    chord[4] = 0.276    #       .
    chord[5] = 0.279    #       .
    chord[6] = 0.269    #       .
    chord[7] = 0.241    #       .
    chord[8] = 0.184    # chord at x_R = 0.9
    # Pitch distribution along blade length relative to pitch at R = 0.7
    relative_pitch_distribution = np.zeros_like(chord)
    relative_pitch_distribution[1] = 0.888  # relative pitch at x_R = 0.2
    relative_pitch_distribution[2] = 1.008  # relative pitch at x_R = 0.3
    relative_pitch_distribution[3] = 1.055  #           .
    relative_pitch_distribution[4] = 1.060  #           .
    relative_pitch_distribution[5] = 1.039  #           .
    relative_pitch_distribution[6] = 1.000  #           .
    relative_pitch_distribution[7] = 0.948  #           .
    relative_pitch_distribution[8] = 0.888  # relative pitch at x_R = 0.9
    # blade parameters
    lift_curve_slope = 0.97 # lift curve slope 
    root_thickness = 0.045  # blade thickness at root
    
    TT,TQ,TE = BEMT(J,P_D,BA_ratio,N_blades,chord,relative_pitch_distribution,
                        lift_curve_slope,root_thickness)
    assert np.allclose(TT,0.060,atol = 1e-3)
    assert np.allclose(TQ,0.011,atol = 1e-3)
    assert np.allclose(TE,0.595,atol = 1e-2)

def test_fortran3():
    J=0.3       #Advance Ratio
    P_D=0.6      # Pitch/Diamter Ratio
    BA_ratio=0.5   #  Blade Area Ratio
    N_blades=5    # number of blades
    
    #Set chord at spanwise positions
    chord = np.zeros(10)
    chord[1] = 0.208    # chord at x_R = 0.2
    chord[2] = 0.241    # chord at x_R = 0.3
    chord[3] = 0.263    #       .
    chord[4] = 0.276    #       .
    chord[5] = 0.279    #       .
    chord[6] = 0.269    #       .
    chord[7] = 0.241    #       .
    chord[8] = 0.184    # chord at x_R = 0.9
    # Pitch distribution along blade length relative to pitch at R = 0.7
    relative_pitch_distribution = np.zeros_like(chord)
    relative_pitch_distribution[1] = 0.888  # relative pitch at x_R = 0.2
    relative_pitch_distribution[2] = 1.008  # relative pitch at x_R = 0.3
    relative_pitch_distribution[3] = 1.055  #           .
    relative_pitch_distribution[4] = 1.060  #           .
    relative_pitch_distribution[5] = 1.039  #           .
    relative_pitch_distribution[6] = 1.000  #           .
    relative_pitch_distribution[7] = 0.948  #           .
    relative_pitch_distribution[8] = 0.888  # relative pitch at x_R = 0.9
    # blade parameters
    lift_curve_slope = 0.97 # lift curve slope 
    root_thickness = 0.045  # blade thickness at root
    
    TT,TQ,TE = BEMT(J,P_D,BA_ratio,N_blades,chord,relative_pitch_distribution,
                        lift_curve_slope,root_thickness)
    
    assert np.allclose(TT,0.175,atol = 1e-3)
    assert np.allclose(TQ,0.018,atol = 1e-3)
    assert np.allclose(TE,0.457,atol = 1e-2)
    
def test_fortran4():
    J=0.7      #Advance Ratio
    P_D=0.95      # Pitch/Diamter Ratio
    BA_ratio=0.4   #  Blade Area Ratio
    N_blades=4    # number of blades
    
    #Set chord at spanwise positions
    chord = np.zeros(10)
    chord[1] = 0.208    # chord at x_R = 0.2
    chord[2] = 0.241    # chord at x_R = 0.3
    chord[3] = 0.263    #       .
    chord[4] = 0.276    #       .
    chord[5] = 0.279    #       .
    chord[6] = 0.269    #       .
    chord[7] = 0.241    #       .
    chord[8] = 0.184    # chord at x_R = 0.9
    # Pitch distribution along blade length relative to pitch at R = 0.7
    relative_pitch_distribution = np.zeros_like(chord)
    relative_pitch_distribution[1] = 0.888  # relative pitch at x_R = 0.2
    relative_pitch_distribution[2] = 1.008  # relative pitch at x_R = 0.3
    relative_pitch_distribution[3] = 1.055  #           .
    relative_pitch_distribution[4] = 1.060  #           .
    relative_pitch_distribution[5] = 1.039  #           .
    relative_pitch_distribution[6] = 1.000  #           .
    relative_pitch_distribution[7] = 0.948  #           .
    relative_pitch_distribution[8] = 0.888  # relative pitch at x_R = 0.9
    # blade parameters
    lift_curve_slope = 0.97 # lift curve slope 
    root_thickness = 0.045  # blade thickness at root
    
    TT,TQ,TE = BEMT(J,P_D,BA_ratio,N_blades,chord,relative_pitch_distribution,
                        lift_curve_slope,root_thickness)
    
    assert np.allclose(TT,0.157,atol = 1e-3)
    assert np.allclose(TQ,0.025,atol = 1e-3)
    assert np.allclose(TE,0.704,atol = 1e-2)
