%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROPLIN3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Matlab implimentation of Blade Element Momentum Propeller Code PROPLIN3
%Orginal Fortran Code By Molland and Turnock
%Matlab Implementation AB Phillips 2010
%
%Modifications
%
%02/02/2012 Started Code Conversion
%03/02/2012 Code Validated against fortran version!
%
%Validation
%
%           Parmeters              |    Proplin3.f            |       Matlab      | 
%   J   P/D     BAR     NosBlades       KT      KQ     eta0    KT      KQ     eta0
%   0.8 0.95    0.4     4               0.113   0.019  0.746   0.1133  0.0193 0.7464  
%   0.7 0.8     0.8     3               0.060   0.011  0.595   0.0605  0.0113 0.5951
%   0.3 0.6     0.5     5               0.175   0.018  0.457   0.1749  0.0183 0.4575
%   0.7 0.95    0.4     4               0.157   0.025  0.704   0.1574  0.0249 0.7044  
t = cputime;
clear all
close all
t = cputime;
fprintf('\n')
fprintf('PROPLIN3\n')
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variable Details
%
%FORTRAN    MATLAB      DESCRIPTION
%
%P          P_D         Pitch/Diameter Ratio at 0.7 Radius
%J          J           Advance Ratio
%AR         BARatio     Blade BARatio Ratio
%RZ         NosBlades   Number of Blades



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0.8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Enter Propeller Parameters
J=0.8;          %Advance Ratio
P_D=0.95;        %Pitch/Diamter Ratio
BARatio=0.4;      %Blade Area Ratio
NosBlades=4;    %Number of Blades



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check Propeller PBARatioameters BARatioe Valid
if or(J>1.1, J<0.3)
    fprintf 'Advance Ratio - J outside limits of 0.3 to 1.1'
    return
else
    fprintf ('Advance Ratio - J=%1.2f \n', J)
end

if or(P_D<J+0.03,or(P_D<0.6,P_D>1.2))
    fprintf 'Pitch/Diameter Ratio at 0.7R is outside limits for this J value'
    return
else
    fprintf ('Pitch/Diameter Ratio at 0.7R - P/D=%1.2f\n', P_D)
end

if or(BARatio<0.4, BARatio>0.8)
    fprintf 'Blade Area Ratio - BAR outside limits of 0.4 to 0.8'
    return
else
    fprintf ('Blade Area Ratio - BAR=%1.2f \n', BARatio)
end

if or(NosBlades<3, NosBlades>5)
    fprintf 'Number of Blades - Number of Blades outside limits of 3 to 5'
    return
else
    fprintf ('Number of Blades - NosBlades=%i \n', NosBlades)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C2(2)=0.208;
C2(3)=0.241;
C2(4)=0.263;
C2(5)=0.276;
C2(6)=0.279;
C2(7)=0.269;
C2(8)=0.241;
C2(9)=0.184;


%Pitch Distribution along blade length relative to pitch distribution at
%R=0.7
PA(2)=0.888;
PA(3)=1.008;
PA(4)=1.055;
PA(5)=1.060;
PA(6)=1.039;
PA(7)=1.000;
PA(8)=0.948;
PA(9)=0.888;


VV=0;
LA=0.97;
RT=0.045;
LM=0.8*11.3;
LI=0.8*12.67;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Iteration Loop (for each blade element!)

for i=2:9  
    
    %assuming hub radius of 0.2R iterate through blade elements to 0.9R assume zero thrust and torque at tip
    
    C1(i)=C2(i)*BARatio*4.0/(NosBlades*0.5);    
    CP=NosBlades*C1(i);                         %Nc/D
    X(i)=i/10;    
    
    %x/R
    
    AI(i)=0;                                    %Assume initial Angle of Attack=0
    PR(i)=P_D*PA(i); 
    
    
    %Pitch/diameter ratio at blade element    
	TD=RT-(RT*0.935)*X(i);
	TC(i)=TD/C1(i);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sub loop iterates on angle of attack
    counter=1;
    SubConverged=0;
    while and(counter <2000, SubConverged==0)
        counter=counter+1;
        
        
        HA=atan(PR(i)/(pi*X(i)))/pi*180-AI(i);  %Calculate local hydrodynamic pitch angle - phi (deg)
        
        H1=tan(HA/180*pi);                      %tan phi
        AA=J/(pi*X(i));                         %tan psi (undisturbed flow angle?)
        
        EI=AA/H1;                               %Ideal Efficiency

        EA(i)=0.9*EI; 
        
        %Local Efficiency Estimate
        AE=0.0;                                 %tan(gamma) =atan(Cd(i)/Cl(i)) assume =0.0
        KF=X(i)*H1;                             %Lambda parameter for Goldstein correction

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculate Goldtein Correction, K
        SF=NosBlades/(2.0*KF)-0.5;
        F1=cosh(SF);
        F2=cosh(SF*X(i));
        F3=F2/F1;
        F4=acos(F3);
        K=2.0*F4/pi;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        GK(i)=K;                                %Store local Goldstien Correction to an array
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %subsubloop - iterates local efficency
        counter2=0;
        SubSubConverged=0;
        while and(counter2<2000, SubSubConverged==0)
            counter2=counter2+1;
           
            FA(i)=(1-EI)/(EI+AA*AA/EA(i)); %axial flow factor
            

            KT(i)=pi*J^2*X(i)*K*FA(i)*(1+FA(i));                %local dKT/dx
            

            AE=tan(AE/180*pi);                                  %tan(gamma)
            FT(i)=1-EI*(1+FA(i)); 
            %a' tangential flow factor
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate required lift coefficient - CL
            CL(i)=KT(i)*cos(HA/180*pi);
            CL(i)=CL(i)*4/(pi^2)/(X(i)^2);
           
            CL(i)=CL(i)/((1-FT(i))^2*(1-H1*AE)*CP);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate Drag Coefficient CD
            F6=0.0107+(AI(i)+1.0)*(-0.0015+AI(i)*(0.0015+0.000965*(AI(i)-1.0)));
            F7=-0.03833+(AI(i)+1.0)*(0.0133+AI(i)*(-0.015-0.01166*(AI(i)-1.0)));
            F8=0.8193+(AI(i)+1.0)*(-0.0138+AI(i)*(0.0903+0.079*(AI(i)-1.0)));
            F9=-3.076+(AI(i)+1.0)*(-0.0728+AI(i)*(-0.3162-0.2437*(AI(i)-1.0)));
            CD(i)=F6+TC(i)*(F7+(TC(i)-0.06)*(F8+F9*(TC(i)-0.12)));
            
            
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            AG=CD(i)/CL(i);                                     
            AE=atan(AG)/pi*180;  
            
            EZ=AA/tan((HA+AE)/180*pi);
             
            %%gamma
            if (abs(EZ-EA(i)) < 0.001)  
                %Initial estimate of local efficiency EA(I) is sufficiently close to calculated value
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Camber Correction, CC
                disp(EA(i))
                pause 
                K2(2)=1.0+2.857*(BARatio-0.4)*(BARatio-0.6)*(BARatio-0.9);
                
                K2(3)=1.19+(BARatio-0.4)*(BARatio-0.6)*(0.267+(BARatio-0.9)*0.1665);
	            K2(4)=1.3+(BARatio-0.4)*(0.1+(BARatio-0.6)*(1+(BARatio-0.9)*0.43));
	            K2(5)=1.54+(BARatio-0.4)*(0.15+(BARatio-0.6)*(1.1+0.286*(BARatio-0.9)));
	            K2(6)=1.67+(BARatio-0.4)*(0.55+(BARatio-0.6)*(0.1667+(BARatio-0.9)*2.9465));
	            K2(7)=1.8+(BARatio-0.4)*(0.75+(BARatio-0.6)*(0.3+(BARatio-0.9)*2.835));
	            K2(8)=1.8+(BARatio-0.4)*(1.0+(BARatio-0.6)*(1.333+(BARatio-0.9)*1.905));
	            K2(9)=1.75+(BARatio-0.4)*(1.25+(BARatio-0.6)*(1.5+(BARatio-0.9)*3.55));
                
	            U1=-0.65*KF*KF+1.1*KF+0.664;
	            U2=(0.85+(KF-0.3)*(-4.0+(KF-0.4)*(15.42-47.95*(KF-0.5))));
	            U2=-0.09+(KF-0.2)*U2;
	            U3=(1.375+(KF-0.3)*(-3.75+(KF-0.4)*(20.85-75.7875*(KF-0.5))));
	            U3=-0.2+(KF-0.2)*U3;
	            K1=U1+(BARatio-0.4)*(U2+U3*(BARatio-0.8));
	            CC(i)=K1*K2(i);
                
                if VV==1
                    MC(i)=MT(i)*TC(i);
                   
                    AC=(CC(i)*CL(i)/LI-MC(i))*LM/0.1097+(MC(i)*LI/LA);
                elseif VV==5
                    MC(i)=0.5*TC(i);
                    AC=(CC(i)*CL(i)/LI-MC(i))*LM/0.1097+(CL(i)/LA);
                    MT(i)=MC(i)/TC(i);
                    
                else
                    AC=CL(i)/LA;
                    
                    MC(i)=CL(i)/LI;
                    MC(i)=MC(i)*CC(i);
                     
                    
                    if (MC(i)<0.5*TC(i))
                        MT(i)=MC(i)/TC(i);
                        
                    else
                        VV=5;
                        MC(i)=0.5*TC(i);
                        AC=(CC(i)*CL(i)/LI-MC(i))*LM/0.1097+(CL(i)/LA);
                        MT(i)=MC(i)/TC(i);
                        
                        
                    end              
                
                end
                
               SubSubConverged=1;                                %Set subsub loop converged flag to 1
                
            else
                EA(i)=EZ;
                %update local efficiency estimate
            end          
        end %subsubloop
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        if (abs((AC-AI(i))/(AC))<0.1) 
            
            %subloop converged difference between angle of attack assumed AI and calculated angle of attack AC <0.1*AC
            KQ(i)=4.935*J*X(i)*X(i)*X(i)*K*FT(i)*(1+FA(i)); 
            %calculate dkq/dx
            
            SubConverged=1;                                      %set subloop converged flag to 1
        else
            AI(i)=(AC+AI(i))/2;
            disp(AI(i))
            pause%New Guess of Angle of Attack AI(i)
        end
    end %subloop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end     %End Main Loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assume 0 thrust and torque production at the tip
 X(10)=1.0;
 KT(10)=0.0;
 KQ(10)=0.0;
 FA(10)=0.0;
 FT(10)=0.0;
 GK(10)=0.0;
 
% MT(10)=MT(9);
% PR(10)=PR(9);
 AI(10)=AI(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Numerical Integration of Total Values
TT=(KT(2)+2.0*(KT(4)+KT(6)+KT(8))+4.0*(KT(3)+KT(5)+KT(7)+KT(9)));   %Integrate local dkt/dx along blade using simpsons rule
TT=0.1*TT/3.0;                                                      %Integrate local dkt/dx along blade using simpsons rule
TQ=(KQ(2)+2.0*(KQ(4)+KQ(6)+KQ(8))+4.0*(KQ(3)+KQ(5)+KQ(7)+KQ(9)));   %Integrate local dkq/dx along blade using simpsons rule
TQ=0.1*TQ/3.0;                                                      %Integrate local dkq/dx along blade using simpsons rule
TE=J*TT/(2.0*3.1416*TQ);                                            %Calculate open water efficiency
    
disp(TC(:))
fprintf ('Thrust Co-efficient, KT=%1.4f\n', TT)
fprintf ('Torque Co-efficient, KQ=%1.4f\n', TQ)
fprintf ('Open Water Efficiency, eta0=%1.4f\n', TE)

e = cputime - t



