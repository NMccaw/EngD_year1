PROGRAM Wind Turbine Performance Analysis by BEM Method
C     Author: J.Wellicome, modified by S.Turnock for sess6058 coursework
C       ******* Computes Power & thrust data and Radial Distributions
C               of Bending Moments for Turbines of specified geometry.
C                      WITH COMP CORRECTIONS.
C       **Declarations**
C
        INTEGER SECTIONS
        REAL A(20),CL(20),CD(20),SIG(20),TH(20),DIFCP(20),DIFCT(20)
        REAL ANGLE(20),LIFT(20),DRAG(20),INTCP(20),INTCT(20),INTCPX(20)
        REAL INTXCT(20),MT(20),MQ(20),MASS(20),TENSION(20),CENFUG(20)
        REAL MXX(20),MYY(20),MCRIT(20),MCRITA,MACH,NF,K,MU
        CHARACTER REPEAT
C
C       **Input and Setup**
C
        PI=3.14159265
        DTOR=0.017453292
        RTOD=1.0/DTOR
        VSOUND=342.0
C
      OPEN(UNIT=1,FILE='AERO.DAT',STATUS='OLD')
C       Open data file of Cl, Cd vs incidence
      OPEN(UNIT=2,FILE='BLADE.DAT',STATUS='OLD')
C       Open data file of Blade Geometry and mass
      OPEN(UNIT=3,FILE='OUTFILE.DAT',STATUS='OLD')
C       Open output file of detailed section performance estimates
      OPEN(UNIT=4,FILE='SUMFILE.DAT',STATUS='OLD')
C       Open output file of integrated turbine performance data.
C       All files must exist prior to running the program
C
C       Read two dimensional section aerodynamic data
C
C       Program assumes one data set applies for all sections
C       IAER=number of data points in set. IAER<=20.
C       MCRIT=critical Mach Number for section drag.
C       ALPHA=incidence angle in degrees
        WRITE(3,1015)
        WRITE(4,1015)
1015 FORMAT(//////,' CWIND TURBINE PERFORMANCE PROGRAM')
        WRITE(3,1020)
1020   FORMAT(/,'    ALPHA    CL      CD     MCRIT')
        READ(1,*)IAER
        DO 600 I=1,IAER
        READ(1,*)A(I),CL(I),CD(I),MCRIT(I)
        WRITE(3,1030)A(I),CL(I),CD(I),MCRIT(I)
1030   FORMAT(4F8.3)
        A(I)=DTOR*A(I)
600  CONTINUE
C
C       Read Blade Geometry
C
C       DIAM=rotor diameter (m)
C       SECTIONS=number of blade sections. SECTIONS<=20
C       XHUB=hub radius fraction. Sections assumed equally spaced from hub to tip.
C       B=number of blades.
C       For each section:
C       WIDTH=blade chord (m).
C       THETA=datum pitch angle to no-lift axis (degrees).
C       SECMASS= blade mass per unit span (Kg/m)
        WRITE(3,1040)
1040   FORMAT(/,'   X     SOLIDITY  BLADE ANGLE  BLADEMASS')
        READ(2,*)DIAM,B,SECTIONS,XHUB
      READ(2,*)VWIND,BPITCH
        DX=(1.0-XHUB)/SECTIONS
        X=XHUB
        DO 700 I=1,SECTIONS
        READ(2,*)WIDTH,THETA,SECMASS
        SIG(I)=B*WIDTH
        TH(I)=THETA*DTOR
        MASS(I)=4.0*SECMASS/DIAM/DIAM
        WRITE(3,1050)X,SIG(I),RTOD*TH(I),MASS(I)
1050   FORMAT(F6.2,F9.4,2F12.4)
        X=X+DX
700  CONTINUE
C
C       Begining of performance loop for one case. New cases are run
C       on answering 'Y'to 'ANOTHER CASE?'
C
C       BPITCH=blade rotation angle from datum pitch (degrees)
C       TSR=rotor tipspeed to windspeed ratio
1    WRITE(*,1010)
1010 FORMAT(' WINDSPEED,RPM,PITCH ANGLE=')
        READ(*,*)VWIND,RPM,BPITCH
        TSR=DIAM*PI*RPM/(VWIND*60.0)
        BPITCH=DTOR*BPITCH
C
C       Copy case data to output files OUTFILE and SUMFILE
C
        WRITE(3,1051)DIAM
        WRITE(4,1051)DIAM
1051 FORMAT(/,'        DIAMETER=',F8.4,'  m')
        WRITE(3,1052)VWIND
        WRITE(4,1052)VWIND
        WRITE(*,1052)VWIND
1052 FORMAT('       WINDSPEED=',F8.4,'  m/s')
        WRITE(3,1053)RPM
        WRITE(4,1053)RPM
        WRITE(*,1053)RPM
1053 FORMAT('             RPM=',F8.2)
        WRITE(3,1054)TSR
1054 FORMAT(' TIP SPEED RATIO=',F8.3)
        WRITE(3,1056)RTOD*BPITCH
        WRITE(4,1056)RTOD*BPITCH
        WRITE(*,1056)RTOD*BPITCH
1056 FORMAT('     PITCH ANGLE=',F8.3,'  deg')
C
        X=XHUB-DX
C
C       ** X loop - All sections:Hub to Tip **
C
C       Section performance estimates using Blade Element-Momentum theory
C
C       Ouput column headings to OUTFILE
        WRITE(3,1060)
1060 FORMAT(/,'    X    PHI   A-AXIAL  A-TANG   GOLD   MACH    MCRIT')
C
C       Set initial relaxation factor (MU) for updating ALPHA in iteration loop
        MU=1.0
C       IMCRIT=Count of secs MACH>MCRITA
        IMCRIT=0
C
        DO 300 IXCOUNT=1,SECTIONS
C
        X=X+DX
C       X=Next Radius Fraction
        THETA=TH(IXCOUNT)+BPITCH
        SIGMA=SIG(IXCOUNT)
        ICOUNT=0
        MACH=VWIND*SQRT(1.0+(TSR*X)**2)/VSOUND
        IF (MACH.LE.0.97) THEN
        COMPFAC=1.0/(SQRT(1.0-MACH*MACH))
        ELSE
        COPMFAC=4.0
        END IF
C       COPFAC=Glauert compressibiity correction factor for Cl,Cd.
C       Valid for moderate Mach Nos.
C
C       ** Alpha loop start **
C
        ALPHA=ATAN(0.8/(TSR*X))-THETA
C       First approximation for ALPHA
C
100    PHI=ALPHA+THETA
C       PHI Calculated from blade angle and angle of attack
        TPHI=TAN(PHI)
        CLIFT=COMPFAC*COEFT(IAER,A,CL,ALPHA)
        CDRAG=COMPFAC*COEFT(IAER,A,CD,ALPHA)
        GAMMA=ATAN(CDRAG/CLIFT)
        K=GOLDK(B,X,TPHI)
C       K=Goldstein    momentum averaging (tip loss) factor-see subroutine
        AA=4.0*PI*X*K*SIN(PHI)*TPHI/(1.0+TPHI*CDRAG/CLIFT)
        AB=SIGMA*CLIFT
        AA=AB/(AB+AA)
C       AA=Axial Inflow Factor
        AT=AA*TAN(PHI-GAMMA)/TSR/X
C       AT=Tangential Inflow Factor
        DCPDX=2.0*X*TSR*SIGMA*CLIFT*(1.0-TAN(GAMMA)/TPHI)/SIN(PHI)/PI
        AB=1.0-AA
        DCPDX=DCPDX*AB*AB
C       DCPDX=Derivative of Power Coefficient calculated from element properties
        PHI=ATAN((1.0-AA)/(1.0+AT)/TSR/X)
C       PHI Calculated from velocity vector diagram
        A0=PHI-THETA
C       A0=Angle of attack derived from velocity vector diagram
        ICOUNT=ICOUNT+1
        DIFF=ABS(A0-ALPHA)
        A2=A1
        A1=ALPHA
        ALPHA=(1.0-MU)*ALPHA+MU*A0
        IF ((ICOUNT.EQ.5).OR.(ICOUNT.EQ.10)) THEN
C       Apply Aitken Extrapolator to esimate a likely final ALPHA value as an aid
C       to selecting a suitable relaxation factor. See notes.
C       MU re-estimated only at 5th and 10th iterations
        A3=(ALPHA*A2-A1*A1)/(ALPHA-2.0*A1+A2)
C       Reset relaxation factor, limiting the allowable range of values.
        AMU=(A3-A1)/(A0-A1)
        IF (AMU.GT.2.0) AMU=2.0
        IF (AMU.LT.0.4) AMU=0.4
        MU=AMU
C       Reset next ALPHA value
        ALPHA=A3
        END IF
C       Updated angle of attack
        IF((DIFF.GE.0.00002).AND.(ICOUNT.LT.25)) GOTO 100
C       Iterate until alpha within about 0.001 deg or max 25 iterations
C
C       ** End of alpha loop **
C
        DCTDX=8.0*X*K*AA*(1.0-AA)
C       DCTDX=Derivative of Thrust coefficient from momentum
C
C       Store section aerodynamic performance estimates
        DIFCP(IXCOUNT)=DCPDX
        DIFCT(IXCOUNT)=DCTDX
        ANGLE(IXCOUNT)=ALPHA
        LIFT(IXCOUNT)=CLIFT
        MCRITA=COEFT(IAER,A,MCRIT,ALPHA)
        IF (MACH.GT.MCRITA) IMCRIT=IMCRIT+1
        DRAG(IXCOUNT)=CDRAG
C
C       Write section performance estimates to OUTFILE
        WRITE(3,1000)X,RTOD*PHI,AA,AT,K,MACH,MCRITA
1000 FORMAT(F5.2,F8.2,5F8.4)
300    CONTINUE
C
        WRITE(3,3000)
3000   FORMAT(/,'   X    ALPHA     CL      CD     DCPDX   DCTDX ')
        X=XHUB
        DO 800 IX=1,SECTIONS
        AINC=RTOD*ANGLE(IX)
        WRITE(3,3010)X,AINC,LIFT(IX),DRAG(IX),DIFCP(IX),DIFCT(IX)
3010 FORMAT(F5.2,F8.2,4F8.4)
        X=X+DX
800  CONTINUE
C
C       **  End of xloop **
C
C       ** Start of Radial Integration routine **
C
C       Trapeziodal integration of power coeft,thrust coeft and bending moment data
C
        X=1.0-DX
        INTCP(SECTIONS)=DIFCP(SECTIONS)
        INTCT(SECTIONS)=DIFCT(SECTIONS)
        INTCPX(SECTIONS)=DIFCP(SECTIONS)/X
        INTXCT(SECTIONS)=DIFCT(SECTIONS)*X
        CENFUG(SECTIONS)=X*MASS(SECTIONS)
        IXCOUNT=SECTIONS
C
400  IPLUS=IXCOUNT
      IXCOUNT=IXCOUNT-1
      INTCP(IXCOUNT)=INTCP(IPLUS)+DIFCP(IPLUS)+DIFCP(IXCOUNT)
      INTCT(IXCOUNT)=INTCT(IPLUS)+DIFCT(IPLUS)+DIFCT(IXCOUNT)
      INTCPX(IXCOUNT)=INTCPX(IPLUS)+DIFCP(IPLUS)/X+DIFCP(IXCOUNT)/(X-DX)
      INTXCT(IXCOUNT)=INTXCT(IPLUS)+DIFCT(IPLUS)*X+DIFCT(IXCOUNT)*(X-DX)
      CENFUG(IXCOUNT)=CENFUG(IPLUS)+MASS(IPLUS)*X+MASS(IXCOUNT)*(X-DX)
     X=X-DX

      IF(IXCOUNT.GT.1) GOTO 400
C
C       Convert integrals to thrust and torque moments
        X=XHUB
C
        WRITE(3,3100)
3100 FORMAT(/,'    X      MOM-XX     MOM-YY    TENSION')
        DO 500 IXCOUNT=1,SECTIONS
C
        MT(IXCOUNT)=(INTXCT(IXCOUNT)-X*INTCT(IXCOUNT))*DX/2.0
        MQ(IXCOUNT)=(INTCP(IXCOUNT)-X*INTCPX(IXCOUNT))*DX/2.0
        CS=COS(TH(IXCOUNT))
        SN=SIN(TH(IXCOUNT))
        MXX(IXCOUNT)=MT(IXCOUNT)*CS+MQ(IXCOUNT)*SN
        MYY(IXCOUNT)=MQ(IXCOUNT)*CS-MT(IXCOUNT)*SN
        TENSION(IXCOUNT)=CENFUG(IXCOUNT)*DX/2.0
        WRITE(3,3110)X,MXX(IXCOUNT),MYY(IXCOUNT),TENSION(IXCOUNT)
3110 FORMAT(F8.4,2F10.6,F10.4)
        X=X+DX
C
500    CONTINUE
C
C       Set total Power and Thrust Coefficients
        CPOWER=INTCP(1)*DX/2.0
        CTHRUST=INTCT(1)*DX/2.0
C
C       ** End of Radial Integration **
C
C       ** Output Non-dimensional Data to OUTFILE **
C
        WRITE(3,3200)
3200 FORMAT(/,' OVERALL NON-DIMENSIONAL PARAMETERS')
        WRITE(3,3210)CPOWER
3210 FORMAT(/,'             POWER COEFFICIENT=',F8.5)
        WRITE(3,3220)CTHRUST
3220 FORMAT('            THRUST COEFFICIENT=',F8.5)
        WRITE(3,3230)MT(1)
3230 FORMAT(' HUB THRUST MOMENT COEFFICIENT=',F8.5)
        WRITE(3,3240)MQ(1)
3240 FORMAT(' HUB TORQUE MOMENT COEFFICIENT=',F8.5)
        WRITE(3,3245)TENSION(1)
3245 FORMAT('       HUB TENSION COEFFICIENT=',F8.5)
C
C       ** Output Real Values to OUTFILE and to SUMFILE **
C
        FCTOR1=1.23*PI*DIAM*DIAM*VWIND*VWIND/8.0
        FCTOR2=FCTOR1*DIAM/2.0
        WRITE(3,3250)
3250 FORMAT(/,' OVERALL PERFORMANCE DATA')
        WRITE(3,3260)CPOWER*FCTOR1*VWIND/1000.0
        WRITE(4,3260)CPOWER*FCTOR1*VWIND/1000.0
        WRITE(*,3260)CPOWER*FCTOR1*VWIND/1000.0
3260 FORMAT(/,'    ABSORBED POWER=',F10.4,'  KW')
        WRITE(3,3270)CTHRUST*FCTOR1
        WRITE(4,3270)CTHRUST*FCTOR1
        WRITE(*,3270)CTHRUST*FCTOR1
3270 FORMAT('            THRUST=',F10.4,'  N')
        WRITE(3,3280)MXX(1)*FCTOR2
        WRITE(4,3280)MXX(1)*FCTOR2
        WRITE(*,3280)MXX(1)*FCTOR2
3280 FORMAT('     HUB XX MOMENT=',F10.4,'  Nm')
        WRITE(3,3290)MYY(1)*FCTOR2/TSR
        WRITE(4,3290)MYY(1)*FCTOR2/TSR
        WRITE(*,3290)MYY(1)*FCTOR2/TSR
3290 FORMAT('     HUB YY MOMENT=',F10.4,'  Nm')
        WRITE(3,3295)TENSION(1)*TSR*TSR*VWIND*VWIND*DIAM/2000.0
        WRITE(4,3295)TENSION(1)*TSR*TSR*VWIND*VWIND*DIAM/2000.0
        WRITE(*,3295)TENSION(1)*TSR*TSR*VWIND*VWIND*DIAM/2000.0
3295 FORMAT('  HUB TENSILE LOAD=',F10.4,'  kN')
        IF (IMCRIT.GT.0)THEN
        WRITE(3,3297)IMCRIT
        WRITE(4,3297)IMCRIT
        WRITE(*,3297)IMCRIT
        END IF
3297 FORMAT(/,' MACH>CRITICAL FOR',I3,'  SECTIONS')
        WRITE(*,3300)
C
C       Select another case or close program
3300 FORMAT(/,' ANOTHER CASE ? (Y/N) ')
3310   READ(*,'(A)')REPEAT
        IF ((REPEAT.NE.'Y').AND.(REPEAT.NE.'N')) GOTO 3310
        IF (REPEAT.EQ.'Y') GOTO 1
       CLOSE(1)
        CLOSE(2)
        CLOSE(3)
C
        STOP
        END
C
C       ******* END OF MAIN PROGRAM *********
C
        REAL FUNCTION COEFT(N,A,C,ALPHA)
C
        INTEGER N,M
        REAL ALPHA,U,A(20),C(20),S(10),V(10)
C
C       Data C(A) defined by N element arrays A,C
C       COEFT=Lagrange M point interpolation at angle ALPHA
C       M=EVEN-Set initially at M=4
C
        M=4
        DTOR=0.017453292
        RTOD=1.0/DTOR
C
C       Scan for interval along curve
        I=0
C
10   I=I+1
        IF ((A(I).LE.ALPHA).AND.(I.LT.N)) GOTO 10
C       Set end value if ALPHA<=A(1) or >=A(N) else interpolate
        IF ((I.EQ.1).OR.(ALPHA.GE.A(N))) THEN
        COEFT=C(I)
        ELSE
C       ** Lagrange interpolation routine-see notes **
C
C       ** Generate S,V arrays **
        I=I-M/2-1
        IF(I.LT.0) I=0
        IF(I.GT.(N-M)) I=N-M
C
        DO 20 J=1,M
        V(J)=C(J+I)
        S(J)=A(J+I)
20   CONTINUE
C
        DO 30 J=1,(M-1)
        DO 40 I=J+1,M
        V(I)=(V(I)-V(J))/(S(I)-S(J))
  40  CONTINUE
  30  CONTINUE
C
C       ** Evaluate COEFT for given ALPHA **
C
        U=V(M)
        I=M-1
50   U=(ALPHA-S(I))*U+V(I)
        I=I-1
        IF(I.GT.0) GOTO 50
C
        COEFT=U
C
        END IF
        END
C
C       ******* END OF COEFT ********
C
        REAL FUNCTION GOLDK(B,X,TPHI)
C       B=number of blades  X=radius fraction  TPHI=tan(phi)
        REAL F,PI,B,X,TPHI
C
        PI=3.14159265
        IF (ABS(TPHI).LT.0.001) THEN
        GOLDK=1.0
        ELSE
        F=B/2.0/X/TPHI-0.50
        IF (ABS(F).GT.85.0) THEN
        GOLDK=1.0
        ELSE
        GOLDK=2.0*ACOS(COSH(X*F)/COSH(F))/PI
        END IF
        END IF
C
        END
C       *****END OF GOLDK *****
