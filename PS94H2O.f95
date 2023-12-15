!	Code "borrowed" from TD website.
!	H2O tested against Thermocalc output and gives exact value at 10000 C, 20,000 bars
!
!===============================================================
!===============================================================
!--------------------------------
!--------------------------------
      SUBROUTINE PS94H2O(P,T,GRE,VRE)
      IMPLICIT NONE
!-- Subroutines by CdC and Erik Düsterhoeft
!-- K.S. Pitzer and S.M. Sterner, 1994
!-- J. Chem. Phys. 101(4),1994,3111-3116
!-- critical point H2O: 373.946 C, 220.64 Bar
!-- critical density H2O: 0.017873993893977 mol/ccm
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
!--
      REAL*8 T,P,R,ZZ,MI,MA,LIM1,LIM2,RHINI, &
      VOL1,VOL2,VOL0,VDP1,VDP2,PDV,F1,RH0,RH1,RH2,VRE,GRE
      REAL*8 ARES,LOF
      INTEGER*4 I,I0,I1,I2,IM,HOW
!--
      REAL*8 GF,K1,K2,K3,K4,CPRDT,CPRTDT,H0,S0,TT,SQT,TT0,SQT0,T0
!-H2O
      DATA H0,S0/-241.81D0,188.80D-3/
      DATA K1,K2,K3,K4/0.0401D0,0.8656D-5,487.5D0,-0.2512D0/
      DATA T0,TT0/298.15D0,88893.4225D0/
!-CO2
!      DATA H0,S0/-393.51D0,213.70D-3/
!      DATA K1,K2,K3,K4/0.0878D0,-0.2644D-5,706.4D0,-0.9989D0/
!      DATA T0,TT0/298.15D0,88893.4225D0/
!--
      R=8.3143D0
      SQT0=DSQRT(T0)
      TT=T*T
      SQT=DSQRT(T)
      RT=R*T
!
      CPRDT=K1*(T-T0)+K2*(TT-TT0)/2.0D0 &
      -K3*(1.0D0/T-1.0D0/T0)+2D0*K4*(SQT-SQT0)
      CPRTDT=K1*DLOG(T/T0)+K2*(T-T0) &
      -K3*(1.0D0/(TT)-1.0D0/(TT0))/2.0D0 &
      -2D0*K4*(1.0D0/SQT-1.0D0/SQT0)
      GF=1.0D3*(H0+CPRDT-T*(S0+CPRTDT))
!
      I0=0
      I1=0
      I2=0
      VRE=0.0D0
      GRE=0.0D0
!----
      CALL PISTH2ODEF(T)
!      CALL PISTCO2DEF(T)
!=====
!----- (0) low RH at 1 Bar. If fail, use ideal gas
      ZZ=1.0D0/RT/10.0D0
      LIM1=-1D20
      LIM2=0.03D0
      RHINI=1.0D0/RT
      CALL PISTFINDRH(ZZ,LIM1,LIM2,RHINI,RH0)
      IF (RH0.EQ.0.0D0) THEN
      VOL0=RT/P
      ELSE
      VOL0=1.0D0/RH0/10.0D0
      I0=1
      END IF
!
!----- (1) try low RH If fail, use ideal gas
      ZZ=P/RT/10.0D0
      LIM1=-1D20
      LIM2=0.03D0
      RHINI=1D-5
      CALL PISTFINDRH(ZZ,LIM1,LIM2,RHINI,RH1)
      IF (RH1.EQ.0.0D0) THEN
      VOL1=RT/P
      RH1=1.0D0/VOL1/10.0D0
      ELSE
      VOL1=1.0D0/RH1/10.0D0
      I1=1
      END IF
!
!----- (2) try high RH If fail, use ideal gas
!      ZZ=P/RT/10.0D0
      LIM1=0.01D0
      LIM2=1D20
      RHINI=6D-2
      CALL PISTFINDRH(ZZ,LIM1,LIM2,RHINI,RH2)
      IF (RH2.EQ.0.0D0) THEN
      VOL2=RT/P
      RH2=1.0D0/VOL2/10.0D0
      ELSE
      VOL2=1.0D0/RH2/10.0D0
      I2=1
      END IF
!
!===== if all fails (RH1=0 and RH2=0) try interval
      IF (I1.EQ.0.AND.I2.EQ.0) THEN
      RH2=0.0D0
      ZZ=P/R/T/10.0D0
      I=0
      MI=1D-7
      MA=10.0D0
      CALL INTERV(ZZ,MI,MA,I,RH2)
      I2=1
      VOL2=1.0D0/RH2/10.0D0
      END IF
!=====
      VDP1=0.0D0
      VDP2=0.0D0
!
      HOW=1
!==========================================================
!==========================================================
!----- IF HOW=1: brute force method
!==========================================================
!==========================================================
      IF (HOW.EQ.1) THEN
!=====
      IF (I1.GT.0) THEN
      CALL INTEGBET(VOL1,VOL0,PDV,I,IM)
      F1=-(VOL0-VOL1)+(P-1.0D0)*VOL1
      VDP1=PDV+F1+GF
      END IF
!
      IF (I2.GT.0) THEN
      CALL INTEGBET(VOL2,VOL0,PDV,I,IM)
      F1=-(VOL0-VOL2)+(P-1.0D0)*VOL2
      VDP2=PDV+F1+GF
      END IF
!=====
      END IF
!==========================================================
!==========================================================
!----- IF HOW=2: calulate ln(f), gas equation (no liquid)
!----- assumes that standard state and Cp are for ideal gas
!==========================================================
!==========================================================
      IF (HOW.EQ.2) THEN
!=====
!
      IF (I1.GT.0) THEN
      CALL PISTARES(RH1,ARES)
      LOF=(DLOG(RH2*10.0D0)+ARES+P/RH2/RT/10.0D0)+DLOG(RT)-1.0D0
      VDP1=GF+RT*LOF
      END IF
!
      IF (I2.GT.0) THEN
      CALL PISTARES(RH2,ARES)
      LOF=(DLOG(RH2*10.0D0)+ARES+P/RH2/RT/10.0D0)+DLOG(RT)-1.0D0
      VDP2=GF+RT*LOF
      END IF
!=====
      END IF
!==========================================================
!==========================================================
      IF (VDP1.LT.VDP2) THEN
      VRE=VOL1
      GRE=VDP1
      ELSE
      VRE=VOL2
      GRE=VDP2
      END IF
!-----
      RETURN
      END
!-----
!--------------------------------
!--------------------------------
      SUBROUTINE PS94CO2(P,T,GRE,VRE)
      IMPLICIT NONE
!-- Subroutines by CdC and Erik Düsterhoeft
!-- K.S. Pitzer and S.M. Sterner, 1994
!-- J. Chem. Phys. 101(4),1994,3111-3116
!-- critical point CO2: 31.04 C, 73.8 Bar
!-- critical density CO2: 0.010656668938878 mol/ccm
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
!--
      REAL*8 T,P,R,ZZ,MI,MA,LIM1,LIM2,RHINI, &
      VOL1,VOL2,VOL0,VDP1,VDP2,PDV,F1,RH0,RH1,RH2,VRE,GRE
      REAL*8 ARES,LOF
      INTEGER*4 I,I0,I1,I2,IM,HOW
!--
      REAL*8 GF,K1,K2,K3,K4,CPRDT,CPRTDT,H0,S0,TT,SQT,TT0,SQT0,T0
!-H2O
!      DATA H0,S0/-241.81D0,188.80D-3/
!      DATA K1,K2,K3,K4/0.0401D0,0.8656D-5,487.5D0,-0.2512D0/
!      DATA T0,TT0/298.15D0,88893.4225D0/
!-CO2
      DATA H0,S0/-393.51D0,213.70D-3/
      DATA K1,K2,K3,K4/0.0878D0,-0.2644D-5,706.4D0,-0.9989D0/
      DATA T0,TT0/298.15D0,88893.4225D0/
!--
      R=8.3143D0
      SQT0=DSQRT(T0)
      TT=T*T
      SQT=DSQRT(T)
      RT=R*T
!
      CPRDT=K1*(T-T0)+K2*(TT-TT0)/2.0D0 &
      -K3*(1.0D0/T-1.0D0/T0)+2D0*K4*(SQT-SQT0)
      CPRTDT=K1*DLOG(T/T0)+K2*(T-T0) &
      -K3*(1.0D0/(TT)-1.0D0/(TT0))/2.0D0 &
      -2D0*K4*(1.0D0/SQT-1.0D0/SQT0)
      GF=1.0D3*(H0+CPRDT-T*(S0+CPRTDT))
!
      I0=0
      I1=0
      I2=0
!----
!      CALL PISTH2ODEF(T)
      CALL PISTCO2DEF(T)
!
!----- (0) low RH at 1 Bar. If fail, use ideal gas
      ZZ=1.0D0/RT/10.0D0
      LIM1=-1D20
      LIM2=0.03D0
      RHINI=1.0D0/RT
      CALL PISTFINDRH(ZZ,LIM1,LIM2,RHINI,RH0)
      IF (RH0.EQ.0.0D0) THEN
      VOL0=RT/P
      ELSE
      VOL0=1.0D0/RH0/10.0D0
      I0=1
      END IF
!
!----- (1) try low RH If fail, use ideal gas
      ZZ=P/RT/10.0D0
      LIM1=-1D20
      LIM2=0.03D0
      RHINI=1D-5
      CALL PISTFINDRH(ZZ,LIM1,LIM2,RHINI,RH1)
      IF (RH1.EQ.0.0D0) THEN
      VOL1=RT/P
      RH1=1.0D0/VOL1/10.0D0
      ELSE
      VOL1=1.0D0/RH1/10.0D0
      I1=1
      END IF
!
!----- (2) try high RH If fail, use ideal gas
!      ZZ=P/RT/10.0D0
      LIM1=0.01D0
      LIM2=1D20
      RHINI=6D-2
      CALL PISTFINDRH(ZZ,LIM1,LIM2,RHINI,RH2)
      IF (RH2.EQ.0.0D0) THEN
      VOL2=RT/P
      RH2=1.0D0/VOL2/10.0D0
      ELSE
      VOL2=1.0D0/RH2/10.0D0
      I2=1
      END IF
!
!===== if all fails try interval
      IF (I1.EQ.0.AND.I2.EQ.0) THEN
      RH2=0.0D0
      ZZ=P/R/T/10.0D0
      I=0
      MI=1D-7
      MA=10.0D0
      CALL INTERV(ZZ,MI,MA,I,RH2)
      I2=1
      VOL2=1.0D0/RH2/10.0D0
      END IF
!=====
!=====
      VDP1=0.0D0
      VDP2=0.0D0
!
      HOW=1
!==========================================================
!==========================================================
!----- IF HOW=1: brute force method
!==========================================================
!==========================================================
      IF (HOW.EQ.1) THEN
!=====
      IF (I1.GT.0) THEN
      CALL INTEGBET(VOL1,VOL0,PDV,I,IM)
      F1=-(VOL0-VOL1)+(P-1.0D0)*VOL1
      VDP1=PDV+F1+GF
      END IF
!
      IF (I2.GT.0) THEN
      CALL INTEGBET(VOL2,VOL0,PDV,I,IM)
      F1=-(VOL0-VOL2)+(P-1.0D0)*VOL2
      VDP2=PDV+F1+GF
      END IF
!=====
      END IF
!==========================================================
!==========================================================
!----- IF HOW=2: calulate ln(f), gas equation (no liquid)
!----- assumes that standard state and Cp are for ideal gas
!==========================================================
!==========================================================
      IF (HOW.EQ.2) THEN
!=====
!
      IF (I1.GT.0) THEN
      CALL PISTARES(RH1,ARES)
      LOF=(DLOG(RH2*10.0D0)+ARES+P/RH2/RT/10.0D0)+DLOG(RT)-1.0D0
      VDP1=GF+RT*LOF
      END IF
!
      IF (I2.GT.0) THEN
      CALL PISTARES(RH2,ARES)
      LOF=(DLOG(RH2*10.0D0)+ARES+P/RH2/RT/10.0D0)+DLOG(RT)-1.0D0
      VDP2=GF+RT*LOF
      END IF
!=====
      END IF
!==========================================================
!==========================================================
      IF (VDP1.LT.VDP2) THEN
      VRE=VOL1
      GRE=VDP1
      ELSE
      VRE=VOL2
      GRE=VDP2
      END IF
!=====
!-----
      RETURN
      END
!-----
!--------------------------------
!--------------------------------
      SUBROUTINE INTERV(ZZ,MI,MA,I,RH)
      IMPLICIT NONE
      REAL*8 ZZ,MI,MA,RH,RH1,RH2,RH3,Z1,Z2,Z3,IST,DX,X,Y
      INTEGER*4 I
!----
      DX=(MA-MI)/10000.0D0
      X=MI
      DO WHILE(X.LE.MA)
       CALL PISTZ(X,Y)
      X=X+DX
      END DO
!-----
      I=0
      RH1=MI
      CALL PISTZ(RH1,Z1)
      RH2=MA
      CALL PISTZ(RH2,Z2)
      IST=0.0D0
      IF (Z1.LT.ZZ.AND.Z2.GT.ZZ) IST=1.0D0
      IF (Z1.GT.ZZ.AND.Z2.LT.ZZ) IST=-1.0D0
      IF (IST.EQ.0.0D0) THEN
       RH=0.0D0
       RETURN
      END IF
!----
      DO WHILE (DABS(Z2-Z1).GT.ZZ/1D6)
      RH=(RH1+RH2)/2.0D0
!----
       I=I+1
       IF (I.GT.100) RETURN
       RH3=(RH1+RH2)/2.0D0
       CALL PISTZ(RH3,Z3)
        IF ((Z3-ZZ)*IST.GT.0.0D0) THEN
         RH2=RH3
         Z2=Z3
        ELSE
         RH1=RH3
         Z1=Z3
        END IF
      END DO
!----
      RETURN
      END
!----
!--------------------------------
!--------------------------------
      SUBROUTINE PISTFINDRH(ZZ,LIM1,LIM2,RH1,RH)
      IMPLICIT NONE
!--
!-- initial guess for 1 Bar, and gas: V=RT, RH=1/RT
!-- falls Ableitung negativ: fuer gas:RH=RH/2
!-- optimistic limits: 10 > RH > 10-7
!-- LIM1: lower limit OF RH, LIM2: upper limit OF RH
!-- RH1 : initial; guess
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
      REAL*8 RH,RH1,RH2,Z1,Z2,ZS,DX,ZZ,DXTEST,ZS2,LIM1,LIM2
      INTEGER*4 I
!-----
      DX=100.0D0
      DXTEST=ZZ/1D6
      CALL PISTZ(RH1,Z1)
      I=0
      DO WHILE(DX.GT.DXTEST.AND.I.LT.100)
      I=I+1
       IF (DX.LT.DXTEST) GOTO 90
       CALL PISTZS(RH1,ZS)
       RH2=RH1+(ZZ-Z1)/ZS
       CALL PISTZS(RH2,ZS2)
       IF (ZS2.LT.0.0D0.OR.RH2.LT.LIM1.OR.RH2.GT.LIM2) THEN
        RH=0.0D0
        RETURN
       END IF
       DX=DABS(ZZ-Z2)
       CALL PISTZ(RH2,Z2)
       DX=DABS(ZZ-Z2)
       RH1=RH2
       Z1=Z2
      END DO
   90 CONTINUE
      RH=RH1
!=====
!----
      RETURN
      END
!----
!--------------------------------
!--------------------------------
      SUBROUTINE PISTARES(RH,ARES)
      IMPLICIT NONE
!--
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
      REAL*8 RH,RH2,RH3,RH4,F1,ARES
!-----
      RH2=RH*RH
      RH3=RH2*RH
      RH4=RH3*RH
      F1=CC(2)+CC(3)*RH+CC(4)*RH2+CC(5)*RH3+CC(6)*RH4
      ARES=CC(1)*RH+(1.0D0/F1-1.0D0/CC(2)) &
      -(CC(7)/CC(8))*(DEXP(-CC(8)*RH)-1.0D0) &
      -(CC(9)/CC(10))*(DEXP(-CC(10)*RH)-1.0D0)
!=====
!----
      RETURN
      END
!----
!--------------------------------
!--------------------------------
      SUBROUTINE PISTZ(RH,Z)
      IMPLICIT NONE
!--
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
      REAL*8 Z,RH,RH2,RH3,RH4,F1,F2
!-----
      RH2=RH*RH
      RH3=RH2*RH
      RH4=RH3*RH
      F1=CC(3)+2.0D0*CC(4)*RH+3.0D0*CC(5)*RH2+4.0D0*CC(6)*RH3
      F2=CC(2)+CC(3)*RH+CC(4)*RH2+CC(5)*RH3+CC(6)*RH4
      Z=RH+CC(1)*RH2-RH2*(F1/F2**2)+CC(7)*RH2*DEXP(-CC(8)*RH)+ &
      CC(9)*RH2*DEXP(-CC(10)*RH)
!=====
!----
      RETURN
      END
!----
!--------------------------------
!--------------------------------
      SUBROUTINE INTEGBET(VMIN,VMAX,G,I,IM)
      IMPLICIT NONE
!--
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
      REAL*8 VMIN,VMAX,G,DX,DXMIN,XU,XO,YU,YO, &
      X(100),Y(100),YPREC,F1,XD
      INTEGER*4 I,IMEGA,II,IC,IM
!----- primary interval e.g. (max-min)/100 DX
!----- minimal interval e.g. (max-min)/10000 DXMIN
!----- Y precision e.g. 1D-6
!----- COMMENT Apr. 2014: DXMIN=(VMAX-VMIN)/1.0D6 is slower but slightly better
!----- COMMENT Apr. 2014: DXMIN=(VMAX-VMIN)/1.0D4 is faster but not good at very high pressures
!----- X(10) and Y(10) probably sufficient
      G=0.0D0
      DX=(VMAX-VMIN)/100.0D0
      DXMIN=(VMAX-VMIN)/1.0D5
      YPREC=1D-5
!-- outer loop
      XU=VMIN
      CALL PISTPV(XU,YU)
      I=0
!-- begin outer loop
      DO IMEGA=1,100
      XO=XU+DX
      CALL PISTPV(XO,YO)
!====
      X(1)=XU
      Y(1)=YU
      X(2)=XO
      Y(2)=YO
      IC=2
      IM=2
!-- begin inner loop
      DO WHILE (IC.GT.1)
!-add a point
      DO II=IC,2,-1
      X(II+1)=X(II)
      Y(II+1)=Y(II)
      END DO
      X(2)=(X(1)+X(3))/2.0D0
      CALL PISTPV(X(2),Y(2))
      IC=IC+1
      I=I+1
      IF (IC.GT.IM) IM=IC
      XD=X(2)-X(1)
      F1=Y(2)-(Y(1)+Y(3))/2.0D0
      IF (DABS(F1).LT.YPREC.OR.DABS(XD).LT.DXMIN) THEN
!
!      G=G+XD*(Y(1)+Y(2))/2.0D0
!      G=G+XD*(Y(2)+Y(3))/2.0D0
!Simpson
       G=G+(X(3)-X(1))/6.0D0*(Y(1)+4.0D0*Y(2)+Y(3))
!-substract 2 points
      DO II=1,IC-2
      X(II)=X(II+2)
      Y(II)=Y(II+2)
      END DO
      IC=IC-2
      END IF
!-- end inner loop: back to add a point
      END DO
 !====
       XU=XO
       YU=YO
 !-- end outer loop
     END DO
!=====
!----
      RETURN
      END
!----
!--------------------------------
!--------------------------------
      SUBROUTINE PISTPV(V,P)
      IMPLICIT NONE
!--
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
      REAL*8 Z,RH,RH2,RH3,RH4,F1,F2,V,P
!-----
      RH=1.0D0/V/10.0D0
      RH2=RH*RH
      RH3=RH2*RH
      RH4=RH3*RH
      F1=CC(3)+2.0D0*CC(4)*RH+3.0D0*CC(5)*RH2+4.0D0*CC(6)*RH3
      F2=CC(2)+CC(3)*RH+CC(4)*RH2+CC(5)*RH3+CC(6)*RH4
      Z=RH+CC(1)*RH2-RH2*(F1/F2**2)+CC(7)*RH2*DEXP(-CC(8)*RH)+ &
      CC(9)*RH2*DEXP(-CC(10)*RH)
      P=Z*RT*10.0D0
!=====
      RETURN
      END
!----
!--------------------------------
!--------------------------------
      SUBROUTINE PISTZS(RH,ZS)
      IMPLICIT NONE
!--
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
      REAL*8 RH,ZS
!----- mit maxima
      ZS=-CC(7)*CC(8)*RH**2*DEXP(-CC(8)*RH)+2.0D0*CC(7)*RH* &
      DEXP(-CC(8)*RH)-CC(10)*CC(9)*RH**2*DEXP(-CC(10)*RH)+ &
      2.0D0*CC(9)*RH*DEXP(-CC(10)*RH)-2.0D0*RH*(4.0D0*CC(6)*RH**3+ &
      3.0D0*CC(5)*RH**2+2.0D0*CC(4)*RH+CC(3))/(CC(6)*RH**4+ &
      CC(5)*RH**3+CC(4)*RH**2+CC(3)*RH+CC(2))**2-RH**2*(12.0D0* &
      CC(6)*RH**2+6.0D0*CC(5)*RH+2.0D0*CC(4))/(CC(6)*RH**4+ &
      CC(5)*RH**3+CC(4)*RH**2+CC(3)*RH+CC(2))**2+2.0D0*RH**2* &
      (4.0D0*CC(6)*RH**3+3.0D0*CC(5)*RH**2+2.0D0*CC(4)*RH+CC(3))* &
      (4.0D0*CC(6)*RH**3+3.0D0*CC(5)*RH**2+2.0D0*CC(4)*RH+CC(3))/ &
      (CC(6)*RH**4+CC(5)*RH**3+CC(4)*RH**2+CC(3)*RH+CC(2))**3+ &
      2.0D0*CC(1)*RH+1.0D0
!=====
!----
      RETURN
      END
!----
!--------------------------------
!--------------------------------
      SUBROUTINE PISTH2ODEF(T)
      IMPLICIT NONE
!--
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
      INTEGER*4 I
      REAL*8 T,C(10,6)
!--
      C(1,1)=0.0D0
      C(1,2)=0.0D0
      C(1,3)=0.24657688D+6
      C(1,4)=0.51359951D+2
      C(1,5)=0.0D0
      C(1,6)=0.0D0
!--
      C(2,1)=0.0D0
      C(2,2)=0.0D0
      C(2,3)=0.58638965D+0
      C(2,4)=-0.28646939D-2
      C(2,5)=0.31375577D-4
      C(2,6)=0.0D0
!--
      C(3,1)=0.0D0
      C(3,2)=0.0D0
      C(3,3)=-0.62783840D+1
      C(3,4)=0.14791599D-1
      C(3,5)=0.35779579D-3
      C(3,6)=0.15432925D-7
!--
      C(4,1)=0.0D0
      C(4,2)=0.0D0
      C(4,3)=0.0D0
      C(4,4)=-0.42719875D+0
      C(4,5)=-0.16325155D-4
      C(4,6)=0.0D0
!--
      C(5,1)= 0.0D0
      C(5,2)= 0.0D0
      C(5,3)= 0.56654978D+4
      C(5,4)=-0.16580167D+2
      C(5,5)= 0.76560762D-1
      C(5,6)= 0.0D0
!--
      C(6,1)= 0.0D0
      C(6,2)= 0.0D0
      C(6,3)= 0.0D0
      C(6,4)= 0.10917883D+0
      C(6,5)= 0.0D0
      C(6,6)= 0.0D0
!--
      C(7,1)= 0.38878656D+13
      C(7,2)=-0.13494878D+9
      C(7,3)= 0.30916564D+6
      C(7,4)= 0.75591105D+1
      C(7,5)= 0.0D0
      C(7,6)= 0.0D0
!--
      C(8,1)= 0.0D0
      C(8,2)= 0.0D0
      C(8,3)=-0.65537898D+5
      C(8,4)= 0.18810675D+3
      C(8,5)= 0.0D0
      C(8,6)= 0.0D0
!--
      C(9,1)=-0.14182435D+14
      C(9,2)= 0.18165390D+9
      C(9,3)=-0.19769068D+6
      C(9,4)=-0.23530318D+2
      C(9,5)= 0.0D0
      C(9,6)= 0.0D0
!--
      C(10,1)= 0.0D0
      C(10,2)= 0.0D0
      C(10,3)= 0.92093375D+5
      C(10,4)= 0.12246777D+3
      C(10,5)= 0.0D0
      C(10,6)= 0.0D0
!====
      DO I=1,10
      CC(I)=C(I,1)/(T**4)+C(I,2)/(T**2)+C(I,3)/T+C(I,4)+ &
      C(I,5)*T+C(I,6)*T**2
      END DO
!====
      RETURN
      END
!----
!--------------------------------
!--------------------------------
      SUBROUTINE PISTCO2DEF(T)
      IMPLICIT NONE
!--
      REAL*8 CC(10),RT
      COMMON /PISTRE/ CC,RT
      INTEGER*4 I
      REAL*8 T,C(10,6)
!--
      C(1,1)= 0.0D0
      C(1,2)= 0.0D0
      C(1,3)=+0.18261340D+7
      C(1,4)=+0.79224365D+2
      C(1,5)= 0.0D0
      C(1,6)= 0.0D0
      C(2,1)= 0.0D0
      C(2,2)= 0.0D0
      C(2,3)= 0.0D0
      C(2,4)=+0.66560660D-4
      C(2,5)=+0.57152798D-5
      C(2,6)=+0.30222363D-9
!--
      C(3,1)= 0.0D0
      C(3,2)= 0.0D0
      C(3,3)= 0.0D0
      C(3,4)=+0.59957845D-2
      C(3,5)=+0.71669631D-4
      C(3,6)=+0.62416103D-8
!--
      C(4,1)= 0.0D0
      C(4,2)= 0.0D0
      C(4,3)=-0.13270279D+1
      C(4,4)=-0.15210731D+0
      C(4,5)=+0.53654244D-3
      C(4,6)=-0.71115142D-7
!--
      C(5,1)= 0.0D0
      C(5,2)= 0.0D0
      C(5,3)=+0.12456776D+0
      C(5,4)=+0.49045367D+1
      C(5,5)=+0.98220560D-2
      C(5,6)=+0.55962121D-5
!--
      C(6,1)= 0.0D0
      C(6,2)= 0.0D0
      C(6,3)= 0.0D0
      C(6,4)=+0.75522299D+0
      C(6,5)= 0.0D0
      C(6,6)= 0.0D0
!--
      C(7,1)=-0.39344644D+12
      C(7,2)=+0.90918237D+8
      C(7,3)=+0.42776716D+6
      C(7,4)=-0.22347856D+2
      C(7,5)= 0.0D0
      C(7,6)= 0.0D0
!--
      C(8,1)= 0.0D0
      C(8,2)= 0.0D0
      C(8,3)=+0.40282608D+3
      C(8,4)=+0.11971627D+3
      C(8,5)= 0.0D0
      C(8,6)= 0.0D0
!--
      C(9,1)= 0.0D0
      C(9,2)=+0.22995650D+8
      C(9,3)=-0.78971817D+5
      C(9,4)=-0.63376456D+2
      C(9,5)= 0.0D0
      C(9,6)= 0.0D0
!--
      C(10,1)= 0.0D0
      C(10,2)= 0.0D0
      C(10,3)=+0.95029765D+5
      C(10,4)=+0.18038071D+2
      C(10,5)= 0.0D0
      C(10,6)= 0.0D0
!====
      DO I=1,10
      CC(I)=C(I,1)/(T**4)+C(I,2)/(T**2)+C(I,3)/T+C(I,4)+ &
      C(I,5)*T+C(I,6)*T**2
      END DO
!====
      RETURN
      END
!---
