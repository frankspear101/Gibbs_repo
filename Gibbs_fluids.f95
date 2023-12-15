!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	subroutine Graphite_sat(Tk,Pb,Gh2ographite)
!	calculate fluid fugacities assuming graphite saturation
!	following Connolly & Cesare 1993 JMG 
!	!&! use Kerrick and Jacobs HSMRK as reference state for pure fluids
!	so this routine calculates activity = fugC&!/fugK&J
!	Gibbs uses Haar equation for water, so G(graphite) is calculated as
!	G(graphite-Haar) = G(pure-Haar) + R*T*ln(activity)
!	S and V are calculated using finite derivatives
!	H is calculated as G + TK*S

!	P in bars; T in Kelvins

	implicit real*8 (a-h,o-z)
	COMMON /FUGACITY/FH2OPUREkj,FCO2PUREkj
	data r /8.31441d0/

	xco2 = .0001	!K&J routine does not work if xco2 = 0 - this value is close to pure H2O
	MIXPIN = 1
	call WFLU3(TK,PB,XCO2,VCO2,VH2O,VMIX,ACO2,AH2O,DACO2X,&
     &      DAH2OX,DACO2T,DAH2OT,DACO2P,DAH2OP,MIXPIN,IZERO)
!	write(*,*)'K&J   ',FH2OPUREkj
        call cohfit (TK,Pb,fh2o,fco2,fo2)
	lnfh2o_cc = fh2o
	fh2o_cc = dexp(fh2o)
!	write(*,*)'!&!   ',fh2o_cc
	activity = fh2o_cc/fh2opurekj
	if(activity.gt.1.0d0)activity = 1.0d0
!	write(*,*)'activity  ',activity
	CALL WHAARx(Pb, TK, GH2O, VOLUMEstraight)
        IF (Pb .lt. 10000.) THEN          
		CALL WHAARx(Pb, TK, GH2O, VOLUMEstraight)
		else
!       	compute extra G if needed for P>10,000 
        	tenK=10000.d0
        	CALL WHAARx(tenk, Tk, GH2O, VOLUMEtenk)
        	GDIF = 0.0
        	CALL WDH78(Tk, Pb, GDIF, Vextra)
        	GH2O = GH2O + GDIF
        	END IF 

        Gh2o = Gh2o - 306674.4      !converto from G-G(1,298) to G at P and T

	gh2ographite = gh2o + R*TK*dlog(activity)
!	write(*,*)'Haar Gpure,Ggraphite, ',gh2o,gh2ographite,gh2o-gh2ographite

	return
	end
	



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine WHAAR2(T2,P2,G2,H2,S2,V2,iok)
!  -------------------------------------------------------------- 

!     Returns G,H,S and V of H2O at P and T of interest
!     units:
!     G joules/mole
!     H joules/mole
!     S joules/mole-K
!     V joules/bar

!     T2 in Kelvin
!     P2 in bars
!
!     iok is an error switch 0 = no error, 1 = error
        
      real*8 T2,P2,G2,H2,S2,V2,GH2O,VOLUME,GDIF,Vextra,Ghigh,Glow,tenK,&
     &volumestraight,volumetenk,T2plus,T2minus
      integer iok
	character*255 ppp      

!
!       Check to see that P and T are in bounds -- abort if not
!
!       note that T is in Kelvins
        iok=0
        if(P2.lt.0.5d0.or.p2.gt.40000.0.or.T2.lt.100.0.or.T2.gt.1800.0)then
!        if(P2.lt.25.or.p2.gt.30000.or.T2.lt.100.or.T2.gt.1800.)then
        	write(*,1000)T2,P2
1000	format(' T or P is out of bounds in Sub WHAAR2.  TK = ',F10.1,'  Pb = ',F10.1)
             iok=1
             return
             endif








!     calculate volume for all Pressures
!     and G for P<10000
      CALL WHAARx(P2, T2, GH2O, VOLUMEstraight)
!      write(*,*)' Volume from first call ',volume

      IF (P2 .gt. 10000.) THEN          
!     compute extra G if needed for P>10,000 
         tenK=10000.d0
         CALL WHAARx(tenk, T2, GH2O, VOLUMEtenk)
         GDIF = 0.0
         CALL WDH78(T2, P2, GDIF, Vextra)
         GH2O = GH2O + GDIF
         END IF 

      G2=Gh2o


!     calculate volume as a finite difference derivative G/P = V
!     for pressures greater than 10 kbar
!      volume_fd=VOLUMEstraight/10.0d0
!         pinc=1.0d0
!      IF (P2 .le. 10000.) THEN          
!         call WHAARx(P2+pinc,T2,Ghigh,Volume)
!         call WHAARx(P2-pinc,T2,Glow,Volume)
!         Volume_fd = (Ghigh-Glow)/(2.D0*Pinc)

!         else
!         tenK=10000.d0
!         CALL WHAARx(tenk, T2, Gtenk, VOLUMEtenk)

!         GDIF = 0.0
!         CALL WDH78(T2, P2+pinc, GDIF, Vextra)
!         Ghigh = Gtenk + GDIF

!         GDIF = 0.0
!         CALL WDH78(T2, P2-pinc, GDIF, Vextra)
!         Glow = Gtenk + GDIF
!         Volume_fd = (Ghigh-Glow)/(2.D0*Pinc)

!         END IF 

!      if(ifd.eq.1)then
!            v2=volume_fd
!            else
            v2 = volumestraight/10.0d0
!            endif

!     calculate entropy
	T2plus  = T2 + 0.1D0
	T2minus = T2 - 0.1D0
      IF (P2 .le. 10000.) THEN          
!         call WHAARx(P2,T2+.1,Ghigh,Volume)
!         call WHAARx(P2,T2-.1,Glow,Volume)
         call WHAARx(P2,T2plus,Ghigh,Volume)
         call WHAARx(P2,T2minus,Glow,Volume)

         else
         tenK=10000.d0
!         CALL WHAARx(tenk, T2+.1, Ghigh, VOLUME)
         CALL WHAARx(tenk, T2plus, Ghigh, VOLUME)
         GDIF = 0.0
!         CALL WDH78(T2+.1, P2, GDIF, Vextra)
         CALL WDH78(T2plus, P2, GDIF, Vextra)
         Ghigh = Ghigh + GDIF

!         CALL WHAARx(tenk, T2-.1, Glow, VOLUME)
         CALL WHAARx(tenk, T2minus, Glow, VOLUME)
         GDIF = 0.0
!         CALL WDH78(T2-.1, P2, GDIF, Vextra)
         CALL WDH78(T2minus, P2, GDIF, Vextra)
         Glow = Glow + GDIF
         END IF 
        
      S2 = -(Ghigh-Glow)/0.2D0

      G2 = G2 - 306674.4      !converto from G-G(1,298) to G at P and T

      H2=G2 + S2*T2

      return
      end








!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


      SUBROUTINE WHAARx(P, T, GH2O, VH2O)
! 
!     FSS Received this code from Rob Berman, August, 1993

!     Calculates and returns the Gibbs free energy of Water and 
!            its volume using the Haar equation of State.  Written 
!            by Christian De Capitan, Dept. Geological Sciences, 
!            UBC, Vancouver, BC, Canada 
! 
!     This routine uses a Redlich-Kwong to obtain a first guess 
!            for the density of water, thus speeding things up 
! 
      IMPLICIT DOUBLEPRECISION(A - H,O - Z)
      DIMENSION TAUI(0:6), ERMI(0:9), GI(40), KI(40), LI(40)
      DIMENSION CI(18), RHOI(37:40), TTTI(37:40), ALPI(37:40)
      DIMENSION BETI(37:40)

! 
! -----GI ARE IN (bar cc / g)  =  10 * (J / g) 
! 
      DATA GI /-.53062968529023D4, .22744901424408D5, .78779333020687D4, &
      -.69830527374994D3, .17863832875422D6, -.39514731563338D6, & 
      .33803884280753D6, -.13855050202703D6, -.25637436613260D7, &
          .48212575981415D7, -.34183016969660D7, .12223156417448D7,& 
          .11797433655832D8, -.21734810110373D8, .10829952168620D8, &
          -.25441998064049D7, -.31377774947767D8, .52911910757704D8, &
          -.13802577177877D8, -.25109914369001D7, .46561826115608D8, &
          -.72752773275387D8, .41774246148294D7, .14016358244614D8, &
          -.31555231392127D8, .47929666384584D8, .40912664781209D7, &
          -.13626369388386D8, .69625220862664D7, -.10834900096447D8, &
          -.22722827401688D7, .38365486000660D7, .68833257944332D5, &
          .21757245522644D6, -.26627944829770D5, -.70730418082074D6, &
          -.225D1, -1.68D1, .055D1, -93.0D1/ 
      DATA KI /4*1, 4*2, 4*3, 4*4, 4*5, 4*6, 4*7, 4*9, 2*3, 1, 5, 3*2,  4/ 

      DATA LI /1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, & 
          6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 0, 3*3, 0, 2, 0, 0/ 
      DATA CI /.19730271018D2, .209662681977D2, -.483429455355D0, &
         .605743189245D1, 22.56023885D0, -9.87532442D0, &
          -.43135538513D1, .458155781D0, -.47754901883D-1, &
          .41238460633D-2, -.27929052852D-3, .14481695261D-4, & 
          -.56473658748D-6, .16200446D-7, -.3303822796D-9, &
          .451916067368D-11, -.370734122708D-13, .137546068238D-15/ 
      DATA RHOI /0.319D0, 0.310D0, 0.310D0, 1.55D0/ 
      DATA TTTI /640.0D0, 640.0D0, 641.6D0, 270.0D0/ 
      DATA ALPI /34.0D0, 40.0D0, 30.0D0, 1050.0D0/ 
      DATA BETI /2.0D4, 2.0D4, 4.0D4, 25.0D0/ 
! 
      R = 4.6152D0
      RT = R * T
      NLOW = 40
      NHIGH = 20
      IF (T .LT. 449.35D0) NHIGH = 40
! 
! ----GREF CALCULATED WITH THIS ROUTINE AT 25 ! AND 1 BAR 
! 
      GREF = -54955.2356146121147D0
      T0 = 647.073D0
! 
! -----The values (T/T0)**i are stored in the array TAUI(i) 
! 
      TAUI(0) = 1.D0
      TAUI(1) = T / T0
      DO 10 I = 2, 6
        TAUI(I) = TAUI(I - 1) * TAUI(1)
   10 CONTINUE 
! 
      B = -0.3540782D0 * DLOG(TAUI(1)) + 0.7478629D0+0.007159876D0 /TAUI(3) - 0.003528426D0 / TAUI(5)
      BB = 1.1278334D0-0.5944001D0 / TAUI(1) - 5.010996D0 / TAUI(2) + 0.63684256D0 / TAUI(4)
! 


      PS = 220.55D0
      IF (T .LE. 647.25) CALL PSAT2_B(T, PS)
! 
! -----SET INITIAL GUESS FOR RHO USING THB-FIT TO REDLICH-KWONG 
! 
      ARK = 1.279186D8-2.241415D4 * T
      BRK = 1.428062D1+6.092237D-4 * T
      RR = 8.31441D0
      OFT = ARK / (P*DSQRT(T))
      BUK = -10.0D0 * RR * T / P
      CUK = OFT - BRK * BRK + BRK * BUK
      DUK = -BRK * OFT
      CALL KUBIK(BUK, CUK, DUK, X1, X2, X2I, X3)
      IF (X2I .NE. 0.0D0) THEN 
        VOL = X1
      ELSE 
        IF (P .LT. PS) THEN 
          VOL = DMAX1(X1,X2,X3)
        ELSE 
          VOL = DMIN1(X1,X2,X3)
        END IF 
      END IF 
      IF (VOL .LE. 0.0D0) THEN 
        RHN = 1.9D0
      ELSE 
        RHN = 1.0D0 / VOL * 18.0152D0
      END IF 
! 
! -----FIND THE TRUE(?) RH(T,P) 
! -----NOTE: PR = PRESSURE CORRESPONDING TO GUESSED RH 
!            DPR = (dP / dRH) 
!            the values (1-EXP(-RH))**i are stored in the array ERMI(i) 
! 
      DO 50 LOO = 1, 100
        RH = RHN
        IF (RH .LE. 0.0D0) RH = 1.D-8
        IF (RH .GT. 1.9D0) RH = 1.9D0
        RH2 = RH ** 2
        Y = RH * B / 4.D0
        ER = DEXP(-RH)
        Y3 = (1.0D0-Y) ** 3
        ALY = 11.D0 * Y
        BETY = 44.33333333333333D0 * Y * Y
        F1 = (1.D0+ALY + BETY) / Y3
        F2 = 4.D0 * Y * (BB/B - 3.5D0)
        ERMI(0) = 1.0D0
        ERMI(1) = 1.0D0-ER
        DO 20 I = 2, 9
          ERMI(I) = ERMI(I - 1) * ERMI(1)
   20   CONTINUE 
        PR = 0.0D0
        DPR = 0.0D0
        DO 30 I = 1, 36
          S = GI(I) / TAUI(LI(I)) * ERMI(KI(I) - 1)
          PR = PR + S
          DPR = DPR + (2D0+RH*(KI(I)*ER - 1D0)/ERMI(1)) * S
   30   CONTINUE 
        DO 40 I = NLOW, NHIGH
          DEL = RH / RHOI(I) - 1.0D0
          RHOI2 = RHOI(I) * RHOI(I)
          TAU = T / TTTI(I) - 1.0D0
          ABC = -ALPI(I) * DEL ** KI(I) - BETI(I) * TAU ** 2
          IF (ABC .GT. - 100.0D00) THEN 
            Q10 = GI(I) * DEL ** LI(I) * DEXP(ABC)
          ELSE 
            Q10 = 0.0D00
          END IF 
          QM = LI(I) / DEL - KI(I) * ALPI(I) * DEL ** (KI(I) - 1)
          S = Q10 * QM * RH2 / RHOI(I)
          PR = PR + S
          DPR = DPR + S * (2.0D0/RH + QM/RHOI(I)) - RH2 / RHOI2 * Q10 * (LI(I)/DEL/DEL + KI(I)*(KI(I) - 1)*ALPI(I)*DEL**(KI(I) - 2))
   40   CONTINUE 
        PR = RH * (RH*ER*PR + RT*(F1 + F2))
        DPR = RH * ER * DPR + RT * ((1.0D0+2D0*ALY + 3D0*BETY)/Y3 + 3D0*Y*F1/(1.0D0-Y) + 2D0*F2)
! 
! ----- 
! 


        IF (DPR .LE. 0.0D0) THEN 
          IF (P .LE. PS) THEN 
            RHN = RHN * 0.95D0
          ELSE 
            RHN = RHN * 1.05D0
          END IF 
        ELSE 
          IF (DPR .LT. 0.01D0) DPR = 0.01D0
          X = (P - PR) / DPR
          IF (DABS(X) .GT. 0.1D0) X = 0.1D0 * X / DABS(X)
          RHN = RH + X
        END IF 
        DP = DABS(1.0D0-PR/P)
        DR = DABS(1.0D0-RHN/RH)
        IF (DP .LT. 5.D-2 .AND. DR .LT. 5.D-2) GO TO 60
   50 CONTINUE 
   60 RH = RHN
! 
! ----- 
! 
      Y = RH * B / 4.D0
      X = 1.0D0-Y
      ER = DEXP(-RH)
      ERMI(0) = 1.0D0
      ERMI(1) = 1.0D0-ER
      DO 70 I = 2, 9
        ERMI(I) = ERMI(I - 1) * ERMI(1)
   70 CONTINUE 
! 
! -----CALCULATE BASE FUNCTION 
! 
      AA = RT * (-DLOG(X) - 43.33333333333333D0/X + 28.16666666666667D0/X/X &
           + 4D0*Y*(BB/B - 3.5D0) + 15.16666666666667D0+DLOG(RH*RT/1.01325D0))
! 
! -----CALCULATE RESIDUAL FUNCTION 
! 
      DO 80 I = 1, 36
        AA = AA + GI(I) / KI(I) / TAUI(LI(I)) * ERMI(KI(I))
   80 CONTINUE 
      DO 90 I = NLOW, NHIGH
        DEL = RH / RHOI(I) - 1.0D0
        TAU = T / TTTI(I) - 1.0D0
        ABC = -ALPI(I) * DEL ** KI(I) - BETI(I) * TAU ** 2
        IF (ABC .GT. - 100.0D00) THEN 
          AA = AA + GI(I) * DEL ** LI(I) * DEXP(ABC)
        ELSE 
          AA = AA
        END IF 
   90 CONTINUE 
! 
! -----CALCULATE IDEAL GAS FUNCTION 
! 
      TR = T / 1.0D2
      W = TR ** (-3)
      AID = 1.D0+(CI(1)/TR + CI(2)) * DLOG(TR)
      DO 100 I = 3, 18
        AID = AID + CI(I) * W
        W = W * TR
  100 CONTINUE 
      AA = AA - RT * AID

! 
! -----CALCULATE G = AA + P/RH  AND  V = 1/RH 
! 
      GH2O = ((AA + P/RH)*1.80152D0-GREF)
      VH2O = (1.0D0/RH) * 18.0152D0
!     WRITE (91,911) P,T,GH2O,VH2O 
!a11  FORMAT (4F15.3) 

	
      RETURN 
      END 

      SUBROUTINE WDH78(TK, P, GDIF, VDH)
      IMPLICIT DOUBLEPRECISION(A - H,O - Z)
! 
!      RETURNS DIFFERENCE IN FREE ENERGY OF WATER 
!      BETWEEN P,T AND 10KB,T 
! 
      DOUBLE PRECISION A(5,5)
! 
      A(1,1) = -5.6130073D+04
      A(1,2) = 3.8101798D-01
      A(1,3) = -2.1167697D-06
      A(1,4) = 2.0266445D-11
      A(1,5) = -8.3225572D-17
      A(2,1) = -1.5285559D+01
      A(2,2) = 1.3752390D-04
      A(2,3) = -1.5586868D-09
      A(2,4) = 6.6329577D-15
      A(3,1) = -2.6092451D-02
      A(3,2) = 3.5988857D-08
      A(3,3) = -2.7916588D-14
      A(4,1) = 1.7140501D-05
      A(4,2) = -1.6860893D-11
      A(5,1) = -6.0126987D-09
! 
      T = TK - 273.15
      II = 0
      P1 = 10000.0D0
! 
   10 GH2O = 0.0D0
      VH2O = 0.0D0
      DO 30 J = 1, 5
        DO 20 L = 1, (6 - J)
          JJ = J - 1
          LL = L - 1
          IF (JJ .EQ. 0 .AND. LL .EQ. 0) GH2O = GH2O + A(J,L)
          IF (JJ .GT. 0 .AND. LL .GT. 0) GH2O = GH2O + (A(J,L)*T**(JJ)*P1**(LL))

          IF (JJ .EQ. 0 .AND. LL .GT. 0) GH2O = GH2O + (A(J,L)*P1**(LL))
          IF (JJ .GT. 0 .AND. LL .EQ. 0) GH2O = GH2O + (A(J,L)*T**(JJ))
!              IF (LL.NE.0) VH2O = VH2O + (A(J,L) * LL * T**JJ 
!    .                                                     * P**(LL-1)) 
   20   CONTINUE 
   30 CONTINUE 
! 
      VDH = VH2O / .0239
! 
      IF (II .EQ. 0) THEN 
        G10KB = GH2O
        P1 = P
        II = 1
        GO TO 10
      END IF 
! 
   40 GDIF = (GH2O - G10KB) * 4.184
! 
      RETURN 
      END 
 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  -------------------------------------------------------------- 
      SUBROUTINE PSAT2_B(T, PS)
! 
!      Written by Christian De Capitan, Dept. Geological Sciences, 
!            UBC, Vancouver, BC, Canada 
! 
      DOUBLE PRECISION T, PS, A(8), W, WSQ, V, FF
! 
      DATA A /-7.8889166D0, 2.5514255D0, -6.716169D0, 33.239495D0,  & 
             -105.38479D0, 174.35319D0, -148.39348D0, 48.631602D0/ 
! 
      IF (T .LE. 314.00D0) THEN 
        PS = DEXP(6.3573118D0-8858.843D0/T + 607.56335D0/(T**0.6D0))
      ELSE 
        V = T / 647.25D0
        W = DABS(1.0D0-V)
        WSQ = DSQRT(W)
        FF = 0.0D0
        DO 10 I = 1, 8
          FF = FF + A(I) * W
          W = W * WSQ
   10   CONTINUE 
        PS = 220.93D0 * DEXP(FF/V)
      END IF 
      RETURN 
      END 
 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  -------------------------------------------------------------- 
      SUBROUTINE KUBIK(B, C, D, X1, X2, X2I, X3)
! 
!      Written by Christian De Capitan, Dept. Geological Sciences, 
!            UBC, Vancouver, BC, Canada 
! 
      DOUBLE PRECISION B, C, D, Q, P, R, PI, PHI3, FF
      DOUBLE PRECISION X1, X2, X2I, X3
! 
      PI = 3.14159263538979D0
      X2 = 0.0D0
      X2I = 0.0D0
      X3 = 0.0D0
      IF (C .EQ. 0.0D0 .AND. D .EQ. 0.0D0) THEN 
        X1 = -B
        RETURN 
      END IF 
      Q = ((2.D0*B**3)/(27.D0) - (B*C)/(3.D0) + D) / 2.D0
      P = (3.D0*C - B**2) / (9.D0)
      FF = DABS(P)
      R = DSQRT(FF)
      FF = R * Q
      IF (FF .LT. 0.0D0) R = -R
      FF = Q / (R**3)
! 
      IF (P .GT. 0.0D0) THEN 
        PHI3 = DLOG(FF + DSQRT(FF**2 + 1.D0)) / 3.D0
        X1 = -R * (DEXP(PHI3) - DEXP(-PHI3)) - B / 3.D0
        X2I = 1.0D0
      ELSE 
        IF (Q**2 + P**3 .GT. 0.0D0) THEN 
          PHI3 = DLOG(FF + DSQRT(FF**2 - 1.D0)) / 3.D0
          X1 = -R * (DEXP(PHI3) + DEXP(-PHI3)) - B / 3.D0
          X2I = 1.0D0
        ELSE 
          PHI3 = DATAN(DSQRT(1.D0-FF**2)/FF) / 3.D0
          X1 = -2.D0 * R * DCOS(PHI3) - B / 3.D0
          X2 = 2.D0 * R * DCOS(PI/3.D0-PHI3) - B / 3.D0
          X2I = 0.0D0
          X3 = 2.D0 * R * DCOS(PI/3.D0+PHI3) - B / 3.D0
        END IF 
      END IF 
      RETURN 
      END 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 
	SUBROUTINE WFLU3(TK,PB,XCO2,VCO2,VH2O,VMIX,ACO2,AH2O,DACO2X, &
     &      DAH2OX,DACO2T,DAH2OT,DACO2P,DAH2OP,MIXPIN,IZERO)
!	  SUBROUTINE to handle mixed H2O-CO2 fluids and make calls
!	     to Kerrick and Jacob's (1981) subroutine 
!	  	
!	  LIST OF VARIABLES
!	TK,
!	PB,
!	XCO2,
!	VCO2,   Molar volume of pure CO2 at P,T
!	VH2O,   Molar volume of pure H2O at P,T
!	ACO2	
!	AH2O,
!	VMIX,
!	DAH2OX, DERIVATIVE OF ln a W R T COMPOSITION
!       DACO2X,
! 	DAH2OT, DERIVATIVE OF ln a W R T TEMPERATURE
! 	DACO2T,
! 	DAH2OP, DERIVATIVE OF ln a W R T PRESSURE
! 	DACO2P,

!     DAH2OT = log(aH2O)/T
!     DACO2T = log(aCO2)/T
!     DAH2OP = log(aH2O)/P
!     DACO2P = log(aCO2)/P
!     DAH2OX = log(aH2O)/XCO2
!     DACO2X = log(aCO2)/XCO2
!	
!               NOTE!!!! These derivatives are changed slightly from those used
!       in the Gibbs'90 program, which multiplies the derivatives by R*T in this
!       subroutine.  Here, we do the R*T multiplication in the subroutine ASV
!



! 	MIXPIN, FLAG
! 	IZERO,  FAILURE FLAG IN sub COMPUTE IN GIBBS
!	
!	TEMP, TEMPORARY VARIABLE SOMETIMES IS VOLUME IN MIXED FLUID
!	
!
!
!
!	
	
!	THIS SUBROUTINE IS IN DOUBLE PRECISION TO MATCH PROGRAM GIBBS	
 
 
!	LOCAL VARIABLES

	IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 VH2O,TEMP,TC
	REAL*8 TK,PB,VCO2,XCO2,XH2O,ACO2,AH2O,VMIX,DAH2OX,DACO2X
	REAL*8 DAH2OT,DACO2T,DAH2OP,DACO2P
!
	DATA R/8.31441D0/
!
	XH2O=1.D0 - XCO2
	TC=TK - 273.15D0
	
!--------------
!	THIS SECTION IS FOR CALCULATING WATER AND MIXED FLUID VOLUMES
!	   ACTIVITIES AND PARTIAL DERIVATIVES OF ACTIVITY
3100	CONTINUE
	IF (XCO2.LT.0.OR.XH2O.LT.0) THEN
	   WRITE(*,*)' A FLUID COMPONENT IS NEGATIVE IN sub FLU'
	   WRITE(*,*)'         END CALCULATION'
	   IZERO=1
	   RETURN
	ENDIF
	
	     
	  
	IF (XCO2.NE.0) THEN	
	
!	   HERE FOR MIXED FLUIDS
	   MIXPIN=3
!	     temp is here mixed volume
	   call MIXFLU(TC,PB,XCO2,TEMP,VCO2,VH2O,ACO2,AH2O,MIXPIN)
	
	
!
	ENDIF

!	entropy in joules/mole-deg
!	VOLUME WAS RETURNED IN joule/BAR
!

!	ASSIGNMENT OF PASSED VARIABLES AND COMPUTATION OF D (ACTIVITY)
!	 CASE WHEN PURE ENDS OR NEARLY PURE ENDS

	IF(XCO2.LT.0.00005D0.or.XCO2.GT.0.99995D0)THEN
!	 NEARLY PURE WATER	
        
          
	   write(*,*)'*******************************************'
	   write(*,*)'Error in subroutine WFLU'
	   write(*,*)'Calculations fail when fluid is nearly pure:'
	   write(*,*)'XCO2 =',xco2
	   write(*,*)'*******************************************'
	   pause 'Hit return to continue'
	   izero=1
	   RETURN
	ENDIF
!
!       Composition derivatives
!       IF (XCO2.GE.0.010D0. AND .XCO2.LE.0.990D0)THEN
       IF (XCO2.GE.0.0050D0.AND.XCO2.LE.0.9950D0)THEN
!	  THE INTERMEDIATE RANGE OF COMPOSITION 
!	  TEMP=XCO2+0.005D0
	  TEMP=XCO2+0.001D0
	  call MIXFLU(TC,PB,TEMP,VMIX,VCO2,VH2O,ACO2A,AH2OA,MIXPIN)
!	  TEMP=XCO2-0.005D0
	  TEMP=XCO2-0.001D0
	  call MIXFLU(TC,PB,TEMP,VMIX,VCO2,VH2O,ACO2B,AH2OB,MIXPIN)
!	
!	  SLOPES  
!	  DACO2X= (dlog(ACO2A)-dlog(ACO2B))/0.01D0
!	  DAH2OX= (dlog(AH2OB)-dlog(AH2OA))/0.01D0
	  DACO2X= (dlog(ACO2A)-dlog(ACO2B))/0.002D0
	  DAH2OX= (dlog(AH2OA)-dlog(AH2OB))/0.002D0
	  
!	  DACO2X=R*TK*DACO2X
!	  DAH2OX=R*TK*DAH2OX
	
       ELSE
	  IF(XCO2.LT.0.005D0)THEN
!	    TEMP=XCO2+0.005D0
	    TEMP=XCO2+0.001D0
	    call MIXFLU(TC,PB,TEMP,VMIX,VCO2,VH2O,ACO2A,AH2OA,MIXPIN)

!	    DACO2X=(dlog(ACO2A)-dlog(ACO2))/0.005D0
!	    DAH2OX=(dlog(AH2O)-dlog(AH2OA))/0.005D0
	    DACO2X=(dlog(ACO2A)-dlog(ACO2))/0.001D0
	    DAH2OX=(dlog(AH2OA)-dlog(AH2O))/0.001D0
	    
!	    DACO2X=R*TK*DACO2X
!	    DAH2OX=R*TK*DAH2OX
	  
	  ELSE
!	    TEMP=XCO2-0.005D0
	    TEMP=XCO2-0.001D0
	    call MIXFLU(TC,PB,TEMP,VMIX,VCO2,VH2O,ACO2B,AH2OB,MIXPIN)

!	    DACO2X=(dlog(ACO2)-dlog(ACO2B))/0.005D0
!	    DAH2OX=(dlog(AH2OB)-dlog(AH2O))/0.005D0
	    DACO2X=(dlog(ACO2)-dlog(ACO2B))/0.001D0
	    DAH2OX=(dlog(AH2O)-dlog(AH2OB))/0.001D0
	    
!	    DACO2X=R*TK*DACO2X
!	    DAH2OX=R*TK*DAH2OX
	  ENDIF
       ENDIF
       
!
!	calculate pressure derivatives       
        TEMP=PB+10.D0
	call MIXFLU(TC,TEMP,XCO2,VMIX,VCO2,VH2O,ACO2A,AH2OA,MIXPIN)

	DACO2P=(dlog(ACO2A)-dlog(ACO2))/10.D0
	DAH2OP=(dlog(AH2OA)-dlog(AH2O))/10.D0
	 
!	calculate temperature derivatives
        TEMP=TC+1.D0
	call MIXFLU(TEMP,PB,XCO2,VMIX,VCO2,VH2O,ACO2A,AH2OA,MIXPIN)

	DACO2T=(dlog(ACO2A)-dlog(ACO2))/1.D0
	DAH2OT=(dlog(AH2OA)-dlog(AH2O))/1.D0
	
! 	DAH2OT=R*TK*DAH2OT
!	DAH2OP=R*TK*DAH2OP
! 	DACO2T=R*TK*DACO2T
! 	DACO2P=R*TK*DACO2P
	
	END
!	
	
	
	
	
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!    FILE FLU.FOR
!
!	 SUB MIXFLU
!	     SUB ZPURE
!	     SUB ZMIX
!	     SUB NEWRAP_B
!	     SUB FPURE
!	     SUB FMIX
!
!
       SUBROUTINE MIXFLU (TC,PBAR,XC,VM,VCO2,VH2O,ACO2,AH2O,MIXPIN)

!
!	MODIFIED FROM:
!	KERRICK, D. M. AND JACOBS, G. K., 1981.  A MODIFIED
!   REDLICH-KWONG EQUATION FOR H2O, CO2, AND H20-CO2 MIXTURES AT
!   ELEVATED PRESSURES AND TEMPURATURES, AM. JOUR. SCI., 1981.
!		      BY T. MENARD, APRIL, 1988
!
!	Subroutines to compute fugacity of H2O-CO2 mixtures from program
!	of Kerrick and Jacobs (Computers in Geosciences, 1981)
	IMPLICIT REAL*8(A-H,P,R-Z)
	COMMON /MVAR/BC,CC,DC,EC,BW,CW,DW,EW,J,L,M
	COMMON /MNRVAR/PK,P,T,R,T12
	COMMON /FUGACITY/FH2OPUREkj,FCO2PUREkj
!
	
!
!	UNITS ARE CALCULATED IN JOULES AND CUBIC CM
!	   VOLUME IS FINALLY CONVERTED TO JOULES/BAR AND RETURNED
!
!	VARIABLES AND RULES ARE AS FOLLOWS:
!
!	 ACO2,AH2O....ACTIVITY OF CO2 AND H2O, RESPECTIVELY
!	 BC,BW,BM.....COVOLUME OF CO2,H2O AND MIXTURE; CC/MOLE
!	 CC,CW,CM.....ATTRACTIVE TERM FOR CO2, H2O, AND MIXTURE IN
!		      MRK EQUATION ; BAR*(CC**2)**2DSQRT(T)/MOLE**2
!	 DC, DW,,DM...ATTRACTIVE TERM FOR CO2, H2O, AND MIXTURE IN
!		      MRK EQUATION; BAR*(CC**3)*DSQRT(T)/MOLE**3
!	 EC,EW,EM.....ATTRACTIVE TERM FOR CO2, H2O, AND MIXTURE IN
!		      MRK EQUATION; BAR*(CC**4)*DSQRT(T)/MOLE**4
!	 CIJ,DIJ,EIJ..CROSS COEFFICIENTS OF !,D,E
!	 FKCM,FKWM....FUGACITY COEFFICIENT OF CO2 AND H2O IN
!		      THE FLUID MIXTURE
!	 FKCP,FKWP....FUGACITY COEFFICIENTS OF PURE CO2 AND PURE
!		      H2O, RESPECTIVELY
!	 PK,PBAR.........PRESSURE; KBARS,BARS, RESPECTIVELY
!	 R............GAS CONSTANT; 83.14 CC*BARS/MOLE*K
!	 TC,T.........TEMPERATURE; CELSIUS,KELVIN, RESPECTIVELY
!	 VC,VW,VM.....MOLAR VOLUME OF CO2, H2O, AND MIXTURE; CC/MOLE
!	 XC,XCO2......MOLE FRACTION OF CO2 IN THE FLUID MIXTURE
!	 XW,XH2O......MOLE FRACTION OF H2O IN THE FLUID MIXTURE
!	 Y............B/V4; VARIABLE IN HARD SPHERE-EQUATION
!	 ZC,ZW,ZM.....COMPRESSIBILITY OF CO2, H2O, AND MIXTURE
!	     XC AND XW	    ARE INTERNAL TO THIS SUBROUTINE
!
!	 MIXPIN IS A FLAG TO TELL WHERE THIS SUBROUTINE WAS CALLED FROM
!	    MIXPIN = 1	FUGACITY CALC	  NOT USED
!		     2	ACTIVITY	  NOT USED
!		     3	S		  sub GSV  (also gets mixed volume)
!		     4	V MIXED 	  sub VOLUME (IN sub COMPUTE)
!		     5	WATER VOLUME	  sub GSV (WATER SECTION)
!
!
!	 DEFINITION OF CONSTANRS AND INPUT FOR THE CALCULATION OF
!    THE P, T, XCO2 VALUES WHICH ARE USED THROUGHOUT THE PROGRAM.
!
!
	IF(TC.GT.1050.0D0)THEN
	   WRITE(*, 29)
29	 FORMAT(1X,'*****WARNING***1050 DEG ! IS MAX TEMP FOR MIXTURES')
	   WRITE(*,*)' TC=',TC
	ENDIF
	R=83.14D0
	BW=29.00D0
	BC=58.00D0
!	CALCULATION OF PARAMETERS USED IN PROGRAM
	XW=1.0D0-XC
	P=PBAR
	PK=P/1000.0D0
	T=TC+273.15D0
	T15=DSQRT(T**3)
	T12=DSQRT(T)
	RT=R*T15
	CC=(28.31D0+0.10721D0*T-0.00000881D0*T*T)*1000000.0D0
	DC=(9380.0D0-8.53D0*T+0.001189D0*T*T)*1000000.0D0
	EC=(-368654.0D0+715.9D0*T+0.1534D0*T*T)*1000000.0D0
	CW=(290.78D0-0.30276D0*T+0.00014774D0*T*T)*1000000.0D0
	DW=(-8374.0D0+19.437D0*T-0.008148D0*T*T)*1000000.0D0
	EW=(76600.0D0-133.9D0*T+0.1071D0*T*T)*1000000.0D0
	IF(TC.GT.1050.0D0)GO TO 299
	BM=(BC*XC)+(BW*XW)
	CIJ=DSQRT(CC*CW)
	DIJ=DSQRT(DC*DW)
	EIJ=DSQRT(EC*EW)
	CM=(CC*XC*XC)+(CW*XW*XW)+(2.0D0*XC*XW*CIJ)
	DM=(DC*XC*XC)+(DW*XW*XW)+(2.0D0*XC*XW*DIJ)
	EM=(EC*XC*XC)+(EW*XW*XW)+(2.0D0*XC*XW*EIJ)
!
!	  SUBROUTINES ZPURE, ZMIX, FPURE, AND FMIX ARE CALLED TO
!    CALCULATE:
!
!	 1) Z OF PURE CO2 AND H2O RESPECTIVELY
!	 2) Z OF CO2-H2O MIXTURES
!	 3) FUGACITY COEFFICIENTS OF CO2 AND H2O RESPECTIVELY
!	 4) FUGACITY COEFFICIENTS OF CO2 AN H2O IN THE MIXTURE
!	 DELTA S MIX IS THEN CALCULATED
!
!    AT EACH P, T, XCO2 CONDITION.
!
!    NOTE THAT VOLUME IS USED IN FUGACITY CALCULAITON
!	       FUGACITY IS USED IN ACTIVITY CALCULATION
!	       ACTIVITY IS USED TO CALCULATE DELTA S MIX
!
       IF (MIXPIN.EQ.5)THEN
	    J=2
	 CALL ZPURE(PK,ZC,VC,ZW,VW,TC)
	 VW=VM
	 GO TO 810
      ENDIF
299	DO 300 J=1,2
	   CALL ZPURE(PK,ZC,VC,ZW,VW,TC)
300	CONTINUE
	VH2O=VW
	VCO2=VC
	IF(TC.GT.1050.0D0)GO TO 499
	CALL ZMIX(XC,XW,VM,ZM,VC,VW,BM,CM,DM,EM)
	IF (MIXPIN.EQ.4)GO TO 810
!	   FOR MIXED VOULUME CALCULATION
499	DO 500 L=1,2
	   CALL FPURE(RT,FKWP,FKCP,FCP,ZC,VC,ZW,VW)
500	CONTINUE
	IF(TC.GT.1050.0D0)GO TO 810
	DO 700 M=1,2
	   CALL FMIX(RT,CIJ,DIJ,EIJ,XC,XW,BM,FCCM,FCWM,FCM,VM,ZM,CM,DM,EM)

	   ACO2  = FCCM * XC / FKCP
	   AH2O = FCWM * XW / FKWP
700	CONTINUE
810	CONTINUE
!
	
!	CONVERT from cm**3 TO joule/BAR
	VM=VM/10.d0
	vco2=vco2/10.d0
	vh2o=vh2o/10.d0
!
	FH2OPUREkj = FKWP*P
	FCO2PUREkj = FKCP*P
!	write(*,*)'Fugacity coeff-pure  ',FKWP
	RETURN
	END
!
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
	SUBROUTINE ZPURE(PK,ZC,VC,ZW,VW,TC)
	IMPLICIT REAL*8(A-H,P,R-Z)
	COMMON /MVAR/BC,CC,DC,EC,BW,CW,DW,EW,J,L,M
!	SOLVE FOR VOLUME OF CO2 AND H2O.AN INITIAL GUESS (VI) IS MADE
!	AND THE SUBROUTINE NEWRAP SOLVES FOR THE EXACT VOLUME
!	BY A NEWTON-RAPHSON TECHNIQUE
	IF(J.NE.1)GO TO 210
	B=BC
	!=CC
	D=DC
	E=EC
	IF(PK.GE.1.0D0)VI=35.00
	IF(PK.GE.0.10D0.AND.PK.LT.1.0D0)VI=100.0
	IF(PK.GE.0.005D0.AND.PK.LT.0.10D0)VI=500.0
	IF(PK.LT.0.005D0)VI=50000.0
	GO TO 250
210	B=BW
	!=CW
	D=DW
	E=EW
	IF(PK.GE.1.0D0)VI=15.0
	IF(PK.GE.0.6D0.AND.PK.LT.1.0D0)VI=22.50
	IF(PK.GE.0.21D0.AND.PK.LT.0.60D0.AND.TC.GE.550.0D0)VI=75.0
	IF(PK.GE.0.21D0.AND.PK.LT.0.60D0.AND.TC.LT.550.0D0)VI=35.0
	IF(PK.GE.0.1D0.AND.PK.LT.0.21D0.AND.TC.LT.400.0D0)VI=15.0
	IF(PK.GE.0.1D0.AND.PK.LT.0.21D0.AND.TC.GE.400.0D0)VI=100.0
	IF(PK.GE.0.005D0.AND.PK.LT.0.10D0)VI=500.0
	IF(PK.LT.0.005D0)VI=7000.00
250	CALL NEWRAP_B(B,C,D,E,Z,V,VI)
	IF(J.NE.1)GO TO 260
	ZC=Z
	VC=V
	GO TO 299
260	ZW=Z
	VW=V
299	RETURN
	END
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
	SUBROUTINE ZMIX(XC,XW,VM,ZM,VC,VW,BM,CM,DM,EM)
	IMPLICIT REAL*8(A-H,P,R-Z)
	COMMON /MVAR/BC,CC,DC,EC,BW,CW,DW,EW,J,L,M
!	SUBROUTINE TO CALCULATE THE VOLUME OF A MIXTURE OF CO2 AND H2O
!	AT T AND P AND XCO2. THE MOLAR VOLUME OF CO2 AND H2O CALCULATED
!	IN ZPURE ARE USED TO DEFINE AN INITIAL ESTIMATE. THEN NEWRAP
!	IS USED TO CALCULATE THE VOL OF THE MIXTURE
	VI=(VC*XC)+(VW*XW)
	B=BM
	!=CM
	D=DM
	E=EM
	CALL NEWRAP_B(B,C,D,E,Z,V,VI)
	VM=V
	ZM=Z
	RETURN
	END
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
	SUBROUTINE NEWRAP_B(B,C,D,E,Z,V,VI)
	IMPLICIT REAL*8(A-H,P,R-Z)
	COMMON /MNRVAR/PK,P,T,R,T12
!	NEWTON-RAPHSON ROUTINE
!
!	INITIALIZE PARAMETERS
	DO 350 K=1,50
	   Y=B/(4.0D0*VI)
	   X=(1.0D0-Y)
	   BI=VI+B
	   BI2=(VI+B)**2
!	   DEFINITION OF F(X)
	   PN=1.0D0+Y+(Y**2)-Y**3
	   PR=(PN/(VI*(X**3)))*R*T
	   PA1=C+(D/VI)+(E/(VI*VI))
	   PA2=PA1/(T12*VI*BI)
	   F=PR-PA2-P
!	   DEFINITION OF DIFFERENTIAL OF F(X)
	   D1=(-3.0D0*B)/(4.0D0*(VI**3)*X**4)
	   D2=-1.0D0/((VI**2)*(X**3))
	   D3=1.0D0/(VI*(X**3))
	   D4=-B/(4.0D0*VI**2)
	   D5=-2.0D0*(B**2)/(16.0D0*(VI**3))
	   D6=3.0D0*(B**3)/(64.0D0*(VI**4))
	   DPR=((PN*(D1+D2))+(D3*(D4+D5+D6)))*R*T
	   D7=(-1.0D0/(VI*BI2))+(-1.0D0/(VI**2*BI))
	   D8=1.0D0/(VI*BI)
	   D9=-D/VI**2+(-2.0D0*E)/(VI**3)
	   DPA=(PA1*D7+D8*D9)/T12
	   DF=DPR-DPA
!	   CALCULATION OF V(K+1) AND CHECK FOR CONVERGENGE
	   V=VI-(F/DF)
	   DIFF=DABS(V-VI)
	   IF(DIFF.LT.0.01D0)GO TO 360
	   VI=V
350	CONTINUE
	WRITE(*,*)' NO CONVERGENCE AFTER 50 ITERATIONS'
360	Z=(V*P)/(R*T)
	RETURN
	END
!
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
	SUBROUTINE FPURE(RT,FKWP,FKCP,FCP,ZC,VC,ZW,VW)
	IMPLICIT REAL*8(A-H,P,R-Z)
	COMMON /MVAR/BC,CC,DC,EC,BW,CW,DW,EW,J,L,M
!	CALCULATE FUGACITY COEFFS OF PURE CO2 AND H2O AT T AND P
	IF(L.NE.1)GO TO 450
	   B=BC
	   !=CC
	   D=DC
	   E=EC
	   V=VC
	   Z=ZC
	GO TO 460
450	B=BW
	!=CW
	D=DW
	E=EW
	V=VW
	Z=ZW
	IF (Z.LE.0.OR.V.LE.0)THEN
	    WRITE(*,4000)
4000	    FORMAT(' OUT OF BOUNDS ERROR FOR ALOG IN sub MIXFLU')
	    RETURN
	ENDIF
460	Y=B/(4.0D0*V)
	FCP=((8.0D0*Y-9.0D0*Y*Y+3.0D0*Y**3)/((1.0D0-Y)**3))-(DLOG(Z))
	FCP=FCP-(C/(RT*(V+B)))-(D/(RT*V*(V+B)))
	FCP=FCP-(E/(RT*V*V*(V+B)))+((C/(RT*B))*(DLOG(V/(V+B))))
	FCP=FCP-(D/(RT*B*V))+((D/(RT*B*B))*(DLOG((V+B)/V)))
	FCP=FCP-(E/(RT*2.0D0*B*V*V))+(E/(RT*B*B*V))
	FCP=FCP-((E/(RT*B**3))*(DLOG((V+B)/V)))
	FCP=DEXP(FCP)
	IF(L.NE.1)GO TO 470
	FKCP=FCP
	GO TO 499
470	FKWP=FCP
499	RETURN
	END
!
!
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
	SUBROUTINE FMIX(RT,CIJ,DIJ,EIJ,XC,XW,BM,FCCM,FCWM,FCM,VM,ZM,CM,DM,EM)

	IMPLICIT REAL*8(A-H,P,R-Z)
	COMMON /MVAR/BC,CC,DC,EC,BW,CW,DW,EW,J,L,M
!	CALCULATE FUG COEFFS OF H2O AND CO2 IN MIXTURE
	B=BM
	V=VM
	Z=ZM
	!=CM
	D=DM
	E=EM
	Y=B/(4.0D0*V)
	IF(M.NE.1)GO TO 650
	B1=BC
	C1=CC
	D1=DC
	E1=EC
	X1=XC
	X2=XW
	GO TO 660
650	B1=BW
	C1=CW
	D1=DW
	E1=EW
	X1=XW
	X2=XC
	IF (Z.LE.0.OR.V.LE.0)THEN
	    WRITE(*,4001)
4001	    FORMAT(' OUT OF BOUNDS ERROR FOR ALOG IN sub MIXFLU')
	    RETURN
	ENDIF
660	FCM=(4.0D0*Y-3.0D0*Y*Y)/((1.0D0-Y)**2)
	FCM=FCM+((B1/B)*((4.0D0*Y-2.0D0*Y*Y)/((1.0D0-Y)**3)))
	FCM=FCM-(((2.0D0*C1*X1+2.0D0*CIJ*X2)/(RT*B))*(DLOG((V+B)/V)))
	FCM=FCM-((C*B1)/(RT*B*(V+B)))
	FCM=FCM+(((C*B1)/(RT*B*B))*(DLOG((V+B)/V)))
	FCM=FCM-((2.0D0*D1*X1+2.0D0*DIJ*X2+D)/(RT*B*V))
	FCM=FCM+(((2.0D0*X1*D1+2.0D0*DIJ*X2+D)/(RT*B*B))*(DLOG((V+B)/V)))

	FCM=FCM+((D*B1)/(RT*V*B*(V+B)))
	FCM=FCM+((2.0D0*B1*D)/(RT*B*B*(V+B)))
	FCM=FCM-(((2.0D0*B1*D)/(RT*(B**3)))*(DLOG((V+B)/V)))
	FCM=FCM-((2.0D0*E1*X1+2.0D0*EIJ*X2+2.0D0*E)/(RT*2.0D0*B*V*V))
	FCM=FCM+((2.0D0*E1*X1+2.0D0*EIJ*X2+2.0D0*E)/(RT*B*B*V))
	FCM=FCM-(((2.0D0*E1*X1+2.0D0*EIJ*X2+2.0D0*E)/(RT*(B**3)))*(DLOG((V+B)/V)))

	FCM=FCM+((E*B1)/(RT*2.0D0*B*V*V*(V+B)))
	FCM=FCM-((3.0D0*E*B1)/(RT*2.0D0*B*B*V*(V+B)))
	FCM=FCM+(((3.0D0*E*B1)/(RT*(B**4)))*(DLOG((V+B)/V)))
	FCM=FCM-((3.0D0*E*B1)/(RT*(B**3)*(V+B)))
	FCM=FCM-(DLOG(Z))
	FCM=DEXP(FCM)
	IF(M.NE.1)GO TO 670
	FCCM=FCM
	GO TO 699
670	FCWM=FCM
699	RETURN
	END

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      subroutine cohfit (T,P,fh2o,fco2,fo2)
!----------------------------------------------------------------------
      implicit double precision (a-z)

! subroutine which returns ln(f2o,fco2) from functions fit to the
! values obtained from cohhyp for fo2 => xh2o max.

! the functions are of the form:

! ln(f) = a*t + b*p + !/t/t + d*p*t + e*p*p + f*t*t + g*sqrt(p*t)
!       + h*t**3 + i*p**3 + j*p/t/t + k ln t + l ln p 
!       + m/p**2 + n*p*p/t + o*p/t + pp*t/p + q*t*t/p 
!       + rp*t*ln(t) + s*p*ln(t) + tp

!      common / cst11 /fh2o,fco2/ cst5 /p,t,xc,u1,u2,tr,pr,rj,ps
 
      save ah,bh,ch,dh,eh,fh,gh,hh,ih,jh,kh,lh,mh,nh,oh,ph,qh,rh,sh,th

      save ac,bc,cc,dc,ec,fc,gc,hc,ic,jc,kc,lc,mc,nc,oc,pc,qc,rc,sc,tc

      save ao,bo,co,do,eo,fo,go,ho,io,jo,ko,lo,mo,no,oo,po,qo,ro,so,to

      data th,ah,bh,ch,dh,eh,fh,gh,hh,ih,jh,kh,lh,mh,nh,oh,ph,qh,rh,sh/ &
     &-121.5649    , -.5498738d-01, -.8052569d-02, -169883.7    , &
     &-.5150895d-06, 0.2909217d-08, 0.1484028d-04, 0.1595388d-02, &
     &-.2286740d-08, -.4611758d-14, -118.8630    ,  24.90956    , &
     &-1.283717    ,  126339.8    , -.1875289d-05,  1.141895    , &
     &-1.661033    , 0.7550874d-03, 0.9051986d-03, 0.1093571d-02/

      data tc,ac,bc,cc,dc,ec,fc,gc,hc,ic,jc,kc,lc,mc,nc,oc,pc,qc,rc,sc/ &
     -68.24622    , -.2797826d-01, -.5658539d-02, -221752.4, &
     -.2227993d-06, -.4785067d-08, -.2949820d-05, -.3942711d-02, &
     0.8136084d-09, 0.6607593d-13, -126.2944    ,  12.42835, &
     -.1328584    , -168530.0    , -.1849930d-05,  1.058393 ,&
      2.131351    , -.9849674d-03, 0.3118428d-02, 0.8233771d-03/
 
      data to,ao,bo,co,do,eo,fo,go,ho,io,jo,ko,lo,mo,no,oo,po,qo,ro,so/ &
     &-804.2316    , -.1652445    , -.5376252d-02, -4037433.d0 ,&
     &-.2091203d-06, -.4638105d-08, 0.3753368d-04, -.3853404d-02,&
     &-.5442896d-08, 0.6484263d-13, -121.6754    ,  127.5998    ,&
     &-.1486220    , -164866.6    , -.1863209d-05, 0.9622612    ,&
     &  2.097447   , -.9838123d-03, 0.3077560d-02, 0.7829503d-03/

      lp = dlog(p)
      lt = dlog(t)
      spt = dsqrt(p*t)
      p2 = p*p
      t2 = t*t
  
      fh2o = th + t*(ah + dh*p + t*(fh + hh*t) + (ph+qh*t)/p + rh*lp) &
     &          + p*(bh + p*(eh + ih*p) + sh*lt) &
     &          + p/t*(jh/t + nh*p + oh) &
     &          + kh*lt + lh*lp + ch/t2 + gh*spt + mh/p2


      fco2 = tc + t*(ac + dc*p + t*(fc + hc*t) + (pc+qc*t)/p + rc*lp) &
     &          + p*(bc + p*(ec + ic*p) + sc*lt) &
     &          + p/t*(jc/t + nc*p + oc) &
     &          + kc*lt + lc*lp + cc/t2 + gc*spt + mc/p2


      fo2  = to + t*(ao + do*p + t*(fo + ho*t) + (po+qo*t)/p + ro*lp) &
     &          + p*(bo + p*(eo + io*p) + so*lt) &
     &          + p/t*(jo/t + no*p + oo)&
     &          + ko*lt + lo*lp + co/t2 + go*spt + mo/p2


      end

