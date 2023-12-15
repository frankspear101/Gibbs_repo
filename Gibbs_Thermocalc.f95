! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE AllKduTPX(ivol,izero)
!	This routine is for backward compatibility
!	All it does is calculate all thermodynamic parameters and derivatives for all phases in the assemblage
	implicit none
	include "Assemb.inc"
	integer*4 ivol,izero,kcur,k
	TK = TC + 273.15d0
	call CalculateCPTemp(TK)
	Do 9000 kCur = 1,numPh
	k = asmCurrent(kCur)
	call duTPX(k,ivol,izero)
	if(izero.gt.0)then
		write(*,*)' izero>0 in Sub duTPX called from AllKduTPX'
		return
		endif
	call dlnAdX(k,izero)
	if(izero.gt.0)then
		write(*,*)' izero>0 in Sub dlnAdX called from Sub AllKduTPX'
		return
		endif
9000	continue
	return
	end	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE CalculateCPTemp(TK)
	implicit none
! ***********/CPTemp/************************************************
      REAL*8 TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS
      REAL*8 TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
      common /CPTemp/TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS,&
     &               TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
! *******************************************************************
	real*8 TK,TR
	data TR/298.15D0/
!     Temperature terms for integration of heat capacity for entropy
	TACPS = DLOG(TK/TR)
	TBCPS = (1.D0/Dsqrt(TK)-1.D0/Dsqrt(TR))
	TCCPS = (1.D0/(TK*TK) - 1.D0/(TR*TR))
	TDCPS = (1.D0/(TK*TK*TK)-1.D0/(TR*TR*TR))
	TECPS = (1.D0/TK - 1.D0/TR)
	TFCPS = (TK-TR)
	TGCPS = (TK*TK - TR*TR)
!     Temperature terms for integration of heat capacity for enthalpy
	TACPH = (TK-TR)
	TBCPH = (Dsqrt(TK)-Dsqrt(TR))
	TCCPH = ((1.D0/TK) - (1.D0/TR))
	TDCPH = (1.D0/(TK*TK) - 1.D0/(TR*TR))
	TECPH = DLOG(TK/TR)
	TFCPH = (TK*TK - TR*TR)
	TGCPH = (TK*TK*TK - TR*TR*TR)
	return
	end	

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE duTPX(k,ivol,izero)
      implicit none
!     SUBROUTINE TO COMPUTE G, H, S, V and Act, FOR PHASE components

!     ivol is a switch:  ivol = 0 compute everything (v, s, dudTPX)
!                        ivol = 1 compute volumes only (this avoids problems with phase components with negative mole fractions)

!     izero is an error flag:  izero = 0 no error
!                              izero = 1 returns indicates that a mole fraction is <= 0

!     Heat Capacity expression is that of Berman and Brown:
!     CP = A + B*T**-0.5 + C*T**-2 + D*T**-3
!            + E*T**-1   + F*T     + G*T**2
!     UBC data base uses first 4 terms (A,B,C,and D)
!     Maier-Kelly equation uses A, F and C terms Cp= A + F*T + C*T**-2
!     G term is rarely used
!     HP Cp expression used a subset

!     INTEGRATED form :  integral (Cp/T)dT:
!     S(TK,1) = Szero + 
!             + a       *ln(TK/Tr)
!             - B*2*    (1/TK**0.5 - 1/Tr**0.5)
!             - (C/2)*  (1./(TK*TK) - 1./(TR*TR))
!             - (D/3.)* (1./(TK*TK*TK)-1./(TR*TR*TR))
!             - E*      (1/TK - 1/Tr)
!             + F*      (TK-TR)
!             + (G/2)*  (TK*TK - Tr*Tr)

!     INTEGRATED form :  integral CpdT:
!     H(TK,1)) = hPhCoZero + 
!             + a*      (TK-TR)
!             + B*2*    (TK**0.5 - Tr**0.5)
!             - C*      (1/TK - 1/Tr)
!             - (D/2)*  (1./(TK*TK) - 1./(TR*TR))
!             + e*      ln(TK/Tr)
!             + (F/2)*  (TK*TK - Tr*Tr)
!             + (G/3)*  (TK*TK*TK - Tr*Tr*Tr)
!                     
!      Volume at T and P is computed from equation of Berman (1988)
!     Volume terms in Gibbs are       v1, v2, v3 and v4 
!     These correspond to
!     in Berman's program             V1, V2, V3 and V4
!     in Berman (1988) (the paper)
!       the correspondence is         V3, V4, V1 and V2
!   VOLUME EQUATION: (Bermans program)                                                                                  
!    V(P,T)/V(1,298) = 1 + v1PhCo(T-298) + v2PhCo(T-298)**2 + v3PhCo(P-1) + v4PhCo(P-1)**2                           

!         vattp(TK,P) = vzero * (1.D0 + v1 *(TK-298.15D0) +
!     &       v2 *(TK-298.15D0)**2 + v3 *(PB-1.0D0) +
!     &       v4 *(PB-1.0D0)**2)

!     The integral of VdP (dG = VdP) is 
!      int_VdP = vzero ((v3/2.0) - v4)*(PB**2 - PR**2) + (v4/3.0)*(PB**3 - 1.**3)
!             + ((1. - v3 + v4 + v1PhCo(TK - Tr) + v2*(TK - Tr)**2) * (PB - 1.0))

!     to correct enthalpy and entropy up to pressure (at TK) we need the
!     integral of (V*alpha)dP     (alpha is coeff of thermal expansion)
!     int_Valpha_dP = V*(alpha)dP
!     int_Valpha_dP = Vzero [v1 + 2*v2*(TK - TR)]*(PB - 1.)


!     Pressure correction to the entropy 
!      dS = - V*(alpha)dP = -int_Valpha_dP

!     Pressure correction for enthalpy
!     dH = VdP - V*(alpha)*T dP = int_VdP - TK * int_Valpha_dP

!+++++++++++++++++++++++++++++++++++++++++++
!     definitions of some of the variable arrays

!     lnAct(j) = ln(aj)

!     This subroutine underwent major modifications in January 1993 by F. Spear to
!     incorporate a new way of storing data that should make the bookkeeping easier
!     It also will permit multisite solutions to be modeled

!     The basic idea is this.  The governing thermodynamic equations are contained in
!     a matrix AA in subroutine COMPUTE.  Each cell in AA contains the sum of the
!     partial derivative of the chemical potentials of the phase components multiplied
!     by the stoichiometric coefficient of the phase component in the reaction.
!     
!     The purpose of the present subroutine is to generate the matrix for the
!     partial derivatives for each chemical potential with respect to T, P and each
!     independent mole fraction.  The resulting matrix has dimensions
!     npct x nvar (number of phase components x number of independent variables) 
!     
!     The partial derivatives are:

!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.

!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.

!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X

!     Note that the partial derivatives of a component with respect to the mole fractions of components in
!     other phases is zero so the matrix is sparsly populated.

!     The exact form of the partial derivatives depends on the type of solution model.  

!     For crystals, we follow the Wood and Nicholls formulation with Berman's generalized
!     margules model for each site.  
!	All partial derivatives in Gibbs3 are calculated using finite difference approximations
!	   for dµ/dT, dµ/dP and dµ/dX
!     Variable names:
!     dudTPX(25,15) contains the partial derivatives for each chemical potential (first subscript)
!                 with respect to each independent variable (second subscript).


! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
! ****************************************
! ***********/CPTemp/************************************************
      REAL*8 TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS
      REAL*8 TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
      common /CPTemp/TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS,&
     &               TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
! *******************************************************************

!     LOCAL VARIABLES
	integer*4 i,j,k,L,izero,iflag,ivol
	REAL*8 R,RTK,TR,PR,TT,DISORDS,DISORDH,ppp,gKyanite,gSillimanite,gAndalusite,	&
		GH2O,HH2O,SH2O,VH2O,gplus,gminus, 					&
		Tplus,Tminus,Pplus,Pminus
	real*8 temp,xH2O,xCO2
	real*8 Xtemp(PhCoMax),uplus(PhCoMax),uminus(PhCoMax)


! *******************************************************************
! 	HP98 CORK fluid variables
!      integer*4 H2Oswitch	! for HoPo fluid 1=pure H2O; 2=pure CO2; 3=H2O-CO2 mix
!	real*8 fh2opure,fco2pure
	real*8 fco2plus,fco2minus,fh2oplus,fh2ominus,fh2o,fco2,VH2Oplus,VH2Ominus

! *********************************************************************    
!     ab Quartz variables
      integer*4 AorB
      Real*8 GLambda,Hlambda,SLambda,GLambda_plus,GLambda_minus,VLambda,	&
      	HLambda_plus,HLambda_minus,SLambda_plus,SLambda_minus

!     Kspar disorder terms
      Real*8 KSPD0,KSPD1,KSPD2,KSPD3,KSPD4,KSPD5,KSPTmin,KSPTmax
      data KSPD0,KSPD1,KSPD2,KSPD3,KSPD4,KSPD5,KSPTmin,KSPTmax			&
       	/282.98291D0,-4831.37500D0,3620706.00000D0,0.00000,-0.15733D0,		&
      	0.00003D0,298.15000D0,1436.15000D0/
! *********************************************************************    
      DATA R/8.3144D0/
      data TR,PR/298.15D0,1.0D0/
! *********************************************************************    

      TK=TC+273.15D0
      RTK = R*TK
      
!     Check for temperature between reasonable bounds
      if(TK.LT.200.or.TK.GT.2000)then
         WRITE(*,*)' *****************************************'
         WRITE(*,*)' Temperature is out of bounds.  T(K) =',TK
         WRITE(*,*)' Temperature limits are 200-2000(K)'
         WRITE(*,*)' Check made in subroutine dudTPX'

         izero=1
         PAUSE 'Hit return to continue...'
         return
         endif

!     Check to see what type of mineral this is

!	These values are stored in array PhaseType(k)
!			1 = 1 component SPaC
!			2 = alpha-beta quartz (SPaC)
!			3 = multicomponent SPaC
!			4 = H2O SPaC (Haar)
!			5 = Kspar (order-disorder)
!			8 = fixed Al2SiO5 minerals (KAS) - should work for all 3 dataset types
!			10 = Graphite saturated fluid (Connolly_and_Cesare, 1993 JMG)
!			101 = 1 component HP98
!			103 = multicomponent HP98 (includes SF input)
!			104 = Fluid HP98 (CorK)
!			105 = Melt HP98
!			201 = 1 component HP11
!			203 = multicomponent HP11 (includes SF input)
!			204 = Fluid HP11 (Pitzer eq)
!			205 = Melt HP11


	select case(PhaseType(k))
	case(1,101,201)
! ----------------------------------------------
!     		1 component phases
		call GHSV(K,ivol,iflag)
		if(iflag.eq.1) go to 9900      
		if(ivol.eq.1)  go to 9000
      		J=1
		dudTPX(k,j,1) = - SATTP(k,j)
		dudTPX(k,j,2) =   VATTP(k,j)
		GO TO 9000
! --------------------------------------------------
	Case(2)
!     		THIS SECTION FOR alpha-beta QUARTZ from Berman's database
      		J=1
!     		call subroutine LabQtz to get extra H, S and G from lambda transition
!     		Parameter AorB = 1 if in alpha field
!               AorB = 2 if in beta field
!     		Calculate implicit B=Vlambda from finite difference approximation
      		call Labqtz(AorB,TK,PB+1.0D0,Hlambda_plus,Slambda_plus,Glambda_plus)
      		call Labqtz(AorB,TK,PB-1.0D0,Hlambda_minus,Slambda_minus,Glambda_minus)
      		Vlambda =    ((Hlambda_plus - Hlambda_minus) - TK*(SLambda_plus - SLambda_minus))/2.0d0
!     		now calculate values at the real T and P
      		call Labqtz(AorB,TK,PB,Hlambda,Slambda,Glambda)
!     		This code is executed only for the alpha-beta quartz mineral
!     		It sets the values of the thermodynamic properties based on the type of mineral
      		L = 1       ! default is alpha field
!      		IF(AorB.eq.2) L = 2     !in Beta field (According to R Berman, we use alpha values throughout
!     		Set reference thermodynamic properties to appropriate values
		hPhCoZero(k,j)= HZEROF(L)
		sPhCoZero(k,j)= SZEROF(L)
		vPhCoZero(k,j)= VZEROF(L)
		aPhCoCp(k,j)  = ACPF(L)
		bPhCoCp(k,j)  = BCPF(L)
		cPhCoCp(k,j)  = CCPF(L)
		dPhCoCp(k,j)  = DCPF(L)
		ePhCoCp(k,j)  = ECPF(L)
		fPhCoCp(k,j)  = FCPF(L)
		gPhCoCp(k,j)  = GCPF(L)
		v1PhCo(k,j)   = v1F(L)
		v2PhCo(k,j)   = v2F(L)
		v3PhCo(k,j)   = v3F(L)
		v4PhCo(k,j)   = v4F(L)

		call GHSV(K,ivol,iflag)
		if(iflag.eq.1) go to 9900      
		! if(ivol.eq.1)  go to 9000
!     		activity
      		lnAct(k,j)=0.0D0
!     		correct for lambda transition
		VATTP(k,j) = VATTP(k,j) + VLambda
		SATTP(k,j) = SATTP(k,j) + SLambda
		HatTP(k,j) = HatTP(k,j) + HLambda
		GatTP(k,j) = HatTP(k,j) - TK*SatTP(k,j)
!		write(*,*)' inThermocalc ',TC,PB,phName(k),GatTP(1,1)
!     		temperature derivative at T and P
		dudTPX(k,j,1) = -SATTP(k,j)    
!     		pressure derivative of  at T and P
		dudTPX(k,j,2) =  VATTP(k,j)
!     		no composition derivatives for a pure phase such as quartz
		GO TO 9000
! -----------------------------------------------------
	case (3,103,105,203,205,206)
!		Multicomponent, multisite phases
	!     THIS SECTION FOR User input generalized solution models
		call GHSV(K,ivol,iflag)
		if(iflag.eq.1)then
			write(*,*)' Error in routine GHSV -- aborting calculations'
			izero = 1
			return
			endif      
		if(ivol.eq.1)go to 9000
		do 3913 J=1,numPhCo(K)
		Xtemp(j) = xPhCo(k,j)
3913    	continue
		!     T and P derivatives
		!     temperature derivatives
		call Activity(K,uplus ,Xtemp,TK+.1D0,PB,R,iflag)
		if(iflag.eq.1)go to 9900
		call Activity(K,uminus,Xtemp,TK-.1D0,PB,R,iflag)
		if(iflag.eq.1)go to 9900
		do 3930 j=1,numPhCo(K)
		temp = (uplus(j)-uminus(j))/0.2D0
		dudTPX(k,j,1) = -SATTP(k,j) + temp
3930    	continue
		!     pressure derivatives
		call Activity(K,uplus ,Xtemp,TK,PB+1.0D0,R,iflag)
		if(iflag.eq.1)go to 9900
		call Activity(K,uminus,Xtemp,TK,PB-1.0D0,R,iflag)
		if(iflag.eq.1)go to 9900
		do 3931 j=1,numPhCo(K)
		temp = (uplus(j)-uminus(j))/2.0D0
		dudTPX(k,j,2) = VATTP(k,j) + temp
3931		continue
		GO TO 9000
! ------------------------------------------------------------------

! -------------------water------------------------------------------------
!     THIS SECTION IS FOR CALCULATING WATER fugacities, activities
!     AND ENTROPIES - HAAR equation
	case(4)
		j = 1
		PPP=PB
		if(PFluidSwitch.ne.0)PPP=PFLUID
	!     	THIS FOR PURE WATER
		call WHAAR2 (TK,PPP,GH2O,HH2O,SH2O,VH2O,izero)
		if(izero.eq.1)return
		VATTP(k,j) = VH2O
		SATTP(k,j) = SH2O
		HATTP(k,j) = HH2O
		GATTP(k,j) = GH2O
	!     	VOLUME IS RETURNED IN joule/BAR
		VMOL(K)=VH2O
	!     	activity
		lnAct(k,j) = 0.0D0
	!     	temperature derivative at T and X (Presumably the entropy is at P and T from this routine)
		dudTPX(k,j,1)= -SATTP(k,j)    
	!     	pressure derivative of  at T and P (note, ideal solution assumes no excess volume of mixing)
		dudTPX(k,j,2)=VATTP(k,j)
	!     	no composition derivatives for a pure phase like pure water
		GO TO 9000
! ------------------------------------------------------------------
	case(5)
!     		this section for K-feldspar, which includes disorder terms
		J=1
		L=1
		call GHSV(K,ivol,iflag)
		if(iflag.eq.1) go to 9900      
		if(ivol.eq.1)  go to 9000
!     		UBC DISORDER equation
!      		Cp(disorder) = D0 + D1/SQRT(T) + D2/T/T +            D3/T + D4*T + D5*T*T
!      		Cp =           K0 + K1/SQRT(T) + K2/T/T + K3/T/T/T + K4/T + K5*T + K6*T*T                                
      		TT=TK
      		IF(TK.GT.KSPTmax) TT = KSPTmax
!     		KSPD0,KSPD1,KSPD2,KSPD3,KSPD4,KSPD5,KSPTmin,KSPTmax
		DISORDH =       KspD0*(TT-TR) 			&
		         + 2.D0*KspD1*(Dsqrt(TT) - Dsqrt(TR))   &
			      - KspD2*((1.D0/TT) - (1.D0/TR))   &
			      + KspD3*DLOG(TT/TR) 		&
			     + (KspD4/2.0D0)*(TT*TT - TR*TR) 	&
			     + (KspD5/3.D0) *(TT*TT*TT - TR*TR*TR)
	      	Hattp(k,j) = Hattp(k,j) + disordH
      		DISORDS       =   KspD0*DLOG(TT/TR) 			&
                          - 2.D0*(KspD1)*(1.D0/sqrt(TT)-1.D0/sqrt(TR)) 	&
                         - 0.5D0*(KspD2)*(1.D0/(TT*TT) - 1.D0/(TR*TR)) 	&
                               - (KspD3/3.D0)*(1.D0/TT - 1.D0/TR) 	&
                               + (KspD4)* (TT-TR) 			&
                               + (KspD5/2.0D0)*  (TT*TT - TR*TR)
		sattp(k,j)= sattp(k,j) + disordS
		GatTP(k,j) = GatTP(k,j) + DISORDH - TK * disords
!     		activity
      		lnAct(k,j)=0.D0
!     		temperature derivative at T and X (note, no pressure correction made to entropy)
      		dudTPX(k,j,1)= -SATTP(k,j)    
!     		pressure derivative of  at T and P (note, ideal solution assumes no excell volume of mising)
      		dudTPX(k,j,2)=VATTP(k,j)
!      		no composition derivatives for a pure phase
      		GO TO 9000
! -------------------------------------------------------------
	case(8,108,208)
	!     	This section for kyanite-sillimanite-andalusite polymorphs
	! 	Find the Al2SiO5 polymorph with the lowest G at the given T and P
	! 	Then use that phase in all the calculations
	! 	Kyanite
		call SetFixedMineral(k,'kyanite')
		call GHSV(K,ivol,iflag)
		gKyanite = gattp(k,1)
		call SetFixedMineral(k,'sillimanite')
		call GHSV(K,ivol,iflag)
		gSillimanite = gattp(k,1)
		call SetFixedMineral(k,'andalusite')
		call GHSV(K,ivol,iflag)
		gAndalusite = gattp(k,1)
	! 	Now pick phase with lowest G
		if(gAndalusite.lt.gKyanite.and.gAndalusite.lt.gSillimanite)then
			! andalusite has lowest G
			call SetFixedMineral(k,'andalusite')
	! 		write(12,*)'Stable phase is andalusite'
		else if(gkyanite.lt.gAndalusite.and.gKyanite.lt.gSillimanite)then
			! kyanite has lowest G
			call SetFixedMineral(k,'kyanite')
	! 		write(12,*)'Stable phase is kyanite'
		else
			! sillimanite has lowest G
			call SetFixedMineral(k,'sillimanite')
	! 		write(12,*)'Stable phase sillimanite'
		endif
		call GHSV(K,ivol,iflag)
	!     	activity
		lnAct(k,1)=0.D0
	!     	temperature derivative at T and X
		dudTPX(k,1,1)= -SATTP(k,1)    
	!     	pressure derivative of  at T and P (note, ideal solution assumes no excess volume of mixing)
		dudTPX(k,1,2)=VATTP(k,1)
	!     	no composition derivatives for a pure phase
		go to 9000
! ------------------------------------------------------------------
! ------------------------------------------------------------------
! ------Graphite-COH fluid H:O = 2:1 Connolly and Cesare, 1993 JMG ------------------------------------------------------------
	case(10)
		j = 1
		call graphite_sat(TK,pb,Gh2o)
		! 	finite difference entropy
		call graphite_sat(TK+.1,pb,Gplus)
		call graphite_sat(TK-.1,pb,Gminus)
		sh2o = -(gplus - gminus)/0.2d0	
		! 	finite difference volume
		call graphite_sat(tk,pb+1,gplus)
		call graphite_sat(tk,pb-1,gminus)
		vh2o = (gplus - gminus)/2.d0	
		hh2o = gh2o + tk*sh2o	
		! 	write(*,*)'G,H,S,V ',Gh2o,Hh2o,Sh2o,Vh2o
		GatTP(k,j) = Gh2o
		SATTP(k,j) = Sh2o
		VATTP(k,j) = vh2o
		VMOL(K)    = vh2o
		HATTP(k,j) = hh2o
		!     	activity
		lnAct(k,j) = 0.0D0
		!     	temperature derivative at T and P
		dudTPX(k,j,1) = -SATTP(k,j)    
		!     	pressure derivative of  at T and P
		dudTPX(k,j,2) = VATTP(k,j)
		!     	no composition derivatives because we are treatin water as a pure phase here
		! 	This is not strictly correct, because we have a COH fluid
		! 	By doing it this way, we are forcing the variance to be correct, and not worrying
		! 	about the composition of the fluid
		! 	Hence, this can only be used for a carbonate-free system (e.g. a pelite)
		GO TO 9000
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
	case(104)
! 	CORK fluid
		j = 1
		PPP=PB
		if(PFluidSwitch.ne.0)PPP=PFLUID
		! 	we calculate volume as dG/dP, because CORK returns only fugacity
		! 	so we set ivol = 0 here to force GHSV to return G at (1,T).
		call GHSV(K,0,iflag)	
		! 	call GHSV(K,ivol,iflag)
		if(iflag.eq.1)return      
		if(numPhCo(K).eq.2)then
			XCO2=xPhCo(k,2)
			XH2O=1.0D0-XCO2
	! 		H2Oswitch = 3
		else        !either pure H2O or pure CO2
			if(UsePhCo(k,1).eq.1)then    			! we are using pure H2O
				XH2O = 1.0d0
				XCO2 = 0.0d0
		! 		H2Oswitch = 1
				j = 1
				call Gat1T(TK,k,j,GH2O,HH2O,SH2O,iflag)
				call hprk(ppp,tk,xco2,fh2o,fco2)	! returns ln fugacity
				GatTP(k,j) = GH2O + RTK*fh2o
				lnAct(k,j) = 0.0d0
				Tplus = TK+0.1d0
				call Gat1T(Tplus,k,j,gplus,HH2O,SH2O,iflag)
				call hprk(ppp,Tplus,xco2,fh2oplus,fco2plus)
				gplus = gplus + R*Tplus*fh2oplus
				Tminus = TK - 0.1d0
				call Gat1T(Tminus,k,j,gminus,HH2O,SH2O,iflag)
				call hprk(ppp,Tminus,xco2,fh2ominus,fco2minus)
				gminus = gminus + R*Tminus*fh2ominus
				SatTP(k,j) = -(gplus-gminus)/0.2d0
				HatTP(k,j) = GatTP(k,j) + TK*SatTP(k,j)
				Pplus = ppp + 1.0d0
				call hprk(Pplus,tk,xco2,fh2oplus,fco2plus)
				gplus = GatTP(k,j) + RTK*fh2oplus
				Pminus = ppp-1.0d0
				call hprk(Pminus,tk,xco2,fh2ominus,fco2minus)
				gminus = GatTP(k,j) + RTK*fh2ominus
				VatTP(k,j) = (gplus-gminus)/2.0d0
				VMOL(K)    = VatTP(k,j)
				dudTPX(k,j,1) = -SatTP(k,j)
				dudTPX(k,j,2) =  VatTP(k,j)
				elseif(UsePhCo(k,2).eq.1)then  		! we are using pure CO2
				XH2O = 0.0d0
				XCO2 = 1.0d0
		! 		H2Oswitch = 2
				endif
		endif
		go to 9000
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
	case(204)
! 	Pitzer and Sterner 1994 EOS for H2O and CO2 model
!	Used in Holland and Powell 2011 dataset
!	Special case = 15
!	Only tested for H2O so far
		j = 1
		PPP=PB
		if(PFluidSwitch.ne.0)PPP=PFLUID
	! 	we calculate volume as dG/dP, because CORK returns only fugacity
	! 	so we set ivol = 0 here to force GHSV to return G at (1,T).
		call GHSV(K,0,iflag)	
	! 	call GHSV(K,ivol,iflag)
		if(iflag.eq.1)return      
		if(numPhCo(K).eq.2)then
			XCO2=xPhCo(k,2)
			XH2O=1.0D0-XCO2
			write(*,*)' Mixed volatile H2O-CO2 not yet implemented --- aborting'
			pause 'aborting program'
			stop

			else        !either pure H2O or pure CO2
			if(UsePhCo(k,1).eq.1)then    			! we are using pure H2O
				XH2O = 1.0d0
				XCO2 = 0.0d0
		! 		H2Oswitch = 1
				j = 1
				if(Pb.lt.100)then
					write(*,*)' Pressure is too low ... abort. Pb = ',Pb
					izero = 1
					go to 9000
					endif
				call PS94H2O(Pb,TK,GH2O,VH2O)
				GatTP(k,j) = GH2O		! Gibbs uses J (not kJ)
				VatTP(k,j) = VH2O
				lnAct(k,j) = 0.0d0
	!			Temperature derivative (to get SH2O)
				Tplus = TK+0.1d0
				call PS94H2O(Pb,Tplus,Gplus,VH2Oplus)
				Tminus = TK - 0.1d0
				call PS94H2O(Pb,Tminus,Gminus,VH2Ominus)
				SatTP(k,j) = -(gplus-gminus)/0.2d0
				HatTP(k,j) = GatTP(k,j) + TK*SatTP(k,j)
				VMOL(K)    = VatTP(k,j)
				dudTPX(k,j,1) = -SatTP(k,j)
				dudTPX(k,j,2) =  VatTP(k,j)

				elseif(UsePhCo(k,2).eq.1)then  		! we are using pure CO2
				XH2O = 0.0d0
				XCO2 = 1.0d0
				write(*,*)' CO2 not yet implemented --- aborting'
				pause 'aborting program'
				stop
				endif
			endif
		go to 9000

	case default
		write(*,*)' The PhaseType(k) was not correct (no code found). In routine duTPX'
		write(*,*)' In subroutine duTPX '
		write(*,*)'K,PhaseName(k),PhaseType(k) = ',k,PhName(k),PhaseType(k)
		pause 'Hit return to continue'
	end select
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
9000  CONTINUE
!     Now the thermodynamic derivatives for all phases are computed
      if(iLong(3).eq.1) then
	      write(12,*)' '
	      write(12,*) 'Phase component partial derivatives - in duTPX'
	      write(12,*)'                      dudT          dudP          dudX2         dudX3         dudX4...'
!		Do 9010 kCur = 1,numPh
!		k = asmCurrent(kCur)
       		DO 9010 j=1,numPhCo(K)
      		write(12,9020)j,phCoName(k,j),(dudTPX(k,j,L),L=1,numPhCo(K)+1)
9010    	continue
9020    	FORMAT(i4,1x,A8,4x,15e14.4)
       		endif
      RETURN
! -----------------
!     error trap when components of phases have negative values
9900  continue
	if(iLong(1).eq.1)then

	         WRITE(*,*)' *****************************************'
	         WRITE(*,*)' A phase component has a zero or negative value'
	         WRITE(*,*)k,phname(k)
	         write(*,'(6A12,2x)')(phCoName(k,i),i=1,numPhCo(k))
	         write(*,'(6E14.4)')(xPhCo(k,i),i=1,numPhCo(k))
		 write(*,*)'Tc,Pb= ',TC,PB
	         PAUSE 'Hit return to continue...'
	endif
         izero=1
         return
      END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine GHSV(K,ivol,iflag)

!     This routine computer G, H, S and V at T and P for each phase component
!     and also computes the molar volume of the phase

      implicit none
      
! ****************************************
	include "Assemb.inc"
! ****************************************
! ***********/CPTemp/************************************************
      REAL*8 TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS
      REAL*8 TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
      common /CPTemp/TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS,&
     &               TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
! *******************************************************************
!     LOCAL VARIABLES
      integer*4 j,k,iflag,ivol
      REAL*8 R,TR,PR
      REAL*8 H_temp,S_temp,V_temp,katT,onethird,dkdT,V_at_1bar_T
      real*8 Vmax,Smax,Tcrit,TcatP,Q4298,Q298,Q4,Q,Glandau,H298,S298,VexatT,intVexatTdP,int_VdP_temp		! variables for Landau calculations
      real*8 Gex,Hex,Sex,Vex
      real*8 ao,ko,kprime,theta,aa,bb,cc,mu,muo,psi,psio,Pth,alphahp11,oneminusCC,X,Y2,xHi,xLo,y2Hi,y2Lo
      real*8 Kpp,Gam,yy,Qosquared,Qsquared,VLand,H_tot,S_tot,dH,dV,W,Wv,Hdis,Vdis,Sdis,Xa1,Xa2,Xs1,Xs2,n,sf
! *******************************************************************

      DATA R/8.3144D0/
      data TR,PR/298.15D0,1.0D0/
      data onethird/0.33333333333333333d0/
! *******************************************************************

      iflag=0           ! 0 = no problems
!     compute volume of the phase at T and P
!     Note that this subroutine also computes the volumes of each phase component at P and T

!	DatasetKey = 0							! This will tell whether the mineral type is SPaC=1, HP98=2, HP11 = 3
!	if(PhaseType(K).gt.0.and.PhaseType(K).LT.100)DatasetKey = 1
!	if(PhaseType(K).ge.100.and.PhaseType(K).LT.200)DatasetKey = 2
!	if(PhaseType(K).ge.200.and.PhaseType(K).LT.300)DatasetKey = 3
!	if(DatasetKey.eq.0)then
!		write(*,*)' Problem with DatasetKey in subroutine ReadMinData'
!		pause ' Hit return to continue'
!		endif

!     MOLAR  VOLUME CALCULATION
      	VMOL(K)=0.D0

!	T terms are calculated in subroutine CalculateCPTemp and passed through in
!	common block CPTemp

	do 10 j = 1,numPhCo(k)		! loop through all phase components
	
         H_temp = hPhCoZero(k,j) &
                     + aPhCoCp(k,j)        *TACPH &
                     + 2.D0*bPhCoCp(k,j)   *TBCPH &
                     - cPhCoCp(k,j)        *TCCPH &
                     - (dPhCoCp(k,j)/2.D0) *TDCPH &
                     + ePhCoCp(k,j)        *TECPH &
                     + (fPhCoCp(k,j)/2.0D0)*TFCPH &
                     + (gPhCoCp(k,j)/3.D0) *TGCPH

!         Temperature correction for entropy
          S_temp  =  sPhCoZero(k,j) 		 &
                     + aPhCoCp(k,j)       *TACPS &
                     - 2.D0*bPhCoCp(k,j)  *TBCPS &
                     - 0.5D0*cPhCoCp(k,j) *TCCPS &
                     - (dPhCoCp(k,j)/3.D0)*TDCPS &
                     - ePhCoCp(k,j)       *TECPS &
                     + fPhCoCp(k,j)       *TFCPS &
                     + (gPhCoCp(k,j)/2.D0)*TGCPS



	select case (DatasetKey(k))
	
	case(1)			! SPaC dataset

		! Berman volume equation
         	vattp(k,j)=vPhCoZero(k,j)*(1.D0 + v1PhCo(k,j)*(TK-TR) +   &
            		v2PhCo(k,j)*(TK-TR)**2 + v3PhCo(k,j)*(PB-PR) +    &
            		v4PhCo(k,j)*(PB-PR)**2)

		int_VdP(k,j) = vPhCoZero(k,j) 					&
		 * (( (v3PhCo(k,j)/2.0D0) - v4PhCo(k,j) )*(PB**2 - PR**2) 	&
		 + (v4PhCo(k,j)/3.0)*(PB**3 - PR**3) 				&
		 + ((1. - v3PhCo(k,j) + v4PhCo(k,j) + v1PhCo(k,j)*(TK - Tr) 	&
				     + v2PhCo(k,j)*(TK - Tr)**2) * (PB - PR)))

		!      V*(alpha)dP
		int_Valpha_dP(k,j) = vPhCoZero(k,j) * 				&
		      ((v1PhCo(k,j) + 2.0D0*v2PhCo(k,j)*(TK - TR))*(PB - PR))

		!     entropy at T and P 
		SatTP(k,j) = S_temp - int_Valpha_dP(k,j)
		!     enthalpy at T and P
		HatTP(k,j) = H_temp + int_VdP(k,j) - TK * int_Valpha_dP(k,j)
		!     Gibbs free energy at T and P 
		GatTP(k,j) = H_temp - TK*S_temp + int_VdP(k,j)



	case(2)			! HoPo98 dataset
		Gex = 0.0d0
		Hex = 0.0d0
		Sex = 0.0d0
		Vex = 0.0d0
		! HoPo volume equation
		katT = v2PhCo(k,j)*(1.d0 - 1.5d-4*(TK - TR))			! HoPo 98 equation
		! Thermal expansion -- gibes V at 1bar, T
		V_temp =  vPhCoZero(k,j)* (1.0d0 + v1PhCo(k,j)* ((TK - TR) - 20.0d0*(dsqrt(TK)-dsqrt(TR))))

!		Note that when KatT = 0 (e.g. for melt components with no compressibility data)
!		this code will result in vattp = 0
!		What I want is vattp = vzero (i.e. no T or P changes. Hence the if statement here	
		if (katT.lt.1.d-10) then
			vattp(k,j) = V_temp
			else
			vattp(k,j) = V_temp*(1.0d0 - 4.0d0*PB/(katT + 4.0D0*Pb))**0.25	! form in H&P 1998
			endif

! 		code for Landau ordering - we'll do it here and add terms in where needed below
		if(LandauTc(k,j).gt.0)then
			Vmax    = LandauVmax(k,j)
			Smax    = LandauSmax(k,j)
			Tcrit   = LandauTc(k,j)
	
			TcatP   = Tcrit + (Vmax/Smax)*PB
			Q4298   = 1 - TR/Tcrit
			Q298    = Q4298**0.25d0
			Q4      = 1 - TK/TcatP

			H298 = Smax * Tcrit * (Q298*Q298 - onethird*Q298**6)
			S298 = Smax * Q298*Q298
			VexatT =  Vmax * Q298*Q298 * (1.0d0+v1PhCo(k,j)* (  (TK-TR) - 20.0d0*(dsqrt(TK)-dsqrt(TR))))
			intVexatTdP = onethird*VexatT*katT*((1.0d0+4.0d0*(PB)/katT)**0.75d0 - 1.0d0)

			if(Q4.le.0)then
				Q = 0
				Glandau = 0.0d0
				Vex = VexatT / (1.0d0 + 4.0d0*PB/katT)**0.25
				else
				Q = Q4**0.25d0
				Glandau = Smax*((TK-TcatP)*Q*Q + onethird*TcatP*Q**6)
				Vex = -1.0d0*(Vmax*(Smax*(TK+2.0d0*Tcrit)+2.0d0*PB*Vmax)*(1.0d0-(Smax*TK)/	&
						(Smax*Tcrit+PB*Vmax))**0.5d0)
				Vex = Vex/(3.0d0*(Smax*Tcrit+PB*Vmax))
		! 		Vex(J) = Vex(J) + VexatT / (1.0d0 - 4.0d0*PB/(katT + 4.0D0*Pb))**0.25
				Vex = Vex + VexatT / (1.0d0 + 4.0d0*PB/katT)**0.25
				endif
			Gex = H298 - TK*S298 + intVexatTdP + Glandau
			Sex = S298 - Smax * Q * Q
			Hex = Gex + TK*Sex
		! 	Vex(J) = VexatT * (1.0d0 - 4.0d0*PB/(katT + 4.0D0*Pb))**0.25
		! 	write(12,*)'Gex   ,   Glandau   , Hex,   Sex   , Vex'
		! 	write(12,*)Gex(J),Glandau,Hex(J),Sex(J),Vex(J)
			vattp(k,j) = vattp(k,j) + Vex
			endif

		!     VdP
		int_VdP(k,j) = (V_temp*katT/3.0d0)*((1.0d0+4.0d0*(PB)/katT)**0.75d0 - 1.0d0)
		!      V*(alpha)dP
		int_Valpha_dP(k,j) = int_VdP(k,j) * v1PhCo(k,j) * (1.0d0 - 10.0d0/sqrt(TK))
		!     entropy at T and P 
		SatTP(k,j) = S_temp - int_Valpha_dP(k,j) + Sex
		!     enthalpy at T and P
		HatTP(k,j) = H_temp + int_VdP(k,j) - TK * int_Valpha_dP(k,j) + Hex
		!     Gibbs free energy at T and P 
		GatTP(k,j) = H_temp - TK*S_temp + int_VdP(k,j)+ Gex



	case(3)		! HoPo2011

		Select case (HP11ModelSwitch(k,j))
		case (0)	! no ordering model (default)
			ao     = v1PhCo(k,j)
			ko     = v2PhCo(k,j)
			kprime = v3PhCo(k,j)
			theta  = v4PhCo(k,j)
			!Pkb    = Pb/1000.0d0

			mu     = theta/TK
			muo    = theta/TR
			psi    = (mu**2)*dexp(mu)/((dexp(mu) - 1.0d0)**2)
			psio   = (muo**2)*dexp(muo)/((dexp(muo) - 1.0d0)**2)
			! note that ao is assumed constant
			Pth    = (ao*theta*ko/psio)*(1.0d0/(dexp(mu) - 1.0d0) - 1.0d0/(dexp(muo) - 1.0d0))

!			aa     = 1.0d0 + kprime
!			bb     = (kprime/ko)*(2.0d0 + kprime)/(1.0d0 + kprime)
!			cc     = 1.0d0/(2.0d0*kprime + kprime**2)
			Kpp    = -Kprime/Ko		! K"
			aa  = (1.0d0 + Kprime) / (1.0d0 + Kprime + Kpp * Ko)
			bb  = Kprime / Ko - Kpp/(1.0d0 + Kprime)
			cc  = (1.0d0 + Kprime + Kpp * Ko) / (Kprime*(1.0d0 + Kprime) - Kpp * Ko)
			yy = 1.0d0 + bb * (Pb - Pth)
!			print *,'////////////'
!			print *,tk,pb,pth
!			print *,aa,bb,cc,yy
!			print *,'////////////'
			! calculate V at P, T
			V_temp =  vPhCoZero(k,j) * ( 1.0 - aa * (1.0-yy**(-cc)) ) 

!			VatTP(k,j) = vPhCoZero(k,j)*(1.0d0 - aa*(1.0d0 - (1.0d0 + bb*(Pb - Pth))**(-cc)))

			oneminuscc = 1.0d0 - cc
		!	int_VdP(k,j) =   Pb*VPhCoZero(k,j)*(1.0d0-aa+(aa*((1.0d0 - bb*Pth)**oneminuscc - (1.0d0 + bb*(Pb - Pth))**oneminuscc))/(bb*(cc - 1.0d0)*Pb))
			int_VdP_temp  = vPhCoZero(k,j)*(Pb*(1.0d0-aa)+aa*((1.0d0-bb*Pth)**oneminuscc-		&
					(1.0d0+bb*(Pb-Pth))**oneminuscc)/(bb*(cc-1.0d0))) 
			!	Assume alpha is independent of P and use the 1 bar expression in HP11
			!	I could also just use alpha at TR,PR, which assumes alpha is independent of both P and T
			alphahp11 = ao*(psi/psio)*(1.0d0/((1.0d0 - bb*Pth)*(aa + (1.0d0-aa)*(1.0d0 - bb*Pth)**cc)))
!			int_Valpha_dP(k,j) = alphaHP11*int_VdP(k,j)		! which one is correct???
			VatTP(k,j) = V_temp
			int_VdP(k,j) = int_VdP_temp
			int_Valpha_dP(k,j) = ao*int_VdP_temp			! which one is correct???
			!     entropy at T and P 
			SatTP(k,j) = S_temp - int_Valpha_dP(k,j)
			!     enthalpy at T and P
			HatTP(k,j) = H_temp + int_VdP_temp - TK * int_Valpha_dP(k,j)
			!     Gibbs free energy at T and P 
			GatTP(k,j) = H_temp - TK*S_temp + int_VdP_temp

		case (1)			! Landau model
			Tcrit   = HP11Model(k,j,1)
			Smax    = HP11Model(k,j,2)
			Vmax    = HP11Model(k,j,3)
			TcatP   = Tcrit + (Vmax/Smax)*PB
			ao     = v1PhCo(k,j)		! thermal expansion
			ko     = v2PhCo(k,j)		! compressibility
			kprime = v3PhCo(k,j)		! K'
			theta  = v4PhCo(k,j)		! Einstein temperature
			Kpp    = -Kprime/Ko		! K"
			!Pkb    = Pb/1000.0d0
			mu     = theta/TK
			muo    = theta/TR
			Gam    = theta*(1.0d0/(dexp(mu) - 1.0d0) - 1.0d0/(dexp(muo) - 1.0d0) )
			psi    = (mu**2)*dexp(mu)/((dexp(mu) - 1.0d0)**2)
			psio   = (muo**2)*dexp(muo)/((dexp(muo) - 1.0d0)**2)
			! note that ao is assumed constant
			Pth    = (ao*ko*Gam)/psio

			aa     = (1.0d0 + kprime)/(1 + Kprime + Kpp*Ko)
			bb     = (kprime/ko) - Kpp/(1.0d0 + Kprime)
			cc     = (1.0d0 + Kprime + Kpp*Ko)/(Kprime*(1.0d0 + Kprime) - Kpp*Ko)
!	check these - should we be using Pb or Pkb???????
			yy = 1.0d0 + bb * (Pb - Pth);
			V_temp =  vPhCoZero(k,j) * ( 1.0 - aa * (1.0-yy**(-cc)) ) 

			oneminuscc = 1.0d0 - cc
			int_VdP_temp  = vPhCoZero(k,j)*(Pb*(1.0d0-aa)+aa*((1.0d0-bb*Pth)**oneminuscc-		&
					(1.0d0+bb*(Pb-Pth))**oneminuscc)/(bb*(cc-1.0d0))) 
!			int_VdP_temp  = vPhCoZero(k,j)*(Pkb*(1-aa)+aa*((1-bb*Pth)**(1-cc)-(1+bb*(Pkb-Pth))**(1-cc))/(bb*(cc-1))) 
			Qosquared   = dsqrt((1.0d0 - TR/Tcrit))
			if(TK.lt.tcatP) then			! in the alpha field
				Qsquared = dsqrt(((TcatP - TK)/Tcrit))
				else		! in beta field
				Qsquared = 0.0d0
				endif

			Sex = smax*(Qosquared-Qsquared)
			Hex = smax*(Tcrit*(Qosquared-  Qosquared*Qosquared**2/3.0d0 + Qsquared*Qsquared**2/3.0d0) - tcatP*Qsquared) 
			Hex = Hex + vmax*Qosquared*Pb
			Vland = Sex*Vmax/SMax
			VatTP(k,j) = V_temp + Vland
			int_VdP(k,j) = int_VdP_temp

			int_Valpha_dP(k,j) = ao*int_VdP_temp

      			H_tot = H_temp + Hex
      			S_tot = S_temp + Sex
      			HatTP(k,j) = H_tot + int_VdP(k,j) - TK * int_Valpha_dP(k,j)
      			SatTP(k,j) = S_tot - int_Valpha_dP(k,j)
			GatTP(k,j) = H_tot - TK*S_tot + int_VdP_temp

		!	Print *, ''
		!	print *, 'G    =', GatTP(k,j),HatTP(k,j),SatTP(k,j),VatTP(k,j)  



		case(2)		! Bragg-Williams (SF) model
			dH = HP11Model(k,j,1)
			dV = HP11Model(k,j,2)
			W  = HP11Model(k,j,3)
			Wv = HP11Model(k,j,4)
			n  = HP11Model(k,j,5)
			sf = HP11Model(k,j,6)
			!Pkb    = Pb/1000.0d0

			ao     = v1PhCo(k,j)
			ko     = v2PhCo(k,j)
			kprime = v3PhCo(k,j)
			theta  = v4PhCo(k,j)
			mu     = theta/TK
			muo    = theta/TR
			psi    = (mu**2)*dexp(mu)/((dexp(mu) - 1.0d0)**2)
			psio   = (muo**2)*dexp(muo)/((dexp(muo) - 1.0d0)**2)
			! ao is assumed constant
			Pth    = (ao*theta*ko/psio)*(1.0d0/(dexp(mu) - 1.0d0) - 1.0d0/(dexp(muo) - 1.0d0))
!		check - Pb or Pkb???
			! Tim Holland's code
			Kpp    = -Kprime/Ko		! K"
			aa  = (1.0d0 + Kprime) / (1.0d0 + Kprime + Kpp * Ko)
			bb  = Kprime / Ko - Kpp/(1.0d0 + Kprime)
			cc  = (1.0d0 + Kprime + Kpp * Ko) / (Kprime*(1.0d0 + Kprime) - Kpp * Ko)

!			aa     = 1.0d0 + kprime
!			bb     = (kprime/ko)*(2.0d0 + kprime)/(1.0d0 + kprime)
!			cc     = 1.0d0/(2.0d0*kprime + kprime**2)


			yy = 1.0d0 + bb * (Pb - Pth);
			V_temp =  vPhCoZero(k,j) * ( 1.0 - aa * (1.0-yy**(-cc)) ) 
			int_VdP_temp  = vPhCoZero(k,j)*(Pb*(1-aa)+aa*((1-bb*Pth)**(1-cc)-(1+bb*(Pb-Pth))**(1-cc))/(bb*(cc-1))) 

		!	print *,phName(k)
		!	print *,dH,dV,W,Wv,n,sf

			dH = dH + Pb*dV
			W = W + Pb*Wv


              
		!	print *,phName(k)
		!	print *,dH,dV,W,Wv,n,sf
		!	print *, 'x, dx',x,dx  
			xHi = .99999
			xLo = 1.0d-7
			x = (xHi+xLo)/2.0d0
!			print *,'dH, W,n,sf',dH,W,n,sf

15			continue
			if (sf.lt.0.0) then
				Y2=   (n/(n+1))*R*TK*(Dlog(n-n*x)  -sf*Dlog(1-x)  -dlog(1+n*x)  +sf*dlog(n+x))  +W*(2*x-1)+dH
				Y2Hi= (n/(n+1))*R*TK*(dlog(n-n*xHi)-sf*dlog(1-xHi)-dlog(1+n*xHi)+sf*dlog(n+xHi))+W*(2*xHi-1)+dH
				Y2Lo= (n/(n+1))*R*TK*(dlog(n-n*xLo)-sf*dlog(1-xLo)-dlog(1+n*xLo)+sf*dlog(n+xLo))+W*(2*xLo-1)+dH
!				dy2 =((n*R*TK)/(n+1))*((-n/(n-n*x))+(SF/(1-X))-(n/(1+n*x))+(SF/(n+x)))+2*W         
				else
				Y2=   sf*(n/(n+1.0d0))*R*TK*(dlog(n -n*x)  +dlog(1.0d0-x)  -dlog(1.0d0+n*x)  -dlog(n+x))  +	&
						W*(2.0d0*x-1.0d0)+dH
				Y2Hi= sf*(n/(n+1.0d0))*R*TK*(dlog(n -n*xHi)+dlog(1.0d0-xHi)-dlog(1.0d0+n*xHi)-dlog(n+xHi))+	&
						W*(2.0d0*xHi-1.0d0)+dH
				Y2Lo= sf*(n/(n+1.0d0))*R*TK*(dlog(n -n*xLo)+dlog(1.0d0-xLo)-dlog(1.0d0+n*xLo)-dlog(n+xLo))+	&
						W*(2.0d0*xLo-1.0d0)+dH
!				dy2 =((sf*n*R*TK)/(n+1))*((-n/(n-n*x))+(-SF/(1-x))-(n/(1+n*x))+(SF/(n+x)))+2*W 
				endif

!			print *,Xlo,X,Xhi
!			print *,Y2lo,Y2,Y2hi

!			Interval halving code
			if(dabs(y2).lt.1.D-5)go to 16
			if(Y2Lo*Y2.lt.0.0d0)then
				xHi = x
				else		! assumes there is a zero root between 0 and 1
				xLo = x
				endif
			x = (xHi+xLo)/2.0d0
!			Pause '   '
			go to 15
16			continue

!		This is code for Newton's method. It doesn't converge to the correct value if the initial guess is off
!			IF (DABS(Y2).GT.0.1d-7) then
!				dx = -(Y2/dy2) 
!				x = x + dx
		!		print *,'y2,dy2',y2,dy2
		!		print *, 'x, dx',x,dx  
!				Goto 15
!				endif

			Q = x
!			print *, 'Q', x
!			PAUSE "PAUSE"



			Hdis = dH + Q * (W - dH) - Q * Q * W
			Vdis = (1 - Q) * dV + Wv * Q * (1 - Q)
			v_temp = V_temp + Vdis
			xa1 = (1 + n * Q) / (n+1);
			xa2 = (1 - Q) / (n+1);
			xs1 = (n - n * Q) / (n+1);
			xs2 = (n + Q) / (n+1);
			if (sf.lt.0.0)then
				Sdis = - R * (xa1 * dlog(xa1) + xs1 * dlog(xs1) - n * sf * (xa2 * dlog(xa2) + xs2 * dlog(xs2)))
				else 
				Sdis = - R*sf*(xa1 * dlog(xa1) + n * xa2 * dlog(xa2) + xs1 * dlog(xs1) + n * xs2 * dlog(xs2))
				endif

	!		print*, 'Hdis=', Hdis
	!		print*, 'Vdis=', Vdis
	!		print*, 'Sdis=', Sdis
			H_tot = H_temp + Hdis
			S_tot = S_temp + Sdis
			VatTP(k,j) = V_temp
			int_VdP(k,j) = int_VdP_temp
			int_Valpha_dP(k,j) = ao*int_VdP_temp
			HatTP(k,j) = H_tot + int_VdP_temp - TK * int_Valpha_dP(k,j)
			SatTP(k,j) = S_tot - int_Valpha_dP(k,j)
			GatTP(k,j) = H_tot - TK*S_tot + int_VdP_temp
	!		print*, 'HatTP=', HatTP(k,j)
	!		print*, 'SatTP=', SatTP(k,j)
	!		print*, 'VatTP=', VatTP(k,j)
	!		print*, 'GatTP=', GatTP(k,j)


		case(4)		! HP11 melt model
			dkdT  = HP11Model(k,j,1)		! derivative of compressibility
			ao     = v1PhCo(k,j)
			ko     = v2PhCo(k,j)
			kprime = v3PhCo(k,j)
			theta  = v4PhCo(k,j)
			!Pkb    = Pb/1000.0d0

			katT   = ko + dkdT*(TK-TR)
			mu     = theta/TK
			muo    = theta/TR
			psi    = (mu**2)*dexp(mu)/((dexp(mu) - 1.0d0)**2)
			psio   = (muo**2)*dexp(muo)/((dexp(muo) - 1.0d0)**2)
			!Pth    = (ao*theta*ko/psio)*(1.0d0/(dexp(mu) - 1.0d0) - 1.0d0/(dexp(muo) - 1.0d0))
			Pth = 0.0d0
!			Theriak code notes 		
! 			Volume at T and P (hopo 11). No thermal Pressure term, PTH=0
!			aa     = 1.0d0 + kprime
!			bb     = (kprime/ko)*(2.0d0 + kprime)/(1.0d0 + kprime)
!			cc     = 1.0d0/(2.0d0*kprime + kprime**2)
			Kpp    = -Kprime/Ko		! K"
			aa  = (1.0d0 + Kprime) / (1.0d0 + Kprime + Kpp * KatT)
			bb  = Kprime / KatT - Kpp/(1.0d0 + Kprime)
			cc  = (1.0d0 + Kprime + Kpp * KatT) / (Kprime*(1.0d0 + Kprime) - Kpp * KatT)
!			yy = 1.0d0 + bb * (Pb - Pth)
			yy = (1.0d0 + bb * Pb)
!			print *,'////////////'
!			print *,tk,pb,pth
!			print *,aa,bb,cc,yy
!			print *,'////////////'
			!	From Theriak code
			! Volume at T and 1 Bar (hopo 98) Thermal expansion
			!       FV1=V0R*(1.0D0+A0*(T-T0)-20.0D0*A0*(SQT-SQT0))
			! this is pure guesswork
			!!      FV1=V0R*(1.0D0+A0*(T-T0))
!			FV1=V0R*DEXP(A0*(T-T0))
			V_at_1bar_T = vPhCoZero(k,j) * DEXP(ao*(TK-TR))
!			pressure correction at T
			v_Temp = V_at_1bar_T * (1.0 - aa * (1.0 - yy**(-cc)) )		
!			V_temp =  vPhCoZero(k,j) * (1.0 - aa * (1.0-yy**(-cc)) ) 

!			VatTP(k,j) = vPhCoZero(k,j)*(1.0d0 - aa*(1.0d0 - (1.0d0 + bb*(Pb - Pth))**(-cc)))
			!print *,'////////////'
			!print *,tk,pb,pth
			!print *,aa,bb,cc,yy
			!print *,vPhCoZero(k,j),v_at_1bar_T,v_Temp
			!print *,'////////////'

			oneminuscc = 1.0d0 - cc
			! T-D code uses v_at_1bar_T rather than vPhCoZero - this makes sense for the pressure derivative
			! This next line is correct, but since Pth = 0 it can be simplified
!			int_VdP_temp  = v_at_1bar_T*(Pb*(1.0d0-aa)+aa*((1.0d0-bb*Pth)**oneminuscc-		&
!					(1.0d0+bb*(Pb-Pth))**oneminuscc)/(bb*(cc-1.0d0))) 

			int_VdP_temp  = v_at_1bar_T*(Pb*(1.0d0-aa) + aa*(1.0d0 - (1.0d0+bb*Pb)**oneminuscc)/(bb*(cc-1.0d0))) 

			alphahp11 = ao*(psi/psio)*(1.0d0/((1.0d0 - bb*Pth)*(aa + (1.0d0-aa)*(1.0d0 - bb*Pth)**cc)))
!			int_Valpha_dP(k,j) = alphaHP11*int_VdP(k,j)		! which one is correct???
			VatTP(k,j) = V_temp
			int_VdP(k,j) = int_VdP_temp
			int_Valpha_dP(k,j) = ao*int_VdP_temp			! which one is correct???
			!     entropy at T and P 
			SatTP(k,j) = S_temp - int_Valpha_dP(k,j)
			!     enthalpy at T and P
			HatTP(k,j) = H_temp + int_VdP_temp - TK * int_Valpha_dP(k,j)
			!     Gibbs free energy at T and P 
			GatTP(k,j) = H_temp - TK*S_temp + int_VdP_temp

		case(5)		! HP11 phase components with no EOS values
				! so VatTP = Vreference etc
			! if we want to add a thermal expansion term, we could incorporate this
			!V_at_1bar_T = vPhCoZero(k,j) * DEXP(ao*(TK-TR))

			V_temp = vPhCoZero(k,j)
			VatTP(k,j) = V_temp
			int_VdP_temp  = V_temp*(Pb - 1.0d0)
			int_VdP(k,j) = int_VdP_temp
			int_Valpha_dP(k,j) = 0.0d0
			!     entropy at T and P  (no pressure correction to entropy because alpha = 0)
			SatTP(k,j) = S_temp
			!     enthalpy at T and P (no -TK*int_Valpha_dP(k,j) term because alpha = 0)
			HatTP(k,j) = H_temp + int_VdP_temp
			!     Gibbs free energy at T and P 
			GatTP(k,j) = H_temp - TK*S_temp + int_VdP_temp

		
		
		
		

		case default
			write(*,*)'*********************************'
			write(*,*)' Code for HP11ModelSwitch is not 0, 1, 2, 3, or 4 ---- fix the Thermo data file'
			write(*,*)'Phase, PhaseComponent,ModelSwitch = ',K,j,HP11ModelSwitch(k,j)
			write(*,*)'*********************************'
			pause 'hit return to abort'
			stop
		end select

	! HP11 DQF calculations
		if(DQFswitch(k,j).eq.1)then
			HatTP(k,j) = HatTP(k,j) + DQFH(k,j) + DQFV(k,j)*Pb  ! --- why did I add in the V term here?
!			HatTP(k,j) = HatTP(k,j) + DQFH(k,j)   ????  This one seemed to always make melt Gdiff < 0 even when melt was stable. It's wrong (don't use)
			SatTP(k,j) = SatTP(k,j) + DQFS(k,j)
			VatTP(k,j) = VatTP(k,j) + DQFV(k,j)
			GatTP(k,j) = GatTP(k,j) + DQFH(k,j) - DQFS(k,j)*TK + DQFV(k,j)*Pb
			endif

	
	case default
		write(*,*)'*********************************'
		write(*,*)' Code for DatasetKey is not 1, 2 or 3 ---- fix the Thermo data file'
		write(*,*)'DatasetKey = ',Datasetkey(k)
		write(*,*)'*********************************'
		pause 'Hit return to abort'
		stop
	end select


      VMOL(K)=VMOL(K)+(xPhCo(k,j)*VatTP(k,j))

10    CONTINUE
      


	if(grainBoundaryProblem.and.K.eq.grainBoundaryPhase)then
		!	Code for Grain Boundaries
		! nphGB is the grain boundary phase
		! Grain boundary enthalpies are modeled as 
		! H = H(se)/M + H(de)*M where (se) is strain energy and (de) is disorder energy
		! this results in a minimum at some value of M (which is proportional to width of grain boundary)
		! The derivative is
		! dH = ( -H(se)/M^2 + H(de) )dM
			! the dertivatives are added into the dM columns in subroutine Compute2.for
			! loop do 1220 in Compute2.for
		! intital values are set in main calling program. Change in subroutine ThermoData
		! HseGB = 1000000.  	! I'm not sure what these numbers should be (or even what the units are)
		! HdeGB = 100000.
		! miminum is at M = sqrt(Gse/Gde)
		! So for these values, min = sqrt(10) = 3.16
		
		K = grainBoundaryPhase
		do j=1,numPhCo(K)
			HatTP(k,j) = HatTP(k,j) + (HseGB/MP0(K) + HdeGB*MP0(K))
			GatTP(k,j) = GatTP(k,j) + (HseGB/MP0(K) + HdeGB*MP0(K))
			end do

		endif


      return
      end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine Gat1T(TK1,k,j,GTEMP,H_temp,S_temp,iflag)
!     This routine computes G at 1 bar and T for CORK fluids
      implicit none
! ****************************************
	include "Assemb.inc"
! ****************************************
!     LOCAL VARIABLES
      integer*4 j,iflag,k
      REAL*8 R,TR,PR,TK1
      REAL*8 H_temp,S_temp,GTEMP
! *******************************************************************
      DATA R/8.3144D0/
      data TR,PR/298.15D0,1.0D0/
! *******************************************************************
      iflag=0           !0 = no problems
!     COMPUTE S and H AT 1 bar and T

!     Temperature correction for enthalpy
         H_temp = hPhCoZero(k,j) &
                     + aPhCoCp(k,j)        *(TK1-TR) &
                     + 2.D0*bPhCoCp(k,j)   *(Dsqrt(TK1)-Dsqrt(TR))&
                     - cPhCoCp(k,j)        *((1.D0/TK) - (1.D0/TR))&
                     - (dPhCoCp(k,j)/2.D0) *(1.D0/(TK1*TK1) - 1.D0/(TR*TR))&
                     + ePhCoCp(k,j)        *DLOG(TK1/TR)&
                     + (fPhCoCp(k,j)/2.0D0)*(TK1*TK1 - TR*TR) &
                     + (gPhCoCp(k,j)/3.D0) *(TK1*TK1*TK1 - TR*TR*TR)
!         Temperature correction for entropy
          S_temp  =  sPhCoZero(k,j)&
                     + aPhCoCp(k,j)       *DLOG(TK1/TR)&
                     - 2.D0*bPhCoCp(k,j)  *(1.D0/Dsqrt(TK1)-1.D0/Dsqrt(TR))&
                     - 0.5D0*cPhCoCp(k,j) *(1.D0/(TK1*TK1) - 1.D0/(TR*TR))&
                     - (dPhCoCp(k,j)/3.D0)*(1.D0/(TK1*TK1*TK1)-1.D0/(TR*TR*TR))&
                     - ePhCoCp(k,j)       *(1.D0/TK1 - 1.D0/TR)&
                     + fPhCoCp(k,j)       *(TK1-TR)&
                     + (gPhCoCp(k,j)/2.D0)*(TK1*TK1 - TR*TR)
	GTEMP = H_temp - TK1*S_temp           
      return
      end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	SUBROUTINE dlnAdX(K,izero)
	implicit none
!	subroutine to compute derivatives of ln(activity)/Xj for independent components  of a phase K
!	this code is extracted from subroutine dudTPX - everything is removed from that except the composition derivatives

!     izero is an error flag:  izero = 0 no error
!                              izero = 1 returns indicates that a mole fraction is <= 0
! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
! ****************************************
! ***********/CPTemp/************************************************
!	REAL*8 TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS
!	REAL*8 TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
!	common /CPTemp/TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS,&
!		     TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
! *******************************************************************
!     LOCAL VARIABLES
	integer*4 i,j,jjj,k,L,izero,iflag
	integer*4 idep,indep,usePlus,useMinus
	REAL*8 R,RTK,TR,PR
	REAL*8 XH2O,XCO2,TEMP
	real*8 Xtemp(PhCoMax),uMin(PhCoMax),uplus(PhCoMax),uminus(PhCoMax),Xinc,xTempPlus(PhCoMax),xTempMinus(PhCoMax)
! *********************************************************************    
!     Kspar disorder terms
!      Real*8 KSPD0,KSPD1,KSPD2,KSPD3,KSPD4,KSPD5,KSPTmin,KSPTmax
!      data KSPD0,KSPD1,KSPD2,KSPD3,KSPD4,KSPD5,KSPTmin,KSPTmax    &
!          /282.98291D0,-4831.37500D0,3620706.00000D0,0.00000,-0.15733D0,0.00003D0,298.15000D0,1436.15000D0/
! *********************************************************************    
      DATA R/8.3144D0/
      data TR,PR/298.15D0,1.0D0/
! *********************************************************************    

!     definitions of some of the variable arrays

!     lnAct(j) = ln(aj)

!	Feb 7, 2010. 
!	we use this same code and storage for partial derivatives except we skip the temperature and pressure derivatives
!	The activity derivatives are exactly the same, so we use the same code
	TK=TC+273.15D0
	RTK = R*TK
!     Check for temperature between reasonable bounds
	if(TK.LT.200.or.TK.GT.2000)then
		 WRITE(*,*)' *****************************************'
		 WRITE(*,*)' Temperature is out of bounds.  T(K) =',TK
		 WRITE(*,*)' Temperature limits are 200-2000(K)'
		 WRITE(*,*)' Check made in subroutine dlnAdX'
		 izero=20	! temperature out of bounds error
		 PAUSE 'Hit return to continue...'
		 return
		 endif

!     Check to see what type of mineral this is
!	These values are stored in array PhaseType(k)
!			1 = 1 component SPaC
!			2 = alpha-beta quartz (SPaC)
!			3 = multicomponent SPaC
!			4 = H2O SPaC (Haar)
!			5 = Kspar (order-disorder)
!			8 = fixed Al2SiO5 minerals (KAS) - should work for all 3 dataset types
!			10 = Graphite saturated fluid (Connolly_and_Cesare, 1993 JMG)
!			101 = 1 component HP98
!			103 = multicomponent HP98 (includes SF input)
!			104 = Fluid HP98 (CorK)
!			105 = Melt HP98
!			201 = 1 component HP11
!			203 = multicomponent HP11 (includes SF input)
!			204 = Fluid HP11 (Pitzer eq)
!			205 = Melt HP11


! ----------------------------------------------
	select case(PhaseType(k))
	case(1,101,201)
!     		1 component phases
	!     	activity
		J=1
        	lnAct(k,j)=0
!     		No composition derivatives for phase with 1 component
		GO TO 9000
! --------------------------------------------------
	case(2)
		!     THIS SECTION FOR alpha-beta QUARTZ from Berman's database
		!	This code should never be executed because quartz has no composition derivatives
		J=1
		!     activity
		lnAct(k,j)=0.0D0
		!     no composition derivatives for a pure phase such as quartz
		GO TO 9000
! ------------------------------------------------------------------
! ------------------------------------------------------------------
	case(3,103,105,203,205,206)
!     		THIS SECTION FOR User input generalized solution models
		do 3913 J=1,numPhCo(K)
		Xtemp(j) = xPhCo(k,j)
		XtempPlus(j) = xPhCo(k,j)
		XtempMinus(j) = xPhCo(k,j)
3913    continue
!     		Activities
		iflag=0
		call Activity(K,uMin,Xtemp,TK,PB,R,iflag)
		if(iflag.eq.1)then
		!	write(*,*) 'Activity out of bounds in first call min = ',k,phName(k),TC,TK,PB
			write(95,*)'Activity out of bounds in first call min = ',k,phName(k),TC,TK,PB
	
			!call DumpEverything()
			
			go to 9900
			endif
		do 3920 J=1,numPhCo(K)
		lnAct(k,j) = uMin(j)/(R*TK)
3920    continue
!     		composition derivative
		Xinc=1.0D-6
!		Xinc=1.0D-12
		idep = 1
		do 3940 indep = 2,numPhCo(K)
		if(abs(Xtemp(indep)).lt.1.0d-4)then
			Xinc = 0.1D0*Xtemp(indep)	! make the increment 10% of the value of the component. This should eliminate overruns.
			else
			Xinc=1.0D-6
			endif
		useMinus = 1
		usePlus = 1
		XtempPlus(indep) = Xtemp(indep) + Xinc
		XtempPlus(idep)  = Xtemp(idep) - Xinc
		iflag=0
		call Activity(K,uplus,XtempPlus,TK,PB,R,iflag)
		if(iflag.eq.1)then
			if(iDoingGrid.eq.0)then
				write(*,*)' Activity out of bounds in plus derivative'
				write(*,*)'Xinc = ',Xinc
				go to 9900  	! only try to keep going if we are doing a pseudosection - otherwise flag and exit
				endif				
			usePlus = 0		!The derivative of the activity failed -- probably too close to a subsystem (some cation = 0)
							! Try using only the negative increment for derivative	
		!	go to 9900
			endif
		XtempPlus(indep) = Xtemp(indep)
		XtempPlus(idep)  = Xtemp(idep)
		XtempMinus(indep) = Xtemp(indep) - Xinc
		XtempMinus(idep)  = Xtemp(idep) + Xinc
		iflag=0
		call Activity(K,uminus,XtempMinus,TK,PB,R,iflag)
		if(iflag.eq.1)then
			if(iDoingGrid.eq.0)then
				write(*,*)' Activity out of bounds in minus derivative'
				write(*,*)'Xinc = ',Xinc
				go to 9900  	! only try to keep going if we are doing a pseudosection - otherwise flag and exit
				endif				
			useMinus = 0		!The derivative of the activity failed -- probably too close to a subsystem (some cation = 0)
		!	go to 9900
			endif
		XtempMinus(indep) = Xtemp(indep)
		XtempMinus(idep)  = Xtemp(idep)
		if(usePlus.eq.0.and.useMinus.eq.0)then
			write(95,*)'Derivative in dlnAdX failed (3900 block) k = ',k,phName(k),TC,Pb
			 write(95,'(6A12,2x)')(phCoName(k,i),i=1,numPhCo(k))
			 write(95,'(6E14.4)')(xPhCo(k,i),i=1,numPhCo(k))
			go to 9900
			endif
		if(usePlus.eq.1.and.useMinus.eq.1)then
			do 3945 j=1,numPhCo(K)
			temp = (uplus(j)-uminus(j))/(2.0D0*Xinc)
			dudTPX(k,j,1+indep) = temp
3945  			continue
			go to 3940
			endif
		if(usePlus.eq.1.and.useMinus.eq.0)then
			do 3946 j=1,numPhCo(K)
			temp = (uplus(j)-umin(j))/Xinc
			dudTPX(k,j,1+indep) = temp
3946  			continue
			go to 3940
			endif
		if(usePlus.eq.0.and.useMinus.eq.1)then
			do 3947 j=1,numPhCo(K)
			temp = (umin(j)-uminus(j))/Xinc
			dudTPX(k,j,1+indep) = temp
3947  			continue
			go to 3940
			endif
3940    	continue
		GO TO 9000
! -----------------------------------------------------
! ------------------------------------------------------------------
	case(4)
! 	water------------------------------------------------
!     	THIS SECTION IS FOR CALCULATING WATER fugacities, activities
!     	AND ENTROPIES
!     		activity
		j = 1
		lnAct(k,j) = 0.0D0
!     		no composition derivatives for a pure phase like pure water
      		GO TO 9000
! -------------------------------------------------------------
! ------------------------------------------------------------------
	case(5)
		!     this section for K-feldspar, which includes disorder terms
		!     activity
		j = 1
		lnAct(k,j)=0.D0
		!      no composition derivatives for a pure phase
		GO TO 9000
! -------------------------------------------------------------
! ------------------------------------------------------------------
	case(8,108,208)
!	Al2SiO5 polymorphs
		!     	activity
		j = 1
		lnAct(k,j)=0.D0
		!     	no composition derivatives for a pure phase
		go to 9000
! ------------------------------------------------------------------
! ------------------------------------------------------------------
	case(10)
! 		Graphite-COH fluid H:O = 2:1 Connolly and Cesare, 1993 JMG ------------------------------------------------------------
!     		activity
		j = 1
		lnAct(k,j) = 0.0D0
		GO TO 9000
! --------------------------------------------------------------
! --------------------------------------------------------------
	case(104)
! 	CORK fluid
		if(numPhCo(K).eq.2)then
			XCO2=xPhCo(k,2)
			XH2O=1.0D0-XCO2
		! 		H2Oswitch = 3
			else        !either pure H2O or pure CO2
			if(UsePhCo(k,1).eq.1)then    			! we are using pure H2O
				XH2O = 1.0d0
				XCO2 = 0.0d0
		! 		H2Oswitch = 1
				jjj = 1
				lnAct(k,jjj) = 0.0d0
				elseif(UsePhCo(k,2).eq.1)then  		! we are using pure CO2
				XH2O = 0.0d0
				XCO2 = 1.0d0
				jjj = 1
				lnAct(k,jjj) = 0.0d0
		! 		H2Oswitch = 2
				endif
			endif
		go to 9000
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
	case(204)
		! 	Pitzer and Sterner 1994 EOS for H2O and CO2 model
		!	Used in Holland and Powell 2011 dataset
		!	Special case = 15
		!	Only tested for H2O so far
		if(numPhCo(K).eq.2)then
			XCO2=xPhCo(k,2)
			XH2O=1.0D0-XCO2
			write(*,*)' Mixed volatile H2O-CO2 not yet implemented --- aborting'
			pause 'aborting program'
			stop

			else        !either pure H2O or pure CO2
			if(UsePhCo(k,1).eq.1)then    			! we are using pure H2O
				XH2O = 1.0d0
				XCO2 = 0.0d0
				j = 1
				lnAct(k,j) = 0.0d0
				elseif(UsePhCo(k,2).eq.1)then  		! we are using pure CO2
				XH2O = 0.0d0
				XCO2 = 1.0d0
				j = 1
				lnAct(k,j) = 0.0d0
				write(*,*)' CO2 not yet implemented --- aborting'
				pause 'aborting program'
				stop
				endif
			endif
		go to 9000
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
	case default
		write(*,*)' The PhaseType(k) was not correct (no code found). In routine dlnAdX'
		write(*,*)' In subroutine dlnAdX '
		write(*,*)'K,PhaseName(k),PhaseType(k) = ',k,PhName(k),PhaseType(k)
		pause 'Hit return to continue'
	end select
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------


9000  CONTINUE
!     Now the thermodynamic derivatives for all phases are computed
      if(iLong(3).eq.1) then
      write(12,*)' '
      write(12,*) 'Phase component partial derivatives - in dlnAdX'
      write(12,*)'                      dudT          dudP          dudX2         dudX3         dudX4...'
       DO 9010 j=1,numPhCo(K)
      write(12,9020)j,phCoName(k,j),(dudTPX(k,j,L),L=1,numPhCo(K)+1)
9010    continue
9020    FORMAT(i4,1x,A8,4x,15e14.4)
       endif
      RETURN
! -----------------
!     error trap when components of phases have negative values
9900  continue
	if(iLong(1).eq.1)then
	         WRITE(*,*)' *****************************************'
	         WRITE(*,*)' A phase component has a zero or negative value'
	         WRITE(*,*)k,phname(k)
	         write(*,'(6A12,2x)')(phCoName(k,i),i=1,numPhCo(k))
	         write(*,'(6E14.4)')(xPhCo(k,i),i=1,numPhCo(k))
		 write(*,*)'Tc,Pb= ',TC,PB
	         PAUSE 'Hit return to continue...'
	endif
         izero=21
         return
      END




! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine Labqtz(AorB,TK,Pbars,Hlambda,Slambda,Glambda)

!     Computes the change in thermodynamic properties of alpha
!     quartz owing to the lambda transition
!     Follows discussion in Berman (1988, J. Pet)
!	called from case(2) in dudTPX
      implicit none
      integer*4 AorB
      real*8 Tref,TL1bar,DTDPslope,L1,L2,TLatP,Pbars,TK,T9
      real*8 x1,x2,x3,x4,Td,Tr,Hlambda,Slambda,Glambda
      data Tref,TL1bar,DTDPslope,L1,L2/373.0D0,848.0D0,0.023743D0,-0.09186959,0.00024607D0/

      TLatP = TL1bar + DTDPslope*(Pbars-1.0D0)
      Td = TL1bar - TLatP
      Tr = Tref - Td

      x1 = L1*L1*Td + 2.0D0*L1*L2*Td*Td + L2*L2*Td*Td*Td
      x2 = L1*L1 + 4.0D0*L1*L2*Td + 3.0D0*L2*L2*Td*Td
      x3 = 2.0D0*L1*L2 + 3.0D0*L2*L2*Td
      x4 = L2*L2

!     Only integrate from Tr up to the temperature of the transition
      if(TK.gt.TLatP)then
            T9 = TLatP
            AorB=2      ! in beta field
            else
            T9 = TK
            AorB=1      ! in alpha field
            endif
      Hlambda = x1*(T9-Tr) + (x2/2.0D0)*(T9*T9 - Tr*Tr) 	&
                           + (x3/3.0d0)*(T9*T9*T9 - Tr*Tr*Tr) 	&
                           + (x4/4.0d0)*(T9*T9*T9*T9 - Tr*Tr*Tr*Tr)
      Slambda = x1*(Dlog(T9/Tr)) 				&
                     + (x2)*(T9 - Tr) 				&
                     + (x3/2.0d0)*(T9*T9 - Tr*Tr) 		&
                     + (x4/3.0d0)*(T9*T9*T9 - Tr*Tr*Tr)
      Glambda = Hlambda - T9*Slambda
      return
      end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine SmolePhases
!     Calculate molar entropy of all phases in assemblage at T,P and X
! 	Also calculates the G of each phase
      implicit none
! ****************************************
	include "Assemb.inc"
! ****************************************
      	integer*4 k,j,kCur
	real*8 R
      	DATA R/8.3144D0/

	Do 100 kCur = 1,numPh
	k = asmCurrent(kCur)
	SMOL(K) = 0.0d0
	gPhase(K) = 0.0d0
	do 200 j = 1,numPhCo(K)
	SMOL(K) = SMOL(K) + dudtpx(k,j,1)*xPhCo(k,j)
	gPhase(k) = gPhase(k) + xPhCo(k,j)*(gatTP(k,j) + R*TK*lnAct(k,j))
200	continue
	SMOL(K) = -SMOL(K)
100	continue
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      	Subroutine SetFixedMineral(k,phase)
	implicit none
!	Routine to set the values for "fixed" minerals into the appropriate array for calculation in GHSV
!	Called from case(8) in dudTPX
! ****************************************
	include "Assemb.inc"
! ****************************************
      	character*(*) phase
	integer L,k,i
	if(trim(phase).eq.'kyanite')then
		L = 3		! kyanite
		elseif(trim(phase).eq.'sillimanite')then
		L = 4		! sillimanite
		elseif(trim(phase).eq.'andalusite')then
		L = 5		! andalusite
		else
		write(12,*)' In SetFixedMineral did not find phase'
		write(12,*)phase
		pause ' Hit return to continue'
		endif
      	hPhCoZero(k,1)= HZEROF(L)
      	sPhCoZero(k,1)= SZEROF(L)
      	vPhCoZero(k,1)= VZEROF(L)
      	aPhCoCp(k,1)  = ACPF(L)
      	bPhCoCp(k,1)  = BCPF(L)
      	cPhCoCp(k,1)  = CCPF(L)
      	dPhCoCp(k,1)  = DCPF(L)
      	ePhCoCp(k,1)  = ECPF(L)
      	fPhCoCp(k,1)  = FCPF(L)
      	gPhCoCp(k,1)  = GCPF(L)
      	v1PhCo(k,1)   = v1F(L)
      	v2PhCo(k,1)   = v2F(L)
      	v3PhCo(k,1)   = v3F(L)
      	v4PhCo(k,1)   = v4F(L)

	LandauVmax(k,1) = LandauVmaxF(L)
	LandauSmax(k,1) = LandauSmaxF(L)
	LandauTc(k,1)   = LandauTcF(L)

	HP11ModelSwitch(k,1) = HP11ModelSwitchF(L)
	do 10 i = 1,6
	HP11Model(k,1,i) = HP11ModelF(L,i)
10	continue

	DQFswitch(k,1) = DQFswitchF(L)
	DQFH(k,1) = DQFHF(L)
	DQFS(k,1) = DQFSF(L)
	DQFV(k,1) = DQFVF(L)
	
	return
	end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE uZeroAtTP(k,ivol,izero)
      implicit none
!     SUBROUTINE TO COMPUTE G, H, S, V ONLY, FOR PHASE components
!	The same code as in duTPX except all calls to activity are removed
!     ivol is a switch:  ivol = 0 compute everything (v, s, dudTPX)
!                        ivol = 1 compute volumes only (this avoids problems with phase components with negative mole fractions)

!     izero is an error flag:  izero = 0 no error
!                              izero = 1 returns indicates that a mole fraction is <= 0

!     Heat Capacity expression is that of Berman and Brown:
!     CP = A + B*T**-0.5 + C*T**-2 + D*T**-3
!            + E*T**-1   + F*T     + G*T**2
!     UBC data base uses first 4 terms (A,B,C,and D)
!     Maier-Kelly equation uses A, F and C terms Cp= A + F*T + C*T**-2
!     G term is rarely used
!     HP Cp expression used a subset

!     INTEGRATED form :  integral (Cp/T)dT:
!     S(TK,1) = Szero + 
!             + a       *ln(TK/Tr)
!             - B*2*    (1/TK**0.5 - 1/Tr**0.5)
!             - (C/2)*  (1./(TK*TK) - 1./(TR*TR))
!             - (D/3.)* (1./(TK*TK*TK)-1./(TR*TR*TR))
!             - E*      (1/TK - 1/Tr)
!             + F*      (TK-TR)
!             + (G/2)*  (TK*TK - Tr*Tr)

!     INTEGRATED form :  integral CpdT:
!     H(TK,1)) = hPhCoZero + 
!             + a*      (TK-TR)
!             + B*2*    (TK**0.5 - Tr**0.5)
!             - C*      (1/TK - 1/Tr)
!             - (D/2)*  (1./(TK*TK) - 1./(TR*TR))
!             + e*      ln(TK/Tr)
!             + (F/2)*  (TK*TK - Tr*Tr)
!             + (G/3)*  (TK*TK*TK - Tr*Tr*Tr)
!                     
!      Volume at T and P is computed from equation of Berman (1988)
!     Volume terms in Gibbs are       v1, v2, v3 and v4 
!     These correspond to
!     in Berman's program             V1, V2, V3 and V4
!     in Berman (1988) (the paper)
!       the correspondence is         V3, V4, V1 and V2
!   VOLUME EQUATION: (Bermans program)                                                                                  
!    V(P,T)/V(1,298) = 1 + v1PhCo(T-298) + v2PhCo(T-298)**2 + v3PhCo(P-1) + v4PhCo(P-1)**2                           

!         vattp(TK,P) = vzero * (1.D0 + v1 *(TK-298.15D0) +
!     &       v2 *(TK-298.15D0)**2 + v3 *(PB-1.0D0) +
!     &       v4 *(PB-1.0D0)**2)

!     The integral of VdP (dG = VdP) is 
!      int_VdP = vzero ((v3/2.0) - v4)*(PB**2 - PR**2) + (v4/3.0)*(PB**3 - 1.**3)
!             + ((1. - v3 + v4 + v1PhCo(TK - Tr) + v2*(TK - Tr)**2) * (PB - 1.0))

!     to correct enthalpy and entropy up to pressure (at TK) we need the
!     integral of (V*alpha)dP     (alpha is coeff of thermal expansion)
!     int_Valpha_dP = V*(alpha)dP
!     int_Valpha_dP = Vzero [v1 + 2*v2*(TK - TR)]*(PB - 1.)


!     Pressure correction to the entropy 
!      dS = - V*(alpha)dP = -int_Valpha_dP

!     Pressure correction for enthalpy
!     dH = VdP - V*(alpha)*T dP = int_VdP - TK * int_Valpha_dP

!+++++++++++++++++++++++++++++++++++++++++++
!     definitions of some of the variable arrays

!     lnAct(j) = ln(aj)

!     This subroutine underwent major modifications in January 1993 by F. Spear to
!     incorporate a new way of storing data that should make the bookkeeping easier
!     It also will permit multisite solutions to be modeled

!     The basic idea is this.  The governing thermodynamic equations are contained in
!     a matrix AA in subroutine COMPUTE.  Each cell in AA contains the sum of the
!     partial derivative of the chemical potentials of the phase components multiplied
!     by the stoichiometric coefficient of the phase component in the reaction.
!     
!     The purpose of the present subroutine is to generate the matrix for the
!     partial derivatives for each chemical potential with respect to T, P and each
!     independent mole fraction.  The resulting matrix has dimensions
!     npct x nvar (number of phase components x number of independent variables) 
!     
!     The partial derivatives are:

!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.

!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.

!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X

!     Note that the partial derivatives of a component with respect to the mole fractions of components in
!     other phases is zero so the matrix is sparsly populated.

!     The exact form of the partial derivatives depends on the type of solution model.  

!     For crystals, we follow the Wood and Nicholls formulation with Berman's generalized
!     margules model for each site.  
!	All partial derivatives in Gibbs3 are calculated using finite difference approximations
!	   for dµ/dT, dµ/dP and dµ/dX
!     Variable names:
!     dudTPX(25,15) contains the partial derivatives for each chemical potential (first subscript)
!                 with respect to each independent variable (second subscript).


! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
! ****************************************
! ***********/CPTemp/************************************************
      REAL*8 TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS
      REAL*8 TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
      common /CPTemp/TACPS,TBCPS,TCCPS,TDCPS,TECPS,TFCPS,TGCPS,&
     &               TACPH,TBCPH,TCCPH,TDCPH,TECPH,TFCPH,TGCPH
! *******************************************************************

!     LOCAL VARIABLES
	integer*4 i,j,k,L,izero,iflag,ivol
	REAL*8 R,RTK,TR,PR,TT,DISORDS,DISORDH,ppp,gKyanite,gSillimanite,gAndalusite,	&
		GH2O,HH2O,SH2O,VH2O,gplus,gminus, 					&
		Tplus,Tminus,Pplus,Pminus
	real*8 temp,xH2O,xCO2
	real*8 Xtemp(PhCoMax),uplus(PhCoMax),uminus(PhCoMax)


! *******************************************************************
! 	HP98 CORK fluid variables
!      integer*4 H2Oswitch	! for HoPo fluid 1=pure H2O; 2=pure CO2; 3=H2O-CO2 mix
!	real*8 fh2opure,fco2pure
	real*8 fco2plus,fco2minus,fh2oplus,fh2ominus,fh2o,fco2,VH2Oplus,VH2Ominus

! *********************************************************************    
!     ab Quartz variables
      integer*4 AorB
      Real*8 GLambda,Hlambda,SLambda,GLambda_plus,GLambda_minus,VLambda,	&
      	HLambda_plus,HLambda_minus,SLambda_plus,SLambda_minus

!     Kspar disorder terms
      Real*8 KSPD0,KSPD1,KSPD2,KSPD3,KSPD4,KSPD5,KSPTmin,KSPTmax
      data KSPD0,KSPD1,KSPD2,KSPD3,KSPD4,KSPD5,KSPTmin,KSPTmax			&
       	/282.98291D0,-4831.37500D0,3620706.00000D0,0.00000,-0.15733D0,		&
      	0.00003D0,298.15000D0,1436.15000D0/
! *********************************************************************    
      DATA R/8.3144D0/
      data TR,PR/298.15D0,1.0D0/
! *********************************************************************    

      TK=TC+273.15D0
      RTK = R*TK
      
!     Check for temperature between reasonable bounds
      if(TK.LT.200.or.TK.GT.2000)then
         WRITE(*,*)' *****************************************'
         WRITE(*,*)' Temperature is out of bounds.  T(K) =',TK
         WRITE(*,*)' Temperature limits are 200-2000(K)'
         WRITE(*,*)' Check made in subroutine dudTPX'

         izero=1
         PAUSE 'Hit return to continue...'
         return
         endif

!     Check to see what type of mineral this is

!	These values are stored in array PhaseType(k)
!			1 = 1 component SPaC
!			2 = alpha-beta quartz (SPaC)
!			3 = multicomponent SPaC
!			4 = H2O SPaC (Haar)
!			5 = Kspar (order-disorder)
!			8 = fixed Al2SiO5 minerals (KAS) - should work for all 3 dataset types
!			10 = Graphite saturated fluid (Connolly_and_Cesare, 1993 JMG)
!			101 = 1 component HP98
!			103 = multicomponent HP98 (includes SF input)
!			104 = Fluid HP98 (CorK)
!			105 = Melt HP98
!			201 = 1 component HP11
!			203 = multicomponent HP11 (includes SF input)
!			204 = Fluid HP11 (Pitzer eq)
!			205 = Melt HP11


	select case(PhaseType(k))
	case(1,101,201)
! ----------------------------------------------
!     		1 component phases
		call GHSV(K,ivol,iflag)
		if(iflag.eq.1) go to 9900      
		if(ivol.eq.1)  go to 9000
		GO TO 9000
! --------------------------------------------------
	Case(2)
!     		THIS SECTION FOR alpha-beta QUARTZ from Berman's database
      		J=1
!     		call subroutine LabQtz to get extra H, S and G from lambda transition
!     		Parameter AorB = 1 if in alpha field
!               AorB = 2 if in beta field
!     		Calculate implicit B=Vlambda from finite difference approximation
      		call Labqtz(AorB,TK,PB+1.0D0,Hlambda_plus,Slambda_plus,Glambda_plus)
      		call Labqtz(AorB,TK,PB-1.0D0,Hlambda_minus,Slambda_minus,Glambda_minus)
      		Vlambda =    ((Hlambda_plus - Hlambda_minus) - TK*(SLambda_plus - SLambda_minus))/2.0d0
!     		now calculate values at the real T and P
      		call Labqtz(AorB,TK,PB,Hlambda,Slambda,Glambda)
!     		This code is executed only for the alpha-beta quartz mineral
!     		It sets the values of the thermodynamic properties based on the type of mineral
      		L = 1       ! default is alpha field
!      		IF(AorB.eq.2) L = 2     !in Beta field (According to R Berman, we use alpha values throughout
!     		Set reference thermodynamic properties to appropriate values
		hPhCoZero(k,j)= HZEROF(L)
		sPhCoZero(k,j)= SZEROF(L)
		vPhCoZero(k,j)= VZEROF(L)
		aPhCoCp(k,j)  = ACPF(L)
		bPhCoCp(k,j)  = BCPF(L)
		cPhCoCp(k,j)  = CCPF(L)
		dPhCoCp(k,j)  = DCPF(L)
		ePhCoCp(k,j)  = ECPF(L)
		fPhCoCp(k,j)  = FCPF(L)
		gPhCoCp(k,j)  = GCPF(L)
		v1PhCo(k,j)   = v1F(L)
		v2PhCo(k,j)   = v2F(L)
		v3PhCo(k,j)   = v3F(L)
		v4PhCo(k,j)   = v4F(L)

		call GHSV(K,ivol,iflag)
		if(iflag.eq.1) go to 9900      
		if(ivol.eq.1)  go to 9000
!     		activity
      		lnAct(k,j)=0.0D0
!     		correct for lambda transition
		VATTP(k,j) = VATTP(k,j) + VLambda
		SATTP(k,j) = SATTP(k,j) + SLambda
		HatTP(k,j) = HatTP(k,j) + HLambda
		GatTP(k,j) = HatTP(k,j) - TK*SatTP(k,j)
	!if(k.eq.1)then
	!write(*,*)'Quartz in uZeroatTP ',GatTP(k,1)
	!pause 'hit return to continue'
	!endif
		GO TO 9000
! -----------------------------------------------------
	case (3,103,105,203,205,206)
!		Multicomponent, multisite phases
	!     THIS SECTION FOR User input generalized solution models
		call GHSV(K,ivol,iflag)
		if(iflag.eq.1)then
			write(*,*)' Error in routine GHSV -- aborting calculations'
			izero = 1
			return
			endif      
		if(ivol.eq.1)go to 9000
		GO TO 9000
! ------------------------------------------------------------------

! -------------------water------------------------------------------------
!     THIS SECTION IS FOR CALCULATING WATER fugacities, activities
!     AND ENTROPIES - HAAR equation
	case(4)
		j = 1
		PPP=PB
		if(PFluidSwitch.ne.0)PPP=PFLUID
	!     	THIS FOR PURE WATER
		call WHAAR2 (TK,PPP,GH2O,HH2O,SH2O,VH2O,izero)
		if(izero.eq.1)return
		VATTP(k,j) = VH2O
		SATTP(k,j) = SH2O
		HATTP(k,j) = HH2O
		GATTP(k,j) = GH2O
	!     	VOLUME IS RETURNED IN joule/BAR
		VMOL(K)=VH2O
		GO TO 9000
! ------------------------------------------------------------------
	case(5)
!     		this section for K-feldspar, which includes disorder terms
		J=1
		L=1
		call GHSV(K,ivol,iflag)
		if(iflag.eq.1) go to 9900      
		if(ivol.eq.1)  go to 9000
!     		UBC DISORDER equation
!      		Cp(disorder) = D0 + D1/SQRT(T) + D2/T/T +            D3/T + D4*T + D5*T*T
!      		Cp =           K0 + K1/SQRT(T) + K2/T/T + K3/T/T/T + K4/T + K5*T + K6*T*T                                
      		TT=TK
      		IF(TK.GT.KSPTmax) TT = KSPTmax
!     		KSPD0,KSPD1,KSPD2,KSPD3,KSPD4,KSPD5,KSPTmin,KSPTmax
		DISORDH =       KspD0*(TT-TR) 			&
		         + 2.D0*KspD1*(Dsqrt(TT) - Dsqrt(TR))   &
			      - KspD2*((1.D0/TT) - (1.D0/TR))   &
			      + KspD3*DLOG(TT/TR) 		&
			     + (KspD4/2.0D0)*(TT*TT - TR*TR) 	&
			     + (KspD5/3.D0) *(TT*TT*TT - TR*TR*TR)
	      	Hattp(k,j) = Hattp(k,j) + disordH
      		DISORDS       =   KspD0*DLOG(TT/TR) 			&
                          - 2.D0*(KspD1)*(1.D0/sqrt(TT)-1.D0/sqrt(TR)) 	&
                         - 0.5D0*(KspD2)*(1.D0/(TT*TT) - 1.D0/(TR*TR)) 	&
                               - (KspD3/3.D0)*(1.D0/TT - 1.D0/TR) 	&
                               + (KspD4)* (TT-TR) 			&
                               + (KspD5/2.0D0)*  (TT*TT - TR*TR)
		sattp(k,j)= sattp(k,j) + disordS
		GatTP(k,j) = GatTP(k,j) + DISORDH - TK * disords
      		GO TO 9000
! -------------------------------------------------------------
	case(8,108,208)
	!     	This section for kyanite-sillimanite-andalusite polymorphs
	! 	Find the Al2SiO5 polymorph with the lowest G at the given T and P
	! 	Then use that phase in all the calculations
	! 	Kyanite
		call SetFixedMineral(k,'kyanite')
		call GHSV(K,ivol,iflag)
		gKyanite = gattp(k,1)
		call SetFixedMineral(k,'sillimanite')
		call GHSV(K,ivol,iflag)
		gSillimanite = gattp(k,1)
		call SetFixedMineral(k,'andalusite')
		call GHSV(K,ivol,iflag)
		gAndalusite = gattp(k,1)
	! 	Now pick phase with lowest G
		if(gAndalusite.lt.gKyanite.and.gAndalusite.lt.gSillimanite)then
			! andalusite has lowest G
			call SetFixedMineral(k,'andalusite')
	! 		write(12,*)'Stable phase is andalusite'
		else if(gkyanite.lt.gAndalusite.and.gKyanite.lt.gSillimanite)then
			! kyanite has lowest G
			call SetFixedMineral(k,'kyanite')
	! 		write(12,*)'Stable phase is kyanite'
		else
			! sillimanite has lowest G
			call SetFixedMineral(k,'sillimanite')
	! 		write(12,*)'Stable phase sillimanite'
		endif
		call GHSV(K,ivol,iflag)
		go to 9000
! ------------------------------------------------------------------
! ------------------------------------------------------------------
! ------Graphite-COH fluid H:O = 2:1 Connolly and Cesare, 1993 JMG ------------------------------------------------------------
	case(10)
		j = 1
		call graphite_sat(TK,pb,Gh2o)
		! 	finite difference entropy
		call graphite_sat(TK+.1,pb,Gplus)
		call graphite_sat(TK-.1,pb,Gminus)
		sh2o = -(gplus - gminus)/0.2d0	
		! 	finite difference volume
		call graphite_sat(tk,pb+1,gplus)
		call graphite_sat(tk,pb-1,gminus)
		vh2o = (gplus - gminus)/2.d0	
		hh2o = gh2o + tk*sh2o	
		! 	write(*,*)'G,H,S,V ',Gh2o,Hh2o,Sh2o,Vh2o
		GatTP(k,j) = Gh2o
		SATTP(k,j) = Sh2o
		VATTP(k,j) = vh2o
		VMOL(K)    = vh2o
		HATTP(k,j) = hh2o
		GO TO 9000
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
	case(104)
! 	CORK fluid
		j = 1
		PPP=PB
		if(PFluidSwitch.ne.0)PPP=PFLUID
		! 	we calculate volume as dG/dP, because CORK returns only fugacity
		! 	so we set ivol = 0 here to force GHSV to return G at (1,T).
		call GHSV(K,0,iflag)	
		! 	call GHSV(K,ivol,iflag)
		if(iflag.eq.1)return      
		if(numPhCo(K).eq.2)then
			XCO2=xPhCo(k,2)
			XH2O=1.0D0-XCO2
	! 		H2Oswitch = 3
			else        !either pure H2O or pure CO2
			if(UsePhCo(k,1).eq.1)then    			! we are using pure H2O
				XH2O = 1.0d0
				XCO2 = 0.0d0
		! 		H2Oswitch = 1
				j = 1
				call Gat1T(TK,k,j,GH2O,HH2O,SH2O,iflag)
				call hprk(ppp,tk,xco2,fh2o,fco2)	! returns ln fugacity
				GatTP(k,j) = GH2O + RTK*fh2o
				Tplus = TK+0.1d0
				call Gat1T(Tplus,k,j,gplus,HH2O,SH2O,iflag)
				call hprk(ppp,Tplus,xco2,fh2oplus,fco2plus)
				gplus = gplus + R*Tplus*fh2oplus
				Tminus = TK - 0.1d0
				call Gat1T(Tminus,k,j,gminus,HH2O,SH2O,iflag)
				call hprk(ppp,Tminus,xco2,fh2ominus,fco2minus)
				gminus = gminus + R*Tminus*fh2ominus
				SatTP(k,j) = -(gplus-gminus)/0.2d0
				HatTP(k,j) = GatTP(k,j) + TK*SatTP(k,j)
				Pplus = ppp + 1.0d0
				call hprk(Pplus,tk,xco2,fh2oplus,fco2plus)
				gplus = GatTP(k,j) + RTK*fh2oplus
				Pminus = ppp-1.0d0
				call hprk(Pminus,tk,xco2,fh2ominus,fco2minus)
				gminus = GatTP(k,j) + RTK*fh2ominus
				VatTP(k,j) = (gplus-gminus)/2.0d0
				VMOL(K)    = VatTP(k,j)
				elseif(UsePhCo(k,2).eq.1)then  		! we are using pure CO2
				XH2O = 0.0d0
				XCO2 = 1.0d0
		! 		H2Oswitch = 2
				endif
			endif
		go to 9000
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
	case(204)
! 	Pitzer and Sterner 1994 EOS for H2O and CO2 model
!	Used in Holland and Powell 2011 dataset
!	Special case = 15
!	Only tested for H2O so far
		j = 1
		PPP=PB
		if(PFluidSwitch.ne.0)PPP=PFLUID
	! 	we calculate volume as dG/dP, because CORK returns only fugacity
	! 	so we set ivol = 0 here to force GHSV to return G at (1,T).
		call GHSV(K,0,iflag)	
	! 	call GHSV(K,ivol,iflag)
		if(iflag.eq.1)return      
		if(numPhCo(K).eq.2)then
			XCO2=xPhCo(k,2)
			XH2O=1.0D0-XCO2
			write(*,*)' Mixed volatile H2O-CO2 not yet implemented --- aborting'
			pause 'aborting program'
			stop

			else        !either pure H2O or pure CO2
			if(UsePhCo(k,1).eq.1)then    			! we are using pure H2O
				XH2O = 1.0d0
				XCO2 = 0.0d0
		! 		H2Oswitch = 1
				j = 1
				if(Pb.lt.100)then
					write(*,*)' Pressure is too low ... abort. Pb = ',Pb
					izero = 1
					go to 9000
					endif
				call PS94H2O(Pb,TK,GH2O,VH2O)
				GatTP(k,j) = GH2O		! Gibbs uses J (not kJ)
				VatTP(k,j) = VH2O
	!			Temperature derivative (to get SH2O)
				Tplus = TK+0.1d0
				call PS94H2O(Pb,Tplus,Gplus,VH2Oplus)
				Tminus = TK - 0.1d0
				call PS94H2O(Pb,Tminus,Gminus,VH2Ominus)
				SatTP(k,j) = -(gplus-gminus)/0.2d0
				HatTP(k,j) = GatTP(k,j) + TK*SatTP(k,j)
				VMOL(K)    = VatTP(k,j)

				elseif(UsePhCo(k,2).eq.1)then  		! we are using pure CO2
				XH2O = 0.0d0
				XCO2 = 1.0d0
				write(*,*)' CO2 not yet implemented --- aborting'
				pause 'aborting program'
				stop
				endif
			endif
		go to 9000

	case default
		write(*,*)' The PhaseType(k) was not correct (no code found). In routine uZeroAtTP'
		write(*,*)' In subroutine uZeroAtTP line 2130   '
		write(*,*)'K,PhaseName(k),PhaseType(k) = ',k,PhName(k),PhaseType(k)
		pause 'Hit return to continue'
	end select
! -----------------------------------------------------------------------
! -----------------------------------------------------------------------
9000  CONTINUE
!     Now the thermodynamic derivatives for all phases are computed
      if(iLong(3).eq.1) then
	      write(12,*)' '
	      write(12,*) 'Phase component partial derivatives - in duTPX'
	      write(12,*)'                      dudT          dudP          dudX2         dudX3         dudX4...'
!		Do 9010 kCur = 1,numPh
!		k = asmCurrent(kCur)
       		DO 9010 j=1,numPhCo(K)
      		write(12,9020)j,phCoName(k,j),(dudTPX(k,j,L),L=1,numPhCo(K)+1)
9010    	continue
9020    	FORMAT(i4,1x,A8,4x,15e14.4)
       		endif
      RETURN
! -----------------
!     error trap when components of phases have negative values
9900  continue
	if(iLong(1).eq.1)then

	         WRITE(*,*)' *****************************************'
	         WRITE(*,*)' A phase component has a zero or negative value'
	         WRITE(*,*)k,phname(k)
	         write(*,'(6A12,2x)')(phCoName(k,i),i=1,numPhCo(k))
	         write(*,'(6E14.4)')(xPhCo(k,i),i=1,numPhCo(k))
		 write(*,*)'Tc,Pb= ',TC,PB
	         PAUSE 'Hit return to continue...'
	endif
         izero=1
         return
      END
