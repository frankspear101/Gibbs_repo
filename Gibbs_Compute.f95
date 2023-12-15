! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE C2OMPUTE(izero)
	use MatrixArrays
      implicit none

! c     version modified June 17, 1988
! C     ROUTINE TO COMPUTE SLOPES AND DELTAS
! C
! C
! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "Singular.inc"
!	include "Solute.inc"
	include "Tangent.inc"
	include "MatrixSingular.inc"
! c*****************************************
!     local variables
      integer*4 i,ii,ii_1,j,jj,k,L,LL,Lcont,jcol,irow,ier,j1,kCur,&
     & izero,imiss,nsolve,newtoncounter,ivol,ALLXrow,iok
      REAL*8 DETA,R,SUM,TEMP,Pi,TEMPT,TEMPP,TEMPX,TEMPH,TEMPS,TEMPK,TempPFluid,xTemp(20),allXTemp(100)
      CHARACTER vn51(100)*8
      DATA R,Pi/8.3144D0,3.141592654D0/


! --------------------------------------------

      IF(iLong(3).eq.1)then

         write(12,*)' '
         write(12,*)'------------------------------------------------------'
         write(12,*)'------------------------------------------------------'
         write(12,*)' Variable list '
         do 11 i = 1,nvar
         write(12,15)I,VN1(I),VN2(I)
15       FORMAT(I3,2X,A2,A16)
11       continue
         write(12,*)'NEQ=',NEQ,'  NVAR=',NVAR
         write(12,*)
         write(12,*)'VARIANCE =',NVAR-NEQ
         write(12,*)
         write(12,*)'MONITOR PARAMETERS    Delta X for monitors'
         if(nvar-neq.gt.0)then
         do 12 i = 1,nvar-neq
         write(12,16)MON(I),VN1(MON(I)),VN2(MON(I)),Deltax(i)*NSTEP
16       FORMAT(I3,2X,A2,A16,F12.3)
12       continue
!         write(12,*)' DELTA Xi FOR MONITOR PHASES'
!         write(12,16)(DELTAX(I)*NSTEP,I=1,NVAR-NEQ)
!16       FORMAT(10F12.3)
         else
         write(12,*)' No monitor parameters (variance <= 0)'
         endif
         write(12,*)'NUMBER OF FINITE DIFF. ITERATIONS=',NSTEP
         write(12,*)
         write(12,*)'INITIAL MINERAL DATA'
         write(12,17)TC,PB
17       format(' T START (DEG C)=',F10.2,'   P START=',F10.2)
	if(PFluidSwitch.gt.0)write(12,18)PFluid
18	format(' PFluid = ',F10.2)
	endif
90    CONTINUE
!     LCONT IS COUNTER IN CONTOUR BLOCK
      LCONT=0
      newtoncounter=0

! -----------------------------------------------------------
!                             
100   CONTINUE


      TK=TC+273.15D0

!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST

!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      izero=0
      ivol=0            !compute everything
      K=0
! x       write(*,*)' In c2ompute #3-calling duTPX'
! x       pause
!	This call will callculate all thermodynamic properties - the TP dependent ones and the composition dependent ones
!	We need this when we're using T or P as dependent variables. But not when they're independent (because they're already set)
      CALL ALLKduTPX(ivol,izero)
! x       write(*,*)' In c2ompute #4-back from duTPX'
! x        pause
      if(izero.gt.0)then
      		write(*,*)' izero > 0 in sub C2OMPUTE'
      		go to 9999
		endif
!     COMPUTE PHASE VOLUMES
      IF(IMASS.EQ.0)GO TO 342
!         DO 340 k = 1,numPh
	Do 340 kCur = 1,numPh
	k = asmCurrent(kCur)
            VP0(K)=MP0(K)*VMOL(K)*10.D0
340         continue

342   CONTINUE

!     OUTPUT THERMO DATA TO DISK

      IF(iLong(3).eq.1)then
350       CONTINUE
	Do 390 kCur = 1,numPh
	k = asmCurrent(kCur)
          write(12,*)
          write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       FORMAT(' ',I4,F8.1,4X,A32)
          write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          IF(IMASS.EQ.1) then
             write(12,507)MP0(K),VP0(K)
          endif
507       FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          write(12,501)(phCoName(k,j),J=1,numPhCo(K))
501       FORMAT('          ',12(A4,6X))
          write(12,504)(xPhCo(k,j),J=1,numPhCo(K))
504       FORMAT('  COMP ',12(F10.4))
          write(12,505)(Dexp(lnAct(k,j)),J=1,numPhCo(K))
505       FORMAT(' Activ  ',12(E10.3))
          write(12,506)(lnAct(k,j),J=1,numPhCo(K))
506       FORMAT('  Ln(a) ',12(F10.5))
          write(12,514)(hPhCoZero(k,j),J=1,numPhCo(K))
514       Format(' hPhCoZero',12E15.7)
          write(12,515)(HATTP(k,j),J=1,numPhCo(K))
515       Format(' HATTP',12E15.7)
          write(12,509)(sPhCoZero(k,j),J=1,numPhCo(K))
509       Format(' SZERO',12F15.5)
          write(12,510)(SATTP(k,j),J=1,numPhCo(K))
510       Format(' SATTP',12F15.5)
          write(12,508)(vPhCoZero(k,j),J=1,numPhCo(K))
508       Format(' VZERO',12F15.5)
          write(12,513)(VATTP(k,j),J=1,numPhCo(K))
513       Format(' VATTP',12F15.5)
          write(12,516)(GATTP(k,j),J=1,numPhCo(K))
516       Format(' GATTP',12E15.7)
          write(12,*)'dudTPX ARRAY:'
          write(12,*)'     dudT   .    dudP   .    dudX2  .    dudX3  .    dudX4...'
            do 512 j = 1,numPhCo(K)
           	write(12,511)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
511       	FORMAT(15E12.5)
512   	continue
!         WRITE(iout,408)' ACP',(aPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' BCP',(bPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' CCP',(cPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' DCP',(dPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' ECP',(ePhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' FCP',(fPhCoCp(k,j),J=1,numPhCo(K))
!         WRITE(iout,408)' GCP',(gPhCoCp(k,j),J=1,numPhCo(K))
408       format(' ',a5,6E10.3)
551       CONTINUE
390       CONTINUE

	endif
395       CONTINUE


!     THIS PART OF THE PROGRAM SETS UP THE MASTER MATRIX

!     ZERO MASTER ARRAY

      DO 1000 J=1,NVAR
      DO 1000 I=1,NEQ
1000  AA(I,J)=0.0D0

!     SET UP equations for the NRX reactions
!     These equations are generated by multiplication of the matrix ARX(I,J)
!     by the matrix dudTPX:
 
!      AA = ARX(nrx,np) x dudTPX(np,nvar)

!     The partial derivatives in matrix dudTPX are:

!     In storage space 1: dudTPX(jj,1):
!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.

!     In storage space 2: dudTPX(jj,2):
!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.

!     In storage space 3 to (numPhCo(K)+1): dudTPX(jj,3), dudTPX(jj,3), dudTPX(jj,3) etc:
!           Note that there are only numPhCo(K)-1 partial derivatives and the subscripts range
!           from 3 to numPhCo(K)+1.
!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X

!       
!     The matrix ARX(nrx,npc) contains the stoichiometric coefficients for each phase component
!     in each of the reactions (subscript 1 refers to reaction number, subscript
!     2 refers to the phase component)

!     The matrix AA therefore contains in each cell terms such as
!     
!       (n(i)*u(i))/T
!       (n(i)*u(i))/P
!       (n(i)*u(i))/X2
!       (n(i)*u(i))/X3

!     for each independent variable (where n(j,i) are the stoichiometric
!     coefficients for phase component i in reaction j)
!       Each row of AA contains a different set of reaction coefficients.

!     The only catch to the calculation of AA is that the matrix dudTPX is
!     not stored in true matrix form.  That is, to save space the independent variables
!     are listed for each phase only.
!     This simply requires a bit of extra bookkeeping during the matrix multiplication

! c-------------code that works with no fluid pressure option--------
!      ii = np-nrx      !note that reaction coefficients are stored in the bottom of array ARX
!      DO 1010 i=1,NRX
!      TEMPT=0.D0
!      TEMPP=0.D0
!      do 1015 J = 1,NP
!      dT coefficient
!      TEMPT =TEMPT + ARX(i+ii,j)*dudTPX(j,1)	!S reaction
!  !     dP coefficient
!       TEMPP =TEMPP + ARX(i+ii,j)*dudTPX(j,2)	!V reaction
! c1015  continue      
!       AA(i,1)=TEMPT
!       AA(i,2)=TEMPP
! c1010  continue
! c--------------------------------------------

! c	New code to implement fluid pressure as independent variable
	ii = np-nrx      !note that reaction coefficients are stored in the bottom of array ARX
	DO 1010 i=1,NRX
	TEMPT=0.D0
	TEMPP=0.D0
	TEMPPFluid = 0.0D0
	j = 0
	Do 1015 kCur = 1,numPh
	k = asmCurrent(kCur)
	do 1015 jj = 1,numPhCo(k)
	j = j + 1
!     	dT coefficient
      	TEMPT =TEMPT + ARX(i+ii,j)*dudTPX(k,jj,1)				! S reaction
!     	dP coefficient
! 	if a fluid is present and PFluidSwitch>0 we must put PFluid*Vfluid into an independent variable
	if(PFluidSwitch.eq.0)then					! Fluid pressure is dependent
		TEMPP =TEMPP + ARX(i+ii,j)*dudTPX(k,jj,2)			! V reaction
		else							! Fluid pressure is independent
! 		now we must determine if "j" is the index for a fluid component
!		jj = Minpoint(KFlu) - 1
!	I'm not sure this code is correct....
		if(k.eq.kflu)then
			do 1016 LL = 1,numPhCo(Kflu)
!		if(j.eq.jj+LL)then					! This is a fluid component if true
			TEMPPFluid =TEMPPFluid + ARX(i+ii,j)*dudTPX(k,jj,2)	! V reaction
1016			continue
			else						! not a fluid component
			TEMPP =TEMPP + ARX(i+ii,j)*dudTPX(k,jj,2)		! V reaction
			endif
		endif
1015  	continue      
      	AA(i,1)=TEMPT
      	AA(i,2)=TEMPP
      	if(PFluidSwitch.gt.0)AA(i,3) = TEMPPFluid
1010  	continue

!     	dX(j,k) coefficients
	DO 1020 i=1,NRX
	L=TANDP
	jj = 0
!	Do 1021 k = 1,numPh
	Do 1021 kCur = 1,numPh
	k = asmCurrent(kCur)
	if(numPhCo(K).eq.1)then
		jj = jj + 1
		go to 1021	! there are no dX terms for a phase of fixed composition
		endif
	Do 1022 LL = 3,numPhCo(K)+1     	! this loops on the number of independent composition derivatives
	tempX = 0.D0
	do 1023 j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
	TEMPX=TEMPX + ARX(i+ii,jj+j)*dudTPX(k,j,LL)
1023  	continue
      	L=L+1
      	AA(i,L)=TEMPX
1022  	continue
	jj = jj + numPhCo(k)
1021  	continue
1020  	continue

!     	SET UP MASS BALANCE EQUATIONS
      	if(imass.eq.0)go to 1100
!     	FIRST SET UP dM Columns (M=MOLE FRACTION OF PHASE)
!     	JCOL COUNTS WHERE dM TERMS ARE
      	JCOL = TANDP + NX
!     	loop through each phase (there is a dM term for each phase)
	Do 1110 kCur = 1,numPh
	k = asmCurrent(kCur)
      	JCOL = JCOL + 1
!     	LOOP THROUGH ALL mass balance equations for each SYSTEM COMPONENT
      	call bulkcomp				! get total moles of each system component (moles(i) in common block)
      	DO 1105 i = 1,nc
!      	if(noOxygenMassBalance.eq.1)then
!		go to 1105			! for grain boundary problems, skip mass balance for oxygen
!	      endif
!     	IROW counts the row number
      	iRow = i + nRx
      	sum=0.D0
!     	LOOP THROUGH ALL PHASE COMPONENTS for this phase
      	DO 1120 j=1,numPhCo(K)
1120  	sum = sum + comp(k,J,I)*xPhCo(k,J)
      	AA(IROW,JCOL)=SUM
      	YY(IROW) = -(moles(i) - openFlux(i) - molesStart(i))	! molesStart(i) is the bulk composition we want to fit
      								! moles(i) is the current bulk comp based on phase comp and M(k)
								! when they are equal, we have the correct X and M for the desired bulk comp
								! openFlux are open system components either dependent or independent variables
1105  	CONTINUE
1110  	CONTINUE	!loop for each phase

!     	NOW SET UP COEFFICIENTS FOR dX TERMS in mass balance equations
!     	Note that the dependent dX term is removed as above
!     	JCOL is the column for the dX term
      	JCOL=TANDP
!     	FIRST LOOP THROUGH ALL PHASES TO FIND independent dX TERMS
	Do 1130 kCur = 1,numPh
	k = asmCurrent(kCur)
      	IF(numPhCo(K).EQ.1) GO TO 1130
!     	IF HERE THEN WE HAVE FOUND A PHASE WITH dX TERMS (I.E. MORE THAN 1 COMPONENT)
!     	LOOP THROUGH ALL dX TERMS FOR THIS PHASE
      	DO 1140 J=2,numPhCo(K)
      	JCOL=JCOL+1
!     	LOOP THROUGH ALL MASS BALANCE EQUATIONS FOR EACH dX TERM
      	DO 1150 I=1,NC
!      	if(noOxygenMassBalance.eq.1)then
!	      go to 1150	! for grain boundary problems, skip mass balance for oxygen
!	      endif
      	IROW= I + NRX
	AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))
1150	continue
1140  	CONTINUE
1130  	CONTINUE

!	Check to see if one of the phases is a GrainBoundary. If so, then add the dM terms to the energy equations.
	if(grainBoundaryProblem)then	
	! nphGB is the grain boundary phase
	! Grain boundary enthalpies are modeled as 
	! H = H(se)/M + H(de)*M where (se) is strain energy and (de) is disorder energy
	! this results in a minimum at some value of M (which is proportional to width of grain boundary)
	! The derivative is
	! dH = ( -H(se)/M^2 + H(de) )dM
       	ii = np-nrx      			! note that reaction coefficients are stored in the bottom of array ARX
      	JCOL = TANDP + NX + grainBoundaryPhase	! JCOL is the column where M(GB) is
      	DO 1220 i=1,NRX
	K = grainBoundaryPhase		   	! is the grain boundary phase number in the input list
							! note code will need to be written if I change assemblages
      	tempX = 0.D0
      	do 1223 j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
      	TEMPX=TEMPX + ARX(i+ii,jj+j)*(-HseGB/MP0(K)**2 + HdeGB)
1223  	continue
      	AA(i,JCOL)=TEMPX
1220  	continue
	endif

1100  continue   ! end of mass balance equations


!     	Set up open system equations
      	if(iopen.NE.0)then
!           (note if iopen0 then imass must=1)
            jcol = TANDP + NX + numPh
            irow=nrx
            do 1210 i=1,iopen
            AA(irow+iopena(i),jcol+i) = -1.D0
1210        continue
            endif

!     Set up EXTRA equations (there are numNew of them)
      	if(numNew.eq.0)go to 1200
      	JCOL=TANDP + NX
      	IROW=NRX		! iRow is the row at the end of the reactions (+ mas balance, if needed) and the start of the newVar equations
      	IF(IMASS.EQ.1)THEN
!		jCol is the column for the start of the new variables in the A matrix
         	JCOL=JCOL+numPh+iopen
         	IROW=IROW+NC
         	ENDIF
!     	Code for setting up generic new variable
!     	Form of new variable is a linear combination of the phase components
!     	for the phase of interest, as specified in the thermodynamic data file
!     	0 = New - (A + aj*Xj)/(B + bj*Xj)
!     	Differential of this expression is
!     	0 = dNew -  (  aj/(B + bj*Xj) - (A + aj*Xj)*bj/(B + bj*Xj)^2 ) dXj
!     	Q: How do we check to ensure that all the components necessary to form the
!        	new variable are actually in use?   -- This should probably be done when new
!        	minerals are read in.
!     	The New variable definitions are written in terms of independent + dependent
!     	components, so we will proceed with the generalized derivative including the
!     	dependent component and then remove it at the end.
!     	NewNumA(k,i) = A constant
!     	NewNumaj(k,i,j) = aj coeffs
!     	NewDenB(k,i) = B constant
!     	NewDenbj(k,i,j) = bj coeffs
      	ALLXrow = TandP + NX
!      	if(imass.eq.1)ALLXrow = ALLXrow + numPh
	if(imass.eq.1)ALLXrow = ALLXrow + numPh + iopen	! should this be the way for open systems?
!	numNew = 0		! total number of new variables
      	ii = TandP
	do 1305 kCur = 1,numPh
	k = asmCurrent(kCur)
	if(includeNew(k).eq.1)then	
		do 1310 i = 1,numNewInPhK(k)
!		do 1310 LL=1,numNew
!		i = Jnew(LL)
!		K = Knew(i)
! 		compute numerator
		numerator = NewNumA(k,i)
		do 1311 j = 1,numPhCo(K)
1311  		numerator = numerator + newNumAj(k,i,j)*xPhCo(k,j)
		! compute denominator
		denominator = NewDenB(k,i)
		do 1312 j = 1,numPhCo(K)
1312  		denominator = denominator + NewDenBj(k,i,j)*xPhCo(k,j)
      		denominator2 = denominator*denominator          ! denominator squared
      		IROW=IROW+1
!     		calculate the coefficient for the dependent component (component 1)
      		temp = NewNumAj(k,i,1)/denominator - numerator*NewDenBj(k,i,1)/denominator2
!     		now we can loop on only the independent variables
!		This equation is wrong - need to calculate where the values go...
!		ii is the columnwhere the phase component derivatives start for the new variables
!			not all phases have new variables	
!      		ii = TandP + NewPoint(i) - (Knew(i) - 1) - 2
		ii_1 = ii-1
      		DO 1315 j = 2,numPhCo(K)
      		AA(IROW,ii_1 + j) = -((NewNumAj(k,i,j)/denominator - numerator*NewDenBj(k,i,j)/denominator2) - temp)   ! temp is the dependent component coefficient
1315  		continue
!		jCol is the column for this new variable
		JCOL = JCOL + 1			! this should be correct
      		AA(IROW,JCOL) = 1.D0          ! new variable derivative
      		ALLXrow = ALLXrow+1      	
      		YY(irow) = -(ALLX(1,ALLXrow) - numerator/denominator)		! (for NEWTON) Current value of new variable is stored in ALLX(1,irow)
1310  		continue
!		ii = ii + numPhCo(k) - 1	! increment ii to skip these phase components
!		else
!		ii = ii + numPhCo(k) - 1	! just increment ii to skip these phase components
		endif 
	ii = ii + numPhCo(k) - 1	! increment ii to skip these phase components
1305  	continue
1200  	continue
!     Finished with extra equations

! 	Set up equation for Pfluid = hydrostatic, if needed
! 	Note: this equation does not incorporate fluiddensity as a derivative of T and P
! 	Also note: there is no derivative of fluid density wrt fluid composition (for mixed fluids)
	if(PFluidSwitch.eq.2)then
		call GetFluidDensity(iok)
		if(iok.eq.1)then
			call fss_alert('ALERT!!','P or T out of range for fluid equation of state')
			izero = 1
			go to 9999
			endif
		AA(NEQ,2) = -1.0d0*(FluidDensity/RockDensity)
		AA(NEQ,3) = 1.0D0
! 			   pfluid	Prock
		YY(NEQ) = -(ALLX(1,3) - AllX(1,2)*FluidDensity/RockDensity)  ! (for NEWTON) Current value of new variable is stored in ALLX(1,irow)
		endif	

!     If we are in MODE=NEWTON then we need to set up a data vector that
!     contains -(H - TS + R T lnKeq) for each reaction
!     (note there is no V term here because H and S are calculated at T and P
      	if(inewton.eq.1) then
		ii = np-nrx      ! note that reaction coefficients are stored in the end of array ARX
		DO 1410 i=1,NRX
		TEMPH=0.D0
		TEMPS=0.D0
		TEMPK=0.D0
		jj = 0
!		do 1415 k = 1,numPh
		Do 1415 kCur = 1,numPh
		k = asmCurrent(kCur)
		do 1415 J = 1,numPhCo(k)
		jj = jj + 1
		TEMPH =TEMPH + ARX(i+ii,jj)*HatTP(k,j)	!        H reaction
		TEMPS =TEMPS + ARX(i+ii,jj)*SatTP(k,j)	!        S reaction
		TEMPK = TEMPK + ARX(i+ii,jj)*lnAct(k,j)	!        LnK reaction
1415     	continue      
         	YY(i)=-(TEMPH - TK * TEMPS + R*TK*TEMPK)
1410     	continue
         	endif      

!     	Rearrange MASTER matrix AA
!     	MOVE MONITOR PARAMETERS TO RIGHT SIDE OF EQUNS IN ARRAY A
	nsolve=nvar
	if(inewton.eq.1)nsolve=nvar+1
	if(noOxygenMassBalance.eq.1)then
		NEQ = NEQ - 1		! this is reset after matrix inversion
		write(*,*)' NEQ RESET = ',NEQ
	      endif
	DO 2100 I=1,NVAR-NEQ
	DO 2100 J=1,NEQ
  	A(J,NEQ+I)=-AA(J,MON(I))
2100	continue
!     	FILL LEFT SIDE OF EQUATIONS IN MATRIX A FROM AA
      	DO 2110 I=1,NEQ
      	DO 2110 J=1,NEQ
	A(J,I)=AA(J,IPOINT(I))
2110	continue
      	if(inewton.eq.1)then
!           Fill Nvar+1=nsolve array element in A with the Y vector      
            DO 2115 J=1,NEQ
2115        A(J,Nsolve)=YY(J)
            endif
            

!     	Output matrix and modified matrix, if desired
	if(iLong(5).eq.1)then
2330  		CONTINUE
      		write(12,*)' '
      		write(12,*)' '
      		write(12,*)' MASTER MATRIX'
	        write(12,2321)(VN1(j),VN2(j),J=1,nvar),'  YY(i)'
2321  		FORMAT('     ',30(A2,A4,6x))
	      	DO 2310 I=1,NEQ
!		write(12,*)'--------------------'
		write(12,2320)(AA(I,J),J=1,NVAR),YY(I)
2310  		continue
2320  		FORMAT(' ',30E12.4)
		do 2311 i = 1,nvar-neq
		vn51(neq+i) = trim(vn1(mon(i)))//vn2(mon(i))
2311		continue
		do 2312 i = 1,neq
		vn51(i) = trim(vn1(ipoint(i)))//vn2(ipoint(i))
2312		continue
		write(12,*)' '
		write(12,*)' '
		write(12,*)' MODIFIED MASTER MATRIX'
	        write(12,2313)(VN51(j),J=1,nvar),' YY(i)'
2313  		FORMAT('     ',30(A6,6x))
		DO 2315 I=1,NEQ
	!         write(12,*)'--------------------'
		write(12,2320)(A(I,J),J=1,nsolve)
2315  		continue
		endif

!     Find solution
      	DETA=0.D0
      	IER=0
      	CALL REDUCE (NEQ,NSOLVE,DETA,IER)
      	IF (IER.EQ.1) then
		izero=1
		if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)go to 2509
		WRITE(*,*)' ************ ERROR **************************'
		WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
		WRITE(*,*)' The most likely cause is that the variables'
		WRITE(*,*)' specified as monitor parameters are not linearly'
		WRITE(*,*)' independent.  Note that only the number of'
		WRITE(*,*)' intensive variables equal to the phase rule'
		WRITE(*,*)' variance of the assemblage are independent.'
		WRITE(*,*)' Return to previous menu and try again.'
		write(*,*)' '
		pause 'Hit return to continue...'
		if(noOxygenMassBalance.eq.1)then
			NEQ = NEQ + 1		! reset
			write(*,*)' NEQ RESET = ',NEQ
			endif
		go to 9999
		endif

2509	continue

      	if(noOxygenMassBalance.eq.1)then
		NEQ = NEQ + 1		! reset to original value
		write(*,*)' NEQ RESET = ',NEQ
		izero = 1
		endif

	if(matrixIsSingular.eq.1)go to 9999

	if(iLong(2).eq.1)then
      		write(12,*)'   '
      		write(12,*)' RESULTS OF THIS FINITE DIFFERENCE ITERATION'
      		write(12,*)'DETERMINANT=',DETA
		endif
        SLOPE=XX(ndep,1)	! slope is used for contouring operations
	if(icont.eq.3)go to 9999	! return if only a slope is calculated
! --------------------CONTOUR BLOCK--------------
      	IF (ICONT.EQ.0) GO TO 8000
!       	IF PERFORMING A CONTOUR OPERATION
!     	check to see if user escaped with escape key
!      	call trapit(izero)
      	if(izero.gt.0)then
		izero = 2
		go to 9999
		endif
!     	ndep is the dependent plot variable=NYPLT
      	GO TO (8100,8200) ICONT
8100  	CONTINUE
!      	FIRST CALCULATION
!      	XINC AND YINC ARE SPECIFIED IN CONTOUR
      	IF (DABS(SLOPE).GT.ONE)THEN
!         	STEEP SLOPE
!         	USE NYPLT AS INDEP TO CALCULATE DELTA NXPLT
          	DELTAX(1)=YINC/SLOPE
       		ELSE
          	DELTAX(1)=XINC
       		ENDIF
       	GO TO 8016
8200   CONTINUE
!     	convert is the tolerance we are using for accepting an angle as ok
!     	convert = (5 degrees)*per cent of axes for precision (in radians)
!     	convert is calculated in subroutine contour
!      	Determine angle between the current slope and the previous slopw
!     	(note that angles are all in radians)
!     	one is the scaling factor for the plot of interest
      	if (dabs(slope).gt.one)then
            ANGLE = dabs(DATAN(one/PSLOP) - DATAN(one/SLOPE))
            else
            ANGLE = dabs(DATAN(PSLOP/one) - DATAN(SLOPE/one))
            endif
!      	ANGLE=DABS(ANGLE)
!      	IF(ANGLE.GT.Pi)ANGLE=2.d0*Pi - ANGLE
!       TEMP=ANGLE*180.D0/3.1415D0
!       DISPLAY IN DEGREES
!       WRITE(*,*)' LCONT=',LCONT,' Slope= ',slope,'  ANGLE =',temp
       	IF(ANGLE.GT. CONVERT)THEN
!       	The ANGLE is greater than our tolerance (5deg * %)
!       	we will LOOP BACK THROUGH compute decreasing the size of XINC and YINC
!       	until the angle is greater than 175 degrees
          	XINC=XINC/2.D0
          	YINC=YINC/2.D0
          	LCONT=LCONT+1
          	IF(LCONT.GE.5)GO TO 8053
!           	ONLY DO 5 TIMES
          	CALL RESET(6,4)
!          	DO NOT RESET PSLOP
          	GO TO 100
       		ENDIF
!     	if here then angle is less than convert.  First check to see
!     	if XINC and YINC can be made larger
!     	find the fraction of the angle as a fraction of convert and
!     	scale the increase in Xinc and yinc to that fraction
!     	the limit is a doubling of xinc and yinc if the angle is 0
       	XINC=DABS(XINC)
       	YINC=DABS(YINC)
      	temp = 1.d0 - (angle)/CONVERT
      	xinc = xinc + xinc*temp
      	yinc = yinc + yinc*temp
!     	check to be sure we have not exceeded the maximum values
        IF(XINC.GT.XINCST)XINC=XINCST
        IF(YINC.GT.YINCST)YINC=YINCST
8053  	continue
       	LCONT=0
       	IF (DABS(SLOPE).GT.ONE)THEN
!          	STEEP SLOPE
!          	USE (y axis=NYPLT) AS independent variable TO CALCULATE DELTA NXPLT
!          	NOTE THIS USES NYPLT BECAUSE OF STEEP SLOPE
!          	CDELY is the last change in Y.  If this is negative then
!           	the next change must also be negative
          	IF(CDELY.LT.0.D0)YINC=-Dabs(YINC)
          	DELTAX(1)=YINC/SLOPE
       		ELSE
!          	SHALLOW SLOPE; USE x axis=NXPLT
!            	NOTE CDELX IS LAST INCREMENT OF x axis
          	IF (CDELX.LT.0.0D0)XINC=-Dabs(XINC)
          	DELTAX(1)=XINC
       		ENDIF
8016  	CONTINUE
      	PSLOP=SLOPE
8000  	CONTINUE

! -------------------------------
!     	Set ALLX(5,L) equal to current values of ALLX(1,L)
!     	Note that below ALLX(1,L) is incremented by DELX
      	DO 3000 L=1,NVAR
3000  	ALLX(5,L)=ALLX(1,L)

! ____________________Newtons method___________________________
      	if(inewton.eq.1)then
!     		XX IS A NEQx(NVAR-NEQ+1) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		The last column (NVAR-NEQ+1) is the solution to the x for Newtons method
!     		Calculate the change in the variables
!     		New variables are computed as
!     		Xnew = Xold + X
! 		Note that we are doing these calculations with the values of the "independent" variables constant
! 		i.e. dxi = 0. For example, if it is univariant and we are doing calculations at constant P then dP = 0
! 		This way, we only need to solve for the values of the "dependent" variables (dT, dXi, etc.)
!

		! skip this newton step stuff -- it seems to be causing trouble
!		go to 5010

!       	Composition variables can only range from 0-1. 
!		If the solution vector component (xx) is too large, it may result in the change in mole fraction (xPhCo(i) exceeding this range
!		This block of code checks this and adjusts newtonStep to ensure this doesn't happen
!		newtonStep = 1.0d0

		newtonStep = newtonStepMax
		jCol = NVAR-NEQ+1	! index of solution vector in array XX
		do 5002 i = 1,nvar
		allXtemp(i) = allX(1,i)	! load up allXTemp
5002		continue
! ----		loop to here for newtonStep iteration --------
5001		continue
		do 5004 i = 1,NEQ
		ii=ipoint(i)
!     		adjust the values of AllXtemp that we have just calculated
!		note that this changes all dependent values, not just the composition ones
!		we could make this a wee bit faster if we just calculated the composition ones but we're not sure which of them are dependent
!		jcol is the column where the solution vector resides in XX
		ALLXtemp(ii)=ALLX(1,ii) + newtonstep*XX(i,jCol)
5004		continue
!		now pick out the composition values and store in xTemp
		l = TandP		! this is where the independent composition derivatives start
		Do 5010 kCur = 1,numPh
		k = asmCurrent(kCur)
		if(numPhCo(k).eq.1)go to 5010	! no composition derivatives for 1 component phases
		sum = 1.0d0		! this is for the dependent composition derivative (the first one)
		DO 5020 J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
		l = l + 1		! index of independent composition variable in AllX
!		j indexes the phase components of phase K
		xTemp(j) = allXTemp(l)
		sum = sum - xTemp(j)
5020  		continue
      		xTemp(1) = sum
		ier = 0
		call CheckStoichiometry(k,xTemp,ier)
		if(ier.eq.1)then			! failed stoichiometry test
			newtonStep = newtonStep/2.0d0
			!if(newtonStep.gt.1.0d-6)go to 5001			! try again with this smaller stepsize
			if(newtonStep.gt.1.0d-50)go to 5001			! try again with this smaller stepsize
!			stepsize is too small. Abort this phase
			if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)then
				izero = 1
				return
				endif
			write(*,*)'stepSize = ',newtonStep
			write(*,*)'This is too small. Something must be wrong in the code. Line 733 of C2OMPUTE'
			write(*,*)'Phase = ',k
			write(*,*)phName(k)
			do 5026 j1 = 1,numPhCo(K)
			write(*,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
5026			continue
			pause ' hit return to continue'
			izero = 1		! set error flag
			return
			endif

5010  		continue
!		write(*,*)'NewtonStep = ',newtonStep
!     		check for convergence
!     		for convergence we check to see whether the x value is < tolerance*X value      
      		jj = NVAR-NEQ+1		! index of solution vector in array XX
      		imiss=0
      		do 5110 i = 1,NEQ
      		ii=ipoint(i)
!     		Newtonstep scales the solution vector so we can step towards the solution in smaller steps
      		ALLX(1,ii)=ALLX(1,ii) + newtonstep*XX(i,jj)
!		AllX(1,ii) is the current value of the dependent variables. xx(i,jj) is the solution vector (variable).
!		When variable is less than a fraction(TolNewton) of the dependent variable, we have converged
      		if(dabs(xx(i,jj)).gt.dabs(TolNewton*Allx(1,ii)))then
      			if(allFractl(ii).eq.0)then		! only count imiss if the phase is not fractionating (seeGibbs_FileInput.f line 520 or so)
      		      		imiss=imiss+1
				endif
			endif
5110  		continue
 		if(iopen.eq.1)then
			do 5120 i = 1,iopen
			openFlux(iOpena(i)) = allX(1,openPtr - 1 + i)
5120			continue
 		endif
!     		if imiss>0 then one or more values failed tolerance test.  Rearrange
!     		things and go do another iteration
      		go to 3200
      		endif
! ________________________________________________________________
      
! _______________________Gibbs Method__________________________
      	if(inewton.eq.0)then
!     		XX IS A NEQx(NVAR-NEQ) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		COMPUTE new values FOR EACH VARIABLE
!     		note that the old variables DELX and SDELX are now defined as
!           	DELX(i)  = ALLX(1,i) - ALLX(5,i)  :current minus last finite difference
!           	SDELX(i) = ALLX(1,i) - ALLX(2,i)  :current minus starting
!     		Definition of ALLX variable:
!     		1,50   current value
!     		2,50     starting value
!     		3,50   ref at start of contour
!     		4,50   user selected reference
!     		5,50   previous finite diff point
!     		6,50     (not used)
!     		Store new values for each independent variable (monitor parameter)
      		DO 3100 J=1,NVAR-NEQ
      		ii=mon(j)
      		ALLX(1,ii)=ALLX(1,ii) + DELTAX(J)
3100  		continue
!     		compute new values for each dependent variable
      		Do 3110 i = 1,NEQ
      		ii=ipoint(i)
      		DO 3110 J=1,NVAR-NEQ
      		ALLX(1,ii)=ALLX(1,ii) + XX(i,J) * DELTAX(J)
3110  		continue    
      		endif       ! gibbs/newton

! ________________________________________________________________
3200  	continue

! ______________BOTH_____________________________________________
!     OUTPUT DERIVATIVES FOR THIS FINITE DIFFERENCE ITERATION
	if(iLong(2).eq.1)then
		call WriteOutJacobian()
		write(12,*)' NewtonStep = ',NewtonStep
	   	endif

! ----------------------------
      	call setTPX
!     	WRITE TO TERMINAL AND DISK
      	IF(iLong(4).eq.1)then
         	write(12,*)
         	write(12,*)
         	write(12,*)'Calculated Deltas'
         	write(12,*)'                      .       Value  .This increment.   Running sum'
         	do 4018 i = 1,nvar
         	write(12,4060)I,VN1(I),VN2(I),ALLX(1,i),ALLX(1,i)-ALLX(5,i),ALLX(1,i)-ALLX(2,i)
4060     	FORMAT(I3,1X,A2,A16,3F15.4)
4018     	continue
		endif

	if(iLong(2)+iLong(3)+iLong(4)+iLong(5)+iLong(6).ge.1)then
		write(*,*) "Type 0 (enter) to continue with next step"
		write(*,*) "Type 1 (enter) to abort"
		write(*,*) "Type -1 (enter) to continue without extended output"
!		pause "Hit return to continue with next step"
		read(*,*)izero
		if(izero.eq.1)then
			go to 9999
		    	endif
		if(izero.lt.0)then
			iLong(2) = 0
			iLong(3) = 0
			iLong(4) = 0
			iLong(5) = 0
			iLong(6) = 0
			endif
		endif

!     DONE WITH THIS STEP
      	if(inewton.eq.1.and.newtonStepsOutput.eq.1)then
		call NewtonStepSave
		write(*,*)newtonCounter,imiss
		endif

      	if(inewton.eq.1.and.imiss.gt.0)then
 		newtoncounter = newtoncounter+1
		if(newtoncounter.gt.1000)then
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			write(*,*)'Newton Counter > 1000 iterations. Aborting....'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			Write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
			izero = 1
			go to 9999
			endif		
        	go to 100	! loop for another iteration
		endif
	
9999	continue		! exit routine

!	if(noOxygenMassBalance.eq.1)then
!		NEQ = NEQ + 1
!		endif
      	RETURN
      	END
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Compute4MMAD(izero)
	use MatrixArrays
      implicit none

! c     version modified June 8, 2015
! C     ROUTINE TO COMPUTE SLOPES AND DELTAS for use in MMAD routines. 
! C	This is the same code as in C2OMPUTE except
!		(1) everything except Newton is eliminated (no differential thermo or contouring)
!		(2) only iterate on composition because everything is done at the same T
!			It is necessary that a call to dudTPX is made previously (not made here)
! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "Singular.inc"
!	include "Solute.inc"
	include "Tangent.inc"
	include "MatrixSingular.inc"
! c*****************************************
!     local variables
      integer*4 i,ii,j,jj,k,L,LL,Lcont,jcol,irow,ier,j1,kCur,&
     & izero,imiss,nsolve,newtoncounter,ivol,iok
      REAL*8 DETA,R,SUM,Pi,TEMPT,TEMPP,TEMPX,TEMPH,TEMPS,TEMPK,TempPFluid,xTemp(20),allXTemp(100)
      CHARACTER vn51(100)*8
      DATA R,Pi/8.3144D0,3.141592654D0/


! --------------------------------------------

      IF(iLong(3).eq.1)then
        	 write(12,*)' '
        	 write(12,*)'------------------------------------------------------'
        	 write(12,*)'------------------------------------------------------'
        	 write(12,*)' Variable list '
        	 do 11 i = 1,nvar
        	 write(12,15)I,VN1(I),VN2(I)
15      	 FORMAT(I3,2X,A2,A16,5x,e15.4)
11      	 continue
        	 write(12,*)'NEQ=',NEQ,'  NVAR=',NVAR
        	 write(12,*)
        	 write(12,*)'VARIANCE =',NVAR-NEQ
        	 write(12,*)
        	 write(12,*)'MONITOR PARAMETERS    Delta X for monitors'
         	if(nvar-neq.gt.0)then
         		do 12 i = 1,nvar-neq
         		write(12,16)MON(I),VN1(MON(I)),VN2(MON(I)),Deltax(i)*NSTEP
16       		FORMAT(I3,2X,A2,A16,F12.3)
12       		continue
!        		 write(12,*)' DELTA Xi FOR MONITOR PHASES'
!        		 write(12,16)(DELTAX(I)*NSTEP,I=1,NVAR-NEQ)
!16      		 FORMAT(10F12.3)
         		else
         		write(12,*)' No monitor parameters (variance <= 0)'
         		endif
        	 write(12,*)'NUMBER OF FINITE DIFF. ITERATIONS=',NSTEP
        	 write(12,*)
        	 write(12,*)'INITIAL MINERAL DATA'
        	 write(12,17)TC,PB
17      	 format(' T START (DEG C)=',F10.2,'   P START=',F10.2)
		if(PFluidSwitch.gt.0)write(12,18)PFluid
18		format(' PFluid = ',F10.2)
		endif
90    CONTINUE
!     LCONT IS COUNTER IN CONTOUR BLOCK
      LCONT=0
      newtoncounter=0
! -----------------------------------------------------------
!	Main loop starts here
100   CONTINUE
      TK=TC+273.15D0
!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST

!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      izero=0
      ivol=0            !compute everything
!	This call will callculate all thermodynamic properties - the TP dependent ones and the composition dependent ones
!	We need this when we're using T or P as dependent variables. But not when they're independent (because they're already set)
!	call DumpEverything('Before dlnAdX .......')

	Do 110 kCur = 1,numPh
	k = asmCurrent(kCur)
	call dlnAdX(k,izero)
      	if(izero.gt.0)then
		write(95,*)'NewtonCounter = ',newtoncounter
		write(95,*)'NewtonStep    = ',newtonstep
	!	call DumpEverything('After dlnAdX')
!		pause
      		go to 9999
		endif
110	continue
!      CALL ALLKduTPX(ivol,izero)
!     COMPUTE PHASE VOLUMES
	Do 340 kCur = 1,numPh
	k = asmCurrent(kCur)
	VP0(K)=MP0(K)*VMOL(K)*10.D0
340	continue
342   	CONTINUE
!     	OUTPUT THERMO DATA TO DISK
      	IF(iLong(3).eq.1)then
350       	CONTINUE
		Do 390 kCur = 1,numPh
		k = asmCurrent(kCur)
          	write(12,*)
          	write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       	FORMAT(' ',I4,F8.1,4X,A32)
          	write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          	IF(IMASS.EQ.1) then
          	   write(12,507)MP0(K),VP0(K)
          	   endif
507       	FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          	write(12,501)(phCoName(k,j),J=1,numPhCo(K))
501       	FORMAT('          ',12(A4,6X))
          	write(12,504)(xPhCo(k,j),J=1,numPhCo(K))
504       	FORMAT('  COMP ',12(F10.4))
          	write(12,505)(Dexp(lnAct(k,j)),J=1,numPhCo(K))
505       	FORMAT(' Activ  ',12(E10.3))
          	write(12,506)(lnAct(k,j),J=1,numPhCo(K))
506       	FORMAT('  Ln(a) ',12(F10.5))
          	write(12,514)(hPhCoZero(k,j),J=1,numPhCo(K))
514       	Format(' hPhCoZero',12E15.7)
          	write(12,515)(HATTP(k,j),J=1,numPhCo(K))
515       	Format(' HATTP',12E15.7)
          	write(12,509)(sPhCoZero(k,j),J=1,numPhCo(K))
509       	Format(' SZERO',12F15.5)
          	write(12,510)(SATTP(k,j),J=1,numPhCo(K))
510       	Format(' SATTP',12F15.5)
          	write(12,508)(vPhCoZero(k,j),J=1,numPhCo(K))
508       	Format(' VZERO',12F15.5)
          	write(12,513)(VATTP(k,j),J=1,numPhCo(K))
513       	Format(' VATTP',12F15.5)
          	write(12,516)(GATTP(k,j),J=1,numPhCo(K))
516       	Format(' GATTP',12E15.7)
          	write(12,*)'dudTPX ARRAY:'
          	write(12,*)'     dudT   .    dudP   .    dudX2  .    dudX3  .    dudX4...'
          	do 512 j = 1,numPhCo(K)
           	write(12,511)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
511       	FORMAT(15E12.5)
512   		continue
408       	format(' ',a5,6E10.3)
551       	CONTINUE
390       	CONTINUE
		endif
395       CONTINUE


!     THIS PART OF THE PROGRAM SETS UP THE MASTER MATRIX
!     ZERO MASTER ARRAY
      DO 1000 J=1,NVAR
      DO 1000 I=1,NEQ
1000  AA(I,J)=0.0D0

!     SET UP equations for the NRX reactions
!     These equations are generated by multiplication of the matrix ARX(I,J)
!     by the matrix dudTPX:
!      AA = ARX(nrx,np) x dudTPX(np,nvar)
!     The partial derivatives in matrix dudTPX are:
!     In storage space 1: dudTPX(jj,1):
!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.
!     In storage space 2: dudTPX(jj,2):
!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.
!     In storage space 3 to (numPhCo(K)+1): dudTPX(jj,3), dudTPX(jj,3), dudTPX(jj,3) etc:
!           Note that there are only numPhCo(K)-1 partial derivatives and the subscripts range
!           from 3 to numPhCo(K)+1.
!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X
!       
!     The matrix ARX(nrx,npc) contains the stoichiometric coefficients for each phase component
!     in each of the reactions (subscript 1 refers to reaction number, subscript
!     2 refers to the phase component)
!     The matrix AA therefore contains in each cell terms such as
!     
!       (n(i)*u(i))/T
!       (n(i)*u(i))/P
!       (n(i)*u(i))/X2
!       (n(i)*u(i))/X3
!     for each independent variable (where n(j,i) are the stoichiometric
!     coefficients for phase component i in reaction j)
!       Each row of AA contains a different set of reaction coefficients.
!     The only catch to the calculation of AA is that the matrix dudTPX is
!     not stored in true matrix form.  That is, to save space the independent variables
!     are listed for each phase only.
!     This simply requires a bit of extra bookkeeping during the matrix multiplication
! c	New code to implement fluid pressure as independent variable
	ii = np-nrx      !note that reaction coefficients are stored in the bottom of array ARX
	DO 1010 i=1,NRX
	TEMPT=0.D0
	TEMPP=0.D0
	TEMPPFluid = 0.0D0
	j = 0
!	do 1015 k = 1,numPh
	Do 1015 kCur = 1,numPh
	k = asmCurrent(kCur)
	do 1015 jj = 1,numPhCo(k)
	j = j + 1
!     	dT coefficient
      	TEMPT =TEMPT + ARX(i+ii,j)*dudTPX(k,jj,1)				! S reaction
!     	dP coefficient
! 	if a fluid is present and PFluidSwitch>0 we must put PFluid*Vfluid into an independent variable
	if(PFluidSwitch.eq.0)then					! Fluid pressure is dependent
		TEMPP =TEMPP + ARX(i+ii,j)*dudTPX(k,jj,2)			! V reaction
		else							! Fluid pressure is independent
! 		now we must determine if "j" is the index for a fluid component
!		jj = Minpoint(KFlu) - 1
!	I'm not sure this code is correct....
		if(k.eq.kflu)then
			do 1016 LL = 1,numPhCo(Kflu)
!		if(j.eq.jj+LL)then					! This is a fluid component if true
			TEMPPFluid =TEMPPFluid + ARX(i+ii,j)*dudTPX(k,jj,2)	! V reaction
1016			continue
			else						! not a fluid component
			TEMPP =TEMPP + ARX(i+ii,j)*dudTPX(k,jj,2)		! V reaction
			endif
		endif
1015  	continue      
      	AA(i,1)=TEMPT
      	AA(i,2)=TEMPP
      	if(PFluidSwitch.gt.0)AA(i,3) = TEMPPFluid
1010  	continue

!     	dX(j,k) coefficients
	DO 1020 i=1,NRX
	L=TANDP
	jj = 0
!	Do 1021 k = 1,numPh
	Do 1021 kCur = 1,numPh
	k = asmCurrent(kCur)
	if(numPhCo(K).eq.1)then
		jj = jj + 1
		go to 1021	! there are no dX terms for a phase of fixed composition
		endif
	Do 1022 LL = 3,numPhCo(K)+1     	! this loops on the number of independent composition derivatives
	tempX = 0.D0
	do 1023 j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
	TEMPX=TEMPX + ARX(i+ii,jj+j)*dudTPX(k,j,LL)
1023  	continue
      	L=L+1
      	AA(i,L)=TEMPX
1022  	continue
	jj = jj + numPhCo(k)
1021  	continue
1020  	continue

!     	SET UP MASS BALANCE EQUATIONS
!      	if(imass.eq.0)go to 1100
!     	FIRST SET UP dM Columns (M=MOLE FRACTION OF PHASE)
!     	JCOL COUNTS WHERE dM TERMS ARE
      	JCOL = TANDP + NX
!     	loop through each phase (there is a dM term for each phase)
	Do 1110 kCur = 1,numPh
	k = asmCurrent(kCur)
      	JCOL = JCOL + 1
!     	LOOP THROUGH ALL mass balance equations for each SYSTEM COMPONENT
      	call bulkcomp				! get total moles of each system component (moles(i) in common block)
      	DO 1105 i = 1,nc
!      	if(noOxygenMassBalance.eq.1)then
!		go to 1105			! for grain boundary problems, skip mass balance for oxygen
!	      endif
!     	IROW counts the row number
      	iRow = i + nRx
      	sum=0.D0
!     	LOOP THROUGH ALL PHASE COMPONENTS for this phase
      	DO 1120 j=1,numPhCo(K)
1120  	sum = sum + comp(k,J,I)*xPhCo(k,J)
      	AA(IROW,JCOL)=SUM
      	YY(IROW) = -(moles(i) - openFlux(i) - molesStart(i))	! molesStart(i) is the bulk composition we want to fit
      								! moles(i) is the current bulk comp based on phase comp and M(k)
								! when they are equal, we have the correct X and M for the desired bulk comp
								! openFlux are open system components either dependent or independent variables
1105  	CONTINUE
1110  	CONTINUE	!loop for each phase

!     	NOW SET UP COEFFICIENTS FOR dX TERMS in mass balance equations
!     	Note that the dependent dX term is removed as above
!     	JCOL is the column for the dX term
      	JCOL=TANDP
!     	FIRST LOOP THROUGH ALL PHASES TO FIND independent dX TERMS
	Do 1130 kCur = 1,numPh
	k = asmCurrent(kCur)
      	IF(numPhCo(K).EQ.1) GO TO 1130
!     	IF HERE THEN WE HAVE FOUND A PHASE WITH dX TERMS (I.E. MORE THAN 1 COMPONENT)
!     	LOOP THROUGH ALL dX TERMS FOR THIS PHASE
      	DO 1140 J=2,numPhCo(K)
      	JCOL=JCOL+1
!     	LOOP THROUGH ALL MASS BALANCE EQUATIONS FOR EACH dX TERM
      	DO 1150 I=1,NC
!      	if(noOxygenMassBalance.eq.1)then
!	      go to 1150	! for grain boundary problems, skip mass balance for oxygen
!	      endif
      	IROW= I + NRX
	AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))
1150	continue
1140  	CONTINUE
1130  	CONTINUE

!	Check to see if one of the phases is a GrainBoundary. If so, then add the dM terms to the energy equations.
	if(grainBoundaryProblem)then	
	! nphGB is the grain boundary phase
	! Grain boundary enthalpies are modeled as 
	! H = H(se)/M + H(de)*M where (se) is strain energy and (de) is disorder energy
	! this results in a minimum at some value of M (which is proportional to width of grain boundary)
	! The derivative is
	! dH = ( -H(se)/M^2 + H(de) )dM
       	ii = np-nrx      			! note that reaction coefficients are stored in the bottom of array ARX
      	JCOL = TANDP + NX + grainBoundaryPhase	! JCOL is the column where M(GB) is
      	DO 1220 i=1,NRX
	K = grainBoundaryPhase		   	! is the grain boundary phase number in the input list
							! note code will need to be written if I change assemblages
      	tempX = 0.D0
      	do 1223 j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
      	TEMPX=TEMPX + ARX(i+ii,jj+j)*(-HseGB/MP0(K)**2 + HdeGB)
1223  	continue
      	AA(i,JCOL)=TEMPX
1220  	continue
	endif

!1100  continue   ! end of mass balance equations


!     	Set up open system equations
!      	if(iopen.NE.0)then
!           (note if iopen0 then imass must=1)
!            jcol = TANDP + NX + numPh
!            irow=nrx
!            do 1210 i=1,iopen
!            AA(irow+iopena(i),jcol+i) = -1.D0
!1210        continue
!            endif
! 	Set up equation for Pfluid = hydrostatic, if needed
! 	Note: this equation does not incorporate fluiddensity as a derivative of T and P
! 	Also note: there is no derivative of fluid density wrt fluid composition (for mixed fluids)
	if(PFluidSwitch.eq.2)then
		call GetFluidDensity(iok)
		if(iok.eq.1)then
			call fss_alert('ALERT!!','P or T out of range for fluid equation of state')
			izero = 10
			go to 9999
			endif
		AA(NEQ,2) = -1.0d0*(FluidDensity/RockDensity)
		AA(NEQ,3) = 1.0D0
! 			   pfluid	Prock
		YY(NEQ) = -(ALLX(1,3) - AllX(1,2)*FluidDensity/RockDensity)  ! (for NEWTON) Current value of new variable is stored in ALLX(1,irow)
		endif	

!     If we are in MODE=NEWTON then we need to set up a data vector that
!     contains -(H - TS + R T lnKeq) for each reaction
!     (note there is no V term here because H and S are calculated at T and P
!      	if(inewton.eq.1) then
		ii = np-nrx      ! note that reaction coefficients are stored in the end of array ARX
		DO 1410 i=1,NRX
		TEMPH=0.D0
		TEMPS=0.D0
		TEMPK=0.D0
		jj = 0
!		do 1415 k = 1,numPh
		Do 1415 kCur = 1,numPh
		k = asmCurrent(kCur)
		do 1415 J = 1,numPhCo(k)
		jj = jj + 1
		TEMPH =TEMPH + ARX(i+ii,jj)*HatTP(k,j)	!        H reaction
		TEMPS =TEMPS + ARX(i+ii,jj)*SatTP(k,j)	!        S reaction
		TEMPK = TEMPK + ARX(i+ii,jj)*lnAct(k,j)	!        LnK reaction
1415     	continue      
         	YY(i)=-(TEMPH - TK * TEMPS + R*TK*TEMPK)
1410     	continue
!         	endif      

!     	Rearrange MASTER matrix AA
!     	MOVE MONITOR PARAMETERS TO RIGHT SIDE OF EQUNS IN ARRAY A
!	nsolve=nvar
	nsolve=nvar+1
!	if(inewton.eq.1)nsolve=nvar+1
!	if(noOxygenMassBalance.eq.1)then
!		NEQ = NEQ - 1		! this is reset after matrix inversion
!		write(*,*)' NEQ RESET = ',NEQ
!	      endif
	DO 2100 I=1,NVAR-NEQ
	DO 2100 J=1,NEQ
  	A(J,NEQ+I)=-AA(J,MON(I))
2100	continue
!     	FILL LEFT SIDE OF EQUATIONS IN MATRIX A FROM AA
      	DO 2110 I=1,NEQ
      	DO 2110 J=1,NEQ
	A(J,I)=AA(J,IPOINT(I))
2110	continue
!      	if(inewton.eq.1)then
!           Fill Nvar+1=nsolve array element in A with the Y vector      
            DO 2115 J=1,NEQ
2115        A(J,Nsolve)=YY(J)
!            endif
            

!     	Output matrix and modified matrix, if desired
	if(iLong(5).eq.1)then
2330  		CONTINUE
      		write(12,*)' '
      		write(12,*)' '
      		write(12,*)' MASTER MATRIX'
	        write(12,2321)(VN1(j),VN2(j),J=1,nvar),'  YY(i)'
2321  		FORMAT('     ',30(A2,A4,6x))
	      	DO 2310 I=1,NEQ
!		write(12,*)'--------------------'
		write(12,2320)(AA(I,J),J=1,NVAR),YY(I)
2310  		continue
2320  		FORMAT(' ',30E12.4)
		do 2311 i = 1,nvar-neq
		vn51(neq+i) = trim(vn1(mon(i)))//vn2(mon(i))
2311		continue
		do 2312 i = 1,neq
		vn51(i) = trim(vn1(ipoint(i)))//vn2(ipoint(i))
2312		continue
		write(12,*)' '
		write(12,*)' '
		write(12,*)' MODIFIED MASTER MATRIX'
	        write(12,2313)(VN51(j),J=1,nvar),' YY(i)'
2313  		FORMAT('     ',30(A6,6x))
		DO 2315 I=1,NEQ
	!         write(12,*)'--------------------'
		write(12,2320)(A(I,J),J=1,nsolve)
2315  		continue
		endif

!     Find solution
      	DETA=0.D0
      	IER=0
      	CALL REDUCE (NEQ,NSOLVE,DETA,IER)
      	IF (IER.EQ.1) then
		izero=11
		if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)go to 2509
		WRITE(*,*)' ************ ERROR **************************'
		WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
		WRITE(*,*)' The most likely cause is that the variables'
		WRITE(*,*)' specified as monitor parameters are not linearly'
		WRITE(*,*)' independent.  Note that only the number of'
		WRITE(*,*)' intensive variables equal to the phase rule'
		WRITE(*,*)' variance of the assemblage are independent.'
		WRITE(*,*)' Return to previous menu and try again.'
		write(*,*)' '
		pause 'Hit return to continue...'
		if(noOxygenMassBalance.eq.1)then
			NEQ = NEQ + 1		! reset
			write(*,*)' NEQ RESET = ',NEQ
			endif
		go to 9999
		endif

2509	continue
	if(matrixIsSingular.eq.1)go to 9999
	if(iLong(2).eq.1)then
      		write(12,*)'   '
      		write(12,*)' RESULTS OF THIS FINITE DIFFERENCE ITERATION'
      		write(12,*)'DETERMINANT=',DETA
		endif
        SLOPE=XX(ndep,1)	! slope is used for contouring operations
	if(icont.eq.3)go to 9999	! return if only a slope is calculated
! -------------------------------
!     	Set ALLX(5,L) equal to current values of ALLX(1,L)
!     	Note that below ALLX(1,L) is incremented by DELX
      	DO 3000 L=1,NVAR
3000  	ALLX(5,L)=ALLX(1,L)
! ____________________Newtons method___________________________
!      	if(inewton.eq.1)then
!     		XX IS A NEQx(NVAR-NEQ+1) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		The last column (NVAR-NEQ+1) is the solution to the x for Newtons method
!     		Calculate the change in the variables
!     		New variables are computed as
!     		Xnew = Xold + X
! 		Note that we are doing these calculations with the values of the "independent" variables constant
! 		i.e. dxi = 0. For example, if it is univariant and we are doing calculations at constant P then dP = 0
! 		This way, we only need to solve for the values of the "dependent" variables (dT, dXi, etc.)
!
!       	Composition variables can only range from 0-1. 
!		If the solution vector component (xx) is too large, it may result in the change in mole fraction (xPhCo(i) exceeding this range
!		This block of code checks this and adjusts newtonStep to ensure this doesn't happen
!		newtonStep = 1.0d0
		newtonStep = newtonStepMax
		jCol = NVAR-NEQ+1	! index of solution vector in array XX
		do 5002 i = 1,nvar
		allXtemp(i) = allX(1,i)	! load up allXTemp
5002		continue
! ----		loop to here for newtonStep iteration --------
5001		continue
		do 5004 i = 1,NEQ
		ii=ipoint(i)
!     		adjust the values of AllXtemp that we have just calculated
!		note that this changes all dependent values, not just the composition ones
!		we could make this a wee bit faster if we just calculated the composition ones but we're not sure which of them are dependent
!		jcol is the column where the solution vector resides in XX
		ALLXtemp(ii)=ALLX(1,ii) + newtonstep*XX(i,jCol)
5004		continue
!		now pick out the composition values and store in xTemp
		l = TandP		! this is where the independent composition derivatives start
!		DO 5010 k = 1,numPh
		Do 5010 kCur = 1,numPh
		k = asmCurrent(kCur)
		if(numPhCo(k).eq.1)go to 5010	! no composition derivatives for 1 component phases
		sum = 1.0d0		! this is for the dependent composition derivative (the first one)
		DO 5020 J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
		l = l + 1		! index of independent composition variable in AllX
!		j indexes the phase components of phase K
		xTemp(j) = allXTemp(l)
		sum = sum - xTemp(j)
5020  		continue
      		xTemp(1) = sum
		ier = 0
		go to 5010
!	Skip the stoichiometry check because it slows things to a crawl
		call CheckStoichiometry(k,xTemp,ier)
		if(ier.eq.1)then			! failed stoichiometry test
			newtonStep = newtonStep/2.0d0
			!if(newtonStep.gt.1.0d-7)go to 5001			! try again with this smaller stepsize
			if(newtonStep.gt.1.0d-50)go to 5001			! try again with this smaller stepsize
!			stepsize is too small. Abort this phase
			if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)then
				izero = 12
				write(95,*)'NewtonStep = ',newtonStep
				return
				endif
			write(*,*)'stepSize = ',newtonStep
			write(*,*)'This is too small. Something must be wrong in the code. Line 1384 of Compute4MMAD'
			write(*,*)'Phase = ',k
			write(*,*)phName(k)
			do 5026 j1 = 1,numPhCo(K)
			write(*,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
5026			continue
			pause ' hit return to continue'
			izero = 12		! set error flag
			return
			endif

5010  		continue
!		write(*,*)'NewtonStep = ',newtonStep
!     		check for convergence
!     		for convergence we check to see whether the x value is < tolerance*X value      
      		jj = NVAR-NEQ+1		! index of solution vector in array XX
      		imiss=0
      		do 5110 i = 1,NEQ
      		ii=ipoint(i)
!     		Newtonstep scales the solution vector so we can step towards the solution in smaller steps
      		ALLX(1,ii)=ALLX(1,ii) + newtonstep*XX(i,jj)
!		AllX(1,ii) is the current value of the dependent variables. xx(i,jj) is the solution vector (variable).
!		When variable is less than a fraction(TolNewton) of the dependent variable, we have converged
      		if(dabs(xx(i,jj)).gt.dabs(TolNewton*Allx(1,ii)))then
      			if(allFractl(ii).eq.0)then		! only count imiss if the phase is not fractionating (seeGibbs_FileInput.f line 520 or so)
      		      		imiss=imiss+1
				endif
			endif
5110  		continue
 		if(iopen.eq.1)then
			do 5120 i = 1,iopen
			openFlux(iOpena(i)) = allX(1,openPtr - 1 + i)
5120			continue
 		endif
!     		if imiss>0 then one or more values failed tolerance test.  Rearrange
!     		things and go do another iteration
      		go to 3200
!      		endif
! ________________________________________________________________
! ________________________________________________________________
3200  	continue

! ______________BOTH_____________________________________________
!     OUTPUT DERIVATIVES FOR THIS FINITE DIFFERENCE ITERATION
	if(iLong(2).eq.1)then
		call WriteOutJacobian()
	   	endif

! ----------------------------
      	call setTPX
!     	WRITE TO TERMINAL AND DISK
      	IF(iLong(4).eq.1)then
         	write(12,*)
         	write(12,*)
         	write(12,*)'Calculated Deltas'
         	write(12,*)'                      .       Value  .This increment.   Running sum'
         	do 4018 i = 1,nvar
         	write(12,4060)I,VN1(I),VN2(I),ALLX(1,i),ALLX(1,i)-ALLX(5,i),ALLX(1,i)-ALLX(2,i)
4060     	FORMAT(I3,1X,A2,A16,3F15.4)
4018     	continue
		endif

	if(iLong(2)+iLong(3)+iLong(4)+iLong(5)+iLong(6).ge.1)then
		write(*,*) "Type 0 (enter) to continue with next step"
		write(*,*) "Type 1 (enter) to abort"
		write(*,*) "Type -1 (enter) to continue without extended output"
!		pause "Hit return to continue with next step"
		read(*,*)izero
		if(izero.eq.1)then
		    go to 9999
		    endif
		if(izero.lt.0)then
			iLong(2) = 0
			iLong(3) = 0
			iLong(4) = 0
			iLong(5) = 0
			iLong(6) = 0
			endif
		endif

!     DONE WITH THIS STEP
      	if(inewton.eq.1.and.newtonStepsOutput.eq.1)then
		call NewtonStepSave
		write(*,*)newtonCounter,imiss
		endif

      	if(inewton.eq.1.and.imiss.gt.0)then
		newtoncounter = newtoncounter+1
		if(newtoncounter.gt.1000)then
			write(95,*)'Newton Counter > 1000 iterations. Aborting....'
			izero = 13
!			pause 'Hit return to continue - Line 1500 of Sub Compute4MMAD'
			go to 9999
			endif		
        	go to 100	! loop for another iteration
		endif
	
9999	continue		! exit routine
      	RETURN
      	END

! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Compute4MDF(izero,logfileon)
	use MatrixArrays
      implicit none

! 	Routine to calculate Jacobian for MDF (OS) model phases
!	The method is this:
!	(1) Set up the equilibrium relations for the base assemblage (same as for any equilibrium calculation)
!	(2) Set up the MDF (parallel tangent) equations for the OS phase
!		Note that there are npc - 1 such independent equations
!	(3) Set up the mass balance equations including all phases
!
!	The variance of this system should be 2 + #OS phases so we can solve for changes in all phases given the growth of the OS phase
!		For example dMph/dMOSphase
!	version written June 8, 2018 by F. Spear

! 	Note that this is the same code as in C2OMPUTE 
!		(1) everything except Newton is eliminated (no differential thermo or contouring)
!		(2) only iterate on composition because everything is done at the same T
!			It is necessary that a call to dudTPX is made previously (not made here)
!	It is copied from COMPUTE4MAD and the OS equations have been added
! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "Singular.inc"
!	include "Solute.inc"
	include "Tangent.inc"
	include "MatrixSingular.inc"
! c*****************************************
!     local variables
      integer*4 i,ii,j,jj,k,L,LL,Lcont,jcol,irow,ier,j1,kCur,izero,imiss,nsolve,newtoncounter,ivol,iok
      integer*4 numOSdXterms,numdXterms,OSrow,OScol,jdep,numOSEQ,jjj,logfileon
      integer*4 numVarEqSave,numEqEqSave,numMonits
!      real*8 LastYY(  ),LastSlope(   )
      REAL*8 DETA,R,SUM,Pi,TEMPT,TEMPP,TEMPX,TEMPH,TEMPS,TEMPK,TempPFluid,xTemp(20),allXTemp(100)
      CHARACTER vn51(100)*8
      DATA R,Pi/8.3144D0,3.141592654D0/


!	Each OS phase adds 1 M variable plus numPhCo-1 composition variables
!	Each OS phase adds numPhCo-1 equations
	numOSdXterms = 0
	do 1001 kcur = 1,numOSPhases
	k = OSPhase(kcur)
	numOSdXterms = numOSdXterms + numPhCo(K)-1
1001	continue
	! we are saving the values of neq and nvar for the equilibrium assemblage here then resetting them when we exit
	numVarEQSave = nvar
	numEqEQSave = nEQ
	! If we're growing only garnet then 3 OSdXterms and 1 OsPhase
	nVar = nvar + numOsdXterms + numOSPhases
	nEQ = nEQ + numOsdXterms	
!	nVar = nvar + 4
!	nEQ = nEQ + 3	
	numMonits = nvar - neq


! --------------------------------------------

      IF(iLong(3).eq.1)then
        	 write(12,*)' '
        	 write(12,*)'------------------------------------------------------'
        	 write(12,*)'------------------------------------------------------'
        	 write(12,*)' Variable list '
        	 do 11 i = 1,nvar
        	 write(12,15)I,VN1(I),VN2(I)
15      	 FORMAT(I3,2X,A2,A16)
11      	 continue
        	 write(12,*)'NEQ=',NEQ,'  NVAR=',NVAR
        	 write(12,*)
        	 write(12,*)'VARIANCE =',NVAR-NEQ
        	 write(12,*)
        	 write(12,*)'MONITOR PARAMETERS    Delta X for monitors'
         	if(nvar-neq.gt.0)then
         		do 12 i = 1,nvar-neq
         		write(12,16)MON(I),VN1(MON(I)),VN2(MON(I)),Deltax(i)*NSTEP
16       		FORMAT(I3,2X,A2,A16,F12.3)
12       		continue
!        		 write(12,*)' DELTA Xi FOR MONITOR PHASES'
!        		 write(12,16)(DELTAX(I)*NSTEP,I=1,NVAR-NEQ)
!16      		 FORMAT(10F12.3)
         		else
         		write(12,*)' No monitor parameters (variance <= 0)'
         		endif
        	 write(12,*)'NUMBER OF FINITE DIFF. ITERATIONS=',NSTEP
        	 write(12,*)
        	 write(12,*)'INITIAL MINERAL DATA'
        	 write(12,17)TC,PB
17      	 format(' T START (DEG C)=',F10.2,'   P START=',F10.2)
		if(PFluidSwitch.gt.0)write(12,18)PFluid
18		format(' PFluid = ',F10.2)
		endif
90    CONTINUE
!     LCONT IS COUNTER IN CONTOUR BLOCK
      LCONT=0
      newtoncounter=0
! -----------------------------------------------------------
!	Main loop starts here
100   CONTINUE
      TK=TC+273.15D0
!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST

!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      izero=0
      ivol=0            !compute everything
!	This call will callculate all thermodynamic properties - the TP dependent ones and the composition dependent ones
!	We need this when we're using T or P as dependent variables. But not when they're independent (because they're already set)
	Do 110 kCur = 1,numPh
	k = asmCurrent(kCur)
	call dlnAdX(k,izero)
      	if(izero.gt.0)go to 9999
110	continue
!      CALL ALLKduTPX(ivol,izero) -- this call must be made outside of this subroutine
!     COMPUTE PHASE VOLUMES
	Do 340 kCur = 1,numPh
	k = asmCurrent(kCur)
	VP0(K)=MP0(K)*VMOL(K)*10.D0
340	continue
342   	CONTINUE
!     	OUTPUT THERMO DATA TO DISK
      	IF(iLong(3).eq.1)then
350       	CONTINUE
		Do 390 kCur = 1,numPh
		k = asmCurrent(kCur)
          	write(12,*)
          	write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       	FORMAT(' ',I4,F8.1,4X,A32)
          	write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          	IF(IMASS.EQ.1) then
          	   write(12,507)MP0(K),VP0(K)
          	   endif
507       	FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          	write(12,501)(phCoName(k,j),J=1,numPhCo(K))
501       	FORMAT('          ',12(A4,6X))
          	write(12,504)(xPhCo(k,j),J=1,numPhCo(K))
504       	FORMAT('  COMP ',12(F10.4))
          	write(12,505)(Dexp(lnAct(k,j)),J=1,numPhCo(K))
505       	FORMAT(' Activ  ',12(E10.3))
          	write(12,506)(lnAct(k,j),J=1,numPhCo(K))
506       	FORMAT('  Ln(a) ',12(F10.5))
          	write(12,514)(hPhCoZero(k,j),J=1,numPhCo(K))
514       	Format(' hPhCoZero',12E15.7)
          	write(12,515)(HATTP(k,j),J=1,numPhCo(K))
515       	Format(' HATTP',12E15.7)
          	write(12,509)(sPhCoZero(k,j),J=1,numPhCo(K))
509       	Format(' SZERO',12F15.5)
          	write(12,510)(SATTP(k,j),J=1,numPhCo(K))
510       	Format(' SATTP',12F15.5)
          	write(12,508)(vPhCoZero(k,j),J=1,numPhCo(K))
508       	Format(' VZERO',12F15.5)
          	write(12,513)(VATTP(k,j),J=1,numPhCo(K))
513       	Format(' VATTP',12F15.5)
          	write(12,516)(GATTP(k,j),J=1,numPhCo(K))
516       	Format(' GATTP',12E15.7)
          	write(12,*)'dudTPX ARRAY:'
          	write(12,*)'     dudT   .    dudP   .    dudX2  .    dudX3  .    dudX4...'
          	do 512 j = 1,numPhCo(K)
           	write(12,511)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
511       	FORMAT(15E12.5)
512   		continue
408       	format(' ',a5,6E10.3)
551       	CONTINUE
390       	CONTINUE
		endif
395       CONTINUE


!     THIS PART OF THE PROGRAM SETS UP THE MASTER MATRIX
!     ZERO MASTER ARRAY
!************** needs to be expanded to include OS phases
	DO 1000 J=1,NVAR
	DO 1000 I=1,NEQ
1000	AA(I,J)=0.0D0

!     SET UP equations for the NRX reactions
!     These equations are generated by multiplication of the matrix ARX(I,J)
!     by the matrix dudTPX:
!      AA = ARX(nrx,np) x dudTPX(np,nvar)
!     The partial derivatives in matrix dudTPX are:
!     In storage space 1: dudTPX(jj,1):
!     Temperature:      /T = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.
!     In storage space 2: dudTPX(jj,2):
!     Pressure   :      /P =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.
!     In storage space 3 to (numPhCo(K)+1): dudTPX(jj,3), dudTPX(jj,3), dudTPX(jj,3) etc:
!           Note that there are only numPhCo(K)-1 partial derivatives and the subscripts range
!           from 3 to numPhCo(K)+1.
!     Composition:      /X = (ideal)/X + (margules)/X + (reciprocal)/X
!       
!     The matrix ARX(nrx,npc) contains the stoichiometric coefficients for each phase component
!     in each of the reactions (subscript 1 refers to reaction number, subscript
!     2 refers to the phase component)
!     The matrix AA therefore contains in each cell terms such as
!     
!       (n(i)*u(i))/T
!       (n(i)*u(i))/P
!       (n(i)*u(i))/X2
!       (n(i)*u(i))/X3
!     for each independent variable (where n(j,i) are the stoichiometric
!     coefficients for phase component i in reaction j)
!       Each row of AA contains a different set of reaction coefficients.
!     The only catch to the calculation of AA is that the matrix dudTPX is
!     not stored in true matrix form.  That is, to save space the independent variables
!     are listed for each phase only.
!     This simply requires a bit of extra bookkeeping during the matrix multiplication
! c	New code to implement fluid pressure as independent variable
	ii = np-nrx      !note that reaction coefficients are stored in the bottom of array ARX
	DO 1010 i=1,NRX
	TEMPT=0.D0
	TEMPP=0.D0
	TEMPPFluid = 0.0D0
	j = 0
!	do 1015 k = 1,numPh
	Do 1015 kCur = 1,numPh
	k = asmCurrent(kCur)
	do 1015 jj = 1,numPhCo(k)
	j = j + 1
!     	dT coefficient
      	TEMPT =TEMPT + ARX(i+ii,j)*dudTPX(k,jj,1)				! S reaction
!     	dP coefficient
! 	if a fluid is present and PFluidSwitch>0 we must put PFluid*Vfluid into an independent variable
	if(PFluidSwitch.eq.0)then					! Fluid pressure is dependent
		TEMPP =TEMPP + ARX(i+ii,j)*dudTPX(k,jj,2)			! V reaction
		else							! Fluid pressure is independent
! 		now we must determine if "j" is the index for a fluid component
!		jj = Minpoint(KFlu) - 1
!	I'm not sure this code is correct....
		if(k.eq.kflu)then
			do 1016 LL = 1,numPhCo(Kflu)
!		if(j.eq.jj+LL)then					! This is a fluid component if true
			TEMPPFluid =TEMPPFluid + ARX(i+ii,j)*dudTPX(k,jj,2)	! V reaction
1016			continue
			else						! not a fluid component
			TEMPP =TEMPP + ARX(i+ii,j)*dudTPX(k,jj,2)		! V reaction
			endif
		endif
1015  	continue      
      	AA(i,1)=TEMPT
      	AA(i,2)=TEMPP
      	if(PFluidSwitch.gt.0)AA(i,3) = TEMPPFluid
1010  	continue

!     	dX(j,k) coefficients
	DO 1020 i=1,NRX
	L=TANDP
	jj = 0
!	Do 1021 k = 1,numPh
	Do 1021 kCur = 1,numPh
	k = asmCurrent(kCur)
	if(numPhCo(K).eq.1)then
		jj = jj + 1
		go to 1021	! there are no dX terms for a phase of fixed composition
		endif
	Do 1022 LL = 3,numPhCo(K)+1     	! this loops on the number of independent composition derivatives
	tempX = 0.D0
	do 1023 j=1,numPhCo(K)          	! this loops on the number of phase components in each K phase
	TEMPX=TEMPX + ARX(i+ii,jj+j)*dudTPX(k,j,LL)
1023  	continue
      	L=L+1
      	AA(i,L)=TEMPX
1022  	continue
	jj = jj + numPhCo(k)
1021  	continue
1020  	continue
	numdXterms = L


!     	SET UP MASS BALANCE EQUATIONS for equilibrium assemblage
!     	FIRST SET UP dM Columns (M=MOLE FRACTION OF PHASE)
!     	JCOL COUNTS WHERE dM TERMS start
      	JCOL = TANDP + NX
!     	loop through each phase (there is a dM term for each phase)
	Do 1110 kCur = 1,numPh
	k = asmCurrent(kCur)
      	JCOL = JCOL + 1
!     	LOOP THROUGH ALL mass balance equations for each SYSTEM COMPONENT
      	call bulkcomp				! get total moles of each system component (moles(i) in common block)
      	DO 1105 i = 1,nc
!     	IROW counts the row number
      	iRow = i + nRx
      	sum=0.D0
!     	LOOP THROUGH ALL PHASE COMPONENTS for this phase
      	DO 1120 j=1,numPhCo(K)
1120  	sum = sum + comp(k,J,I)*xPhCo(k,J)
      	AA(IROW,JCOL)=SUM
      	YY(IROW) = -(moles(i) - openFlux(i) - molesStart(i))	! molesStart(i) is the bulk composition we want to fit
      								! moles(i) is the current bulk comp based on phase comp and M(k)
								! when they are equal, we have the correct X and M for the desired bulk comp
								! openFlux are open system components either dependent or independent variables
1105  	CONTINUE
1110  	CONTINUE	!loop for each phase


!     	NOW SET UP COEFFICIENTS FOR dX TERMS in mass balance equations
!     	Note that the dependent dX term is removed as above
!     	JCOL is the column for the dX term
      	JCOL=TANDP
!     	FIRST LOOP THROUGH ALL PHASES TO FIND independent dX TERMS
	Do 1130 kCur = 1,numPh
	k = asmCurrent(kCur)
      	IF(numPhCo(K).EQ.1) GO TO 1130
!     	IF HERE THEN WE HAVE FOUND A PHASE WITH dX TERMS (I.E. MORE THAN 1 COMPONENT)
!     	LOOP THROUGH ALL dX TERMS FOR THIS PHASE
      	DO 1140 J=2,numPhCo(K)
      	JCOL=JCOL+1
!     	LOOP THROUGH ALL MASS BALANCE EQUATIONS FOR EACH dX TERM
      	DO 1150 I=1,NC
      	IROW= I + NRX
	AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))
1150	continue
1140  	CONTINUE
1130  	CONTINUE



!	Set up the MDF equations
!	These take the form 
!	µalm(grt) - uprp(Grt) = ualm(Matrix) - uprp(Matrix)
!	µsps(grt) - uprp(Grt) = usps(Matrix) - uprp(Matrix)
!	µgrs(grt) - uprp(Grt) = ugrs(Matrix) - uprp(Matrix)

!	That is find F = (µ(j)-µ(1)) - (µOnTan(j)-µOnTan(1)) = 0
!	where µ(j)-µ(1) is the tangent to the Gibbs surface for the phase and
!	µOnTan(j)-µOnTan(1) is the slope of the system tangent in the same composition direction
!
!	I'm going to set this up to solve for everything at once
!	(a) composition and amounts of equilibrium assemblage
!	(b) composition of OS phase(s)
!	There should be 3 degrees of freedom so I can specify ∆P=∆T=0 and set ∆MOS phase to solve for everything else

!	We want the Xj in the OS phase where they are equal
!
!	Using Newton's method we set up a system of equations of the form
!
!	Fj(Xo) + d(Fj(atXo))/dXj*Xj = 0
!
!	there is one equation for each independent composition variable (Xj).
!	All necessary values and derivatives are computed in subroutine dlnAdX ( taken from dudTPX)
!	
!	After a call to dlnAdX, we have for each phase component (even dependent ones)
!	lnAct(jj) = (jj)/R*T		(for ideal solutions)
!	duTPX(jj,2+j) = d/dXj where the derivatives are only for each independent phase component (there are numPhCo terms and numPhCo-1 derivatives)
!		(the first to array elements are the T and P derivatives of each  but we don't need those here so we don't calculate them)
!
!	So our equations are calculated as
!	F(Xoj) = (uo(j)-uo(1)) + R*T*(lnAct(j) - lnAct(1)) - (µOnTan(j)-µOnTan(1))
!	d(FXoj)/dXj = dudTPX(j,2+j) - dudTPX(1,2+j)
!
!	so our matrix is
!	d(F2Xo)/dX2 + d(F2Xo)/dX3 + d(F2Xo)/dX4 + ...   * X2 = -F2(Xo)
!	d(F3Xo)/dX2 + d(F3Xo)/dX3 + d(F3Xo)/dX4 + ...   * X3 = -F3(Xo)
!	d(F4Xo)/dX2 + d(F4Xo)/dX3 + d(F4Xo)/dX4 + ...   * X4 = -F4(Xo)
!

!	Note -- all OS phase code has 8000 labels to make it easier to identify
	OSRow = NRX + nc		! this is the start of the rows containing the OS equations: NRX is the number of independent reactions, nc is the number of mass balance equations
	OSCol = TandP + NX + numPh	! column where MDF equations start: NX is the number of dX terms, nph the number of Mphase terms
	do 8010 kcur = 1,numOSPhases	! we might eventually have more than 1 OS phase
	k = OSPhase(kcur)

	if(numPhCo(K).eq.1)then
!		no OS composition equations for a 1 component phase
!		gOnTanAU(K)   = uOnTanAU(k,1) 				not needed
!		gOfPhaseAU(K) = uZero(k,1)/atomNorm(k,1)		not needed
!		gDifferenceAU(K) = gOfPhaseAU(K) - gOnTanAU(K)		not needed
		go to 8010				! done with this phase
		endif

!	This is code for phases with more than 1 phase component
	jdep = 1
	numOSEq = numPhCo(K) - 1

!	These store the last values incase there is a change of slope for a phase with a solvus
!	I'm not sure how they will work, if at all
!	They may not be necessary, but they shouldn't hurt at this point
!	do 8020 j = 1,numOSEq
!	LastYY(j) = 0.0d0
!	LastYY(OSRow+j) = 0.0d0
!	do 8020 i = 1,numOSEq
!	LastSlope(OSRow+j,i) = 0.0d0
!8020	continue

	izero = 0
	! note that this routine doesn't work on MIF arrays
	call dLnAdX(K,izero)				! update the activities based on the current phase composition
	if(izero.gt.0)then
      		write(12,*)' error in Subroutine dLnAdX called from subroutine Compute4MDF in file Gibbs_Compute.f95'
		write(12,*)' Phase, comp',K,phName(k),(xPhCo(k,j),j=1,numPhCo(k))
      		pause ' Hit return to continue -- Line 1906 in Sub Compute4MDF'
      		endif
	if(logFileOn.eq.1.and.izero.eq.1)then
		write(95,*)' Error in Compute4MDF -> Call dLnAdX for OS phase'
		write(95,*)' Phase, comp',K,phName(k),(xPhCo(k,j),j=1,numPhCo(k))
		endif

!	So our equations are calculated as
!	F(Xoj) = (uo(j)-uo(1)) + R*T*(lnAct(j) - lnAct(1)) - (µOnTan(j)-µOnTan(1))
!	d(FXoj)/dXj = dudTPX(j,2+j) - dudTPX(1,2+j)

	! calculate 
!	This loop goes through each equation, row by row
	do 8030 j = 1,numOSEq   ! numEq = numPhCo(K) - 1
	jjj = j + 1	! we skip the first phase component
!	YY(j) = -(uzeroDeltaAU(k,jjj) + Rjoules*TK*((lnAct(k,jjj)/atomNorm(k,jjj)) - (lnAct(k,jdep)/atomNorm(k,jdep))) &
!                   - uOnTanDeltaAU(k,jjj))
	YY(OSRow+j) = -(uzeroDeltaAU(k,jjj) + R*TK*((lnAct(k,jjj)/atomNorm(k,jjj)) - (lnAct(k,jdep)/atomNorm(k,jdep))) &
                   - uOnTanDeltaAU(k,jjj))
	if(minRec(k).eq.118)then
		write(96,*)xPhCo(k,jjj),uzeroDeltaAU(k,jjj),uOnTanDeltaAU(k,jjj),   &
     			R*TK*((lnAct(k,jjj)/atomNorm(k,jjj)) - (lnAct(k,jdep)/atomNorm(k,jdep))),YY(j)
     		endif
	
	! remember that uzeroDelta(K,j) are differenes - the first j=1 is a dummy
!	A(j,numEq+1) = YY(j)
!	A(OSRow,numEq+1) = YY(OSRow)	! actually, I don't need this one here now because these are set in DO 2115 loop
!	This loop goes through each column (dX terms)
	do 8040 i = 1,numOSEq
	slope = dudTPX(k,jjj,2+i)/atomNorm(k,jjj) - dudTPX(k,jdep,2+i)/atomNorm(k,jdep)
!	If this phase may have a solvus, then check whether there is a solution
!	if(minRec(k).eq.118.or.minRec(k).eq.96)then		! these are paragonite and alkali feldspar
!		if(LastSlope(OSRow+j,i)*slope.lt.0.0d0)then		! the slope has changed sign
!			if(LastYY(OSRow+j)*YY(OSRow+j).gt.0.0d0)then		! the function has not changed sign
				! there is no solution - we need to flag this and go on
!				gDifferenceAU(K) = 1.0d4		! Make = 10000. just to move on
!				go to 8010
!				endif
!			endif			
!		endif
!	LastSlope(OSRow+j,OSCol+i) = slope
!	AA(j,i) = slope
	AA(OSRow+j,OSCol+i) = slope
!	A(j,i) = AA(j,i)   ------- done below in loop Do 2100	
8040	continue
!	LastYY(OSRow) = YY(OSRow)
8030	continue
	OSRow = OSRow + numOSEq
	OSCol = OSCol + numOSEq		! update the column we are working on

8010	continue	!end of loop for all OS phases





!	Mass balance for OS phases
!     	loop through each OS phase (there is a dM term for each phase)
!	The dMOSphase column is the last variable(s)
	Do 8110 kCur = 1,numOSPhases
	k = OSPhase(kCur)
!      	JCOL = JCOL + 1
	JCOL = OSCol + 1
!     	LOOP THROUGH ALL mass balance equations for each SYSTEM COMPONENT
!      	call bulkcomp				! get total moles of each system component (moles(i) in common block)
      	DO 8105 i = 1,nc
!     	IROW counts the row number
      	iRow = i + nRx
      	sum=0.D0
!     	LOOP THROUGH ALL PHASE COMPONENTS for this phase
      	DO 8120 j=1,numPhCo(K)
8120  	sum = sum + comp(k,J,I)*xPhCo(k,J)
      	AA(IROW,JCOL)=SUM
      	YY(IROW) = -(moles(i) - openFlux(i) - molesStart(i))	! molesStart(i) is the bulk composition we want to fit
      								! moles(i) is the current bulk comp based on phase comp and M(k)
								! when they are equal, we have the correct X and M for the desired bulk comp
								! openFlux are open system components either dependent or independent variables
8105  	CONTINUE
8110  	CONTINUE	!loop for each phase


!     	NOW SET UP COEFFICIENTS FOR dX TERMS in mass balance equations for OS phases
!     	Note that the dependent dX term is removed as above
!     	JCOL is the column for the dX term
      	JCOL=TANDP + NX + numPh 	! we start with the column after those with the equilibrium phases (TandP + dXterms + dMterms)
!     	FIRST LOOP THROUGH ALL PHASES TO FIND independent dX TERMS
	Do 8130 kCur = 1,numOSPhases
	k = OSPhase(kCur)
      	IF(numPhCo(K).EQ.1) GO TO 8130
!     	IF HERE THEN WE HAVE FOUND A PHASE WITH dX TERMS (I.E. MORE THAN 1 COMPONENT)
!     	LOOP THROUGH ALL dX TERMS FOR THIS PHASE
      	DO 8140 J=2,numPhCo(K)
      	JCOL=JCOL+1
!     	LOOP THROUGH ALL MASS BALANCE EQUATIONS FOR EACH dX TERM
      	DO 8150 I=1,NC
      	IROW= I + NRX
!	I'm not sure what MP0 for an OS phase would or should be? 
!	Presumably, we'll want to solve for some value here if it is the independent component, so it will be set
!	otherwise, is there any value other than 0? Seems OS phases, by definition, are given moles = 0	
	AA(IROW,JCOL)=MP0(K)*(COMP(k,J,I)-COMP(k,1,I))
8150	continue
8140  	CONTINUE
8130  	CONTINUE



! end of mass balance equations


!     	Set up open system equations
!      	if(iopen.NE.0)then
!           (note if iopen0 then imass must=1)
!            jcol = TANDP + NX + numPh
!            irow=nrx
!            do 1210 i=1,iopen
!            AA(irow+iopena(i),jcol+i) = -1.D0
!1210        continue
!            endif
! 	Set up equation for Pfluid = hydrostatic, if needed
! 	Note: this equation does not incorporate fluiddensity as a derivative of T and P
! 	Also note: there is no derivative of fluid density wrt fluid composition (for mixed fluids)
	if(PFluidSwitch.eq.2)then
		call GetFluidDensity(iok)
		if(iok.eq.1)then
			call fss_alert('ALERT!!','P or T out of range for fluid equation of state')
			izero = 1
			go to 9999
			endif
		AA(NEQ,2) = -1.0d0*(FluidDensity/RockDensity)
		AA(NEQ,3) = 1.0D0
! 			   pfluid	Prock
		YY(NEQ) = -(ALLX(1,3) - AllX(1,2)*FluidDensity/RockDensity)  ! (for NEWTON) Current value of new variable is stored in ALLX(1,irow)
		endif	

!     If we are in MODE=NEWTON then we need to set up a data vector that
!     contains -(H - TS + R T lnKeq) for each reaction
!     (note there is no V term here because H and S are calculated at T and P
	ii = np-nrx      ! note that reaction coefficients are stored in the end of array ARX
	DO 1410 i=1,NRX
	TEMPH=0.D0
	TEMPS=0.D0
	TEMPK=0.D0
	jj = 0
	Do 1415 kCur = 1,numPh
	k = asmCurrent(kCur)
	do 1415 J = 1,numPhCo(k)
	jj = jj + 1
	TEMPH =TEMPH + ARX(i+ii,jj)*HatTP(k,j)	!        H reaction
	TEMPS =TEMPS + ARX(i+ii,jj)*SatTP(k,j)	!        S reaction
	TEMPK = TEMPK + ARX(i+ii,jj)*lnAct(k,j)	!        LnK reaction
1415     continue      
         YY(i)=-(TEMPH - TK * TEMPS + R*TK*TEMPK)
1410     continue
!
!	Note: Y vector for OS phases is set up above 

!     	Rearrange MASTER matrix AA
!     	MOVE MONITOR PARAMETERS TO RIGHT SIDE OF EQUNS IN ARRAY A
	nsolve=nvar+1
	DO 2100 I=1,numMonits		!NVAR-NEQ
	DO 2100 J=1,nEq		!NEQ
  	A(J,NEQ+I)=-AA(J,MON(I))
2100	continue
!     	FILL LEFT SIDE OF EQUATIONS IN MATRIX A FROM AA
      	DO 2110 I=1,nEQ		!NEQ
      	DO 2110 J=1,nEQ		!NEQ
	A(J,I)=AA(J,IPOINT(I))
2110	continue
!      	if(inewton.eq.1)then
!           Fill Nvar+1=nsolve array element in A with the Y vector      
            DO 2115 J=1,nEQ		!NEQ
2115        A(J,Nsolve)=YY(J)
!            endif
            

!     	Output matrix and modified matrix, if desired
	if(iLong(5).eq.1)then
2330  		CONTINUE
      		write(12,*)' '
      		write(12,*)' '
      		write(12,*)' MASTER MATRIX'
	        write(12,2321)(VN1(j),VN2(j),J=1,nvar),'  YY(i)'
2321  		FORMAT('     ',30(A2,A4,6x))
	      	DO 2310 I=1,NEQ
!		write(12,*)'--------------------'
		write(12,2320)(AA(I,J),J=1,NVAR),YY(I)
2310  		continue
2320  		FORMAT(' ',30E12.4)
		do 2311 i = 1,nvar-neq
		vn51(neq+i) = trim(vn1(mon(i)))//vn2(mon(i))
2311		continue
		do 2312 i = 1,neq
		vn51(i) = trim(vn1(ipoint(i)))//vn2(ipoint(i))
2312		continue
		write(12,*)' '
		write(12,*)' '
		write(12,*)' MODIFIED MASTER MATRIX'
	        write(12,2313)(VN51(j),J=1,nvar),' YY(i)'
2313  		FORMAT('     ',30(A6,6x))
		DO 2315 I=1,NEQ
	!         write(12,*)'--------------------'
		write(12,2320)(A(I,J),J=1,nsolve)
2315  		continue
		endif

!     Find solution
      	DETA=0.D0
      	IER=0
      	CALL REDUCE (NEQ,NSOLVE,DETA,IER)
      	IF (IER.EQ.1) then
		izero=1
!		if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)go to 2509
		WRITE(*,*)' ************ ERROR **************************'
		WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
		WRITE(*,*)' The most likely cause is that the variables'
		WRITE(*,*)' specified as monitor parameters are not linearly'
		WRITE(*,*)' independent.  Note that only the number of'
		WRITE(*,*)' intensive variables equal to the phase rule'
		WRITE(*,*)' variance of the assemblage are independent.'
		WRITE(*,*)' Return to previous menu and try again.'
		write(*,*)' '
		pause 'Hit return to continue...'
		go to 9999
		endif

2509	continue
	if(matrixIsSingular.eq.1)go to 9999
	if(iLong(2).eq.1)then
      		write(12,*)'   '
      		write(12,*)' RESULTS OF THIS FINITE DIFFERENCE ITERATION'
      		write(12,*)'DETERMINANT=',DETA
		endif
        SLOPE=XX(ndep,1)	! slope is used for contouring operations
	if(icont.eq.3)go to 9999	! return if only a slope is calculated
! -------------------------------
!     	Set ALLX(5,L) equal to current values of ALLX(1,L)
!     	Note that below ALLX(1,L) is incremented by DELX
      	DO 3000 L=1,NVAR
3000  	ALLX(5,L)=ALLX(1,L)
! ____________________Newtons method___________________________
!      	if(inewton.eq.1)then
!     		XX IS A NEQx(NVAR-NEQ+1) MATRIX WITH THE SOLUTIONS OF THE
!     		SLOPE EQUATIONS FOR EACH MONITOR PARAMETER
!     		The last column (NVAR-NEQ+1) is the solution to the x for Newtons method
!     		Calculate the change in the variables
!     		New variables are computed as
!     		Xnew = Xold + X
! 		Note that we are doing these calculations with the values of the "independent" variables constant
! 		i.e. dxi = 0. For example, if it is univariant and we are doing calculations at constant P then dP = 0
! 		This way, we only need to solve for the values of the "dependent" variables (dT, dXi, etc.)
!
!       	Composition variables can only range from 0-1. 
!		If the solution vector component (xx) is too large, it may result in the change in mole fraction (xPhCo(i) exceeding this range
!		This block of code checks this and adjusts newtonStep to ensure this doesn't happen
!		newtonStep = 1.0d0
		newtonStep = newtonStepMax
		jCol = NVAR-NEQ+1	! index of solution vector in array XX
		do 5002 i = 1,nvar
		allXtemp(i) = allX(1,i)	! load up allXTemp
5002		continue
! ----		loop to here for newtonStep iteration --------
5001		continue
		do 5004 i = 1,NEQ
		ii=ipoint(i)
!     		adjust the values of AllXtemp that we have just calculated
!		note that this changes all dependent values, not just the composition ones
!		we could make this a wee bit faster if we just calculated the composition ones but we're not sure which of them are dependent
!		jcol is the column where the solution vector resides in XX
		ALLXtemp(ii)=ALLX(1,ii) + newtonstep*XX(i,jCol)
5004		continue
!		now pick out the composition values and store in xTemp
		l = TandP		! this is where the independent composition derivatives start
!		DO 5010 k = 1,numPh
		Do 5010 kCur = 1,numPh
		k = asmCurrent(kCur)
		if(numPhCo(k).eq.1)go to 5010	! no composition derivatives for 1 component phases
		sum = 1.0d0		! this is for the dependent composition derivative (the first one)
		DO 5020 J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
		l = l + 1		! index of independent composition variable in AllX
!		j indexes the phase components of phase K
		xTemp(j) = allXTemp(l)
		sum = sum - xTemp(j)
5020  		continue
      		xTemp(1) = sum
		ier = 0
		call CheckStoichiometry(k,xTemp,ier)
		if(ier.eq.1)then			! failed stoichiometry test
			newtonStep = newtonStep/2.0d0
			!if(newtonStep.gt.1.0d-7)go to 5001			! try again with this smaller stepsize
			if(newtonStep.gt.1.0d-50)go to 5001			! try again with this smaller stepsize
!			stepsize is too small. Abort this phase
			if(idoingpseudo.eq.1.and.ioutputpseudo.eq.0)then
				izero = 1
				go to 9999
!				return
				endif
			write(*,*)'stepSize = ',newtonStep
			write(*,*)'This is too small. Something must be wrong in the code.  Line 2202 of Compute4MDF'
			write(*,*)'Phase = ',k
			write(*,*)phName(k)
			do 5026 j1 = 1,numPhCo(K)
			write(12,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
5026			continue
			pause ' hit return to continue'
			izero = 1		! set error flag
			go to 9999
			!return
			endif

5010  		continue
!		write(*,*)'NewtonStep = ',newtonStep
!     		check for convergence
!     		for convergence we check to see whether the x value is < tolerance*X value      
      		jj = NVAR-NEQ+1		! index of solution vector in array XX
      		imiss=0
      		do 5110 i = 1,NEQ
      		ii=ipoint(i)
!     		Newtonstep scales the solution vector so we can step towards the solution in smaller steps
      		ALLX(1,ii)=ALLX(1,ii) + newtonstep*XX(i,jj)
!		AllX(1,ii) is the current value of the dependent variables. xx(i,jj) is the solution vector (variable).
!		When variable is less than a fraction(TolNewton) of the dependent variable, we have converged
      		if(dabs(xx(i,jj)).gt.dabs(TolNewton*Allx(1,ii)))then
      			if(allFractl(ii).eq.0)then		! only count imiss if the phase is not fractionating (seeGibbs_FileInput.f line 520 or so)
      		      		imiss=imiss+1
				endif
			endif
5110  		continue
! 		if(iopen.eq.1)then
!			do 5120 i = 1,iopen
!			openFlux(iOpena(i)) = allX(1,openPtr - 1 + i)
!5120			continue
!	 		endif
!     		if imiss>0 then one or more values failed tolerance test.  Rearrange
!     		things and go do another iteration
      		go to 3200
!      		endif
! ________________________________________________________________
! ________________________________________________________________
3200  	continue

! ______________BOTH_____________________________________________
!     OUTPUT DERIVATIVES FOR THIS FINITE DIFFERENCE ITERATION
	if(iLong(2).eq.1)then
		call WriteOutJacobian()
	   	endif

! ----------------------------
!      	call setTPX		! not sure what getting rid of this will do... may need to fix
!     	WRITE TO TERMINAL AND DISK
      	IF(iLong(4).eq.1)then
         	write(12,*)
         	write(12,*)
         	write(12,*)'Calculated Deltas'
         	write(12,*)'                      .       Value  .This increment.   Running sum'
         	do 4018 i = 1,nvar
         	write(12,4060)I,VN1(I),VN2(I),ALLX(1,i),ALLX(1,i)-ALLX(5,i),ALLX(1,i)-ALLX(2,i)
4060     	FORMAT(I3,1X,A2,A16,3F15.4)
4018     	continue
		endif

	if(iLong(2)+iLong(3)+iLong(4)+iLong(5)+iLong(6).ge.1)then
		write(*,*) "Type 0 (enter) to continue with next step"
		write(*,*) "Type 1 (enter) to abort this iteration"
		write(*,*) "Type 2 (enter) to abort entirely"
		write(*,*) "Type -1 (enter) to continue without extended output"
!		pause "Hit return to continue with next step"
		read(*,*)izero
		if(izero.eq.1)then
		    go to 9999
		    endif
		if(izero.lt.0)then
			iLong(2) = 0
			iLong(3) = 0
			iLong(4) = 0
			iLong(5) = 0
			iLong(6) = 0
			endif
		endif

!     DONE WITH THIS STEP
      	if(inewton.eq.1.and.newtonStepsOutput.eq.1)then
		call NewtonStepSave
		write(*,*)newtonCounter,imiss
		endif

      	if(inewton.eq.1.and.imiss.gt.0)then
		newtoncounter = newtoncounter+1
		if(newtoncounter.gt.1000)then
		!	write(*,*)'Newton Counter > 1000 iterations. Aborting....'
			izero = 1
			go to 9999
			endif		
        	go to 100	! loop for another iteration
		endif
	
9999	continue		! exit routine
	nvar = numVarEQSave
	nEQ = numEqEQSave
      	RETURN
      	END




