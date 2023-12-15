! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ComputePTXM(TC4,PB4,izero)
      implicit none
! --------------------------------------------------------------------
	include "Assemb.inc"
	include "Monit.inc"
! --------------------------------------------------------------------
!     local variables
      	integer*4 i,j,l,izero
      	REAL*4 TC4,PB4
! --------------------------------------------------------------------
      	TC=TC4
      	Pb=Pb4
      	ALLX(1,1)=TC
      	ALLX(1,2)=PB
      	ALLX(2,1)=TC
      	ALLX(2,2)=PB
! 	This is a real kluge.....
! 	To ensure that Mwater is always positive, I'm going to add some moles to H2O
! 	After calculating the stable phase assemblage, I will NOT include MH2O in the free energy calculation
! 	molesStart(6) = molesStart(6) + 10.
! ----------------------------
! Make system open to H2O
! Solve for m(H2O) necessary to make M(Water) = 0. 
! m(H2O) will be a dependent variable
! 	iopen = 1		! number of open system components
! 	nvar = nvar + 1		! increases variance by 1
! 	iopena(1) = 6 		! H2O is component 6 in the KFMASH system (this will need to be modified for other systems
! 	do 8110 i = 1,nvar
! 	if(trim(vn2(i)).eq.'Water') then
! 		smon(3) = i	! The third monitor parameter is M(H2O)
! 		sdel(3) = 0.0d0	! M(H2O) will be forced to stay zero
! 		go to 8111
! 	endif
! 8110	continue
! 	call fss_alert('Did not find water in phase list')
! 	return
! 8111	continue		! found water - continue
! ----------------------------
!      	set monitor parameters
      	do 3401 i=1,nvar-neq
      	deltax(i)=sdel(i)
3401   	mon(i)=smon(i)
!      	Set up array IPOINT to contain pointers to non-monitor parameters
      	J=0
      	Do 3410 i=1,NVAR
      	DO 3411 L=1,NVAR-NEQ
      	IF(MON(L).eq.i)go to 3410
3411   	continue
      	J=J+1
      	IPOINT(J) = I
3410   	CONTINUE
!     	THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
      	DO 3403 I=1,nvar
!        ALLX(3,I) IS WHERE calculations STARTED-save
         ALLX(3,I)=ALLX(1,I)
3403   	CONTINUE
! 	if(inewton.eq.1)then
!       	Store new values for each independent variable (monitor parameter)
        	DO 3481 J=1,NVAR-NEQ
        	i=mon(j)
        	ALLX(1,i)=ALLX(1,i) + DELTAX(J)
3481    	continue
        	call SetTPX
!         endif
!     THIS SETS REFERENCES FOR resetting if there is a Newtonstep decrease
      DO 3404 I=1,nvar
         ALLX(4,I)=ALLX(1,I)		! ALLX(4,i) is user-defined. 
3404   CONTINUE
	if(NEQ.eq.0) go to 3422		! in some cases (e.g. sil+grt+qtz+kfs+h20) there are no reactions
					! so we must skip calculations and just plot tie lines
        IZERO=0
!	call dumpeverything('In ComputePTXM just before call to compute....')
        call Compute4MMAD(izero)
3422	continue
	return
	end
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ComputePTXGrid(TC4,PB4,izero,isubtend)
      implicit none
!	This routine is never called!!! Why is it here??? Who knows???
! --------------------------------------------------------------------
	include "Assemb.inc"
	include "Monit.inc"
! --------------------------------------------------------------------
!     local variables
      	integer*4 i,j,l,izero,isubtend
      	REAL*4 TC4,PB4
! --------------------------------------------------------------------
      	TC=TC4
      	Pb=Pb4
      	ALLX(1,1)=TC
      	ALLX(1,2)=PB
      	ALLX(2,1)=TC
      	ALLX(2,2)=PB
! ----------------------------
!      	set monitor parameters
      	do 3401 i=1,nvar-neq
      	deltax(i)=sdel(i)
3401   	mon(i)=smon(i)
!      	Set up array IPOINT to contain pointers to non-monitor parameters
      	J=0
      	Do 3410 i=1,NVAR
      	DO 3411 L=1,NVAR-NEQ
      	IF(MON(L).eq.i)go to 3410
3411   	continue
      	J=J+1
      	IPOINT(J) = I
3410   	CONTINUE
!     	THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
      	DO 3403 I=1,nvar
!        ALLX(3,I) IS WHERE calculations STARTED-save
         ALLX(3,I)=ALLX(1,I)
3403   	CONTINUE
! 	if(inewton.eq.1)then
!       	Store new values for each independent variable (monitor parameter)
        	DO 3481 J=1,NVAR-NEQ
        	i=mon(j)
        	ALLX(1,i)=ALLX(1,i) + DELTAX(J)
3481    	continue
        	call SetTPX
!         endif
!     THIS SETS REFERENCES FOR resetting if there is a Newtonstep decrease
      DO 3404 I=1,nvar
         ALLX(4,I)=ALLX(1,I)		! ALLX(4,i) is user-defined. 
3404   CONTINUE
	if(NEQ.eq.0) go to 3422		! in some cases (e.g. sil+grt+qtz+kfs+h20) there are no reactions
					! so we must skip calculations and just plot tie lines
        IZERO=0
        call C2ompute(izero)
!	no need to reset anything because it is done in Tangent routine when izero=1
3422	continue
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ResetTangent()
      implicit none
!     SUBROUTINE TO RESET THE MINERAL COMPOSITIONS TO THEIR STARTING
!     VALUES
! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
!	include "Output.inc"
!	include "Plotted.inc"
! ****************************************
!     local variables
      integer*4 j,k,L,kCur
! real*4 XXX,YYY
! character*128 ppp
! --------------------------------------------------
!
      DO 1490 L=1,nvar
            ALLX(1,L)=ALLX(3,L)
            ALLX(5,L)=ALLX(1,L)     !set previous finite difference point equal to current point
1490  CONTINUE
      L=TANDP
	Do 25 kCur = 1,numPh
	k = asmCurrent(kCur)
         IF(numPhCo(K).EQ.1)GO TO 424
         xPhCo(k,1)=1.0d0
         DO 24 J=2,numPhCo(K)
            L=L+1
            xPhCo(k,J)=ALLX(1,L)
            xPhCo(k,1)=xPhCo(k,1)-xPhCo(k,J)
24       CONTINUE
424      CONTINUE
25    CONTINUE
!     RESET MINERAL ABUNDANCE MOLES
      IF(IMASS.EQ.0) GO TO 30
         L=TANDP+NX
	Do 35 kCur = 1,numPh
	k = asmCurrent(kCur)
            L=L+1
            MP0(K)=ALLX(1,L)
            MP1(K)=MP0(K)
            MP2(K)=MP0(K)
35       CONTINUE
30    CONTINUE

      TC=ALLX(1,1)
      PB=ALLX(1,2)
      IF(KFLU.NE.0)PFLUID=ALLX(1,3)

1000  continue
!      if(isub.eq.3)go to 1020
!      IF(IRESET.NE.6)THEN
!     calls to plot only to move plotter pen back to starting conditions
!         xxx=ALLX(1,NXPLT)
!         YYY=ALLX(1,NYPLT)
!         IF (NYPLT.EQ.2)yyy=YYY/1000.
!         call plot(XXX,YYY,0)
!	Call FSS_SetPort(3)
!         call symb(XXX,YYY,6,10)    !draw a symbol on the screen just to show we got there
!      ENDIF
!1020  continue

!      if(imass.eq.1.and.ireset.ne.6)then
!        PHASE VOLUMES RECOMPUTED HERE
!         izero=0
!         call duTPX(1,izero)
!        returns in J/bar; display in cc/mole
!        moles per cm**3 (note VMOL is in joule/bar =CC/10)
!        note: if this code is changed, it must also be changed in routine CHANGE
!         do 300 k = 1,numPh
!	Do 300 kCur = 1,numPh
!	k = asmCurrent(kCur)
!         vp0(k)=mp0(K)*vmol(k)*10.D0
!         vp1(k)=vp0(k)
!         vp2(k)=vp0(k)
!         if(FRACTL(k).eq.1. or .FRACTL(k).eq.2)then
!            ppp = 'Warning: The phase '//Trim(phname(K))//' is fractionating. Volumes and moles may not be reset correctly.'
!            call FSS_alert(ppp)
!            endif
! 300     continue
!	 masserroralert = 0  	!switch for massbalance error alert (reset to start conditions)
!      endif

      RETURN
      END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine CalcMuZero(izero)
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 k,j,ivol,izero

      ! compute thermodyanmic values of all phases at the P&T of interest
	! This routine will return GatTP = µo for every end member phase component (not the phase total)
	! With P&T constant we do not need to recompute the values for the end member phase components - once is enough
!
	TK = TC + 273.15d0
	call CalculateCPTemp(TK)
!	numPh = numPhMIF		! I need to set this because duTPX used to use numPh to do all phases. But I changed this
	do 10 K = 1,numPhMIF
	ivol = 0	!computer everything
	izero = 0
!      	CALL duTPX(k,ivol,izero) 		! all calculations in duTPX are done with molar units (not atom units)
      	CALL uZeroAtTP(k,ivol,izero) 		! all calculations in duTPX are done with molar units (not atom units)
      	if(izero.eq.1)then
      		write(*,*)' error in Subroutine uZeroAtTP called from Subroutine CalcMuZero in TangentRoutines.f'
      		pause ' Hit return to continue - Line 250 in file TangentRoutines'
      		return
      		endif
	! the j here refers to the independent components in an assemblage for each phase. 
	!	probably should be subscripted that way
	!	note that gattp is only computed for the assemblage under consideration, yet we need all values
	!	for subroutine paralleltotangent
10	continue
!	write(12,*)' &&&&&&&&&&&&&&&&&&&&&&&'
!	write(12,*)' CalcMuZero'
!	write(12,*)TK,Pb
	do 25 K = 1,numPhMIF
	do 20 j = 1,numPhCo(k)
	uZero(k,j) = gattp(k,j)
	uZeroAU(k,j) = gattp(k,j)/atomNorm(k,j)
!	write(12,*)phName(k),phCoName(k,j),gattp(k,j),sattp(k,j),vattp(k,j)
!	µZero stores the µzero value at this T and P for each phase component
	uZeroDelta(k,j) = uzero(k,j) - uzero(k,1)
	uZeroDeltaAU(k,j) = uzeroAU(k,j) - uzeroAU(k,1)
! 	uZeroDelta = (independent phase component - dependent phase component)
20	continue
25	continue
!	write(12,*)' &&&&&&&&&&&&&&&&&&&&&&&'

	return
	end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine CalcMuOnTangent
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 j,l,k
	real*8 sum
!----------------------------------------------------
! 	Given the current tangent plane, 
!	calculate the value of µ for each phase component on the tangent
!	These calculations are all done in atom units
	do k = 1,numPhMIF
		do j = 1,numPhCo(k)
			sum = 0.0d0
			do l = 1,nc
				sum = sum + tanPlane(l)*comp(k,j,l)/atomNorm(k,j)		! comp() is in molar units
				end do
			uOnTanAU(k,j) = sum
			end do
		end do

	do K = 1,numphMIF
		do j = 1,numPhCo(k)
			uOnTanDeltaAU(k,j) = uOnTanAU(k,j) - uOnTanAU(k,1)
			! uOnTanDelta = (independent phase component - dependent phase component)
			end do
		end do
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine ParallelToTangent(logFileOn)
!	routine to find composition of the phase where the tangent plane is parallel to the system tangent plane
!	That is find F = (µ(j)-µ(1)) - (µOnTan(j)-µOnTan(1)) = 0
!	where µ(j)-µ(1) is the tangent to the Gibbs surface for the phase and
!	µOnTan(j)-µOnTan(1) is the slope of the system tangent in the same composition direction
!
!	We want the Xj where they are equal
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
!	which we solve for Xj . then compute
!	Xnewj = Xoj + Xj
!	and repeat until we converge
!
!	The compositions of each phase at the minimum are stored in phCoatGMin(noPh,numPhCo)
!	The G of the phase is stored in gOfPhase(k)
!	The G on the tangent at this composition is stored in gOnTan(k)
!	The difference is stored in gDifference(k)
!	
	use MatrixArrays
	implicit none
! *****************************************
	include "Assemb.inc"
! 	include "Solute.inc"
	include "Tangent.inc"
	include "Newton.inc"
!	include "MIFarrays.inc"
! --------------------------------------------------------------------
	integer*4 i,k,ivol,j,jjj,izero,nsolve,ier,iflag,notconverged,jdep,j1
	integer*4 loopcount,logFileOn,looping,fixTheStepsize
	real*8 Xnew(phMax,phCoMax),deta,sum,stepsize,xTemp(20),LastYY(phCoMax),LastSlope(phCoMax,phCoMax),slope
	real*8 Rjoules
      	DATA Rjoules/8.3144D0/
! --------------------------------------------------------------------

	TK = TC+273.15
	Pfluid = Pb
	izero=0
	ivol=0            

	newtonStepMax = 0.1

!	If a phase displays a solvus, there are situations where there is no solution to the parallel tangent
!	This occurs because the equation to solve is a cubic and there are either 1 or 3 roots. 
!		When there is only 1 root there is no solution for the other phase
!		Two examples are muscovite-paragonite and alkali feldspars
!	We need to be able to identify these situations and stop trying to find a root
!	We will do this by examining the slope to the equations
!		F(Xoj) = (uo(j)-uo(1)) + R*T*(lnAct(j) - lnAct(1)) - (µOnTan(j)-µOnTan(1))
!		d(FXoj)/dXj = dudTPX(j,2+j) - dudTPX(1,2+j)
!	If the slope changes sign (eg negative to positive) and the function does not change sign, then
!		we don't have a root (remember we're trying to find the composition where the function = 0)
!		if the function is never zero we can't have a root.
!	I suppose there are situations where this algorithm won't work and we'll probably need to look at the matrix of second derivatives
!		to determine if it's positive or negative definite... etc.
!		For example, this might only work on a binary.
!	A more robust solution might be to find the minimum and maximum of the function and see if the function has different signs at these points.
!		That would ensure that there was a root.
!		

	do 5000 k = 1,numPhMIF
	if(tanOut(2).eq.1)then		
		write(12,*)'+++++++++++++++'
		write(12,*)'Phase,PhaseType = ',k,phName(k),minRec(k),PhaseType(k)
		write(95,*)'+++++++++++++++'
		write(95,*)'Phase,PhaseType = ',k,phName(k),minRec(k),PhaseType(k)
		endif

!	if(minRec(k).eq.118)then
!		stepsize = 0.01
!		fixTheStepsize = 1	! this will keep the stepsize from being set to 1.0 if we don't converge on the first iteration
!		else
!		stepsize = 1.0		! this will be adjusted below so we stay in bounds
!		fixTheStepsize = 0	! this will allow the stepsize to be set to 1 each new iteration
!		endif
	stepsize = 1.0		! this will be adjusted below so we stay in bounds
!	stepsize = 0.1		! this will be adjusted below so we stay in bounds
	fixTheStepsize = 0	! this will allow the stepsize to be set to 1 each new iteration
		
!	For a phase of fixed composition we need only to calculate the G on the tangent
!		at the same composition as the phase(k)
!	This is already done in Subroutine CalcMuOnTangent
!		Result is stored in uOnTan(j)
!	if(k.eq.tanPhase)then
!		gDifferenceAU(K) = 0.0d0		! this is the tangent plane - by definition, the difference is 0
!		go to 999				! done with this phase
!		endif	
	if(numPhCo(K).eq.1)then
		gOnTanAU(K)   = uOnTanAU(k,1)
		gOfPhaseAU(K) = uZero(k,1)/atomNorm(k,1)
		gDifferenceAU(K) = gOfPhaseAU(K) - gOnTanAU(K)
		!IF(k.eq.1)then		! quartz
		!	write(*,*)TC,PB,TanPlane(1),uOnTanAU(k,1),uZero(k,1),atomNorm(k,1),gofphaseAU(k),gDifferenceAU(k)
		!	pause ' this is for quartz -- return to continue (in sub Parallel to tangent)'
		!	endif
		go to 999				! done with this phase
		endif

!	This is code for phases with more than 1 phase component
	jdep = 1
	numEq = numPhCo(K) - 1
	loopcount = 0
	looping = 0
	do 3100 j = 1,numEq
	LastYY(j) = 0.0d0
	do 3100 i = 1,numEq
	LastSlope(j,i) = 0.0d0
3100	continue

!	loop on Newton's method begins here
4000	continue
	loopcount = loopcount+1
	izero = 0
	! note that this routine doesn't work on MIF arrays
	call dLnAdX(K,izero)				! update the activities based on the current phase composition
	if(izero.eq.1.and.iOutputPseudo.eq.1)then
      		write(12,*)' error in Subroutine dLnAdX called from subroutine ParallelToTangent in file TangentRoutines.f'
		write(12,*)' Loop count - ',loopcount,TC,PB
		write(12,*)' Phase, comp',K,phName(k),(xPhCo(k,j),j=1,numPhCo(k))
		call listCurrentAsm
		write(12,*)' Hit return to continue'
      		pause ' Hit return to continue - Line 452 in Sub ParallelToTangent'
      		endif
	if(logFileOn.eq.1.and.izero.eq.1)then
		write(95,*)' Error in ParallelToTangent -> Call dLnAdX'
		write(95,*)' Phase, comp',K,phName(k),(xPhCo(k,j),j=1,numPhCo(k))
		endif

!	So our equations are calculated as
!	F(Xoj) = (uo(j)-uo(1)) + R*T*(lnAct(j) - lnAct(1)) - (µOnTan(j)-µOnTan(1))
!	d(FXoj)/dXj = dudTPX(j,2+j) - dudTPX(1,2+j)

6001	continue	
	! calculate 
	do 4010 j = 1,numEq   ! numEq = numPhCo(K) - 1
	jjj = j + 1	! we skip the first phase component
!	YY(j) = -((uzeroDeltaMIF(jjj) + Rjoules*TK*(lnAct(jjj)-lnAct(jdep)))/atomNormMIF(jjj) - uOnTanDeltaMIFAU(jjj))
!	YY(j) = -(uzeroDeltaMIFAU(jjj) + Rjoules*TK*((lnAct(jjj)/atomNormMIF(jjj)) - (lnAct(jdep)/atomNormMIF(jdep))) 
!     &              - uOnTanDeltaMIFAU(jjj))
	YY(j) = -(uzeroDeltaAU(k,jjj) + Rjoules*TK*((lnAct(k,jjj)/atomNorm(k,jjj)) - (lnAct(k,jdep)/atomNorm(k,jdep))) &
     &              - uOnTanDeltaAU(k,jjj))
	if(minRec(k).eq.118)then
		write(96,*)xPhCo(k,jjj),uzeroDeltaAU(k,jjj),uOnTanDeltaAU(k,jjj),   &
     &			Rjoules*TK*((lnAct(k,jjj)/atomNorm(k,jjj)) - (lnAct(k,jdep)/atomNorm(k,jdep))),YY(j)
     		endif
	
	! remember that uzeroDelta(K,j) are differenes - the first j=1 is a dummy
	A(j,numEq+1) = YY(j)
	do 4020 i = 1,numEq
	slope = dudTPX(k,jjj,2+i)/atomNorm(k,jjj) - dudTPX(k,jdep,2+i)/atomNorm(k,jdep)
!	If this phase may have a solvus, then check whether there is a solution
	if(minRec(k).eq.118.or.minRec(k).eq.96)then		! these are paragonite and alkali feldspar
		if(LastSlope(j,i)*slope.lt.0.0d0)then		! the slope has changed sign
			if(LastYY(j)*YY(j).gt.0.0d0)then		! the function has not changed sign
				! there is no solution - we need to flag this and go on
				gDifferenceAU(K) = 1.0d4		! Make = 10000. just to move on
				go to 999
				endif
			endif			
		endif
	LastSlope(j,i) = slope
	AA(j,i) = slope
	A(j,i) = AA(j,i)	
4020	continue
	LastYY(j) = YY(j)
4010	continue

	if(tanOut(2).eq.1)then		
		write(12,*)' '
		write(12,*)'Y = ',(YY(j),j=1,numEq),(xPhCo(k,j),j=1,numPhCo(k))
		write(12,*)'A matrix'
		do 4040 i = 1,numEq
		write(12,80)(A(i,j),j=1,numEq+1)
80		format(50F15.5)
4040		continue
		endif		

!     Find solution
	nsolve = numEq + 1
	DETA=0.D0
	IER=0
	CALL REDUCE (numEq,NSOLVE,DETA,IER)
	IF (IER.EQ.1) then
	      WRITE(12,*)' ************ ERROR **************************'
	      WRITE(12,*)' Matrix failed to invert in SUBROUTINE REDUCE'
	      write(12,*)' We are in Subroutine ParallalToTangent in file TangentRoutines.f '
	      write(12,*)' Phase = ',K,phName(k)
	      izero=1
!	      write(12,*) 'Hit return to continue...'
!	      pause 'Hit return to continue...'
	      return
	      endif

	if(tanOut(2).eq.1)then		
		write(12,*)'Solution xx'
		write(12,80)(xx(j,1),j=1,numEq)
		endif		

!	Check to see if the stepsize is too large.
!	In Newton's method, the solution finds the roots of the equation. 
!		If the curvature of the function is too great, this can result in the calculated composition being
!			outside of the range 0-1
!		Check this and reduce stepsize until we stay within range.
!		The same code is used in Subroutine C2ompute

!	make initial calculation with original stepSize
5001	continue
	sum = 1.0d0
	do 5004 j = 1,numEq
	jjj = jdep+j		! this is component number 2, 3, 4, etc
	j1 = j+1
	xTemp(j1) = xPhCo(k,jjj) + stepsize*xx(j,1)			
	sum = sum - xTemp(j1)
5004	continue
	xTemp(1) = sum	!dependent composition variable (number 1)
	if(tanOut(2).eq.1)then
		write(12,*)'stepsize = ',stepsize		
		do 50 j = 1,numPhCo(K)
		write(12,*)phCoName(k,j),xPhCo(k,j),xTemp(j)
50		continue
		endif

	iflag = 0
	call CheckStoichiometry(k,xTemp,iFlag)
	if(iflag.eq.1)then			! failed stoichiometry test
		stepSize = stepSize/2.0d0
		if(stepSize.gt.1.0d-10)go to 5001			! try again with this smaller stepsize
!		stepsize is too small. Abort this phase
		write(12,*)'stepSize = ',stepSize
		write(12,*)'This is too small. Something must be wrong in the code SUB ParallelToTangent - line 533'
		write(12,*)'Phase = ',k
		write(12,*)phName(k)
		do 5026 j1 = 1,numPhCo(K)
		write(12,*)phCoName(k,j1),xPhCo(k,j1),xTemp(j1)
5026		continue
		write(12,*)' '
		write(12,*)'Y = ',(YY(j),j=1,numEq),(xPhCo(k,j),j=1,numPhCo(k))
		write(12,*)'A matrix'
		do 5040 i = 1,numEq
		write(12,80)(A(i,j),j=1,numEq+1)
!80		format(50F15.5)
5040		continue
		write(12,*)'Solution xx'
		write(12,80)(xx(j1,1),j1=1,numEq)
		write(12,*) ' hit return to continue with the next phase in the list'
	!	pause ' hit return to continue with the next phase in the list'
		if(tanOut(4).eq.1)call thermodata   ! to check how things are  stored
		go to 4999	! jump to the next phase in the assemblage and try that
		endif
!	if here, then iflag = 0 and we are OK
! 	calculate new compositions with this stepSize
	sum = 1.0d0
	do 4100 j = 1,numEq		!number of independent components
	jjj = jdep+j		! this is component number 2, 3, 4, etc
	xnew(k,jjj) = xPhCo(k,jjj) + stepsize*xx(j,1)
	sum = sum - xnew(k,jjj)
4100	continue
	xnew(k,jdep) = sum	!dependent composition variable (number 1)
	notConverged = 0
	do 4130 j = 1,numPhCo(K)
	jjj = j
	! This convergence criterion was designed for major elements
	! trace elements can have concentrations much smaller (e.g. 10-10). 
!	if(xPhCo(k,jjj).gt.1.d-4)then
!		if(Dabs(xnew(k,jjj)-xPhCo(k,jjj)).gt.1.d-5)then
!		if(Dabs(xnew(k,jjj)-xPhCo(k,jjj)).gt.1.d-5)then
!			notConverged = 1
!			endif
!	else
!		if(Dabs(xnew(k,jjj)-xPhCo(k,jjj)).gt.1.D-15)then
!			notConverged = 1
!			endif
!	endif
!	if(Dabs(xnew(k,jjj)-xPhCo(k,jjj)).gt.TolNewton*xPhCo(k,jjj))then	! TolNewton is the fraction of X to scale convergence criterion
	if(Dabs(xnew(k,jjj)-xPhCo(k,jjj)).gt.0.00001)then	! TolNewton is the fraction of X to scale convergence criterion
		notConverged = 1
		endif


4130	continue
!	if(minRec(k).eq.118)then	This was for paragonite
!		write(95,*)loopcount,xnew(k,1),xnew(k,2)
!		endif

	if(tanOut(2).eq.1)then
		if(notConverged.eq.0)then			! notConverged = 0 means we have converged
			write(12,*)'Converged     '
			do 4131 j = 1,numPhCo(k)
			write(12,83)phCoName(k,j),xPhCo(k,j),xnew(k,j),xnew(k,j)-xPhCo(k,j)
83			format(A8,T10,20F15.5)
4131			continue
			else
			write(12,*)'NOT converged '
			do 4132 j = 1,numPhCo(k)
			write(12,83)phCoName(k,j),xPhCo(k,j),xnew(k,j),xnew(k,j)-xPhCo(k,j)
4132			continue
			endif	
		write(12,*) 'hit return to continue'
		pause 'hit return to continue - Line 629 in Sub ParallelToTangent'
		endif
	if(notConverged.eq.1)then 		! we have not converged
		if(loopCount.gt.1000)then	! something is wrong - we should have converged but we havent
						! this can happen when paragonite bounces to muscovite and back (across the solvus)
			write(95,*)'         LoopCount = 1000  Stepsize,Phase = ',stepSize,k
			looping = looping + 1
			if(looping.gt.3)then
				write(95,*)'         Looping = 3 -- abort this calculation ',k
				gDifferenceAU(K) = 1.0d4		! Make = 10000. just to move on
				go to 999
				endif
			do 4136 j = 1,numPhCo(K)
!			xPhCo(k,j) = xPhCoLast(k,j)	! start over with last good compositions (this will probably oscillate)
			xPhCo(k,j) = xPhCoInitial(k,j)	! start over with initial compositions
4136			continue
			loopCount = 0			! reset loopcount		
			if(looping.eq.1)then
				stepsize = 0.01			! make stepsize small this time
				else
				stepsize = .001
				endif
			else				! Loopcount < 100
			do 4135 j = 1,numPhCo(K)
			xPhCo(k,j) = xnew(k,j)
4135			continue
			endif
		if(fixTheStepSize.eq.0)then		! only make stepsize large if we are not looping
			stepsize = 1.0
!			stepsize = 0.1
			endif
		go to 4000
		endif

	do 4140 j = 1,numPhCo(K)
	xPhCo(k,j)    = xnew(k,j)
4140	continue
!	Converged - now calculate activities one more time to get gPhase
	izero = 0
	call dLnAdX(K,izero)
	! this call should calculate the activities for phase K provided it was loaded as an assemblage
      if(izero.eq.1.and.iOutputPseudo.eq.1)then
	      	write(12,*)' error in Subroutine dLnAdX called from subroutine ParallelToTangent in file TangentRoutines.f'
		write(12,*)' Loop count - ',loopcount
      		write(12,*) ' This is second call - Hit return to continue'
      		pause ' This is second call - Hit return to continue'
      		return
      		endif
	if(logFileOn.eq.1.and.izero.eq.1)then
		write(95,*)' Error in ParallelToTangent -> Call dLnAdX after convergence (this should never be a problem) '
		write(95,*)' Phase, comp',K,phName(k),(xPhCo(k,j),j=1,numPhCo(k))
		endif
		
!	for phases with 1 component we skip down to here		
4220	continue

	! Calculate G of the phase at this composition and compare with G on the tangent
	gOnTanAU(K) = 0.0d0
	gOfPhaseAU(K) = 0.0D0
	do 4210 j = 1,numPhCo(K)
	gOnTanAU(K)   = gOnTanAU(K)   + xPhCo(k,j)*uOnTanAU(k,j)
	gOfPhaseAU(K) = gOfPhaseAU(K) + xPhCo(k,j)*(uzero(k,j) + Rjoules*TK*lnAct(k,j))/atomNorm(k,j)
	gDifferenceAU(K) = (gOfPhaseAU(K) - gOnTanAU(K))
4210	continue

!	End routine for this phase
999	continue
	if(tanOut(2).eq.1)then
		write(12,81)phName(K),gOfPhaseAU(K),gOnTanAU(K),gDifferenceAU(k),(xPhCo(k,j),j=1,numPhCo(K))
81		format(A10,T15,3F15.2,15F12.5)
		endif

4999	continue
	if(tanOut(2).eq.1)then
		pause 'Look at results for this phase and hit return when ready'
		endif
5000	continue
!	newtonStepMax = 1.
	newtonStepMax = .4




	return
	end
	

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine AdjustTangentToAsm
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
! --------------------------------------------------------------------
	integer*4 k,j,l,jjj,kCur
	real*8 sum,uPhCoCurrent(100)			! temporary storage
	real*8 Rjoules
	include "CoTransCommon.inc"
	DATA Rjoules/8.3144D0/
! --------------------------------------------------------------------
!	Calculate the transformation from phase components to system components and
!	calculate the  of each system component from this matrix using the µ
!	 of each phase comonent
!	note that we have just done a COMPUT calculation
	call SysCoToPhCoTrans()
	! matrix B(np,nc) now contains transformation matrix
	if(tanOut(5).eq.1)then
		write(12,*)' Transformation matrix (B)'
!		write(12,4211)(phCoName(j),j=1,np) 
4211		format(T10,30(2x,A6))
		do 412 l=1,nc
		write(12,4212)coName(l),(B(l,j),j=1,np) 
4212		format(A12,30F8.3)
412		continue
		endif	
! Compute the chemical potential of each phase component in this assemblage
	jjj = 0
	do 10 kCur = 1,numPh		! loop on all phases in current assemblage
	k = asmCurrent(kCur)
	do 11 j = 1,numPhCo(K)
	jjj = jjj + 1
!	uPhCoCurrentAU(jjj) = (gattp(jjj) + Rjoules*TK*lnAct(jjj))/atomNormMIF(jjjMIF)	! convert to atom units
	uPhCoCurrent(jjj) = (gattp(k,j) + Rjoules*TK*lnAct(k,j))
	!IF(K.EQ.1)THEN
	!	write(*,*)' debug uquartz AdjustTangenttoASM ',Rjoules,TK,phCoName(k,1),uPhCoCurrent(jjj),gattp(k,1),lnAct(k,1)
	!	pause ' hit return to continue'
	!	endif
	! note that we don't need AU because each µ is defined as one AU
	! NOTE: uPhCo is only computed for the equilibrium assemblage we just computed, and NOT for the entire MIF
11	continue
10	continue
	! calculate coefficients of tangent plane for each system component
	! this will be the energy of the system component calculated from the linear combination of the phase component energies
	! calculated using the transformation matrix
	if(np.ne.jjj)then
		write(*,*)' In sub AdjustTangentToAsm -- np not equal to jjj (np,jjj):  ',np,jjj
		Write(*,*)' TP ',TC,PB
		write(*,'(20A12,2x)')(phName(asmCurrent(k)),k=1,numCurrent)
		
		pause 'Hit return to continue -- Line 765 of Sub AdjustTangentToAsm'
		endif
	do 425 L = 1,nc
	sum = 0.0d0
	do 420 j = 1,np		! np should equal jjj
!	actually, it doesn't make any sense not to use all phase components - check this during testing
!	Looking at the transformation matrix (B) it appears that some columns are all 0 and there are nc columns with non-zero values
!	By looping through all np phase components we don't have to worry as to which columns are all zero (i.e. no sorting)
	sum = sum + B(L,j)*uPhCoCurrent(j)
420	continue
	tanPlane(L) = sum
425	continue
	return
	end
	






	
	
		
