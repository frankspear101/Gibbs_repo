! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ComputeJacobian(izero)
	use MatrixArrays
      implicit none

!     	May 5, 2020 by F. S. Spear
!	ROUTINE TO COMPUTE the Jacobian based on the tangent plane (apple tree) approach
!	Based on code from COMPUTE4MAD
!  *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "Singular.inc"
!	include "Solute.inc"
	include "Tangent.inc"
	include "MatrixSingular.inc"
!  *****************************************
!     local variables
      	integer*4 i,ii,j,jj,k,L,LL,jcol,jcol2,krow,ier,kCur,numPhCoAsm,				&
           izero,nsolve,newtoncounter,ivol,iok,k1,j2,neqJ,nvarJ,variance
	integer*4 isteps,icount
      	REAL*8 DETA,R,AllXJ(5,50),sum,NormalMoles(30)
      	CHARACTER vn51(100)*8
      	DATA R/8.3144D0/

! --------------------------------------------
	numPh = numCurrent
!	normalize moles to 1.0
	sum = 0.0d0
	do l = 1,nc
		sum = sum + molesStart(l)
		end do
	do l = 1,nc
		NormalMoles(l) = molesStart(l)/sum
		end do

!	Set the tangent plane based on the current assemblage and mineral compositions
!	call AdjustTangentToAsm			why am I making this call? It changes the tangent plane

!	Set initial values in array ALLXJ
	sum = 0
	krow = 0
	do l = 1,nc
		sum = sum + TanPlane(l) * NormalMoles(l)
		end do
	krow = krow + 1
	ALLXJ(1,krow) = sum
	vn1(krow) = 'dG'
	vn2(krow) = ' '
	do l = 1,nc
		krow = krow + 1
		ALLXJ(1,krow) = TanPlane(l)
		vn1(krow) = 'SYu'
		vn2(krow) = CoName(l)
		end do
	do kcur = 1,numPh
		k = asmCurrent(kcur)
		do jj = 1, numPhCo(k)
			krow = krow + 1
			ALLXJ(1,krow) = gattp(k,jj)/atomNorm(k,jj)	 
			vn1(krow) = 'u'
			vn2(krow) = PhCoName(k,jj)
			end do
		end do
	do kcur = 1,numPh
		k = asmCurrent(kcur)
		do jj = 2, numPhCo(k)		! this should skip the first phase component
			krow = krow + 1
			ALLXJ(1,krow) = XPhCo(k,jj)	 
			vn1(krow) = 'X'
			vn2(krow) = PhCoName(k,jj)
			end do
		end do
	write(*,*)'ALLXJ array'
	do i = 1,krow
!		ALLXJ(2,i) contains the starting values. Set equal to the initial values
		ALLXJ(2,i) = ALLXJ(1,i)
		write(*,*)i,vn1(i),vn2(i),AllXJ(1,i),ALLXJ(2,i)
		end do
	write(*,*)'   '
	write(*,*)'   '

90    CONTINUE
      	newtoncounter=0
	variance = nc - numph
	write(*,*)' Variance = ',variance
	write(*,*)' Input monitor parameters'
	do i = 1,variance
		read(*,*)mon(i)
		end do
	write(*,*)' Input delta values for monitor parameters'
	do i = 1,variance
		write(*,*)i
		read(*,*)Deltax(i)
		end do

! -----------------------------------------------------------
!	Main loop starts here
	isteps = 1
100   CONTINUE
!	icount = 1
	do 200 icount = 1,isteps

	TK=TC+273.15D0
!     	THESE CALLS CALCULATE THE THERMODYNAMIC
!     	PROPERTIES OF THE PHASES OF INTEREST
!     	SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      	izero=0
      	ivol=0            !compute everything
!	This call will callculate the activities of each phase component
!	We don't need the µo terms because they were already calculated and we're working at constant T&P

	Do kCur = 1,numPh
		k = asmCurrent(kCur)
		call dlnAdX(k,izero)
		if(izero.gt.0)then
			write(95,*)'NewtonCounter = ',newtoncounter
			write(95,*)'NewtonStep    = ',newtonstep
		!	call DumpEverything('After dlnAdX')
	!		pause
			go to 9999
			endif
		end do

!     THIS PART OF THE PROGRAM SETS UP THE MASTER MATRIX
!     ZERO MASTER ARRAY
      	DO J=1,NVAR
      		DO I=1,NEQ
			AA(I,J)=0.0D0
			end do
		end do
!	Variables are these
!	dG (system) - we will increment this through each iteration
!	dµ of each system component (nc total)
!	dµ of each phase component in the current "assemblage". There are numPhCoAsm total of these
!	dX of each phase component in the current "assemblage". There are numPhCoAsm total of these
!
!	Equations are these:
!	Definition of dGsys in terms of dµ of system components
!		These take the form 
!		dG	dµSiO2	dµMgO	dµFeO	
!		-1	.6	.1	.3	
!		where the coefficients for the dµ terms are the moles in the bulk composition
!	Definition of each phase component in the current "assemblage" in terms of system components
!		dG	dµSiO2	dµMgO	dµFeO	dµQ	dµFo	dµFa	dXFa
!		0	1	0	0	-1	0	0	0
!		0	.33333	.66667	0	0	-1	0	0
!		0	.33333	0     .66667	0	0	-1	0	
!
!	A Gibbs-Duhem equation for each phase of the form
!		dG	dµSiO2	dµMgO	dµFeO	dµQ	dµFo	dµFa	dXFa
!		0	0	0	0	0	.25	.75	0		G-D for olivine
!		0	0	0	0	1	0	0	0		G-D quartz
!		Where the coefficients are the mole fractions of the phase components
!		Note that the G-D equation for quartz ensures that µSiO2 is fixed
!
!	Differentials of the phase components: the derivative of
!		µ = µo + RTln(a) or du = (RT/a)dX
!		These have the form
!		dG	dµSiO2	dµMgO	dµFeO	dµQ	dµFo	dµFa	dXFa
!		0	0	0	0	0	0	-1	26009.5176	0		du/dxFa Olivine
!		Note that we only use derivatives of the independent phase components (i.e. Xfa)
!			The G-D equation takes care of this dependency
!		The coefficients are dµ/dX, which are typically calculated by finite difference for complex activities
!
!	As in the C2OMPUTE routine, the derivatives are stored in the array duTPX(k,j,l)
!	The partial derivatives in matrix dudTPX are:
!	In storage space 1: dudTPX(jj,1):
!		Temperature:     dµ/dT = - S(P,T,X)  -- Partial molar Entropy at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state entropy, plus the entropy of mixing plus any exess entropy
!                                              terms.
!	In storage space 2: dudTPX(jj,2):
!		Pressure   :      dµ/dP =   V(P,T,X)  -- Partial molar Volume at P, T and composition so it incorporates all of the T and P corrections
!                                              to the standard state volume, plus any exess volume
!                                              terms.
!	In storage space 3 to (numPhCo(K)+1): dudTPX(jj,3), dudTPX(jj,3), dudTPX(jj,3) etc:
!           Note that there are only numPhCo(K)-1 partial derivatives and the subscripts range
!           from 3 to numPhCo(K)+1.
!     Composition:      dµ/dX = (ideal)/X + (margules)/X + (reciprocal)/X
!       
!	Example: SiO2 - MgO - FeO with the phases quartz + olivine
!		dG	dµSiO2	dµMgO	dµFeO	dµQ	dµFo	dµFa	dXFa
!		-1	.6	.1	.3	0	0	0	0		G sys definition
!		0	1	0	0	-1	0	0	0		uQtz definition
!		0	.33333	.66667	0	0	-1	0	0		µFo definition
!		0	.33333	0     .66667	0	0	-1	0		µFa definition
!		0	0	0	0	0	.2	.8	0		G-D for olivine
!		0	0	0	0	1	0	0	0		G-D quartz
!		0	0	0	0	0	0	-1	24383		du/dxFa Olivine




	ll =10
!	Definition of Gsystem
	AA(1,1) = -1.0d0
	do l = 1,nc			! loop through every system component
		AA(1,l+1) = NormalMoles(l)		! I'm not sure if these need to be normalized -- probably they do
		end do
!	write(*,29)'Gsys   ',(AA(1,l),l=1,ll)

!	definitions of µphase components
	jcol = nc + 1			! column index of the last system component µ
	krow = 1			! row index just before the first µPhCo definition
	numPhCoAsm = 0			! counts total number of phase components
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
		do jj = 1,numPhCo(k)
			numPhCoAsm = numPhCoAsm + 1
			krow = krow + 1
			jcol = jcol + 1
			AA(krow,jcol) = -1.0d0
			do l = 1,nc
				AA(krow,l+1) = comp(k,jj,l)
				end do
	!		write(*,29)'mu PhCo',(AA(krow,l),l=1,ll)
			end do
		end do

!	Gibbs-Duhem equations
	jcol = nc + 1			! index of the last system component µ
	krow = numPhCoAsm + 1		! index of the last equation before the first G-D equation
	do kCur = 1,numPh
		k = asmCurrent(kCur)
		krow = krow + 1
		do jj = 1,numPhCo(k)
			jcol = jcol + 1
			AA(krow,jcol) = XPhCo(k,jj)
			end do
!		write(*,29)'GD EQun',(AA(krow,l),l=1,ll)
		end do
	
!     	dµ/dX(k,j) equations
	jcol = nc + 1			! index of the last system component µ
	jcol2 = nc + numPhCoAsm + 1	! index of the last duPhCo term (one column before the first dX column)
	krow = numPh + numPhCoAsm + 1		! index of the equation before first dµ/dX equation
	do kCur = 1,numPh
		k = asmCurrent(kCur)
		jcol = jcol + 1			! column counter for dµPhCo
		if(numPhCo(k).eq.1)then
			goto 25
			endif
		do jj = 2,numPhCo(k)
			krow = krow + 1			! equation counter
			jcol = jcol + 1			! column counter for duPhCo
!			write(*,*)k,jj,krow,jcol,jj,jcol2
			AA(krow,jcol) = -1.0d0
			jcol2 = jcol2 + 1
			AA(krow,jcol2) = dudTPX(k,jj,jj+1)/atomNorm(k,jj)	! convert to atom units
			end do
!		write(*,29)'du/dX  ',(AA(krow,l),l=1,ll)
29		format(A7,2x,30E12.4)
25		continue
		end do

	NEQJ  = krow
	NVARJ = Jcol2
	nsolve = nvarj + 1
	variance = nvarj-neqj
!	write(*,*)' Variance = ',nvarj-neqj
	k = 0
	do i = 1,nvarj
		do j = 1,variance
			if(i.eq.mon(j))go to 31
			end do
		k = k + 1
		ipoint(k) = i
31		continue
		end do
!	write(*,*)'Monitors ',(mon(i),i=1,variance)
!	write(*,*)'Non-Mons ',(ipoint(i),i=1,neqj)
	!pause 'hit return to continue'

!     	Rearrange MASTER matrix AA
!     	MOVE MONITOR PARAMETERS TO RIGHT SIDE OF EQUNS IN ARRAY A
	DO I=1,variance
		DO k=1,NEQj
  			A(k,NEQj+I)=-AA(k,MON(I))
			end do
		end do

!     	FILL LEFT SIDE OF EQUATIONS IN MATRIX A FROM AA
      	DO I=1,NEQj
	      	DO k=1,NEQj
			A(k,I)=AA(k,IPOINT(I))
			end do
		end do

!     	Output matrix and modified matrix, if desired
	if(icount.eq.isteps)then
      		write(12,*)' '
      		write(12,*)' '
      		write(12,*)' MASTER MATRIX'
	        !write(12,2321)(VN1(j),VN2(j),J=1,nvarj),'  YY(i)'
2321  		FORMAT('     ',30(A2,A4,6x))
	      	DO I=1,NEQj
!			write(12,*)'--------------------'
			write(12,2320)(AA(I,J),J=1,NVARj)
2320  			FORMAT(' ',30E12.4)
			end do
!		do 2311 i = 1,variance
		!vn51(neq+i) = trim(vn1(mon(i)))//vn2(mon(i))
!2311		continue
!		do 2312 i = 1,neqj
!		!vn51(i) = trim(vn1(ipoint(i)))//vn2(ipoint(i))
!2312		continue
		write(12,*)' '
		write(12,*)' '
		write(12,*)' MODIFIED MASTER MATRIX'
	        !write(12,2313)(VN51(j),J=1,nvarj),' YY(i)'
2313  		FORMAT('     ',30(A6,6x))
		DO I=1,neqj
		!       write(12,*)'--------------------'
			write(12,2320)(A(I,J),J=1,nvarj)
			end do
		endif

	!pause 'Hit return to continue'
!     Find solution
      	DETA=0.D0
      	IER=0
      	CALL REDUCE (NEQJ,NSOLVE,DETA,IER)
      	IF (IER.EQ.1) then
		izero=11
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

	if(matrixIsSingular.eq.1)go to 9999

	if(icount.eq.isteps)then
      		write(12,*)'   '
      		write(12,*)' RESULTS OF THIS FINITE DIFFERENCE ITERATION'
      		write(12,*)'DETERMINANT=',DETA
		endif


!     OUTPUT DERIVATIVES FOR THIS FINITE DIFFERENCE ITERATION
	if(icount.eq.isteps)then
		write(12,*)
		write(12,3507)
3507    	FORMAT('    Jacobian matrix:')
		write(12,3502)(VN1(MON(I)),VN2(MON(I)),i=1,variance)
3502    	FORMAT(' Monitors:',5x,15(A2,A4,9x))
		write(12,3506)(deltax(I),i=1,variance)
3506    	format('      ',15(f15.4))
		JJ=0
		DO I=1,nvarj
			DO J=1,variance
				IF(I.EQ.MON(J))GO TO 3520
				end do
			JJ=JJ+1
			write(12,3501)VN1(I),VN2(I),(XX(JJ,LL),ll=1,variance)	! writes out soln in last XX
3501    		FORMAT(' ',A2,A4,' = ',15E15.5)
3520    		CONTINUE
			end do
	   	endif

!	I need to calculate the changes in all of the dependent variables.
!	This will include the compositions of the phases in the "assemblage" as well as the 
!	chemical potentials of the tangent plane

! -------------------------------
!     	Set ALLXJ(5,L) equal to current values of ALLXJ(1,L)
!     	Note that below ALLXJ(1,L) is incremented by DELX
      	DO L=1,NVARJ
		ALLXJ(5,L)=ALLXJ(1,L)
		end do
!	First calculate new compositions in array ALLXJ
!	These are dependent variables
      	DO j=1,variance
	      	ii=mon(j)
	      	ALLXJ(1,ii)=ALLXJ(1,ii) + DELTAX(j)
		end do

!     	compute new values for each dependent variable
      	Do i = 1,NEQJ
	      	ii=ipoint(i)
      		DO j=1,variance
      			ALLXJ(1,ii)=ALLXJ(1,ii) + XX(i,j) * DELTAX(j)
			end do
		end do

!	Change the tangent plane to the new values
	do l = 1,nc
		TanPlane(l) = ALLXJ(1,l+1)
		end do

	if(icount.eq.isteps)Call listtangentplane()		

!     set new values of mineral composition in phase component array X
!	This code is the same as in subroutine SETTPX
	L = 1 + nc + numPhCoAsm			! Last index before composition derivatives
	Do kCur = 1,numPh
		k = asmCurrent(kCur)
	      	sum=0.0d0
	      	DO j = 2,numPhCo(K) ! loop from second to last phase component for this mineral
		      	L=L+1
		      	xPhCo(k,j) = ALLXJ(1,L)
		      	sum = sum + xPhCo(k,j)
			end do
      		xPhCo(k,1) = 1.0D0-sum		! dependent compositional variable
		end do
	
! ----------------------------
!     	WRITE TO TERMINAL AND DISK
	if(icount.eq.isteps)then
         	write(12,*)
         	write(12,*)
         	write(12,*)'Calculated Deltas'
         	write(12,*)'                      .       Value  .This increment.   Running sum'
         	do i = 1,nvarj
         		write(12,4060)I,VN1(I),VN2(I),ALLXJ(1,i),ALLXJ(1,i)-ALLXJ(5,i),ALLXJ(1,i)-ALLXJ(2,i)
4060     		FORMAT(I3,1X,A2,A16,3F15.4)
			end do
		endif


! ________________________________________________________________
! ________________________________________________________________
!	Then run Parallel to tangent

!	Calculate parallel to tangent for current assemblage
!		(b) After just resetting the phases that define the tangent
!			This is done with routine AdjustTangentToAsm

!		write(12,*)' -------------------'
!		write(12,*)'Tangent before call to AdjustTangentToAsm '
!		tanOut(7) = 1
!		if(tanOut(7).eq.1)call ListTangentPlane
!		write(12,*)' -------------------'
!		call AdjustTangentToAsm
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
!		write(12,*)' -------------------'
!		write(12,*)'Tangent after call to AdjustTangentToAsm '
!		if(tanOut(7).eq.1)call ListTangentPlane
!		write(12,*)' -------------------'
		i = 0			! log file is off
		call ParalleltoTangent(i)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!		if(tanOut(3).eq.1)call ListGDifference
!		write(*,*)' quartz before PRINTT call ',gattp(1,1)
!	         CALL PRINTT(1)
!		write(*,*)' quartz after PRINTT call ',gattp(1,1)

200	continue

		call WriteAllStuffToTerminal(12)		! write results to output window

 
42	continue
!	if(iLong(2)+iLong(3)+iLong(4)+iLong(5)+iLong(6).ge.1)then
		write(*,*) "Type 0 (enter) to continue with next step"
		write(*,*) "Type 1 (enter) to specify the number of steps"

		write(*,*) "Type -1 (enter) to abort"
		write(*,*) "Type -2 (enter) to change delta value"
		write(*,*) "Type -3 (enter) to continue without extended output"
		read(*,*)izero
		select case(izero)
		case(0)
			isteps = 1
			go to 100
		case(1)
			write(*,*)' How many steps? (0 to abort)'
			read(*,*) isteps
			if(isteps.eq.0)go to 42
			go to 100
		case(-1)
		    	go to 9999
		case(-2)
			write(*,*)' Current monitor and delta values'
			do i = 1,variance
				write(*,*)i,mon(i),deltax(i)
				end do
			write(*,*)' Input new delta values for monitor parameters'
			do i = 1,variance
				write(*,*)i
				read(*,*)Deltax(i)
				end do
			go to 42
		case(-3)
			iLong(2) = 0
			iLong(3) = 0
			iLong(4) = 0
			iLong(5) = 0
			iLong(6) = 0
			go to 42
		case default
			go to 42
		end select

!		endif

!     DONE WITH THIS STEP
	go to 42

9999	continue		! exit routine
      	RETURN
      	END

