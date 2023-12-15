! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
!	routine to read the next record in file 87
!	if we hit EOF then rewind 87 and try again
	implicit none
	integer*4 i,numcurrent,logFileOn,iExclude,EOFflag,asmcurrent(*)
10	continue
	read(87,*,end=99)numCurrent,(asmCurrent(i),i=1,numCurrent)
	if(logFileOn.eq.1)then
		write(95,*)'Read new asm   ',(asmCurrent(i),i=1,numCurrent)
		endif
	if(iExclude.ne.0)then
		do 12 i=1,numCurrent
		if(iExclude.eq.asmCurrent(i))go to 10		! we need to skip this assemblage
12		continue
		endif
!	if we get here then iExclude is not in this new assemblage
	return
99	continue
	EOFflag = 1
	return
!	call RewindTheFile(87,logFileOn)
	!	rewind(87)
!	go to 10
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine RewindTheFile(iunit,logFileOn)
	implicit none
	integer*4 iunit,logFileOn
	rewind(iunit)
	if(logFileOn.eq.1)then
		write(95,*)'Rewind file = ',iunit
		endif
	return
	end		
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine CalculateGSystem
	implicit none
	include "Assemb.inc"
	include "Tangent.inc"
	integer*4 l

	gSystem = 0.0d0
	do 10 l = 1,nc
	gSystem = gSystem + tanPlane(l)*sysCoMoleFraction(l)	
10	continue
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine CalculateMode(sum)
	implicit none
	include "Assemb.inc"
	include "Tangent.inc"
	real*8 sum
	integer*4 kcur,k
!	sum = 0.0d0		! set to 0 in calling routine
! 	calculate mode
	Do 10 kCur = 1,numPh
	k = asmCurrent(kCur)
	sum = sum + vp0(k)
10	continue
!	write(*,*)' Mode -------- moles(k),volume(k),mode(k)'
	Do 20 kCur = 1,numPh
	k = asmCurrent(kCur)
	mode(k) = vp0(k)*100.d0/sum
!	write(*,44)k,phName(k),mp0(k),vp0(k),mode(k)
!44	format(I5,A12,3E20.10)
20	continue
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine ListTangentPlane
	implicit none
	include "Assemb.inc"
	include "Tangent.inc"
	integer*4 i
	write(12,*)'   '
	write(12,*)' Current tangent Plane ',TC,PB
	write(12,80)(coname(i),i=1,nc)
80	format(T10,20(A5,10x))
	write(12,81)(tanPlane(i),i=1,nc) !tangent plane is from system components
81	format(20F15.3)
	write(12,*)' Gsystem = ',gSystem
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine ListGDifference
	implicit none
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 k,j
	write(12,*)' '
	write(12,*)' Phase        G of Phase     G on Tangent  ignore G diff(GPh-GTan)       Phase composition'
	do 10 k = 1,numPhMIF
!	jj = minpointMIF(k)-1
	write(12,80)phName(K),gOfPhaseAU(K),gOnTanAU(K),ignorePhase(k),gDifferenceAU(k),(xPhCo(k,j),j=1,numPhCo(K))
80	format(A10,T12,2F15.2,I3,F15.2,15F12.5)
10	continue
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine ListCurrentAsm
	implicit none
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 k
!	write(12,*)' '
!	write(12,*)' Current assemblage:'
	write(12,101)(phName(asmCurrent(k)),k=1,numCurrent)
101	format(30(A8,2x))
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine CalculateSystemMoleFractions
	implicit none	
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
! --------------------------------------------------------------------
	integer*4 l
	real*8 sum
! ---------------------------
!	calculate system component mole fractions
	sum = 0.0d0
	do 10 l = 1,nc
	sum = sum + bulkCompMoles(l)
10	continue
	do 11 l = 1,nc
	sysCoMoleFraction(l) = bulkCompMoles(l)/sum
11	continue
	write(12,*)' '
	write(12,*)'Bulk composition = ', bulkCompTitle
	write(12,600)(coname(l),l=1,nc)
600	format(T12,24a10)
      write(12,602)(bulkCompWt(l),l=1,nc)
602   format(' Wt%',T12,24F10.3)
	write(12,601)(1000.*bulkCompMoles(l),l=1,nc)
601   format(' mMoles',T12,24F10.3)
	write(12,603)(100.*sysCoMoleFraction(l),l=1,nc)
603   format(' mole%',T12,24F10.3)
	
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine CalculateAtomUnits
!	Routine to calculate the conversion factors from molar units to atom units
!	All calculations of equilibria are done with molar units
!	The tangent plane is in atom units
!	The conversion factors are used in Subroutine ParallelToTangent
!		where the phases are compared to the tangent plane
!	and in Subroutine AdjustTangentToAsm
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
! --------------------------------------------------------------------
	integer*4 k,j,i
	real*8 sum
! ************************************************************************
	do 5 k = 1,numPhMIF
!	jj = minpointMIF(K)-1
	do 10 j = 1,numPhCo(K)
	sum = 0.0d0
!	jjj = jj + j
	do 11 i = 1,nc
	sum = sum + comp(k,j,i)
11	continue
	atomNorm(k,j) = sum		! there is one conversion factor for each phase component
10	continue
5	continue
	return
	end

! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE LoadCurrentAsm(TC4,PB4,izero,isubtend)
      implicit none
! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
! --------------------------------------------------------------------
	include "Tangent.inc"
! --------------------------------------------------------------------
!     local variables
      	integer*4 k,izero,isubtend,kCur,i,L,j
      	REAL*4 TC4,PB4
!	This code should be exactly similar to the exit code in Subroutine Change
!	put the phases into arrays for computation
	np = 0
	do 20 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	np = np + numPhCo(k)
!	by commenting these out, we should be using the moles from the previous solution
!	mp0(k) = 0.1		! reset the moles of phases to .1
!	mp1(k) = .1
!	mp2(k) = .1

20	continue
	numPh = numCurrent
!	xPhCo = xPhCoInitial 	! set all compositions back to initial (input) values
	NX=NP-numPh		!	compute number of independent compositional (dX) variables (before SetAllX)
	TandP = 2
	numNew = 0

!	Set ALLX array
      	DO 40 I=1,6
      	ALLX(I,1)=TC
      	ALLX(I,2)=PB
!     	THESE ARE THE X'S
      	L=TANDP
	Do 41 kCur = 1,numPh
	k = asmCurrent(kCur)
      	DO 42 J=2,numPhCo(K)
      	L=L+1
      	ALLX(I,L)=xPhCo(k,j)
42   	CONTINUE
41	CONTINUE
!     THESE ARE MOLES
      	L=TANDP+NX
	Do 44 kCur = 1,numPh
	k = asmCurrent(kCur)
	L=L+1
	ALLX(I,L)=MP0(K)
44	CONTINUE
40   	CONTINUE
!	done setting ALLX array
	call NAMES()
	ipseudo = 1	! supresses output in REXN
	idegen = 0	! returns 1 if assemblage is degenerate in REXN
	call REXN()
	if(idegen.eq.1)then	! assemblage is degenerate
		isubtend = 0
!		if(ioutputPseudo.eq.1)write(12,*)' Assemblage is degenerate'
		return
		endif
!	nvar = TandP + np	!	nvar = TandP + np - nPh + nPh
	NVAR= TANDP + NX	!	compute total number of variables
	NVAR = NVAR + numPh
! 	Here we need to check to be sure the assemblage subtends the bulk composition
	NEQ = NRX + NC		! total number of equations
	isubtend = 1		! = yes (where we start)
	call subtend(isubtend)
	if(isubtend.eq.0)then
!		write(12,*)' Assemblage does not subtend bulk composition'
!		pause 'Hit return to continue'
		isubtend = 0		! no, it doesn't subtend bulk composition
		return
		endif
	return
	end

! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE LoadCurrentAsmNoMassBal(TC4,PB4,izero,isubtend)
      implicit none
! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
! --------------------------------------------------------------------
	include "Tangent.inc"
! --------------------------------------------------------------------
!     local variables
      	integer*4 k,izero,kCur,isubtend
      	REAL*4 TC4,PB4

!	put the phases into arrays for computation
	np = 0
	do 20 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	np = np + numPhCo(k)
20	continue
	numPh = numCurrent
	NX=NP-numPh
	call SetAllX()
	call NAMES()
	ipseudo = 1
	idegen = 0
        CALL REXN()
	if(idegen.eq.1)then	! assemblage is degenerate
		isubtend = 0
!		if(ioutputPseudo.eq.1)write(12,*)' Assemblage is degenerate'
		return
		endif
 	TandP = 2
	NVAR= TANDP + NX	!	compute total number of variables
	NEQ = NRX		! total number of equations
     	TC=TC4
      	Pb=Pb4
      	ALLX(1,1)=TC
      	ALLX(1,2)=PB
      	ALLX(2,1)=TC
      	ALLX(2,2)=PB
        call names
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	Subroutine AdjustALLX_Array()
!	Routine to set the ALLX array after loading a new assemblage from the MIF
!	This is necessary because it is the ALLX array that is changed during the COMPUTE cycle
!	Note: this code is copied from Subroutine BEGIN (lines 487-545 of file Gibbs_FileInput.f)
!	implicit none
! *****************************************
!	include "Assemb.inc"
!	include "Monit.inc"
!	include "Tangent.inc"
!	include "MIFarrays.inc"
!	integer*4 i,k,j,l,kCur

!        ARRAY ALLX CONTAINS VALUES OF THE FOLLOWING
!     set all columns of ALLX(i,j) to starting values
!     Definition of ALLX variable:
!     1,50   current value
!     2,50     starting value
!     3,50   ref at start of contour
!     4,50   user selected reference
!     5,50   previous finite diff point
!     6,50     (not used)

!      	DO 444 I=1,6
!      	ALLX(I,1)=TSTART
!      	ALLX(I,2)=PSTART
!      	IF (KFLU.NE.0)ALLX(I,3)=PFStart
!     	THESE ARE THE X'S
!      	L=TANDP
!      	DO 441 kCur =1,numPh
!	k = asmCurrent(kCur)
!      	DO 555 J=2,numPhCo(K)
!      	L=L+1
!      	ALLX(I,L)=xPhCo(k,j)
!555   	CONTINUE
!441   	CONTINUE
!     	THESE ARE MOLES
!      	L=TANDP+NX
!      	IF(IMASS.EQ.1) THEN
!            DO 442 K=1,numPh
!            L=L+1
!            ALLX(I,L)=MP0(K)
!            if(fractl(k).eq.1.or.fractl(k).eq.2)then
!            	!This is done as a flag for COMPUT2 where it checks for Newton convergence (around line 722)
!            	! If FRACTL = 1 or 2 then NEWTON doesn't converge properly.
!            	! This flag will omit moles of fractionating phases from Newton convergence
!            	allx(6,L)=1
!            	else
!            	allx(6,L)=0
!            	endif
!442         CONTINUE
!            ENDIF


!     Set allFractl switch if necessary
!      IF(IMASS.EQ.1) THEN
!      L=TANDP+NX
!            DO 443 K=1,numPh
!            L=L+1
!            if(fractl(k).eq.1.or.fractl(k).eq.2)then
            	!This is done as a flag for COMPUT2 where it checks for Newton convergence (around line 722)
            	! If FRACTL = 1 or 2 then NEWTON doesn't converge properly.
            	! This flag will omit moles of fractionating phases from Newton convergence
!            	allFractl(L)=1
!            	else
!            	allFractl(L)=0
!            	endif
!443         CONTINUE
!        	ENDIF
!444   CONTINUE

!	return
!	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	Subroutine LoadXfromMIF
!	Routine to load compositions of phases from array xMIF into working array x
!	implicit none
! *****************************************
!	include "Assemb.inc"
!	include "Tangent.inc"
!	include "MIFarrays.inc"
!	integer*4 k,j,jj,jjj
! --------------------------------------------------------------------
! Load the mineral compositions from the array xMIF into the array x
!	do 10 k = 1,nph
!	jj = minPoint(K) - 1
!	do 20 j = 1,noPhCo(k)
!	jjj = jj+j
!	x(jjj) = xMIF(jjj)
!20	continue
!10	continue
!	return
!	end	



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	Subroutine StoreXinMIF()
!	Routine to load compositions of phases from array x into MIF array xMIF
!	implicit none
! *****************************************
!	include "Assemb.inc"
!	include "Tangent.inc"
!	include "MIFarrays.inc"
!	integer*4 k,j,jjMIF,jj,kMIF
! --------------------------------------------------------------------
! Load the mineral compositions from the array xMIF into the array x
!	do 10 k = 1,nph	! loop through the phases in this assemblage only
				! note numCurrent = nph
!	kMIF = asmCurrent(k)	! index of phase(k) in MIF
!	molesPhaseMIF(kMIF) = MP0(k)	! store the moles of the phase in MIF file
!	jjMIF = minPointMIF(kMIF)-1
!	jj = minPoint(K) - 1
!	do 20 j = 1,noPhCo(k)
!	xMIF(jjMIF + j) = x(jj+j)
!20	continue
!10	continue
!	return
!	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine OpenOutputFile(status)
	implicit none
	integer*4 status
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
! --------------------------------------------------------------------
!	open(73,FILE='',status='UNKNOWN',iostat=status)
	open(73,FILE='',status='NEW',iostat=status,action='WRITE')
	if(status.ne.0)return
	inquire(73,NAME=tangentOutputFile)
	write(73,*)' This file is intentionally left empty'
	write(73,*)' Feel free to delete it'
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteOutputHeader(iunit)
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 i,iunit
! --------------------------------------------------------------------

! write some stuff to the header of file = iunit = output file
	write(iunit,*)'Bulk composition = ', bulkCompTitle
	write(iunit,*)' Bulk composition - input wt% and moles'
	write(iunit,600)(coname(i),i=1,nc)
600	format(10x,24a8)
	write(iunit,601)(1000.*bulkCompMoles(i),i=1,nc)
601	format(' mMoles   ',24F8.3)
	write(iunit,602)(bulkCompWt(i),i=1,nc)
602	format(' Wt %     ',24F8.3)
	write(iunit,*)' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'				! extra line to preserve compatibility with PseudoForwardModel routine
	write(iunit,*)'***************************************'
	write(iunit,*)'MasterInputFile = ',FileIn
!	write(iunit,7100)noStable
!7100	format(I8,4x,'Number of different stable mineral assemblages')	
	! I should list all considered minerals here
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteGoodStuff			! formerly ListGoodStuff
	implicit none
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 k,j
	write(72,*)TC,Pb
	write(72,*)'         Phase          mMoles phase       G diff         Phase composition'
	do 10 k = 1,numPhMIF
!	jj = minpointMIF(k)-1
	write(72,80)k,minRec(K),phName(K),mp0(k)*1000.,gDifferenceAU(k),(xPhCo(k,j),j=1,numPhCo(K))
80	format(I5,I5,2x,A10,T25,F12.3,F15.3,15F12.5)
10	continue
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteOutputAsm
	implicit none
	integer*4 k,j,i
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
! --------------------------------------------------------------------
	write(73,80)TC,PB,gSystem,numCurrent
80	format(2F10.1,F15.2,I5,' = numPhases')
	do 10 i = 1,numCurrent
	k = asmCurrent(i)
	write(73,82)minRec(k),phName(k),mp0(k),vP0(k),numPhCo(k),(phCoName(k,j),xPhCo(k,j),j=1,numPhCo(k))
82	format(T40,I5,2x,A8,F15.5,4x,F15.5,4x,I5,20(5x,A8,F12.5))
10	continue	
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	subroutine WriteAllOutHeader(iunit,Tinc,Pinc)
	subroutine WriteAllOutHeader(iunit)
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 i,iunit
!	real*4 Tinc,Pinc
! --------------------------------------------------------------------

! write some stuff to the header of file = iunit = output file
	write(iunit,*)'Bulk composition = ', bulkCompTitle
	write(iunit,*)' Bulk composition - input wt% and moles'
	write(iunit,600)(coname(i),i=1,nc)
600	format(10x,24a8)
	write(iunit,601)(1000.*bulkCompMoles(i),i=1,nc)
601	format(' mMoles   ',24F8.3)
	write(iunit,602)(bulkCompWt(i),i=1,nc)
602	format(' Wt %     ',24F8.3)
	write(iunit,*)' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'				! extra line to preserve compatibility with PseudoForwardModel routine
	write(iunit,*)'***************************************'
	write(iunit,*)'MasterInputFile = ',FileIn
!	write(iunit,7100)noStable
!7100	format(I8,4x,'Number of different stable mineral assemblages')	
	! I should list all considered minerals here
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	write(iunit,*)Tinc,Pinc,'         Tincrement   Pincrement in file'
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	write(iunit,*)numPhMIF,'           Total number of phases in MIF'
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteAllStuff			!Routine to write out information to unit 74
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
! --------------------------------------------------------------------
	integer*4 k,j
	real*4 xPhCoTemp(phCoMax)
! --------------------------------------------------------------------
	write(74,80)TC,PB,gSystem,numCurrent
80	format(2F10.1,F20.5,I5,' = numPhases')
	write(74,*)'                 Phase           GDiff         mMoles phase           V phase              Phase composition'
	do 10 k = 1,numPhMIF
	if(Dabs(gDifferenceAU(k)).gt.0.1d0)then
		vP0(k) = 0.0d0
		mp0(k) = 0.0d0
		endif
	do 20 j = 1,numPhCo(k)
	if(Dabs(gDifferenceAU(k)).eq.1.0d4)then
		xPhCoTemp(j) = 0.0		! we zero these for output when there is no solution to the paralleltotangent
		else
		xPhCoTemp(j) = xPhCo(k,j)
		endif
20	continue
!	write(74,82)minRec(k),phName(k),gDifferenceAU(k),mp0(k)*1000,vP0(k),numPhCo(k),(phCoName(k,j),xPhCo(k,j),j=1,numPhCo(k))
	write(74,82)minRec(k),phName(k),gDifferenceAU(k),mp0(k)*1000.,vP0(k),numPhCo(k),(phCoName(k,j),xPhCoTemp(j),j=1,numPhCo(k))
!82	format(T12,I5,2x,A8,F15.5,4x,F15.5,4x,F15.5,4x,I5,20(5x,A8,E12.5))
82	format(T12,I5,2x,A8,F15.5,4x,E15.5,4x,E15.5,4x,I5,20(5x,A8,E15.5))
10	continue	
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteAllStuffToTerminal(iunit)		!Routine to write out information to unit 74
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
! --------------------------------------------------------------------
	integer*4 k,j,iunit
	real*4 xPhCoTemp(phCoMax)
! --------------------------------------------------------------------
	write(iunit,80)TC,PB,gSystem,numCurrent
80	format(2F10.1,F20.5,I5,' = numPhases')
	write(iunit,*)'                 Phase           GDiff         mMoles phase           V phase              Phase composition'
	do 10 k = 1,numPhMIF
	if(Dabs(gDifferenceAU(k)).gt.0.1d0)then
		vP0(k) = 0.0d0
		mp0(k) = 0.0d0
		endif
	do 20 j = 1,numPhCo(k)
	if(Dabs(gDifferenceAU(k)).eq.1.0d4)then
		xPhCoTemp(j) = 0.0		! we zero these for output when there is no solution to the paralleltotangent
		else
		xPhCoTemp(j) = xPhCo(k,j)
		endif
20	continue
!	write(74,82)minRec(k),phName(k),gDifferenceAU(k),mp0(k)*1000,vP0(k),numPhCo(k),(phCoName(k,j),xPhCo(k,j),j=1,numPhCo(k))
	write(iunit,82)k,minRec(k),phName(k),gDifferenceAU(k),mp0(k)*1000.,vP0(k),numPhCo(k),	&
			(phCoName(k,j),xPhCoTemp(j),j=1,numPhCo(k))
!82	format(T12,I5,2x,A8,F15.5,4x,F15.5,4x,F15.5,4x,I5,20(5x,A8,E12.5))
82	format(T5,I5,2x,I5,2x,A8,F15.5,4x,E15.5,4x,E15.5,4x,I5,20(5x,A8,E15.5))
10	continue	
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteAllStuffToLogFile()		!Routine to write out information to unit 74
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
! --------------------------------------------------------------------
	integer*4 k,j
	real*4 xPhCoTemp(phCoMax)
! --------------------------------------------------------------------
	write(95,80)TC,PB,gSystem,numCurrent
80	format(2F10.1,F20.5,I5,' = numPhases')
	write(95,*)'                 Phase           GDiff         mMoles phase           V phase              Phase composition'
	do 10 k = 1,numPhMIF
	if(Dabs(gDifferenceAU(k)).gt.0.1d0)then
		vP0(k) = 0.0d0
		mp0(k) = 0.0d0
		endif
	do 20 j = 1,numPhCo(k)
	if(Dabs(gDifferenceAU(k)).eq.1.0d4)then
		xPhCoTemp(j) = 0.0		! we zero these for output when there is no solution to the paralleltotangent
		else
		xPhCoTemp(j) = xPhCo(k,j)
		endif
20	continue
!	write(74,82)minRec(k),phName(k),gDifferenceAU(k),mp0(k)*1000,vP0(k),numPhCo(k),(phCoName(k,j),xPhCo(k,j),j=1,numPhCo(k))
	write(95,82)minRec(k),phName(k),gDifferenceAU(k),mp0(k)*1000.,vP0(k),numPhCo(k),(phCoName(k,j),xPhCoTemp(j),j=1,numPhCo(k))
!82	format(T12,I5,2x,A8,F15.5,4x,F15.5,4x,F15.5,4x,I5,20(5x,A8,E12.5))
82	format(T12,I5,2x,A8,F15.5,4x,E15.5,4x,E15.5,4x,I5,20(5x,A8,E15.5))
10	continue	
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteBulkCompHeader(iunit,itype)
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 i,iunit,itype
!	real*4 Tinc,Pinc
! --------------------------------------------------------------------

! write some stuff to the header of file = iunit = output file
	write(iunit,*)'Bulk composition = ', bulkCompTitle
	write(iunit,*)' Bulk composition - input wt% and moles'
	write(iunit,600)(coname(i),i=1,nc)
600	format(10x,24a8)
	write(iunit,601)(1000.*bulkCompMoles(i),i=1,nc)
601	format(' mMoles   ',24F8.3)
	write(iunit,602)(bulkCompWt(i),i=1,nc)
602	format(' Wt %     ',24F8.3)
	write(iunit,*)' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'				! extra line to preserve compatibility with PseudoForwardModel routine
	write(iunit,*)'***************************************'
	write(iunit,*)'MasterInputFile = ',FileIn
!	write(iunit,7100)noStable
!7100	format(I8,4x,'Number of different stable mineral assemblages')	
	! I should list all considered minerals here
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	write(iunit,*)Tinc,Pinc,'         Tincrement   Pincrement in file'
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	write(iunit,*)nc,'        Number of system components'
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	select case(itype)
	case(1)
!	Writing out the bulk composition wt
		write(iunit,603)(coname(i),i=1,nc)
603		format('       T        P       BC_m_moles    rockVolume         ',24(a8,7x))
	case(2)
!	bulk composition moles
		write(iunit,604)(coname(i),i=1,nc)
604		format('       T        P         BCmass      rockVolume         ',24(a8,7x))
	case(3)
!	Melt composition
		write(iunit,605)(coname(i),i=1,nc)
605		format('       T        P      Melt_m_moles     MeltVolume         ',24(a8,7x))
	case default
	end select
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteBCMoles			!Routine to write out information to unit 74
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
! --------------------------------------------------------------------
	integer*4 i
! --------------------------------------------------------------------
	write(75,80)TC,PB,BCmoles,BCvol,(1000.*bulkCompMoles(i),i=1,nc)
80	format(2F10.1,24F15.8)
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteBCWt			!Routine to write out information to unit 74
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
! --------------------------------------------------------------------
	integer*4 i
! --------------------------------------------------------------------
	write(76,80)TC,PB,BCmass,BCvol,(bulkCompWt(i),i=1,nc)
80	format(2F10.1,24F15.8)
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteMeltWt			!Routine to write out information to unit 74
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
! --------------------------------------------------------------------
	integer*4 i
! --------------------------------------------------------------------
	write(78,80)TC,PB,MeltMoles,MeltVol,(MeltWt(i),i=1,nc)
80	format(2F10.1,24F15.8)
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	subroutine WriteAllOutHeader(iunit,Tinc,Pinc)
	subroutine WriteFractHeader(iunit)
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Tangent.inc"
!	include "MIFarrays.inc"
	integer*4 i,iunit
!	real*4 Tinc,Pinc
! --------------------------------------------------------------------

! write some stuff to the header of file = iunit = output file
	write(iunit,*)'Bulk composition = ', bulkCompTitle
	write(iunit,*)' Bulk composition - input wt% and moles'
	write(iunit,600)(coname(i),i=1,nc)
600	format(10x,24a8)
	write(iunit,601)(1000.*bulkCompMoles(i),i=1,nc)
601	format(' mMoles   ',24F8.3)
	write(iunit,602)(bulkCompWt(i),i=1,nc)
602	format(' Wt %     ',24F8.3)
	write(iunit,*)' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'				! extra line to preserve compatibility with PseudoForwardModel routine
	write(iunit,*)'***************************************'
	write(iunit,*)'MasterInputFile = ',FileIn
!	write(iunit,7100)noStable
!7100	format(I8,4x,'Number of different stable mineral assemblages')	
	! I should list all considered minerals here
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	write(iunit,*)Tinc,Pinc,'         Tincrement   Pincrement in file'
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	write(iunit,*)numFractionate,'           Total number of phases that are fractionating'
	write(iunit,*)' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine WriteFractionate			!Routine to write out information to unit 74
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
! --------------------------------------------------------------------
	integer*4 k,j,iFrac
! --------------------------------------------------------------------
	write(77,80)TC,PB,gSystem,numFractionate
80	format(2F10.1,F20.5,I5,' = numFractionating Phases')
	write(77,*)'                 Phase           GDiff         mMoles removed         Vol removed           Phase composition'
	do 10 iFrac = 1,numFractionate
	k = toFractionate(iFrac)
	write(77,82)minRec(k),phName(k),gDifferenceAU(k),removeMoles(iFrac)*1000.,removeVol(iFrac),	&
			numPhCo(k),(phCoName(k,j),xPhCo(k,j),j=1,numPhCo(k))
!	write(74,82)minRec(k),phName(k),gDifferenceAU(k),mp0(k)*1000.,vP0(k),numPhCo(k),(phCoName(k,j),xPhCoTemp(j),j=1,numPhCo(k))
82	format(T12,I5,2x,A8,F15.5,4x,F15.5,4x,F15.5,4x,I5,20(5x,A8,E15.5))
	
10	continue
	return
	end
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine SetUpFractionation()
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Tangent.inc"
! --------------------------------------------------------------------
	integer*4 i
! --------------------------------------------------------------------
	numFractionate = 0
	do 4411 i = 1,numPhMIF
4411	write(*,*)i,phName(i)
4412	continue
	write(*,*)' Input the number of the phase you want to fractionate.'
	write(*,*)'   Enter 0 to end the list'
	read(*,*)i
	if(i.gt.0)then
		numFractionate = numFractionate + 1
		toFractionate(numFractionate) = i
		write(*,*)' Input the residual amount (lower bound) of the phase to remain in the rock'
		write(*,*)'   Units are volume % (mode) that will not be removed'
		read(*,*)residual(numFractionate)
		write(*,*)' Input threshold (upper bound) to remove phase'
		write(*,*)'   This is the volume % (mode) at which the phase will be removed'
		write(*,*)' To maintain a constant amount (e.g. porosity) make the threshold the same as the residual'
		read(*,*)threshold(numFractionate)
		write(*,*)numFractionate,toFractionate(numFractionate),residual(numFractionate),threshold(numFractionate)
		go to 4412
		endif

	return	
	end
	
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine MakeParserFile
!	Routine to make a file that can be read by program MAD_Parser
!	The file will contain the GFWs and other information about the system components
!	as well as the composition of each phase component in each phase of a MIF in terms of the systme components
	implicit none
! ****************************************
	include "Assemb.inc"
! ****************************************
!     local variables
	integer*4 i,j,k
! *******************************************************************
!	Open the parser file for output
	open(23,file='',status='NEW')

!	write file header
	write(23,*)nc		! number of system components
	write(23,105)(coname(i),i=1,nc)
105	format(' Names',5x,25(A4,8x))
	write(23,106)(molwt(i),i=1,nc)
106	format(' GFWs ',5x,25(F12.3))
	write(23,107)(numCatInOxide(i),i=1,nc)
107	format(' Cats ',5x,25(F12.3))
	write(23,108)(OxToElWt(i),i=1,nc)
108	format(' ElOx ',5x,25(F12.5))
!	write out information for each phase and phase component
	Do 41 k = 1,numPhMIF
	write(23,101)minRec(k),phName(k),numPhCo(k)
101	format(I5,2x,A8,I5)
	do 42 j=1,numPhCo(K)
	write(23,102)PhCoName(k,j),(comp(k,j,i),i=1,nc)
102	format(T5,A8,5x,15(F12.3))
42   	continue
41   	continue
	
	close(23)
	return
	end
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Subtend (iyes)
      	implicit none
!	routine to determine if an assemblage subtends the bulk composition
!	Actually, whether all system components are in the collection of phase components
! c*****************************************
	include "Assemb.inc"
	include "Tangent.inc"
! c*****************************************
	integer*4 j,l,iyes,k,kCur
	real*8 sum
!	Does this assemblage subtend the bulk composition?
!	several criteria must be met
!	FIrst, does the chemical variability of the phases span the bulk composition?
!	That is, are the number of phase components at least as large as the number of system components
!	e.g. System Al2O3-SiO3 has 2 components. Al2SiO5 has all the components, but no chemical variability, so it does not subtend
	j = 0
	do 20 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	j = j + numPhCo(k)
20	continue
	if(j.lt.nc)then
		iyes = 0	!  0 = NO
		return
		endif
! 	Add up colums for each phase component
! 	If a column = 0 then this system component isn't in the phase assemblage
!	put the phases into arrays for computation
        DO 10 L=1,NC		! loop for each column = a system component
	sum = 0.0d0
	do 5 kCur = 1,numCurrent	! loop on all phases
	k = asmCurrent(kCur)
      	DO 12 j=1,numPhCo(k)		! loop on all phase components in the phase
            sum = sum + COMP(k,j,L)
12	continue
5	continue
	if(sum.lt.1.e-8)then			! sum is for a single element - if 0 then there ain't enough
		iyes = 0			!  0 means NO, it does not subtend
		return
		endif
10    	CONTINUE
	iyes = 1			! means yes it does subtend
	return
	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine ChooseConstantPorosity()
	implicit none
	integer*4 i
	include "Assemb.inc"
	include "Tangent.inc"
! *****************************************
	write(*,*)' Do you want to keep porosity for fluid constant?'
	write(*,*)' 0 = no '
	write(*,*)' 1 = yes'
	read(*,*)constantPorosity
	if(constantPorosity.eq.1)then
		do 250 i = 1,numPhMIF
250		write(*,*)i,phName(i)
		write(*,*)'Which phase is the fluid?'
		read(*,*)fluidIndex
		write(*,*)'What value of porosity do you want (in %)'
		read(*,*)porosity
		endif
	return
	end



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine ChooseInitialAssemblage()
	implicit none
	integer*4 i
	include "Assemb.inc"
	include "Tangent.inc"
! *****************************************

	open(87,file='Fort_87_ComboScratch',status = 'UNKNOWN')
	iParagonite = 0
	do 228 i = 1,numPhMIF
	if(minRec(i).eq.118)then
		iParagonite = i			! this is a switch so we can initialize Xparag to a large value
		endif
228	write(*,*)i,phName(i)
	write(*,*)'Please provide your best guess for the initial assemblage.'
	write(*,*)'Input numbers for these phases (0 to end list)'
	numCurrent = 0
226	continue
	read(*,*)i
	if(i.eq.0)then
		go to 229
		else
		numCurrent = numCurrent + 1
		asmCurrent(numCurrent) = i
		endif
	go to 226
229	continue
	if(numCurrent.gt.0)then		! only write out if an assemblage was specified
		write(87,*)numCurrent,(asmCurrent(i),i=1,numCurrent)
		endif
	write(*,*)'iParagonite = ',iParagonite
!	iall = 0
!		now sort through and reorder possible assemblages with those with likely phases first
	rewind(85)
!240		continue
!		read(85,*,end=248)numCurrent,(asmCurrent(i),i=1,numCurrent)
!		do 244 j = 1,iAll
!		do 242 i = 1,numCurrent
!		if(asmCurrent(i).eq.iAllAsm(j))go to 244		! we have a match
!242		continue
!		if we get here, then we do not have a match for this iAllAsm(j)
!		write to unit 87 then go get the next assemblage
!		write(86,*)numCurrent,(asmCurrent(i),i=1,numCurrent)
!		go to 240	
!244		continue
!		if we get here then iAll phases match. Write this out
!		write(87,*)numCurrent,(asmCurrent(i),i=1,numCurrent)
!		go to 240
!248		continue			! we looked through entire file unit = 84
!		write other posible assemblages to unit 85
!		rewind(86)
246	continue
	read(85,*,end=249)numCurrent,(asmCurrent(i),i=1,numCurrent)
	write(87,*)numCurrent,(asmCurrent(i),i=1,numCurrent)
	go to 246
249	continue
	close(84)	! temp file with all possible assemblages
	close(85)	! temp file with subtended assemblages
!	close(86)	! temp file with "preferred" assemblages
	rewind(87)	! temp file with preferred assemblages followed by the rest
!	write(12,*)'=========Initial setup done============================='
	write(*,*)'=========Initial setup done============================='

	return
	end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine ZeroTanOut()
	implicit none
	integer*4 i
	include "Assemb.inc"
	include "Tangent.inc"
! *****************************************
	do 3 i = 1,20
	tanOut(i) = 0
3	continue
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine SetTanOut()
	implicit none
	integer*4 i
	include "Assemb.inc"
	include "Tangent.inc"
! *****************************************
1	continue
!	do 3 i = 1,7
!	write(*,*)i,tanOut(i)
	write(*,*)'1 ',tanOut(1),'  Output results of ParallelToTangent (summary)'
	write(*,*)'2 ',tanOut(2),'  Output every iteration of ParallelToTangent'
	write(*,*)'3 ',tanOut(3),'  Pause after every phase ParallelToTangent'
	write(*,*)'4 ',tanOut(4),'  Call to ThermoData with new assemblage'
	write(*,*)'5 ',tanOut(5),'  Output transformation matrix (B) - phase to system components'
	write(*,*)'6 ',tanOut(6),'  Call to PRINTT'
	write(*,*)'7 ',tanOut(7),'  Call ListTangentPlane'
	write(*,*)'8 ',tanOut(8),'  Call DumpEverything()'
!3	continue
	write(*,*)'Input number to change, 0 to exit'
	read(*,*)i
	if(i.eq.0)return
	if(tanOut(i).eq.1)then
		tanOut(i) = 0
		else
		tanOut(i) = 1
		endif	
	go to 1	
	end
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SetInitialTangent
      implicit none
!     ROUTINE TO find an initial tangent plane that subtends the bulk composition
!	given a system of phase components and a chemical system
!	The algorithm does this by computing a set of linearly independent reactions among
!	the phase components and the apexes of the tangent plane (e.g. SiO2, MgO, AlO3/2, etc)
!	This code uses the same approach as in Subroutine REXN except that a "fictive" tangent phase is
!	added to the right side of the input matrix. After the Gauss reduction, this results in each
!	"component" of the tangent "phase" being defined in terms of the other phase components

! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
! c*****************************************
! --------------------------------------------------------------------
	include "Tangent.inc"
! --------------------------------------------------------------------

!     Local arrays
!     DIMENSIONED FOR 12 SYSTEM COMPONENTS AND 25 PHASE COMPONENTS
      INTEGER*4 L,j,jj,k,i,ll,kCur,npTemp,npTotal
      REAL*8 A(100,100),TEMP,FACT
      common /REXNcommon/A,Temp,Fact
      	Character VN3(100)*4
	common/VN3names/VN3
! c*****************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	This part is exactly like the REXN routine
!     Set up array of phase component names
	L=0
	npTemp = 0
	jj = 0
	Do 10 kCur = 1,numPh
	k = asmCurrent(kCur)
	DO 10 J=1,numPhCo(k)
	jj = jj + 1
        VN3(jj)=phCoName(k,J)
	npTemp = npTemp + 1
10  CONTINUE
!	npMIF = npTemp			! total number of phase components in the MIF
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Add on "components" for the tangent phase
	l = npTemp
	do 15 i = 1,nc
	l = l + 1
	vn3(l) = coname(i)
15	continue
!++++++++++++++++++++++++++++++++++++++++++++++++++
!     THIS Routine IS DESIGNED TO COMPUTE A LINEARLY INDEPENDENT
!     SET OF REACTIONS AMONG PHASE COMPONENTS IN A SYSTEM OF
!     NC SYSTEM COMPONENTS AND NP PHASE COMPONENTS
!     USES METHOD OF GAUSS-JORDAN REDUCTION OPERATING ON AN IDENTITY
!     MATRIX

      IF (NP.GT.50) THEN
         WRITE(*,*) 'Dimension of array ARX in subroutine REXN is too small -----stop'
         pause 'Hit return to continue'
         stop
      endif
	L = 0

!	asmCurrent should include all phases in the input file
!		although it could just include a set of phases that will subtend the bulk composition
	Do 20 kCur = 1,numPh
	k = asmCurrent(kCur)
	do 20 j = 1,numPhCo(k)
	L = L+1
        DO 20 ll=1,NC
        A(L,ll)=COMP(k,j,ll)
20    CONTINUE
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Add equations for the tangent phase
	L = npTemp			! total number of phase components
	do 25 i = 1,nc		! do each row
	L = L + 1
	do 26 j = 1,nc		! do each column
	A(l,j) = 0.0d0
26	continue
	A(l,i) = 1.0d0
25	continue
	npTotal = npTemp + nc
!++++++++++++++++++++++++++++++++++++++++++++++++++
!     ECHO INPUT DATA
      IF(iLong(7).eq.1)then
		write(12,*)
		write(12,*)
		write(12,*)'INPUT DATA'
		write(12,*)
		write(12,*)
		write(12,2000)(CONAME(I),I=1,NC)
2000    	FORMAT(' ',10X,20(4X,A4))
		DO 30 I=1,npTotal
		write(12,2001)VN3(I),(A(I,J),J=1,NC)
2001    	FORMAT(' ',A8,2X,20(F8.3))
30      	CONTINUE
		endif

	DO 50 J=1,npTotal
	DO 50 I=1,npTotal
50	ARX(I,J)=0.
	DO 51 J=1,npTotal
51	ARX(J,J)=1.

      	DO 100 J=1,NC
	JJ=J+1
	IF(DABS(A(J,J)).GT..00001)GO TO 115
	DO 110 K=J+1,npTotal
110     IF(DABS(A(K,J)).GT.0.00001)GO TO 120

! 	If here, then the assemblage is degenerate. Sometimes this does not cause a problem
! 	but sometimes it does.
! 	If we are doing a "PseudoFromScratch", then we want to avoid these degeneracies.
	If (ipseudo.eq.1)then
		idegen = 1
		return
		endif
         WRITE(*,*)' WARNING----'
         WRITE(*,*)' The number of system components specified is'
         WRITE(*,*)' greater than the number required to describe'
         WRITE(*,*)' the chemical variability of the phases'
         WRITE(*,*)' (the last column of array A in SUBROUTINE REXN'
         WRITE(*,*)' is all zero. COL=',J
         WRITE(*,*)' NOTE: usually, this does not cause a problem'
         write(*,*) 'Hit return to continue...'
         read(*,*)
         GO TO 100
120     CONTINUE
!        SWITCH ROW K WITH ROW J
      	DO 125 I=1,NC
      	TEMP=A(J,I)
      	A(J,I)=A(K,I)
125    	A(K,I)=TEMP
      	DO 130 I=1,npTotal
      	TEMP=ARX(J,I)
	ARX(J,I)=ARX(K,I)
130   	ARX(K,I)=TEMP
115    CONTINUE
      	DO 140 I=JJ,npTotal
	FACT=A(I,J)/A(J,J)
	DO 145 K=1,npTotal
145	ARX(I,K)=ARX(I,K)-ARX(J,K)*FACT
	DO 150 K=1,NC
150	A(I,K)=A(I,K)-A(J,K)*FACT
140     CONTINUE
100   	CONTINUE

!     ALL DONE WITH ROW REDUCTION
!     CHECK TO SEE HOW MANY ROWS ARE ALL ZERO--START FROM BOTTOM UP
      DO 160 I=npTotal,1,-1
         DO 160 J=1,NC
            IF(DABS(A(I,J)).GT.0.0001)GO TO 410
160  CONTINUE
410   CONTINUE

!     	THE ITH ROW HAS NON ZERO ELEMENTS
!     	CALCULATE # INDEPENDENT REACTIONS=NRX
      	NRX = npTotal - I
!     	CALCULATE TOTAL # EQUATIONS=NEQ
!      	NEQ = NRX
!      	IF(IMASS.EQ.1)NEQ = NEQ + NC
!      	if(numNew.gt.0)NEQ = NEQ + numNew
!      	if(PFluidSwitch.eq.2)NEQ = NEQ + 1

!     	OUTPUT SECTION

      IF(iLong(7).eq.1)then
420      CONTINUE
         write(12,*)
         write(12,*)'NC=',NC,'  npTotal=',npTotal,'  # REACTIONS=',NRX
         write(12,*)
         write(12,*)'REDUCED COMPOSITION MATRIX'
         write(12,2100)(CONAME(I),I=1,NC)
 2100    FORMAT(' ',10X,30(4X,A4))
         DO 449 I=1,npTotal
         write(12,2101)VN3(I),(A(I,J),J=1,NC)
2101     FORMAT(' ',A8,2X,30(F8.3))
449      CONTINUE
         write(12,*)
         write(12,*)
         write(12,*)'REDUCED PHASE MATRIX'
         write(12,2110)(VN3(J),J=1,npTotal)
2110     FORMAT(' ',4X,30(A4,4X))
         DO 429 I=1,npTotal
         write(12,2120)(ARX(I,J),J=1,npTotal)
         IF(I.eq.npTotal-NRX)then
           write(12,*)' '
           write(12,*)' Linearly independent reactions:'
           write(12,2110)(VN3(J),J=1,npTotal)
!           write(12,*)' '
           endif
2120     FORMAT(' ',30F8.3)
429      CONTINUE
	endif      

      return
      end
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      	SUBROUTINE BuildJacobian(izero)
 	use MatrixArrays
     	implicit none
!     	ROUTINE TO build the matrix for calculation of the Jacobian for the chemical potentials
!	The matrix is set up as follows:
!	Row 1: definition of the dG for the bulk composition
!	Rows 2 to numph+1: definitions of the phase components in terms of the system components (i.e. tangent plane)
!	Rows numph+2 to numph+2+numph: Gibbs-Duhem equations for each phase
!
!	Example where the only phase is olivine (Fo-Fa)
!	      dÂµFo       dÂµFa       dÂµSiO      dÂµMgO      dÂµFe       dG    
!	   0.0000E+00  0.0000E+00  0.6000E+00  0.1000E+00  0.3000E+00 -0.1000E+01
!	  -0.1000E+01  0.0000E+00  0.3333E+00  0.6667E+00  0.0000E+00  0.0000E+00
!	   0.0000E+00 -0.1000E+01  0.3333E+00  0.0000E+00  0.6667E+00  0.0000E+00
!	   0.8000E+00  0.2000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00


! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
	include "Tangent.inc"
! c*****************************************
! --------------------------------------------------------------------
! --------------------------------------------------------------------

!     Local arrays
!     DIMENSIONED FOR 12 SYSTEM COMPONENTS AND 25 PHASE COMPONENTS
      INTEGER*4 L,j,jj,k,i,kCur,npTemp,izero,irow,jcol,nrows,ncols,variance,ier
      integer*4 indy(15),depVarNumber(100)
      real*8 deta,dG
!      REAL*8 A(50,50),TEMP,FACT
!      common /REXNcommon/A,Temp,Fact
      	Character*4 VN3(100),depvarname(100),indyvarname(15)
	common/VN3names/VN3
! c*****************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	This part is exactly like the REXN routine
!     Set up array of phase component names
	L=0
	npTemp = 0
	jcol = 0
!	asmCurrent should include all phases in the current assemblage
	Do 10 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	DO 10 J=1,numPhCo(k)
	jcol = jcol + 1
        VN3(jcol)=phCoName(k,J)
	npTemp = npTemp + 1
10  CONTINUE
	np = npTemp			! total number of phase components in the assemblage
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Add on "components" for the tangent phase
	do 15 i = 1,nc
	npTemp = npTemp + 1
	vn3(npTemp) = coname(i)
15	continue
!	Add the dG variable
	npTemp = npTemp + 1
	vn3(npTemp) = 'dG'
	ncols = npTemp
	nrows = 1 + np + numCurrent
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	zero out matrix
	do 16 i = 1,nrows
	do 16 j = 1,ncols
	AA(i,j) = 0.0D0
16	continue
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Define the system free energy
	iRow = 1
	jcol = np		! number of phase components
        DO 22 i=1,NC
        jcol = jcol + 1
        AA(irow,jcol) = sysCoMoleFraction(i)
22    CONTINUE
	AA(irow,npTemp) = -1
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Add equations that define the phase components in terms of the system components
	jcol = 0
	do 30 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	do 30 j = 1,numPhCo(k)
	jcol = jcol + 1
	irow = irow + 1
	AA(irow,jcol) = -1.0D0
30	continue
	irow = 1
	jcol = np
	do 32 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	do 32 j = 1,numPhCo(k)
	irow = irow + 1
	do 33 i = 1,nc
	AA(irow,jcol+i) = comp(k,j,i)/atomNorm(k,j)
33	continue
32	continue
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	build Gibbs-Duhem equations
	jcol = 0
!	irow = continue from last value
	do 40 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	irow = irow + 1
	do 40 j = 1,numPhCo(k)
	jcol = jcol + 1
	AA(irow,jcol) = XPhCo(k,j)
40	continue

!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Pick the independent variables
	variance = ncols - nrows
	indy(1) = ncols			! the first independent variable is always dG of the system
	jcol = np
	do 50 i = 2,variance
	jcol = jcol + 1
	indy(i) = jcol
50	continue
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Generate modified master matrix
	jcol = 0
	do 110 j = 1,ncols
	do 112 k = 1,variance
	if(j.eq.indy(k))go to 110
112	continue
!	if here, then this is a dependent variable
	jcol = jcol + 1
	do 120 i = 1,nrows
	A(i,jcol) = AA(i,j)
	depVarName(jcol) = vn3(j)
	depVarNumber(jcol) = j
120	continue
110	continue
!	now do dependent variables
	do 125 k = 1,variance
	jcol = jcol + 1
	do 126 i = 1,nrows
	A(i,jcol) = -AA(i,indy(k))
	indyVarName(k) = vn3(indy(k))
126	continue	
125	CONTINUE

!++++++++++++++++++++++++++++++++++++++++++++++++++
!     ECHO INPUT DATA
      IF(iLong(7).eq.2)then
		write(12,*)
		write(12,*)
		write(12,*)'Original Matrix'
		write(12,2000)(vn3(j),j=1,ncols)
		DO 150 i=1,nrows
		write(12,2001)(AA(i,j),j=1,ncols)
150      	CONTINUE
		write(12,*)
		write(12,*)' Rows, columns = ',nrows,ncols
		write(12,*)'      Variance = ',variance
		write(12,*)
		write(12,*)'Modified Matrix'
		write(12,2000)(depVarName(i),i=1,nrows),(indyVarName(i),i=1,variance)
2000    	FORMAT(' ',20(4X,A4))
		DO 155 I=1,nrows
		write(12,2001)(A(I,J),J=1,ncols)
2001    	FORMAT(' ',2X,20(F8.3))
155      	CONTINUE
		endif
	
!++++++++++++++++++++++++++++++++++++++++++++++++++
	deta = 0.0d0
	ier = 0
	call reduce(nrows,ncols,deta,ier)
	if(ier.eq.1)then
		write(12,*)'Problem inverting matrix'
      		write(12,*)'--------------------'
	      	write(12,*)'Reduced matrix'
		write(12,2000)(adjustl(depVarName(j)),j=1,neq),(adjustl(indyVarName(j)),j=1,variance)
      		DO 231 I=1,NEQ
      		write(12,2001)(A(I,J),J=1,NVAR)
231  		continue
		pause 'hit return to continue'
		return
		endif
      IF(iLong(7).eq.1.or.iLong(7).eq.2)then
		write(12,*)''
		write(12,*)'Determinant   ',deta
		write(12,*)' Jacobian'
		write(12,3502)(vn3(INDY(K)),K=1,VARIANCE)
3502		format(10x,15(A15))
		DO 252 I=1,nrows
      		write(12,3501)adjustl(depVarName(i)),(XX(i,j),j=1,variance)	! writes out soln in last XX
3501  		FORMAT(A6,15E15.5)
252  		CONTINUE
		write(12,*)'Determinant   ',deta
		write(12,*)''
		ENDIF
!	Multiply Jacobian by dG and calculate new values of tangent plane [tangentPlane(i)]
!	The Jacobian is in array XX(i,j) where i=rows (each dependent variable) and j = columns (the first is dG)
!	The first np rows contain the dependent phase components for the current assemblage
!	ncols (total) = np + nc + 1
!	The next nc -
!	np = total number of phase components (all are dependent variables)
!	nc = total number of system components (some are dependent, some are independent)
!	ncols = total number of variables = np + nc + 1
!	nrows = total number of equations = 1 + np + numPh
!	variance = ncols - nrows = np + nc + 1 - 1 - np - numPh = nc - nph (i.e. the phase rule at constant T and P)
	dG = 1.0d0
	j = np
	jj = variance - 1
	do 310 i = 1,nc
	irow = i + jj		! this should be the index of the tangentPlane element
	j = np + i
	tanPlane(irow) = tanPlane(irow) + dG*xx(j,1)
310	continue
      IF(iLong(7).eq.1.or.iLong(7).eq.2)then
		write(12,*)
		write(12,*)'Tangent plane'
		write(12,*)(tanPlane(i),i=1,nc)
		endif

      	return
      	end

! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      	SUBROUTINE BuildJacobian2(izero)
 	use MatrixArrays
     	implicit none
!     	ROUTINE TO build the matrix for calculation of the Jacobian for the chemical potentials
!	The matrix is set up as follows:
!	Row 1: definition of the dG for the bulk composition
!	Rows 2 to numph+1: definitions of the phase components in terms of the system components (i.e. tangent plane)
!	Rows numph+2 to numph+2+numph: Gibbs-Duhem equations for each phase
!
!	Example where the only phase is olivine (Fo-Fa)
!	      dÂµFo       dÂµFa       dÂµSiO      dÂµMgO      dÂµFe       dG    
!	   0.0000E+00  0.0000E+00  0.6000E+00  0.1000E+00  0.3000E+00 -0.1000E+01
!	  -0.1000E+01  0.0000E+00  0.3333E+00  0.6667E+00  0.0000E+00  0.0000E+00
!	   0.0000E+00 -0.1000E+01  0.3333E+00  0.0000E+00  0.6667E+00  0.0000E+00
!	   0.8000E+00  0.2000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00


! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
	include "Tangent.inc"
! c*****************************************
! --------------------------------------------------------------------
! --------------------------------------------------------------------

!     Local arrays
!     DIMENSIONED FOR 12 SYSTEM COMPONENTS AND 25 PHASE COMPONENTS
      INTEGER*4 L,j,jj,k,i,kCur,npTemp,izero,irow,jcol,nrows,ncols,variance,ier
      integer*4 indy(15),depVarNumber(100)
      real*8 deta,xTemp,dX(15)
!      REAL*8 A(50,50),TEMP,FACT
!      common /REXNcommon/A,Temp,Fact
      	Character*4 VN3(100),depvarname(100),indyvarname(15)
	common/VN3names/VN3
! c*****************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	This part is exactly like the REXN routine
!     Set up array of phase component names
	L=0
	npTemp = 0
	jcol = 0
!	asmCurrent should include all phases in the current assemblage
	Do 10 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	DO 10 J=1,numPhCo(k)
	jcol = jcol + 1
        VN3(jcol)=phCoName(k,J)
	npTemp = npTemp + 1
10  	CONTINUE
	np = npTemp			! total number of phase components in the assemblage
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Add on "components" for the tangent phase
	do 15 i = 1,nc
	npTemp = npTemp + 1
	vn3(npTemp) = coname(i)
15	continue
!	Add the dG variable
	npTemp = npTemp + 1
	vn3(npTemp) = 'dG'
	ncols = npTemp
	nrows = 1 + np + numCurrent
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	zero out matrix
	do 16 i = 1,nrows
	do 16 j = 1,ncols
	AA(i,j) = 0.0D0
16	continue
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Define the system free energy
	iRow = 1
	jcol = np		! number of phase components
        DO 22 i=1,NC
        jcol = jcol + 1
        AA(irow,jcol) = sysCoMoleFraction(i)
22    	CONTINUE
	AA(irow,npTemp) = -1
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Add equations that define the phase components in terms of the system components
	jcol = 0
	do 30 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	do 30 j = 1,numPhCo(k)
	jcol = jcol + 1
	irow = irow + 1
	AA(irow,jcol) = -1.0D0
30	continue
	irow = 1
	jcol = np
	do 32 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	do 32 j = 1,numPhCo(k)
	irow = irow + 1
	do 33 i = 1,nc
	AA(irow,jcol+i) = comp(k,j,i)/atomNorm(k,j)
33	continue
32	continue
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	build Gibbs-Duhem equations
	jcol = 0
!	irow = continue from last value
	do 40 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	irow = irow + 1
	do 40 j = 1,numPhCo(k)
	jcol = jcol + 1
	AA(irow,jcol) = XPhCo(k,j)
40	continue

!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Pick the independent variables
	variance = ncols - nrows
!	The independent variables will be the first "variance" of the system components
!	indy(1) = ncols			! the first independent variable is always dG of the system
!	jcol = np
!	jcol = np+1		!See if a small mod will render matrix non-singular
	jcol = ncols+1		! start at the end (dGsys) and move down to pick independent variables
	!				This code might need improving if matrix becomes singular
	do 50 i = 1,variance
	jcol = jcol - 1
	indy(i) = jcol
50	continue
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Generate modified master matrix
	jcol = 0
	do 110 j = 1,ncols
	do 112 k = 1,variance
	if(j.eq.indy(k))go to 110
112	continue
!	if here, then this is a dependent variable
	jcol = jcol + 1
	do 120 i = 1,nrows
	A(i,jcol) = AA(i,j)
	depVarName(jcol) = vn3(j)
	depVarNumber(jcol) = j
120	continue
110	continue
!	now do dependent variables
	do 125 k = 1,variance
	jcol = jcol + 1
	do 126 i = 1,nrows
	A(i,jcol) = -AA(i,indy(k))
	indyVarName(k) = vn3(indy(k))
126	continue	
125	CONTINUE

!++++++++++++++++++++++++++++++++++++++++++++++++++
!     ECHO INPUT DATA
      IF(iLong(7).eq.2)then
		write(12,*)
		write(12,*)
		write(12,*)'Original Matrix'
		write(12,2000)(vn3(j),j=1,ncols)
		DO 150 i=1,nrows
		write(12,2001)(AA(i,j),j=1,ncols)
150      	CONTINUE
		write(12,*)
		write(12,*)' Rows, columns = ',nrows,ncols
		write(12,*)'      Variance = ',variance
		write(12,*)
		write(12,*)'Modified Matrix'
		write(12,2000)(depVarName(i),i=1,nrows),(indyVarName(i),i=1,variance)
2000    	FORMAT(' ',50(4X,A4))
		DO 155 I=1,nrows
		write(12,2001)(A(I,J),J=1,ncols)
2001    	FORMAT(' ',2X,50(F8.3))
155      	CONTINUE
		pause 'pausing'
		endif
	
!++++++++++++++++++++++++++++++++++++++++++++++++++
	deta = 0.0d0
	ier = 0
	call reduce(nrows,ncols,deta,ier)
	if(ier.eq.1)then
		write(12,*)'Problem inverting matrix in Subroutine BuildJacobian2'
      		write(12,*)'--------------------'
	      	write(12,*)'Reduced matrix'
		write(12,2000)(adjustl(depVarName(j)),j=1,nrows),(adjustl(indyVarName(j)),j=1,variance)
      		DO 231 I=1,NROWS
      		write(12,2001)(A(I,J),J=1,NCOLS)
231  		continue
		pause 'hit return to continue'
		return
		endif
      IF(iLong(7).eq.1.or.iLong(7).eq.2)then
		write(12,*)''
		!write(12,*)'Determinant   ',deta
		write(12,*)' Jacobian'
		write(12,3502)(vn3(INDY(K)),K=1,VARIANCE)
3502		format(10x,50(A15))
		DO 252 I=1,nrows
      		write(12,3501)adjustl(depVarName(i)),(XX(i,j),j=1,variance)	! writes out soln in last XX
3501  		FORMAT(A6,50E15.5)
252  		CONTINUE
		write(12,*)'Determinant   ',deta
		write(12,*)''
		ENDIF
!	Multiply Jacobian by dG and calculate new values of tangent plane [tangentPlane(i)]
!	The Jacobian is in array XX(i,j) where i=rows (each dependent variable) and j = columns (the first is dG)
!	The first np rows contain the dependent phase components for the current assemblage
!	ncols (total) = np + nc + 1
!	The next nc -
!	np = total number of phase components (all are dependent variables)
!	nc = total number of system components (some are dependent, some are independent)
!	ncols = total number of variables = np + nc + 1
!	nrows = total number of equations = 1 + np + numPh
!	variance = ncols - nrows = np + nc + 1 - 1 - np - numPh = nc - nph (i.e. the phase rule at constant T and P)
!
!	The last row of xx(,) contains the partial derivatives for the dGsystem
!	The amount the independent variables change will be the amount necessary to change dGsys by 1

!	dG = 1.0d0
	do 310 i = 1,variance 	! loop on the independent variables
	xTemp = xx(nrows,i)
	if(dabs(xTemp).lt.1.0d-2)then
		dX(i) = 10.0d0		! this is the maximum allowable change
		else
		dX(i) = 1.0d0/xx(nrows,i)
		endif
	tanPlane(i) = tanPlane(i) + dX(i)
310	continue
!	Now calculate changes in the dependent tangent plane variables
	jj = np
	do 320 j = variance+1,nc
	jj = jj+1
	do 325 i = 1,variance
	xTemp = xTemp + dX(i)*xx(jj,i)
325	continue
	tanPlane(j) = tanPlane(j) + xTemp	
320	continue

      IF(iLong(7).eq.1.or.iLong(7).eq.2)then
		write(12,*)
		write(12,*)'Tangent plane'
		write(12,*)(tanPlane(i),i=1,nc)
		endif

      	return
      	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Combinatorial3(np,npt,index,inum)
! 	Note that this subroutine is recursive
	integer*4 np,npt,index,inum
	integer*4 i1
	integer*4 in(20)
	common /combo/in
	integer combos(10000,15),comboCounter
	common /comboarray/combos,comboCounter
	
	do 10 i1 = inum,npt
	in(index) = i1
	if(np.eq.index)then
!		write(iunit,*)(in(i),i=1,np)
		comboCounter = comboCounter + 1
		do 12 i = 1,np
		combos(comboCounter,i) = in(i)
12		continue
		go to 10
		endif
	call Combinatorial3(np,npt,index+1,i1+1)
10	continue
	return
	end


