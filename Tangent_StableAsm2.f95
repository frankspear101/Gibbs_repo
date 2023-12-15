! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine StableASM2(logFileOn,TC4,PB4,ierr)
    	use AWE_Interfaces

!	This routine is designed to calculate the stable phase assemblage for a specified bulk composition 
!		at the specified P and T
!	The routine runs through all possible assemblages

!! 	Here's how I implement it
!
!	For each P-T point - Do this:
!		(A) Calculate Go (o) for each phase component at P&T
!		(B) Read the first assemblage from 87
!		(C) Try to solve for this assemblage (in Compute)
!			(i) If successful then
!				(a) Check to see if all moles are positive. 
!					- If so then proceed to (D)
!					- If Not then remove phase with negative moles and go back to (C) and try again
!			(ii) If not successful (any error code from Compute) then go back to (B) and read next assemblage and try again
!		(D) Calculate tangent plane and determine whether
!			(i) All phases lie on the tangent (i.e. in the equilibrium assemblage) or above.
!				(a) If so, then we are done. 
!					- Write information to unit 73 (output file) then
!					- Go to (4) and get an new PT point
!			(ii) One or more phases lie below the tangent. 
!				(a) Find the phase that lies furthest below the tangent
!					- If the current assemblage has v > 2 (trivariant or higher) then add this phase to the assemblage and go back to (C) to try and solve
!					- If the current assemblage has v = 2 then iterate through the combinations of v = 2 and go back to (C) to try and solve
!		(E) If all else fails, go to (B), read in a new assemblage and try again
!		
!
	use ieee_arithmetic
	use MatrixArrays
      implicit none
!
! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
! 	include "Solute.inc"
! *****************************************
! --------------------------------------------------------------------
	include "Tangent.inc"
	include "MatrixSingular.inc"	
! --------------------------------------------------------------------
!      local variables
      	Character VN3(maxPC)*4
	common/VN3names/VN3
	integer combos(10000,15),comboCounter,numCurrentStore,icombo
	common /comboarray/combos,comboCounter

	integer*4 i,izero,isubtend
	integer*4 ierr,WriteAll
	real*4 Tc4,Pb4
	integer*4 logFileOn
!-----
!	these only needed for testing in block 46
	real*8 R
	DATA R/8.3144D0/

! ---------------------

	ierr = 0
!	tanPhase = 0
	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
!	iLong(1) = 1		! turn on error reporting
!	iLong(3) = 1		! turn on output thermo data
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	iDoingGrid = 1
!	ioutputpseudo = 1			! allows error output in subroutine compute
	WriteAll = 0	
	call ZeroTanOut()

	write(*,*)'   '
	write(*,*)' Phases in this MIF'
	do 228 i = 1,numPhMIF
228	write(*,*)i,phName(i)
	write(*,*)'   '
	write(*,*)'Please provide your best guess for the initial assemblage.'
	write(*,*)'Input numbers for these phases (0 to end list)'
	write(*,*)' From the MIF:'
	write(*,*)AssembTitle
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


	call SetLast()
	call SetNewP()
	write(*,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
402	format(2F12.2,20F12.5)

436	continue
	write(*,*)' Starting TP calculations'

	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	if(izero.gt.0)then
		write(*,*)' Error calculating MU zero for this T&P... abort the run'
		write(95,*)' Error calculating MU zero for this T&P... abort the run'
		pause 'take a look and hit return when ready'
		ierr = 1	! error code
		return
		endif


	call ResetToStart()	! start with initial compositions when examining a new assemblage

	write(*,*)numCurrent,(asmCurrent(i),i=1,numCurrent)
420	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
!	call dumpeverything('Before Compute')
	if(isubtend.eq.0)then
		write(*,*)'      Assemblage does not subtend BC (line 168 in Sub StableAsm) - get new asm'
		pause 'hit return to continue'
		go to 426
		endif
	izero = 0
	call ComputePTXM(TC4,PB4,izero)
!	call dumpeverything('After Compute')
	if(izero.gt.0)then
		write(*,*)'      Error in ComputePTXM (Line 175 of Sub StableASM2)'
		pause 'hit return to continue'
!		xPhCo = xPhCoNewP		! reset to the initial compositions
!		MP0 = MPnewP
!		call ResetToLast()
!	do 20 i = 1,numPhMIF
!	MP0(i) = MPLast(i)
!	do 24 j = 1,numPhCo(i)
!	xPhCo(i,j) = xPhCoLast(i,j)
!24	continue
!20	continue

!		xPhCo = xPhCoLast		! reset to the the last compositions that worked
!		MP0 = MPLast
		go to 426
		endif
425	continue
!	If here, then we have found a solution using this assemblage

	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	

426	continue
	write(*,*)' Current assemblage'
	write(*,*)numCurrent,(asmCurrent(i),i=1,numCurrent)

	call WriteAllStuffToTerminal(6)
	write(*,*)' 0 = Change assemblage; 1 = Accept these results and exit; 2 = abort; 3 = Go to Steps'
	write(*,*)' 4 = SetLast array; 5 = ResetToLast; 6 = ResetToStart; 7 = Change T & P'
	read(*,*)WriteAll
	select case (WriteAll)
	case(0)
		numPh = numCurrent
		call Change()
		numCurrent = numPh
		write(*,*)numCurrent,(asmCurrent(i),i=1,numCurrent)
		go to 420
	case(1)
		ierr = 0	! error code
		return
	case(2)
		ierr = 1	! error code
		return
	case(3)
		call Steps()
		write(*,*)' If P or T has changed, you need to set this using option 6 before continuing'
		go to 426
	case(4)
		call SetLast()
		numPh = numCurrent
		go to 426
	case(5)
		call ReSetToLast()
		numPh = numCurrent
		go to 426
	case(6)
		call ReSetToStart()
		go to 426
	case(7)
		write(*,*)' Specify new T and P for calculations'
		read(*,*)TC,PB
		TC4 = TC
		Pb4 = PB
		go to 436
	case default
		go to 426
	end select





	call SetLast()
	call SetNewP()

400	continue			! loop back to get the next temperature
	ierr = 0		! no error
	return
	end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PickSecondBulkComposition()
      implicit none
! 	Routine to pick a bulk rock composition from the file "BulkRockAnalyses" in Gibbs_Essentials folder
!	This routine is only used from Sub AutoTan_BC in order to pick a second bulk composition to use in calculations
!	After selection, a list of bulk compositions is created and stored in a temporary file that the program uses for sequential execution

! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
! *****************************************
!      	local variables
      	INTEGER i,j,iok
      	CHARACTER TempTitle*80,dummy*4
	Character*80 TheList(200),WhatToPick,BulkCompTitle2
	Integer*4 TheOne,ListSize,iselect
	REAL*8 molesStart2(sysCoMax),bulkCompMoles2(sysCoMax),bulkCompWt2(sysCoMax),delta(sysCoMax),newMoles(sysCoMax),increments

! -----------------------------------------------

	open(unit=24,file=BulkRockCompositionFile,status="OLD")

!      	Build a list of analysis titles for listmanager
      	READ(24,'(A)')dummy			! read the header
	j = 0
     	DO
		READ(24,'(A)',END=10)dummy          ! read the row of dashes
		READ(24,'(A)',END=10)TempTitle
		j = j+1
		if(j.gt.200)go to 10			! maximum dimensions of TheList(200)
		TheList(j) = Trim(TempTitle)
		DO I=1,3
			READ(24,'(A)',END=10)dummy	! read 3 lines
			end do
		REPEAT

10     CONTINUE

!     Put the list up and let the user pick one plot type
	iSelect = 1
	ListSize = j
	WhatToPick = 'Select a bulk composition '

!	Sub  AWE_Pick_one(TheOne,ListLength,TheList,WhatToPick,iok,iselect)
      	call AWE_Pick_One(TheOne,ListSize,TheList,WhatToPick,iok,iSelect)
      	if(iok.ne.0)then
		bulkCompSwitch = .false.
		return
		endif

	j = TheOne
	REWIND(24)
	READ(24,'(A)')dummy			! read the header
	DO i=1,(j-1)*5
		READ(24,'(A)')dummy     		! ignore the unwanted plot definitions
		end do
	READ(24,'(A)')dummy     		! Read line of dashes
	READ(24,'(A)')bulkCompTitle2
	read(24,*)nc				! number of components
	READ(24,*)dummy				! names of oxides
	read(24,*)(bulkCompWt2(i),i=1,nc)	! bulk composition in wt %
	close(24)

!	Convert to moles
	do i=1,nc
		bulkCompMoles2(i)=bulkCompWt2(i)*numCatInOxide(i)/molwt(i)
!		bulkCompWtSum2=bulkCompWtSum2 + bulkCompWt(i)
		molesStart2(i) = bulkCompMoles2(i)	! molesStart is the array used in Compute
!		bulkCompMolesStart2(i) = molesStart2(i)		! bulkCompMolesStart is the initial input composition in case we need to reset to starting values
!		bulkCompWtStart(i) = bulkCompWt(i)		! save the initial composition in case we need it later
		end do

	! interpolate between the 2 bulk compositions and make a file with the intermediate values
	open(unit=25,file="TempBulkRockMolesFile",status="UNKNOWN")
	write(25,*)(molesStart(i),i=1,nc)
	write(*,*)' How many increments do you want (e.g. 10 or so)?'
	read(*,*)increments
!	increments = 10.
	do j = 1,increments
		do i = 1,nc
			delta(i) = (molesStart(i) - molesStart2(i))/increments
			newMoles(i) = molesStart(i) + delta(i)*float(j)
			end do
		write(25,*)(newMoles(i),i=1,nc)
		end do
	rewind(25)

!	bulkCompSwitch = .true.			! a bulk composition is now open
	return
	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
