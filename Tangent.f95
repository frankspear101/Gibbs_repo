! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Tangent
    	use AWE_Interfaces

!	This routine is designed to calculate the stable phase assemblage for a specified bulk composition 
!! 	Here's how I implement it
!
!	1) Pick bulk composition to model
!	2) Pick list of phases to considerc
!		- read in thermodynamic data for all phases as well as compositions
!		- Run the combinatorial formula to generate a file (unit 84) with all possible assemblages up to v=2
!		- Sort through to weed out assemblages that do not subtend the bulk composition
!			- store result in unit 85
!		- Ask the user if there are "preferred phases"
!		- Sort through assemblages in unit 85 to make a list of assemblages with "preferred phases" (if any)
!			- Write to unit 87
!			- Write others to unit 86
!			- After all "preferred phase" assemblages are in 87 write the contents of 86 to the end of 87
!				- This ensures that all assemblages are represented and the "preferred phase" assemblages are listed first
!	3) Pick array of P-T points to consider
!	4) For each P-T point - Do this:
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

	integer*4 iTan,i,ii,j,jj,L,iTP,autorun,iok,status,k,izero,isubtend,maxNegativeMolesPhase,maxgDifferencePhase,kGrt
	integer*4 asmCurrentStore(phMax),jRemove,rewind87,iAllAsm(phMax),iAll,inew,ibegin,kCur,newP,countTheMisses,debug
	real*4 Tc4,Pb4,Pb4old,Ttemp,Ptemp,gridset,longhit
	real*8 maxgDifference,maxNegativeMoles,xTemp(5),xvalue
	real*8 molesSyComp,molesFluidToChange,molesFluidNew,molesChange,weightChange,sum,mGrt,mGrtTotal,radius,shiftG
	integer*4 logFileOn,iGrid,iremove,ierr
	integer*4 niterate,EOFflag,doesGarnetFractionate
	character*64 PTinputFile
	character*132 aline
!-----
	real*8 Xalmcore,Xspscore,Xgrscore,XAnCore,XalmCalc,XspsCalc,XgrsCalc,XAnCalc,						&
	XalmCalcPlus,XspsCalcPlus,XgrsCalcPlus,XAnCalcPlus,XalmCalcMinus,XspsCalcMinus,tol,mincFD,  	&
	DalmDmolesFe,DalmDmolesMn,DalmDmolesCa,DalmDmolesNa,DspsDmolesFe,DspsDmolesMn,DspsDmolesCa,DspsDmolesNa,		&
	DgrsDmolesFe,DgrsDmolesMn,DgrsDmolesCa,DgrsDmolesNa,DAnDmolesFe,DAnDmolesMn,DAnDmolesCa,DAnDmolesNa,  			&
	Deta,damp
	integer*4 iFe,iMn,iCa,iNa,iCount,nsolve,ier,iFrac,kPlg
!	these only needed for testing in block 46
	real*8 R,Rjoules
	DATA R/8.3144D0/

! ---------------------

	Rjoules = R
!	ioutputPseudo = 1		! default output is summary only
	ioutputPseudo = 3		! default output is summary only
	ioutputPseudo = 4		! lots of output
	iDoingGrid = 1				! keep going through activity calculations if we are doing a grid

	autorun = 0			! default is auto run
!	autorun = 1			! default is auto run
	numPTpoints = 1		! default to do 1 point
	Tgrid(1) = 500
	Pgrid(1) = 4000		! default P&T for testing

!	Tcstart = 450
!	TCend = 750
!	Tinc = 2
!	Pbstart = 2000
!	Pbend = 12000
!	Pinc = 100
	gridset = 0
	TC = tcstart
	Pb = pbstart
	Pfluid = pbstart
!	tanPhase = 0

	tolnewton = 1.0d-9		! we need to set this small to get within 1 joule for G asm
!	tolnewton = 1.0d-3		! Convergence is on compositions so .1% should be fine
	iNewton = 1				! this only works if Newton's method is switched to ON

	call ZeroTanOut()
	logFileOn = 1		! default is to make a log file
!	open(15,file='./Gibbs_Essentials/ASMCombos.txt',status = 'OLD')
!10	continue
!	read(15,*,end=100)j,jj
!	numASMCombos(j) = jj
!	do 15 i = 1,jj
!	read(15,*,end=100)ASMCombo(j,i,1),(ASMCombo(j,i,k),k=2,ASMCombo(j,i,1)+1)
!15	continue
!	go to 10
100	continue
	iLong(1) = 1		! turn on error reporting
	iExclude = 0



	write(*,*)'====================================='
	write(*,*)'====================================='
	write(*,*)'====================================='
     	If(inewton.eq.1.and.imass.eq.1)write(*,*)'Newton''s method; Mass balance = ON'
     	If(inewton.eq.1.and.imass.eq.0)write(*,*)'Newton''s method; Mass balance = OFF'
     	If(inewton.eq.0.and.imass.eq.1)write(*,*)'Gibbs'' method; Mass balance = ON'
     	If(inewton.eq.0.and.imass.eq.0)write(*,*)'Gibbs'' method; Mass balance = OFF'
	if(imass.eq.1)then
		if(bulkCompSwitch)then
			write(*,*)'Bulk composition = ',bulkCompTitle
		else
			write(*,*)'Bulk composition = (none - bulk comp may change)'
		endif
	endif
	write(*,*)'====================================='
	write(*,*)'Tangent menu '
	write(*,*)'  0 = return'
	write(*,*)'---------------------------'
! 	write(*,*)'Pick P and T endpoints for PT-X diagram '
	write(*,*)'  1 = Pick a bulk composition'
 	write(*,*)'  2 = Choose minerals to consider'
 	write(*,*)'  3 = Pick P&T grid to examine'
	write(*,*)'---------------------------'
!	write(*,*)' 4 = Auto calc  --- this is no longer supported. Please use option 42'
!	write(*,*)'41 = Find stable assemblage at T&P -- Please use option 410 instead'
	write(*,*)'410 = Find stable assemblage at T&P -- manual code'
!	write(*,*)'411 = Find stable assemblage at new bulk composition -- manual code'
	write(*,*)' 42 = Auto calc (new algorithm)'
!	write(*,*)'421 = Auto calc (new algorithm - Takes tiny steps if solution is not found)'
	write(*,*)' 44 = Auto calc with fractionation - Equilibrium only'
	write(*,*)' 45 = Auto calc with fractionation - OS model'
! 	write(*,*)' 5 = Pick a P&T and calc (debug mode)'
! 	write(*,*)' 25 = New -- Pick a P&T and calc (debug mode)'
	write(*,*)'---------------------------'
	write(*,*)'Other MAD routines'
	write(*,*)'  6 = Turn log file on or off'
 	write(*,*)'  7 = Call Printt'
 	write(*,*)'  8 = Call GLOBAL'
 	write(*,*)'  9 = Set plotting axes'
	write(*,*)' 10 = Gibbs/Newton switch'
 	write(*,*)'------Debug----------------'
! 	write(*,*)' 11 = Calculate affinity file'
! 	write(*,*)'      (this replaces the OverstepGarnet routines by calculating all affinities on the P-T grid)'
! 	write(*,*)' 12 = Grow garnet in overstepped environment with fixed assemblage and garnet fractionation'
! 	write(*,*)' 13 = Affinity node nucleation calculations'
! 	write(*,*)' 14 = Grow garnet in equilibrium with fractionation along a path'
! 	write(*,*)' 15 = Check out activity calculation for HoPo dataset'
! 	write(*,*)' 16 = Speed test'
! 	write(*,*)' 17 = Grow garnet in OS environment (PT path) with assemblage evolution and garnet fractionation'
! 	write(*,*)' 18 = Calculate tangent function for debugging'
! 	write(*,*)' 19 = Grow garnet in OS environment (constant P&T) with assemblage evolution and garnet fractionation'
! 	write(*,*)'---------------------------'
! 	write(*,*)' 20 = Jacob"s ladder'
! 	write(*,*)' 21 = Calculate EBC for garnet core composition OS model -- requires P&T of nucleation (3 elements)'
! 	write(*,*)' 22 = Calculate EBC for garnet core composition EQ or OS model -- requires P&T of nucleation (2 elements)'
! 	write(*,*)' 23 = Calculate EBC for garnet core composition EQ or OS model -- requires P&T of nucleation (4 elements)'
! 	write(*,*)' 25 = Debug routine'
! 	write(*,*)' 26 = Jacobi''s ladder version 2'
	write(*,*)'---------------------------'
 	write(*,*)' 30 = Make parser conversion file for this MIF and thermo data file'
! 	write(*,*)' 50 = Tangent AFM routine'
	write(*,*)' -2 = Set global output length'

	read(*,*)iTan

	if(iTan.eq.0)go to 999
! ---------------------------------------
	if(iTan.eq.30)then
		call MakeParserFile
		goto 100
		endif
! ---------------------------------------
	if(iTan.eq.-2)then
		call AWE_setoutput
		goto 100
		endif
! ------Close output file-------------
	if(iTan.eq.6)then
		call SetTanOut()
		write(*,*)'logFileOn = ',logFileOn
		write(*,*)'0 = turn logFile OFF'
		write(*,*)'1 = turn logFile ON'
		write(*,*)'2 = Send log file to window'
		read(*,*)logFileOn
		go to 100
		endif
! ------Pick bulk composition-------------
	if(iTan.eq.1)then
		call PickBulkComp
!		do 110 i = 1,nc
!		bulkCompWtStart(i) = bulkCompWt(i)				! save the starting bulk composition - we can change bulkCompWt if we want to fractionate
!		actually, this isn't correct. the Wt% arrays aren't used again
!		molesStart(i) is the array used in COMPUTE to do the mass balance
!		bulkCompositionMoles (calculated in FileInput) contains the same information
!110		continue
!	this is done in FileInput because we used to have to read in different assemblages with the same BC. No longer done that way
!		 - I should moce this code to PickBulkComp
		go to 100
		endif
! ------Choose minerals to consider-------------
	if(iTan.eq.2)then
! 		open MIF and use all phases in file
		write(12,*)'-------Open MIF file --------------------'
		INEW=0
		numNew=0
		iBegin = 1
		CALL BEGIN(IBEGIN,INEW)
	        write(*,*)' Input file name: ',FILEIN
	        write(12,*)' Input file name: ',FILEIN

		write(*,*)' NC = ',nc
!		call CheckSystemComponents(iok)	! this doesn't work because the system components are Si, Al vs SiO2, Al2O3 etc
!		if(iok.eq.1)then
!			write(*,*)'The system components in the bulk composition file do not'
!			write(*,*)'    match those in the MIF file'
!			write(*,*)' Please fix!!'
!			go to 100
!			endif
		numPhMIF = numPh				! Total number of phases in the input MIF file
		IF(INEW.EQ.1)then       			! if true, we actually started a new problem
			isfileopen = 1
			else
			write(*,*)'File did not open -- Perhaps user hit cancel'
			pause 'Hit return to continue'
			go to 100
			endif
!!		printMIF = 0			! this will supress printing monitor parameters (causes a bomb)
!		replace with call Printt(0) as needed to supress printing monitor parameters
		call PrintMinAssemblage
!		call LoadAllToMIF
		call CalculateAtomUnits
!		call SetInitialTangent		! I don't think this is needed anymore because the tangent isn't used until we have an initial solution
		call CalculateSystemMoleFractions
!		call ListTangentPlane
!		Save the initial compositions read in from the file to an array
		xPhCoInitial = xPhCo			! this should assign the entire array - store initial compositions in case we need them
!		Generate file with codes for assemblages to test
!		do 222 i = nc,1,-1				! this will find all possibilities from numPh = 1 to numPh = nc (divariant)
!		do 222 i = nc-1,1,-1			! this will find all possibilities from numPh = nc-1 (trivariant) to numPh = 1 


		go to 255	! skip the combinatorial routines -- see what happens


		do 222 i = 1,nc,1				! this will find all possibilities from numPh = 1 to numPh = nc (divariant)
		write(*,*)'Finding combinations:   ',i
		call Combinatorial2(84,i,numPhMIF,1,1)
222		continue
		rewind(84)
225		continue
		read(84,*,end=227)numCurrent,(asmCurrent(i),i=1,numCurrent)
		isubtend = 1
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)	! don't need to load everything, just subtend
!		call subtend(isubtend)
!		LoadCurrentAsm also checks for degeneracies, which subtend doesn't
		if(isubtend.eq.0)then
			go to 225
			else
			write(85,*)numCurrent,(asmCurrent(i),i=1,numCurrent)
			!call ListCurrentAsm
			go to 225
			endif
227		continue
		close(84)


		go to 255



255		continue
		write(12,*)'=========Initial setup done============================='
		write(*,*)'=========Initial setup done============================='
		npMIF = np		! total number of phase components in the MIF
		go to 100
		endif

! ------ pick P-T grid for calculations
	if(iTan.eq.3)then
	write(*,*)' 0 = return'
	write(*,*)' 1 = Pick PT grid from dialog box -- isobaric mode'
	write(*,*)' 2 = Open file with P-T points or paths'
	write(*,*)' 3 = Pick PT grid from dialog box -- isothermal mode (does not allow constant porosity or fractionation)'
	write(*,*)' 4 = list PT array (after doing 1,2, or 3)'
	read(*,*)iGrid
	if(iGrid.eq.0)go to 100
	if(iGrid.eq.4) then
		do i = 1,numPTpoints
			write(93,335)Tgrid(i),Pgrid(i)
335			format(2F12.1)
			end do
		close(93)
		call fss_Alert('Alert','Your PT array is in file fort.93')
		endif
	if(iGrid.eq.1)then
		call PickPT_grid(Tcstart,Pbstart,TCend,PBend,Tinc,Pinc,longhit)
		if(longhit.eq.2)go to 100		! cancel
! Not sure why this code is here. Probably to keep from bombing if Pinc=0. But we should try negative Pinc
!		if(Pinc.le.0.)then		
!			Pinc = .1
!			endif
		! Fill the arrays Tgrid, Pgrid with the P-T points to examine
		i = 0
		Ttemp = TCstart - Tinc
		Ptemp = PBstart
310		continue
		i = i+1
		Ttemp = Ttemp + Tinc
!		if(Ttemp.gt.TCend.or.Ttemp.lt.TCstart)then

		if(Tinc.ge.0.0d0)then
			if(Ttemp.gt.TCend.or.Ttemp.lt.TCstart)then
				Ptemp = Ptemp + Pinc
	!			Ttemp = Ttemp - Tinc
				Ttemp = TCstart
				if(Pinc.gt.0.)then
					if(Ptemp.gt.Pbend)then	! we are increasing P on each increment
						go to 320
						endif
					else			! we are decreasing pressure on each increment
					if(Ptemp.lt.Pbend)then
						go to 320
						endif
					endif				
	!			Tinc = -Tinc
				endif
			else	! Tinc must be negative
			if(Ttemp.lt.TCend)then	! if true, then increment P
				Ptemp = Ptemp + Pinc
	!			Ttemp = Ttemp - Tinc
				Ttemp = TCstart
				if(Pinc.gt.0.)then
					if(Ptemp.gt.Pbend)then	! we are increasing P on each increment
						go to 320
						endif
					else			! we are decreasing pressure on each increment
					if(Ptemp.lt.Pbend)then
						go to 320
						endif
					endif				
	!			Tinc = -Tinc
				endif
			endif
		Tgrid(i) = Ttemp
		Pgrid(i) = Ptemp
		go to 310
320		continue
		numPTpoints = i-1
		gridset = 1
		PTinputFile = 'PTgrid'
		endif	! end dialog box routine
	
	if(iGrid.eq.2)then
		write(*,*)' Open file with P-T points'
		open(78,FILE='',status='OLD',iostat=status)
		if(status.ne.0)go to 100
		inquire(78,NAME=PTinputFile)
		i = 0
		read(78,*)aline		! title line in file
322		continue
		i = i + 1
		read(78,*,end=325)Tgrid(i),Pgrid(i)
		go to 322
325		continue
		numPTpoints = i - 1
		gridset = 1
		close (78)
		endif			!end read PT file
	if(iGrid.eq.3)then
		call PickPT_grid(Tcstart,Pbstart,TCend,PBend,Tinc,Pinc,longhit)
		if(longhit.eq.2)go to 100		! cancel
! Not sure why this code is here. Probably to keep from bombing if Pinc=0. But we should try negative Pinc
!		if(Pinc.le.0.)then		
!			Pinc = .1
!			endif
		! Fill the arrays Tgrid, Pgrid with the P-T points to examine
		i = 0
		Ttemp = TCstart
		Ptemp = PBstart - Pinc
330		continue
		i = i+1
!		Ttemp = Ttemp + Tinc
		Ptemp = Ptemp + Pinc
!		if(Ttemp.gt.TCend.or.Ttemp.lt.TCstart)then

		if(Pinc.ge.0.0d0)then
			if(Ptemp.gt.Pbend.or.Ptemp.lt.Pbstart)then
				Ttemp = Ttemp + Tinc
	!			Ttemp = Ttemp - Tinc
				Ptemp = Pbstart
				if(Tinc.gt.0.)then
					if(Ttemp.gt.TCend)then	! we are increasing T on each increment
						go to 333
						endif
					else			! we are decreasing pressure on each increment
					if(Ttemp.lt.TCend)then
						go to 333
						endif
					endif				
	!			Tinc = -Tinc
				endif
			else	! Pinc must be negative
			if(Ptemp.lt.Pbend)then	! if true, then increment P
				Ttemp = Ttemp + Tinc
	!			Ttemp = Ttemp - Tinc
				Ptemp = Pbstart
				if(Tinc.gt.0.)then
					if(Ttemp.gt.TCend)then	! we are increasing P on each increment
						go to 333
						endif
					else			! we are decreasing pressure on each increment
					if(Ttemp.lt.TCend)then
						go to 333
						endif
					endif				
	!			Tinc = -Tinc
				endif
			endif
		Tgrid(i) = Ttemp
		Pgrid(i) = Ptemp
		go to 330
333		continue
		numPTpoints = i-1
		gridset = 1
		PTinputFile = 'PTgrid'
		endif	! end dialog box routine

!	do this for either type of input
	write(*,*)'The number of P-T points in your grid is = ',numPTpoints
	go to 100
	endif

! ------   Auto calc ---  Auto calculate from "Find initial tangent" from assemblage-------------
	if(iTan.eq.4)then
!	tanPhase = 0
	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
!	call SetTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort
	do 412 i = 1,numPhMIF
412	write(*,*)i,phName(i)
	write(*,*)'Is there any phase to exclude from equilibrim assemblages (for calculation of affinities)?'
	write(*,*)'Input a number from the list or 0 for none (i.e. include all phases)'
	read(*,*)iExclude
	if(Tgrid(1).lt.Tgrid(numPTpoints))then
		call ChooseConstantPorosity()
		endif
!	stuffFile = trim(tangentOutputFile)//'.stuff'
!	open(72,FILE=stuffFile,status = 'UNKNOWN')
	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
	select case (logFileOn)
	case(1)
!	if(logFileOn.eq.1)then
		logFile = trim(tangentOutputFile)//'.log'
		open(95,FILE=logFile,status = 'UNKNOWN')
	case(2)
		OPEN(95,FILE = 'MAD LOG',ACCESS = 'window, 800, 1200')
	case default
	end select
!		endif
!	call WriteOutputHeader(73)
	call WriteAllOutHeader(74,Tinc,Pinc)

	call ChooseInitialAssemblage()

	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the first assemblage to check
	Call RewindTheFile(87,logFileOn)
!	rewind87 = 1
!	numCurrent = nc
!	read(87,*)(asmCurrent(i),i=1,numCurrent)
	!write(12,*)'GNA = 339'
	call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)		! get the initial assemblage to test
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	MPLast = MP0
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	MPNewP = MP0
	PB4old = 0
	if(logFileOn.eq.1)then
		write(95,*)' Starting TP loop'
		endif
	do 400 iTP = 1,numPTpoints
	Call RewindTheFile(87,logFileOn)	! rewind the file at every new PT point so we only go through once
	EOFflag = 0
!	RewindCount = 0
	newP = 0
	TC = Tgrid(iTP)
	PB = Pgrid(iTP)
	tc4 = tc
	pb4 = pb
!	write(12,*)TC4,PB4	! for debugging only
	write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
	if(Dabs(PB4-PB4old).gt.1.)then		! this should print every time P changes - just to show progress
!		If we are holding the porosity constant, then we need to reset the H2O content to have it in excess here
		if(constantPorosity.eq.0)then
!			write(12,*)TC4,PB4,(asmCurrent(iok),iok=1,numCurrent)
			write(*,*)TC4,PB4
			else
			if(iGrid.eq.1)then
			! reset the moles of H2O back to the starting value for the next row of calculations
			!	but only if we are doing a grid (not a P-T path)
				bulkCompMoles = bulkCompMolesStart
				molesStart = bulkCompMolesStart
				endif
			write(*,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)
402			format(2F12.2,20F12.5)
			endif
		PB4old = Pb4
		newP = 1
		endif
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	if(izero.gt.0)then
		write(*,*)' Error calculating MU zero for this T&P... abort the run'
		write(95,*)' Error calculating MU zero for this T&P... abort the run'
		pause 'take a look and hit return when ready'
		go to 100
		endif
	if(newP.eq.1)then
!		If we are starting the next pressure, use the results from the previous P at the same T as a starting point
		asmCurrent = asmCurrentNewP		! set the assemblage to the low T, new P values
		numCurrent = numCurrentNewP
		xPhCo = xPhCoNewP
		MP0 = MPNewP
		else
!		Otherwise, use the last assemblage as the starting point
		asmCurrent = asmCurrentLast		! set the assemblage to the last good one
		numCurrent = numCurrentLast
		MP0 = MPLast
		endif
	countTheMisses = 0
410	continue
	if(iParagonite.gt.0)then
!		xPhCo(iParagonite,1) = 0.05
!		xPhCo(iParagonite,2) = 0.95
		endif		
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
!		if(logFileOn.eq.1)then
		write(95,*)' Assemblage does not subtend BC (line 468 in Sub Tangent) - get new asm'
!			endif
		numCurrent = nc
		!if(rewind87.eq.0)then
		!	Call RewindTheFile(87,logFileOn)    !rewind(87)
		!	rewind87 = 1
		!	RewindCount = RewindCount + 1
		!	write(95,*)' Rewind on line 503'
		!	if(RewindCount.ge.2)then
		!		write(*,*)' Every assemblage has been examined and no solution found'
		!		write(*,*)' Skipping this P,T point and moving on to the next'
		!		RewindCount = 0
!				go to 400
!				endif
!			endif
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
		if(EOFflag.eq.1)then
			!write(*,*)' Every assemblage has been examined and no solution found'
			!write(*,*)' Skipping this P,T point and moving on to the next'
			write(*,*)'   No solution..Skipping PT point...T,P ',TC4,PB4
			go to 400
			endif		
		xPhCo = xPhCoLast			! Reset to last good compositions
		if(iParagonite.gt.0)then
		!	xPhCo(iParagonite,1) = 0.05
		!	xPhCo(iParagonite,2) = 0.95
			endif		
		go to 410
		endif
	izero = 0
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
!		if(logFileOn.eq.1)then
		write(95,*)' Error in ComputePTXM (Line 485 of Sub Tangent) - get new asm'
		!write(*,*)'Problem...trying to solve...',TC4,PB4
		!call Printt(1)
		!pause ' hit return'
!			endif
!		numCurrent = nc			! Why this? divariant? This statement is irrelevent because numcurrent is read from file 87
!		if(rewind87.eq.0)then
!			Call RewindTheFile(87,logFileOn)               !rewind(87)
!			rewind87 = 1
!			RewindCount = RewindCount + 1
!			write(95,*)' Rewind on line 503'
!			if(RewindCount.ge.2)then
!				write(*,*)' Every assemblage has been examined and no solution found'
!				write(*,*)' Skipping this P,T point and moving on to the next'
!				RewindCount = 0
!				go to 400
!				endif
!			endif
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
		if(EOFflag.eq.1)then
			write(*,*)' Every assemblage has been examined and no solution found'
			write(*,*)' Skipping this P,T point and moving on to the next'
			write(*,*)'T,P ',TC4,PB4
			go to 400
			endif		
		xPhCo = xPhCoLast			! Reset to last good compositions
		if(iParagonite.gt.0)then
		!	xPhCo(iParagonite,1) = 0.05
		!	xPhCo(iParagonite,2) = 0.95
			endif		
		go to 410
		endif
425	continue
!	If here, then we have found a solution using this assemblage

!	Check to see if any phases have negative moles - 
!		if so then find the largest negative moles and remove that phase
	maxNegativeMoles = 0.0d0
	maxNegativeMolesPhase = 0
	do 430 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	if(mp0(k)*atomNorm(k,1).lt.maxNegativeMoles)then
!		A phase has negative moles - record and loop back
		maxNegativeMoles = mp0(k)*atomNorm(k,1)		! we are storing moles/atom (not per molecule)
		maxNegativeMolesPhase = kCur
		endif
430	continue
	if(maxNegativeMolesPhase.ne.0)then
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
!			if(rewind87.eq.0)then
!				Call RewindTheFile(87,logFileOn)           !rewind(87)
!				RewindCount = RewindCount + 1
!				if(RewindCount.ge.2)then
!					write(*,*)' Every assemblage has been examined and no solution found'
!					write(*,*)' Skipping this P,T point and moving on to the next'
!					RewindCount = 0
!					go to 400
!					endif
!				rewind87 = 1
!				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			if(EOFflag.eq.1)then
				write(*,*)' Every assemblage has been examined and no solution found'
				write(*,*)' Skipping this P,T point and moving on to the next'
				write(*,*)'T,P ',TC4,PB4
				go to 400
				endif		
			xPhCo = xPhCoLast			! Reset to last good compositions
			if(iParagonite.gt.0)then
			!	xPhCo(iParagonite,1) = 0.05
			!	xPhCo(iParagonite,2) = 0.95
				endif		
			write(95,*)'Missed>20 (Negative moles)... loop back'
			go to 410
			endif
		iremove = asmCurrent(maxNegativeMolesPhase)
		mp0(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		do 432 j = maxNegativeMolesPhase,numPh
		asmCurrent(j) = asmCurrent(j+1)
432		continue
		numCurrent = numCurrent - 1
		write(95,*)'Removed phase= ',iremove,'asm ',(MinRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		countTheMisses = countTheMisses + 1
		go to 410
		endif


	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	if here, then no phase has negative moles
!	Check to see if any phases are below the tangent plane
435	continue
	maxgDifference = -1.0d-2
!	If a phase is more than 0.01 J below the tangent, then flag it.
	maxgDifferencePhase = 0
	do 440 k = 1,numPhMIF
	if(k.eq.iExclude)go to 440			! skip the excluded phase even if it is below the tangent
	if(gDifferenceAU(k).lt.maxgDifference)then			
		maxgDifference = gDifferenceAU(k)
		maxgDifferencePhase = k
		endif
440	continue
	if(maxgDifferencePhase.ne.0)then		
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
!			if(rewind87.eq.0)then
!				Call RewindTheFile(87,logFileOn)           !rewind(87)
!				rewind87 = 1
!				RewindCount = RewindCount + 1
!				if(RewindCount.ge.2)then
!					write(*,*)' Every assemblage has been examined and no solution found'
!					write(*,*)' Skipping this P,T point and moving on to the next'
!					RewindCount = 0
!					go to 400
!					endif
!				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			if(EOFflag.eq.1)then
				write(*,*)' Every assemblage has been examined and no solution found'
				write(*,*)' Skipping this P,T point and moving on to the next'
				write(*,*)'T,P ',TC4,PB4
				go to 400
				endif		
			xPhCo = xPhCoLast			! Reset to last good compositions
			if(iParagonite.gt.0)then
			!	xPhCo(iParagonite,1) = 0.05
			!	xPhCo(iParagonite,2) = 0.95
				endif		
			write(95,*)'Missed>20 (maxGdifference)... loop back'
			go to 410
			endif


!		A phase is below the tangent plane - add to the assemblage
!		This can work fine so long as the new assemblage is not nearly univariant.
!		Two cases identified are
!			1) When kspar enters assemblage at the expense of muscovite (univariant in KASH subsystem)
!			2) Ctoic + chl + Kyan + Sta (univariant in KFMASH subsystem)
!		The solution we try is
!			1) Try the calculation with the new phase to see if it works. If it does work then we don't have this degeneracy
!			2) Go through sub assemblages to see if we can find ones that work. The indicies for these are stored in array ASMCombo(numCurrent,indexes)
445		continue
!		if (numCurrent.lt.nc)then	! the variance < 2 - add this phase to the assemblage and continue
		numCurrent = numCurrent + 1
		asmCurrent(numCurrent) = maxgDifferencePhase
		write(95,*)'Phase added    ',(minRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
!		write(95,*)'Phase added    ',(asmCurrent(iok),iok=1,numCurrent)
		countTheMisses = countTheMisses + 1
		isubtend = 1
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)		! this one MUST subtend because the last one did and we added a phase
		izero = 0
		call ComputePTXM(TC4,PB4,izero)
		if(izero.eq.0)then
			go to 425			! this assemblage worked - go back and check for negative moles
			endif
!		If we get here, then the last assemblage didn't work (izero not = 0)
		write(95,*)'	Calculation failed in ComputePTXM (line 722 in Sub Tangent) - Work on subassemblages'
!		Begin trying sub assemblages
		do 451 k = 1,numCurrent
		asmCurrentStore(k) = asmCurrent(k)
451		continue
		numCurrentStore = numCurrent		! this is the number of phases in the assemblage that didn't work
!		The current assemblage has failed after adding a phase that sat below the tangent
!		This might be the result of a degeneracy (e.g. And + Sil but not on the reaction curve)
!		The trouble is, we don't know what phase to remove
!		So we will try every sub assemblage starting with numCurrent-1 phases
452		continue						! loop back to here if a subassemblage has another phase below tangent. Try all over
		do 480 icombo = numCurrent-1,1,-1	! Loop on all
!		Call Combinatorial to get an array (Combos) that contains the subassemblages with icombo phases
		combocounter = 0
		call Combinatorial3(icombo,numCurrentStore,1,1)
		write(95,*)' icombo, numCurrentStore, comboCounter ',icombo, numCurrentStore, comboCounter
		
		do 450 i = 1,comboCounter		
		xPhCo = xPhCoLast				! Reset to last good compositions
!		if(iParagonite.gt.0)then
		!	xPhCo(iParagonite,1) = 0.05
		!	xPhCo(iParagonite,2) = 0.95
!			endif		
		numCurrent = icombo		          !the number of phases in the assemblage
		write(95,*)i,numCurrent,(Combos(i,k),k=1,numCurrent)
		do 455 k = 1,numCurrent
		asmCurrent(k) = asmCurrentStore(Combos(i,k))
455		continue
		write(95,*)'	 Subassemblage  ',numCurrent,(MinRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		isubtend = 1
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)
		if(isubtend.eq.0)then
			write(95,*)'		Assemblage does not subtend BC (line 605 in Sub Tangent) - get new asm'
			go to 450
			endif
		izero = 0
		call ComputePTXM(TC4,PB4,izero)
		if(izero.gt.0)then
			write(95,*)'		   Error in ComputePTXM (Line 609 of Sub Tangent) - get new asm'
			go to 450			! get the next assemblage
			endif
!		If we get to here, then we have found a subassemblage that works.
!		Check for negative moles, etc.
		do 460 kCur = 1,numCurrent
		k = asmCurrent(kCur)
		if(mp0(k)*atomNorm(k,1).le.-1.0d-10)then
!			A phase has negative moles - go check out next assemblage
			write(95,*)'		   Phase has negative moles = ',k
			go to 450
			endif
460		continue

!		if here, then no phase has negative moles
!		Check to see if any phases are below the tangent plane
		call AdjustTangentToAsm
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
		do 462 k = 1,numPhMIF
		if(gDifferenceAU(k).lt.-1.0d-4)then			
			write(95,*)'		   Phase lies below tangent = ',k
!	The code fails at times when a new phase (one not in the original asmCurrent) lies below the tangent
!		first check to see that this phase isn't in the original assemblage
			do 463 j = 1,numCurrentStore
			if(k.eq.asmCurrentStore(j))go to 450		! this phase is already being considered. Don't reset anything
463			continue
!	It appears that this is a new phase. Try adding that phase in and going through the asmCombos again
!	Make sure this isn't the excluded phase
			if(iExclude.eq.k)go to 450
			numCurrentStore = numCurrentStore + 1
			asmCurrentStore(numCurrentStore) = k
			numCurrent = numCurrentStore
			go to 452		! Start the subassemblage routine all over
			endif
462		continue
!		if we get here, then 
!			(a) we have a solution
!			(b) no phase has negative moles
!			(c) no phase sits below the tangent
		go to 475

450		continue

480	continue	! finished looping on all combos
!	If we get to here then we have gone through all of the subassemblages in the ASMCompo array and none has worked
!	Loop back and try again reading assemblages from the main Compo file
		xPhCo = xPhCoLast			! Reset to last good compositions
		if(iParagonite.gt.0)then
		!	xPhCo(iParagonite,1) = 0.05
		!	xPhCo(iParagonite,2) = 0.95
			endif		
		go to 410

		endif

475	continue
!	If we get to here, then no phases sit below the tangent and we should be done
	sum = 0.0d0
	call CalculateMode(sum)		! sum is the sum of the volumes of phases in cm^3
	call CalculateGSystem
!	It seems that some assemblages result in NaN for a number of variables.
!	We need to check this and cycle back if we have NaN rather than a valid solution
	if(ieee_is_nan(gSystem))then
	!	write(*,*)'NaN detected'
		write(95,*)' Nan detected'
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
		if(EOFflag.eq.1)then
			write(*,*)' Every assemblage has been examined and no solution found'
			write(*,*)' Skipping this P,T point and moving on to the next'
			write(*,*)'T,P ',TC4,PB4
			go to 400
			endif		
		xPhCo = xPhCoLast			! Reset to last good compositions
		go to 410
		endif
!	call WriteGoodStuff		! output G values and comps to unit 72
!	call WriteOutputAsm		! write results to output file 73
	call WriteAllStuff		! write results to output file 74
!	if we are holding the porosity constant, then here is where we remove sufficient H2O to maintain constant porosity
	if(constantPorosity.eq.1)then
!		check to see if the fluid is part of the current assemblage
		do 472 k = 1,numCurrent
		if(fluidIndex.eq.asmCurrent(k))then
			molesFluidNew = (porosity*(sum - vp0(fluidIndex))/(100.0d0 - porosity))/(vmol(fluidIndex)*10)
!			molesFluidToChange = (vp0(fluidIndex) - porosity)/vatTP(fluidIndex,1)
!			molesFluidToChange = (mp0(fluidIndex)*vmol(fluidIndex) - porosity)/vatTP(fluidIndex,1)
			molesFluidToChange = mp0(fluidIndex) - molesFluidNew
			if(molesFluidToChange.le.0.0d0)then
				molesFluidToChange = 0.0d0
				go to 473		! do not allow H2O to be added back into the rock
				endif
			do 470 i = 1,nc
			molesSyComp = 0.0d0
!			for pure H2O I really don't need this loop 471, but I will need it for garnet fractionation
			do 471 j = 1,numPhCo(fluidIndex)
			molesSyComp = molesSyComp + xPhCo(fluidIndex,j)*comp(fluidIndex,j,i)		
471			continue
			molesChange = molesSyComp*molesFluidToChange
			bulkCompMoles(i) = bulkCompMoles(i) - molesChange		!bulkCompositionMoles is the same as molesStart, but used elsewhere (I should consolidate these)
			molesStart(i) = bulkCompMoles(i)					!molesStart is the bulk comp we want to fit in C2ompute
!			wtChange = molesChange*molwt(i)/numCatInOxide(i)
!			bulkCompWt(i) = bulkCompWt(i) - wtChange
470			continue
			endif
472		continue
		endif
473	continue
	rewind87 = 0			! reset the rewind switch - this says we have not rewound unit 87 for the new assemblage
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	MPLast = MP0
	if(iParagonite.gt.0)then
	!	xPhCo(iParagonite,1) = 0.05
	!	xPhCo(iParagonite,2) = 0.95
		endif		
	if(newP.eq.1)then		! We have a solution for the initial T for the new pressure - save
		asmCurrentNewP = asmCurrent	! Store the last good assemblage
		numCurrentNewP = numCurrent
		xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
		MPNewP = MP0
		endif
400	continue			! loop back to get the next temperature
!	close (72)
	close (73)
	close (74)
	if(logFileOn.eq.1)close (95)
	write(*,*)'All done with calculations.....'
	write(*,*)' The stable assemblages are in the file:            ',tangentOutputFile
	write(*,*)' All phases (stable and metastable) are in the file:',trim(tangentOutputFile)//'.ALL'
	write(*,*)' Debugging information is in file:                  ',trim(tangentOutputFile)//'.log'
	write(*,*)' '
	write(*,*)' Run program MADPlotter3 to plot the MAD diagram'
	write(*,*)' '
	write(*,*)' '
	pause ' hit return to continue'

	endif

! -------  Calculate the stable assemblage at the specified P&T --------------
	if(iTan.eq.41)then
	write(*,*)' Specify T and P for calculations'
	read(*,*)TC,PB
	TC4 = TC
	Pb4 = PB
	logFile = trim(tangentOutputFile)//'.log'
	logFile = 'TempLogFile.log'
	open(95,FILE=logFile,status = 'UNKNOWN')
	write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
	call StableASM(logfileon,TC4,PB4,ierr)
	close(87)	! close combo file
	close(95)	! close log file
	if(ierr.eq.0)then
		write(*,*)'Stable assemblage:'
		call WriteAllStuffToTerminal(6)
		write(*,*)' Do you want to save this file as a MIF?'
		write(*,*)' 0 = no'
		write(*,*)' 1 = yes'
		read(*,*)i
		if(i.eq.1)then
			call SaveMIF(TC4,PB4)
			endif
		pause ' hit return to continue'
		else
		write(*,*)' No solution found'
		pause ' hit return to continue'
		endif
	go to 100
	endif
! -------  Calculate the stable assemblage at the specified P&T --------------
	if(iTan.eq.410)then
	write(*,*)' Specify T and P for calculations'
	read(*,*)TC,PB
	TC4 = TC
	Pb4 = PB
	logFile = trim(tangentOutputFile)//'.log'
	logFile = 'TempLogFile.log'
	open(95,FILE=logFile,status = 'UNKNOWN')
	write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
	call StableASM2(logfileon,TC4,PB4,ierr)
	close(87)	! close combo file
	close(95)	! close log file
	if(ierr.eq.0)then
		write(*,*)'Stable assemblage:'
		call WriteAllStuffToTerminal(6)
		write(*,*)' Do you want to save this file as a MIF?'
		write(*,*)' 0 = no'
		write(*,*)' 1 = yes'
		read(*,*)i
		if(i.eq.1)then
			call SaveMIF(TC4,PB4)
			endif
		pause ' hit return to continue'
		else
		write(*,*)' No solution found'
		pause ' hit return to continue'
		endif
	go to 100
	endif

! -------  Calculate the stable assemblage at the specified P&T for new bulk compositions --------------
	if(iTan.eq.411)then
	call PickSecondBulkComposition()


	write(*,*)' Specify T and P for calculations'
	read(*,*)TC,PB
	TC4 = TC
	Pb4 = PB
	logFile = trim(tangentOutputFile)//'.log'
	logFile = 'TempLogFile.log'
	open(95,FILE=logFile,status = 'UNKNOWN')

4115	continue
	read(25,*,end=4110)(bulkCompmoles(i),i=1,nc)
	write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
	call StableASM2(logfileon,TC4,PB4,ierr)
	go to 4115

4110	continue
	close(87)	! close combo file
	close(95)	! close log file
	close(25)	! close temporary moles file
	if(ierr.eq.0)then
		write(*,*)'Stable assemblage:'
		call WriteAllStuffToTerminal(6)
		write(*,*)' Do you want to save this file as a MIF?'
		write(*,*)' 0 = no'
		write(*,*)' 1 = yes'
		read(*,*)i
		if(i.eq.1)then
			call SaveMIF(TC4,PB4)
			endif
		pause ' hit return to continue'
		else
		write(*,*)' No solution found'
		pause ' hit return to continue'
		endif
	go to 100
	endif


! -------  Autocalc (new algorithm) --------------
	if(iTan.eq.42.or.iTan.eq.421)then

	TC = Tgrid(1)
	PB = Pgrid(1)
	tc4 = tc
	pb4 = pb
!	tanPhase = 0
	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	iDoingGrid = 1				! keep going through activity calculations if we are doing a grid
	call ZeroTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort
	logFile = trim(tangentOutputFile)//'.log'
!	logFile = 'TempLogFile.log'
	open(95,FILE=logFile,status = 'UNKNOWN')
	write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
	do 4212 i = 1,numPhMIF
4212	write(*,*)i,phName(i)
	write(*,*)'Is there any phase to exclude from equilibrim assemblages (for calculation of affinities)?'
	write(*,*)'Input a number from the list or 0 for none (i.e. include all phases)'
	read(*,*)iExclude
	if(Tgrid(1).lt.Tgrid(numPTpoints))then		! we must be going up temperature
		if(iGrid.eq.1.or.iGrid.eq.2)then	! we cannot be going isothermally
			call ChooseConstantPorosity()
			endif
		endif
!	stuffFile = trim(tangentOutputFile)//'.stuff'
!	open(72,FILE=stuffFile,status = 'UNKNOWN')
	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
	call WriteAllOutHeader(74,Tinc,Pinc)


	call StableASM(logfileon,TC4,PB4,ierr)
	if(ierr.eq.0)then
		write(*,*)'Initial stable assemblage found:',TC4,PB4
		call WriteAllStuff()
!		pause ' hit return to continue'
		else
		write(*,*)' No solution for initial assemblage found'
		pause ' hit return to continue'
		go to 100
		endif

	if(iGrid.eq.1.or.iGrid.eq.2)then
		if(iTan.eq.42)then
			Call AutoTan(logFileOn,iGrid,ierr)
			endif
		if(iTan.eq.421)then
			Call AutoTan2(logFileOn,iGrid,ierr)
			endif
		endif
	if(iGrid.eq.3)Call AutoTanIsothermal(logFileOn,iGrid,ierr)
	if(ierr.eq.1)then
		write(*,*)' Some error occurred.... hit return to continue'
		pause
		go to 100
		endif

	go to 100
	endif
! ------   Auto calc with fractionation ---  Auto calculate from "Find initial tangent" from assemblage-------------
	if(iTan.eq.44)then

	call SetUpFractionation()
	TC = Tgrid(1)
	PB = Pgrid(1)
	tc4 = tc
	pb4 = pb
	
!	tanPhase = 0
	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
!	call SetTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort
!	do 4413 i = 1,numPhMIF
!4413	write(*,*)i,phName(i)
!	write(*,*)'Is there any phase to exclude from equilibrim assemblages (for calculation of affinities)?'
!	write(*,*)'Input a number from the list or 0 for none (i.e. include all phases)'
!	read(*,*)iExclude
	iExclude = 0
!	stuffFile = trim(tangentOutputFile)//'.stuff'
!	open(72,FILE=stuffFile,status = 'UNKNOWN')
	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
	BCMolesFile = trim(tangentOutputFile)//'.BCMoles'
	open(75,FILE=BCMolesFile,status = 'UNKNOWN')
	BCWtFile = trim(tangentOutputFile)//'.BCWt'
	open(76,FILE=BCWtFile,status = 'UNKNOWN')
	FractFile = trim(tangentOutputFile)//'.FRACT'
	open(77,FILE=FRACTFile,status = 'UNKNOWN')
	MeltFile = trim(tangentOutputFile)//'.Melt'
	open(78,FILE=MeltFile,status = 'UNKNOWN')
	logFile = trim(tangentOutputFile)//'.log'
	open(95,FILE=logFile,status = 'UNKNOWN')

	call WriteAllOutHeader(74,Tinc,Pinc)
	call WriteBulkCompHeader(75,1)		! Bulk composition moles
	call WriteBulkCompHeader(76,2)		! Bulk composition weight
	call WriteBulkCompHeader(78,3)		! Melt extraction
	call WriteFractHeader(77)

	call StableASM(logfileon,TC4,PB4,ierr)
	if(ierr.eq.0)then
		write(*,*)'Initial stable assemblage found:',TC4,PB4
		call WriteAllStuff()
!		pause ' hit return to continue'
		else
		write(*,*)' No solution for initial assemblage found'
		pause ' hit return to continue'
		go to 100
		endif
	Call AutoTan(logFileOn,iGrid,ierr)
	if(ierr.eq.1)then
		write(*,*)' Some error occurred.... hit return to continue'
		pause
		go to 100
		endif

	close (75)
	close (76)
	close (77)
	close (78)


	go to 100
	endif

! ------   Auto calc with fractionation ---  Auto calculate from "Find initial tangent" from assemblage-------------
! -------   This is the OS model -- (use option 44 for equilibrium model)
! -------   You can use assemblage evolution by using the full MIF
! -------   If you want a constant assemblage (probably more realistic) then make a MIF with only the desired phases
	if(iTan.eq.45)then

	!call SetUpFractionation()-- this only works for equilibrium when we are actually growing phases. 
	!In OS model, we need to specify how much to grow at each step
	TC = Tgrid(1)
	PB = Pgrid(1)
	tc4 = tc
	pb4 = pb
	
!	tanPhase = 0
	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
!	call SetTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort

	do i = 1,numPhMIF
		write(*,*)i,phName(i)
		end do
	write(*,*)'Input number for garnet (for fractionation)'
	read(*,*)kGrt
	iExclude = kGrt
	write(*,*)'Input millimoles of garnet to make on every iteration (e.g. .1 or .01 etc)'
	read(*,*)mGrt
! 	write(*,*)' Do you want to fractionate garnet? 0 = no, 1 = yes'
! 	read(*,*)doesGarnetFractionate
	doesGarnetFractionate = 1
	mGrt = mGrt/1000.0d0		! convert from millimoles to moles
	mGrtTotal = 0.0d0

! 	do i = 1,numPhMIF
!  		write(*,*)i,phName(i)
! 		end do
! 	write(*,*)'Is there any phase to exclude from equilibrim assemblages (for calculation of affinities)?'
! 	write(*,*)'Input a number from the list or 0 for none (i.e. include all phases)'
! 	read(*,*)iExclude
!	iExclude = 0

	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
! 	BCMolesFile = trim(tangentOutputFile)//'.BCMoles'
! 	open(75,FILE=BCMolesFile,status = 'UNKNOWN')
! 	BCWtFile = trim(tangentOutputFile)//'.BCWt'
! 	open(76,FILE=BCWtFile,status = 'UNKNOWN')
! 	FractFile = trim(tangentOutputFile)//'.FRACT'
! 	open(77,FILE=FRACTFile,status = 'UNKNOWN')
! 	MeltFile = trim(tangentOutputFile)//'.Melt'
! 	open(78,FILE=MeltFile,status = 'UNKNOWN')
	logFile = trim(tangentOutputFile)//'.log'
	open(95,FILE=logFile,status = 'UNKNOWN')

	call WriteAllOutHeader(74,Tinc,Pinc)
! 	call WriteBulkCompHeader(75,1)		! Bulk composition moles
! 	call WriteBulkCompHeader(76,2)		! Bulk composition weight
! 	call WriteBulkCompHeader(78,3)		! Melt extraction
! 	call WriteFractHeader(77)

	call StableASM(logfileon,TC4,PB4,ierr)
	if(ierr.eq.0)then
		write(*,*)'Initial stable assemblage found:',TC4,PB4
		call WriteAllStuff()
!		pause ' hit return to continue'
		else
		write(*,*)' No solution for initial assemblage found'
		pause ' hit return to continue'
		go to 4599
		endif
	write(*,*)'iTP,   TC4,     Pb4,       gDifferenceAU(kGrt),        mGrtTotal,        radius(µm)'

	Call AutoTanOSFract(logFileOn,iGrid,ierr,kGrt,mGrt,MgrtTotal,Radius)
	if(ierr.eq.1)then
		write(*,*)' Some error occurred.... hit return to continue'
		pause
		go to 4599
		endif
4599	continue
	close (74)
! 	close (75)
! 	close (76)
! 	close (77)
! 	close (78)
	close (95)


	go to 100
	endif



! --------------------------------
! ------Auto calculate from "Find initial tangent" from assemblage-------------
	if(iTan.eq.5)then

500   CONTINUE

	iCareifMIS = 1		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 1		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 1			! supresses error output in subroutine compute
	tanOut(1) = 1
!	call SetTanOut()

	write(*,*)'Input new T and P'
	read(*,*)TC,PB
	if(TC.eq.0.)go to 100
	tc4 = tc
	pb4 = pb
      	ALLX(1,1)=TC
      	ALLX(1,2)=PB
      	ALLX(2,1)=TC
      	ALLX(2,2)=PB
	write(12,*)' ===================================================================================='
	write(12,*)' ===================================================================================='
	write(12,*)tc4,pb4
	call CalcMuZero(izero)			! calculates for all phases at T&P - only needs to be done once at T&P

	write(*,*)'Quartz after calling CalcMuZero ',uzero(1,1),gattp(1,1)


	do 528 i = 1,numPhMIF
528	write(*,*)i,phName(i)
	write(*,*)'   '
	write(*,*)'Please provide your best guess for the initial assemblage.'
	write(*,*)'Input numbers for these phases (0 to end list)'
	write(*,*)' From the MIF:'
	write(*,*)AssembTitle
	numCurrent = 0
526	continue
	read(*,*)i
	if(i.eq.0)then
		go to 529
		else
		numCurrent = numCurrent + 1
		asmCurrent(numCurrent) = i
		endif
	go to 526
529	continue

	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
!      	SET DEFAULTS FOR MONITOR AND DELTAX
	SNstep = 1
	!ndep = 1
      	NSTEP=1
      	SNSTEP=1
	do 501 i=1,nvar-neq   ! should always be 2
	smon(i) = i
	mon(i) = i
	sdel(i)=0.d0
	deltax(i)= 0.00
501	continue
!   	Set up array IPOINT to contain pointers to non-monitor parameters
	J=0
	Do 510 i=1,NVAR
	DO 511 L=1,NVAR-NEQ
	IF(MON(L).eq.i)go to 510
511	continue
	J=J+1
	IPOINT(J) = I
510	CONTINUE
!      	IF(iLong(4).eq.1) write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)
      	write(12,*)' Monitors are:    ',(Mon(I),I=1,nvar-NEQ)
      	write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)

!      THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
	DO 503 I=1,nvar
!       ALLX(3,I) IS WHERE calculations STARTED-save
        ALLX(3,I)=ALLX(1,I)
503	CONTINUE

!	do 599 istep=1,snstep
!         Store new values for each independent variable (monitor parameter)
! 	This is where Neton's method increments composition, T or P
	DO 520 J=1,NVAR-NEQ
	i=mon(j)
	ALLX(1,i)=ALLX(1,i) + DELTAX(J)
520	continue
        call SetTPX
         
	write(*,561)TC,PB
561   	format (' T = ',f10.1,'  P = ',f10.1)
!      	write(*,*)' NSTEP = ',snstep,' Counter = ',istep

        IZERO=0
        call c2ompute(izero)
        IF(IZERO.gt.0)then
        	write(*,*)' izero error in call to c2ompute. izero =  ',izero
        	pause 'hit return to continue'
        	GO TO 100
		endif
!599   	continue

	
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	if(tanOut(7).eq.1)call ListTangentPlane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	if(tanOut(3).eq.1)call ListGDifference

	write(*,*)' quartz before PRINTT call ',gattp(1,1)
         CALL PRINTT(1)
	write(*,*)' quartz after PRINTT call ',gattp(1,1)
	call WriteAllStuffToTerminal(12)		! write results to output window

	pause 'Hit return to continue'	

	go to 100



	endif
! ------List assemblage information-------------
	if(iTan.eq.7)then
		call Printt(1)
		go to 100
		endif
! ------Call global routines -----------------------
	if(iTan.eq.8)then
		CALL GLOBAL(iTan)
		go to 100
		endif
! ------Adjust plotting axes-----------------------
	if(iTan.eq.9)then
	      CALL PLOTIN
	      go to 100
		endif
! ------Gibbs or Newton-----------------------
	if(iTan.eq.10)then
		call GibbsorNewton
	      	go to 100
		endif

! ------Calculate Affinity File-------------
	if(iTan.eq.11)then
!	This routine will generate an "ALL" file that gives the value of Gdiff (Gphase - Gtangent) for all phases in the input MIF file
!		but only compared to a single assemblage (not the most stable assemblage as is done in the  MMAD routine)
!		For example, assume the rock is chlorite + biotite + muscovite + quartz + plagioclase + H2O
!			 and calculate the affinity for the nucleation of every phase
!	It uses much of the same code as in option 4
!	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort
!	stuffFile = trim(adjustL(tangentOutputFile))//'.stuff'
!	open(72,FILE=stuffFile,status = 'UNKNOWN')
	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
	if(logFileOn.eq.1)then
		logFile = trim(adjustL(tangentOutputFile))//'.log'
		open(95,FILE=logFile,status = 'UNKNOWN')
		endif
!	call WriteOutputHeader(73)
	call WriteAllOutHeader(74,Tinc,Pinc)
	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the base assemblage from user
	do 1105 i = 1,numPhMIF
1105	write(*,*)i,phName(i)
	write(*,*)' Input numbers for phases in your base assemblage(0 to end list)'
	numCurrent = 0
1106	continue
	read(*,*)k
	if(k.eq.0)go to 1107
	numCurrent = numCurrent + 1
	asmCurrent(numCurrent) = k
	go to 1106
1107	continue
	numPh = numCurrent
	write(*,*)(asmCurrent(k),k=1,numCurrent)
	write(12,*)(asmCurrent(k),k=1,numCurrent)
	call PrinTT(0)
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	PB4old = 0
!	----------------- LOOP starts here -----------------
	do 1100 iTP = 1,numPTpoints
	newP = 0
	TC = Tgrid(iTP)
	PB = Pgrid(iTP)
	tc4 = tc
	pb4 = pb
	if(logFileOn.eq.1)then
		write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)
		endif
	if(Dabs(PB4-PB4old).gt.1.)then		! this should print every time P changes - just to show progress
		write(*,*)TC4,PB4
		PB4old = Pb4
		newP = 1
		endif
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	if(newP.eq.1)then
		asmCurrent = asmCurrentNewP		! set the assemblage to the low T, new P values
		numCurrent = numCurrentNewP
		xPhCo = xPhCoNewP
		else
		asmCurrent = asmCurrentLast		! set the assemblage to the last good one
		numCurrent = numCurrentLast
		endif
1150	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
		write(*,*)'This assemblage doesnot subtend the bulk composition - pick another assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 100
		endif
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 100
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	If we get to here, then no phases sit below the tangent and we should be done
	sum = 0.0d0
	call CalculateMode(sum)		! sum is the sum of the volumes of phases in cm^3
	call CalculateGSystem
!	call WriteGoodStuff		! output G values and comps to unit 72
!	call WriteOutputAsm		! write results to output file
	call WriteAllStuff		! write results to output file 74
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	if(newP.eq.1)then		! We have a solution for the initial T for the new pressure - save
		asmCurrentNewP = asmCurrent	! Store the last good assemblage
		numCurrentNewP = numCurrent
		xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
		endif
1100	continue			! loop back to get the next temperature
!	close (72)
	close (73)
	close (74)
	go to 100

	endif

! ------Grow garnet in overstepped environment-------------
	if(iTan.eq.12)then
!	This routine will calculate the zoning in a garnet at constant T&P after nucleating from a "base" assemblage via overstepping
!		The first part is exactly the same as option = 11 except T is set to the desired conditions
!	When the affinity of garnet reaches a prescribed value (e.g. 300 j/mole of O) then a garnet nucleates with the composition given
!		by the maximum G change (actually, I didn't implement this part - just picked the P&T to calculate
!	After nucleation garnet grows with the composition determined by the maximum G
!		after every growth increment, the bulk composition is fractionated by removing the grown garnet and
!			the tangent to the matrix is recomputed and a new garnet composition determined by the parallel tangent method
!		This continues until the growth increment is completed for this temperature.
!	The affinity of garnet is computed after every growth increment. 
!		It is expected that the affinity for producing garnet will decrease as Mn is removed the bulk composition
!			This will supress more nucleation until the temperature increases 
!	This routine will generate a modified "stuff" file that also includes the bulk composition 
!		and perhaps also the garnet zoning profile
!	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort
!	stuffFile = trim(adjustl(tangentOutputFile))//'.stuff'
!	open(72,FILE=stuffFile,status = 'UNKNOWN')
	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
	if(logFileOn.eq.1)then
		logFile = trim(adjustl(tangentOutputFile))//'.log'
		open(95,FILE=logFile,status = 'UNKNOWN')
		endif
!	call WriteOutputHeader(73)
	call WriteAllOutHeader(74,Tinc,Pinc)
	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the base assemblage from user
	do 1205 i = 1,numPhMIF
1205	write(*,*)i,phName(i)
	write(*,*)' Input numbers for phases in your base assemblage(0 to end list)'
	numCurrent = 0
1206	continue
	read(*,*)k
	if(k.eq.0)go to 1207
	numCurrent = numCurrent + 1
	asmCurrent(numCurrent) = k
	go to 1206
1207	continue
	write(*,*)'Input number for garnet'
	read(*,*)kGrt
	numPh = numCurrent
	write(*,*)(asmCurrent(k),k=1,numCurrent)
	write(12,*)(asmCurrent(k),k=1,numCurrent)
	call PrinTT(0)
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	PB4old = 0
	write(*,*)'Input T and P for calculations'
	read(*,*)Tgrid(1),Pgrid(1)
	numPTpoints = 1
!	----------------- LOOP starts here -----------------
!	do 1200 iTP = 1,numPTpoints
!	mGrt = 1.0d-5		! this is .1 millimole
	write(*,*)'Input millimoles of garnet to make on every iteration (e.g. .1 or .01 etc)'
	read(*,*)mGrt
	mGrt = mGrt/1000.0d0		! convert from millimoles to moles
	mGrtTotal = 0.0d0
	newP = 0
	iTP = 1
	TC = Tgrid(iTP)
	PB = Pgrid(iTP)
	tc4 = tc
	pb4 = pb
	if(logFileOn.eq.1)then
		write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)
		endif
	if(Dabs(PB4-PB4old).gt.1.)then		! this should print every time P changes - just to show progress
		write(12,*)TC4,PB4
		PB4old = Pb4
		newP = 1
		endif
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	if(newP.eq.1)then
		asmCurrent = asmCurrentNewP		! set the assemblage to the low T, new P values
		numCurrent = numCurrentNewP
		xPhCo = xPhCoNewP
		else
		asmCurrent = asmCurrentLast		! set the assemblage to the last good one
		numCurrent = numCurrentLast
		endif
1250	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
		write(*,*)'This assemblage doesnot subtend the bulk composition - pick another assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 100
		endif
	write(12,*)'                   BULK                                                                                Garnet'
	write(12,1260)(coname(i),i=1,nc),(coname(i),i=1,nc)
1260	format(40a15)

1200	continue		! loop back to here from below for incremental garnet growth
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 100
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	If we get to here, then no phases sit below the tangent and we should be done
	sum = 0.0d0
	call CalculateMode(sum)		! sum is the sum of the volumes of phases in cm^3
	call CalculateGSystem
	mp0(kGrt) = mGrt
!	call WriteGoodStuff		! output G values and comps to unit 72
!	call WriteOutputAsm		! write results to output file
	call WriteAllStuff		! write results to output file 74
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	if(gDifferenceAU(kGrt).ge.0)go to 1299		! we are now at equilibrium so exit
	if(newP.eq.1)then		! We have a solution for the initial T for the new pressure - save
		asmCurrentNewP = asmCurrent	! Store the last good assemblage
		numCurrentNewP = numCurrent
		xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
		endif
!	Now grow some garnet and subtract this from the bulk composition
!	I don't know how many moles to grow
!		this would depend on a number of factors in reality
!		For this routine, just make the mGrt small
! 	calculates the composition of each phase in moles and wt% units	
	k = kGrt
	do 1232 i = 1,NC
	PhComp(K,i) = 0.
	do 1230 j=1,numPhCo(K)
	PhComp(K,i) = PhComp(K,i) + xPhCo(k,j)*comp(k,j,i)
1230	continue
1232	continue
!	Now subtract from the bulk composition
	do 1240 i = 1,nc
	molesStart(i) = molesStart(i) - mGrt*phComp(k,i)		! molesStart is the bulk composition we fit in C2OMPUTE()
1240	continue
	write(12,1261)(molesStart(i),i=1,nc),(phComp(k,i),i=1,nc)
1261	format(40E15.7)
	write(*,*)gDifferenceAU(kGrt)
!	mGrt = mGrt*2.d0 				! see how this scaling works - the goal is to have increments small in the core then get bigger
!	pause ' Hit return to continue'
	go to 1200
!1200	continue			! loop back to get the next temperature
1299	continue			! end routine
!	close (72)
	close (73)
	close (74)
	write(*,*)'All done with calculations.....'
	write(*,*)' The stable assemblages are in the file:            ',tangentOutputFile
	write(*,*)' All phases (stable and metastable) are in the file:',trim(tangentOutputFile)//'.ALL'
	write(*,*)' Debugging information is in file:                  ',trim(tangentOutputFile)//'.log'
	write(*,*)' '
	write(*,*)' Run program MADPlotter3 to plot the MAD diagram'
	write(*,*)' '
	write(*,*)' '
	pause ' hit return to continue'
	go to 100

	endif
! ------Affinity node nucleation calculations-------------
	if(iTan.eq.13)then
!	This routine will calculate the affinity for nucleating garnet at a "node"
!		where a node is defined as 2, 3, 4, 5, etc grains intersecting at a grain boundary
!	The local bulk composition is comprised of the phases in contact only
!	The routine does this:
!	(1) Determine the compositions of phases in equilibrium given the whole-rock bulk composition and the P&T
!	(2) Check the combination file (contains all possible nodes) for ones that subtend the bulk composition
!		e.g. quartz + quartz doesn't but chlorite + quartz + ilmenite does
!	(3) Calculate the local bulk composition from the phases at the node
!	(4) Calculate the affinity for garnet nucleation
!	(5) Output the results in a "node" file that gives the value of Gdiff (Gphase - Gtangent) for all phases in the input MIF file
!	It uses much of the same code as in option 11
!	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort
	stuffFile = trim(adjustl(tangentOutputFile))//'.node'
	open(72,FILE=stuffFile,status = 'UNKNOWN')
	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
	if(logFileOn.eq.1)then
		logFile = trim(tangentOutputFile)//'.log'
		open(95,FILE=logFile,status = 'UNKNOWN')
		endif
!	call WriteOutputHeader(73)
	call WriteOutputHeader(72)
	call WriteAllOutHeader(74,Tinc,Pinc)
	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the base assemblage from user
	do 1305 i = 1,numPhMIF
1305	write(*,*)i,phName(i)
	write(*,*)' Input numbers for phases in your base assemblage(0 to end list)'
	numCurrent = 0
1306	continue
	read(*,*)k
	if(k.eq.0)go to 1307
	numCurrent = numCurrent + 1
	asmCurrent(numCurrent) = k
	go to 1306
1307	continue
	numPh = numCurrent
	write(*,*)(asmCurrent(k),k=1,numCurrent)
	write(12,*)(asmCurrent(k),k=1,numCurrent)
	call PrinTT(0)
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	PB4old = 0
	numPTpoints = 1
	Tgrid(1) = 600
	Pgrid(1) = 6000
!	----------------- LOOP starts here -----------------
	do 1300 iTP = 1,numPTpoints
	newP = 0
	TC = Tgrid(iTP)
	PB = Pgrid(iTP)
	tc4 = tc
	pb4 = pb
	if(logFileOn.eq.1)then
		write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)
		endif
	if(Dabs(PB4-PB4old).gt.1.)then		! this should print every time P changes - just to show progress
		write(12,*)TC4,PB4
		PB4old = Pb4
		newP = 1
		endif
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	if(newP.eq.1)then
		asmCurrent = asmCurrentNewP		! set the assemblage to the low T, new P values
		numCurrent = numCurrentNewP
		xPhCo = xPhCoNewP
		else
		asmCurrent = asmCurrentLast		! set the assemblage to the last good one
		numCurrent = numCurrentLast
		endif
1350	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
		write(*,*)'This assemblage doesnot subtend the bulk composition - pick another assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 100
		endif
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 100
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	If we get to here, then no phases sit below the tangent and we should be done
	sum = 0.0d0
	call CalculateMode(sum)		! sum is the sum of the volumes of phases in cm^3
	call CalculateGSystem
!	call WriteGoodStuff		! output G values and comps to unit 72
!	call WriteOutputAsm		! write results to output file
	call WriteAllStuff		! write results to output file 74
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	if(newP.eq.1)then		! We have a solution for the initial T for the new pressure - save
		asmCurrentNewP = asmCurrent	! Store the last good assemblage
		numCurrentNewP = numCurrent
		xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
		endif
1300	continue			! loop back to get the next temperature
	close (72)
	close (73)
	close (74)
	go to 100

	endif
! ------   Calculate garnet zoning along a EQUILIBRIUM P-T path with fractionation -------------
	if(iTan.eq.14)then
!	tanPhase = 0
	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort
!	stuffFile = trim(adjustl(tangentOutputFile))//'.stuff'
!	open(72,FILE=stuffFile,status = 'UNKNOWN')
	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
	if(logFileOn.eq.1)then
		logFile = trim(adjustl(tangentOutputFile))//'.log'
		open(95,FILE=logFile,status = 'UNKNOWN')
		endif
!	call WriteOutputHeader(73)
	call WriteAllOutHeader(74,Tinc,Pinc)
	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the first assemblage to check
!	get the base assemblage from user
	do 1405 i = 1,numPhMIF
1405	write(*,*)i,phName(i)
	write(*,*)'Input number for garnet (for growing/fractionation)'
	read(*,*)kGrt
	write(*,*)' Do you want to fractionate garnet? 0 = no, 1 = yes'
	read(*,*)doesGarnetFractionate
	Call RewindTheFile(87,logFileOn)
	rewind87 = 1
!	numCurrent = nc
!	read(87,*)(asmCurrent(i),i=1,numCurrent)
	!write(12,*)'GNA = 339'
	call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)		! get the initial assemblage to test
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	PB4old = 0
	if(logFileOn.eq.1)then
		write(95,*)' Starting TP loop'
		endif
	do 1400 iTP = 1,numPTpoints
	newP = 0
	TC = Tgrid(iTP)
	PB = Pgrid(iTP)
	tc4 = tc
	pb4 = pb
	if(logFileOn.eq.1)then
		write(95,1402)TC4,PB4,(bulkCompmoles(i),i=1,nc)
		endif
	if(Dabs(PB4-PB4old).gt.1.)then		! this should print every time P changes - just to show progress
!		If we are holding the porosity constant, then we need to reset the H2O content to have it in excess here
!		if(constantPorosity.eq.0)then
!			write(12,*)TC4,PB4
!			else
!			write(12,1402)TC4,PB4,(bulkCompmoles(i),i=1,nc)
1402			format(2F12.2,20F12.5)
!			if(iGrid.eq.1)then
			! reset the moles of H2O back to the starting value for the next row of calculations
			!	but only if we are doing a grid (not a P-T path)
!				bulkCompMoles = bulkCompMolesStart
!				molesStart = bulkCompMolesStart
!				endif
!			endif
		PB4old = Pb4
		newP = 1
		endif
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	if(izero.gt.0)then
		write(*,*)' Error calculating MU zero for this T&P... abort the run'
		write(95,*)' Error calculating MU zero for this T&P... abort the run'
		pause 'take a look and hit return when ready'
		go to 100
		endif
	if(newP.eq.1)then
!		If we are starting the next pressure, use the results from the previous P at the same T as a starting point
		asmCurrent = asmCurrentNewP		! set the assemblage to the low T, new P values
		numCurrent = numCurrentNewP
		xPhCo = xPhCoNewP
		else
!		Otherwise, use the last assemblage as the starting point
		asmCurrent = asmCurrentLast		! set the assemblage to the last good one
		numCurrent = numCurrentLast
		endif
	countTheMisses = 0
1450	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
!		write(*,*)'Problem in ComputePTXM - line 1020 of Tangent.f'
!		pause 'Hit return to continue'
		if(logFileOn.eq.1)then
			write(95,*)' Error in LoadCurrentAsm - get new asm'
			endif
		numCurrent = nc
		if(rewind87.eq.0)then
			Call RewindTheFile(87,logFileOn)    !rewind(87)
			rewind87 = 1
			endif
		!write(12,*)'GNA = 361'
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
! I may want to reset mineral compositions to starting values here
		xPhCo = xPhCoLast			! Reset to last good compositions
		go to 1450
		endif
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
!		write(*,*)'Problem in ComputePTXM - line 1020 of Tangent.f'
!		pause 'Hit return to continue'
		if(logFileOn.eq.1)then
			write(95,*)' Error in ComputePTXM - get new asm'
			endif
		numCurrent = nc
		if(rewind87.eq.0)then
			Call RewindTheFile(87,logFileOn)               !rewind(87)
			rewind87 = 1
			endif
		!write(12,*)'GNA = 378'
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
! 		I may want to reset mineral compositions to starting values here
		xPhCo = xPhCoLast			! Reset to last good compositions
		go to 1450
		endif

!	call StoreXinMIF()
1425	continue
!	Check to see if any phases have negative moles - 
!		if so then find the largest negative moles and remove that phase
	maxNegativeMoles = 0.0d0
	maxNegativeMolesPhase = 0
	do 1430 kCur = 1,numCurrent
	k = asmCurrent(kCur)
!	iMIF = asmCurrent(k)			! index of phase in MIF
	if(mp0(k)*atomNorm(k,1).lt.maxNegativeMoles)then
!		A phase has negative moles - remove from assemblage and loop back
		maxNegativeMoles = mp0(k)*atomNorm(k,1)		! we are storing moles/atom (not per molecule)
		maxNegativeMolesPhase = kCur
		endif
1430	continue
	if(maxNegativeMolesPhase.ne.0)then
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
			if(rewind87.eq.0)then
				Call RewindTheFile(87,logFileOn)           !rewind(87)
				rewind87 = 1
				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			xPhCo = xPhCoLast			! Reset to last good compositions
			go to 1450
			endif
!		molesPhaseMIF(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		mp0(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		do 1432 j = maxNegativeMolesPhase,numPh
		asmCurrent(j) = asmCurrent(j+1)
1432		continue
		numCurrent = numCurrent - 1
		if(logFileOn.eq.1)then
			write(95,*)(asmCurrent(iok),iok=1,numCurrent)
			endif
		countTheMisses = countTheMisses + 1
!		pause 'Take a look then hit return'
		go to 1450
		endif


	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	if here, then no phase has negative moles
!	Check to see if any phases are below the tangent plane
1435	continue
	maxgDifference = -1.0d-2
	maxgDifferencePhase = 0
	do 1440 k = 1,numPhMIF
	if(gDifferenceAU(k).lt.maxgDifference)then			
		maxgDifference = gDifferenceAU(k)
		maxgDifferencePhase = k
		endif
1440	continue
	if(maxgDifferencePhase.ne.0)then		
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
			if(rewind87.eq.0)then
				Call RewindTheFile(87,logFileOn)           !rewind(87)
				rewind87 = 1
				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			xPhCo = xPhCoLast			! Reset to last good compositions
			go to 1450
			endif


!		A phase is below the tangent plane - add to the assemblage

!		There is an issue when K-spar is added to a muscovite-bearing assemblage - it takes a long time to recover
!			this code will (maybe) fix that
		if(minRec(maxGDifferencePhase).eq.92)then  ! this is K-feldspar
			do 1441 kCur = 1,numCurrent
			k = asmCurrent(kCur)
			if(minRec(k).eq.17)then		! this is muscovite - take it out and add the K-feldspar
!				First try this new assemblage (sometimes it works)
				numCurrent = numCurrent + 1
				asmCurrent(numCurrent) = maxGDifferencePhase
				if(logFileOn.eq.1)then
					write(95,*)(asmCurrent(iok),iok=1,numCurrent)
					endif
				isubtend = 1
				call LoadCurrentAsm(TC4,PB4,izero,isubtend)
				izero = 0
				call ComputePTXM(TC4,PB4,izero)
				if(izero.eq.0)go to 1425	! ASM with Mx+Kfs didn't bomb - go check it out
!				izero must = 1. Remove muscovite and try again
				do 1442 j = kCur,numPh
				asmCurrent(j) = asmCurrent(j+1)
1442				continue
				numCurrent = numCurrent - 1
!				asmCurrent(numPh) = maxGDifferencePhase
				if(logFileOn.eq.1)then
					write(95,*)(asmCurrent(iok),iok=1,numCurrent)
					endif
				countTheMisses = countTheMisses + 1
				go to 1450	
				endif				
1441			continue
!			muscovite is not in the current assemblage - just keep going
			endif
!		if K-spar is not the new phase or if muscovite is not in the current assemblage, then just do regular recovery
1445		continue
		if (numCurrent.lt.nc)then	! the variance < 2 - add this phase to the assemblage and continue
			numCurrent = numCurrent + 1
			asmCurrent(numCurrent) = maxgDifferencePhase
			if(logFileOn.eq.1)then
				write(95,*)(asmCurrent(iok),iok=1,numCurrent)
				endif
			countTheMisses = countTheMisses + 1
			go to 1450
			else		! skip univariant stuff as an experiment
			numCurrent = nc
			if(rewind87.eq.0)then
				Call RewindTheFile(87,logFileOn)           !rewind(87)
				rewind87 = 1
				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			xPhCo = xPhCoLast			! Reset to last good compositions
!			write(12,*)' '
!			write(12,*)' Current assemblage:'
!			call ListCurrentAsm
			go to 1450
			endif
		endif
!	If we get to here, then no phases sit below the tangent and we should be done
	sum = 0.0d0
	call CalculateMode(sum)		! sum is the sum of the volumes of phases in cm^3
	call CalculateGSystem
!	call WriteGoodStuff		! output G values and comps to unit 72
!	call WriteOutputAsm		! write results to output file
	call WriteAllStuff		! write results to output file 74
	if(doesGarnetFractionate.eq.1)then
!		We are fractionating garnet so see if garnet is part of this assemblage
		do 1460 k = 1,numCurrent
		if(kGrt.eq.asmCurrent(k))then			! if true, then subtract garnet from the bulk composition
			mGrt = mp0(kGrt)
			do 1462 i = 1,NC
			PhComp(kGrt,i) = 0.
			do 1461 j=1,numPhCo(kGrt)
			PhComp(kGrt,i) = PhComp(kGrt,i) + xPhCo(kGrt,j)*comp(kGrt,j,i)
1461			continue
1462			continue
!			Now subtract from the bulk composition
			do 1465 i = 1,nc
			molesStart(i) = molesStart(i) - mGrt*phComp(kGrt,i)		! molesStart is the bulk composition we fit in C2OMPUTE()
1465			continue
			write(12,1466)(molesStart(i),i=1,nc),(phComp(kGrt,i),i=1,nc)
1466			format(40E15.7)
			pause 'hit return to continue'
			go to 1469 
			endif
1460		continue
1469		continue
		endif
	write(*,*)TC4,Pb4,gDifferenceAU(kGrt)

!	if we are holding the porosity constant, then here is where we remove sufficient H2O to maintain constant porosity
	if(constantPorosity.eq.1)then
!		check to see if the fluid is part of the current assemblage
		do 1472 k = 1,numCurrent
		if(fluidIndex.eq.asmCurrent(k))then
			molesFluidNew = (porosity*(sum - vp0(fluidIndex))/(100.0d0 - porosity))/(vmol(fluidIndex)*10)
!			molesFluidToChange = (vp0(fluidIndex) - porosity)/vatTP(fluidIndex,1)
!			molesFluidToChange = (mp0(fluidIndex)*vmol(fluidIndex) - porosity)/vatTP(fluidIndex,1)
			molesFluidToChange = mp0(fluidIndex) - molesFluidNew
			if(molesFluidToChange.le.0.0d0)then
				molesFluidToChange = 0.0d0
				go to 1473		! do not allow H2O to be added back into the rock
				endif
			do 1470 i = 1,nc
			molesSyComp = 0.0d0
!			for pure H2O I really don't need this loop 1471, but I will need it for garnet fractionation
			do 1471 j = 1,numPhCo(fluidIndex)
			molesSyComp = molesSyComp + xPhCo(fluidIndex,j)*comp(fluidIndex,j,i)		
1471			continue
			molesChange = molesSyComp*molesFluidToChange
			bulkCompMoles(i) = bulkCompMoles(i) - molesChange		!bulkCompositionMoles is the same as molesStart, but used elsewhere (I should consolidate these)
			molesStart(i) = bulkCompMoles(i)					!molesStart is the bulk comp we want to fit in C2ompute
!			wtChange = molesChange*molwt(i)/numCatInOxide(i)
!			bulkCompWt(i) = bulkCompWt(i) - wtChange
1470			continue
			endif
1472		continue
		endif
1473	continue
	rewind87 = 0			! reset the rewind switch - this says we have not rewound unit 87 for the new assemblage
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	if(newP.eq.1)then		! We have a solution for the initial T for the new pressure - save
		asmCurrentNewP = asmCurrent	! Store the last good assemblage
		numCurrentNewP = numCurrent
		xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
		endif
1400	continue			! loop back to get the next temperature
!	close (72)
	close (73)
	close (74)
	close (95)
	endif


! ------Check out activity calculation for HoPo dataset
	if(iTan.eq.15)then
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 1			! supresses error output in subroutine compute
	call SetTanOut()
	write(*,*)'Input new T and P'
	read(*,*)TC,PB
	if(TC.eq.0.)then
		go to 100
		endif
	tc4 = tc
	pb4 = pb
	write(12,*)' ===================================================================================='
	write(12,*)' ===================================================================================='
	write(12,*)tc4,pb4
	call CalcMuZero(izero)			! calculates for all phases at T&P - only needs to be done once at T&P
1500	continue
	do 1501 i = 1,numPhMIF
	write(*,*)i,phName(i)
1501	continue
	write(*,*)'Pick phases for the test assemblages (0 to end list)'
	numCurrent = 0
1510	continue
	if(numCurrent.eq.nc)then
		write(*,*)'You cannot pick more phases (v = 2)'
		go to 1511
		endif
	numCurrent = numCurrent + 1
	read(*,*)asmCurrent(numCurrent)
	if(asmCurrent(numCurrent).ne.0)go to 1510
	numCurrent = numCurrent - 1
1511	continue
	if(numCurrent.eq.0)go to 100			! bail out
1550	continue
	write(12,*)' This is the new assemblage...'
	call ListCurrentAsm
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(tanOut(4).eq.1)call thermodata   ! to check how things are  stored
	go to 1500
	endif
! ------Speed test------------------------------
	if(iTan.eq.16)then
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 1			! supresses error output in subroutine compute
	call SetTanOut()
1600	continue
	write(*,*)' Choose..'
	write(*,*)' 0 = return'
	write(*,*)' 1 = test ComputePTXM'
	write(*,*)' 2 = test ParallelToTangent'
	read(*,*)i
	if(i.eq.0)go to 100
	rewind(87)
	rewind87 = 1
	call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
	write(12,*)' '
	write(12,*)' Current assemblage:'
	call ListCurrentAsm
	write(*,*)'Input new T and P'
	read(*,*)TC,PB
	if(TC.eq.0.)then
		go to 100
		endif
	tc4 = tc
	pb4 = pb
	write(*,*)' How many iterations (1-lots)'
	read(*,*)j
	izero = 0
	isubtend = 0
	if(i.eq.1)then
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)
		do 1610 k = 1,j
		call ComputePTXM(TC4,PB4,izero)
1610		continue
		go to 1600
		endif
	if(i.eq.2)then
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)
		do 1620 k = 1,j
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
1620		continue
		go to 1600
		endif
	go to 1600
	endif

! ------   Calculate garnet zoning along a P-T path with fractionation -------------
	if(iTan.eq.17)then
!	This routine will calculate the zoning in a garnet along a PT path (set using option 3) after nucleating from an initial assemblage via overstepping
!	It is nearly exactly like option 12 except in option 12 the base assemblage can't change. 
!	In this routine, the stable paleoassemblage is computed as the bulk composition is modified by garnet fractionation
!	It was written specifically to examine the evolution of inclusion suites within garnet under isothermal, isobaric growth
!	tanPhase = 0
	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort
!	stuffFile = trim(adjustl(tangentOutputFile))//'.stuff'
!	open(72,FILE=stuffFile,status = 'UNKNOWN')
	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
	if(logFileOn.eq.1)then
		logFile = trim(adjustl(tangentOutputFile))//'.log'
		open(95,FILE=logFile,status = 'UNKNOWN')
		endif
!	call WriteOutputHeader(73)
	call WriteAllOutHeader(74,Tinc,Pinc)
	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the first assemblage to check
!	get the base assemblage from user
	do 1705 i = 1,numPhMIF
1705	write(*,*)i,phName(i)
	write(*,*)'Input number for garnet (for fractionation)'
	read(*,*)kGrt
	iExclude = kGrt
	Call RewindTheFile(87,logFileOn)
	rewind87 = 1
!	numCurrent = nc
!	read(87,*)(asmCurrent(i),i=1,numCurrent)
	!write(12,*)'GNA = 339'
	call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)		! get the initial assemblage to test
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	PB4old = 0
!	write(*,*)'Input T and P for calculations'
!	read(*,*)Tgrid(1),Pgrid(1)
!	numPTpoints = 1



	if(logFileOn.eq.1)then
		write(95,*)' Starting TP loop'
		endif

	write(*,*)'Input millimoles of garnet to make on every iteration (e.g. .1 or .01 etc)'
	read(*,*)mGrt
	write(*,*)' Do you want to fractionate garnet? 0 = no, 1 = yes'
	read(*,*)doesGarnetFractionate
	mGrt = mGrt/1000.0d0		! convert from millimoles to moles
	mGrtTotal = 0.0d0
!	iTP = 1
	do 1700 iTP = 1,numPTpoints
	newP = 0
	TC = Tgrid(iTP)
	PB = Pgrid(iTP)
	tc4 = tc
	pb4 = pb
	if(logFileOn.eq.1)then
		write(95,1702)TC4,PB4,(bulkCompmoles(i),i=1,nc)
		endif
	if(Dabs(PB4-PB4old).gt.1.)then		! this should print every time P changes - just to show progress
!		If we are holding the porosity constant, then we need to reset the H2O content to have it in excess here
!		if(constantPorosity.eq.0)then
!			write(12,*)TC4,PB4
!			else
!			write(12,1702)TC4,PB4,(bulkCompmoles(i),i=1,nc)
1702			format(2F12.2,20F12.5)
!			if(iGrid.eq.1)then
			! reset the moles of H2O back to the starting value for the next row of calculations
			!	but only if we are doing a grid (not a P-T path)
!				bulkCompMoles = bulkCompMolesStart
!				molesStart = bulkCompMolesStart
!				endif
!			endif
		PB4old = Pb4
		newP = 1
		endif
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	if(izero.gt.0)then
		write(*,*)' Error calculating MU zero for this T&P... abort the run'
		write(95,*)' Error calculating MU zero for this T&P... abort the run'
		pause 'take a look and hit return when ready'
		go to 100
		endif
	if(newP.eq.1)then
!		If we are starting the next pressure, use the results from the previous P at the same T as a starting point
		asmCurrent = asmCurrentNewP		! set the assemblage to the low T, new P values
		numCurrent = numCurrentNewP
		xPhCo = xPhCoNewP
		else
!		Otherwise, use the last assemblage as the starting point
		asmCurrent = asmCurrentLast		! set the assemblage to the last good one
		numCurrent = numCurrentLast
		endif
	countTheMisses = 0
!1700	continue

1750	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
!		write(*,*)'Problem in ComputePTXM - line 1020 of Tangent.f'
!		pause 'Hit return to continue'
		if(logFileOn.eq.1)then
			write(95,*)' Error in LoadCurrentAsm - get new asm'
			endif
		numCurrent = nc
		if(rewind87.eq.0)then
			Call RewindTheFile(87,logFileOn)    !rewind(87)
			rewind87 = 1
			endif
		!write(12,*)'GNA = 361'
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
! I may want to reset mineral compositions to starting values here
		xPhCo = xPhCoLast			! Reset to last good compositions
		go to 1750
		endif
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
!		write(*,*)'Problem in ComputePTXM - line 1020 of Tangent.f'
!		pause 'Hit return to continue'
		if(logFileOn.eq.1)then
			write(95,*)' Error in ComputePTXM - get new asm'
			endif
		numCurrent = nc
		if(rewind87.eq.0)then
			Call RewindTheFile(87,logFileOn)               !rewind(87)
			rewind87 = 1
			endif
		!write(12,*)'GNA = 378'
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
! 		I may want to reset mineral compositions to starting values here
		xPhCo = xPhCoLast			! Reset to last good compositions
		go to 1750
		endif

!	call StoreXinMIF()
1725	continue
!	Check to see if any phases have negative moles - 
!		if so then find the largest negative moles and remove that phase
	maxNegativeMoles = 0.0d0
	maxNegativeMolesPhase = 0
	do 1730 kCur = 1,numCurrent
	k = asmCurrent(kCur)
!	iMIF = asmCurrent(k)			! index of phase in MIF
	if(mp0(k)*atomNorm(k,1).lt.maxNegativeMoles)then
!		A phase has negative moles - remove from assemblage and loop back
		maxNegativeMoles = mp0(k)*atomNorm(k,1)		! we are storing moles/atom (not per molecule)
		maxNegativeMolesPhase = kCur
		endif
1730	continue
	if(maxNegativeMolesPhase.ne.0)then
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
			if(rewind87.eq.0)then
				Call RewindTheFile(87,logFileOn)           !rewind(87)
				rewind87 = 1
				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			xPhCo = xPhCoLast			! Reset to last good compositions
			go to 1750
			endif
!		molesPhaseMIF(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		mp0(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		do 1732 j = maxNegativeMolesPhase,numPh
		asmCurrent(j) = asmCurrent(j+1)
1732		continue
		numCurrent = numCurrent - 1
		if(logFileOn.eq.1)then
			write(95,*)(asmCurrent(iok),iok=1,numCurrent)
			endif
		countTheMisses = countTheMisses + 1
!		pause 'Take a look then hit return'
		go to 1750
		endif


	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	if here, then no phase has negative moles
!	Check to see if any phases are below the tangent plane
1735	continue
	maxgDifference = -1.0d-2
	maxgDifferencePhase = 0
	do 1740 k = 1,numPhMIF
	if(k.eq.iExclude)go to 1740			! skip the excluded phase even if it is below the tangent
	if(gDifferenceAU(k).lt.maxgDifference)then			
		maxgDifference = gDifferenceAU(k)
		maxgDifferencePhase = k
		endif
1740	continue
	if(maxgDifferencePhase.ne.0)then		
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
			if(rewind87.eq.0)then
				Call RewindTheFile(87,logFileOn)           !rewind(87)
				rewind87 = 1
				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			xPhCo = xPhCoLast			! Reset to last good compositions
			go to 1750
			endif


!		A phase is below the tangent plane - add to the assemblage

!		There is an issue when K-spar is added to a muscovite-bearing assemblage - it takes a long time to recover
!			this code will (maybe) fix that
		if(minRec(maxGDifferencePhase).eq.92)then  ! this is K-feldspar
			do 1741 kCur = 1,numCurrent
			k = asmCurrent(kCur)
			if(minRec(k).eq.17)then		! this is muscovite - take it out and add the K-feldspar
!				First try this new assemblage (sometimes it works)
				numCurrent = numCurrent + 1
				asmCurrent(numCurrent) = maxGDifferencePhase
				if(logFileOn.eq.1)then
					write(95,*)(asmCurrent(iok),iok=1,numCurrent)
					endif
				isubtend = 1
				call LoadCurrentAsm(TC4,PB4,izero,isubtend)
				izero = 0
				call ComputePTXM(TC4,PB4,izero)
				if(izero.eq.0)go to 1725	! ASM with Mx+Kfs didn't bomb - go check it out
!				izero must = 1. Remove muscovite and try again
				do 1742 j = kCur,numPh
				asmCurrent(j) = asmCurrent(j+1)
1742				continue
				numCurrent = numCurrent - 1
!				asmCurrent(numPh) = maxGDifferencePhase
				if(logFileOn.eq.1)then
					write(95,*)(asmCurrent(iok),iok=1,numCurrent)
					endif
				countTheMisses = countTheMisses + 1
				go to 1750	
				endif				
1741			continue
!			muscovite is not in the current assemblage - just keep going
			endif
!		if K-spar is not the new phase or if muscovite is not in the current assemblage, then just do regular recovery
1745		continue
		if (numCurrent.lt.nc)then	! the variance < 2 - add this phase to the assemblage and continue
			numCurrent = numCurrent + 1
			asmCurrent(numCurrent) = maxgDifferencePhase
			if(logFileOn.eq.1)then
				write(95,*)(asmCurrent(iok),iok=1,numCurrent)
				endif
			countTheMisses = countTheMisses + 1
			go to 1750
			else		! skip univariant stuff as an experiment
			numCurrent = nc
			if(rewind87.eq.0)then
				Call RewindTheFile(87,logFileOn)           !rewind(87)
				rewind87 = 1
				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			xPhCo = xPhCoLast			! Reset to last good compositions
!			write(12,*)' '
!			write(12,*)' Current assemblage:'
!			call ListCurrentAsm
			go to 1750
			endif
		endif
!	If we get to here, then no phases sit below the tangent and we should be done
	sum = 0.0d0
	call CalculateMode(sum)		! sum is the sum of the volumes of phases in cm^3
	call CalculateGSystem
	mp0(kGrt) = mGrt
!	call WriteGoodStuff		! output G values and comps to unit 72
!	call WriteOutputAsm		! write results to output file
	call WriteAllStuff		! write results to output file 74

	if(doesGarnetFractionate.eq.1)then
	!	Now grow some garnet and subtract this from the bulk composition
	!	I don't know how many moles to grow
	!		this would depend on a number of factors in reality
	!		For this routine, just make the mGrt small
	! 	calculates the composition of each phase in moles and wt% units	
	!	do 1760 k = 1,numCurrent
	!	if(kGrt.eq.asmCurrent(k))then			! if true, then subtract garnet from the bulk composition
			!mGrt = mp0(kGrt)
			do 1762 i = 1,NC
			PhComp(kGrt,i) = 0.
			do 1761 j=1,numPhCo(kGrt)
			PhComp(kGrt,i) = PhComp(kGrt,i) + xPhCo(kGrt,j)*comp(kGrt,j,i)
1761			continue
1762			continue
!			Now subtract from the bulk composition
			do 1765 i = 1,nc
			molesStart(i) = molesStart(i) - mGrt*phComp(kGrt,i)		! molesStart is the bulk composition we fit in C2OMPUTE()
1765			continue
!			write(12,1766)(molesStart(i),i=1,nc),(phComp(kGrt,i),i=1,nc)
1766			format(40E15.7)
			write(95,1704)TC4,PB4,(molesStart(i)*1000,i=1,nc)
1704			format(2F12.2,20F12.5)
			write(95,1703)(molwt(i)*molesStart(i)/numCatInOxide(i),i=1,nc)
1703			format(24x,20F12.5)
			write(95,1703)(xPhCo(kGrt,j),j=1,numPhCo(kGrt))
!			pause 'hit return to continue'
			go to 1769 
!			endif
!1760		continue
1769		continue
		endif		! end garnet fractionation
	mGrtTotal = mGrtTotal + mGrt
	radius = (((3./(4.*3.14159))*(mGrtTotal*117.))**0.333333333)*10000.	! radius in microns
	write(*,*)iTP,TC4,Pb4,gDifferenceAU(kGrt),mGrtTotal,radius


!	if we are holding the porosity constant, then here is where we remove sufficient H2O to maintain constant porosity
	if(constantPorosity.eq.1)then
!		check to see if the fluid is part of the current assemblage
		do 1772 k = 1,numCurrent
		if(fluidIndex.eq.asmCurrent(k))then
			molesFluidNew = (porosity*(sum - vp0(fluidIndex))/(100.0d0 - porosity))/(vmol(fluidIndex)*10)
!			molesFluidToChange = (vp0(fluidIndex) - porosity)/vatTP(fluidIndex,1)
!			molesFluidToChange = (mp0(fluidIndex)*vmol(fluidIndex) - porosity)/vatTP(fluidIndex,1)
			molesFluidToChange = mp0(fluidIndex) - molesFluidNew
			if(molesFluidToChange.le.0.0d0)then
				molesFluidToChange = 0.0d0
				go to 1773		! do not allow H2O to be added back into the rock
				endif
			do 1770 i = 1,nc
			molesSyComp = 0.0d0
!			for pure H2O I really don't need this loop 1771, but I will need it for garnet fractionation
			do 1771 j = 1,numPhCo(fluidIndex)
			molesSyComp = molesSyComp + xPhCo(fluidIndex,j)*comp(fluidIndex,j,i)		
1771			continue
			molesChange = molesSyComp*molesFluidToChange
			bulkCompMoles(i) = bulkCompMoles(i) - molesChange		!bulkCompositionMoles is the same as molesStart, but used elsewhere (I should consolidate these)
			molesStart(i) = bulkCompMoles(i)					!molesStart is the bulk comp we want to fit in C2ompute
!			wtChange = molesChange*molwt(i)/numCatInOxide(i)
!			bulkCompWt(i) = bulkCompWt(i) - wtChange
1770			continue
			endif
1772		continue
		endif
1773	continue
	rewind87 = 0			! reset the rewind switch - this says we have not rewound unit 87 for the new assemblage
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	if(newP.eq.1)then		! We have a solution for the initial T for the new pressure - save
		asmCurrentNewP = asmCurrent	! Store the last good assemblage
		numCurrentNewP = numCurrent
		xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
		endif
!	Rather than stopping, just keep going until the end of the PT file
!	if(gDifferenceAU(kGrt).ge.0)go to 1799		! we are now at equilibrium so exit

1700	continue
!	go to 1700		! do next iteration
!	close (72)
1799	continue
	close (73)
	close (74)
	close (95)
	endif



! -------  Debug the tangent function
	if(iTan.eq.18)then
	do 1802 i = 1,numPhMIF
1802	write(*,*)i,phName(i)
	write(*,*)' Which phase?'
	read(*,*)k
	if(k.eq.0)go to 100
	xtemp(1) = xPhCo(k,1)
	xtemp(2) = xPhCo(k,2)

	xPhCo(k,1) = 0.0d0
	do 1810 i = 1,99
	xPhCo(k,1) = xPhCo(k,1) + 0.01d0	
	xPhCo(k,2) = 1.0d0 - xPhCo(k,1)	
	iOK = 0
	call dLnAdX(K,iOK)				! update the activities based on the current phase composition
	if(iok.eq.1)then
		write(*,*)'iOK = 1 .... abort'
		pause
		go to 100
		endif
!	YY(j) = -(uzeroDeltaAU(k,jjj) + Rjoules*TK*((lnAct(k,jjj)/atomNorm(k,jjj)) - (lnAct(k,jdep)/atomNorm(k,jdep))) 
!     &              - uOnTanDeltaAU(k,jjj))
	xvalue = lnAct(k,2)/atomNorm(k,2) - lnAct(k,1)/atomNorm(k,1)
	write(96,*)xPhCo(k,1),xvalue
1810	continue
	close(96)
	go to 100
	endif



! ------   Calculate garnet zoning along a P-T path with fractionation -------------
	if(iTan.eq.19)then
!	This routine will calculate the zoning in a garnet at constant T&P after nucleating from an initial assemblage via overstepping
!	It is nearly exactly like option 12 except in option 12 the base assemblage can't change. 
!	In this routine, the stable paleoassemblage is computed as the bulk composition is modified by garnet fractionation
!	It was written specifically to examine the evolution of inclusion suites within garnet under isothermal, isobaric growth
!	tanPhase = 0
	iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
	call OpenOutputFile(iok)
	if(iok.ne.0)go to 100			! abort
!	stuffFile = trim(adjustl(tangentOutputFile))//'.stuff'
!	open(72,FILE=stuffFile,status = 'UNKNOWN')
	AllFile = trim(tangentOutputFile)//'.All'
	open(74,FILE=AllFile,status = 'UNKNOWN')
	if(logFileOn.eq.1)then
		logFile = trim(adjustl(tangentOutputFile))//'.log'
		open(95,FILE=logFile,status = 'UNKNOWN')
		endif
!	call WriteOutputHeader(73)
	call WriteAllOutHeader(74,Tinc,Pinc)
	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the first assemblage to check
!	get the base assemblage from user
	do 1905 i = 1,numPhMIF
1905	write(*,*)i,phName(i)
	write(*,*)'Input number for garnet (for fractionation)'
	read(*,*)kGrt
	iExclude = kGrt
	Call RewindTheFile(87,logFileOn)
	rewind87 = 1
!	numCurrent = nc
!	read(87,*)(asmCurrent(i),i=1,numCurrent)
	!write(12,*)'GNA = 339'
	call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)		! get the initial assemblage to test
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	PB4old = 0
	write(*,*)'Input T and P for calculations'
	read(*,*)Tgrid(1),Pgrid(1)
	numPTpoints = 1



	if(logFileOn.eq.1)then
		write(95,*)' Starting TP loop'
		endif

!	do 1900 iTP = 1,numPTpoints

	write(*,*)'Input millimoles of garnet to make on every iteration (e.g. .1 or .01 etc)'
	read(*,*)mGrt
	mGrt = mGrt/1000.0d0		! convert from millimoles to moles
	newP = 0
	iTP = 1
	TC = Tgrid(iTP)
	PB = Pgrid(iTP)
	tc4 = tc
	pb4 = pb
	if(logFileOn.eq.1)then
		write(95,1902)TC4,PB4,(bulkCompmoles(i),i=1,nc)
		endif
	if(Dabs(PB4-PB4old).gt.1.)then		! this should print every time P changes - just to show progress
!		If we are holding the porosity constant, then we need to reset the H2O content to have it in excess here
!		if(constantPorosity.eq.0)then
!			write(12,*)TC4,PB4
!			else
!			write(12,1902)TC4,PB4,(bulkCompmoles(i),i=1,nc)
1902			format(2F12.2,20F12.5)
!			if(iGrid.eq.1)then
			! reset the moles of H2O back to the starting value for the next row of calculations
			!	but only if we are doing a grid (not a P-T path)
!				bulkCompMoles = bulkCompMolesStart
!				molesStart = bulkCompMolesStart
!				endif
!			endif
		PB4old = Pb4
		newP = 1
		endif
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	if(izero.gt.0)then
		write(*,*)' Error calculating MU zero for this T&P... abort the run'
		write(95,*)' Error calculating MU zero for this T&P... abort the run'
		pause 'take a look and hit return when ready'
		go to 100
		endif
	if(newP.eq.1)then
!		If we are starting the next pressure, use the results from the previous P at the same T as a starting point
		asmCurrent = asmCurrentNewP		! set the assemblage to the low T, new P values
		numCurrent = numCurrentNewP
		xPhCo = xPhCoNewP
		else
!		Otherwise, use the last assemblage as the starting point
		asmCurrent = asmCurrentLast		! set the assemblage to the last good one
		numCurrent = numCurrentLast
		endif
	countTheMisses = 0
1900	continue

1950	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
!		write(*,*)'Problem in ComputePTXM - line 1020 of Tangent.f'
!		pause 'Hit return to continue'
		if(logFileOn.eq.1)then
			write(95,*)' Error in LoadCurrentAsm - get new asm'
			endif
		numCurrent = nc
		if(rewind87.eq.0)then
			Call RewindTheFile(87,logFileOn)    !rewind(87)
			rewind87 = 1
			endif
		!write(12,*)'GNA = 361'
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
! I may want to reset mineral compositions to starting values here
		xPhCo = xPhCoLast			! Reset to last good compositions
		go to 1950
		endif
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
!		write(*,*)'Problem in ComputePTXM - line 1020 of Tangent.f'
!		pause 'Hit return to continue'
		if(logFileOn.eq.1)then
			write(95,*)' Error in ComputePTXM - get new asm'
			endif
		numCurrent = nc
		if(rewind87.eq.0)then
			Call RewindTheFile(87,logFileOn)               !rewind(87)
			rewind87 = 1
			endif
		!write(12,*)'GNA = 378'
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
! 		I may want to reset mineral compositions to starting values here
		xPhCo = xPhCoLast			! Reset to last good compositions
		go to 1950
		endif

!	call StoreXinMIF()
1925	continue
!	Check to see if any phases have negative moles - 
!		if so then find the largest negative moles and remove that phase
	maxNegativeMoles = 0.0d0
	maxNegativeMolesPhase = 0
	do 1930 kCur = 1,numCurrent
	k = asmCurrent(kCur)
!	iMIF = asmCurrent(k)			! index of phase in MIF
	if(mp0(k)*atomNorm(k,1).lt.maxNegativeMoles)then
!		A phase has negative moles - remove from assemblage and loop back
		maxNegativeMoles = mp0(k)*atomNorm(k,1)		! we are storing moles/atom (not per molecule)
		maxNegativeMolesPhase = kCur
		endif
1930	continue
	if(maxNegativeMolesPhase.ne.0)then
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
			if(rewind87.eq.0)then
				Call RewindTheFile(87,logFileOn)           !rewind(87)
				rewind87 = 1
				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			xPhCo = xPhCoLast			! Reset to last good compositions
			go to 1950
			endif
!		molesPhaseMIF(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		mp0(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		do 1932 j = maxNegativeMolesPhase,numPh
		asmCurrent(j) = asmCurrent(j+1)
1932		continue
		numCurrent = numCurrent - 1
		if(logFileOn.eq.1)then
			write(95,*)(asmCurrent(iok),iok=1,numCurrent)
			endif
		countTheMisses = countTheMisses + 1
!		pause 'Take a look then hit return'
		go to 1950
		endif


	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	if here, then no phase has negative moles
!	Check to see if any phases are below the tangent plane
1935	continue
	maxgDifference = -1.0d-2
	maxgDifferencePhase = 0
	do 1940 k = 1,numPhMIF
	if(k.eq.iExclude)go to 1940			! skip the excluded phase even if it is below the tangent
	if(gDifferenceAU(k).lt.maxgDifference)then			
		maxgDifference = gDifferenceAU(k)
		maxgDifferencePhase = k
		endif
1940	continue
	if(maxgDifferencePhase.ne.0)then		
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
			if(rewind87.eq.0)then
				Call RewindTheFile(87,logFileOn)           !rewind(87)
				rewind87 = 1
				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			xPhCo = xPhCoLast			! Reset to last good compositions
			go to 1950
			endif


!		A phase is below the tangent plane - add to the assemblage

!		There is an issue when K-spar is added to a muscovite-bearing assemblage - it takes a long time to recover
!			this code will (maybe) fix that
		if(minRec(maxGDifferencePhase).eq.92)then  ! this is K-feldspar
			do 1941 kCur = 1,numCurrent
			k = asmCurrent(kCur)
			if(minRec(k).eq.17)then		! this is muscovite - take it out and add the K-feldspar
!				First try this new assemblage (sometimes it works)
				numCurrent = numCurrent + 1
				asmCurrent(numCurrent) = maxGDifferencePhase
				if(logFileOn.eq.1)then
					write(95,*)(asmCurrent(iok),iok=1,numCurrent)
					endif
				isubtend = 1
				call LoadCurrentAsm(TC4,PB4,izero,isubtend)
				izero = 0
				call ComputePTXM(TC4,PB4,izero)
				if(izero.eq.0)go to 1925	! ASM with Mx+Kfs didn't bomb - go check it out
!				izero must = 1. Remove muscovite and try again
				do 1942 j = kCur,numPh
				asmCurrent(j) = asmCurrent(j+1)
1942				continue
				numCurrent = numCurrent - 1
!				asmCurrent(numPh) = maxGDifferencePhase
				if(logFileOn.eq.1)then
					write(95,*)(asmCurrent(iok),iok=1,numCurrent)
					endif
				countTheMisses = countTheMisses + 1
				go to 1950	
				endif				
1941			continue
!			muscovite is not in the current assemblage - just keep going
			endif
!		if K-spar is not the new phase or if muscovite is not in the current assemblage, then just do regular recovery
1945		continue
		if (numCurrent.lt.nc)then	! the variance < 2 - add this phase to the assemblage and continue
			numCurrent = numCurrent + 1
			asmCurrent(numCurrent) = maxgDifferencePhase
			if(logFileOn.eq.1)then
				write(95,*)(asmCurrent(iok),iok=1,numCurrent)
				endif
			countTheMisses = countTheMisses + 1
			go to 1950
			else		! skip univariant stuff as an experiment
			numCurrent = nc
			if(rewind87.eq.0)then
				Call RewindTheFile(87,logFileOn)           !rewind(87)
				rewind87 = 1
				endif
			call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
			xPhCo = xPhCoLast			! Reset to last good compositions
!			write(12,*)' '
!			write(12,*)' Current assemblage:'
!			call ListCurrentAsm
			go to 1950
			endif
		endif
!	If we get to here, then no phases sit below the tangent and we should be done
	sum = 0.0d0
	call CalculateMode(sum)		! sum is the sum of the volumes of phases in cm^3
	call CalculateGSystem
	mp0(kGrt) = mGrt
!	call WriteGoodStuff		! output G values and comps to unit 72
!	call WriteOutputAsm		! write results to output file
	call WriteAllStuff		! write results to output file 74
!	Now grow some garnet and subtract this from the bulk composition
!	I don't know how many moles to grow
!		this would depend on a number of factors in reality
!		For this routine, just make the mGrt small
! 	calculates the composition of each phase in moles and wt% units	
!	do 1960 k = 1,numCurrent
!	if(kGrt.eq.asmCurrent(k))then			! if true, then subtract garnet from the bulk composition
		!mGrt = mp0(kGrt)
		do 1962 i = 1,NC
		PhComp(kGrt,i) = 0.
		do 1961 j=1,numPhCo(kGrt)
		PhComp(kGrt,i) = PhComp(kGrt,i) + xPhCo(kGrt,j)*comp(kGrt,j,i)
1961		continue
1962		continue
!		Now subtract from the bulk composition
		do 1965 i = 1,nc
		molesStart(i) = molesStart(i) - mGrt*phComp(kGrt,i)		! molesStart is the bulk composition we fit in C2OMPUTE()
1965		continue
!		write(12,1966)(molesStart(i),i=1,nc),(phComp(kGrt,i),i=1,nc)
1966		format(40E15.7)
		write(95,1904)TC4,PB4,(molesStart(i)*1000,i=1,nc)
1904		format(2F12.2,20F12.5)
		write(95,1903)(molwt(i)*molesStart(i)/numCatInOxide(i),i=1,nc)
1903		format(24x,20F12.5)
		write(95,1903)(xPhCo(kGrt,j),j=1,numPhCo(kGrt))
!		pause 'hit return to continue'
		go to 1969 
!		endif
!1960	continue
1969	continue
	write(*,*)TC4,Pb4,gDifferenceAU(kGrt)

!	if we are holding the porosity constant, then here is where we remove sufficient H2O to maintain constant porosity
	if(constantPorosity.eq.1)then
!		check to see if the fluid is part of the current assemblage
		do 1972 k = 1,numCurrent
		if(fluidIndex.eq.asmCurrent(k))then
			molesFluidNew = (porosity*(sum - vp0(fluidIndex))/(100.0d0 - porosity))/(vmol(fluidIndex)*10)
!			molesFluidToChange = (vp0(fluidIndex) - porosity)/vatTP(fluidIndex,1)
!			molesFluidToChange = (mp0(fluidIndex)*vmol(fluidIndex) - porosity)/vatTP(fluidIndex,1)
			molesFluidToChange = mp0(fluidIndex) - molesFluidNew
			if(molesFluidToChange.le.0.0d0)then
				molesFluidToChange = 0.0d0
				go to 1973		! do not allow H2O to be added back into the rock
				endif
			do 1970 i = 1,nc
			molesSyComp = 0.0d0
!			for pure H2O I really don't need this loop 1971, but I will need it for garnet fractionation
			do 1971 j = 1,numPhCo(fluidIndex)
			molesSyComp = molesSyComp + xPhCo(fluidIndex,j)*comp(fluidIndex,j,i)		
1971			continue
			molesChange = molesSyComp*molesFluidToChange
			bulkCompMoles(i) = bulkCompMoles(i) - molesChange		!bulkCompositionMoles is the same as molesStart, but used elsewhere (I should consolidate these)
			molesStart(i) = bulkCompMoles(i)					!molesStart is the bulk comp we want to fit in C2ompute
!			wtChange = molesChange*molwt(i)/numCatInOxide(i)
!			bulkCompWt(i) = bulkCompWt(i) - wtChange
1970			continue
			endif
1972		continue
		endif
1973	continue
	rewind87 = 0			! reset the rewind switch - this says we have not rewound unit 87 for the new assemblage
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	if(newP.eq.1)then		! We have a solution for the initial T for the new pressure - save
		asmCurrentNewP = asmCurrent	! Store the last good assemblage
		numCurrentNewP = numCurrent
		xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
		endif
	if(gDifferenceAU(kGrt).ge.0)go to 1999		! we are now at equilibrium so exit
	go to 1900		! do next iteration
!	close (72)
1999	continue
	close (73)
	close (74)
	close (95)
	endif

! ------ Jacob's ladder routine (Jacobian) ------------------------------
	if(iTan.eq.20)then
	tc4 = 700.
	pb4 = 8000
2001	continue
	do 2005 k = 1,numPhMif
	asmCurrent(k) = k
2005	continue
	iLong(7) = 1
	Call SetInitialTangent
	iLong(7) = 0
	write(12,*)tc4,pb4
	call CalcMuZero(izero)			! calculates for all phases at T&P - only needs to be done once at T&P

!	Assign values to the tangent plane based on the results linear combinations of phase components and GatTP values 
	ii = npMIF	! this is the index for the beginning of the independent reactions in array ARX
!	write(12,*)' ii, npMIF,nc',ii,npmif
	do 2010 i = 1,nc	! loop for every system component on tangent plane
	j = 0
	sum = 0.0d0
	do 2012 k = 1,numPhMIF
!	k = asmCurrent(kCur)		! should be every phase in order
	do 2014 jj = 1,numPhCo(k)
!	write(12,*)k,jj,gattp(k,jj)
	j = j + 1
	sum = sum - ARX(i+ii,j)*Gattp(k,jj)
2014	continue
2012	continue
!	write(12,*)'i, sum = ',i,sum
	tanPlane(i) = sum
2010	continue

	call CalculateGSystem()
	call ListTangentPlane()

!	tc4 = tc
!	pb4 = pb
!	Somewhere around here is where we will loop on T and P
!	numCurrent = 1
!	asmCurrent(numCurrent) = tanPhase	! The tangent phase is always the first phase in the list
	write(12,*)' ===================================================================================='
	write(12,*)' ===================================================================================='
	write(12,*)tc4,pb4
	call CalcMuZero(izero)			! calculates for all phases at T&P - only needs to be done once at T&P
		! needed when we loop on T and P
	if(izero.gt.0)then
		! try and recover
		endif
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
!	if(tanOut(7).eq.1)
	call CalculateGSystem()
	call ListTangentPlane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	if(tanOut(3).eq.1)
	call ListGDifference
!	pause 'pausing...'	
!	Find the phase that sits farthest below the tangent and lower the tangent to this phase
	maxgDifference = 1.0d20
	maxgDifferencePhase = 0
	do 2040 k = 1,numPhMIF
	if(gDifferenceAU(k).lt.maxgDifference)then			
		maxgDifference = gDifferenceAU(k)
		maxgDifferencePhase = k
		endif
2040	continue
	do 2045 i = 1,nc
	tanPlane(i) = tanPlane(i) + maxgDifference
2045	continue
	call CalculateGSystem()
	call ListTangentPlane
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	if(tanOut(3).eq.1)
	call ListGDifference
!	pause 'The first phase should now sit on the tangent plane...'	


!	Add this phase to the assemblage list
	numCurrent = 1		! the first time through numCurrent will = 2
	asmCurrent(numCurrent) = maxgDifferencePhase
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	loop to here when we add a new phase to the asemblage
2050	continue
	write(12,*)' '
	write(12,*)' Current assemblage:'
	call ListCurrentAsm()
!	pause 'Look at Current assemblage ... pausing...'	
!++++++++++++++++++++++++++++++++++++++++++++++++++
!	Check this assemblage to see if it is the equilibrium one
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
		write(12,*)' Assemblage does not subtend bulk composition - line 2032 of Tangent.f'
!		pause 'Hit return to continue'
!		xPhCo = xPhCoLast			! Reset to last good compositions
		call ListCurrentAsm
!		Keep iterating on Jacobian
		go to 2059
		endif
	izero = 0
	isubtend = 0
!	if(tanOut(6).eq.1)call PRINTT(1)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)' izero = 1: Problem in ComputePTXM - line 2045 of Tangent.f'
		xPhCo = xPhCoLast			! Reset to last good compositions
		call ListCurrentAsm
		pause 'Hit return to continue'
!		Keep iterating on Jacobian
		go to 2059
		endif

!	If here, then the calculation of the assemblage succeeded
	call Printt(1)
	call AdjustTangentToAsm
	call ListTangentPlane
!	Check to see if any phases have negative moles - 
	write(*,*)numCurrent,(mp0(asmCurrent(k)),k=1,numCurrent)
	do 2065 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	if(mp0(k).lt.0.0d0)then
!		A phase has negative moles
!		remove this phase from the assemblage
		do 2064 j = kCur,numCurrent
		asmCurrent(j) = asmCurrent(j+1)
2064		continue
		numCurrent = numCurrent - 1
		write(12,*)'Removed phase= ',k,'asm ',(MinRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		call ListCurrentAsm
		call ListGDifference
		write(12,*)' Looping to 2059 to continue Jacobian iteration'
		go to 2059			
		endif
2065	continue
!	The calculation succeeded and no phases have negative moles.
!	This should be the stable assemblage
!		I should maybe check with ParallelToTangent one final time
	call Printt(1)
	call AdjustTangentToAsm
	call ListTangentPlane
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	call ListGDifference
	call CalculateGSystem()
	call PrinTT(1)
	Write(*,*) ' this should be the stable assemblage'
	write(*,*)'Input new T and P'
	read(*,*)TC4,PB4
	if(TC.eq.0.)then
		go to 100
		endif
	TC = TC4
	Pb = Pb4
	go to 2001

!++++++++++++++++++++++++++++++++++++++++++++++++++
	niterate = 0
!++++++++++++++++++++++++++++++++++++++++++++++++++
	ilong(7) = 2
2059	continue
	niterate = niterate + 1
!	we should probably check the bulk composition at this point but...
!	Build the Jacobian
	call BuildJacobian2(izero)
	ilong(7) = 0
!	pause 'Check out Jacobian calculation'

!++++++++++++++++++++++++++++++++++++++++++++++++++
	call CalculateGSystem()
!	call ListTangentPlane
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	write(12,*)'Iteration = ',niterate,'Gsys= ',Gsystem,'Gdiff= ',(gDifferenceAU(i),i=1,numPhMIF)
	do 2060 k = 1,numPhMIF
	do 2061 kcur = 1,numCurrent
	if(k.eq.asmCurrent(kcur))go to 2060
2061	continue
	if(gDifferenceAU(k).lt.0.0d0)then
!		Add this phase to the current assemblage
		numCurrent = numCurrent+1
		asmCurrent(numCurrent) = k
		if(numCurrent.eq.nc)then
!			This is divariant.... take a look
			call ListCurrentAsm()
			pause 'We now have a divariant assemblage .. pausing'
			endif
		call ListCurrentAsm()
		write(12,*) 'Looping to 2050 to check this new assemblage for stability'
		pause 'We just added a phase to the assemblage - pause'	
		iLong(7) = 2
		go to 2050		! Check this assemblage with COMPUTE for stability
		endif

2060	continue
!	pause 'Iterate on tangent plane...'	
	if(iLong(7).eq.2)then
		call ListGDifference
		write(12,*)' Continue Jacobian loop with current assemblage - go to 2059'
		!pause '  - pause'
		endif
	go to 2059

!++++++++++++++++++++++++++++++++++++++++++++++++++
	go to 100
	endif



!---------------------------------------------------------------------
! ------Calculate EBC-------------------------------------------------
	if(iTan.eq.21)then
!	This routine calculates the effective bulk composition to match the composition of the garnet core at nucleation
!	It starts with the starting bulk composition and iterates on Fe, Mn and Ca to match Xalm, Xsps and Xgrs
!	It uses code from option 12
!	It requires 
!	(1) A starting bulk composition
!	(2) The "base" (paleo) assemblage
!	(3) The P and T of nucleation
!	(4) The composition of the garnet core
!	The output is the bulk composition with modified Fe, Mn and Ca concentrations
!	The code uses Newton's method to find the solution
!	Make all statement labels start with 21000
!----------------------
!	This code does exactly the same thing as option 21 except:
!	Option 21 solves for Fe, Mn and Ca	(i.e. 3 element in the EBC)
!	Option 22 solves for only Fe and Mn	(i.e. 2 elements in the EBC)
!	Option 23 solves for Fe, Mn, Ca and Na	(i.e. 4 elements in the EBC)
!	Note that even though the code says "Fe, Mn, Ca etc." one can actually solve for any 2, 3 or 4 elements in the EBC)
!----------------------
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
!	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the base assemblage from user
	do 21205 i = 1,numPhMIF
21205	write(*,*)i,phName(i)
	write(*,*)' Input numbers for phases in your base assemblage(0 to end list)'
	numCurrent = 0
21206	continue
	read(*,*)k
	if(k.eq.0)go to 21207
	numCurrent = numCurrent + 1
	asmCurrent(numCurrent) = k
	go to 21206
21207	continue
	write(*,*)'Input number for garnet'
	read(*,*)kGrt
	numPh = numCurrent
	call PrinTT(0)
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	PB4old = 0
	write(*,*)'Input Xalm, Xsps, Xgrs in the core in mole fractions'
	read(*,*)XalmCore,XspsCore,XgrsCore

21501	continue
	write(*,*)'Input T and P (bars) for calculations (i.e. the T and P of nucleation)(0 to exit)'
	read(*,*)TC,PB
	if(TC.eq.0)go to 100

	tanOut(8) = 0

21500	continue	! Loop here for picking new indicies
	write(*,*)
	do i = 1,nc
	write(*,21600)i,coname(i)
	molesStart(i) = bulkCompMolesStart(i)		! reset the bulk composition to the initial composition
	end do
21600	format(I5,2x,A8)
	write(*,*)' Input indicies for Fe, Mn, and Ca in bulk rock'
	read(*,*)iFe,iMn,iCa
	if(iFe.eq.0)go to 21501
	xPhCo = xPhCoInitial		! Reset to the initial compositions
	write(*,*)' Input 0 = no debugging output; 1 = output each iteration'
	read(*,*)tanOut(8)
	numPTpoints = 1
!	----------------- LOOP starts here -----------------
	newP = 0
	iTP = 1
	tc4 = tc
	pb4 = pb
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)

21250	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
		write(*,*)'This assemblage doesnot subtend the bulk composition - pick another assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 100
		endif
	if(tanOut(8).eq.1)then
		write(12,*)' ---------------------------'
		write(12,*)' Modifying components ',coname(iFe),coname(iMn),coname(iCa)
		write(12,*)'                   BULK                                                                                Garnet'
		write(12,21260)(coname(i),i=1,nc)
21260		format(3x,40a10)
		endif
!===========================================================
!	mincFD = 0.01d0		! increment of moles for finite difference calculations
	mincFD = 0.0001d0		! increment of moles for finite difference calculations
	write(*,*)'Input value for mincFSD  (around 0.001)'
	read(*,*)mincFD
	iCount = 0		! counter for iterations
21200	continue		! loop back to here from below finite difference convergence
	icount = icount + 1
	if(icount.gt.1000)then
		write(*,*)'iCount iterations > 1000 .... something must be wrong'
		pause 'Hit return to continue'
		go to 21500
		endif

!	Get the initial compositions of phases using the current bulk composition
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#1 Calc base: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		write(12,21261)(molesStart(i),i=1,nc)
		pause 'Hit return to continue'
		go to 21500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalc = XphCo(kGrt,2)
	XspsCalc = XphCo(kGrt,3)
	XgrsCalc = XphCo(kGrt,4)

	if(tanOut(8).eq.1)then		
		Write(12,*)' iCount = ',iCount
		write(12,21261)(molesStart(i),i=1,nc)
21261		format(20F10.5)
		write(12,*)'     Xprp      Xalm      Xsps      Xgrs'
		write(12,21261)(1.-XalmCalc - XspsCalc - XgrsCalc),XalmCalc,XspsCalc,XgrsCalc
		endif

!	write(12,21261)(molesStart(i),i=1,nc),XalmCalc,XspsCalc,XgrsCalc

	tol = 0.0001d0
!	Check for convergence
	if(dabs(XalmCore - XalmCalc).gt.tol)go to 21201
	if(dabs(XspsCore - XspsCalc).gt.tol)go to 21201
	if(dabs(XgrsCore - XgrsCalc).gt.tol)go to 21201

!	if we get to here, then we have converged

!	Calculate the Wt% composition from the moles
!	sum = 0.0d0 	! we don't normalize because that changes rock volume we are modeling
	write(12,*)'         '
	write(12,*)' We have converged'
	write(12,*)' Modifying components ',coname(iFe),coname(iMn),coname(iCa)
	call PrinTT(0)
	write(12,21260)(coname(i),i=1,nc)
	write(12,21261)(molesStart(i),i=1,nc)	
	do 21208 i = 1,nc
	wtpct(i) = molesStart(i)*molwt(i)/NumCatInOxide(i)
21208	continue
	write(12,21261)(WtPct(i),i=1,nc)
	write(12,*)'     Xprp      Xalm      Xsps      Xgrs'
	write(12,21261)(1.-XalmCalc - XspsCalc - XgrsCalc),XalmCalc,XspsCalc,XgrsCalc
	write(12,*)' '
	write(*,*)' We have converged -- solution is in the output window'
	pause 'Hit return to continue'
	go to 21500

21201	continue
!	calculate finite difference derivatives
!	Derivative for Fe in bulk composition
	molesStart(iFe) = molesStart(iFe) + mincFD
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#2 Calc Fe FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 21500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcPlus = XphCo(kGrt,2)
	XspsCalcPlus = XphCo(kGrt,3)
	XgrsCalcPlus = XphCo(kGrt,4)
!	Calculate derivatives
	DalmDMolesFe = (XalmCalc - XalmCalcPlus)/mincFD
	DspsDMolesFe = (XspsCalc - XspsCalcPlus)/mincFD
	DgrsDMolesFe = (XgrsCalc - XgrsCalcPlus)/mincFD
	molesStart(iFe) = molesStart(iFe) - mincFD		! set back to starting composition
	
!	Derivative for Mn in bulk composition
	molesStart(iMn) = molesStart(iMn) + mincFD
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#3 Calc Mn FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 21500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcPlus = XphCo(kGrt,2)
	XspsCalcPlus = XphCo(kGrt,3)
	XgrsCalcPlus = XphCo(kGrt,4)
!	Calculate derivatives
	DalmDMolesMn = (XalmCalc - XalmCalcPlus)/mincFD
	DspsDMolesMn = (XspsCalc - XspsCalcPlus)/mincFD
	DgrsDMolesMn = (XgrsCalc - XgrsCalcPlus)/mincFD
	molesStart(iMn) = molesStart(iMn) - mincFD		! set back to starting composition

!	Derivative for Ca in bulk composition
	molesStart(iCa) = molesStart(iCa) + mincFD
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#4 Calc Ca FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 21500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcPlus = XphCo(kGrt,2)
	XspsCalcPlus = XphCo(kGrt,3)
	XgrsCalcPlus = XphCo(kGrt,4)
!	Calculate derivatives
	DalmDMolesCa = (XalmCalc - XalmCalcPlus)/mincFD
	DspsDMolesCa = (XspsCalc - XspsCalcPlus)/mincFD
	DgrsDMolesCa = (XgrsCalc - XgrsCalcPlus)/mincFD
	molesStart(iCa) = molesStart(iCa) - mincFD		! set back to starting composition

!	Set up to solve the equations
!	-(XalmCore - XalmCalc) = 0 = (XalmCore - XalmCalc) + (dXalm/dFe)*ÆFe + (dXalm/dMn)*ÆMn + (dXalm/dCa)*ÆCa
!	-(XspsCore - XspsCalc) = 0 = (XspsCore - XspsCalc) + (dXsps/dFe)*ÆFe + (dXsps/dMn)*ÆMn + (dXsps/dCa)*ÆCa
!	-(XgrsCore - XgrsCalc) = 0 = (XgrsCore - XgrsCalc) + (dXgrs/dFe)*ÆFe + (dXgrs/dMn)*ÆMn + (dXgrs/dCa)*ÆCa

	A(1,4) = -(XalmCore - XalmCalc)
	A(2,4) = -(XspsCore - XspsCalc)
	A(3,4) = -(XgrsCore - XgrsCalc)

	A(1,1) = DalmDMolesFe
	A(1,2) = DalmDMolesMn
	A(1,3) = DalmDMolesCa
	A(2,1) = DspsDMolesFe
	A(2,2) = DspsDMolesMn
	A(2,3) = DspsDMolesCa
	A(3,1) = DgrsDMolesFe
	A(3,2) = DgrsDMolesMn
	A(3,3) = DgrsDMolesCa

	if(tanOut(8).eq.1)then		
		write(12,*)' '
		write(12,*)'A matrix'
		do 21040 i = 1,numEq
		write(12,21080)(A(i,j),j=1,numEq+1)
21080		format(50F15.5)
21040		continue
		endif		

!     Find solution
	numEq = 3
	nsolve = 4
	DETA=0.D0
	IER=0
	CALL REDUCE (numEq,NSOLVE,DETA,IER)
	IF (IER.EQ.1) then
	      WRITE(12,*)' ************ ERROR **************************'
	      WRITE(12,*)' Matrix failed to invert in SUBROUTINE REDUCE'
	      write(12,*)' We are in Subroutine ParallalToTangent in file TangentRoutines.f '
	      izero=1
	      write(12,*) 'Hit return to continue...'
	      pause 'Hit return to continue...'
	      go to 21500
	      endif

	if(tanOut(8).eq.1)then		
		write(12,*)'Solution xx'
		write(12,21080)(xx(j,1),j=1,numEq)
		write(12,*)'Determinant = ',DETA
		endif		


!	Calculate new bulk composition
!	First check to see that no moles go negative
!	Damp = 1.
	Damp = .1D0		!  maximum damp
21041	continue
	if(tanOut(8).eq.1)then		
		write(12,*)' Damp = ',damp
		endif
	if((damp*xx(1,1)+molesStart(iFe)).le.0.0d0.or.	&
	   (damp*xx(2,1)+molesStart(iMn)).le.0.0d0.or.	&
	   (damp*xx(3,1)+molesStart(iCa)).le.0.0d0) then
	   	Damp = Damp/2.0d0
		if(damp.lt.1.0D-8)then
			write(*,*)' Damp = ',damp
			pause 'something must be wrong ... damp is too small'
			go to 21500
			endif
	   	go to 21041
	   	endif
	molesStart(iFe) = molesStart(iFe) + xx(1,1)*damp
	molesStart(iMn) = molesStart(iMn) + xx(2,1)*damp
	molesStart(iCa) = molesStart(iCa) + xx(3,1)*damp
	if(tanOut(8).eq.1)then		
		write(12,21261)(molesStart(i),i=1,nc)
		endif
	
	go to 21200

	endif
!===================================
!---------------------------------------------------------------------
! ------Calculate EBC-------------------------------------------------
	if(iTan.eq.22)then
!	This routine calculates the effective bulk composition to match the composition of the garnet core at nucleation
!	It starts with the starting bulk composition and iterates on Fe, Mn and Ca to match Xalm, Xsps and Xgrs
!	It uses code from option 12
!	It requires 
!	(1) A starting bulk composition
!	(2) The "base" (paleo) assemblage
!	(3) The P and T of nucleation
!	(4) The composition of the garnet core
!	The output is the bulk composition with modified Fe, Mn and Ca concentrations
!	The code uses Newton's method to find the solution
!	Make all statement labels start with 22000
!----------------------
!	This code does exactly the same thing as option 21 except:
!	Option 21 solves for Fe, Mn and Ca	(i.e. 3 element in the EBC)
!	Option 22 solves for only Fe and Mn	(i.e. 2 elements in the EBC)
!	Option 23 solves for Fe, Mn, Ca and Na	(i.e. 4 elements in the EBC)
!	Note that even though the code says "Fe, Mn, Ca etc." one can actually solve for any 2, 3 or 4 elements in the EBC)
!----------------------

	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
!	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the base assemblage from user
	do 22205 i = 1,numPhMIF
22205	write(*,*)i,phName(i)
	write(*,*)' Input numbers for phases in your base assemblage(0 to end list)'
	numCurrent = 0
22206	continue
	read(*,*)k
	if(k.eq.0)go to 22207
	numCurrent = numCurrent + 1
	asmCurrent(numCurrent) = k
	go to 22206
22207	continue
	write(*,*)'Input number for garnet'
	read(*,*)kGrt
	numPh = numCurrent
	call PrinTT(0)
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	PB4old = 0
	write(*,*)'Input Xalm, Xsps in the core in mole fractions'
	read(*,*)XalmCore,XspsCore

22501	continue
	write(*,*)'Input T and P (bars) for calculations (i.e. the T and P of nucleation)(0 to exit)'
	read(*,*)TC,PB
	if(TC.eq.0)go to 100

	tanOut(8) = 0

22500	continue	! Loop here for picking new indicies
	write(*,*)
	do i = 1,nc
	write(*,22600)i,coname(i)
	molesStart(i) = bulkCompMolesStart(i)		! reset the bulk composition to the initial composition
	end do
22600	format(I5,2x,A8)
	write(*,*)' Input indicies for Fe and Mn in bulk rock'
	read(*,*)iFe,iMn
	if(iFe.eq.0)go to 22501
	xPhCo = xPhCoInitial		! Reset to the initial compositions
	write(*,*)' Input 0 = no debugging output; 1 = output each iteration'
	read(*,*)tanOut(8)
	numPTpoints = 1
!	----------------- LOOP starts here -----------------
	newP = 0
	iTP = 1
	tc4 = tc
	pb4 = pb
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)

22250	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
		write(*,*)'This assemblage doesnot subtend the bulk composition - pick another assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 100
		endif
	if(tanOut(8).eq.1)then
		write(12,*)' ---------------------------'
		write(12,*)' Modifying components ',coname(iFe),coname(iMn)
		write(12,*)'                   BULK                                                                                Garnet'
		write(12,22260)(coname(i),i=1,nc)
22260		format(3x,40a10)
		endif
!===========================================================
!	mincFD = 0.01d0		! increment of moles for finite difference calculations
	mincFD = 0.0001d0		! increment of moles for finite difference calculations
	write(*,*)'Input value for mincFSD  (around 0.001)'
	read(*,*)mincFD
	iCount = 0		! counter for iterations
22200	continue		! loop back to here from below finite difference convergence
	icount = icount + 1
	if(icount.gt.1000)then
		write(*,*)'iCount iterations > 1000 .... something must be wrong'
		pause 'Hit return to continue'
		go to 22500
		endif

!	Get the initial compositions of phases using the current bulk composition
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#1 Calc base: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		write(12,22261)(molesStart(i),i=1,nc)
		pause 'Hit return to continue'
		go to 22500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalc = XphCo(kGrt,2)
	XspsCalc = XphCo(kGrt,3)
!	XgrsCalc = XphCo(kGrt,4)

	if(tanOut(8).eq.1)then		
		Write(12,*)' iCount = ',iCount
		write(12,22261)(molesStart(i),i=1,nc)
22261		format(20F10.5)
		write(12,*)'     Xprp      Xalm      Xsps      Xgrs'
!		write(12,22261)(1.-XalmCalc - XspsCalc - XgrsCalc),XalmCalc,XspsCalc,XgrsCalc
		write(12,22261)(1.-XalmCalc - XspsCalc),XalmCalc,XspsCalc
		endif

!	write(12,22261)(molesStart(i),i=1,nc),XalmCalc,XspsCalc,XgrsCalc

	tol = 0.0001d0
!	Check for convergence
	if(dabs(XalmCore - XalmCalc).gt.tol)go to 22201
	if(dabs(XspsCore - XspsCalc).gt.tol)go to 22201
!	if(dabs(XgrsCore - XgrsCalc).gt.tol)go to 22201

!	if we get to here, then we have converged

!	Calculate the Wt% composition from the moles
!	sum = 0.0d0 	! we don't normalize because that changes rock volume we are modeling
	write(12,*)'         '
	write(12,*)' We have converged'
!	write(12,*)' Modifying components ',coname(iFe),coname(iMn),coname(iCa)
	write(12,*)' Modifying components ',coname(iFe),coname(iMn)
	call PrinTT(0)
	write(12,22260)(coname(i),i=1,nc)
	write(12,22261)(molesStart(i),i=1,nc)	
	do 22208 i = 1,nc
	wtpct(i) = molesStart(i)*molwt(i)/NumCatInOxide(i)
22208	continue
	write(12,22261)(WtPct(i),i=1,nc)
	write(12,*)'     Xprp      Xalm      Xsps      Xgrs'
!	write(12,22261)(1.-XalmCalc - XspsCalc - XgrsCalc),XalmCalc,XspsCalc,XgrsCalc
	write(12,22261)(1.-XalmCalc - XspsCalc),XalmCalc,XspsCalc
	write(12,*)' '
	write(*,*)' We have converged -- solution is in the output window'
	pause 'Hit return to continue'
	go to 22500

22201	continue
!	calculate finite difference derivatives
!	Derivative for Fe in bulk composition
	molesStart(iFe) = molesStart(iFe) + mincFD
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#2 Calc Fe FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 22500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcPlus = XphCo(kGrt,2)
	XspsCalcPlus = XphCo(kGrt,3)
!	XgrsCalcPlus = XphCo(kGrt,4)
	molesStart(iFe) = molesStart(iFe) - 2.0d0*mincFD		! calculate minus
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#2 Calc Fe FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 22500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcMinus = XphCo(kGrt,2)
	XspsCalcMinus = XphCo(kGrt,3)
!	XgrsCalcMinus = XphCo(kGrt,4)
	molesStart(iFe) = molesStart(iFe) + mincFD		! set back to starting composition
!	Calculate derivatives
	DalmDMolesFe = (XalmCalcPlus - XalmCalcMinus)/(2.0d0*mincFD)
	DspsDMolesFe = (XspsCalcPlus - XspsCalcMinus)/(2.0d0*mincFD)
!	DgrsDMolesFe = (XgrsCalcPlus - XgrsCalcMinus)/(2.0d0*mincFD)
	
!	Derivative for Mn in bulk composition
	molesStart(iMn) = molesStart(iMn) + mincFD
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#3 Calc Mn FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 22500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcPlus = XphCo(kGrt,2)
	XspsCalcPlus = XphCo(kGrt,3)
!	XgrsCalcPlus = XphCo(kGrt,4)
	molesStart(iMn) = molesStart(iMn) - 2.0d0*mincFD		! set back to starting composition
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#3 Calc Mn FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 22500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcMinus = XphCo(kGrt,2)
	XspsCalcMinus = XphCo(kGrt,3)
!	XgrsCalcMinus = XphCo(kGrt,4)
	molesStart(iMn) = molesStart(iMn) + mincFD		! set back to starting composition
!	Calculate derivatives
	DalmDMolesMn = (XalmCalcPlus - XalmCalcMinus)/(2.0d0*mincFD)
	DspsDMolesMn = (XspsCalcPlus - XspsCalcMinus)/(2.0d0*mincFD)
!	DgrsDMolesMn = (XgrsCalcPlus - XgrsCalcMinus)/(2.0d0*mincFD)

!	Derivative for Ca in bulk composition
!	molesStart(iCa) = molesStart(iCa) + mincFD
!	izero = 0
!	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
!	call ComputePTXM(TC4,PB4,izero)
!	if(izero.gt.0)then
!		write(*,*)'#4 Calc Ca FD: Problem in ComputePTXM trying to calculate the base assemblage'
!		write(*,*)TC4,PB4
!		write(*,*)(asmCurrent(k),k=1,numCurrent)
!		pause 'Hit return to continue'
!		go to 22500
!		endif
!	call AdjustTangentToAsm
!	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
!	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
!	XalmCalcPlus = XphCo(kGrt,2)
!	XspsCalcPlus = XphCo(kGrt,3)
!!	XgrsCalcPlus = XphCo(kGrt,4)
!	molesStart(iCa) = molesStart(iCa) - 2.0D0*mincFD		! set back to starting composition
!	izero = 0
!	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
!	call ComputePTXM(TC4,PB4,izero)
!	if(izero.gt.0)then
!		write(*,*)'#4 Calc Ca FD: Problem in ComputePTXM trying to calculate the base assemblage'
!		write(*,*)TC4,PB4
!		write(*,*)(asmCurrent(k),k=1,numCurrent)
!		pause 'Hit return to continue'
!		go to 22500
!	call AdjustTangentToAsm
!	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
!	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
!	XalmCalcMinus = XphCo(kGrt,2)
!	XspsCalcMinus = XphCo(kGrt,3)
!	XgrsCalcMinus = XphCo(kGrt,4)
!	molesStart(iCa) = molesStart(iCa) + mincFD		! set back to starting composition
!	Calculate derivatives
!	DalmDMolesCa = (XalmCalcPlus - XalmCalcMinus)/(2.0d0*mincFD)
!	DspsDMolesCa = (XspsCalcPlus - XspsCalcMinus)/(2.0d0*mincFD)
!	DgrsDMolesCa = (XgrsCalcPlus - XgrsCalcMinus)/(2.0d0*mincFD)

!	Set up to solve the equations
!	-(XalmCore - XalmCalc) = 0 = (XalmCore - XalmCalc) + (dXalm/dFe)*ÆFe + (dXalm/dMn)*ÆMn + (dXalm/dCa)*ÆCa
!	-(XspsCore - XspsCalc) = 0 = (XspsCore - XspsCalc) + (dXsps/dFe)*ÆFe + (dXsps/dMn)*ÆMn + (dXsps/dCa)*ÆCa
!	-(XgrsCore - XgrsCalc) = 0 = (XgrsCore - XgrsCalc) + (dXgrs/dFe)*ÆFe + (dXgrs/dMn)*ÆMn + (dXgrs/dCa)*ÆCa

	A(1,3) = -(XalmCore - XalmCalc)
	A(2,3) = -(XspsCore - XspsCalc)
!	A(3,4) = -(XgrsCore - XgrsCalc)

	A(1,1) = -DalmDMolesFe
	A(1,2) = -DalmDMolesMn
!	A(1,3) = -DalmDMolesCa
	A(2,1) = -DspsDMolesFe
	A(2,2) = -DspsDMolesMn
!	A(2,3) = -DspsDMolesCa
!	A(3,1) = -DgrsDMolesFe
!	A(3,2) = -DgrsDMolesMn
!	A(3,3) = -DgrsDMolesCa

	numEq = 2
	if(tanOut(8).eq.1)then		
		write(12,*)' '
		write(12,*)'A matrix'
		do 22040 i = 1,numEq
		write(12,22080)(A(i,j),j=1,numEq+1)
22080		format(50F15.5)
22040		continue
		endif		

!     Find solution
	nsolve = 3
	DETA=0.D0
	IER=0
	CALL REDUCE (numEq,NSOLVE,DETA,IER)
	IF (IER.EQ.1) then
	      WRITE(12,*)' ************ ERROR **************************'
	      WRITE(12,*)' Matrix failed to invert in SUBROUTINE REDUCE'
	      write(12,*)' We are in Subroutine ParallalToTangent in file TangentRoutines.f '
	      izero=1
	      write(12,*) 'Hit return to continue...'
	      pause 'Hit return to continue...'
	      go to 22500
	      endif

	if(tanOut(8).eq.1)then		
		write(12,*)'Solution xx'
		write(12,22080)(xx(j,1),j=1,numEq)
		write(12,*)'Determinant = ',DETA
		endif		


!	Calculate new bulk composition
!	First check to see that no moles go negative
!	Damp = 1.
	Damp = .1D0		!  maximum damp
22041	continue
	if(tanOut(8).eq.1)then		
		write(12,*)' Damp = ',damp
		endif
!	if((damp*xx(1,1)+molesStart(iFe)).le.0.0d0.or.	&
!	   (damp*xx(2,1)+molesStart(iMn)).le.0.0d0.or.	&
!	   (damp*xx(3,1)+molesStart(iCa)).le.0.0d0) then
	if((damp*xx(1,1)+molesStart(iFe)).le.0.0d0.or.	&
	   (damp*xx(2,1)+molesStart(iMn)).le.0.0d0) then
	   	Damp = Damp/2.0d0
		if(damp.lt.1.0D-8)then
			write(*,*)' Damp = ',damp
			pause 'something must be wrong ... damp is too small'
			go to 22500
			endif
	   	go to 22041
	   	endif
	molesStart(iFe) = molesStart(iFe) + xx(1,1)*damp
	molesStart(iMn) = molesStart(iMn) + xx(2,1)*damp
!	molesStart(iCa) = molesStart(iCa) + xx(3,1)*damp
	if(tanOut(8).eq.1)then		
		write(12,22261)(molesStart(i),i=1,nc)
		endif
	
	go to 22200

	endif
!===================================
!---------------------------------------------------------------------
! ------Calculate EBC-------------------------------------------------
	if(iTan.eq.23)then
!	This routine calculates the effective bulk composition to match the composition of the garnet core at nucleation
!	It starts with the starting bulk composition and iterates on Fe, Mn and Ca to match Xalm, Xsps and Xgrs
!	It uses code from option 12
!	It requires 
!	(1) A starting bulk composition
!	(2) The "base" (paleo) assemblage
!	(3) The P and T of nucleation
!	(4) The composition of the garnet core
!	The output is the bulk composition with modified Fe, Mn and Ca concentrations
!	The code uses Newton's method to find the solution
!	Make all statement labels start with 23000
!----------------------
!	This code does exactly the same thing as option 21 except:
!	Option 21 solves for Fe, Mn and Ca	(i.e. 3 element in the EBC)
!	Option 22 solves for only Fe and Mn	(i.e. 2 elements in the EBC)
!	Option 23 solves for Fe, Mn, Ca and Na	(i.e. 4 elements in the EBC)
!	Note that even though the code says "Fe, Mn, Ca etc." one can actually solve for any 2, 3 or 4 elements in the EBC)
!----------------------
	iLong(1) = 0		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 0			! supresses error output in subroutine compute
	call ZeroTanOut()
!	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	get the base assemblage from user
	do 23205 i = 1,numPhMIF
23205	write(*,*)i,phName(i)
	write(*,*)' Input numbers for phases in your base assemblage(0 to end list)'
	numCurrent = 0
23206	continue
	read(*,*)k
	if(k.eq.0)go to 23207
	numCurrent = numCurrent + 1
	asmCurrent(numCurrent) = k
	go to 23206
23207	continue
	write(*,*)'Input number for garnet'
	read(*,*)kGrt
	write(*,*)'Input number for plagioclase'
	read(*,*)kPlg
	numPh = numCurrent
	call PrinTT(0)
	asmCurrentLast = asmCurrent	! Store the last good assemblage
	numCurrentLast = numCurrent
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	asmCurrentNewP = asmCurrent	! Store the last good assemblage
	numCurrentNewP = numCurrent
	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
	PB4old = 0
	write(*,*)'Input Xalm, Xsps, Xgrs in the core in mole fractions'
	read(*,*)XalmCore,XspsCore,XgrsCore
	write(*,*)'Input Xan in plagioclase in mole fractions'
	read(*,*)XanCore

23501	continue
	write(*,*)'Input T and P (bars) for calculations (i.e. the T and P of nucleation)(0 to exit)'
	read(*,*)TC,PB
	if(TC.eq.0)go to 100

	tanOut(8) = 0

23500	continue	! Loop here for picking new indicies
	write(*,*)
	do i = 1,nc
	write(*,23600)i,coname(i)
	molesStart(i) = bulkCompMolesStart(i)		! reset the bulk composition to the initial composition
	end do
23600	format(I5,2x,A8)
	write(*,*)' Input indicies for Fe, Mn, Ca and Na in bulk rock'
	read(*,*)iFe,iMn,iCa,iNa
	if(iFe.eq.0)go to 23501
	xPhCo = xPhCoInitial		! Reset to the initial compositions
	write(*,*)' Input 0 = no debugging output; 1 = output each iteration'
	read(*,*)tanOut(8)
	numPTpoints = 1
!	----------------- LOOP starts here -----------------
	newP = 0
	iTP = 1
	tc4 = tc
	pb4 = pb
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)

23250	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
		write(*,*)'This assemblage doesnot subtend the bulk composition - pick another assemblage'
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 100
		endif
	if(tanOut(8).eq.1)then
		write(12,*)' ---------------------------'
		write(12,*)' Modifying components ',coname(iFe),coname(iMn),coname(iCa),coname(iNa)
		write(12,*)'                   BULK                                                                      Garnet'
		write(12,23260)(coname(i),i=1,nc)
23260		format(3x,40a10)
		endif
!===========================================================
!	mincFD = 0.01d0		! increment of moles for finite difference calculations
	mincFD = 0.0001d0		! increment of moles for finite difference calculations
	write(*,*)'Input value for mincFSD  (around 0.001)'
	read(*,*)mincFD
	iCount = 0		! counter for iterations
23200	continue		! loop back to here from below finite difference convergence
	icount = icount + 1
	if(icount.gt.1000)then
		write(*,*)'iCount iterations > 1000 .... something must be wrong'
		pause 'Hit return to continue'
		go to 23500
		endif

!	Get the initial compositions of phases using the current bulk composition
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#1 Calc base: Problem in ComputePTXM trying to calculate the base assemblage'
		write(12,*)TC4,PB4
		write(12,*)(asmCurrent(k),k=1,numCurrent)
		write(12,23261)(molesStart(i),i=1,nc)
		write(12,*)'Iteration = ',icount
		pause 'Hit return to continue'
		go to 23500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalc = XphCo(kGrt,2)
	XspsCalc = XphCo(kGrt,3)
	XgrsCalc = XphCo(kGrt,4)
	XAnCalc =  XphCo(kPlg,2)

	if(tanOut(8).eq.1)then		
		Write(12,*)' iCount = ',iCount
		write(12,23261)(molesStart(i),i=1,nc)
23261		format(20F10.5)
		write(12,*)'     Xprp      Xalm      Xsps      Xgrs      XAn'
		write(12,23261)(1.-XalmCalc - XspsCalc - XgrsCalc),XalmCalc,XspsCalc,XgrsCalc,XAnCalc
		endif

!	write(12,21261)(molesStart(i),i=1,nc),XalmCalc,XspsCalc,XgrsCalc

	tol = 0.0001d0
!	Check for convergence
	if(dabs(XalmCore - XalmCalc).gt.tol)go to 23201
	if(dabs(XspsCore - XspsCalc).gt.tol)go to 23201
	if(dabs(XgrsCore - XgrsCalc).gt.tol)go to 23201
	if(dabs(XanCore  - XanCalc).gt.tol)go to 23201

!	if we get to here, then we have converged

!	Calculate the Wt% composition from the moles
!	sum = 0.0d0 	! we don't normalize because that changes rock volume we are modeling
	write(12,*)'         '
	write(12,*)' We have converged'
	write(12,*)' Modifying components ',coname(iFe),coname(iMn),coname(iCa),coname(iNa)
	call PrinTT(0)
	write(12,23260)(coname(i),i=1,nc)
	write(12,23261)(molesStart(i),i=1,nc)	
	do 23208 i = 1,nc
	wtpct(i) = molesStart(i)*molwt(i)/NumCatInOxide(i)
23208	continue
	write(12,23261)(WtPct(i),i=1,nc)
	write(12,*)'     Xprp      Xalm      Xsps      Xgrs      Xan'
	write(12,23261)(1.-XalmCalc - XspsCalc - XgrsCalc),XalmCalc,XspsCalc,XgrsCalc,XanCalc
	write(12,*)' '
	write(*,*)' We have converged -- solution is in the output window'
	pause 'Hit return to continue'
	go to 23500

23201	continue
!	calculate finite difference derivatives
!	Derivative for Fe in bulk composition
	molesStart(iFe) = molesStart(iFe) + mincFD
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#2 Calc Fe FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(12,*)'Iteration = ',icount
		write(*,*)TC4,PB4
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 23500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcPlus = XphCo(kGrt,2)
	XspsCalcPlus = XphCo(kGrt,3)
	XgrsCalcPlus = XphCo(kGrt,4)
	XanCalcPlus = XphCo(kPlg,2)
!	Calculate derivatives
	DalmDMolesFe = (XalmCalc - XalmCalcPlus)/mincFD
	DspsDMolesFe = (XspsCalc - XspsCalcPlus)/mincFD
	DgrsDMolesFe = (XgrsCalc - XgrsCalcPlus)/mincFD
	DanDMolesFe = (XanCalc - XanCalcPlus)/mincFD
	molesStart(iFe) = molesStart(iFe) - mincFD		! set back to starting composition
	
!	Derivative for Mn in bulk composition
	molesStart(iMn) = molesStart(iMn) + mincFD
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#3 Calc Mn FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(12,*)'Iteration = ',icount
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 23500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcPlus = XphCo(kGrt,2)
	XspsCalcPlus = XphCo(kGrt,3)
	XgrsCalcPlus = XphCo(kGrt,4)
	XanCalcPlus = XphCo(kPlg,2)
!	Calculate derivatives
	DalmDMolesMn = (XalmCalc - XalmCalcPlus)/mincFD
	DspsDMolesMn = (XspsCalc - XspsCalcPlus)/mincFD
	DgrsDMolesMn = (XgrsCalc - XgrsCalcPlus)/mincFD
	DanDMolesMn = (XanCalc - XanCalcPlus)/mincFD
	molesStart(iMn) = molesStart(iMn) - mincFD		! set back to starting composition

!	Derivative for Ca in bulk composition
	molesStart(iCa) = molesStart(iCa) + mincFD
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#4 Calc Ca FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(12,*)'Iteration = ',icount
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 23500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcPlus = XphCo(kGrt,2)
	XspsCalcPlus = XphCo(kGrt,3)
	XgrsCalcPlus = XphCo(kGrt,4)
	XanCalcPlus = XphCo(kPlg,2)
!	Calculate derivatives
	DalmDMolesCa = (XalmCalc - XalmCalcPlus)/mincFD
	DspsDMolesCa = (XspsCalc - XspsCalcPlus)/mincFD
	DgrsDMolesCa = (XgrsCalc - XgrsCalcPlus)/mincFD
	DanDMolesCa = (XanCalc - XanCalcPlus)/mincFD
	molesStart(iCa) = molesStart(iCa) - mincFD		! set back to starting composition

!	Derivative for Na in bulk composition
	molesStart(iNa) = molesStart(iNa) + mincFD
	izero = 0
	isubtend = 0		! subtend isn't used in ComputPTXM (should be removed from argument)
	call ComputePTXM(TC4,PB4,izero)
	if(izero.gt.0)then
		write(*,*)'#5 Calc Ca FD: Problem in ComputePTXM trying to calculate the base assemblage'
		write(*,*)TC4,PB4
		write(12,*)'Iteration = ',icount
		write(*,*)(asmCurrent(k),k=1,numCurrent)
		pause 'Hit return to continue'
		go to 23500
		endif
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
	XalmCalcPlus = XphCo(kGrt,2)
	XspsCalcPlus = XphCo(kGrt,3)
	XgrsCalcPlus = XphCo(kGrt,4)
	XanCalcPlus = XphCo(kPlg,2)
!	Calculate derivatives
	DalmDMolesNa = (XalmCalc - XalmCalcPlus)/mincFD
	DspsDMolesNa = (XspsCalc - XspsCalcPlus)/mincFD
	DgrsDMolesNa = (XgrsCalc - XgrsCalcPlus)/mincFD
	DanDMolesNa =  (XanCalc - XanCalcPlus)/mincFD
	molesStart(iNa) = molesStart(iNa) - mincFD		! set back to starting composition


!	Set up to solve the equations
!	-(XalmCore - XalmCalc) = 0 = (XalmCore - XalmCalc) + (dXalm/dFe)*ÆFe + (dXalm/dMn)*ÆMn + (dXalm/dCa)*ÆCa
!	-(XspsCore - XspsCalc) = 0 = (XspsCore - XspsCalc) + (dXsps/dFe)*ÆFe + (dXsps/dMn)*ÆMn + (dXsps/dCa)*ÆCa
!	-(XgrsCore - XgrsCalc) = 0 = (XgrsCore - XgrsCalc) + (dXgrs/dFe)*ÆFe + (dXgrs/dMn)*ÆMn + (dXgrs/dCa)*ÆCa

	A(1,5) = -(XalmCore - XalmCalc)
	A(2,5) = -(XspsCore - XspsCalc)
	A(3,5) = -(XgrsCore - XgrsCalc)
	A(4,5) = -(XanCore  - XanCalc)

	A(1,1) = DalmDMolesFe
	A(1,2) = DalmDMolesMn
	A(1,3) = DalmDMolesCa
	A(1,4) = DalmDMolesNa
	A(2,1) = DspsDMolesFe
	A(2,2) = DspsDMolesMn
	A(2,3) = DspsDMolesCa
	A(2,4) = DspsDMolesNa
	A(3,1) = DgrsDMolesFe
	A(3,2) = DgrsDMolesMn
	A(3,3) = DgrsDMolesCa
	A(3,4) = DgrsDMolesNa
	A(4,1) = DanDMolesFe
	A(4,2) = DanDMolesMn
	A(4,3) = DanDMolesCa
	A(4,4) = DanDMolesNa

	numEq = 4
	if(tanOut(8).eq.1)then		
		write(12,*)' '
		write(12,*)'A matrix'
		do 23040 i = 1,numEq
		write(12,21080)(A(i,j),j=1,numEq+1)
23080		format(50F15.5)
23040		continue
		endif		

!     Find solution
	nsolve = 5
	DETA=0.D0
	IER=0
	CALL REDUCE (numEq,NSOLVE,DETA,IER)
	IF (IER.EQ.1) then
	      WRITE(12,*)' ************ ERROR **************************'
	      WRITE(12,*)' Matrix failed to invert in SUBROUTINE REDUCE'
	      write(12,*)' We are in Subroutine ParallalToTangent in file TangentRoutines.f '
	      izero=1
	      write(12,*) 'Hit return to continue...'
	      pause 'Hit return to continue...'
	      go to 23500
	      endif

	if(tanOut(8).eq.1)then		
		write(12,*)'Solution xx'
		write(12,23080)(xx(j,1),j=1,numEq)
		write(12,*)'Determinant = ',DETA
		endif		


!	Calculate new bulk composition
!	First check to see that no moles go negative
!	Damp = 1.
	Damp = .1D0		!  maximum damp
23041	continue
	if(tanOut(8).eq.1)then		
		write(12,*)' Damp = ',damp
		endif
	if((damp*xx(1,1)+molesStart(iFe)).le.0.0d0.or.	&
	   (damp*xx(2,1)+molesStart(iMn)).le.0.0d0.or.	&
	   (damp*xx(3,1)+molesStart(iCa)).le.0.0d0.or.	&
	   (damp*xx(4,1)+molesStart(iNa)).le.0.0d0) then
	   	Damp = Damp/2.0d0
		if(damp.lt.1.0D-8)then
			write(*,*)' Damp = ',damp
			pause 'something must be wrong ... damp is too small'
			go to 23500
			endif
	   	go to 23041
	   	endif
	molesStart(iFe) = molesStart(iFe) + xx(1,1)*damp
	molesStart(iMn) = molesStart(iMn) + xx(2,1)*damp
	molesStart(iCa) = molesStart(iCa) + xx(3,1)*damp
	molesStart(iNa) = molesStart(iNa) + xx(4,1)*damp
	if(tanOut(8).eq.1)then		
		write(12,23261)(molesStart(i),i=1,nc)
		endif
	
	go to 23200

	endif
!===================================


! ------Auto calculate from "Find initial tangent" from assemblage-------------
	if(iTan.eq.25)then
!	tanPhase = 0
	iCareifMIS = 1		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 1		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 1			! supresses error output in subroutine compute
	tanOut(1) = 1
	tanOut(3) = 1
	tanOut(7) = 1
!	call SetTanOut()
2500	continue
	write(*,*)' Debug options'
	write(*,*)' 0 = return'
	write(*,*)' 1 = Input P and T'
	write(*,*)' 2 = input new assemblage'
	write(*,*)' 3 = Solve using COMPUTPTXM'
	write(*,*)' 4 = Call Parallel to Tangent'
	write(*,*)' 5 = Calculate new tangent with phase compositions'
	write(*,*)' 6 = set initial tangent'
	write(*,*)' 8 = set output options'
	write(*,*)' 7 = PRINTT'
	write(*,*)' 9 = Save as a MIF file'
	read(*,*)debug
	if(debug.eq.0)then
		go to 100
		endif
!------------------------------------------------------
	if(debug.eq.8)then
		call SetTanOut()
		go to 2500
		endif
!------------------------------------------------------
	if(debug.eq.7)then
		call PRINTT(1)
		go to 2500
		endif
!------------------------------------------------------
	if(debug.eq.1)then

		write(*,*)'Input new T and P'
		read(*,*)TC,PB
		if(TC.eq.0.)go to 100
		tc4 = tc
		pb4 = pb
		ALLX(1,1)=TC
		ALLX(1,2)=PB
		ALLX(2,1)=TC
		ALLX(2,2)=PB
		write(12,*)' ===================================================================================='
		write(12,*)' ===================================================================================='
		write(12,*)tc4,pb4
		call CalcMuZero(izero)			! calculates for all phases at T&P - only needs to be done once at T&P
		go to 2500
		endif
!------------------------------------------------------
	if(debug.eq.2)then
!		xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
		write(*,*)'Specify initial assemblage'
		write(*,*)AssembTitle
		do 2510 i = 1,NumPhMIF
2510		write(*,*)i,phName(i)
		numCurrent = 0
2511		continue
		read(*,*)i
		if(i.ne.0)then
			numCurrent = numCurrent + 1
			asmCurrent(numCurrent) = i
			go to 2511
			else
			call listCurrentAsm
		!	call PRINTT(1)
			endif 
!		xPhCo = xPhCoInitial			! see if starting with initial comps works
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)
		go to 2500
		endif
		
!------------------------------------------------------
!	calculate initial tangent using COMPUT4MAD
	if(debug.eq.3)then
!      		SET DEFAULTS FOR MONITOR AND DELTAX
		SNstep = 1
		!ndep = 1
		NSTEP=1
		SNSTEP=1
		do 2501 i=1,nvar-neq   ! should always be 2
		smon(i) = i
		mon(i) = i
		sdel(i)=0.d0
		deltax(i)= 0.00
2501		continue
!   		Set up array IPOINT to contain pointers to non-monitor parameters
		J=0
		Do 2512 i=1,NVAR
		DO 2513 L=1,NVAR-NEQ
		IF(MON(L).eq.i)go to 2512
2513		continue
		J=J+1
		IPOINT(J) = I
2512		CONTINUE
!      		IF(iLong(4).eq.1) write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)
		write(12,*)' Monitors are:    ',(Mon(I),I=1,nvar-NEQ)
		write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)

!      		THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
		DO 2503 I=1,nvar
!       	ALLX(3,I) IS WHERE calculations STARTED-save
	        ALLX(3,I)=ALLX(1,I)
2503		CONTINUE

!		do 599 istep=1,snstep
!       	  Store new values for each independent variable (monitor parameter)
! 		This is where Neton's method increments composition, T or P
		DO 2520 J=1,NVAR-NEQ
		i=mon(j)
		ALLX(1,i)=ALLX(1,i) + DELTAX(J)
2520		continue
        	call SetTPX
         
		write(*,2561)TC,PB
2561   		format (' T = ',f10.1,'  P = ',f10.1)
!      		write(*,*)' NSTEP = ',snstep,' Counter = ',istep

        	IZERO=0
        	call c2ompute(izero)
        	IF(IZERO.gt.0)then
	        	write(*,*)' izero error in call to c2ompute. izero =  ',izero
	        	pause 'hit return to continue'
	        	GO TO 2500
			endif
		endif


!------------------------------------------------------
	if(debug.eq.4)then

!	Calculate parallel to tangent for current assemblage
!	Note that this is the assemblage we put in option=2
!	It can be after either 
!		(a) COMPUT4MAD or 
!		(b) After just resetting the phases that define the tangent
!			This is done with routine AdjustTangentToAsm

		write(12,*)' -------------------'
		write(12,*)'Tangent before call to AdjustTangentToAsm '
		if(tanOut(7).eq.1)call ListTangentPlane
		write(12,*)' -------------------'
		call AdjustTangentToAsm
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
		write(12,*)' -------------------'
		write(12,*)'Tangent after call to AdjustTangentToAsm '
		if(tanOut(7).eq.1)call ListTangentPlane
		write(12,*)' -------------------'
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
		if(tanOut(3).eq.1)call ListGDifference
!		write(*,*)' quartz before PRINTT call ',gattp(1,1)
	         CALL PRINTT(1)
!		write(*,*)' quartz after PRINTT call ',gattp(1,1)
		call WriteAllStuffToTerminal(12)		! write results to output window

		pause 'Hit return to continue'	

		go to 2500
		endif
!------------------------------------------------------
	if(debug.eq.6)then
!		Explore the tangent
!		Set the tangent plane
! Current tangent Plane   700.000000000000  20000.0000000000
!         Si             Al             Ti             Mg             Fe             Mn             Ca             Na             K              H2O  
!    -932679.695    -865184.536    -979100.248    -665108.099    -354753.177    -514965.223    -798943.652    -399801.571    -433680.158    -348012.637

!         Si             Al             Ti             Mg             Fe             Mn             Ca             Na             K              H2O  
!    
		tanPlane(1) = -932679.695   	!Si
		tanPlane(2) = -865184.536	!Al
		tanPlane(3) = -979100.248  	!Ti
		tanPlane(4) = -665108.099   	!Mg
		tanPlane(5) = -354753.177  	!Fe
		tanPlane(6) = -514965.223  	!Mn
		tanPlane(7) = -798943.652  	!Ca
		tanPlane(8) = -399801.571  	!Na
		tanPlane(9) = -433680.158   	!K
		tanPlane(10)= -348012.637	!H2O


!		tanPlane(1) = 	
!		write(12,81)(tanPlane(i),i=1,nc) !tangent plane is from system components
	
		write(12,*)' -------------------'
		write(12,*)'Tangent before call to AdjustTangentToAsm '
		if(tanOut(7).eq.1)call ListTangentPlane
		write(12,*)' -------------------'
!		call AdjustTangentToAsm
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
!		write(12,*)' -------------------'
!		write(12,*)'Tangent after call to AdjustTangentToAsm '
!		if(tanOut(7).eq.1)call ListTangentPlane
		write(12,*)' -------------------'
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!		if(tanOut(3).eq.1)call ListGDifference
!		write(*,*)' quartz before PRINTT call ',gattp(1,1)
!	         CALL PRINTT(1)
!		write(*,*)' quartz after PRINTT call ',gattp(1,1)
		call WriteAllStuffToTerminal(12)		! write results to output window

		pause 'Hit return to continue'	

		write(*,*)' Lowering tangent plane - input value to lower'
		read(*,*)maxGdifference
		do 25601 i = 1,nc
		tanPlane(i) = tanPlane(i) - maxGdifference
25601		continue
		call listTangentPlane
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!		if(tanOut(3).eq.1)call ListGDifference
!		write(*,*)' quartz before PRINTT call ',gattp(1,1)
!	        CALL PRINTT(1)
!		write(*,*)' quartz after PRINTT call ',gattp(1,1)
		call WriteAllStuffToTerminal(12)		! write results to output window
		pause ' hit return to continue'
		go to 2500
		endif
!------------------------------------------------------
	if(debug.eq.5)then

		call PhaseComp()
!		this routine returns the composition of phases in terms of system components
!		in the arran PhComp(K,i) for the current assemblage
!		It also returns Wt% oxides and Wt% elements, but I dont' need those.

! 		Calculate G of the phase at this composition
		gOfPhaseAU(K) = 0.0D0
		do 25210 j = 1,numPhCo(K)
		gOfPhaseAU(K) = gOfPhaseAU(K) + xPhCo(k,j)*(uzero(k,j) + Rjoules*TK*lnAct(k,j))/atomNorm(k,j)
25210		continue

		endif
!------------------------------------------------------
	if(debug.eq.9)then
		call SaveMIF(TC4,PB4)
		go to 2500
		endif
	endif
!------------------------------------------------------
!	Jacobi's ladder routine
	if(iTan.eq.26)then
!	tanPhase = 0
	iCareifMIS = 1		! i do not care if matrix is singular (do not print error message)
	iLong(1) = 1		! turn off error reporting
	idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
	ioutputpseudo = 1			! supresses error output in subroutine compute
	tanOut(1) = 1
	tanOut(3) = 1
	tanOut(7) = 1
!	call SetTanOut()
2600	continue
	write(*,*)' Debug options'
	write(*,*)' 0 = return'
	write(*,*)' 1 = Input P and T'
	write(*,*)' 2 = input new assemblage'
	write(*,*)' 3 = Set initial tangent'
	write(*,*)' 4 = Calculate Jacobian'
	write(*,*)' 8 = set output options'
	write(*,*)' 7 = PRINTT'
	write(*,*)' 9 = Save as a MIF file'
	write(*,*)' 10 = Solve using COMPUTPTXM'
	write(*,*)' 11 = Call Parallel to Tangent'
	read(*,*)debug
	if(debug.eq.0)then
		go to 100
		endif
!------------------------------------------------------
	if(debug.eq.1)then

		write(*,*)'Input new T and P'
		read(*,*)TC,PB
		if(TC.eq.0.)go to 100
		tc4 = tc
		pb4 = pb
		ALLX(1,1)=TC
		ALLX(1,2)=PB
		ALLX(2,1)=TC
		ALLX(2,2)=PB
		write(12,*)' ===================================================================================='
		write(12,*)' ===================================================================================='
		write(12,*)tc4,pb4
		call CalcMuZero(izero)			! calculates for all phases at T&P - only needs to be done once at T&P
		go to 2600
		endif
!------------------------------------------------------
	if(debug.eq.2)then
		TC = 700
		PB = 9000
		write(12,*)' ===================================================================================='
		write(12,*)' ===================================================================================='
		write(12,*)tc4,pb4
		call CalcMuZero(izero)			! calculates for all phases at T&P - only needs to be done once at T&P
!		xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
		write(*,*)'Specify initial assemblage'
		write(*,*)AssembTitle
		do 2610 i = 1,NumPhMIF
2610		write(*,*)i,phName(i)
		numCurrent = 0
2611		continue
		read(*,*)i
		if(i.ne.0)then
			numCurrent = numCurrent + 1
			asmCurrent(numCurrent) = i
			go to 2611
			else
			call listCurrentAsm
		!	call PRINTT(1)
			endif 
!		xPhCo = xPhCoInitial			! see if starting with initial comps works
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)
		go to 2600
		endif
		


!------------------------------------------------------
	if(debug.eq.3)then
!		Explore the tangent
!		Set the tangent plane
! Current tangent Plane   700.000000000000  20000.0000000000
!         Si             Al             Ti             Mg             Fe             Mn             Ca             Na             K              H2O  
!    -932679.695    -865184.536    -979100.248    -665108.099    -354753.177    -514965.223    -798943.652    -399801.571    -433680.158    -348012.637

!         Si             Al             Ti             Mg             Fe             Mn             Ca             Na             K              H2O  
!    
		tanPlane(1) = -932679.695   	!Si
		tanPlane(2) = -865184.536	!Al
		tanPlane(3) = -979100.248  	!Ti
		tanPlane(4) = -665108.099   	!Mg
		tanPlane(5) = -354753.177  	!Fe
		tanPlane(6) = -514965.223  	!Mn
		tanPlane(7) = -798943.652  	!Ca
		tanPlane(8) = -399801.571  	!Na
		tanPlane(9) = -433680.158   	!K
		tanPlane(10)= -348012.637	!H2O

		tanPlane(1) = -932679.695   	!Si
		tanPlane(2) = -678276.874   	!Mg
		tanPlane(3) = -351305.200  	!Fe
 !   -957298.237    -689726.674    -350572.533
!    -957298.237    -678276.874    -351305.200
		call CalcMuZero(izero)			! calculates for all phases at T&P - only needs to be done once at T&P


!		tanPlane(1) = 	
!		write(12,81)(tanPlane(i),i=1,nc) !tangent plane is from system components
	
		write(12,*)' -------------------'
		write(12,*)'Initial Tangent plane'
		if(tanOut(7).eq.1)call ListTangentPlane
		write(12,*)' -------------------'
!		call AdjustTangentToAsm
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
!		write(12,*)' -------------------'
!		write(12,*)'Tangent after call to AdjustTangentToAsm '
!		if(tanOut(7).eq.1)call ListTangentPlane
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!		if(tanOut(3).eq.1)call ListGDifference
!		write(*,*)' quartz before PRINTT call ',gattp(1,1)
!	         CALL PRINTT(1)
!		write(*,*)' quartz after PRINTT call ',gattp(1,1)
		call WriteAllStuffToTerminal(12)		! write results to output window
!		pause 'Hit return to continue'	
!		Find phase that sits farthest below tangent
		maxgDifference = -1.0d-2
!		If a phase is more than 0.01 J below the tangent, then flag it.
		maxgDifferencePhase = 0
		do k = 1,numPhMIF
			if(k.eq.iExclude)go to 26440			! skip the excluded phase even if it is below the tangent
			if(gDifferenceAU(k).lt.maxgDifference)then			
				maxgDifference = gDifferenceAU(k)
				maxgDifferencePhase = k
				endif
26440			continue
			end do


		write(*,*)' Lowering tangent plane by',maxGdifference
		write(12,*)' Maximum Gdifference phase is ',maxGDifferencePhase
		write(12,*)' Lowering tangent plane by    ',maxGdifference
		write(*,*)
!		read(*,*)maxGdifference
		do i = 1,nc
			tanPlane(i) = tanPlane(i) + maxGdifference
			end do
		numCurrent = 1
		asmCurrent(numCurrent) = maxGDifferencePhase
		write(12,*)' -------------------'
		write(12,*)'Adjusted Tangent plane'
		call listTangentPlane
		write(12,*)' -------------------'
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!		if(tanOut(3).eq.1)call ListGDifference
!		write(*,*)' quartz before PRINTT call ',gattp(1,1)
!	        CALL PRINTT(1)
!		write(*,*)' quartz after PRINTT call ',gattp(1,1)
		call WriteAllStuffToTerminal(12)		! write results to output window
		pause ' hit return to continue'
		go to 2600
		endif
!------------------------------------------------------
	if(debug.eq.4)then

		call ComputeJacobian(izero)
		write(*,*)'izero = ',izero
		pause 'Hit return to continue'
		go to 2600
		endif
!------------------------------------------------------
	if(debug.eq.8)then
		call SetTanOut()
		go to 2600
		endif
!------------------------------------------------------
	if(debug.eq.7)then
		call PRINTT(1)
		go to 2600
		endif
!------------------------------------------------------
	if(debug.eq.9)then
		call SaveMIF(TC4,PB4)
		go to 2600
		endif
!------------------------------------------------------
!------------------------------------------------------
!	calculate initial tangent using COMPUT4MAD
	if(debug.eq.10)then
!      		SET DEFAULTS FOR MONITOR AND DELTAX
		SNstep = 1
		!ndep = 1
		NSTEP=1
		SNSTEP=1
		do 2601 i=1,nvar-neq   ! should always be 2
		smon(i) = i
		mon(i) = i
		sdel(i)=0.d0
		deltax(i)= 0.00
2601		continue
!   		Set up array IPOINT to contain pointers to non-monitor parameters
		J=0
		Do 2612 i=1,NVAR
		DO 2613 L=1,NVAR-NEQ
		IF(MON(L).eq.i)go to 2612
2613		continue
		J=J+1
		IPOINT(J) = I
2612		CONTINUE
!      		IF(iLong(4).eq.1) write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)
		write(12,*)' Monitors are:    ',(Mon(I),I=1,nvar-NEQ)
		write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)

!      		THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
		DO 2603 I=1,nvar
!       	ALLX(3,I) IS WHERE calculations STARTED-save
	        ALLX(3,I)=ALLX(1,I)
2603		CONTINUE

!		do 599 istep=1,snstep
!       	  Store new values for each independent variable (monitor parameter)
! 		This is where Neton's method increments composition, T or P
		DO 2620 J=1,NVAR-NEQ
		i=mon(j)
		ALLX(1,i)=ALLX(1,i) + DELTAX(J)
2620		continue
        	call SetTPX
         
		write(*,2661)TC,PB
2661   		format (' T = ',f10.1,'  P = ',f10.1)
!      		write(*,*)' NSTEP = ',snstep,' Counter = ',istep

        	IZERO=0
        	call c2ompute(izero)
        	IF(IZERO.gt.0)then
	        	write(*,*)' izero error in call to c2ompute. izero =  ',izero
	        	pause 'hit return to continue'
	        	GO TO 2600
			endif
!599   		continue
		endif

!------------------------------------------------------
	if(debug.eq.11)then

!	Calculate parallel to tangent for current assemblage
!	Note that this is the assemblage we put in option=2
!	It can be after either 
!		(a) COMPUT4MAD or 
!		(b) After just resetting the phases that define the tangent
!			This is done with routine AdjustTangentToAsm

		write(12,*)' -------------------'
		write(12,*)'Tangent before call to AdjustTangentToAsm '
		if(tanOut(7).eq.1)call ListTangentPlane
		write(12,*)' -------------------'
		call AdjustTangentToAsm
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
		write(12,*)' -------------------'
		write(12,*)'Tangent after call to AdjustTangentToAsm '
		if(tanOut(7).eq.1)call ListTangentPlane
		write(12,*)' -------------------'
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
		if(tanOut(3).eq.1)call ListGDifference
!		write(*,*)' quartz before PRINTT call ',gattp(1,1)
	         CALL PRINTT(1)
!		write(*,*)' quartz after PRINTT call ',gattp(1,1)
		call WriteAllStuffToTerminal(12)		! write results to output window

		!pause 'Hit return to continue'	

		go to 2600
		endif
		
!--------------------------------------------------------------------

	go to 2600

	endif


! -------  Tangent AFM routine --------------
	if(iTan.eq.50)then

		TC = Tgrid(1)
		PB = Pgrid(1)
		tc4 = tc
		pb4 = pb
	!	tanPhase = 0
		iCareifMIS = 0		! i do not care if matrix is singular (do not print error message)
		iLong(1) = 0		! turn off error reporting
		idoingpseudo = 1			! lets subroutine compute know we are calling from a pseudosection routine
		ioutputpseudo = 0			! supresses error output in subroutine compute
		iDoingGrid = 1				! keep going through activity calculations if we are doing a grid
		call ZeroTanOut()
		call FSS_Alert('Alert','Open output file')
		call OpenOutputFile(iok)
		if(iok.ne.0)go to 100			! abort
		logFile = trim(tangentOutputFile)//'.log'
	!	logFile = 'TempLogFile.log'
		open(95,FILE=logFile,status = 'UNKNOWN')
		write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
	! 	do 4212 i = 1,numPhMIF
	! 4212	write(*,*)i,phName(i)
	! 	write(*,*)'Is there any phase to exclude from equilibrim assemblages (for calculation of affinities)?'
	! 	write(*,*)'Input a number from the list or 0 for none (i.e. include all phases)'
	! 	read(*,*)iExclude
	! 	if(Tgrid(1).lt.Tgrid(numPTpoints))then
	! 		call ChooseConstantPorosity()
	! 		endif
	!	stuffFile = trim(tangentOutputFile)//'.stuff'
	!	open(72,FILE=stuffFile,status = 'UNKNOWN')
		AllFile = trim(tangentOutputFile)//'.All'
		open(74,FILE=AllFile,status = 'UNKNOWN')
		call WriteAllOutHeader(74,Tinc,Pinc)

	! 	call StableASM(logfileon,TC4,PB4,ierr)
	! 	if(ierr.eq.0)then
	! 		write(*,*)'Initial stable assemblage found:',TC4,PB4
	! 		call WriteAllStuff()
	! !		pause ' hit return to continue'
	! 		else
	! 		write(*,*)' No solution for initial assemblage found'
	! 		pause ' hit return to continue'
	! 		go to 100
	! 		endif

		Call AFMTangent(logFileOn,iGrid,ierr)
		if(ierr.eq.1)then
			write(*,*)' Some error occurred.... hit return to continue'
			pause
			go to 100
			endif

		go to 100
		endif




!--------------------------------------------------------------------

	go to 100
999	continue
	return
	end
	
