! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine StableASM(logFileOn,TC4,PB4,ierr)
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

	integer*4 i,j,iok,k,izero,iall,isubtend,maxNegativeMolesPhase,maxgDifferencePhase
	integer*4 asmCurrentStore(phMax),kCur,countTheMisses,ierr,WriteAll
	real*4 Tc4,Pb4,Pb4old
	real*8 maxgDifference,maxNegativeMoles,gDiffThreshold
	real*8 sum
	integer*4 logFileOn,iremove
	integer*4 EOFflag
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

	open(87,file='Fort_87_ComboScratch',status = 'UNKNOWN')
	iParagonite = 0
	write(95,*)'   '
	write(95,*)' Phases in this MIF'
	do 228 i = 1,numPhMIF
	if(minRec(i).eq.118)then
		iParagonite = i			! this is a switch so we can initialize Xparag to a large value
		endif
	write(95,*)i,PhName(i)
228	write(*,*)i,phName(i)
	write(95,*)'   '
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
	if(numCurrent.gt.0)then		! only write out if an assemblage was specified
		write(87,*)numCurrent,(asmCurrent(i),i=1,numCurrent)
		endif
	write(*,*)'iParagonite = ',iParagonite
	iall = 0
!	now sort through and reorder possible assemblages with those with likely phases first
	rewind(85)
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

!	xPhCoLast = xPhCoInitial		! Store the Initial set of phase compositions
!	xPhCoNewP = xPhCoInitial		! Store the last set of good phase compositions for resetting if needed
!	MPNewP = MP0
!	MPlast = MP0
	call SetLast()
	call SetNewP()
	write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
402	format(2F12.2,20F12.5)
	write(95,*)' Starting TP calculations'

	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	if(izero.gt.0)then
		write(*,*)' Error calculating MU zero for this T&P... abort the run'
		write(95,*)' Error calculating MU zero for this T&P... abort the run'
		pause 'take a look and hit return when ready'
		ierr = 1	! error code
		return
		endif

	Call RewindTheFile(87,logFileOn)
	EOFflag = 0
	countTheMisses = 0



410	continue
!	get the first assemblage to check
!	call ResetToLast()	! try this 1/11/2020
!	Reset to the compositions in the input file
	call ResetToStart()	! start with initial compositions when examining a new assemblage
	call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)		! get the initial assemblage to test
		if(EOFflag.eq.1)then
			write(*,*)'      Every assemblage has been examined and no solution found'
			write(*,*)'      Aborting the run.....'
			write(95,*)'      Every assemblage has been examined and no solution found'
			write(95,*)'      Aborting the run.....'
			ierr = 1	! error code
			return
			endif		
420	continue
	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
!	call dumpeverything('Before Compute')
	if(isubtend.eq.0)then
		write(95,*)'      Assemblage does not subtend BC (line 168 in Sub StableAsm) - get new asm'
		go to 410
		endif
	izero = 0
	call ComputePTXM(TC4,PB4,izero)
!	call dumpeverything('After Compute')
	if(izero.gt.0)then
		write(95,*)'      Error in ComputePTXM (Line 175 of Sub StableASM) - get new asm'
!		xPhCo = xPhCoNewP		! reset to the initial compositions
!		MP0 = MPnewP
		call ResetToLast()
!		xPhCo = xPhCoLast		! reset to the the last compositions that worked
!		MP0 = MPLast
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
	if(mp0(k)/atomNorm(k,1).lt.maxNegativeMoles)then
!		A phase has negative moles - record and loop back
		maxNegativeMoles = mp0(k)*atomNorm(k,1)		! we are storing moles/atom (not per molecule)
		maxNegativeMolesPhase = kCur
		endif
430	continue
	if(maxNegativeMolesPhase.ne.0)then
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
			call ResetToLast()
!			xPhCo = xPhCoNewP		! reset to the initial compositions
!			MP0 = MPnewP
			write(95,*)'      Missed>20 (Negative moles)... loop back'
			go to 410
			endif
		iremove = asmCurrent(maxNegativeMolesPhase)
		mp0(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		do 432 j = maxNegativeMolesPhase,numPh
		asmCurrent(j) = asmCurrent(j+1)
432		continue
		numCurrent = numCurrent - 1
		write(95,*)'Removed phase= ',iremove,'asm ',(asmCurrent(iok),iok=1,numCurrent)
		countTheMisses = countTheMisses + 1

!		call ResetToLast()	!try this 1/11/2020 (Why reset if we have just removed a phase?)
		go to 420
		endif

!	if here, then no phase has negative moles
!	Check to see if any phases are below the tangent plane
!	write(*,*)'               First call to ATTA'
	call AdjustTangentToAsm
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
	if(WriteAll.eq.0)then
		call WriteAllStuffToTerminal(6)
		write(*,*)' 0 = next with output; 1 = continue w/o output; 2 = abort'
		read(*,*)WriteAll
		if(WriteAll.eq.2)then
			ierr = 1	! error code
			return
			endif
!		pause ' hit return to continue'
		endif
!	xPhCoLast = xPhCo		! store all compositions that are parallel to the tangent
435	continue
	maxgDifference = -1.0d20	! Set very negative to start
!	maxgDifference = -1.0d-2	
!	If a phase is more than 0.01 J below the tangent, then flag it.
! &&&&&&&&&&&&&&&
!	I need to make this "minimum negative G"
	gDiffThreshold = -1.0e-2
	maxgDifferencePhase = 0
	do 440 k = 1,numPhMIF
	if(k.eq.iExclude)go to 440			! skip the excluded phase even if it is below the tangent
	if(gDifferenceAU(k).lt.gDiffThreshold)then	! if true then the phase lies below the tangent
		if(gDifferenceAU(k).gt.maxgDifference)then			
			maxgDifference = gDifferenceAU(k)	! this should find the phase with the smallest GDiff (closest to tangent)
			maxgDifferencePhase = k
			endif

		endif
!	if(gDifferenceAU(k).lt.maxgDifference)then			
!		maxgDifference = gDifferenceAU(k)
!		maxgDifferencePhase = k
!		endif
440	continue
	if(maxgDifferencePhase.ne.0)then		
		if(countTheMisses.gt.20)then		! we can't find a solution by iterating, so just keep reading assemblages from the file
			call ResetToLast()	!try this 1/11/2020
!			if(iParagonite.gt.0)then
			!	xPhCo(iParagonite,1) = 0.05
			!	xPhCo(iParagonite,2) = 0.95
!				endif		
			write(95,*)'      Missed>20 (maxGdifference)... loop back'
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
!		write(95,*)'Phase added  ',(minRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		write(95,*)'Phase added    ',(asmCurrent(iok),iok=1,numCurrent),maxGDifference
		countTheMisses = countTheMisses + 1
		isubtend = 1
	!	call ResetToLast()	!try this 1/11/2020  -- WHY reset since we just added a phase
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)		! this one MUST subtend because the last one did and we added a phase
		izero = 0
		call ComputePTXM(TC4,PB4,izero)
		if(izero.eq.0)then
			go to 425			! this assemblage worked - go back and check for negative moles
			endif
	!call dumpeverything('Before Subassemblages')

!		If we get here, then the last assemblage didn't work (izero not = 0)
		write(95,*)'*****************************************************************'
		write(95,*)' Start subassemblage routine'
		write(95,*)'*****************************************************************'
		write(95,*)'	Calculation failed in ComputePTXM (line 278 in Sub StableAsm) - Work on subassemblages'
		write(95,*)'	Current assemblage: ',(asmCurrent(iok),iok=1,numCurrent)

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
		call ResetToLast()	!try this 1/11/2020
!		if(iParagonite.gt.0)then
		!	xPhCo(iParagonite,1) = 0.05
		!	xPhCo(iParagonite,2) = 0.95
!			endif		
		numCurrent = icombo		          !the number of phases in the assemblage
		write(95,*)i,numCurrent,(Combos(i,k),k=1,numCurrent)
		do 455 k = 1,numCurrent
		asmCurrent(k) = asmCurrentStore(Combos(i,k))
455		continue
!		write(95,*)'	 Subassemblage  ',numCurrent,(MinRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		write(95,*)'	 Subassemblage  ',(asmCurrent(iok),iok=1,numCurrent)
		isubtend = 1
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)
		if(isubtend.eq.0)then
			write(95,*)'		Assemblage does not subtend BC (line 313 in Sub StableAsm) - get new asm'
			go to 450
			endif
		izero = 0
		call ComputePTXM(TC4,PB4,izero)
		if(izero.gt.0)then
			write(95,*)'		   Error in ComputePTXM (Line 319 of Sub StableAsm) - get new asm'
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
!		write(95,*)'               Second call to ATTA'
		call AdjustTangentToAsm
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	

!	If a phase lies below the tangent, don't I just want to skip to the next assemblage?
!	I don't think I want to add a phase at this point because we're working on subassemblages

		do 462 k = 1,numPhMIF
		if(gDifferenceAU(k).lt.-1.0d-4)then			
			write(95,*)'		   Phase lies below tangent = ',k,gDifferenceAU(k)
!	The code fails at times when a new phase (one not in the original asmCurrent) lies below the tangent
!		first check to see that this phase isn't in the original assemblage
			go to 450		! get the next assemblage
			do 463 j = 1,numCurrentStore
			if(k.eq.asmCurrentStore(j))go to 450		! this phase is already being considered. Don't reset anything
463			continue
!	It appears that this is a new phase. Try adding that phase in and going through the asmCombos again
!	Make sure this isn't the excluded phase
			if(iExclude.eq.k)go to 450
			numCurrentStore = numCurrentStore + 1
			asmCurrentStore(numCurrentStore) = k
			numCurrent = numCurrentStore
!			call SetLast()			! 3/28/2020 This assemblage "worked" but one phase was below tangent. Save 
						! -- NOT -- this assemblage didn't work. We just added a phase
			go to 452		! Start the subassemblage routine all over
			endif
462		continue
!		if we get here, then 
!			(a) we have a solution
!			(b) no phase has negative moles
!			(c) no phase sits below the tangent
		go to 475

450		continue


480		continue	! finished looping on all combos
		write(95,*)'*****************************************************************'
		write(95,*)' End of subassemblage routine'
		write(95,*)'*****************************************************************'
	
!	If we get to here then we have gone through all of the subassemblages in the ASMCompo array and none has worked
!	Loop back and try again reading assemblages from the main Compo file
!		xPhCo = xPhCoLast			! Reset to last good compositions
		call ResetToLast()			! added 3/28/2020
!		if(iParagonite.gt.0)then
		!	xPhCo(iParagonite,1) = 0.05
		!	xPhCo(iParagonite,2) = 0.95
!			endif		
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
		write(95,*)'      Nan detected'
		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
		if(EOFflag.eq.1)then
			write(*,*)' Every assemblage has been examined and no solution found'
			write(*,*)' Skipping this P,T point and moving on to the next'
			write(*,*)'T,P ',TC4,PB4
			go to 400
			endif		
!		xPhCo = xPhCoLast			! Reset to last good compositions
		call ResetToLast()			! added 3/28/2020
		go to 410
		endif
473	continue

	call SetLast()
	call SetNewP()
!	if(iParagonite.gt.0)then
	!	xPhCo(iParagonite,1) = 0.05
	!	xPhCo(iParagonite,2) = 0.95
!		endif		
400	continue			! loop back to get the next temperature
	ierr = 0		! no error
	return
	end
