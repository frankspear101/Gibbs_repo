! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine AutoTan2(logFileOn,iGrid,ierr)
    	use AWE_Interfaces
	use MyCanvas

!	************
!	This routine is exactly the same as AutoTan except:
!		When a solution isn't found, the program backs up and tries to find a solution at the new P&T by taking small steps
!		It hasn't been tested very well and doesn't seem to make much difference.....


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
	TYPE(AWE_Canvas) :: MyWindow		! The plotting canvas
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

	integer*4 i,j,iTP,iok,k,izero,isubtend,maxNegativeMolesPhase,maxgDifferencePhase,iTried
	integer*4 asmCurrentStore(phMax),kCur,newP,countTheMisses
	real*4 Tc4,Pb4,Pb4old,TCLast,xplot,yplot
	real*8 maxgDifference,maxNegativeMoles,gDiffThreshold
	real*8 sum
	real*8 TCsmall,PBsmall,deltaT,deltaP
	integer*4 logFileOn,iGrid,iremove
	integer*4 ierr,idumped,ierrsmall
	character*16 TPtext
!-----
!	these only needed for testing in block 46
	real*8 R
	DATA R/8.3144D0/

! ---------------------

!	We are using the solution from Sub StableASM(T,P) for the first TP point in the grid
!		as the starting point for calculations
!	asmCurrentLast = asmCurrent	! Store the last good assemblage
!	numCurrentLast = numCurrent
!	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
!	MPLast = MP0
!	asmCurrentNewP = asmCurrent	! Store the last good assemblage
!	numCurrentNewP = numCurrent
!	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
!	MPNewP = MP0
	call SetLast()
	call SetNewP()
!	PB4old = 0
	PB4old = PGrid(1)
	TC = TGrid(1)			! temporarily set TC to first point
	idumped = 0
!-----------------------------------------------
	write(95,*)' Starting TP loop'
	write(*,*)' Starting TP loop'
	write(95,401)TC4,PB4,(asmCurrent(k),k=1,numCurrent)			! log file

	xplot = Tgrid(numPTpoints) + 100
	yplot = (Pgrid(numPTpoints) + 1000)/1000.

	call OpenCanvasWindow(MyWindow)
	
	do 400 iTP = 2,numPTpoints
!	iLong(1) = 1		! turn on error reporting
!	iLong(2) = 1		! turn on output Jacobian
!	iLong(3) = 1		! turn on output thermo data
!	iLong(5) = 1		! turn on output MasterMatrix
!	iLong(7) = 1		! turn on output REXN
	newP = 0
	TCLast = TC		! store the previous PT point
	TC = Tgrid(iTP)
	PB = Pgrid(iTP)
	tc4 = tc
	pb4 = pb

!	write(*,*)TC4,PB4	! just to check progress
!	SUBROUTINE TextOnPlot(canvas,xplot,yplot,text,textSize)
	write(TPtext,390)TC4,PB4
390	format(2F8.1)
!	call TextOnPlot(XYPlot,xplot,yplot,TPtext,16)
	call TextOnMyWindow(MyWindow,50,50,TPtext,16)

	if(Dabs(PB4-PB4old).gt.1.)then		! this should print every time P changes - just to show progress
!		If we are holding the porosity constant, then we need to reset the H2O content to have it in excess here
		if(constantPorosity.eq.0.and.numFractionate.eq.0)then
!			write(12,*)TC4,PB4,(asmCurrent(iok),iok=1,numCurrent)
			write(*,*)TC4,PB4
			else
			if(iGrid.eq.1)then
			! reset the moles of H2O back to the starting value for the next row of calculations
			!	but only if we are doing a grid (not a P-T path)
				do 403 i=1,nc
				bulkCompMoles(i) = bulkCompMolesStart(i)
				molesStart(i) = bulkCompMolesStart(i)
				bulkCompWt(i)=molesStart(i)*molwt(i)/NumCatInOxide(i)
403				continue
				endif
			write(*,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)
402			format(2F12.2,20F12.5)
			endif
		PB4old = Pb4
		newP = 1
		endif
	if(newP.eq.1)then
!		If we are starting the next pressure, use the results from the previous P at the same T as a starting point
		call ResetToNewP()
!		call DumpEverything()
		else
!		Otherwise, use the last assemblage as the starting point
		call ResetToLast()
		endif
	if(constantPorosity.eq.0.and.numFractionate.eq.0)then
		write(95,401)TC4,PB4,(asmCurrent(k),k=1,numCurrent)			! log file
401		format(2F12.2,20I5)
		else
		write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
!		Calculate the total mass and moles of the bulk composition
		BCmoles = 0.0d0
		BCmass = 0.0d0
		do 413 i = 1,nc
		BCmoles = BCmoles + bulkCompMoles(i)
		BCmass = BCmass + bulkCompWt(i)
413		continue
!		calculate the total volume of the assemblage as the sum of the volumes of each phase
		BCVol = 0.0d0
		do 415 kcur = 1,numCurrent
		k = asmCurrent(kCur)
		BCVol = BCVol + VP0(k)
415		continue
		call WriteBCMoles
		call WriteBCWt
		endif
!	Note that a call to duPTX to get uZero also calls activity (this could/should be changed)
!	So, if compositions are off for some reason, the code fails due to the activity calculation (not the uzero calculation)
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	!write(*,*)uzero(1,1),gattp(1,1)
	if(izero.gt.0)then
		write(95,*)' Error calculating MU zero for this T&P... abort the run'
		write(*,*)' Error calculating MU zero for this T&P... abort the run'
		write(*,*)'izero = ',izero
	!	call DumpEverything()
		close(95)
		close(74)
		return
		endif
	countTheMisses = 0
410	continue

	if(countTheMisses.gt.30) then
		if(newP.eq.0)then
			call ResetToLast()
			else					! this calculation is at the new pressure so reset accordingly
			call ResetToNewP()
			endif
		if(tanout(8).eq.1)call dumpEverything()
		go to 400			! start on a new T&P
		endif	

	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
		write(95,*)' Assemblage does not subtend BC (line 468 in Sub AutoTan) - get new asm'
		Write(*,*)' Assemblage does not subtend BC (line 468 in Sub AutoTan)'
!		Write(*,*)' This should never happen -- (already checked for Combo file)'
!		Write(*,*)' Try restarting the program....'		
!		pause 'hit return to abort '
!		ierr = 1
!		close(95)
!		close(74)
!		return
		endif
	itried = 0
414	continue
	izero = 0
	call ComputePTXM(TC4,PB4,izero)
	!write(*,*)uzero(1,1),gattp(1,1)
	if(izero.gt.0)then
		if(newP.eq.0)then
			call ResetToLast()
			else					! this calculation is at the new pressure so reset accordingly
			call ResetToNewP()
			endif
		write(95,*)'	Calculation failed in ComputePTXM in first call at this T&P '
		if(iTried.lt.1)then
			iTried = itried + 1
			write(95,*)' Try again with a smaller newtonStep'
			newtonStepMax = 0.01
			newtonStep    = 0.01
			go to 414
			endif
		write(95,*)'   Jump to work on subassemblages'
		newtonStepMax = 0.4
		newtonStep    = 0.4
		go to 500		! This last assemblage didn't work...try working on subassemblages

!		go to 400		! This last assemblage didn't work... Try the next temperature
		endif
!		Continue calculations
425	continue
!	If here, then we have found a solution using this assemblage

!	Check to see if any phases have negative moles - 
!		if so then find the largest negative moles and remove that phase
	!write(*,*)uzero(1,1),gattp(1,1)
	!write(*,*)' Check for negative moles'
	maxNegativeMoles = 0.0d0
	maxNegativeMolesPhase = 0
	do 430 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	if(mp0(k)/atomNorm(k,1).lt.maxNegativeMoles)then
!		A phase has negative moles - record and loop back
		maxNegativeMoles = mp0(k)*atomNorm(k,1)		! we are storing moles/atom (not per molecule)
		maxNegativeMolesPhase = kCur
!		write(12,*)'  '
!		write(12,*)' &*&*&*&*&*&*&*&*&*&*&*&*&*'
!		write(12,*)'MaxNegativeMoles = ',maxNegativeMoles,k,kcur,mp0(k),atomNorm(k,1)
		endif
430	continue


	if(maxNegativeMolesPhase.ne.0)then
!		call DumpEverything()
		iremove = asmCurrent(maxNegativeMolesPhase)
!		if(iremove.eq.4)call DumpEverything()			
		mp0(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		do 432 j = maxNegativeMolesPhase,numPh
		asmCurrent(j) = asmCurrent(j+1)
432		continue
		numCurrent = numCurrent - 1
		write(95,*)'Removed phase= ',iremove,'asm ',(asmCurrent(iok),iok=1,numCurrent)	!(MinRec(asmCurrent(iok)),iok=1,numCurrent)  
!		write(*,*)'Removed phase= ',iremove,'asm ',(asmCurrent(iok),iok=1,numCurrent)	!(MinRec(asmCurrent(iok)),iok=1,numCurrent)  
		countTheMisses = countTheMisses + 1
!	Don't reset to last good TP. Rather, use the compositions from the last calculation before the remove
!		if(newP.eq.0)then
!			xPhCo = xPhCoLast			! Reset to last good compositions! these should be the ones from the last T&P
!			MP0 = MPLast
!			else					! this calculation is at the new pressure so reset accordingly
!			xPhCo = xPhCoNewP
!			MP0 = MPNewP
!			endif


		go to 410
		endif
!	if here, then no phase has negative moles
	!write(*,*)' Done checking for negative moles'
	!write(*,*)uzero(1,1),gattp(1,1)

	call AdjustTangentToAsm
!	call listTangentPlane()
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
!	write(*,*)'Quartz uOnTanAU,uOnTanDeltaAU: ',uOnTanAU(1,1),uOnTanDeltaAU(1,1)
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	Check to see if any phases are below the tangent plane
435	continue
	!write(*,*)' Check for a phase below the tangent'
	maxgDifference = -1.0d-2
!	maxgDifference = -1.0d0
!	If a phase is more than 0.01 J below the tangent, then flag it.
	gDiffThreshold = -1.0e-2
	maxgDifferencePhase = 0
	!write(*,*)' --------------------'
	do 440 k = 1,numPhMIF
	!write(*,*)Phname(k),gDifferenceAU(k)
	if(k.eq.iExclude)go to 440			! skip the excluded phase even if it is below the tangent
	if(gDifferenceAU(k).lt.gDiffThreshold)then	! if true then the phase lies below the tangent
		if(gDifferenceAU(k).lt.maxgDifference)then			
			maxgDifference = gDifferenceAU(k)
			maxgDifferencePhase = k
			endif
		endif
440	continue
	if(maxgDifferencePhase.ne.0)then		

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
		mp0(numCurrent) = 0.001		! make the moles of this new phase non-zero to avoid possible singularities
!		write(95,*)'Phase added    ',(minRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		write(95,*)'Phase added    ',(asmCurrent(iok),iok=1,numCurrent),maxGdifference
!		write(*,*)'Phase added    ',(asmCurrent(iok),iok=1,numCurrent)
		countTheMisses = countTheMisses + 1
		isubtend = 1
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)		! this one MUST subtend because the last one did and we added a phase
		izero = 0
		call ComputePTXM(TC4,PB4,izero)
		if(izero.eq.0)then
			go to 425			! this assemblage worked - go back and check for negative moles
			endif
		else
		go to 475			! no phase sits below tangent -- we should be good
		endif

!	This next section should be done in a subroutine working on subassemblages.... but this might work (just a bit of spagetti)
500		continue
!		If we get here, then the last assemblage didn't work (izero not = 0)
		write(95,*)'*****************************************************************'
		write(95,*)' Start subassemblage routine'
		write(95,*)'*****************************************************************'
		write(95,*)'	Calculation failed in ComputePTXM after adding phase - Work on subassemblages'
		write(95,*)'	Current assemblage: ',(asmCurrent(iok),iok=1,numCurrent)

!		write(*,*)'	Calculation failed in ComputePTXM (line 269 in Sub AutoTan) - Work on subassemblages'
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
	!	write(*,*)' icombo, numCurrentStore, comboCounter ',icombo, numCurrentStore, comboCounter
		
		do 450 i = 1,comboCounter		
		if(newP.eq.0)then
			call ResetToLast()
!			xPhCo = xPhCoLast			! Reset to last good compositions
!								! these should be the ones from the last T&P
!			MP0 = MPLast
			else					! this calculation is at the new pressure so reset accordingly
			call ResetToNewP()
!			xPhCo = xPhCoNewP
!			MP0 = MPNewP
			endif
!		xPhCo = xPhCoLast				! Reset to last good compositions
!		MP0 = MPLast
!		if(iParagonite.gt.0)then
		!	xPhCo(iParagonite,1) = 0.05
		!	xPhCo(iParagonite,2) = 0.95
!			endif		
		numCurrent = icombo		          !the number of phases in the assemblage
		write(95,*)i,numCurrent,(Combos(i,k),k=1,numCurrent)
		!write(*,*)i,numCurrent,(Combos(i,k),k=1,numCurrent)
		do k = 1,numCurrent
			asmCurrent(k) = asmCurrentStore(Combos(i,k))
			end do
		write(95,*)'	 Subassemblage   ',(asmCurrent(iok),iok=1,numCurrent)
!		write(95,*)'	 Subassemblage  ',numCurrent,(MinRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		!write(*,*)'	 Subassemblage  ',numCurrent,(MinRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		isubtend = 1
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)
		if(isubtend.eq.0)then
			write(95,*)'		Assemblage does not subtend BC (line 605 in Sub Tangent) - get new asm'
			go to 450
			endif
		izero = 0
		call ComputePTXM(TC4,PB4,izero)
		if(izero.gt.0)then
			write(95,*)'		   Error in ComputePTXM on this subassemblage - get new asm'
			write(95,*)'izero = ',izero
			!write(*,*)'		   Error in ComputePTXM (Line 308 of Sub AutoTan) - get new asm'
			!write(*,*)'izero = ',izero
			go to 450			! get the next assemblage
			endif
!		If we get to here, then we have found a subassemblage that works.
!		Check for negative moles, etc.
		do kCur = 1,numCurrent
			k = asmCurrent(kCur)
			if(mp0(k)*atomNorm(k,1).le.-1.0d-10)then
				!A phase has negative moles - go check out next assemblage
				write(95,*)'		   Phase has negative moles = ',k
				!write(*,*)'		   Phase has negative moles = ',k
				go to 450
				endif
			end do

!		if here, then no phase has negative moles
!		Check to see if any phases are below the tangent plane
		call AdjustTangentToAsm
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
		do k = 1,numPhMIF
			if(gDifferenceAU(k).lt.gDiffThreshold)then			
				write(95,*)        '		   Phase lies below tangent = ',k
				if(iExclude.eq.k)then
					write(95,*)'		   This is the excluded phase--skip'
					go to 462
					endif			
				go to 450		! get the next assemblage - since this is a subassemblage, don't add a new phase
! 				do j = 1,numCurrentStore
! 					if(k.eq.asmCurrentStore(j))go to 450		! this phase is already being considered. Don't reset anything
! 					end do
! 				!It appears that this is a new phase. Try adding that phase in and going through the asmCombos again
! 				!Make sure this isn't the excluded phase
! 				numCurrentStore = numCurrentStore + 1
! 				asmCurrentStore(numCurrentStore) = k
! 				numCurrent = numCurrentStore
! 				go to 452		! Start the subassemblage routine all over
				endif
462			continue
			end do
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

!		If we get to here then we have gone through all of the subassemblages in the ASMCompo array and none has worked
!		Loop back and try again reading assemblages from the main Compo file
		if(newP.eq.0)then
			call ResetToLast()
!			xPhCo = xPhCoLast			! Reset to last good compositions
!								! these should be the ones from the last T&P
!			MP0 = MPLast
!			asmCurrent = asmCurrentLast		! set the assemblage to the last good one
!			numCurrent = numCurrentLast
			else					! this calculation is at the new pressure so reset accordingly
			call ResetToNewP()
!			xPhCo = xPhCoNewP
!			MP0 = MPNewP
!			asmCurrent = asmCurrentNewP		! set the assemblage to the low T, new P values
!			numCurrent = numCurrentNewP
			endif

		! if we get here, then we haven't found a solution at this new T and P
		! Try to find a solution using smaller increments of T or P.
		
		!set values to the previous T and P
		TCsmall = Tgrid(iTP-1)
		PBsmall = Pgrid(iTP-1)
		deltaT = 0.01d0		! Temperature increment (C)
		deltaP = 1.d0		! pressure increment (bars)
		if(newP.eq.0)then
			! this is a new temperature
			TCsmall = TCsmall + deltaT
			else
			PBsmall = PBsmall + deltaP
			endif
		call AutoTanSmall(TCsmall,PBsmall,ierrSmall)
		write(*,*)' Current assemblage'
		write(*,*)numCurrent,(asmCurrent(i),i=1,numCurrent)
		call WriteAllStuffToTerminal(6)

		if(ierrSmall.eq.1)then
			write(*,*)'The small did not work'
			write(*,*)TCsmall,PBsmall
			!pause 'hit return to continue'
			go to 400	! go get the next PT point on the grid
			endif
		if(abs(TCsmall-TC).lt.0.01.or.abs(PBsmall-Pb).lt.1)then
			! we're at the new grid T and P. 
			! This should be a good solution (no error)
			go to 475				
			endif
				
!		go to 400		! go on to the next TP point

!		endif

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
		!pause 'NAN detected... hit return to continue'
!		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
!		if(EOFflag.eq.1)then
!			write(*,*)' Every assemblage has been examined and no solution found'
!			write(*,*)' Skipping this P,T point and moving on to the next'
!			write(*,*)'T,P ',TC4,PB4
!			go to 400
!			endif		
		xPhCo = xPhCoLast			! Reset to last good compositions
		go to 400	! go get a new T,P
		endif
	call WriteAllStuff		! write results to output file 74
	write(95,404)(asmCurrent(k),k=1,numCurrent)			! log file
404	format(24x,20I5)
!	if we are holding the porosity constant, then here is where we remove sufficient H2O to maintain constant porosity
	if(constantPorosity.eq.1)then
		call AdjustWater(sum)		! sum is returned from subroutine CalculateMode. It is the total volume of the current assemblage
!		write(95,4402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
		endif
!----- Fractionation------
	if(numFractionate.gt.0)then
		Call Fractionate(sum)		! sum is returned from subroutine CalculateMode. It is the total volume of the current assemblage
!		write(95,4402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
		endif

	call SetLast()
	if(newP.eq.1)then		! We have a solution for the initial T for the new pressure - save
		call SetNewP()
		endif
400	continue			! loop back to get the next temperature
	close (74)
	if(logFileOn.eq.1)close (95)
	Close(87)		! combo file

	call AWE_closeCanvas(MyWindow)

	write(*,*)' '
	write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	write(*,*)'All done with calculations.....'
	write(*,*)' All phases (stable and metastable) are in the file:',trim(tangentOutputFile)//'.ALL'
	write(*,*)' Debugging information is in file:                  ',trim(tangentOutputFile)//'.log'
	write(*,*)' '
	write(*,*)' Run program MADPlotter_2 to plot the MAD diagram'
	write(*,*)' '
	write(*,*)' '
	pause ' hit return to continue'

	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine AutoTanSmall(TCsmall,PBsmall,ierr)
    	use AWE_Interfaces
	use MyCanvas

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
	TYPE(AWE_Canvas) :: MyWindow		! The plotting canvas
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

	integer*4 i,j,iok,k,izero,isubtend,maxNegativeMolesPhase,maxgDifferencePhase,iTried
	integer*4 asmCurrentStore(phMax),kCur,newP,countTheMisses
	real*4 Tc4,Pb4
	real*8 maxgDifference,maxNegativeMoles,gDiffThreshold
	real*8 TCsmall,PBsmall
	integer*4 logFileOn,iremove
	integer*4 ierr
	character*16 TPtext
!-----
!	these only needed for testing in block 46
	real*8 R
	DATA R/8.3144D0/

! ---------------------

! !	We are using the solution from Sub StableASM(T,P) for the first TP point in the grid
! !		as the starting point for calculations
! !	asmCurrentLast = asmCurrent	! Store the last good assemblage
! !	numCurrentLast = numCurrent
! !	xPhCoLast = xPhCo		! Store the last set of good phase compositions for resetting if needed
! !	MPLast = MP0
! !	asmCurrentNewP = asmCurrent	! Store the last good assemblage
! !	numCurrentNewP = numCurrent
! !	xPhCoNewP = xPhCo		! Store the last set of good phase compositions for resetting if needed
! !	MPNewP = MP0
! 	call SetLast()
! 	call SetNewP()
! !	PB4old = 0
! 	PB4old = PGrid(1)
! 	TC = TGrid(1)			! temporarily set TC to first point
! 	idumped = 0
! !-----------------------------------------------
! 	write(95,*)' Starting TP loop'
! 	write(*,*)' Starting TP loop'
! 	write(95,401)TC4,PB4,(asmCurrent(k),k=1,numCurrent)			! log file
! 
! 	xplot = Tgrid(numPTpoints) + 100
! 	yplot = (Pgrid(numPTpoints) + 1000)/1000.
! 
! 	call OpenCanvasWindow(MyWindow)
! 	
! 	do 400 iTP = 2,numPTpoints
!	iLong(1) = 1		! turn on error reporting
!	iLong(2) = 1		! turn on output Jacobian
!	iLong(3) = 1		! turn on output thermo data
!	iLong(5) = 1		! turn on output MasterMatrix
!	iLong(7) = 1		! turn on output REXN
! 	newP = 0
! 	TCLast = TC		! store the previous PT point
! 	TC = Tgrid(iTP)
! 	PB = Pgrid(iTP)
 	TC = TCsmall
 	PB = PBsmall
	tc4 = tc
	pb4 = pb

!	write(*,*)TC4,PB4	! just to check progress
!	SUBROUTINE TextOnPlot(canvas,xplot,yplot,text,textSize)
	write(TPtext,390)TC4,PB4
390	format(2F8.1)
!	call TextOnPlot(XYPlot,xplot,yplot,TPtext,16)
	call TextOnMyWindow(MyWindow,50,50,TPtext,16)

! 	if(Dabs(PB4-PB4old).gt.1.)then		! this should print every time P changes - just to show progress
! !		If we are holding the porosity constant, then we need to reset the H2O content to have it in excess here
! 		if(constantPorosity.eq.0.and.numFractionate.eq.0)then
! !			write(12,*)TC4,PB4,(asmCurrent(iok),iok=1,numCurrent)
! 			write(*,*)TC4,PB4
! 			else
! 			if(iGrid.eq.1)then
! 			! reset the moles of H2O back to the starting value for the next row of calculations
! 			!	but only if we are doing a grid (not a P-T path)
! 				do 403 i=1,nc
! 				bulkCompMoles(i) = bulkCompMolesStart(i)
! 				molesStart(i) = bulkCompMolesStart(i)
! 				bulkCompWt(i)=molesStart(i)*molwt(i)/NumCatInOxide(i)
! 403				continue
! 				endif
! 			write(*,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)
! 402			format(2F12.2,20F12.5)
! 			endif
! 		PB4old = Pb4
! 		newP = 1
! 		endif
! 	if(newP.eq.1)then
! !		If we are starting the next pressure, use the results from the previous P at the same T as a starting point
! 		call ResetToNewP()
! !		call DumpEverything()
! 		else
! !		Otherwise, use the last assemblage as the starting point
! 		call ResetToLast()
! 		endif
! 	if(constantPorosity.eq.0.and.numFractionate.eq.0)then
! 		write(95,401)TC4,PB4,(asmCurrent(k),k=1,numCurrent)			! log file
! 401		format(2F12.2,20I5)
! 		else
! 		write(95,402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
! !		Calculate the total mass and moles of the bulk composition
! 		BCmoles = 0.0d0
! 		BCmass = 0.0d0
! 		do 413 i = 1,nc
! 		BCmoles = BCmoles + bulkCompMoles(i)
! 		BCmass = BCmass + bulkCompWt(i)
! 413		continue
! !		calculate the total volume of the assemblage as the sum of the volumes of each phase
! 		BCVol = 0.0d0
! 		do 415 kcur = 1,numCurrent
! 		k = asmCurrent(kCur)
! 		BCVol = BCVol + VP0(k)
! 415		continue
! 		call WriteBCMoles
! 		call WriteBCWt
! 		endif
!	Note that a call to duPTX to get uZero also calls activity (this could/should be changed)
!	So, if compositions are off for some reason, the code fails due to the activity calculation (not the uzero calculation)
	call CalcMuZero(izero)				! do only once at T&P (note that this resets asmCurrent to include all phases)
	!write(*,*)uzero(1,1),gattp(1,1)
	if(izero.gt.0)then
		write(95,*)' Error calculating MU zero for this T&P... abort the run'
		write(*,*)' Error calculating MU zero for this T&P... abort the run'
		write(*,*)'izero = ',izero
	!	call DumpEverything()
		close(95)
		close(74)
		return
		endif
	countTheMisses = 0
410	continue

! 	if(countTheMisses.gt.30) then
! 		if(newP.eq.0)then
! 			call ResetToLast()
! 			else					! this calculation is at the new pressure so reset accordingly
! 			call ResetToNewP()
! 			endif
! 		if(tanout(8).eq.1)call dumpEverything()
! 
! 		go to 400			! start on a new T&P
! 
! 		endif	

	isubtend = 1
	call LoadCurrentAsm(TC4,PB4,izero,isubtend)
	if(isubtend.eq.0)then
		write(95,*)' Assemblage does not subtend BC (line 468 in Sub AutoTan) - get new asm'
		Write(*,*)' Assemblage does not subtend BC (line 468 in Sub AutoTan)'
!		Write(*,*)' This should never happen -- (already checked for Combo file)'
!		Write(*,*)' Try restarting the program....'		
!		pause 'hit return to abort '
!		ierr = 1
!		close(95)
!		close(74)
!		return
		endif
	itried = 0
414	continue
	izero = 0
	call ComputePTXM(TC4,PB4,izero)
	!write(*,*)uzero(1,1),gattp(1,1)
	if(izero.gt.0)then
		if(newP.eq.0)then
			call ResetToLast()
			else					! this calculation is at the new pressure so reset accordingly
			call ResetToNewP()
			endif
		write(95,*)'	Calculation failed in ComputePTXM in first call at this T&P '
		if(iTried.lt.1)then
			iTried = itried + 1
			write(95,*)' Try again with a smaller newtonStep'
			newtonStepMax = 0.01
			newtonStep    = 0.01
			go to 414
			endif
		write(95,*)'   Jump to work on subassemblages'
		newtonStepMax = 0.4
		newtonStep    = 0.4
		go to 500		! This last assemblage didn't work...try working on subassemblages

!		go to 400		! This last assemblage didn't work... Try the next temperature
		endif
!		Continue calculations
425	continue
!	If here, then we have found a solution using this assemblage

!	Check to see if any phases have negative moles - 
!		if so then find the largest negative moles and remove that phase
	!write(*,*)uzero(1,1),gattp(1,1)
	!write(*,*)' Check for negative moles'
	maxNegativeMoles = 0.0d0
	maxNegativeMolesPhase = 0
	do 430 kCur = 1,numCurrent
	k = asmCurrent(kCur)
	if(mp0(k)/atomNorm(k,1).lt.maxNegativeMoles)then
!		A phase has negative moles - record and loop back
		maxNegativeMoles = mp0(k)*atomNorm(k,1)		! we are storing moles/atom (not per molecule)
		maxNegativeMolesPhase = kCur
!		write(12,*)'  '
!		write(12,*)' &*&*&*&*&*&*&*&*&*&*&*&*&*'
!		write(12,*)'MaxNegativeMoles = ',maxNegativeMoles,k,kcur,mp0(k),atomNorm(k,1)
		endif
430	continue


	if(maxNegativeMolesPhase.ne.0)then
!		call DumpEverything()
		iremove = asmCurrent(maxNegativeMolesPhase)
!		if(iremove.eq.4)call DumpEverything()			
		mp0(asmCurrent(maxNegativeMolesPhase)) = 0.0d0		
		do 432 j = maxNegativeMolesPhase,numPh
		asmCurrent(j) = asmCurrent(j+1)
432		continue
		numCurrent = numCurrent - 1
		write(95,*)'Removed phase= ',iremove,'asm ',(asmCurrent(iok),iok=1,numCurrent)	!(MinRec(asmCurrent(iok)),iok=1,numCurrent)  
!		write(*,*)'Removed phase= ',iremove,'asm ',(asmCurrent(iok),iok=1,numCurrent)	!(MinRec(asmCurrent(iok)),iok=1,numCurrent)  
		countTheMisses = countTheMisses + 1
!	Don't reset to last good TP. Rather, use the compositions from the last calculation before the remove
!		if(newP.eq.0)then
!			xPhCo = xPhCoLast			! Reset to last good compositions! these should be the ones from the last T&P
!			MP0 = MPLast
!			else					! this calculation is at the new pressure so reset accordingly
!			xPhCo = xPhCoNewP
!			MP0 = MPNewP
!			endif


		go to 410
		endif
!	if here, then no phase has negative moles
	!write(*,*)' Done checking for negative moles'
	!write(*,*)uzero(1,1),gattp(1,1)

	call AdjustTangentToAsm
!	call listTangentPlane()
	call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
!	write(*,*)'Quartz uOnTanAU,uOnTanDeltaAU: ',uOnTanAU(1,1),uOnTanDeltaAU(1,1)
	call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
!	Check to see if any phases are below the tangent plane
435	continue
	!write(*,*)' Check for a phase below the tangent'
	maxgDifference = -1.0d-2
!	maxgDifference = -1.0d0
!	If a phase is more than 0.01 J below the tangent, then flag it.
	gDiffThreshold = -1.0e-2
	maxgDifferencePhase = 0
	!write(*,*)' --------------------'
	do 440 k = 1,numPhMIF
	!write(*,*)Phname(k),gDifferenceAU(k)
	if(k.eq.iExclude)go to 440			! skip the excluded phase even if it is below the tangent
	if(gDifferenceAU(k).lt.gDiffThreshold)then	! if true then the phase lies below the tangent
		if(gDifferenceAU(k).lt.maxgDifference)then			
			maxgDifference = gDifferenceAU(k)
			maxgDifferencePhase = k
			endif
		endif
440	continue
	if(maxgDifferencePhase.ne.0)then		

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
		mp0(numCurrent) = 0.001		! make the moles of this new phase non-zero to avoid possible singularities
!		write(95,*)'Phase added    ',(minRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		write(95,*)'Phase added    ',(asmCurrent(iok),iok=1,numCurrent),maxGdifference
!		write(*,*)'Phase added    ',(asmCurrent(iok),iok=1,numCurrent)
		countTheMisses = countTheMisses + 1
		isubtend = 1
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)		! this one MUST subtend because the last one did and we added a phase
		izero = 0
		call ComputePTXM(TC4,PB4,izero)
		if(izero.eq.0)then
			go to 425			! this assemblage worked - go back and check for negative moles
			endif
		else
		go to 475			! no phase sits below tangent -- we should be good
		endif

!	This next section should be done in a subroutine working on subassemblages.... but this might work (just a bit of spagetti)
500		continue
!		If we get here, then the last assemblage didn't work (izero not = 0)
		write(95,*)'*****************************************************************'
		write(95,*)' Start subassemblage routine'
		write(95,*)'*****************************************************************'
		write(95,*)'	Calculation failed in ComputePTXM after adding phase - Work on subassemblages'
		write(95,*)'	Current assemblage: ',(asmCurrent(iok),iok=1,numCurrent)

!		write(*,*)'	Calculation failed in ComputePTXM (line 269 in Sub AutoTan) - Work on subassemblages'
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
	!	write(*,*)' icombo, numCurrentStore, comboCounter ',icombo, numCurrentStore, comboCounter
		
		do 450 i = 1,comboCounter		
		if(newP.eq.0)then
			call ResetToLast()
!			xPhCo = xPhCoLast			! Reset to last good compositions
!								! these should be the ones from the last T&P
!			MP0 = MPLast
			else					! this calculation is at the new pressure so reset accordingly
			call ResetToNewP()
!			xPhCo = xPhCoNewP
!			MP0 = MPNewP
			endif
!		xPhCo = xPhCoLast				! Reset to last good compositions
!		MP0 = MPLast
!		if(iParagonite.gt.0)then
		!	xPhCo(iParagonite,1) = 0.05
		!	xPhCo(iParagonite,2) = 0.95
!			endif		
		numCurrent = icombo		          !the number of phases in the assemblage
		write(95,*)i,numCurrent,(Combos(i,k),k=1,numCurrent)
		!write(*,*)i,numCurrent,(Combos(i,k),k=1,numCurrent)
		do k = 1,numCurrent
			asmCurrent(k) = asmCurrentStore(Combos(i,k))
			end do
		write(95,*)'	 Subassemblage   ',(asmCurrent(iok),iok=1,numCurrent)
!		write(95,*)'	 Subassemblage  ',numCurrent,(MinRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		!write(*,*)'	 Subassemblage  ',numCurrent,(MinRec(asmCurrent(iok)),iok=1,numCurrent),(asmCurrent(iok),iok=1,numCurrent)
		isubtend = 1
		call LoadCurrentAsm(TC4,PB4,izero,isubtend)
		if(isubtend.eq.0)then
			write(95,*)'		Assemblage does not subtend BC (line 605 in Sub Tangent) - get new asm'
			go to 450
			endif
		izero = 0
		call ComputePTXM(TC4,PB4,izero)
		if(izero.gt.0)then
			write(95,*)'		   Error in ComputePTXM on this subassemblage - get new asm'
			write(95,*)'izero = ',izero
			!write(*,*)'		   Error in ComputePTXM (Line 308 of Sub AutoTan) - get new asm'
			!write(*,*)'izero = ',izero
			go to 450			! get the next assemblage
			endif
!		If we get to here, then we have found a subassemblage that works.
!		Check for negative moles, etc.
		do 460 kCur = 1,numCurrent
		k = asmCurrent(kCur)
		if(mp0(k)*atomNorm(k,1).le.-1.0d-10)then
!			A phase has negative moles - go check out next assemblage
			write(95,*)'		   Phase has negative moles = ',k
			!write(*,*)'		   Phase has negative moles = ',k
			go to 450
			endif
460		continue

!		if here, then no phase has negative moles
!		Check to see if any phases are below the tangent plane
		call AdjustTangentToAsm
		call CalcMuOnTangent		! calculates mu of each phase component on the current tangent plane
		call ParalleltoTangent(logFileOn)	! determines (a) composition closest to tangent (b) difference in G between this phase and tangent	
		do k = 1,numPhMIF
			if(gDifferenceAU(k).lt.gDiffThreshold)then			
				write(95,*)'		   Phase lies below tangent = ',k
				if(iExclude.eq.k)then
					write(95,*)'		   This is the excluded phase--skip'
					go to 462
					endif			
				go to 450		! get the next assemblage - since this is a subassemblage, don't add a new phase
! 				do j = 1,numCurrentStore
! 					if(k.eq.asmCurrentStore(j))go to 450		! this phase is already being considered. Don't reset anything
! 					end do
! !			It appears that this is a new phase. Try adding that phase in and going through the asmCombos again
! !			Make sure this isn't the excluded phase
! 				numCurrentStore = numCurrentStore + 1
! 				asmCurrentStore(numCurrentStore) = k
! 				numCurrent = numCurrentStore
! 				go to 452		! Start the subassemblage routine all over
				endif
462			continue
			end do
		
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

	! If we get to here then we have gone through all of the subassemblages in the ASMCompo array and none has worked

	! return with error = 1
	ierr = 1
	return

475	continue
!	If we get to here, then no phases sit below the tangent and we should be done

	! return with error = 0
	ierr = 0
	return
	
! 	sum = 0.0d0
! 	call CalculateMode(sum)		! sum is the sum of the volumes of phases in cm^3
! 	call CalculateGSystem
! !	It seems that some assemblages result in NaN for a number of variables.
! !	We need to check this and cycle back if we have NaN rather than a valid solution
! 	if(ieee_is_nan(gSystem))then
! 	!	write(*,*)'NaN detected'
! 		write(95,*)' Nan detected'
! 		!pause 'NAN detected... hit return to continue'
! !		call GetNewAssemblage(asmCurrent,numCurrent,logFileOn,iExclude,EOFflag)
! !		if(EOFflag.eq.1)then
! !			write(*,*)' Every assemblage has been examined and no solution found'
! !			write(*,*)' Skipping this P,T point and moving on to the next'
! !			write(*,*)'T,P ',TC4,PB4
! !			go to 400
! !			endif		
! 		xPhCo = xPhCoLast			! Reset to last good compositions
! 		go to 400	! go get a new T,P
! 		endif
! 	call WriteAllStuff		! write results to output file 74
! 	write(95,404)(asmCurrent(k),k=1,numCurrent)			! log file
! 404	format(24x,20I5)
! !	if we are holding the porosity constant, then here is where we remove sufficient H2O to maintain constant porosity
! 	if(constantPorosity.eq.1)then
! 		call AdjustWater(sum)		! sum is returned from subroutine CalculateMode. It is the total volume of the current assemblage
! !		write(95,4402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
! 		endif
! !----- Fractionation------
! 	if(numFractionate.gt.0)then
! 		Call Fractionate(sum)		! sum is returned from subroutine CalculateMode. It is the total volume of the current assemblage
! !		write(95,4402)TC4,PB4,(bulkCompmoles(i),i=1,nc)			! log file
! 		endif
! 
! 	call SetLast()
! 	if(newP.eq.1)then		! We have a solution for the initial T for the new pressure - save
! 		call SetNewP()
! 		endif
! 400	continue			! loop back to get the next temperature
! 	close (74)
! 	if(logFileOn.eq.1)close (95)
! 	Close(87)		! combo file
! 
! 	call AWE_closeCanvas(MyWindow)
! 
! 	write(*,*)' '
! 	write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
! 	write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
! 	write(*,*)'All done with calculations.....'
! 	write(*,*)' All phases (stable and metastable) are in the file:',trim(tangentOutputFile)//'.ALL'
! 	write(*,*)' Debugging information is in file:                  ',trim(tangentOutputFile)//'.log'
! 	write(*,*)' '
! 	write(*,*)' Run program MADPlotter_2 to plot the MAD diagram'
! 	write(*,*)' '
! 	write(*,*)' '
! 	pause ' hit return to continue'
! 
! 	return
	end

