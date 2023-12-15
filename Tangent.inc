	integer*4 tanOut(20)
!	tanOut(1) = Output results of ParallelToTangent (summary)
!	tanOut(2) = Output every iteration of ParallelToTangent
!	tanOut(3) = Pause after every phase ParallelToTangent
!	tanOut(4) = Call to ThermoData with new assemblage
!	tanOut(5) = Output transformation matrix (B) - phase to system components
!	tanOut(6) = 
!	tanOut(7) = 
!	tanOut(8) =  
!	tanOut(9) =  
!	tanOut(10) = 
!	tanOut(11) = 
!	tanOut(12) = 
!	tanOut(13) = 
!	tanOut(14) = 
!
!
	integer*4 maxGridPT
!	parameter (maxGridPT = 20000)
	parameter (maxGridPT = 500000)
	real*8 tanPlane(sysCoMax)
	real*8 uZero(phMax,phCoMax),uZeroDelta(phMax,phCoMax),uZeroAU(phMax,phCoMax),uZeroDeltaAU(phMax,phCoMax), &
     		uOnTanAU(phMax,phCoMax),uOnTanDeltaAU(phMax,phCoMax)
	real*8 gSystem,gOnTanAU(phMax),gOfPhaseAU(phMax),gDifferenceAU(phMax)
	real*8 sysCoMoleFraction(sysCoMax)
	integer*4 numVar,numEq,numCurrent,numCurrentLast,numCurrentNewP,ioutputPseudo,idoingpseudo,ipseudo,idegen
	integer*4 ignorePhase(phMax),iExclude,OSPhase(phMax),numOSPhases
	real*8 molesCurrent(phMax),gofPhaseCurrent(phMax),gofCurrentAsm
	real*8 Tgrid(0:maxGridPT),Pgrid(0:maxGridPT),PFluidgrid(0:maxGridPT)
	integer*4 numPTpoints
	character*16 varNames(phCoMax)
	character*4 equationNames(phCoMax)
	character*255 tangentOutputFile,logFile,stuffFile,allFile,BCMolesFile,BCWtFile,FractFile,MeltFile
!	integer*2 ASMCombo(14,10000,15),numASMCOMBOS(14)
	real*4 TcStart,TcEnd,Tinc,PbStart,PbEnd,Pinc
	real*8 residual(phMax),threshold(phMax),removeVol(phMax),removeMoles(phMax),MeltWt(sysCoMax),     &
		MeltMoles,MeltVol,BCmoles,BCMass,BCVol,porosity
	integer*4 numFractionate,toFractionate(phMax),constantPorosity,fluidIndex,iparagonite
	common/tangentCom/tanOut,  &
     		tanPlane, &
     		uZero,uZeroDelta,uZeroau,uZeroDeltaau,uOnTanau,uOnTanDeltaau, &
     		gSystem,gOnTanau,gOfPhaseau,gDifferenceau, &
     		sysCoMoleFraction, &
     		numVar,numEq,numCurrent,numCurrentLast,numCurrentNewP,ioutputPseudo,idoingpseudo,ipseudo,idegen, &
     	  	ignorePhase,iExclude,OSPhase,numOSPhases,		 	&
     		molesCurrent,gofPhaseCurrent,gofCurrentAsm, 			&
     		Tgrid,Pgrid,PFluidgrid,numPTpoints,				 &
     		varNames,equationNames,tangentOutputFile,logFile,stuffFile,allFile,BCMolesFile,BCWtFile,FractFile,MeltFile,   &
!     		ASMCombo,numAsmCombos,	&
     		TcStart,TcEnd,Tinc,PbStart,PbEnd,Pinc,			&
     		residual,threshold,removeVol,removeMoles,				&
     		numFractionate,toFractionate,MeltWt,constantPorosity,fluidIndex,iparagonite,	&
     		MeltMoles,MeltVol,BCMoles,BCMass,BCVol,porosity
	
	! variables with "au" at the end are in atom units. All others in molar units
	! uZerodelta(J) holds for each phase the difference between the independent and dependent phase components (j-1)
	!      note that storage goes from 2 to numPhCo(K) so that the first one is a dummy (actually 0.0)
	! uZero(k,j) stores the value of zero for each phase component at the T and P
	!		This is fixed after initialization of routine
	! uonTan(j) stores the value of the tangent plane where it hits a phase component composition
	!		This needs to be updated whenever the tangent plane is reset
	!		Storage goes from 1 to numPhCo(K)
	! uonTanDelta(j) is the same as uzerodelta except that it stores the difference between independent and dependent
	!		chemical potentials from the tangent plane
	!      note that storate goes from 2 to numPhCo(K) so that the first one is a dummy (actually 0.0)
	! atomNorms(K) is the factor used to convert all energies to conservative atom units

!	Arrays to store working assemblages
!	AsmCurrent(K)
!		-number is kcurrent
!	asmMolesCurrent(K) - saves the moles for the current assemblage

