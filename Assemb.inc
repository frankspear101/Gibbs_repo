! 	parameters to fix maximum array boundaries
	integer*4 sysCoMax			! maximum system components in file
	integer*4 phCoMax			! maximum phase components in problem
	integer*4 phMax				! maximum phases in problem 
	integer*4 numSiteAtomMax
	integer*4 newVarMax			! maximum new variables in a phase
	parameter (sysCoMax = 30, phCoMax = 20, phMax = 30,numSiteAtomMax = 20, newVarMax = 10)
!	integer maxPC; parameter(maxPC=150)		! maximum total phase components in a problem 
	integer maxPC; parameter(maxPC=100)		! maximum total phase components in a problem 
	! maxPC must be the same value as maxPCMonit in Monit.inc -- if you change it here, you must change it there 
!********/ASSEMB/***************************************************

	character*255 AssembTitle
!     phase arrays (each entry corresponds to one phase in the problem)
!	integer*4 NREC(2000)
      	integer*4 numPh,DatasetKey(PhMax)   		! 1 = SPaC, 2 = HoPo98,   3 = HoPo11
	integer*4 HP11ModelSwitch(phMax,PhCoMax)
	CHARACTER phName(phMax)*32
	CHARACTER phCoName(phMax,phCoMax)*16
	integer*4 MINREC(phMax),numPhCo(phMax),FRACTL(phMax),PhaseType(phMax)
	REAL*8 SITMUL(phMax),MPSTART(phMax),MP0(phMax),MP1(phMax),MP2(phMax),MPLast(phMax),MPNewP(phMax)
        REAL*8 VPSTART(phMax),VP0(phMax),VP1(phMax),VP2(phMax),VMOL(phMax),SMOL(phMax),mode(phMax)
        REAL*8 PhComp(phMax,sysCoMax),PhWtOx(phMax,sysCoMax),PhWtEl(phMax,sysCoMax),gPhase(phMax)
!     phase component arrays (each entry corresponds to one phase component)
        REAL*8 xPhCo(phMax,phCoMax),xPhCoInitial(phMax,phCoMax),xPhCoLast(phMax,phCoMax),xPhCoNewP(phMax,phCoMax)
        REAL*8 hPhCoZero(phMax,phCoMax),sPhCoZero(phMax,phCoMax),vPhCoZero(phMax,phCoMax)
        REAL*8 aPhCoCp(phMax,phCoMax),bPhCoCp(phMax,phCoMax),cPhCoCp(phMax,phCoMax),dPhCoCp(phMax,phCoMax),			&
        	ePhCoCp(phMax,phCoMax),fPhCoCp(phMax,phCoMax),gPhCoCp(phMax,phCoMax), 						&
         	v1PhCo(phMax,phCoMax),v2PhCo(phMax,phCoMax),v3PhCo(phMax,phCoMax),v4PhCo(phMax,phCoMax),			&
     	    	GatTP(phMax,phCoMax),HatTP(phMax,phCoMax),SatTP(phMax,phCoMax),VatTP(phMax,phCoMax),lnAct(phMax,phCoMax),	&
     		COMP(phMax,phCoMax,sysCoMax),												&
     		int_VdP(phMax,phCoMax),int_Valpha_dP(phMax,phCoMax),								&
     		LandauTc(phMax,phCoMax),LandauSmax(phMax,phCoMax),LandauVmax(phMax,phCoMax),					&
		HP11Model(phMax,phCoMax,6),											&
     		atomNorm(phMax,phCoMax)
!     system component arrays (each entry corresponds to one system component)
	REAL*8 molWt(sysCoMax),oxToElWt(sysCoMax),numCatInOxide(sysCoMax),moles(sysCoMax),molesStart(sysCoMax),wtPct(sysCoMax)
	real*8 bulkCompWt(sysCoMax),bulkCompMoles(sysCoMax),bulkCompWtSum,openFlux(sysCoMax),					&
              bulkCompWtStart(sysCoMax),bulkCompMolesStart(sysCoMax)
 	real*8 bulkCompNorm
	logical*4 bulkCompSwitch
	logical*4 grainBoundaryProblem
	integer*4 grainBoundaryPhase
	Real*8 HseGB,HdeGB        !grain boundary energies (strain and disorder)
	INTEGER*4 noEl(sysCoMax),iOpenA(sysCoMax)
	CHARACTER coName(sysCoMax)*8,BCcoName(sysCoMax)*8,bulkCompTitle*64
!	Current assemblage variables
	integer*4 asmCurrent(phMax),asmCurrentLast(phMax),asmCurrentNewP(phMax),numPhMIF
!     misc assemblage variables
	REAL*8 TSTART,PSTART,PFSTART,TC,TK,PB,PFLUID,RockDensity,FluidDensity
	INTEGER*4 NC,NRX,NX,TANDP,NVAR,NEQ,IMASS,IFRACT,ICONS,iOpen,KFlu,PFluidSwitch
	integer*4 np,npMIF			! total number of phase components
	integer*4 gibbsVariance,duhemsVariance
	common/ASSEMB/datasetKey,assembTitle,minRec,numPhCo,fractl,					&
		sitMul,PhaseType,mpStart,mp0,mp1,mp2,vpStart,vp0,vp1,vp2,vMol,sMol,gPhase,mode,		&
		phName,xPhCo,xPhCoInitial,xPhCoLast,xPhCoNewP,MPLast,MPNewP,			&
		hPhCoZero,SPhCoZERO,VPhCoZERO,APhCoCP,BPhCoCP,CPhCoCP,DPhCoCP,EPhCoCP,FPhCoCP,GPhCoCP,	&
		v1PhCo,v2PhCo,v3PhCo,v4PhCo,GatTP,HatTP,SatTP,VatTP,int_VdP,int_Valpha_dP,		&
		LandauTc,LandauSmax,LandauVmax,atomNorm,						&
		HP11ModelSwitch,HP11Model,								&
		lnAct,COMP,phCoName,									&
		molWt,oxToElWt,numCatInOxide,moles,molesStart,wtPct,noel,				&
		iOpenA,openFlux,CONAME,BCcoName,PhComp,PhWtOx,PhWtEl,					&
		bulkCompWt,bulkCompMoles,bulkCompWtSum,bulkCompTitle,bulkCompWtStart,bulkCompMolesStart,&
		asmCurrent,asmCurrentLast,asmCurrentNewP,numPhMIF,					&
		TSTART,PSTART,PFSTART,TC,TK,PB,PFLUID,RockDensity,FluidDensity,				&
		NC,numPh,NP,npMIF,NRX,NX,TANDP,NVAR,NEQ,IMASS,IFRACT,ICONS,				&
		iOpen,KFlu,PFluidSwitch,bulkCompSwitch,bulkCompNorm,					&
		gibbsVariance,duhemsVariance,grainBoundaryProblem,grainBoundaryPhase,HseGB,HdeGB
!********************************************************************
!************/AtomSites/*********************************************
!	Common block to do ideal activity and excess calculations (Margules and Reciprocal -- SPaC only)
!	Ideal activity (config entropy) variables
!	Ideal activity part is made from expression
!	R*TK*(ActivityConstant + alpha1*ln(X1) + alpha2*ln(X2) + alpha3*ln(X3) + .....)
	Integer*4 numSiteAtom(phMax)
	real*8 siteAtomToX(phMax,phCoMax,numSiteAtomMax),xToSiteAtom(phMax,numSiteAtomMax,phCoMax),&
     		siteAtom(phMax,phCoMax),SiteMultiplicity(phMax,phCoMax)
	Character SiteAtomName(phMax,phCoMax)*8
	real*8 	Alpha(phMax,phCoMax,numSiteAtomMax),					&
		ActivityConstant(phMax,phCoMax)					 

!	Variables for reciprocal terms (SPaC only)
	integer*4 numRecipTermsMax		! maximum number of Non-ideal terms in a phase component
	parameter (numRecipTermsMax = 12)			! maximum number of reciprocal terms in a phase
	Integer*4 numReciprocal(phMax),							&
		  RecipConstant(phMax,phCoMax,numRecipTermsMax),			&
		  RecipIndex(phMax,PhCoMax,numRecipTermsMax,numSiteAtomMax)
	real*8 	WHrecip(phMax,numRecipTermsMax),					&
		WSrecip(phMax,numRecipTermsMax),					&
		WVrecip(phMax,numRecipTermsMax)

!	Variables for Margules excess terms (Berman model - SPaC only)
!	Nonideal part is made from a generalized expression of Berman (1990)
	integer*4 numMargulesSitesMax			! maximum number of sites with Margules terms
	parameter (numMargulesSitesMax = 4)
	integer*4 numMargulesSiteCatsMax		! maximum number of cations in each Margules site
	parameter (numMargulesSiteCatsMax = 9)
	integer*4 numMargulesWtermsMax			! maximum number of Margules terms in each Margules site
	parameter (numMargulesWtermsMax = 20)
	integer*4 numMargulesSites(phMax),numMargulesSiteCats(phMax,numMargulesSitesMax),		&
		  MargulesSiteCats(phMax,numMargulesSitesMax,numMargulesSiteCatsMax),			&
		  numMargulesWterms(phMax,numMargulesSitesMax),						&
		  MargulesWindex(phMax,numMargulesSitesMax,numMargulesWtermsMax)
	real*8  MargulesWH(phMax,numMargulesSitesMax,numMargulesWtermsMax),				&
		MargulesWS(phMax,numMargulesSitesMax,numMargulesWtermsMax),				&
		MargulesWV(phMax,numMargulesSitesMax,numMargulesWtermsMax)


      common /AtomSites/numSiteAtom,SiteAtomName,SiteAtomToX,XToSiteAtom,SiteAtom,SiteMultiplicity,		&
		Alpha,ActivityConstant,										&
		numReciprocal,recipConstant,RecipIndex,WHrecip,WSrecip,WVrecip,				&
		numMargulesSites,numMargulesSiteCats,MargulesSiteCats,						&
		numMargulesWterms,MargulesWindex,MargulesWH,MargulesWS,MargulesWV
		
!********************************************************************
!*******/SymmFormalism/*****************************************************
	! Variables specific to symmetric formalism margules modeling
	! ASF = size parameter
	! Note: Maximum SF terms is n*(n-1)/2. PhCoMax = 15, so this would be 14*15/2=105
	!		But I am setting max = 50 anyway
	Integer*4 SFID(phMax,50,2)				! IDs for each W term -- (Phase,Wterm,component 1, component 2)
	Integer*4 numSFterms(phMax)				! number SF terms in each phase
	real*8 SFWH(phMax,50),SFWS(phMax,50),SFWV(phMax,50),ASF(phMax,phCoMax,3)
	integer*4 DQFswitch(phMax,50)
	real*8 DQFH(phMax,50),DQFS(phMax,50),DQFV(phMax,50)
	common/SymmFormalism/SFWH,SFWS,SFWV,ASF,SFID,numSFterms,DQFswitch,DQFH,DQFS,DQFV
!********************************************************************
!*******/NEWVAR/*****************************************************
!     numNew   = Number of new variables in use
!     Knew   = Phase number for the jth new variable
!     numNewInPhK = number of new variables in phase K
!     NewPoint = Pointer to phase component array for start of Knewth phase
!     newVarName  = Name of new variable
!     NewNumA, NewNumaj = constant and coefficients for numerator of jth new variable
!     NewDenB, NewDenbj = constant and coefficients for denominator of jth new variable
!	Dimensions are (phMax,newVarMax) - this allows every phase to have a maximum of newVarMax new variables
!		Each new variable can have a maximum of phCoMax entries (in fact, it MUST have exactly numPhCo(k) for the phase in question
!     numerator, denominator,denominator are temporary storage used in C2pute      

      INTEGER numNew,KNEW(30),numNewInPhK(phMax),NewPoint(30),newVarGlobal
!      character NewVarTitle(30)*32,NewName(30)*16
      character newVarName(phMax,newVarMax)*32
!      real*8 NewNumA(30),NewNumaj(30,10),NewDenB(30),NewDenbj(30,10)
      real*8 NewNumA(phMax,newVarMax),NewNumaj(phMax,newVarMax,phCoMax),NewDenB(phMax,newVarMax),NewDenbj(phMax,newVarMax,phCoMax)
      real*8 numerator,denominator,denominator2
      COMMON /NEWVAR/numNew,Knew,numNewInPhK,newPoint,newVarName,			&
		NewNumA,NewNumaj,NewDenB,NewDenbj,					&
     		numerator,denominator,denominator2,newVarGlobal
!********************************************************************
!******/Fixed_Min_Data/****************************************************
!     Common block to store "fixed" mineral data for kyanite, sillimanite, andalusite, alpha and beta quartz
!     The subscripts for the storage are set at the time the thermodynamic data file is read in
!	---previously these were alpha quartz = 1, beta quartz = 2, Kyanite = 3, sillimanite = 4, andalusite = 5
!     phase component arrays (each entry corresponds to one phase component)
	integer*4 fnMax,fcMax
	parameter (fnMax = 4,fcMax=12)
	integer*4 HP11ModelSwitchF(fcMax),DQFswitchF(fcMax)
! 	fixedNo = number of phases that will have "fixed" status
!	fixedID = phase number of the fixed phase (use in input file to select this phase)
!	fixedIDNoCo = number of "components" in this fixed phase
! 	fixedPhase = numbers from the data file that are the "component" choices of this "fixed" phase
!	e.g. KySilAnd would have 5,6 207 corresponding to Ky, Sil and And in the data file
	integer*4 fixedNo,fixedID(fnMax),fixedIDNoCo(fnMax),fixedPhase(fcMax)
      REAL*8 HZEROF(fcMax),SZEROF(fcMax),VZEROF(fcMax),						&
          ACPF(fcMax),BCPF(fcMax),CCPF(fcMax),DCPF(fcMax),ECPF(fcMax),FCPF(fcMax),GCPF(fcMax),	&
          v1F(fcMax),v2F(fcMax),v3F(fcMax),v4F(fcMax),						&
          LandauTcF(fcMax),LandauSmaxF(fcMax),LandauVmaxF(fcMax),				&
	  HP11ModelF(fcMax,6),									&
	  DQFHF(fcMax),DQFSF(fcMax),DQFVF(fcMax)      
      common/Min_Data/										&
		HZEROF,SZEROF,VZEROF,ACPF,BCPF,CCPF,DCPF,ECPF,FCPF,GCPF,v1F,v2F,v3F,v4F,	&
		landauTcF,landauSmaxF,landauVmaxF,						&
		fixedNo,fixedID,fixedIDNoCo,fixedPhase,						&
		HP11ModelSwitchF,HP11ModelF,DQFswitchF,DQFHF,DQFSF,DQFVF
!********************************************************************
!*******/IOstuff/*****************************************************
!      	storage for Thermodynamic datafile information
!     	and for opening and saving problems
!     	NoSysCoInFile = number of system components in the file
!     	CoNameInFile = Names of system components in the data file
!     	UseSysCo = 0 or 1 depending on whether a system component is used
!     	NoMinInFile = total number of minerals in the data file
!     	MineralNumber = the number of the mineral in the datafile
!     	MineralPosition = line number in the data file on which the mineral starts
!     	MineralName = name of mineral
!     	numPhCoInFile = the number of phase components for a mineral in the data file (used or not)
!     	phCoNameinFile = Names of phase components in the data file
!     	UsePhCo = 0 or 1 depending on whether a phase component from data file is used.
!     	IncludeNew = 0 or 1 depending of whether the "new" variables are to be included for a phase
      	integer*4 NoSysCoInFile,UseSysCo(sysCoMax),NoMinInFile
	integer*2 MineralNumber(512),MineralPosition(3000)
      	Integer*4 numPhCoInFile(phMax),phCoUsedinFile(phMax,phCoMax),UsePhCo(phMax,phCoMax),	&
     		IncludeNew(phMax),MinInFilePoint(phMax),numNewInFile(phMax)
      	character mineralname*32,coNameInFile(sysCoMax)*8
      	CHARACTER phCoNameinFile(phMax,phCoMax)*16
      	real*4 compPlotCoeffInFile(20,sysCoMax)
      	REAL*8 MOLWTinFile(sysCoMax),OxToElWtInFile(sysCoMax),NumCatInOxideInFile(sysCoMax)
      	common /IOstuff/ NoSysCoInFile,UseSysCo,NoMinInFile,MineralNumber,			&
     		MineralPosition,MineralName,numPhCoInFile,phCoUsedinFile,UsePhCo,		&
     		IncludeNew,CoNameInFile,phCoNameinFile,MinInFilePoint,				&
     		numNewInFile,CompPlotCoeffInFile,MolWtinFile,OxToElWtInFile,NumCatInOxideInFile
!********************************************************************
!*******/ThermoRx/*****************************************************
!     	Dimension (phCoMax) is the max number of phase components in a single phase + 2 (for T and P)
!     	This variable contains the partial derivatives of each chemical w/r/t
!     	i/T = dudTPX(i,1)
!     	i/P = dudTPX(i,2)
!     	i/Xj= dudTPX(i,3 to numPhCo(K)+1)
!     	DIMENSION (maxPC,maxPC) REFERS TO THE MAXIMUM TOTAL PHASE COMPONENTS IN A problem
      	REAL*8 dudTPX(phMax,phCoMax,phCoMax)
      	REAL*8 ARX(maxPC,maxPC)
      	COMMON /ThermoRx/ dudTPX,ARX
!********************************************************************
