! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GLOBAL(jopt)
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "Tangent.inc"
! *****************************************


! ********************************************************************
!      local variables
      INTEGER iglobe,imiss,L,i,ii,j,jj,k,iopt3,izero,iok,PFluidSwitchOld,newPtrOld,iOpenOld,jopt,kCur
	character*256 TitleLine
	integer*4 phNumberMIFTT,numPhCoMIFTT,lsite
	character*8 phNameMIFTT,phCoNameMIFTT(PhMax)
	real*8 PhGDiffTT,PhMolesTT,PhVolTT,XPhCoTemp(PhMax)
	real*4 TCMIF,PBMIF
! ********************************************************************

	if(jopt.lt.0)then
		iglobe = abs(jopt)
		go to 2
	endif

1     CONTINUE
	if(jopt.lt.0)return			! exit directly if jopt is negative
      iglobe=0                			! default
      WRITE(*,*)' ******************************'
      WRITE(*,*)' Change GLOBAL variables menu'
      WRITE(*,*)'   0 = Return '
      WRITE(*,*)'   1 = Change mineral assemblage'
      WRITE(*,*)'   2 = Set output length '
      WRITE(*,*)'   3 = Add or drop mass balance'
      WRITE(*,*)'   4 = Display reaction matrix'
      WRITE(*,*)'   5 = Set pfluid'
!      WRITE(*,*)'   6 = Add new variables to list'
      WRITE(*,*)'   7 = Change P & T without changing anything else'
      WRITE(*,*)'   8 = OPEN/CLOSED System behavior'
      WRITE(*,*)'   9 = Change Modal proportions of a phase'
      WRITE(*,*)'  10 = List ALLX variable array for debugging'
      WRITE(*,*)'  11 = Create a MIF from a part of a .ALL file '
      WRITE(*,*)' INPUT CHOICE'
	read(*,*)iglobe


	if(iglobe.eq.0)return

2	continue	! get here directly if jopt < 0

! -------------------------------------------------------
!      CHANGE ASSEMBLAGE
! -------------------------------------------------------
	if(iglobe.eq.1)then
! 51    CONTINUE
!      		if we are changing the assemblage, we must zero out any "new" variables
!      		because there is no way to determine whether they will still be in the 
!      		list after a mineral is removed.
      		NVAR=NVAR-numNew
      		NEQ=NEQ-numNew
      		numNew=0
       		call change
!       		call rexn
		call PrintMinAssemblage
       		goto 1
	endif
! -------------------------------------------------------

	if(iglobe.eq.2)then
! -------------------------------------------------------
!      OUTPUT LENGTH
! -------------------------------------------------------
	call AWE_setoutput
       goto 1

	endif
	
	if(iglobe.eq.3)then
! -------------------------------------------------------
!      IMASS
! -------------------------------------------------------
53    CONTINUE
       WRITE(*,*)' OPTION TO ADD OR REMOVE MASS BALANCE CONSTRAINTS'
1800   WRITE(*,*)' Current value of IMASS=',imass
       WRITE(*,*)' 0 = RETURN'
       WRITE(*,*)' 1 = Drop mass balance constraints'
       WRITE(*,*)' 2 = Add mass balance constraints'
	write(*,*)'3 = Turn off BulkCompositionSwitch (negates bulk rock analysis)'
       READ(*,*)imiss
       if(imiss.eq.0)go to 1
       if(imass.eq.1.and.imiss.eq.1)then
          imass=0
          nvar=nvar-numPh
          neq=neq-nc
          call names
          go to 1
       endif
       if(imass.eq.0.and.imiss.eq.2)then
          imass=1
          nvar=nvar+numPh
          neq=neq+nc
          call names
          go to 1
       endif
	if(imiss.eq.3)then
		bulkCompSwitch = .false.			! a bulk composition is now negated
		go to 1
	endif
       WRITE(*,*)'You cannot do that--- please try again'
       go to 1800

	endif

	if(iglobe.eq.4)then
! -------------------------------------------------------
!      REACTION MATRIX
! -------------------------------------------------------
54    CONTINUE
!          RECOMPUTE REACTION MATRIX
	  iLong(7) = 1
          CALL REXN
	  iLong(7) = 0
        GO TO 1

	endif


	if(iglobe.eq.5)then
! -------------------------------------------------------
!      PFLUID
! -------------------------------------------------------
56    CONTINUE

!      SET UP FOR PFLUID
	call FluidPresent		! returns KFlu>0 if a fluid is present
	if(KFlu.gt.0)then		! ask user how to deal with fluid pressure
		PFluidSwitchold = PFluidSwitch
		call SetPFluid(iok)  ! routine to set P fluid variables 
		if(iok.eq.1)go to 1
		if(PFluidSwitchold.eq.2)NEQ=NEQ-1
		IF(PFluidSwitch.gt.0)then
			PFStart = PFluid
			TandP = 3
		        NVAR = TANDP+NX+numNew
		        IF(IMASS.EQ.1)NVAR=NVAR+numPh
			if(PFluidSwitch.eq.2)NEQ=NEQ+1
			if(PFluidSWitchOld.eq.0)then
				! Adjust ALLX array to make room for PFluid as row 3
				! Only do this if previously there was no PFluid
	           		DO 444 L=NVAR,4,-1
	            		Do 444 J=1,6
				ALLX(J,L)=ALLX(J,L-1)
444     	   		CONTINUE
           			endif
			! Set new value for PFluidStart
			Do 445 J=1,6
445     		ALLX(J,3)=PFSTART
			else	! Pfluid is lithostatic
			TandP=2
		        NVAR = TANDP+NX+numNew
		        IF(IMASS.EQ.1)NVAR=NVAR+numPh
			if(PFluidSwitchOld.gt.0)then
				! Close up AllX array to get rid of PFluid
				DO 446 L=3,NVAR
        			ALLX(1,L)=ALLX(1,L+1)
446     			CONTINUE
				endif
			endif

		else	! no fluid is present
		call fss_alert('ALERT!!','Sorry, but there is no fluid present')
		endif
        CALL NAMES
	GO TO 1
	endif



	if(iglobe.eq.6)then
		go to 1
		endif
	
	
	
	if(iglobe.eq.7)then
! -------------------------------------------------------
!      P AND T
! -------------------------------------------------------
58    CONTINUE
        WRITE(*,*)' Current values of TC, Prock are:'
        WRITE(*,*)TC,PB,PFluid
        WRITE(*,*)'Give new value for T(degrees C)'
        READ(*,*)TC
        ALLX(1,1)=TC
        ALLX(2,1)=TC
        WRITE(*,*)'Give new value for Prock(bars C)'
        READ(*,*)Pb
        ALLX(1,2)=PB
        ALLX(2,2)=PB
        GO TO 1

	endif
	
	
	if(iglobe.eq.8)then
! -------------------------------------------------------
!      Setup open system behavior
! -------------------------------------------------------
700   	continue
!      	NVAR = NVAR - IOPEN
	iOpenOld = iOpen
	newPtrOld = newVarPtr
710   	write(*,*)' '
      	write(*,*)' Current status of open system is:'
      	write(*,*)' IOPEN = ',iopen
      	write(*,704)(CONAME(iopena(i)),i=1,iopen)
704   	format (' ',3x,15(a4,1x))
      	write(*,705)(IOPENA(i),i=1,iopen)
705   	format(' ',15i5)
      	write(*,*)' '
      	write(*,*)' OPTIONS:'
      	write(*,*)' 0 = RETURN to GLOBAL menu'
      	write(*,*)' 1 = MODIFY open system components'
      	read(*,*)i
      	if(i.eq.0)then
!            NVAR = NVAR + IOPEN
!            call names
            go to 1
            endif
      	write(*,*)' Input number of open system components (iOpen)'
      	read(*,*)iOpen
	if(iOpen.eq.0)go to 715
      	WRITE(*,*)' Input integers corresponding to system components from list below:'
      	write(*,704)(CONAME(i),i=1,NC)
      	write(*,705)(i,i=1,NC)
      	read(*,*)(iOpena(i),i=1,iopen)

715	continue
        nVar = nVar + iOpen - iOpenOld
	call Names
! 	Note: a pointer to the first open system component in allX array (openPtr in common/monit) is set up in subroutine names
	jj = newVarPtr - newPtrOld
	if(jj.ne.0)then				! we must adjust allX array
		if(jj.lt.0)then			! move allX up
			do 721 i = 1,numNew
			DO 720 ii=1,5
     			allX(ii,newVarPtr + i - 1) = allX(ii,newPtrOld + i - 1)
720			continue          
721			continue
		else				! move allX down	
			do 726 i = numNew,1,-1
			DO 725 ii=1,5
     			allX(ii,newVarPtr + i - 1) = allX(ii,newPtrOld + i - 1)
725			continue          
726			continue
		endif

		! zero out open system in allX
		do 730 i = 1,iOpen
		do 731 ii = 1,5
		allX(ii,openPtr + i - 1) = 0.0d0
731		continue
730		continue

	endif

	iOpenOld = iOpen
	newPtrOld = newVarPtr

      	go to 710

	endif

	if(iglobe.eq.9)then
! -------------------------------------------------------
!      Change modal propostions of phases
! -------------------------------------------------------
1300  continue
      WRITE(*,1330)(MINREC(asmCurrent(I)),PHNAME(asmCurrent(I)),I=1,numPh)
1330  FORMAT(3(I4,2X,A8,2x))
      WRITE(*,*)' Change modal proportions of a phase'
      WRITE(*,*)' Input number of phase to adjust (0 to exit) '
      READ(*,*)iopt3
      if(iopt3.eq.0)go to 1
!      DO 1310 k = 1,numPh
	Do 1310 kCur = 1,numPh
	k = asmCurrent(kCur)
1310  IF(MINREC(K).EQ.IOPT3)GO TO 1350
!      MINERAL IOPT3 IS NOT IN MINERAL FILE
      WRITE(*,*)' MINERAL NUMBER',IOPT3,' IS NOT IN LIST. TRY AGAIN'
      GO TO 1300
1350  CONTINUE
!      THE MINERAL IS FOUND and its number is K
      WRITE(*,1330)MINREC(K),PHNAME(K)
      izero=0
      call AllKduTPX(1,izero)     ! calculate molar volume of phase K

      VP0(K)=MP0(K)*VMOL(K)
      write(*,*)' Molar volume (j/bar-mol) = ',VMOL(K)
      write(*,*)' Current values for this phase are:'
      write(*,*)'   Moles       Volume (cm**3)'
      write(*,*)MP0(K),VP0(K)*10.D0
      
      WRITE(*,*)'Input new VOLUME for this phase (in cm**3)'
      READ(*,*)VP0(K)
      VP0(K)=VP0(K)/10.D0           ! convert to j/bar
      MP0(K)=VP0(K)/VMOL(K)
      write(*,*)' New values for this phase are:'
      write(*,*)'   Moles       Volume (cm**3)'
      write(*,*)MP0(K),VP0(K)*10.D0
      go to 1300


	endif
	
	
	if(iglobe.eq.10)then
! -----------------------------------------------------------------
!      List ALLX array for debugging
! -----------------------------------------------------------------
1400  continue
!      Definition of ALLX variable:
!      1,50   current value
!      2,50     starting value
!      3,50   ref at start of contour
!      4,50   user selected reference
!      5,50   previous finite diff point
! 
      WRITE(12,*)'   '
      WRITE(12,*)'***** Check variable arrays ********'
      WRITE(12,*)'Name, ALLX(1),ALLX(2),ALLX(3),ALLX(4),ALLX(5),ALLX(6)'
      WRITE(12,*)'  Name                Current     Starting    Contour  User ref    previous fdiff    FRACTL(Y/N)'
      do 1441 I=1,nvar
1441    WRITE(12,1442)I,VN1(I),VN2(I),(ALLX(j,i),J=1,5),allFractl(i)
1442    FORMAT(I3,1X,A2,A16,5F12.4,I5)
!      write(*,*) 'Hit return to continue'
!      pause
      go to 1
	endif

! -----------------------------------------------------------------
!      Translate mole fractions into site occupancies
! -----------------------------------------------------------------
	if(iglobe.eq.11)then
		write(*,*)' This routine requires you read in a single entry from a .ALL file (from the MAD routine)'
		write(*,*)' You need to create a new .ALL file with ONLY the P and T that you want'
		write(*,*)' Do NOT include the header. Only include the Line with the T and P, the column header, and the phases'
		write(*,*)' Here''s how it should look'
		write(*,*)'     500.0    2000.0       -779347.62470    7 = numPhases'
		write(*,*)'                  Phase           GDiff         mMoles phase           V phase              Phase composition'
		write(*,*)'            2001  Quartz-a        0.00000        0.76190E+03        0.17626E+02        1     Quartz      0.10000E+01'
		write(*,*)'            2003  Fluid           0.00000        0.15686E+03        0.40627E+01        1     H2O         0.10000E+01'
		write(*,*)'            2017  Muscovit        0.00000        0.59704E+02        0.84244E+01        4     &
			&Celadoni    0.14291E-01     Fe-celad    0.13768E-01     Muscovit    0.77408E+00     Paragoni    0.19786E+00'
		write(*,*)' etc.'
		write(*,*)' '		
		write(*,*)' The .ALL file MUST be in the same order as the MIF file that created it'
		write(*,*)' AND you MUST have already opened this MIF file'
		write(*,*)' Have you already opened the MIF file and are ready to open the modified .ALL file?'
		write(*,*)' 0 = no (abort); 1 = yes (continue and open modified .ALL file)'
		read(*,*)iok
		if(iok.eq.0)go to 1
		call FSS_Alert('ALERT','Open modified .ALL file for input')
		open(42,file='',status='OLD',iostat=iok)
		if(iok.ne.0)go to 1 	! user hit cancel
!		Here we parse the modified .ALL file
		read(42,*)TCMIF,PBMIF
		write(*,*)TCMIF,PBMIF
		read(42,1001)TitleLine
		write(*,1001)TitleLine
1001		format(A256)
		numCurrent = 0
		do 1105 k = 1,numPhMIF
		read(42,*,end=1115)phNumberMIFTT,phNameMIFTT,phGDiffTT,phMolesTT,phVolTT,numPhCoMIFTT,   &
	                 (phCoNameMIFTT(j),xPhCo(k,j),j=1,numPhCoMIFTT)
		if(phGDiffTT.lt.1e-2)then
			numCurrent = numCurrent + 1
			asmCurrent(numCurrent) = k
			endif
		write(*,1104)phNumberMIFTT,phNameMIFTT
1104		format(I5,5x,A8)
		write(*,1112)(phCoNameMIFTT(j),xPhCo(k,j),j=1,numPhCoMIFTT)
1112		format(20(A8,E15.5,5x))
		if(numSiteAtom(k).gt.0)then 
	      		call XtoSite(K)
      			write(*,1101)(SiteAtomName(k,Lsite),Lsite = 1,numSiteAtom(K))
1101   			format(2x,20A13)
      			write(*,1102)(SiteAtom(k,Lsite),Lsite = 1,numSiteAtom(K))
1102   			format (20E13.6)
			endif
1105		continue	
1115		continue
	call SAVEMIF(TCMIF,PBMIF)

	endif
! --------------------------------------------
      go to 1

      END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE CHANGE
      implicit none
!      ROUTINE TO ALLOW THE MINERAL ASSEMBLAGE TO BE CHANGED

! *****************************************
	include "Assemb.inc"
! *****************************************
! ********************************************************************
!      LOCAL VARIABLES
      	INTEGER*4 iopt3,k,kCur,iAdd,iChanged
	integer*4 minRecTemp(100)
	character*8 phNameTemp(100)	

! ---------------------------------------------------------------
	iChanged = 0		! we haven't changed anything yet
10    CONTINUE

      	WRITE(*,*)' *****CHANGE MINERAL ASSEMBLAGE:*****'
      	write(*,*)' Current mineral assemblage'
      	WRITE(*,133)(MINREC(asmCurrent(K)),PHNAME(asmCurrent(K)),k = 1,numPh)
133   	FORMAT(3(I4,2X,A8,2x))
!	write out minerals not in the assemblage
	iAdd = 0
	do 12 k = 1,numPhMIF
	do 13 kcur = 1,numPh
	if(minRec(k).eq.minRec(asmCurrent(kcur)))go to 12
13	continue
	iAdd = iAdd + 1
	minRecTemp(iAdd) = minRec(k)
	phNameTemp(iAdd) = phName(k)
12	continue
	write(*,*)'Other minerals'
	write(*,133)(minRecTemp(k),phNameTemp(k),k = 1,iAdd)

!      	NOTE  THIS CODE DOES NOT ACCOUNT FOR numNew (EXTRA) VARIABLES AND EQUATIONS

      	iopt3=0           !default
      	write(*,*)
      	WRITE(*,*)' OPTIONS:'
      	WRITE(*,*)' 0 = All done - Return'
      	WRITE(*,*)'+n = ADD mineral to list (note: minerals are added to the end of the list)'
      	WRITE(*,*)'-n = REMOVE mineral number N from the list '
      	WRITE(*,*)'     note: "n" refers to the number in the thermodynamic data file listed above'
      	READ(*,*)iopt3
      	if(iopt3.eq.0)then
!		Calculate total number of phase components
		NP=0
	!	newVarGlobal = 0			! flag that says we have new variables to consider (0 = no, 1 = yes)
      		do 15 kCur = 1,numPh
		k = asmCurrent(kCur)
		NP = numPhCo(K)+NP
		if(iChanged.ne.0)mp0(k) = 0.1		! reset the moles of phases to .1  (100 millimoles) but only if we actually changed the assemblage
!		if(includeNew(k).eq.1)newVarGlobal = 1		! this flag is used in subrotine setAllX
15		CONTINUE
!     		compute number of independent compositional (dX) variables (before SetAllX)
	!	don't reset the compositions .... it's a pain
	!	xPhCo = xPhCoInitial 	! set all compositions back to initial (input) values
      		NX=NP-numPh
		call SetALLX()		! set up ALLX array
!		numNew is calculated in SetALLX
		call Names()		! sets up names of variables for asmCurrent
      		call REXN		! calculate linearly independent reactions
!     		compute total number of variables
      		NVAR= TANDP + NX
      		if(IMASS.EQ.1)NVAR=NVAR + numPh
		if(numNew.gt.0)nvar = nvar + numNew
	      	NEQ = NRX
      		IF(IMASS.EQ.1)NEQ = NEQ + NC
      		if(numNew.gt.0)NEQ = NEQ + numNew
      		if(PFluidSwitch.eq.2)NEQ = NEQ + 1
      		return
		endif
! ------------------------------------------
      IF(IOPT3.gt.0)then
!		write(*,133)(Minrec(k),PhName(k),k=1,numPhMIF)
!		write(*,*)'Input number of phase to add -- 0 to abort'
!		read(*,*)iAdd
		iAdd = iopt3
		do 20 k = 1,numPhMIF
		if(iAdd.eq.Minrec(k))then
			numPh = numPh+1
			asmCurrent(numPh) = k
			iChanged = 1		! we changed the assemblage
			go to 10
			endif
20		continue
		write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	      	WRITE(*,*)' Mineral number',IOPT3,' is not in the list.  Try again'
		write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	      	go to 10
      		endif
! --------------------------------------------------------
!      REMOVE A MINERAL
      if(iopt3.lt.0)then
!      Check to see if this mineral is really in data file
	Do 1010 kCur = 1,numPh
	k = asmCurrent(kCur)
1010    IF(MINREC(K).EQ.-IOPT3)GO TO 1020
!      	MINERAL IOPT3 IS NOT IN MINERAL FILE
		write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      		WRITE(*,*)' Mineral number',-IOPT3,' is not in the list.  Try again'
		write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      	GO TO 10
1020  	CONTINUE
	k = 0
	do 1030 kCur = 1,numPh
        IF(MINREC(asmCurrent(kCur)).NE.-IOPT3)then
		k = k + 1
		asmCurrent(k) = asmCurrent(kCur)
		endif
1030	continue
	iChanged = 1		! we changed the assemblage
	numPh = numPh - 1
      	endif
! --------------------------------------------
      GO TO 10
      END



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE ThermoData
	use AWE_Interfaces
	use MyCanvas
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
! *****************************************
! ********************************************************************
!      	local variables
      	INTEGER ioption,L,i,j,jj,jjj,k,iopt,iopt3,izero,iok,ivol,j1,j2,myColor
      	real*8 TCsave,Pbsave,TKsave,Pfluidsave,Get8,pbPlot
      	real*4 Ttemp,Ptemp
!      	character ppp*128
      	character*64 fileName

	integer*4 iup,yplot,xplot
	real*4 yyy,xxx,xinc
      	real*8 xTemp(20),uMin(20),R,Hsetemp,Hdetemp

 	real*8 pValue(PhCoMax),lnaIdeal(PhCoMax),RTlnGamma(PhCoMax)
 	common /gamma/pValue,lnaIdeal,RTlnGamma


      	DATA R/8.3144D0/


1     CONTINUE
	ioption=0                	! default
        write(*,*)' ***********************************'
        write(*,*)'Thermo file: ',THMOFILE         ! thermodynamic data file
        WRITE(*,*)' ***********************************'
	WRITE(*,*)' List thermodynamic data menu'
	WRITE(*,*)'  0 = Return '
	write(*,*)'  1 = Open new Thermodynamic datafile (in Gibbs_Essentials )'
	WRITE(*,*)'  2 = List Thermo data file info for a phase'
	WRITE(*,*)'  3 = List thermo data at T and P'
	WRITE(*,*)'  4 = Plot thermodynamic data for a phase'
	WRITE(*,*)'  5 = Examine and adjust thermo data for a phase'
	WRITE(*,*)'  6 = Solve for H of phase components'
	WRITE(*,*)'  7 = Save and reopen thermofile'
	WRITE(*,*)'  8 = Make cordierite calculations'
	WRITE(*,*)'  9 = Change Grain Boundary energies'
	WRITE(*,*)' 10 = Table of thermodynamic datafile'
	WRITE(*,*)' 11 = Least Squares calculation of Ho'
	
	WRITE(*,*)' INPUT CHOICE'
	read(*,*)ioption

	if(ioption.eq.0)return

! 	--------Solve for H of phase component ------------------
	if(ioption.eq.11)then
	call LeastSquares_H()
	go to 1
	endif
! 	----------Open new thermo file ----------------------
	if(ioption.eq.1)then

!	call volset(volrefnum)        ! set volume to startup directory


	call ReadThomoFileList(fileName)
	thmoFile = trim(fileName)
	call OpenThermoFile

! 	      close(1)
! 	      Open thermodynamic data file
! 	      fileName = ':Gibbs_Essentials:'//trim(THMOFILE)		!put path in .fig file, not here
! 	      OPEN (1,FILE=fileName,STATUS='OLD')
! 	      iok=0
! 	      call setthermofile(iok)
! 	      if(iok.eq.1)then
! 	            call fss_alert(' There is a problem reading the thermodynamic data file.')
! 	            endif

! 	     Load "fixed" minerals from thermodynamic data file
! 	     1 = alpha quartz, 2 = beta quartz, 3=kyanite,4=sill,5=andalusite
! 	      call loadfx
	
	go to 1
	endif



	if(ioption.eq.10)then
! -------------------------------------------------------
!      List thermo data file info for phases - Table format
! -------------------------------------------------------


	write(12,*)'PhCo  H S  V   A   B   C  D   E   F    G    v1 v2  v3    v4    LandauTc    LandauSmax    LandauVmax '
	do 5015 K = 1, numPhMIF


      write(12,*)
      write(12,5022)Phname(K)
5022	format(A16)
      do 5011 j=1,numPhCo(K)
      jjj=jj+j
	if(minrec(K).ge.1000)then
      write(12,5021)phCoName(k,j),hPhCoZero(k,j),sPhCoZero(k,j),vPhCoZero(k,j),&
     &  aPhCoCp(k,j),bPhCoCp(k,j),cPhCoCp(k,j),dPhCoCp(k,j),ePhCoCp(k,j),fPhCoCp(k,j),gPhCoCp(k,j),&
     &  v1PhCo(k,j),v2PhCo(k,j),v3PhCo(k,j),v4PhCo(k,j),LandauTc(k,j),LandauSmax(k,j),LandauVmax(k,j)
	else
      write(12,5021)phCoName(k,j),hPhCoZero(k,j),sPhCoZero(k,j),vPhCoZero(k,j),&
     &  aPhCoCp(k,j),bPhCoCp(k,j),cPhCoCp(k,j),dPhCoCp(k,j),ePhCoCp(k,j),fPhCoCp(k,j),gPhCoCp(k,j),&
     &  v1PhCo(k,j),v2PhCo(k,j),v3PhCo(k,j),v4PhCo(k,j)
	endif
5021	format(A8,2x,F12.1,F12.3,F12.3,F12.3,F12.3,F15.1,F15.1,3E15.5,4e15.5,3E15.5)

5011   continue
5015	continue
      GO TO 1
	endif		! end ioption = 10
	



	if(ioption.eq.2)then
! -------------------------------------------------------
!      List thermo data file info for a phase
! -------------------------------------------------------
55    CONTINUE
      write(*,*)
      WRITE(*,301)(MINREC(K),PHNAME(K),k = 1,numPhMIF)
301   FORMAT(3(I4,2X,A8,2x))

      iopt=0           		! default
      write(*,*)
      WRITE(*,*)' OPTIONS:'
      WRITE(*,*)' 0 = Return'
      WRITE(*,*)' n = Pick mineral n to list'
      READ(*,*)IOPT
      if(iopt.eq.0)go to 1
!      Check to see if this mineral is really in data file
      DO 302 k = 1,numPhMIF
302   IF(MINREC(K).EQ.IOPT)GO TO 303
!      MINERAL IOPT IS NOT IN MINERAL FILE
      WRITE(*,*)' Mineral number',-IOPT3,' is not in the list.  Try again'
      GO TO 55
303   CONTINUE
      write(12,*)
      write(12,*)Phname(K)
      write(12,*)Minrec(K),numPhCo(K),PhaseType(K),numSiteAtom(K),numReciprocal(k)
      if(numSiteAtom(k).gt.0)then
	      write(12,*)(SiteAtomName(k,L),L=1,numSiteAtom(K))
	      write(12,*)(SiteMultiplicity(k,L),L=1,numSiteAtom(K))
	      endif
      do 310 j=1,numPhCo(K)
 !     jjj=jj+j
      write(12,*)'---------------------------------------'
      write(12,*)phCoName(k,j)
      write(12,*)hPhCoZero(k,j),sPhCoZero(k,j),vPhCoZero(k,j)
      write(12,*)aPhCoCp(k,j),bPhCoCp(k,j),cPhCoCp(k,j),dPhCoCp(k,j),ePhCoCp(k,j),fPhCoCp(k,j),gPhCoCp(k,j)
      write(12,*)v1PhCo(k,j),v2PhCo(k,j),v3PhCo(k,j),v4PhCo(k,j)
	Select Case (DataSetKey(k))
		case(2)		!HP98 dataset
			write(12,*)'Landau values ',LandauTc(k,j),LandauSmax(k,j),LandauVmax(k,j)
		case(3)		!HP11 dataset
			write(12,*)'HP11 model ',HP11ModelSwitch(k,j),(HP11Model(k,j,L),L=1,6)
			write(12,*)'DQF values ',DQFswitch(k,j),DQFH(k,j),DQFS(k,j),DQFV(k,j)
		case default
		end select


      write(12,*)(CoName(i),i=1,NC)
      write(12,*)(Comp(k,j,i),i=1,NC)
      if(numSiteAtom(K).gt.0)then
      write(12,*)'ActivityConstant  ',(SiteAtomName(k,L),L=1,numSiteAtom(K))
      write(12,*)ActivityConstant(k,j),(Alpha(k,j,i),i=1,numSiteAtom(K))
!      write(12,*)numNIDterms(k,j)
!      if(numNIDterms(k,j).gt.0)then
!      write(12,*)' WH   WS  WV   Ca  NXa1  NXa2  NXa3  Cb   NXb1  NXb2    NXb3'
!      do 315 n=1,numNIDterms(k,j)
!      nnn=nn+n
!      write(12,*)NIDWH(k,j,n),NIDWS(k,j,n),NIDWV(k,j,n),numNIDX(k,j,n),			&
!     		(NIDC(k,j,n,i),NIDX(k,j,n,i,1),NIDX(k,j,n,i,2),NIDX(k,j,n,i,3),i=1,numNIDX(k,j,n))
!315   continue
!      endif
      endif
310   continue
       GO TO 55
	endif		! end ioption = 1
	

! -------------------------------------------------------
!      List thermodynamic properties at T and P 
! -------------------------------------------------------
	if(ioption.eq.3)then

      TKsave = TK
      TCsave = TC
      Pbsave = Pb
      Pfluidsave = Pfluid

503   continue


      Ttemp = TC
      Ptemp = Pb
      write(*,*)' TC and Pb = ',Ttemp,Ptemp
      write(*,*)' Specify T in deg C and P in bars: TC,Pb'
      read(*,*)(Ttemp),Ptemp
	if (Ttemp.eq.0)go to 503      


59    CONTINUE
      WRITE(*,630)(MINREC(I),PHNAME(I),I=1,numPhMIF)
      WRITE(*,*)' Options:'
      iopt = 0
      WRITE(*,*)' 0 = Return'
      WRITE(*,*)' +n = Specify mineral number'
      read(*,*)(iopt)
      IF(IOPT.EQ.0)then
        TC = TCsave
        TK = TKsave
        Pb = Pbsave
        Pfluid = Pfluidsave
        GO TO 1
        endif

      DO 510 k = 1,numPhMIF
510   IF(MINREC(K).EQ.IOPT)GO TO 511
!      MINERAL IOPT IS NOT IN MINERAL FILE
      call FSS_alert('ALERT!!','The mineral number you specified is not in the list. Try again') 
      GO TO 503
511   CONTINUE
!      THE MINERAL IS FOUND and its number is K



      TC = Ttemp
      TK = TC+273.15
      Pb = Ptemp
      Pfluid = Pb

      izero=0
      ivol=0            ! compute everything
!       K=0
!      CALL AllKduTPX(ivol,izero)
	call CalculateCPTemp(TK)
	call duTPX(k,ivol,izero)
	call dlnAdX(k,izero)

      if(izero.eq.1)go to 1


!      write out current thermodynamic data
      write(12,*)'************************************************'
      write(12,*)' Option to list thermodynamic properties at T and P'
      write(12,*)'   TC     TK       Pbars'
      write(12,*)TC,TK,Pb
!      K is the phase to write out

      write(12,*)
      write(12,*)'Minrec,sitmul,phName,special'
      write(12,520)MINREC(K),SITMUL(K),PHNAME(K),PhaseType(K)
520   FORMAT(' ',I8,F8.1,4X,A32)


      write(12,522)(phCoName(k,J),J=1,numPhCo(K))
522   FORMAT('               ',15(A8,12X))
      write(12,526)(hPhCoZero(k,J),J=1,numPhCo(K))
526   Format(' hPhCoZero',15E20.9)
      write(12,528)(sPhCoZero(k,J),J=1,numPhCo(K))
528   Format(' SZERO',15F20.5)
      write(12,530)(vPhCoZero(k,J),J=1,numPhCo(K))
530   Format(' VZERO',15F20.5)
	WRITE(12,535)' ACP',(aPhCoCp(k,J),J=1,numPhCo(K))
	WRITE(12,535)' BCP',(bPhCoCp(k,J),J=1,numPhCo(K))
	WRITE(12,535)' CCP',(cPhCoCp(k,J),J=1,numPhCo(K))
	WRITE(12,535)' DCP',(dPhCoCp(k,J),J=1,numPhCo(K))
	WRITE(12,535)' ECP',(ePhCoCp(k,J),J=1,numPhCo(K))
	WRITE(12,535)' FCP',(fPhCoCp(k,J),J=1,numPhCo(K))
	WRITE(12,535)' GCP',(gPhCoCp(k,J),J=1,numPhCo(K))
535     format(' ',a5,15E20.5)
	WRITE(12,535)'  v1',(v1PhCo(k,J),J=1,numPhCo(K))
	WRITE(12,535)'  v2',(v2PhCo(k,J),J=1,numPhCo(K))
	WRITE(12,535)'  v3',(v3PhCo(k,J),J=1,numPhCo(K))
	WRITE(12,535)'  v4',(v4PhCo(k,J),J=1,numPhCo(K))
	Select case(DatasetKey(k))
		case(2)		!HP98 dataset
			write(12,535)'L Tc', (LandauTc(k,j)  ,j=1,numPhCo(K))
			write(12,535)'L Sm', (LandauSmax(k,j),j=1,numPhCo(K))
			write(12,535)'L Vm', (LandauVmax(k,j),j=1,numPhCo(K))
		case(3)		!HP11 dataset
			write(12,*)' HP11 model values'
			write(12,*)' Swch = 0 ; Swch = 1 = Landau model (3 terms); Swch = 2 = Bragg-Williams model (6 terms)'
			write(12,536)'Swch',(HP11ModelSwitch(k,j),j=1,numPhCo(K))
			write(12,535)'    ',(HP11Model(k,j,1),j=1,numPhCo(K))
			write(12,535)'    ',(HP11Model(k,j,2),j=1,numPhCo(K))
			write(12,535)'    ',(HP11Model(k,j,3),j=1,numPhCo(K))
			write(12,535)'    ',(HP11Model(k,j,4),j=1,numPhCo(K))
			write(12,535)'    ',(HP11Model(k,j,5),j=1,numPhCo(K))
			write(12,535)'    ',(HP11Model(k,j,6),j=1,numPhCo(K))
			write(12,*)' HP11 DQF values '
			write(12,590)(DQFswitch(k,j),j=1,numPhCo(K))
590			FORMAT('DQF     ',15(I20))
			write(12,591)(DQFH(k,j),j=1,numPhCo(K))
591			FORMAT('DQFH    ',15(E20.5))
			write(12,592)(DQFS(k,j),j=1,numPhCo(K))
592			FORMAT('DQFS    ',15(E20.5))
			write(12,593)(DQFV(k,j),j=1,numPhCo(K))
593			FORMAT('DQFV    ',15(E20.5))
			write(12,594)((DQFH(k,j) - DQFS(k,j)*TK + DQFV(k,j)*PB),j=1,numPhCo(K))
594			FORMAT('DQFtotal',15(E20.5))
			write(12,585)(Pvalue(j),j=1,numPhCo(K))
		case default
		end select
536     format(' ',a5,15I20)


! 	if(iopt.ge.1000)then
! 	write(12,535)'L Tc', (LandauTc(k,j)  ,j=1,numPhCo(K))
! 	write(12,535)'L Sm', (LandauSmax(k,j),j=1,numPhCo(K))
! 	write(12,535)'L Vm', (LandauVmax(k,j),j=1,numPhCo(K))
! 	endif



      write(12,527)(HATTP(k,J),J=1,numPhCo(K))
527   Format(' HATTP',15E20.8)
      write(12,529)(SATTP(k,J),J=1,numPhCo(K))
529   Format(' SATTP',15F20.8)
      write(12,531)(VATTP(k,J),J=1,numPhCo(K))
531   Format(' VATTP',15F20.7)
      write(12,532)(GATTP(k,J),J=1,numPhCo(K))
532   Format(' GATTP',15E20.9)
      write(12,534)((HATTP(k,J)-TK*SATTP(k,j)),J=1,numPhCo(K))
534   Format(' H-TS ',15E20.9)


      write(12,522)(phCoName(k,J),J=1,numPhCo(K))
! 522   FORMAT('          ',12(A4,6X))
      write(12,523)(xPhCo(k,J),J=1,numPhCo(K))
523   FORMAT('  COMP',15(F20.7))
	Select case(DatasetKey(k))
		case(3)		!HP11 dataset only
585			FORMAT('PValue  ',15(E20.5))
			write(12,589)(exp(lnaIdeal(j)),j=1,numPhCo(K))
589			FORMAT('aIdeal  ',15(E20.5))
			write(12,586)(lnaIdeal(j),j=1,numPhCo(K))
586			FORMAT('lnaIdeal',15(E20.5))
			write(12,587)(RTlnGamma(j),j=1,numPhCo(K))
587			FORMAT('RTlnGama',15(E20.5))
			write(12,588)(dexp(RTlnGamma(j)/(R*TK)),j=1,numPhCo(K))
588			FORMAT('Gamma   ',15(E20.5))
		case default
		! no other special output
		end select

      write(12,524)(Dexp(lnAct(k,J)),J=1,numPhCo(K))
524   FORMAT(' Activ  ',15(E20.5))
      write(12,525)(lnAct(k,J),J=1,numPhCo(K))
525   FORMAT('  Ln(a) ',15(F20.7))

      write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
      IF(IMASS.EQ.1) then
      write(12,521)MP0(K),VP0(K)
521   FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
      endif

      write(12,*)'dudTPX ARRAY:'
      write(12,*)'     dudT           dudP           dudX2          dudX           dudX4...'
      do 550 j = 1,numPhCo(K)
      write(12,533)dudTPX(k,j,1),dudTPX(k,j,2),(dudTPX(k,j,L),L=3,numPhCo(K)+1)
533   FORMAT(15E15.5)
550   continue

	call SmolePhases			! Calculate G of each phase
	write(12,*)' G of phase = ',gPhase(K)
      go to 59

	endif		! end ioption = 3




! -------------------------------------------------------
!      Plot thermodynamic data as a function of P or T
! -------------------------------------------------------
	if(ioption.eq.4)then
! 	save current values so they can be reset after done      
      TKsave = TK
      TCsave = TC
      Pbsave = Pb
      Pfluidsave = Pfluid

3000  CONTINUE
	write(*,*)' Options to plot:'
	write(*,*)' 0 = return'
	write(*,*)' 1 = G'
	write(*,*)' 2 = H'
	write(*,*)' 3 = S'
	write(*,*)' 4 = V'
	write(*,*)' 5 = Cp'
	write(*,*)' 6 = uMineral (= Xi RTlnax)'
	read(*,*)(yPlot)
	IF(yPlot.EQ.0)go to 3099

	select case (yPlot)
	case(0)
		go to 3099
	case(1)
		ymax = 0
		ymin = -10000000
		nystep = 10
		Ylab = 'G'
	case(2)
		ymin = 0
		ymax = 10000000
		nystep = 10
		Ylab = 'H'
	case(3)
		ymin = 0
		ymax = 1000
		nystep = 10
		Ylab = 'S'
	case(4)
		ymin = 0
		ymax = 20
		nystep = 10
		Ylab = 'V'
	case(5)
		ymin = 0
		ymax = 1000
		nystep = 10
		Ylab = 'Cp'
	case(6)
		ymin = -10000.0
		ymax = 0
		nystep = 10
		Ylab = 'mu'
	case default
		go to 3099
	end select


	write(*,*)' Specify plot type'
	write(*,*)' 0 = return'
	write(*,*)' 1 = T versus property'
	write(*,*)' 2 = P versus property'
	write(*,*)' 3 = X versus property'
	read(*,*)xPlot
	select case (XPlot)
	case(0)
		go to 3099
	case(1)
		xmin = 200
		xmax = 1000
		nxstep = 16
		Xlab = 'Tc'
		write(*,*)'Input P in bars'
		read(*,*)PbPlot
	case(2)
		xmin = 0
		xmax = 20
		nxstep = 20
		XLab = 'Pkb'
	case(3)
		xmin = 0
		xmax = 1
		nxstep = 10
		nxdec = 3
		XLab = 'X'
	case default
		go to 3099
	end select
	

! 	Draw plot

	call SetPlot(iok)
!      	call Setplt(XOR,XMIN,XMAX,XLEN,nXstep,nXdec,XLAB,&
!       &YOR,YMIN,YMAX,YLEN,nYstep,nYdec,YLAB,Pltitl,iok)

      	if(iok.eq.1)go to 1  ! cancel button was selected

!      do PostScript file stuff
      	i=0
	myColor = currentColor		! save the previous color
	currentColor = 1		! Black for axes
	CALL PlotAxes(XYPlot)
	currentColor = myColor		! set back to the previous color
      

3003   	continue
      	WRITE(*,630)(MINREC(I),PHNAME(I),I=1,numPhMIF)
      	WRITE(*,*)' Options:'
      	iopt = 0
      	WRITE(*,*)' 0 = Return'
      	WRITE(*,*)' +n = Specify mineral number'
      	read(*,*)(iopt)
      	IF(IOPT.EQ.0)go to 3099

      	DO 3010 k = 1,numPhMIF
3010   	IF(MINREC(K).EQ.IOPT)GO TO 3011
!      	MINERAL IOPT IS NOT IN MINERAL FILE
      	call FSS_alert('ALERT!!',' The mineral number you specified is not in the list. Try again') 
      	GO TO 3003
3011   	CONTINUE
!      	THE MINERAL IS FOUND and its number is K
!      	write out current thermodynamic data

3025   	continue

        WRITE(*,*)
        WRITE(*,3060)MINREC(K),SITMUL(K),PHNAME(K)
3060    FORMAT(' ',I2,F8.1,4X,A32)
	  do 3062 j = 1,numPhCo(K)
        WRITE(*,3061)j,phCoName(k,J)
3062	  continue
3061    FORMAT(I6,2x,A8)


	select case (XPlot)
	case(0)
		go to 3099
	case(1,2)
      		write(*,*)' Input one phase component to plot'
      		j1=0
      		read(*,*)j1
      		if(j1.eq.0)go to 3003
      		WRITE(*,3063)phCoName(k,J1)
3063    		FORMAT(8A8,5x)
	case(3)
      		write(*,*)' Input numbers (two) of phases components to plot'
      		read(*,*)j1,j2
      		if(j1.eq.0)go to 3003
      		WRITE(*,3063)phCoName(k,J1),phCoName(k,J2)
	case default
	end select

! 	now make the plot
	select case (XPlot)
	case(1)			! T versus property
		xinc = 10
		TC = xmin - xinc
		Pb = Pbsave
		Pfluid = Pfluidsave
		Pb = PbPlot
	case(2)			! P versus property
		xinc = 100
		TC = TCsave
		write(*,*)'Input T in degC'
		read(*,*)TC
		Pb = xmin*1000. - xinc + .01		! the 0.01 is to ensure we don't try Pb = 0
		Pfluid = Pfluidsave
	case(3)			! X versus property
		xinc = .01
		TC = TCsave
		Pb = Pfluidsave
		write(*,*)'Input T in degC'
		read(*,*)TC
		write(*,*)'Input P in bars'
		read(*,*)Pb
		xPhCo(k,j1) = 0.001 - xinc
		xPhCo(k,j2) = .999 + xinc
	case default
	end select

	call PickAWEColor(myColor)
	CurrentColor = myColor				! CurrentColor is in PlotStuff.inc common block
	
	iup = 0

3100	continue
	select case (XPlot)
	case(1)			! T versus property
		TC = TC + xinc
		if(TC.gt.xmax)go to 3120
	case(2)			! P versus property
		Pb = Pb + xinc
		if(Pb.gt.xmax*1000.)go to 3120
	case(3)			! X versus property
		xPhCo(k,j1) = xPhCo(k,j1) + xinc
		xPhCo(k,j2) = xPhCo(k,j2) - xinc
		if(xPhCo(k,j1).ge.1)go to 3120
	case default
	end select

	TK = TC + 273.15

	select case(yPlot)
	case(1,2,3,4,5)
		izero=0
		ivol=0            ! compute everything
!		K=0
		!CALL AllKduTPX(ivol,izero)
		TK = TC + 273.15d0
		call CalculateCPTemp(TK)
		call duTPX(k,ivol,izero)
		if(izero.eq.1)then
			write(*,*)'Problem calculating thermodynamic properties. Line 1127 in Global.f95'
			pause
			go to 1
			endif	
	case(6)
	      	do 3913 J=1,numPhCo(K)
      		Xtemp(j) = 1.e-9		! zero out phase components
3913    	continue
      		Xtemp(j1) = xPhCo(k,j1)
      		Xtemp(j2) = xPhCo(k,j2)
!     		Activities
      		iok=0
      		call Activity(K,uMin,Xtemp,TK,PB,R,iok)
      		if(iok.eq.1)go to 1
	case default
	end select



	select case(XPlot)
	case(1)
		xxx = TC
	case(2)
		xxx = Pb/1000.
	case(3)
		xxx = xPhCo(k,j1)
	case default
	end select

		
	select case (Xplot)
	case(1,2)
		select case (yPlot)
		case (1)	! G
			yyy = GATTP(k,j1)		
		case (2)	! H
			yyy = HATTP(k,j1)		
		case (3)	! S
			yyy = SATTP(k,j1)		
		case (4) 	! V
			yyy = VATTP(k,j1)		
		case (5)	! Cp
!      		CP = A + B*T**-0.5 + C*T**-2 + D*T**-3 + E*T**-1   + F*T     + G*T**2
			yyy = aPhCoCp(k,j1) + bPhCoCp(k,j1)/dsqrt(TK) + cPhCoCp(k,j1)/(TK*TK)&
    		    & +dPhCoCp(k,j1)/(TK*TK*TK) + ePhCoCp(k,j1)/TK + fPhCoCp(k,j1)*TK + gPhCoCp(k,j1)*TK*TK
		case(6)
			yyy = uMin(j1)
		case default
		end select
	case(3)	
		select case (yPlot)
		case (1)	! G
			yyy = xPhCo(k,j1)*GATTP(k,j1) + xPhCo(k,j2)*GATTP(k,j2)		
		case (2)	! H
			yyy = xPhCo(k,j1)*HATTP(k,j1) + xPhCo(k,j2)*HATTP(k,j2)		
		case (3)	! S
			yyy = xPhCo(k,j1)*SATTP(k,j1) + xPhCo(k,j2)*SATTP(k,j2)		
		case (4) 	! V
			yyy = xPhCo(k,j1)*VATTP(k,j1) + xPhCo(k,j2)*VATTP(k,j2)		
		case (5)	! Cp
!      		CP = A + B*T**-0.5 + C*T**-2 + D*T**-3 + E*T**-1   + F*T     + G*T**2
			yyy = xPhCo(k,j1)*(aPhCoCp(k,j1) + bPhCoCp(k,j1)/dsqrt(TK) + cPhCoCp(k,j1)/(TK*TK) + &
     &			    dPhCoCp(k,j1)/(TK*TK*TK) + ePhCoCp(k,j1)/TK + fPhCoCp(k,j1)*TK + gPhCoCp(k,j1)*TK*TK)&
     &			+ xPhCo(k,j2)*(aPhCoCp(k,j2) + bPhCoCp(k,j2)/dsqrt(TK) + cPhCoCp(k,j2)/(TK*TK) + &
     &			    dPhCoCp(k,j2)/(TK*TK*TK) + ePhCoCp(k,j2)/TK + fPhCoCp(k,j2)*TK + gPhCoCp(k,j2)*TK*TK)
		case(6)
			yyy = xTemp(j1)*uMin(j1) + xTemp(j2)*uMin(j2)
		case default
		end select
	case default
	end select
	
	write(*,*)xxx,yyy
			
	call plot(XYPlot,xxx,yyy,iup)
!	if(jps.eq.1)call pplot(xxx,yyy,iup)
	iup=1
	go to 3100

3120	continue
!	if(jps.eq.1)call ppenup
	!call !scolor(1)
!	call FSS_ClosePicture(3)

	go to 3003


3099	continue
       	TC = TCsave
       	TK = TKsave
       	Pb = Pbsave
       	Pfluid = Pfluidsave
       	GO TO 1

	endif
	
	




! -------------------------------------------------------
!       Adjust SZERO, hPhCoZero, Vzero value for a phase
! -------------------------------------------------------
	if(ioption.eq.5)then
600   continue
      WRITE(*,630)(MINREC(I),PHNAME(I),I=1,numPhMIF)
630   FORMAT(3(I4,2X,A8,2x))
      WRITE(*,*)' Option to adjust hPhCoZero,Szero or Vzero of a phase'
      WRITE(*,*)' Input number of phase to adjust (0 to exit) '
      iopt3 = 0         !default
      READ(*,*)iopt3
      if(iopt3.eq.0)go to 1
      DO 610 k = 1,numPhMIF
610   IF(MINREC(K).EQ.IOPT3)GO TO 650
!     MINERAL IOPT3 IS NOT IN MINERAL FILE
      WRITE(*,*)' MINERAL NUMBER',IOPT3,' IS NOT IN LIST. TRY AGAIN'
      GO TO 600
650   CONTINUE
!     THE MINERAL IS FOUND and its number is K
!     write out current thermodynamic data

625   continue

          WRITE(*,*)
          WRITE(*,660)MINREC(K),SITMUL(K),PHNAME(K)
660       FORMAT(' ',I2,F8.1,4X,A32)
          WRITE(*,661)(phCoName(k,J),J=1,numPhCo(K))
661       FORMAT('          ',12(A8,7X))
          WRITE(*,662)(hPhCoZero(k,J),J=1,numPhCo(K))
662       Format(' hPhCoZero',12F15.2)
          WRITE(*,663)(sPhCoZero(k,J),J=1,numPhCo(K))
663       Format(' SZERO',12F15.5)
          WRITE(*,664)(vPhCoZero(k,J),J=1,numPhCo(K))
664       Format(' VZERO',12F15.5)
          WRITE(*,665)' ACP',(aPhCoCp(k,J),J=1,numPhCo(K))
          WRITE(*,665)' BCP',(bPhCoCp(k,J),J=1,numPhCo(K))
          WRITE(*,665)' CCP',(cPhCoCp(k,J),J=1,numPhCo(K))
          WRITE(*,665)' DCP',(dPhCoCp(k,J),J=1,numPhCo(K))
          WRITE(*,665)' ECP',(ePhCoCp(k,J),J=1,numPhCo(K))
          WRITE(*,665)' FCP',(fPhCoCp(k,J),J=1,numPhCo(K))
          WRITE(*,665)' GCP',(gPhCoCp(k,J),J=1,numPhCo(K))
665       format(' ',a5,12E15.3)
          WRITE(*,665)'  v1',(v1PhCo(k,J),J=1,numPhCo(K))
          WRITE(*,665)'  v2',(v2PhCo(k,J),J=1,numPhCo(K))
          WRITE(*,665)'  v3',(v3PhCo(k,J),J=1,numPhCo(K))
          WRITE(*,665)'  v4',(v4PhCo(k,J),J=1,numPhCo(K))

      write(*,*)' Inut zero to retain previous value'
      write(*,*)' Input number of phase component to change'
      j=0
      read(*,*)j
      if(j.eq.0)go to 600
      WRITE(*,661)phCoName(k,J)
      
      write(*,*)'Input new value for hPhCoZero(j) (0 to keep old value)'
      WRITE(*,662)hPhCoZero(k,J)
      read(*,*)Get8
      if(Get8.ne.0)hPhCoZero(k,j)=Get8

      write(*,*)'Input new value for sPhCoZero(j) (0 to keep old value)'
      WRITE(*,663)sPhCoZero(k,J)
      read(*,*)Get8
      if(Get8.ne.0)sPhCoZero(k,j)=Get8

      write(*,*)'Input new value for vPhCoZero(j) (0 to keep old value)'
      WRITE(*,664)vPhCoZero(k,J)
      read(*,*)Get8
      if(Get8.ne.0)vPhCoZero(k,j)=get8
      go to 625
      
      endif		!end ioption = 5






! 	--------Solve for H of phase component ------------------
	if(ioption.eq.6)then
	call Solve_for_H(ThmoFile)
	go to 1
	endif


! 	----------Close and open thermo file ---- 12 ------------------
	if(ioption.eq.7)then

! 	      call volset(volrefnum)        !set volume to startup directory
	      close(1)
	      write(*,*)' It is now OK to save Thermofile'
	      Write(*,*)' Hit return to continue with program execution'
	      pause

	      call OpenThermoFile

! c	      Open thermodynamic data file
!  	      FileName = ':Gibbs_Essentials:'//trim(THMOFILE)		!put path in .fig file, not here
! 	      OPEN (1,FILE=FileName,STATUS='OLD')
! 	      iok=0
! 	      call setthermofile(iok)
! 	      if(iok.eq.1)then
! 	            call fss_alert(' There is a problem reading the thermodynamic data file.')
! 	            endif
! c	     Load "fixed" minerals from thermodynamic data file
! c	     1 = alpha quartz, 2 = beta quartz, 3=kyanite,4=sill,5=andalusite
! 	      call loadfx
	
	      go to 1
	      endif


! -------------------------------------------------------
!      Make cordierite calculations at T and P 
! -------------------------------------------------------
	if(ioption.eq.8)then

      TKsave = TK
      TCsave = TC
      Pbsave = Pb
      Pfluidsave = Pfluid


!	Code assumes water = mineral #2 and cordierite = mineral #4 in data file
	write(12,*)'                               Cpdt                          Cp/T dt          Gattp'
	write(12,*)'   Tc        Pb          DryCrd     WetCrd         DryCrd          WetCrd       Water'
79    CONTINUE
      Ttemp = TC
      Ptemp = Pb
      write(*,*)' TC and Pb = ',Ttemp,Ptemp
      write(*,*)' Specify T in deg C and P in bars: TC,Pb'
      read(*,*)Ttemp,Ptemp
	if (Ttemp.eq.0)then
	        TC = TCsave
	        TK = TKsave
	        Pb = Pbsave
	        Pfluid = Pfluidsave
	        GO TO 1
	endif
	
      TC = Ttemp
      TK = TC+273.15
      Pb = Ptemp
      Pfluid = Pb

      izero=0
      ivol=0            ! compute everything
      K=0
      CALL AllKduTPX(ivol,izero)
      if(izero.eq.1)go to 1


	k=4
	write(12,777)TC,Pb,&
     &  (hattp(k,1)-hPhCoZero(k,1)),(hattp(k,2)-hPhCoZero(k,2)),&
     &  (sattp(k,1)-sPhCoZero(k,1)),(sattp(k,2)-sPhCoZero(k,2)),&
     &  GATTP(k,1)
777	Format(2F10.1,12E16.8)

      go to 79

	endif		! end ioption = 8


! ----------------
! -------------------------------------------------------
!      Adjust grain boundary energy 
! -------------------------------------------------------
	if(ioption.eq.9)then

902	continue
	write(*,*)' Current values are (Hse and Hde):'
	write(*,*)HseGB,HdeGB
	write(*,*)'input new values (0,0 to exit)'
	read(*,*)Hsetemp,Hdetemp
	if(Hsetemp.eq.0)then
		go to 1
		else
		HseGB = Hsetemp
		HdeGB = Hdetemp
		go to 902
		endif
	
	endif		! end ioption = 9


	go to 1
	
	end
	
	

	
