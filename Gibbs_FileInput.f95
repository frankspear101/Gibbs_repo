! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE BEGIN (IBEGIN,INEW)
    	use AWE_Interfaces
      implicit none
!     THIS VERSION USES joule AND J/BAR
!     THIS SUBROUTINE IS DESIGNED TO
!     READ IN AND SET UP THE MINERALS FOR THIS PROBLEM

!     ibegin = 1  Open data file from inside this routine
!     ibegin = 2  Assume data file has already been opened
!     ibegin = 3  Read input data from the internal master file
!			This option defunct in Gibbs2 because we now read everything into master file
! 	inew is error code
! 		inew = 0 - did not open a new file successfully
! 		inew = 1 - new file opened successfully
      
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "InputMineral.inc"
	include "Tangent.inc"

      integer*4 ibegin,inew,i,j,k,izero,jj,jjj,nr,L,LL,zz,iok,irec,irecStart,ier,status,kCur
      CHARACTER FileHeader*32,ppp*255,thmoFileTemp*32,plotName*20
	character*256 dummy
	character*256 aline
	integer*4 skipKFlu
	common /skipswitch/skipKFlu
	common /local_begin/i,j,k,izero,jj,jjj,nr,L,LL,zz,iok,irec,irecstart,ier,status,&
     &       fileheader,ppp,thmofiletemp,plotname,dummy	
!     initialize variables



	DO 9051 I=1,50
       	allFractl(i) = 0				! flag keyed to ALLX to tell whether a phase is fractionating - used in COMPUT2
        DO 9051 J=1,6
        ALLX(J,I)=0.D0
9051  	CONTINUE
      do 9052 K=1,15          ! zero phase arrays
         mpstart(k)=0.D0
         vpstart(k)=0.D0
         mp0(k)=0.D0
         mp1(k)=0.D0
         mp2(k)=0.D0
         vp1(k)=0.D0
         vp0(k)=0.D0
         vp2(k)=0.D0
         FRACTL(k)=0
9052     continue
      iopen  = 0
      numNew   = 0

!     Read component names and molecular wts from thermo data file
      rewind(1)

	
!    ========================================================================================
!     Begin input from disk
!    ========================================================================================
      if(IBEGIN.eq.1)then
      	OPEN(5,FILE='',STATUS='OLD',action='read',iostat=status)
      	if(status.ne.0)return
        INQUIRE (5, NAME=filein)
!      	psfilename = trim(filein)//'.ps'
	write(*,*)' You chose...',filein
      ENDIF

!     If Ibegin = 2 we will assume that the data file has already been opened elsewhere (mostly from MakeBuildGrid and MakeMyAFM routines)
!     If ibegin = 3 we are reading mineral info from the InputMineralFile and the mineral record
! 		from the array named InputMineral(16) (int*4) -
!		iBegin = 3 is now defunct in version Gibbs2

!     Read input file from disk

!     SETUP FILES (*.IN) HAVE THE FOLLOWING FORMAT
!        Gibbs v. 4 Input file    (This is a header)
!	if version = 5, this line contains the thermodynamic data file name 
!        Title line
!        TSTART PSTART
!        System COMPONENT names
!        Use System components (0 or 1)
!        IMASS
!  ***************************************  row of stars as a separator
!        PHASE #
!           VOLUME%,FRACTL,UseNewVariables
!           phase component names
!           Use Phase Component (0 or 1)
!           names of phase components or cations on sites (depending on type of mineral)
!           COMPosition (either mole fractions of phase components or cations on sites)
!  ***************************************  row of stars as a separator
!        PHASE #
!           VOLUME%,FRACTL,UseNewVariables
!           phase component names
!           Use Phase Component (0 or 1)
!           names of phase components or cations on sites (depending on type of mineral)
!           COMPosition (either mole fractions of phase components or cations on sites)
!  ***************************************  row of stars as a separator
!        PHASE #
!           VOLUME%,FRACTL,UseNewVariables
!           phase component names
!           Use Phase Component (0 or 1)
!           names of phase components or cations on sites (depending on type of mineral)
!           COMPosition (either mole fractions of phase components or cations on sites)
!  ***************************************  row of stars as a separator


	read(5,'(A32)')fileHeader
10	continue
	write(12,*)fileHeader
	read(5,'(A255)')AssembTitle
	write(12,'(A255)')AssembTitle
	read(5,*)TSTART,PSTART
	write(12,*)TSTART,PSTART
	
      	pfstart=pstart

!     CHOOSE SYSTEM COMPONENTS FOR THIS PROBLEM
	read(5,*)dummy         			! names of the elements
	write(12,*)dummy         			! names of the elements
	read(5,*)(UseSysCo(i),i=1,NoSysCoinFile)
	write(12,*)(UseSysCo(i),i=1,NoSysCoinFile)
	read(5,*)IMASS
	write(12,*)IMASS

!     Sort out the ones we need
      NC = 0
      do 2010 i = 1,NoSysCoinFile
      if(UseSysCo(i).eq.1)then
            NC=NC+1
            NOEL(NC) = i		! NOEL is pointer to system components used.
            endif
2010  continue


!     PICK MINERALS FOR THIS PROBLEM
       K=0
!       LOOP TO READ EACH PHASE
2120   CONTINUE
       K=K+1
2225   continue

	irecStart = 8		! all current files are v 4
	read(5,'(A4)',end=2003)dummy          	! row of stars for separator
	read(5,*,END=2003)NR                	! read a new mineral number
	!if(ilong(11).eq.1)then
		write(12,*)'  '
		write(12,*)'------------------------'
		write(12,*)'Reading Mineral number = ',NR
	!	endif
	IF(NR.EQ.0)GO TO 2003        		! end of list of minerals
	!if(ilong(11).eq.1)
	write(12,*)' Read VPstart(K),FRACTL(K),includeNew(K)'
	read(5,*,END=2003)VPStart(K),FRACTL(K),includeNew(K)
	write(12,*)VPStart(K),FRACTL(K),includeNew(K)

      IF(IMASS.EQ.0)THEN
         VPStart(K)=0.D0
         FRACTL(K)=0
         ENDIF

      call readminheader(K,NR,ier)
      if(ier.eq.1)then
      	inew = 0
		return
        	endif

!     read mineral composition from disk
!      NREC(NR) = K		! FSS 2/26/98 to index mineral number
      MINREC(K)=NR
      xPhCo(k,1)=1.
      xPhCoInitial(k,1)=1.
      zz = MinInFilePoint(K) - 1

	!if(ioutputPseudo.ge.2)
	write(12,*)K,NR,phName(K)
	write(12,*)'Read use PhCo'
	read(5,'(A256)')dummy          			! names of phase components
	write(12,*)dummy          			! names of phase components
	!if(ilong(11).eq.1)
	read(5,*,err=8100)(UsePhCo(k,j),j = 1,numPhCoInFile(K))
	write(12,*)(UsePhCo(k,j),j = 1,numPhCoInFile(K))

!     store a pointer that keeps track of which phase components are actually being
!     used relative to the thermodynamic data file.  This is only needed in
!     subroutine SAVEIN when we write out an input file.
!     This is also used in hard coded, mineral-specific activity subroutines such as
!           Chl1
!           Mus1
!           Bio1
      L = 0
      do 505 j = 1,numPhCoInFile(K)
      if(UsePhCo(k,j).eq.1)then
      	L=L+1
      	phCoUsedinFile(k,L)=j
      	endif
505   continue
      numPhCo(K)=L	! This is the number of phase components we are using

	call readmindata(K,NR)

	write(12,*)'Read Xvalues'
	read(5,'(A256)')aline         ! names of compositions we are reading in
	write(12,*)aline

      if(numSiteAtom(K).eq.0)then
		!if(ilong(11).eq.1)
		read(5,*)(xPhCo(k,j),J=1,numPhCoInFile(K))
		write(12,*)(xPhCo(k,j),J=1,numPhCoInFile(K))
		!strip out unused values
		jjj = 0
		do 501 j = 1,numPhCoInFile(K)
		if(UsePhCo(k,j).eq.1)then
			jjj = jjj+1
			xPhCo(k,jjj) = xPhCo(k,j)
			xPhCoInitial(k,jjj) = xPhCo(k,j)
			endif
501   continue
		write(12,*)(xPhCo(k,j),J=1,jjj)

		else
	!     Input atom site fractions for these atoms:'
		read(5,*)(SiteAtom(k,L),L=1,numSiteAtom(K))
		write(12,*)' Input site atoms from MIF'
		write(12,*)(SiteAtom(k,L),L=1,numSiteAtom(K))
	
	!     calculate mole fractions of components from site atom information
		write(12,*)' '
		write(12,*)'SiteAtomToX matrix (from CoTrans)'
		do 2210 j = 1,numPhCo(K)
		write(12,*)(SiteAtomToX(k,J,L),l=1,numSiteAtom(k))
		xPhCo(k,j) = 0
		do 2211 L = 1,numSiteAtom(K)
		xPhCo(k,j) = xPhCo(k,j) + SiteAtomToX(k,J,L)*SiteAtom(k,L)
2211   		continue
		xPhCoInitial(k,j) = xPhCo(k,j)
2210   		continue
		write(12,*)' '
		write(12,*)'Calculated PhaseComponents from SiteAtomToX matrix'
		write(12,*)(xPhCo(k,j),J=1,numPhCo(K))
	
		endif


      GO TO 2120          !go get another mineral


2003  CONTINUE
!     numPh IS TOTAL # PHASES FOR THIS PROBLEM
	write(*,*)' ------------------'
	write(*,*)' EOF on input file'
	write(*,*)' ------------------'
	write(*,*)' '
      numPh=K-1
!	Assign all phases to the current assemblage to start
	do 2005 k = 1,numPh
	asmCurrent(k) = k
2005	continue
	numPhMIF = numPh		! this is the total number of phases in the MIF file
!     SET INITIAL VALUES
      TC=TSTART
      TK=TC+273.15D0
      PB=PSTART
!
!     SET UP FOR PFLUID
      KFLU=0
	PFStart = PStart
	PFluid = PFStart
	PFluidSwitch = 0
	TandP=2
	call FluidPresent	!returns KFlu>0 if a fluid is present
	if(KFlu.gt.0.and.skipKFlu.eq.0)then	! ask user how to deal with fluid pressure
		!call SetPFluid(iok)		! routine to set P fluid variables --- Needs to be implemented
		if(iok.eq.1)go to 60
		IF(PFluidSwitch.gt.0)then
			TandP = 3
			PFStart = PFluid
			endif
		endif

60	continue
	CLOSE(5)

!    ========================================================================================
!     End input from disk
!    ========================================================================================
!
!     Adjust arrays CONAME(L) and MOLWT(L) so that they only contain
!     system components for this problem
      DO 70 L=1,NC
      MOLWT(L)=MOLWTInFile(NOEL(L))
      OxToElWt(L)=OxToElWtInFile(NOEL(L))
      NumCatInOxide(L)=NumCatInOxideInFile(NOEL(L))
      CONAME(L)=CONameInFile(NOEL(L))
	do 71 i = 1,NoCompPlotVars
	CompPlotCoeff(i,L) = CompPlotCoeffInFile(i,NOEL(L))
71	continue
70    continue

! forget this part - I don't think its the way to go
!	for grain boundary problems we must include oxygen as a system component
!	but we don't want to have an oxygen mass balance equation
	noOxygenMassBalance = 0
!	do 2011 i = 1,nc
!	write(*,*)coName(i)
!      if(trim(coName(i)).eq.'Ox')then
!		noOxygenMassBalance = 1
!		write(*,*)' noOxygenMassBalance = 1 (no mass balance for oxygen) '
!	      endif
!2011	continue
	

!     compute initial volumes of phases in assemblage
!     and compute moles of phases corresponding to that volume
!     set MPSTART,mp0,mp1 AND MP2 to initial values
	izero=0
      call AllKduTPX(1,izero)

      if (izero.eq.1)then
      		call FSS_alert('ALERT!!','There is a problem in subroutine dudTPX.  A phase has a zero or negative value')
            	return
            	endif


	if(bulkCompSwitch.eq..true.)then	! we have set an initial bulk composition
! 	   	Convert bulk composition to moles
		bulkCompWtSum=0.D0
		do i=1,nc
			bulkCompMoles(i)=bulkCompWt(i)*numCatInOxide(i)/molwt(i)
			bulkCompWtSum=bulkCompWtSum + bulkCompWt(i)
			molesStart(i) = bulkCompMoles(i)	! molesStart is the array used in Compute
			bulkCompMolesStart(i) = molesStart(i)		! bulkCompMolesStart is the initial input composition in case we need to reset to starting values
			bulkCompWtStart(i) = bulkCompWt(i)		! save the initial composition in case we need it later
			end do
		iMass = 1	! set mass balance switch to "on"
		nvar=nvar+numPh
		neq=neq+nc
!		do 350 k = 1,numPh
		Do 350 kCur = 1,numPh
		k = asmCurrent(kCur)
		mpStart(k) = .10D0			! set to 1 (must be computed later in compute)
		mp0(k) = .10D0
		mp1(k) = .10D0
		mp2(k) = .10D0
		vpStart(k) = vmol(k)*10.0d0/mpStart(k)	! vmol is j/bar; vpStart is in cm^3
		vp0(K)=VPstart(K)
		vp1(k)=VPstart(K)
		vp2(k)=VPstart(K)
		select case (Fractl(K))
			case(0)	! not fractionating
			case(1)	! fractionating and staying in rock
				mp0(k) = 0.0d0
				vp0(K) = 0.0d0
			case(2)	! fractionating and leaving the rock
				mp0(k) = 0.0d0
				mp1(k) = 0.0d0
				vp0(K) = 0.0d0
				vp1(K) = 0.0d0
			case default
			end select
350	   	continue

		else

!      		do 300 k = 1,numPh
		Do 300 kCur = 1,numPh
		k = asmCurrent(kCur)
!       	set initial volumes for all volume arrays
! 		Adjust volumes for cases where phases are fractionating
		select case (Fractl(K))
			case(0)	! not fractionating
			 vp0(K)=VPstart(K)
			 vp1(k)=VPstart(K)
			 vp2(k)=VPstart(K)
			case(1)	! fractionating and staying in rock
			 vp0(K) = 0.0d0
			 vp1(k)=VPstart(K)
			 vp2(k)=VPstart(K)
			case(2)	! fractionating and leaving the rock
			 vp0(K) = 0.0d0
			 vp1(K) = 0.0d0
			 vp2(k)=VPstart(K)
			case default
			end select
!       	compute initial moles of phases
!       	moles per cm**3 (note VMOL is in joule/bar =CC/10)
!    		note: if this code is changed, it must also be changed in routine CHANGE
      	if(vmol(K).lt.1.e-20)then
        	write(*,*)'There is a problem with vmol(K)'
        	write(*,*)' K, phase = ',K,phName(k)
        	pause 'This problem should be aborted--Hit return to continue'
        	else
        	mp0(K)=vp0(K)/(vmol(K)*10.0D0)
        	mp1(k)=vp1(K)/(vmol(K)*10.0D0)
        	mp2(k)=vp2(K)/(vmol(K)*10.0D0)
        	MPSTART(k)=mp0(k)			! mpStart is the starting moles of each phase
      	  	endif
300   	continue
! 	   	Compute the starting number of moles of each system component
	   	do 310 i=1,nc
	   	molesStart(i)=0.D0			! molesStart is the starting moles of each system component
		Do 312 kCur = 1,numPh
		k = asmCurrent(kCur)
	   	do 315 j=1,numPhCo(K)
	   	molesStart(i) = molesStart(i) + MPStart(K)*xPhCo(k,j)*comp(k,j,i)
315	   	continue
312	   	continue
310	   	continue

		endif		! end bulkCompSwitch true or false


!      	call names		(called below)

	masserroralert = 0		! we have not issued a mass balance error alert yet

!     compute total number of phase components
      	NP=0
	Do 653 kCur = 1,numPh
	k = asmCurrent(kCur)
      	NP=numPhCo(K)+NP
653   	CONTINUE

!     compute number of independent compositional (dX) variables
      	NX=NP-numPh
!     compute total number of variables
      	NVAR= TANDP + NX
      	if(IMASS.EQ.1)NVAR=NVAR + numPh

!     (note: variable NEQ is calculated in subroutine rxn after reactions are figured)

!        ARRAY ALLX CONTAINS VALUES OF THE FOLLOWING
!     set all columns of ALLX(i,j) to starting values
!     Definition of ALLX variable:
!     1,50   current value
!     2,50     starting value
!     3,50   ref at start of contour
!     4,50   user selected reference
!     5,50   previous finite diff point
!     6,50     (not used)

      	DO 444 I=1,6
      	ALLX(I,1)=TSTART
      	ALLX(I,2)=PSTART
      	IF (KFLU.NE.0)ALLX(I,3)=PFStart
!     	THESE ARE THE X'S
      	L=TANDP
	Do 441 kCur = 1,numPh
	k = asmCurrent(kCur)
      	DO 555 J=2,numPhCo(K)
      	L=L+1
      	ALLX(I,L)=xPhCo(k,j)
555   	CONTINUE
441   	CONTINUE
!     	THESE ARE MOLES
      	L=TANDP+NX
      	IF(IMASS.EQ.1) THEN
		Do 442 kCur = 1,numPh
		k = asmCurrent(kCur)
            	L=L+1
            	ALLX(I,L)=MP0(K)
            	if(fractl(k).eq.1.or.fractl(k).eq.2)then
            		!This is done as a flag for COMPUT2 where it checks for Newton convergence (around line 722)
            		! If FRACTL = 1 or 2 then NEWTON doesn't converge properly.
            		! This flag will omit moles of fractionating phases from Newton convergence
            		allx(6,L)=1
            		else
            		allx(6,L)=0
            		endif
442         	CONTINUE
        	ENDIF


!     Set allFractl switch if necessary
      IF(IMASS.EQ.1) THEN
      		L=TANDP+NX
!            	DO 443 k = 1,numPh
		Do 443 kCur = 1,numPh
		k = asmCurrent(kCur)
            	L=L+1
            	if(fractl(k).eq.1.or.fractl(k).eq.2)then
            		!This is done as a flag for COMPUT2 where it checks for Newton convergence (around line 1738)
            		! If FRACTL = 1 or 2 then NEWTON doesn't converge properly.
            		! This flag will omit moles of fractionating phases from Newton convergence
            		allFractl(L)=1
            		else
            		allFractl(L)=0
            		endif
443         	CONTINUE
        	ENDIF

444   CONTINUE

	call SetNewVar()

	CALL NAMES
	NVAR= NVAR + numNew
	INEW=1		! we successfully opened a new file

! 	Sort through the variable list for this new problem and try to find the same variables that were plotted in the last problem
	do 6010 i = 1,noTieMinerals
	do 6015 j = i,nvar
	plotName = trim(VN1(j))//VN2(j)
	if(trim(plotIndexName(i,1)).eq.trim(plotName))then
		plotIndex(i,1) = j
		endif
	if(trim(plotIndexName(i,2)).eq.trim(plotName))then
		plotIndex(i,2) = j
		endif
6015	continue
6010	continue

      	RETURN
      
8100	continue      
      	write(*,*)" Error reading input data file."
   	write(*,*)" Most likely this results from not having the correct components specified"
   	write(*,*)" Mineral number = ",NR
   	write(*,*)" Check the input file against the thermodynamic data file for consistency"
   	pause 'hit return to end program'
	stop
	
      	END


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE readminheader (K,NR,ier)
!     THIS VERSION USES joule AND J/BAR
!     THIS SUBROUTINE IS DESIGNED TO
!     READ read a mineral from the thermodynamic datafile
	implicit none
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
	include "Output.inc"
! ********************************************************************
!      local variables
	integer*4 i,k,nr,kk,newline,deltaline,ier
	CHARACTER aline*64
	common /localreadminheader/ i,kk,newline,deltaline,aline

	ier = 0
	newline = mineralposition(NR) 		! position of mineral in datafile

	rewind(1)
	do 426 i=1,newline-1
	read(1,*)
426     continue

                        
!     	read data for new mineral      
        READ(1,*)KK,PHNAME(K)
        write(*,*)KK,PHNAME(K)
      	if(kk.ne.NR)then
      		write(*,*)' The mineral in the input file could not be found in the thermo file.'
     		write(*,*)' Check that the correct thermo file is loaded.'
 	 	write(*,*)' Check to see the correct thermo file is loaded' 
             	write(*,*)' There is an error reading the thermodynamic data file in Subroutine ReadMin'
	    	write(*,*)'Thermodynamic data file name:'
	    	write(*,*)thmofile
            	write(*,*)' Tried to read mineral'
            	write(*,*)NR
! 134       	FORMAT(I5,i12,5X,A32)
	    	write(*,*)' Actually read number and name'
            	write(*,*)KK,' ',phname(K)
            	write(*,*)
            	write(*,*)' 4 subsequent lines in data file'
		read(1,'(A)')aline
		write(*,*)aline
		read(1,'(A)')aline
		write(*,*)aline
		read(1,'(A)')aline
		write(*,*)aline
		read(1,'(A)')aline
		write(*,*)aline
            	write(*,*)
!             	write(12,*)'Hit return to continue'
!           	read (12,*)
	    	ier = 1	! return the error code
	    	return
	   	! stop
        endif
		
        numNewInFile(K) = 0
        numSiteAtom(K) = 0
        READ(1,*)numPhCoInFile(K),PhaseType(K),NumSiteAtom(K),NumReciprocal(K)
        Write(*,*)numPhCoInFile(K),PhaseType(K),NumSiteAtom(K),NumReciprocal(K)

!	if(special(K).eq.14) then	! this is coded so that special case = 14 is the grain boundary phase

!		if(imass.eq.0)then
!			call fss_alert('ALERT!!','You must be using mass balance to use the Grain Boundary model')
!			endif
!		grainBoundaryProblem = .TRUE.
!		grainBoundaryPhase   = K
!		endif

	ier = 0
	return
	end
	

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE readmindata (K,NR)
!     THIS VERSION USES joule AND J/BAR
!     THIS SUBROUTINE IS DESIGNED TO
!     READ read a mineral from the thermodynamic datafile
      
      implicit none
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
!********************************************************************
!      local variables
      integer*4 i,j,k,jj,jjj,nr,L,LL,zz,n,nn,n1
      CHARACTER dummy*16
      real*8 Xtemp1(20),Xtemp2(20)

	integer*4 SFID1Temp,SFID2Temp,numSFtermsTemp
	real*8 SFWHTemp,SFWSTemp,SFWVTemp

!      REAL*8 A(15,12),B(15,15)
!	common /CoTransCommon/a,b
	include "CoTransCommon.inc"

	common /local6/xtemp1,xtemp2,i,j,jj,jjj,L,LL,zz,n,nn,dummy
      real*8 coValue(sysCoMax)
	common /coValueCommon/covalue

	DatasetKey(k) = 0							! This will tell whether the mineral type is SPaC=1, HP98=2, HP11 = 3
	if(PhaseType(K).gt.0.and.PhaseType(K).LT.100)DatasetKey(k) = 1
	if(PhaseType(K).ge.100.and.PhaseType(K).LT.200)DatasetKey(k) = 2
	if(PhaseType(K).ge.200.and.PhaseType(K).LT.300)DatasetKey(k) = 3
	if(DatasetKey(k).eq.0)then
		write(*,*)' Problem with DatasetKey in subroutine ReadMinData'
		write(*,*)' Phase = ',k,minrec(k),PhName(k),PhaseType(k)
		pause ' Hit return to continue'
		endif

	Select Case(numPhCoInFile(k))
	

      case (1)					! All 1 compontnt minerals minerals   
	jjj = 0
	j = 1
      	READ(1,'(A16)')phCoNameinFile(k,j)
      	write(*,'(A16,I8)')'Non multisite phase  ',phCoNameinFile(k,j),phasetype(K)
      	if(UsePhCo(k,j).eq.1)then	! for a 1 component mineral we should ALWAYS use the phase component
	      	jjj=jjj+1
		phCoName(k,jjj)=phCoNameinFile(k,j)
		write(*,*)'Read HSV component=',phCoName(k,jjj)
		read(1,*)hPhCoZero(k,jjj),sPhCoZero(k,jjj),vPhCoZero(k,jjj)
		write(*,*)hPhCoZero(k,jjj),sPhCoZero(k,jjj),vPhCoZero(k,jjj)
		READ(1,*)aPhCoCp(k,jjj),bPhCoCp(k,jjj),cPhCoCp(k,jjj),dPhCoCp(k,jjj),ePhCoCp(k,jjj),fPhCoCp(k,jjj),gPhCoCp(k,jjj)
		write(*,*)aPhCoCp(k,jjj),bPhCoCp(k,jjj),cPhCoCp(k,jjj),dPhCoCp(k,jjj),ePhCoCp(k,jjj),fPhCoCp(k,jjj),gPhCoCp(k,jjj)
		READ(1,*)v1PhCo(k,jjj),v2PhCo(k,jjj),v3PhCo(k,jjj),v4PhCo(k,jjj)
		write(*,*)v1PhCo(k,jjj),v2PhCo(k,jjj),v3PhCo(k,jjj),v4PhCo(k,jjj)
		select Case (DatasetKey(k))
		   case(2)			! HP98 dataset
!     			read Landau coefficients for Holland&Powell dataset (min number > 1000)
			read(1,*)LandauTc(k,jjj),LandauSmax(k,jjj),LandauVmax(k,jjj)
			write(*,*)LandauTc(k,jjj),LandauSmax(k,jjj),LandauVmax(k,jjj)
!
		   case(3)    ! HP11 dataset
			read(1,*)HP11ModelSwitch(k,jjj),(HP11Model(k,jjj,L),L=1,6)
			write(*,*)HP11ModelSwitch(k,jjj),(HP11Model(k,jjj,L),L=1,6)
			read(1,*)DQFswitch(k,jjj),DQFH(k,jjj),DQFS(k,jjj),DQFV(k,jjj)
			write(*,*)DQFswitch(k,jjj),DQFH(k,jjj),DQFS(k,jjj),DQFV(k,jjj),'   ! DQF switch, H, S, V'
		   case Default     ! SPaC
		   ! no additional terms in this part of SPaC input files
		   end select

! 		routine to read composition of phase component and put into array COMP
		call ReadComp()
      		DO 10 L=1,NC
      		COMP(k,jjj,L)=coValue(NOEL(L))
10   		CONTINUE

      		else					! skip this phase component

		write(*,*)'Why are we here -- we should always use phase components in 1 component minerals'
		pause 'hit return to continue'
	      	read(1,*)dummy				! line of dashes
		do while(dummy(1:5).ne.'-----')
   	   	read(1,*)dummy				! H,S,V
		end do
		endif		! end UsePhCo.eq.1


	Case default 		! all multi-site minerals
 	  read(1,*)dummy			! site atom names
	  read(1,*)dummy			! site atom index
	  read(1,*)(SiteMultiplicity(k,L),L=1,numSiteAtom(K))
	  write(*,*)(SiteMultiplicity(k,L),L=1,numSiteAtom(K))
	  read(1,*)dummy			! line of dashes
	  jjj = 0
	  do 20 j=1,numPhCoInFile(K)
	  READ(1,'(A16)')phCoNameinFile(k,j)
	  write(*,'(A16)')phCoNameinFile(k,j)
	  if(UsePhCo(k,j).eq.1)then
      		jjj=jjj+1
      		phCoName(k,jjj)=phCoNameinFile(k,j)
      		read(1,*)hPhCoZero(k,jjj),sPhCoZero(k,jjj),vPhCoZero(k,jjj)
      		write(*,*)hPhCoZero(k,jjj),sPhCoZero(k,jjj),vPhCoZero(k,jjj)
      		READ(1,*)aPhCoCp(k,jjj),bPhCoCp(k,jjj),cPhCoCp(k,jjj),dPhCoCp(k,jjj),ePhCoCp(k,jjj),fPhCoCp(k,jjj),gPhCoCp(k,jjj)
      		write(*,*)aPhCoCp(k,jjj),bPhCoCp(k,jjj),cPhCoCp(k,jjj),dPhCoCp(k,jjj),ePhCoCp(k,jjj),fPhCoCp(k,jjj),gPhCoCp(k,jjj)
      		READ(1,*)v1PhCo(k,jjj),v2PhCo(k,jjj),v3PhCo(k,jjj),v4PhCo(k,jjj)
      		write(*,*)v1PhCo(k,jjj),v2PhCo(k,jjj),v3PhCo(k,jjj),v4PhCo(k,jjj)

		select Case (DatasetKey(k))
		case(2)			! HP98 dataset
!     			read Landau coefficients for Holland&Powell dataset (min number > 1000)
			read(1,*)LandauTc(k,jjj),LandauSmax(k,jjj),LandauVmax(k,jjj)
			write(*,*)LandauTc(k,jjj),LandauSmax(k,jjj),LandauVmax(k,jjj)
!
		case(3)    ! HP11 dataset
			read(1,*)HP11ModelSwitch(k,jjj),(HP11Model(k,jjj,L),L=1,6)
			write(*,*)HP11ModelSwitch(k,jjj),(HP11Model(k,jjj,L),L=1,6)
			read(1,*)DQFswitch(k,jjj),DQFH(k,jjj),DQFS(k,jjj),DQFV(k,jjj)
			write(*,*)DQFswitch(k,jjj),DQFH(k,jjj),DQFS(k,jjj),DQFV(k,jjj),'   ! DQF switch, H, S, V'

		case Default     ! SPaC
		! no additional terms in this part of SPaC input files
		end select

! 		routine to read composition of phase component and put into array COMP
		call ReadComp()
      		DO 25 L=1,NC
      		COMP(k,jjj,L)=coValue(NOEL(L))
25   		CONTINUE
      		read(1,*)ActivityConstant(k,jjj),(Alpha(k,jjj,L),L=1,numSiteAtom(K))
      		Write(*,*)ActivityConstant(k,jjj),(Alpha(k,jjj,L),L=1,numSiteAtom(K))

		select case (DatasetKey(k))
		case(1)			! SPaC
			if(numReciprocal(K).gt.0)then
				do 26 i = 1,numReciprocal(K)
				read(1,*)RecipConstant(k,jjj,i),(RecipIndex(k,jjj,i,L),L=1,numSiteAtom(K))
				write(*,*)RecipConstant(k,jjj,i),(RecipIndex(k,jjj,i,L),L=1,numSiteAtom(K))
26				continue
				endif
	
		case(3)		!HP11
			read(1,*)ASF(k,jjj,1),ASF(k,jjj,2),ASF(k,jjj,3)			! ASF size parameters
			Write(*,*)ASF(k,jjj,1),ASF(k,jjj,2),ASF(k,jjj,3),   '	! ASF size parameters'
			if(ASF(k,jjj,1).eq.0)then
				write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
				write(*,*)' First ASF parameter of this phase component = 0...that is a no-no'
				write(*,*)' Aborting program'
				write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
				stop
				endif
		case default		! no terms here for HP98 input files
		end select
		
	      	read(1,*)dummy				! line of dashes

      		else						! we don't want this phase co, but we must advance file
!		find the next line of dashes
	      	read(1,*)dummy				! line of dashes
		do while(dummy(1:5).ne.'-----')
   	   	read(1,*)dummy				! H,S,V
		end do
      		endif					! end usephco=1
20   continue

!	We have read the stuff for each individual phase component
!	Now read the phase specific stuff
!		Sites as a function of components
!		Margules terms
!		Reciprocal terms
!		New variables
!     Read the atom to site equations (these are valid for all subsystems)
      	read(1,*)dummy
      	do 70 L = 1,numSiteAtom(K)
      	read(1,*)(Xtemp1(j),j=1,numPhCoInFile(K)),SiteAtomName(k,L)
	write(*,*)(Xtemp1(j),j=1,numPhCoInFile(K)),SiteAtomName(k,L)
      	!strip out unused phase components
      	jj=0
      	do 71 j=1,numPhCoInFile(K)
      	if(UsePhCo(k,j).eq.1)then
      	jj = jj+1
      	XToSiteAtom(k,L,jj)=Xtemp1(j)
      	endif
71   	continue
70   	continue

! Here is where we will insert code to calculate component transformation matrix
! Take the XToSiteAtom array and pass a subset to a routine that will do the reduction
! Return the transformation array, which can be put into the SiteAtomToX array


	do 91 L = 1,numSiteAtom(K)
	do 91 j = 1,numPhCo(K)
	A(L,J) = XToSiteAtom(k,l,j)
91	continue
	call CoTrans(numSiteAtom(K),numPhCo(K))
	do 92 j = 1,numPhCo(K)
	do 92 L = 1,numSiteAtom(K)
	SiteAtomToX(k,j,L) = B(j,L)
92	continue
!	endif

	read(1,*)dummy				! line of dashes
	select case (DatasetKey(k))
	case(1)			!SPaC
!		Read Margules parameters
		read(1,*)numMargulesSites(k)
		write(*,*)' Number of Margules sites = ',numMargulesSites(k)
		if(numMargulesSites(k).gt.0)then
			do 120 jj = 1,numMargulesSites(k)
			read(1,*)numMargulesSiteCats(k,jj),(MargulesSiteCats(k,jj,L),L = 1,numMargulesSiteCats(k,jj))
			read(1,*)numMargulesWterms(k,jj)
			do 121 L=1,numMargulesWterms(k,jj)
			read(1,*)MargulesWIndex(k,jj,L),MargulesWH(k,jj,L),MargulesWS(k,jj,L),MargulesWV(k,jj,L)
121			continue
120			continue
			endif
	read(1,*)dummy			! line of dashes
!	read reciprocal terms
	write(*,*)' Number of reciprocal terms = ',numReciprocal(k)
	if(numReciprocal(k).gt.0)then
		do 130 i = 1,numReciprocal(k)
		read(1,*)WHrecip(k,i),WSrecip(k,1),WVrecip(k,1)
130		continue
		endif
	case(2,3)		!HP98 and HP11

!	Read excess terms for symmetric formalism (SF)
	read(1,*)numSFtermsTemp
	write(*,*)' number of SF terms = ',numSFtermsTemp
	if(numSFtermsTemp.gt.0)then
!		read(1,'(A64)')dummy	!names of nonideal terms
!		write(*,*)dummy	!names of nonideal terms
		n1 = 0
		do 230 n = 1,numSFtermsTemp
		read(1,*) SFWHTemp,SFWSTemp,SFWVTemp,SFID1Temp,SFID2Temp
		write(*,*)SFWHTemp,SFWSTemp,SFWVTemp,SFID1Temp,SFID2Temp
		if(UsePhCo(k,SFID1Temp).eq.1.and.UsePhCo(k,SFID2Temp).eq.1)then
!			we are using both phase components for this excess term
			n1 = n1+1
			SFWH(k,n1) = SFWHTemp
			SFWS(k,n1) = SFWSTemp
			SFWV(k,n1) = SFWVTemp
			do 232 j = 1,numPhCo(k)
			if(phCoUsedInFile(k,j).eq.SFID1Temp)then
				SFID(k,n1,1) = j
				endif
			if(phCoUsedInFile(k,j).eq.SFID2Temp)then
				SFID(k,n1,2) = j
				endif
232			continue
!			SFID(k,n1,2) = SFID2Temp
			endif
230   		continue
		numSFterms(k) = n1
		write(*,*)' SF terms used in this problem'
		do 231 n = 1,numSFterms(k)
		write(*,*)SFWH(k,n),SFWS(k,n),SFWV(k,n),SFID(k,n,1),SFID(k,n,2)
231		continue

		endif					! end reading nonideal terms

	case default
	end select		! done reading excess terms
	
      	read(1,*)dummy				! line of dashes
!	Read new variables
	read(1,*)numNewInFile(k)
	write(*,*)' Number of new variables = ',numNewInFile(k)
!     	read definitions of new variables
      	if(numNewInFile(K).gt.0.and.IncludeNew(k).eq.1)then
	      	do 300 i = 1,numNewInFile(K)
        	if(numNew.eq.31)then
			write(*,*)'There are more than 30 New Variables. Additional new variables will be omitted from the list'
               	 	call FSS_alert('ALERT!!','There are more than 30 New Variables. ')
         		go to 300
       			endif
      		numNew = numNew + 1
      		read(1,'(A)')newVarName(k,i)
!      		NewName(k,i)=NewVarTitle(k,i)(1:16)
      		read(1,*)NewNumA(k,i),(Xtemp1(j),j=1,numPhCoInFile(K))
      		read(1,*)NewDenB(k,i),(Xtemp2(j),j=1,numPhCoInFile(K))
      		! strip out unused phase components
      		jj=0
      		do 302 j=1,numPhCoInFile(K)
      		if(UsePhCo(k,j).eq.1)then
      			jj = jj+1
      			NewNumaj(k,i,jj)=Xtemp1(j)
      			NewDenbj(k,i,jj)=Xtemp2(j)
      			endif
302   		continue
300   		continue
		numNewInPhK(k) = numNewInFile(k)
      		endif

      end select					! end special case multi component minerals
      RETURN
      END

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine ReadComp()
	! program to decode mineral bulk composition
	implicit none
! !*****************************************
 	include "Assemb.inc"
! !*****************************************
	character*8 coNameTemp,coValueChar
	character*16 compChar(sysCoMax),compCharTemp
	character*128 compLine
	integer*4 lineLength,ncomp,i,j,k,iLen
	
	common /localReadComp/ coNameTemp,coValueChar,compChar,compCharTemp,compLine,&
     &             lineLength,ncomp,i,j,k,iLen
     
      real*8 coValue(sysCoMax)
	common /coValueCommon/covalue

! 	write(*,*)(coName(i),i=1,5)
! 	open(1,file='',status='old')
1	continue

	do 2 i = 1,noSysCoInFile
	coValue(i)=0.0d0
2	continue

	read(1,100)compLine
100	format(A128)
	write(*,*)compLine
	! count the number of system components in this line
	! each system component has an '=' sign
	lineLength=len(trim(compLine))
	ncomp=0
	do 5 i = 1,lineLength
	if(compLine(i:i).eq.'=')then
		ncomp=ncomp+1
	endif
5	continue
	
	! now look at each one individually
	read(compLine,*)(compChar(i),i=1,ncomp)	! internal read
101	format(A)

	do 20 k = 1,ncomp
	compCharTemp=compChar(k)
	compCharTemp=adjustL(compCharTemp)
	iLen = len(trim(compCharTemp))
	do 10 i = 1,iLen
	if(compCharTemp(i:i).eq.'=')then
		coNameTemp = compCharTemp(1:i-1)
		! figure out which one this is
		coNameTemp=adjustL(coNameTemp)
		do 15 j = 1,NoSysCoInFile
		if(trim(coNameTemp).eq.trim(coNameInFile(j)))then
			! read value for this component
			coValueChar = compCharTemp(i+1:iLen)
			! read value into array
			read(coValueChar,*)coValue(j)
			coValue(j)=coValue(j)
			go to 10
		endif	
15		continue
	call fss_alert('ALERT!!','Could not decode composition of this mineral')
	endif
10	continue
20	continue
	
! 	write(*,*)compLine
! 	write(*,110)(coValue(j),j=1,noSysCoInFile)
110	format(40F6.2)
! 	pause ' Hit return to continue'
	return
	end
	
	


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SAVEIN (icall,Kskip)
      implicit none
!     ROUTINE TO save the input problem as a *.in file

!     Routine is called from Main menu, where the entire file is saved
!     Routine is also called from CHANGE, where a subset of the file is saved and reread later

!     switches
!     icall = 0 called from Main menue - open new file and save all
!     icall = 1 called from "add mineral" - write entire file to a temp file
!     icall = 2 called from "delete mineral" - write to a temp file and skip mineral Kskip
!     icall = 3 called from "OverstepGarnet" - write to a temp file but skip call to duTPX 

	include "Assemb.inc"
	include "Gibbsfiles.inc"
! ***********************************
!     local variables
      integer*4 icall,Kskip
      integer*4 i,j,k,L,IZERO,ivol,status,kCur
      REAL*4 xx(12)
      character Xsymb*1

	if(icall.eq.3)go to 4 
!        	COMPUTE MOLAR VOLUME OF PHASE FOR NEW COMPOSITION, T, P AND PF
            	izero=0
            	ivol=1
            	CALL AllKduTPX(ivol,izero)
!         	DO 3630 k = 1,numPh
		Do 3630 kCur = 1,numPh
		k = asmCurrent(kCur)
!     		note: vp0, vp1 and vp2 are volumes of phases based on the moles present
!           	these are not modes!!
            	VP0(K)=MP0(K)*VMOL(K)*10.D0   	! note that molar volume of phase (VMOL) is in J/mol
            	vp1(k)=mp1(k)*vmol(k)*10.D0         ! and VP0, VP1, and VP2 are in cm^3
            	vp2(k)=mp2(k)*vmol(k)*10.D0
3630      	CONTINUE

4		continue
!      		WRITE(*,*)' Input name of the output save file (IN-** format)'
!     		if icall = 1 or 2 we will assume the file is already open (unit = 5)
      		if(icall.eq.0)then
      			OPEN(5,FILE='',STATUS='new',iostat=status)
      			if(status.ne.0)return
        		INQUIRE (5, NAME=OutFil)
	      	endif


      write(5,1008)'Gibbs input file'
1008  format(A32)
      write(5,1009)Assembtitle
1009  format(A255)
      write(5,1000)TC,PB
1000  format(2f8.1,9x,'Tstart and Pstart')
      WRITE(5,1010)(CoNameInFile(i),i=1,NoSysCoInFile)
1010  format(40(A4,2x))
      WRITE(5,1001)(UseSysCo(i),i=1,NoSysCoInFile)
1001  format(40(I2,4x))
      WRITE(5,1004)imass
1004  format(I6,'              IMass')
      write(5,1011)
1011  format('***************************************')


!     for each mineral write 
!     (1)  Mineral number and name
!     (2)  volume, nfract and IncludeNew variables
!     (3)  the phase component names
!     (4)  UsePhaseComponent switch (o or 1)
!     (5)  Phase Component or cation on sites names
!     (6)  Mole fractions of components or cations on sites
!     (7)  A row of stars

	Do 10 kCur = 1,numPh
	k = asmCurrent(kCur)
      	if(icall.ge.2.and.K.eq.KSkip)go to 11
      	WRITE(5,1002)MINREC(K),PHNAME(K)
1002  	format(I6,21x,a32)
      	WRITE(5,1006)VP2(K),FRACTL(K),IncludeNew(K)
      ! store VP2 here, even though a phase may be fractionating. Then, when we read in the input file
      ! we will check. If the phase is fractionating with ifract = 1 VP0 = 0; VP1 = VP2 = Stored volume
      ! 						  ifract = 2 VP0 = VP1 = 0; VP2 = stored volume
      ! This gives a means of tracking a fractionating phase volume through an assemblage change
      ! (i.e. staurolite-out or xenotime-out)
1006  format(1x,E15.5,2I8,'     VStart,Fractl,IncludeNew')

      write(5,1012)(phCoNameinFile(k,j),j=1,numPhCoInFile(K))
1012  Format(5x,15(A12,3x))
      write(5,1013)(UsePhCo(k,j),j=1,numPhCoInFile(K))
1013  Format(5x,15(I3,12x))

      if(numSiteAtom(K).eq.0)then

!     If here, we must output mole fractions of components
!     note that the compositions stored in array X refer to the adjusted mineral
!     compositions with components not used in this problem removed.
!     However, the input data file must contain all phase components for that mineral
!     We use the pointer in array phCoUsedinFile(15,12) to sort things out

!     do this loop for each phase component
      do 12 L=1,numPhCoInFile(K)
      xx(L)=0.
12    continue

      do 13 J=1,numPhCo(K)
      XX(phCoUsedinFile(k,j)) = sngl(xPhCo(k,j))
13    continue
      Xsymb='X'
      write(5,1014)(Xsymb,phCoNameinFile(k,j),j=1,numPhCoInFile(K))
1014  Format(4x,15(A1,A12,3x))
      write(5,1003)(xx(L),L=1,numPhCoInFile(K))
1003  format(15(F12.9,3x))


      else
!     If here, we must output atoms on sites
      call XtoSite(K)
      write(5,1020)(SiteAtomName(k,L),L = 1,numSiteAtom(K))
1020  format (3x,20A12)
      write(5,1021)(SiteAtom(k,L),L=1,numSiteAtom(K))
1021  format (20F13.9)


      endif

      write(5,1011)          ! write a row of stars

11    continue
!     We must skip the components of all phases, not just the ones we want to keep
!      ico = ico + numPhCoInFile(K)

10    continue

      if(icall.eq.0) close(5)
!     if icall eq.1 or 2 we will close the file in routine CHANGE

      RETURN
      END

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SAVEMIF(TC4,PB4)
      implicit none
!     ROUTINE TO save the input problem as a MIF file
!	This saves both stable and metastable phases using the compositions from the Routine ParallelToTangent
!	It's called from Routine Tangent only after establishing a valid set of PTX


	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Tangent.inc"
! ***********************************
!     local variables
      integer*4 i,j,k,L,IZERO,ivol,status,kCur,iformat
      REAL*4 xx(12),TC4,PB4
      character Xsymb*1

	CALL FSS_ALERT(' ALERT','Give a name for your MIF file')

	OPEN(5,FILE='',STATUS='new',iostat=status)
	if(status.ne.0)return
        INQUIRE (5, NAME=OutFil)

	write(*,*)' Do you want output in fixed format or E format?'
	write(*,*)' Note that with trace elements, fixed format might not give sufficient digits (only 9 decimal places)'
	write(*,*)' 0 = fixed, 1 = E format'
	read(*,*)iformat

      	write(5,1008)'Gibbs3 MIF file'
1008  	format(A15)

	write(Assembtitle,1022)(asmCurrent(k),k=1,numCurrent)
1022	format(' Indices of stable asm: ',20I5)
!	Assembtitle = 'Indicies of stable assemblage'
      	write(5,1009)Assembtitle
1009  	format(A255)
      	write(5,1000)TC4,PB4
1000  	format(2f8.1,9x,'Tstart and Pstart')
      	WRITE(5,1010)(CoNameInFile(i),i=1,NoSysCoInFile)
1010  	format(40(A4,2x))
      	WRITE(5,1001)(UseSysCo(i),i=1,NoSysCoInFile)
1001  	format(40(I2,4x))
      	WRITE(5,1004)imass
1004  	format(I6,'              IMass')
      	write(5,1011)
1011  	format('***************************************')


!     for each mineral write 
!     (1)  Mineral number and name
!     (2)  volume, nfract and IncludeNew variables
!     (3)  the phase component names
!     (4)  UsePhaseComponent switch (o or 1)
!     (5)  Phase Component or cation on sites names
!     (6)  Mole fractions of components or cations on sites
!     (7)  A row of stars

!	Do 10 kCur = 1,numPh
!	k = asmCurrent(kCur)
	Do 10 k = 1,numPhMIF
      	WRITE(5,1002)MINREC(K),PHNAME(K)
1002  	format(I6,21x,a32)
      	WRITE(5,1006)VP2(K),FRACTL(K),IncludeNew(K)
      ! store VP2 here, even though a phase may be fractionating. Then, when we read in the input file
      ! we will check. If the phase is fractionating with ifract = 1 VP0 = 0; VP1 = VP2 = Stored volume
      ! 						  ifract = 2 VP0 = VP1 = 0; VP2 = stored volume
      ! This gives a means of tracking a fractionating phase volume through an assemblage change
      ! (i.e. staurolite-out or xenotime-out)
1006  	format(1x,E15.5,2I8,'     VStart,Fractl,IncludeNew')

      	write(5,1012)(phCoNameinFile(k,j),j=1,numPhCoInFile(K))
1012  	Format(5x,15(A12,3x))
      	write(5,1013)(UsePhCo(k,j),j=1,numPhCoInFile(K))
1013  	Format(5x,15(I3,12x))

      	if(numSiteAtom(K).eq.0)then

!     		If here, we must output mole fractions of components
!     		note that the compositions stored in array X refer to the adjusted mineral
!     		compositions with components not used in this problem removed.
!     		However, the input data file must contain all phase components for that mineral
!     		We use the pointer in array phCoUsedinFile(15,12) to sort things out

!     		do this loop for each phase component
      		do 12 L=1,numPhCoInFile(K)
      		xx(L)=0.
12    		continue

      		do 13 J=1,numPhCo(K)
      		XX(phCoUsedinFile(k,j)) = sngl(xPhCo(k,j))
13    		continue
      		Xsymb='X'
      		write(5,1014)(Xsymb,phCoNameinFile(k,j),j=1,numPhCoInFile(K))
1014  		Format(4x,15(A1,A12,3x))
      		write(5,1003)(xx(L),L=1,numPhCoInFile(K))
1003  		format(15(F12.9,3x))


      		else
!     		If here, we must output atoms on sites
      		call XtoSite(K)
      		write(5,1020)(SiteAtomName(k,L),L = 1,numSiteAtom(K))
1020  		format (3x,20A12)
 		if(iformat.eq.0)then
	      		write(5,1021)(SiteAtom(k,L),L=1,numSiteAtom(K))
			else
	      		write(5,1221)(SiteAtom(k,L),L=1,numSiteAtom(K))
			endif
1021  		format (20F13.9)
1221		format(20E17.9)

      		endif

      	write(5,1011)          ! write a row of stars
11    	continue

10    	continue

      	close(5)

      	RETURN
      	END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine savfil(iunit,iok)
!     subroutine to save a scratch file
      integer*4 iok,iunit,status
      character*32 prompt,FileName
      character*255 readit
!------------------------------------------
	prompt='Type in filename'      
	FileName='Gibbs.output'
	IOK=0
	OPEN(22,FILE='',STATUS='UNKNOWN',iostat=status)
	if(status.ne.0)return
	INQUIRE (22, NAME=FileName)
	rewind(iunit)
	do
	read(iunit,'(A255)',end=923)readit
	write(22,'(A255)')trim(readit)
	repeat
923   	continue
	backspace(iunit)
	close(22)
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine OpenThermoFile

! 	routine to open a new thermodynamic data file
	implicit none
      	integer*4 i
	character*64 fileName
	character*256 aline
! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "HoPoMeltW.inc"
! *****************************************

	fileName = './Gibbs_Essentials/'//trim(thmoFile)		
	close(2)
	close(1)
	OPEN (2,FILE=fileName,STATUS='OLD')
	OPEN (1,FILE='ThermoTempFile',STATUS='UNKNOWN')

!	Strip all comment lines from datafile
10	continue
	read(2,11,end=20)aline
11	format(A256)
	if(aline(1:1).eq.'!')then
		go to 10
		else
		write(1,11)aline
		endif
	go to 10
20	continue

	rewind(1)

	i=0
	call setthermofile(i)
	if(i.eq.1)then
		write(*,*) ' There is a problem reading the thermodynamic data file.'
		pause ' hit return to continue'
		return
		endif
			
! 	Load "fixed" minerals from thermodynamic data file
! 	1 = alpha quartz, 2 = beta quartz, 3=kyanite,4=sill,5=andalusite
      	call LoadFixed

! 	Open and read HP98 Excess term file
!	This is replaced with generalized SF input
!	open(88,file = './Gibbs_Essentials/HP98MeltWterms',status='OLD')
!	write(*,*)'Reading HP98 Melt excess terms'
!	do 205 i = 1,11		! read file header
!	read(88,*)dummy
!	write(*,*)dummy
!205	continue	
!	do 215 i = 1,30
!	do 215 j = 1,30
!	HoPoWH(i,j) = 0.0d0
!	HoPoWS(i,j) = 0.0d0
!	HoPoWV(i,j) = 0.0d0
!	HoPoWG(i,j) = 0.0d0
!215	continue
!	do 210 i = 1,8
!	read(88,*)dummy 	! read line of dashes
!	do 220 j = i+1,8
!	read(88,*)HoPoWH(i,j),HoPoWS(i,j),HoPoWV(i,j)
!!	write(*,*)i,j,HoPoWH(i,j),HoPoWS(i,j),HoPoWV(i,j)
!220	continue
!210	continue
!	close(88)
	return
	end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine setthermofile(iok)
!    	use AWE_Interfaces
	implicit none
! *****************************************
	include "Assemb.inc"
	include "Newton.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
! *****************************************

      character tempname*32,stuff*64,stars*5,first5*5,dummy*16,dollars*5
      integer*4 linenumber,nummineral,k,iok,i,j
      data stars /'*****'/
      data dollars /'$$$$$'/


!     Read component names and molecular wts from thermo data file
	READ(1,*)NoSysCoInFile
	write(*,*)' Reading component names, mol wts, cations in oxides, oxide to element weights'
      	read(1,*)(CONameInFile(J),J=1,NoSysCoInFile)
      	read(1,*)(MOLWTInFile(J),J=1,NoSysCoInFile)
      	read(1,*)(NumCatInOxideInFile(J),J=1,NoSysCoInFile)
      	read(1,*)(OxToElWtInFile(J),J=1,NoSysCoInFile)


	write(*,*)'Reading component plotting variables'
	read(1,*)dummy				! row of dashes
!	read(1,*)dummy				! label "predefined components"
	read(1,*)NoCompPlotVars
!	read(1,*)dummy				! row of component labels
	do 2 i = 1,NoCompPlotVars
	read(1,*)CompPlotName(i),CompPlotConst(i),(CompPlotCoeffInFile(i,j),j=1,NoSysCoInFile)
2	continue

! 	Read in information about fixed mineral entries
! 	The thermodynamic data for these minerals is read in from subroutine LoadFixed
	write(*,*)'Reading fixed phases'
	read(1,*)dummy				! row of dashes
	read(1,*)fixedNo			! number of "fixed" minerals to read
	if(fixedNo.eq.0)go to 8
!	read(1,*)dummy				! header line
	j = 0					! pointer to start of the mineral in fixedPhase array
	do 20 i = 1,fixedNo
	read(1,*)fixedID(i),fixedIDNoCo(i),(fixedPhase(k+j),k=1,fixedIDNoCo(i))
	j = j+fixedIDNoCo(i)
20	continue


! 	Now find the beginning of the thermodynamic data
! 	We are looking for a row of $$$$$$$$$$$$$$$$$$$$$$$$
8	continue
	rewind(1)				! start over at the beginning to make counting lines easier
	linenumber = 0
	do
	read(1,'(A64)',end=999)stuff
	linenumber = linenumber + 1
	first5=stuff(1:5)
	if(first5.eq.dollars)then        		! we have found the first row of stars in the data file
            go to 9
	endif
	repeat
	
9     continue
!     Read positions of minerals in data files
      write(*,*)'Reading phase headers'
      nummineral=0
10    continue
      linenumber=linenumber+1
      read(1,'(A64)',end=999)stuff
      first5=stuff(1:5)
      if(first5.eq.stars)then        		! we have found the beginning of a new mineral
		nummineral=nummineral+1
		linenumber=linenumber+1
		read(1,'(A64)',end=999)stuff
		write(*,'(A64)')stuff
		read(stuff,*,end=999)K,tempname
		mineralnumber(nummineral)=K
		mineralposition(K)=linenumber
		mineralname=tempname
		write(*,206)nummineral,K,mineralname,mineralposition(K)
206   		format(i5,i5,2x,a32,i7)
            	endif
      go to 10

999   continue
      mineralnumber(nummineral+1)=150
      iok=0
      return
      end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine LoadFixed
      implicit none      
      integer*4 K,j,KK,NR,newline,TnoPhco,iPh,iCo,i,L,TPhaseType
      character*8 Tphname(12)
!      integer*4 imin(5)
! *****************************************
	include "Assemb.inc"
 	include "Output.inc"
! *****************************************
!     load "fixed" minerals from thermodyanmic data file

! 	integer*4 fixedNo,fixedID(fnMax),fixedIDNoCo(fcMax),fixedPhase(fcMax)
! 	fixedNo = number of phases that will have "fixed" status
! 	fixedID = phase number of the fixed phase (use in input file to select this phase)
! 	fixedIDNoCo = number of "components" in this fixed phase
! 	fixedPhase = numbers from the data file that are the "component" choices of this "fixed" phase
! 	e.g. KySilAnd would have 5,6 207 corresponding to Ky, Sil and And in the data file

!     fixed minerals in SPaC are
!      imin(1) = 101           !alpha quartz
!      imin(2) = 102           !beta quartz
!      imin(3) = 5             !kyanite
!      imin(4) = 6             !sillimanite
!      imin(5) = 7             !andalusite
 
 	write(*,*) ' Loading fixed minerals'
	rewind(1)
	K = 0
	do 10 iPh = 1,fixedNo
	write(*,*)
	do 10 iCo = 1,fixedIDNoCo(iPh)
	K = K+1
	NR = fixedPhase(K)      
	newline = mineralposition(NR) ! position of mineral in datafile

        rewind(1)
	do 426 i=1,newline-1		! skip up to where mineral NR is located
        read(1,*)
426     continue
        READ (1,*)KK,Tphname(k)
	write(*,*)kk,Tphname(k)
      	if(kk.ne.NR)then
            write(*,*)' ROUTINE LoadFixed'
            write(*,*)' There is an error reading the thermodynamic data file'
            write(*,*)' Offending record:'
            Write(*,*)' Tried to find, Actually read, Phase name read'
            write(*,134)NR,KK,Tphname(k)
134         format(2I5,5x,8a)
            Write(*,*) 'Hit return to abort'
            pause
            stop
            endif

      	READ(1,*)TNoPHCO,TphaseType
!     	do this loop for each phase component
      	do 20 j=1,Tnophco
      	READ(1,'(5X,A16)')
      	read(1,*)HZEROF(K),SzeroF(K),VzeroF(K)
      	READ(1,*)ACPF(K),BCPF(K),CCPF(K),DCPF(K),ECPF(K),FCPF(K),GCPF(K)
      	READ(1,*)v1F(K),v2F(K),v3F(K),v4F(K)

	select case(TphaseType)
	case(1)				! SPaC (no landau terms)
	case(101)				! HP98 landau terms
		read(1,*)LandauTcF(k),LandauSmaxF(k),LandauVmaxF(k)
	case(201,208)				! HP11 Model terms
		read(1,*)HP11ModelSwitchF(k),(HP11ModelF(k,L),L=1,6)
		read(1,*)DQFswitchF(k),DQFHF(k),DQFSF(k),DQFVF(k)
!		write(*,*)DQFswitch(k,jjj),DQFH(k,jjj),DQFS(k,jjj),DQFV(k,jjj),'   ! DQF switch, H, S, V'
	case default
		write(*,*)' Something went wrong reading the fixed minerals'
		write(*,*)' The phase type should be 1, 101 or 201 (a 1 component SPaC, HP98, or HP11 phase)'
		pause ' Hit return to abort'
		stop
	end select

20    	continue
10    	continue
       write(*,*)' Exit load fixed minerals'
      return
      end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine ReadThomoFileList(thisThmoFile)
! 	routine to pick a mineral assemblage from a Master Input File (MIF)
	implicit none
      	integer*4 i,j
	character*64 thisThmoFile,thmoFileList(20)
	character*1 dummy
! *****************************************
	rewind(3)
      	READ(3,'(A)')dummy          	! read title
      	READ(3,'(A)')dummy          	! read the row of stars
	read(3,*)i 			! Number of thermo files in list
	do 315 j = 1,i
      	READ(3,'(a)')thmoFileList(j)
315	continue
	write(*,*)' Thermodynamic data files in Gibbs.fig file:'
	write(*,*)' '
	do 320 j = 1,i
	write(*,*)j,' ',thmoFileList(j)
320	continue
	write(*,*)' Pick thermodynamic data file to use in this session'
	read(*,*)j
	thisThmoFile = thmoFileList(j)
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PickBulkComp
      implicit none
! Routine to pick a bulk rock composition from the file "BulkRockAnalyses" in Gibbs_Essentials folder
! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
!	include "PlotStuff.inc"
!	include "PlotGibbs.inc"
! *****************************************
!      	local variables
      	INTEGER i,j,iok
      	CHARACTER TempTitle*80,dummy*4
	Character*80 TheList(200),WhatToPick
	Integer*4 TheOne,ListSize,iselect
! -----------------------------------------------

!      	call volget(voltemp)			! Save current folder reference
!      	call volset(volrefnum)			! set to Gibbs  (where program started)
!	open(unit=24,file="./Gibbs_Essentials/BulkRockAnalyses",status="OLD")
	open(unit=24,file=BulkRockCompositionFile,status="OLD")

!      	Build a list of analysis titles for listmanager
      	READ(24,'(A)')dummy			! read the header
	j = 0
     	DO
	  READ(24,'(A)',END=4170)dummy          ! read the row of dashes
        READ(24,'(A)',END=4170)TempTitle
	  j = j+1
	if(j.gt.200)go to 4170			! maximum dimensions of TheList(200)
        TheList(j) = Trim(TempTitle)
        DO 415 I=1,3
            READ(24,'(A)',END=4170)dummy	! read 3 lines
415     CONTINUE
	REPEAT

4170     CONTINUE

!     Put the list up and let the user pick one plot type
      iSelect = 1
      ListSize = j
      TheOne = 1		! always select the first bulk composition in the list
      WhatToPick = 'Select a bulk composition '

!	Sub  AWE_Pick_one(TheOne,ListLength,TheList,WhatToPick,iok,iselect)
      	call AWE_Pick_One(TheOne,ListSize,TheList,WhatToPick,iok,iSelect)
      	if(iok.ne.0)go to 999    			! user hit cancel

	j = TheOne
      REWIND(24)
      READ(24,'(A)')dummy			! read the header


      DO 419 i=1,(j-1)*5
      READ(24,'(A)')dummy     		! ignore the unwanted plot definitions
419   CONTINUE

      READ(24,'(A)')dummy     		! Read line of dashes
      READ(24,'(A)')bulkCompTitle
	read(24,*)nc				! number of components
      READ(24,*)dummy				! names of oxides
	read(24,*)(bulkCompWt(i),i=1,nc)	! bulk composition in wt %

! Bulk composition is converted to moles in subroutine Begin where input data are read
	
999	continue
	close(24)
!     	call volset(voltemp)			! set folder back to where we started

	bulkCompSwitch = .true.			! a bulk composition is now open
	return
	end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      	SUBROUTINE SetALLX()
      	implicit none
!	Routine to set the array ALLX based on the current mineral assemblage
!	Called in Begin (after opening a new file)
!	and from CHANGE after the assemblage is changed
!	Will also be called from TANGENT after assemblage is changed      
	include "Assemb.inc"
	include "Monit.inc"

      	integer*4 i,j,k,L,kCur

!        ARRAY ALLX CONTAINS VALUES OF THE FOLLOWING
!     set all columns of ALLX(i,j) to starting values
!     Definition of ALLX variable:
!     1,50   current value
!     2,50     starting value
!     3,50   ref at start of contour
!     4,50   user selected reference
!     5,50   previous finite diff point
!     6,50     (not used)


      	DO 10 I=1,50
       	allFractl(i) = 0				! flag keyed to ALLX to tell whether a phase is fractionating - used in COMPUT2
        DO 10 J=1,6
        ALLX(J,I)=0.D0
10  	CONTINUE

      	DO 444 I=1,6
      	ALLX(I,1)=TSTART
      	ALLX(I,2)=PSTART
      	IF (KFLU.NE.0)ALLX(I,3)=PFStart
!     	THESE ARE THE X'S
      	L=TANDP
	Do 441 kCur = 1,numPh
	k = asmCurrent(kCur)
      	DO 555 J=2,numPhCo(K)
      	L=L+1
      	ALLX(I,L)=xPhCo(k,j)
555   	CONTINUE
441   CONTINUE
!     THESE ARE MOLES
      	L=TANDP+NX
      	IF(IMASS.EQ.1) THEN
		Do 442 kCur = 1,numPh
		k = asmCurrent(kCur)
            	L=L+1
            	ALLX(I,L)=MP0(K)
            	if(fractl(k).eq.1.or.fractl(k).eq.2)then
            		!This is done as a flag for COMPUT2 where it checks for Newton convergence (around line 722)
            		! If FRACTL = 1 or 2 then NEWTON doesn't converge properly.
            		! This flag will omit moles of fractionating phases from Newton convergence
            		allx(6,L)=1
            		else
            		allx(6,L)=0
            		endif
442         	CONTINUE
        	ENDIF
!     Set allFractl switch if necessary
      IF(IMASS.EQ.1) THEN
      		L=TANDP+NX
		Do 443 kCur = 1,numPh
		k = asmCurrent(kCur)
            	L=L+1
            	if(fractl(k).eq.1.or.fractl(k).eq.2)then
            		!This is done as a flag for COMPUT2 where it checks for Newton convergence (around line 722)
            		! If FRACTL = 1 or 2 then NEWTON doesn't converge properly.
            		! This flag will omit moles of fractionating phases from Newton convergence
            		allFractl(L)=1
            		else
            		allFractl(L)=0
            		endif
443         	CONTINUE
        	ENDIF
444   CONTINUE
	call SetNewVar()
      	CALL NAMES
      	NVAR= NVAR + numNew


      	RETURN
      	END


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine SetNewVar()
      	implicit none
!	Routine to set the new variables into the array ALLX
	include "Assemb.inc"
	include "Monit.inc"

      	integer*4 i,j,jj,k,L,LL,kCur

	LL = TANDP + NX
      	if(IMASS.EQ.1)LL = LL + numPh

	numNew = 0
	!if(newVarGlobal.eq.1)then        		! set up for using all new variables
		jj = 0		
		do 5701 kCur = 1,numPh
		k = asmCurrent(kCur)
		if(includeNew(k).eq.1)then
			numNew = numNew + numNewInPhK(k)
		      	do 5702 i = 1,numNewInPhK(k)
!     			Set ALLX variable
!     			0 = New - (A + aj*Xj)/(B + bj*Xj)
      			!compute numerator
		      	numerator = NewNumA(k,i)
      			do 5751 j = 1,numPhCo(k)
5751  			numerator = numerator + newNumaj(k,i,j)*xPhCo(k,j)
      			!compute denominator
      			denominator = newDenB(k,i)
      			do 5752 j = 1,numPhCo(k)
5752  			denominator = denominator + newDenBj(k,i,j)*xPhCo(k,j)
			jj = jj + 1
      			DO 5750 L=1,5
      			ALLX(L,LL+jj) = numerator/denominator
5750  			continue          
5702  			continue
			endif
5701  		continue
	!	endif
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine CheckSystemComponents(iok)
!	Routine to check if the system components in the bulk composition file are the same as those
!		in the MIF
!	The order of components in the MIF must be the same as the order of elements in the bulkcomposition file
	include "Assemb.inc"
	integer*4 i
!
	write(*,*)(coname(i),i=1,nc)
	write(*,*)(BCconame(i),i=1,nc)
	pause' Take a look'
	iok = 0
	do 10 i = 1,nc
!	do 15 j = 1,nc	
	if(trim(adjustL(BCcoName(i))).ne.trim(adjustL(coName(i))))then
!15	continue
		write(*,*)' System components in the bulk composition file do not match those'
		write(*,*)'    in the MIF'
		write(*,*)'    Please fix'
		iok = 1
		return
		endif
10	continue
	iok = 0
	return
	end
	
	
