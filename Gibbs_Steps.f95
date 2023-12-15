! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE STEPS
    	use AWE_Interfaces
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
	include "TPSave.inc"
! *****************************************
! --------------------------------------------------------------------
	include "Tangent.inc"
! --------------------------------------------------------------------
!      local variables
	integer*4 logfileon,nvar2
      integer*4 ifdplt,i,isteps,j,L,K,istep,izero,irefer,jsymb,isymsize,NeedToPick,iok,myColor
      real*4 Xp(4)
	integer*4 Picked(15)
	real*4 ValuesPicked(15)

      if(noOxygenMassBalance.eq.1)then
		NEQ = NEQ - 1		! this is reset just before exiting the routine at 9999
		write(*,*)' NEQ RESET = ',NEQ
	      endif

	TPCounter = 1
	TPSave(TPCounter,1) = ALLX(1,1)
	TPSave(TPCounter,2) = ALLX(1,2)

	iDoingGrid = 0			! flag to alert activity calculations to stop if finite difference derivative fails at a subsystem
						! Only used in Subroutine dlnAdX in file Gibbs_Thermocalc.f95
      Plotcounter=0
      ifdplt=0
1     CONTINUE
  !    	call FSS_SetPort(1)
		!call !scolor(1)
      WRITE(*,*)' *********************************'
      write(*,461)TC,PB
	call Variance			! write out variance
      if(noOxygenMassBalance.eq.1)then
		write(*,*)' No Oxygen Mass Balance'
	      endif
      WRITE(*,*)' CURRENT MONITOR PARAMETERS ARE'
      if(nvar-neq.gt.0)then
      do 6 i = 1,nvar-neq
      j = smon(i)
      if(VN1(j)(1:2).eq."M_")then
		WRITE(*,5)j,VN1(j),VN2(j),ALLX(1,j)*1000.,SDEL(I)*SNSTEP*1000.
		else
		WRITE(*,5)j,VN1(j),VN2(j),ALLX(1,j),SDEL(I)*SNSTEP
		endif
5     FORMAT(I3,1X,A2,A6,2F13.4)
6     continue
      else
      write(*,*)' No monitor parameters - variance <= 0'
      endif
      WRITE(*,*)'NSTEP =',SNSTEP
      if(idraw.eq.0)write(*,*)'IDRAW = 0 (do not draw the line)'
      if(idraw.eq.1)write(*,*)'IDRAW = 1 (draw the line)'


! ----------------------------------------------------
	if(iternary.eq.0)then			! X-Y plot
		write(*,*)'Current X-Y axes are:'
		write(*,*)'       X             Y'
		do 9 k = 1,NoTieMinerals
		i = PlotIndex(k,1)
		j = PlotIndex(k,2)
9		write(*,7)i,VN1(i),VN2(i),j,VN1(j),VN2(j)
! 9		write(*,7)PlotIndex(k,1),VN1(PlotIndex(k,1)),VN2(PlotIndex(k,1)),PlotIndex(k,2),VN1(PlotIndex(k,2)),VN2(PlotIndex(k,2))
7		FORMAT(4(I4,' ',A2,A8))

	else   					! ternary plot
	      	write(*,*)'Coordinates are:'
	      	write(*,12)(CompPlotName(Apx(i)),i=1,iternary+2)
12		Format(4(A8,8x))
		do 13 k = 1,NoTieMinerals
13		write(*,14)PlotIndex(k,1),PhName(PlotIndex(k,1))
14		FORMAT(4(I4,' ',A12))

	endif
!-----------------------------------------------



      WRITE(*,*)' *********************************'
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
      WRITE(*,*)' *********************************'
      WRITE(*,*)' SINGLE STEPS MENU OPTIONS'
      WRITE(*,*)'    0 = Return'
      WRITE(*,*)'    1 = Line or tie-lines (don''t touch unless you want to change plotting)'
      WRITE(*,*)'    2 = Choose monitors/set deltas'
      WRITE(*,*)'    3'     ! = Set deltas'
      WRITE(*,*)'    4 = Compute one increment'
      WRITE(*,*)'    5 = Unstep one finite difference'
      WRITE(*,*)'    6 = Reference points'
      WRITE(*,*)'    7 = Print assemblage information to Output window'
      WRITE(*,*)'    8 = Go to global menu'
      write(*,*)'    9 = Go to Plotter options'
      write(*,*)'   10 = Gibbs/Newton switch'
      write(*,*)'   11 = Save a digitized reaction'
      write(*,*)'   20 = Begin/Save menu'
      write(*,*)'---------------------------------'
      write(*,*)'   -1 = Change mineral assemblage'
      write(*,*)'   -2 = Set output length'
      write(*,*)'   -3 = NewtonStepsOutput switch'
      write(*,*)'   -4 = Pick line color'
      WRITE(*,*)' CHOOSE OPTION'
      read(*,*)isteps
      

! 	if(isteps.lt.0)call global(isteps)
! ---------------------------------------
	if(isteps.eq.-1)then
!      		if we are changing the assemblage, we must zero out any "new" variables
!      		because there is no way to determine whether they will still be in the 
!      		list after a mineral is removed.
      		NVAR=NVAR-numNew
      		NEQ=NEQ-numNew
      		numNew=0
       		call change
!       	call rexn
		go to 1
		endif
! ---------------------------------------
	if(isteps.eq.-2)then
		isteps = 0
		call AWE_setoutput
		goto 1
		endif
! ---------------------------------------
	if(isteps.eq.-3)then
		write(*,*)'NewtonStepsOutput = ',newtonStepsOutput
		write(*,*)'Input new value (0 or 1)'
		read(*,*)newtonStepsOutput
		if(newtonStepsOutput.eq.0)close(43)
		if(newtonStepsOutput.eq.1)then
			open(43,file='',status='UNKNOWN')
			endif
		goto 1
		endif
! ---------------------------------------
	if(isteps.eq.-4)then
!		myColor = 1		! Keep the same myColor as previously
		call PickAWEColor(myColor)
		CurrentColor = myColor				! CurrentColor is in PlotStuff.inc common block
		goto 1
		endif

	
! ------Begin/Save routine-------------
	if(isteps.eq.20)then
		call GibbsBegin
		go to 1
	endif

	
      GO TO (50,100,200,300,400,500,600,700,800,900,1000,1100)ISTEPS+1
      go to 1

50    CONTINUE
      if(noOxygenMassBalance.eq.1)then
		NEQ = NEQ + 1		! this is reset just before exiting the routine at 9999
		write(*,*)' NEQ RESET = ',NEQ
	      endif

      RETURN
!  ----Reset reaction counter------------------------------------------
1100   CONTINUE
      write(*,*)' 0 = return' 
      write(*,*)' 1 = Reset counter for saving reaction'
      write(*,*)' 2 = Save this reaction (since counter reset)'
      write(*,*)' 3 = List contents of reaction storage array'
	read(*,*)i
	if(i.eq.0)go to 1
	if(i.eq.1)then
		!save current T,P conditions
		TPcounter=1
		TPSave(TPCounter,1) = ALLX(1,1)
		TPSave(TPCounter,2) = ALLX(1,2)
		go to 1
		endif
		
	if(i.eq.2)then
		call savrxn
		go to 1
		endif
	if(i.eq.3)then
		write(*,*)'i        T         P'
		do 1105 i = 1,TPcounter
		write(*,*)i,TPSave(i,1),TPSave(i,2)
1105		continue
		go to 1100
		endif

      GO TO 1100
!  ----inewton------------------------------------------
1000   CONTINUE

	call GibbsorNewton

      GO TO 1
!  ----Set plot stuff------------------------------------------
100   CONTINUE
!      call listsy
      read(*,*)isymb
      isymsize=2
      if(isymb.gt.0)then
            write(*,*)' Input symbol size in pixels (1-10)'
            write(*,*)' (1 is 2x2 pixels, 2 is 4x4 pixels etc.)'
            read(*,*)isymsize
            if(isymsize.le.0)isymsize=1
            endif
      write(*,*)' Do you want to draw lines between points (IDRAW)?'
      write(*,*)' 0 = Do not draw lines'
      write(*,*)' 1 = Draw lines'
      read(*,*)idraw
      GO TO 1


!  ----CHOOSE MONITORS/set deltas-------------------------------
200   CONTINUE

        NeedToPick = Nvar-Neq
	do i = 1,NeedToPick
		Picked(i) = SMon(i)
		ValuesPicked(i) = SDel(i)*Dfloat(SNstep)
		end do
	nvar2 = nvar
	iok = 0
	call PickMonitors(nvar2,NeedToPick,Picked,ValuesPicked,iok)
!	pause 'Back from PickMonitors -- hit return'
        if(iok.eq.0)then
!        SNstep = NStep        !in common block Monitor
	        if(SNStep.le.0)SNStep = 1
	        do i = 1,NeedToPick
			SMon(i) = Picked(i)
			SDel(i) = ValuesPicked(i)/DFloat(SNstep)
! 			write(*,*)i,Smon(i),Sdel(i)
			end do
	        endif

      GO TO 1
!  -----SET DELTAS----------------------------------
300   CONTINUE
	! This routine now done in subroutine Pkmon
      GO TO 1
! ----COMPUTE ONE INCREMENT------------------------
400   CONTINUE

      if(noOxygenMassBalance.eq.1)then
		NEQ = NEQ + 1		! this is reset just before exiting the routine at 9999
		write(*,*)' NEQ RESET = ',NEQ
	      endif


      jsymb=isymb
      if(isymb.eq.0)jsymb=3
!      set monitor parameters
      NSTEP=SNSTEP
      do 401 i=1,nvar-neq
      deltax(i)=sdel(i)
401   mon(i)=smon(i)
!      Set up array IPOINT to contain pointers to non-monitor parameters
      J=0
      Do 410 i=1,NVAR
      DO 411 L=1,NVAR-NEQ
      IF(MON(L).eq.i)go to 410
411   continue
      J=J+1
      IPOINT(J) = I
410   CONTINUE
      IF(iLong(4).eq.1) write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)

!      first draw the starting point on the screen if idraw=1
      Plotcounter=0
!	!call !scolor(icolor)
	call stepplot(1)
	!call !scolor(1)

!      THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
      DO 403 I=1,nvar
!         ALLX(3,I) IS WHERE calculations STARTED-save
         ALLX(3,I)=ALLX(1,I)
403   CONTINUE

      do 499 istep=1,snstep

      if(inewton.eq.1)then
!         Store new values for each independent variable (monitor parameter)
! 	This is where Neton's method increments composition, T or P
         DO 420 J=1,NVAR-NEQ
         i=mon(j)
         ALLX(1,i)=ALLX(1,i) + DELTAX(J)
420      continue
         call SetTPX
         endif
         
      write(*,461)TC,PB
461   format (' T = ',f10.1,'  P = ',f10.1)
      write(*,*)' NSTEP = ',snstep,' Counter = ',istep

        IZERO=0
        call c2ompute(izero)
        IF(IZERO.gt.0)GO TO 470


      if(noOxygenMassBalance.eq.1)then
		NEQ = NEQ - 1		! this is reset just before exiting the routine at 9999
		write(*,*)' NEQ RESET = ',NEQ
	      endif



!      plot this finite difference iteration
	call stepplot(2)

! 	Check for mass balance error
	if(imass.eq.1)then
		call CheckMassBal
		endif

! 	save reaction data for possible output to file
	TPCounter = TPCounter+1
	if(TPCounter.gt.1000)then
		call fss_alert('ALERT!!','TPCounter greater than 1000.  Resetting to 1')
		TPCounter = 1
		endif
	TPSave(TPCounter,1) = ALLX(1,1)
	TPSave(TPCounter,2) = ALLX(1,2)

499   continue

470	continue

	PenSave(PlotCounter+1) = 0	! Just to force a penup in postscript routine
! 	Now replot in the Picture
	do 480 i = 1,PlotCounter
	do 481 j = 1,iternary+2
481	Xp(j) = PlotSave(i,j)
	ipen = PenSave(i)
!	Call GPlot (Xp,ipen,0)
	if(jps.eq.1) then
		call Gplot(Xp,ipen,1)			! 1 means plot in the illustrator file
		if(PenSave(i+1).eq.0)Call Ppenup
		else
		Call GPlot (Xp,ipen,0)			! 0 means no illustrator plot
		endif
480	continue

!      now draw symbol at the endpoint
      if(isymb.gt.0)then
!	!call !scolor(icolor)
         call Gsymb(Xp,isymb,isymsize,0)  ! plot a symbol if needed
!         if(jps.eq.1.and.idraw.eq.1)call Gsymb(Xp,isymb,isymsize,1)
	!call !scolor(1)
         endif

!      done with this set of calculations
      !call !scolor(1)
      go to 1
!  ---UNSTEP---------------------------------------
500   CONTINUE
         CALL RESET(5,1)
         GO TO 1

!  ---REFERENCE POINTS--------------------------
600     CONTINUE
       WRITE(*,*)' REFERENCE POINT OPTIONS'
       WRITE(*,*)'   0 = Return  '
       WRITE(*,*)'   1 = Set reference point'
       WRITE(*,*)'   2 = Return to reference'
       WRITE(*,*)'   3 = Return to starting conditions'
       READ(*,*)IREFER
!      Definition of ALLX variable:
!      1,50   current value
!      2,50     starting value
!      3,50   ref at start of contour
!      4,50   user selected reference
!      5,50   previous finite diff point
!      6,50     (not used)
       IF(IREFER.EQ.1) THEN
        DO 610 I=1,nvar
           ALLX(4,I)=ALLX(1,I)
610      CONTINUE
       ENDIF
       IF(IREFER.EQ.2) CALL RESET(4,1)
       IF(IREFER.EQ.3) CALL RESET(2,1)
       GO TO 1
!  ----PRINT------------------------------------------
700     CONTINUE
         CALL PRINTT(1)
         GO TO 1
!  ----GLOBAL------------------------------------------
800     CONTINUE
         CALL GLOBAL(isteps)
         GO TO 1
! ----------plotter stuff------------------------------
900   continue
        CALL PLOTIN
        go to 1


      END

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine stepplot(ifirst)
    	use AWE_Interfaces

      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
! *****************************************
!      local variables
      integer*4 i,ii,ifirst,iax,j,K,K1,K2
      real*4 Plot1(4),Plot2(4),onethird,pi
	data onethird,pi/0.33333333333333333,3.14159/
! ------------------------------------------


	call phasecomp						! calculate composition of each mineral

! 	Plot single minerals separately from tie line minerals

	if(NoTieMinerals.eq.1)then				! Plot single points

	do 5 iax = 1,iternary+2					! iax loops on axis number

	if(iternary.eq.0)then					! Plotting X-Y diagram
		select case (axisIsPhase(iax,1))
		case(0)						! plot as ALLX variable (whatever it may be)
		Plot1(iax) = ALLX(1,PlotIndex(1,iax))	
		if(PlotIndex(1,iax).eq.2)Plot1(iax)=Plot1(iax)/1000.	! Plot pressure as kbar
		case(1)						! plot as moles of phase (MP2)
		Plot1(iax) = MP2(axisIsPhase(iax,2))
		case(2)						! plot as volume of phase (VP2)
		Plot1(iax) = VP2(axisIsPhase(iax,2))
		case(3)						! plot as radius of phase (VP2)
		Plot1(iax) = ((3./(pi*4.))*VP2(axisIsPhase(iax,2)))**onethird
		case default
		end select

	else								! plotting ternary or tetrahedron
		K = PlotIndex(1,1)
		Plot1(iax) = CompPlotConst(Apx(iax))
		do 7 j = 1,NC
		Plot1(iax) = Plot1(iax) + PhComp(K,j)*CompPlotCoeff(Apx(iax),j)
7		continue

	endif

5	continue

	CALL GSymb(Plot1,3,3,0)					! draw a symbol

        if(idraw.eq.1)then					! draw the line
!     	 	Write X and Y to save arrays
	    	Plotcounter=Plotcounter+1
	    	if(Plotcounter.gt.plotCounterMax)then
	       		call FSS_Alert('ALERT!!','Too many points: redimension SaveRxn')
	       		Plotcounter=1
	       		return
	       		endif
		do 8 iax = 1,iternary+2
8		PlotSave(Plotcounter,iax)=Plot1(iax)

	   	if(ifirst.eq.1)then
	   		PenSave(PlotCounter) = 0		! Do not draw line if this is before calculation
	   		else
	   		PenSave(PlotCounter) = 1		! Draw line, if required
	   		endif
	   
           	endif              				! end idraw = 1

! 	Plot tie lines
	else
	DO 100 i=1,NoTieMinerals-1
        	DO 100 ii=i+1,NoTieMinerals
			do 120 iax = 1,iternary+2		
			
			if(iternary.eq.0)then			! Plotting X-Y diagram
				Plot1(iax) = ALLX(1,PlotIndex(i,iax))	
				if(PlotIndex(i,iax).eq.2)Plot1(iax)=Plot1(iax)/1000.	! Plot pressure as kbar
				Plot2(iax) = ALLX(1,PlotIndex(ii,iax))	
				if(PlotIndex(ii,iax).eq.2)Plot2(iax)=Plot2(iax)/1000.	! Plot pressure as kbar

				else					! plotting ternary or tetrahedron
				K1 = PlotIndex(i,1)
				K2 = PlotIndex(ii,1)
				Plot1(iax) = CompPlotConst(Apx(iax))
				Plot2(iax) = CompPlotConst(Apx(iax))
				do 130 j = 1,NC
				Plot1(iax) = Plot1(iax) + PhComp(K1,j)*CompPlotCoeff(Apx(iax),j)
				Plot2(iax) = Plot2(iax) + PhComp(K2,j)*CompPlotCoeff(Apx(iax),j)
130				continue
				endif	

120			continue

! 		Plot the tie line on the screen
        	CALL GPlot(Plot1,0,0)
		CALL GPlot(Plot2,1,0)

! 		Plot a symbol at each point
!              	CALL GSymb(Plot2,3,3,0)
! 	     	CALL GSymb(Plot1,3,3,0)

!      		Write X and Y to save arrays
	    	Plotcounter=Plotcounter+1
	    	if(Plotcounter+1.gt.plotCounterMax)then
	       		call FSS_Alert('ALERT!!','Too many points: redimension SaveRxn')
	       		Plotcounter=1
	       		return
	       		endif

		do 140 iax = 1,iternary+2	
140	    	Plotsave(Plotcounter,iax) = Plot1(iax)
	    	PenSave(PlotCounter) = 0

	    	Plotcounter=Plotcounter+1
		do 142 iax = 1,iternary+2	
142	    	Plotsave(Plotcounter,iax) = Plot2(iax)
	    	PenSave(PlotCounter) = 1

100		continue

	endif					! end NoTieMinerals 


	return
	end

