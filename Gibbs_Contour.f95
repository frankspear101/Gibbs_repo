! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE CONTOUR
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
	include "Singular.inc"
! *****************************************
!      local variables
      INTEGER*4 i,jopt,irefer,NeedToPick,Picked(10),iok,myColor
!      real*4 zinc4
      real*4 ValuesPicked(10)

!       data ncont/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20/
! --------------------------------------------------------

!   	set some defaults for Contour routine
      	idraw=1
       	ISYMB=0
!       DEFAULT PEN # and symbol
!      	Default increment size
!      	ZINC4=1
!	zinc = zinc4
	zinc = 1.0d0
      	XINCST=(xmax-xmin)*zinc/100.D0
      	YINCST=(ymax-ymin)*zinc/100.D0
!       PROGRAM USES BARS
!       adjust because pressure axis IS  KBARS
      	IF(NYPLT.EQ.2)YINCST=YINCST*1000.D0
!      	set the variable convert, which is the tolerance for angles in routine compute
!      	convert = (5 degrees)*zinc expressed in radians
!            (5 degrees = .087266 radians)
!      	this way if zinc = 2 the tolerance is 10 deg, if zinc = 3 tol = 15 etc. 
      	convert = .087266d0 * zinc

!        WRITE(*,*)' THIS OPTION PLOTS ISOPLETHS OF CONSTANT VALUE'
!        WRITE(*,*)'     OF A VARIABLE ON A GENERAL X-Y DIAGRGAM'
!        WRITE(*,*)' '
!  ---------------------------------------
1     	CONTINUE
   !   	call FSS_SetPort(1)
      	!call scolor(1)
      	WRITE(*,*)' ********************************'
      	write(*,461)TC,PB
461   format (' T = ',f10.1,'  P = ',f10.1)
	call Variance			! write out variance
      	WRITE(*,*)' CURRENT CONSTANT VARIABLES ARE'
      	do 6 I = 2,NVAR-NEQ
      if(VN1(NCONT(I))(1:2).eq."M_")then
      	WRITE(*,5)NCONT(I),VN1(NCONT(I)),VN2(NCONT(I)),ALLX(1,NCONT(I))*1000.
		else
      	WRITE(*,5)NCONT(I),VN1(NCONT(I)),VN2(NCONT(I)),ALLX(1,NCONT(I))
		endif
!      	WRITE(*,5)NCONT(I),VN1(NCONT(I)),VN2(NCONT(I)),ALLX(1,NCONT(I))
5     	FORMAT(' ',I3,2X,A2,A16,F12.4)
6     	continue
      	WRITE(*,7)NXPLT,VN1(NXPLT),VN2(NXPLT)
7     	FORMAT('Current X-Y axes are: X = ',I4,' ',A2,A16)
      	WRITE(*,8)NYPLT,VN1(NYPLT),VN2(NYPLT)
8     	FORMAT('                      Y = ',I4,' ',A2,A16)
      	WRITE(*,*)' *********************************'
! 	inewton=0
!       call FSS_TextSize(14)
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
!       call FSS_TextSize(9)      
      	WRITE(*,*)' *********************************'
      	WRITE(*,*)'  CONTOUR MENU'
      	WRITE(*,*)'   0 = Return'
      	WRITE(*,*)'   1 = '
      	WRITE(*,*)'   2 = Absolute contour routine'
      	WRITE(*,*)'   3 = Choose variable to hold constant'
      	WRITE(*,*)'   4 = Plot a contour'
      	WRITE(*,*)'   5 = Move to a new contour (single steps)'
      	WRITE(*,*)'   6 = Reference points'
      	WRITE(*,*)'   7 = Print assemblage information to Output window'
      	WRITE(*,*)'   8 = Go to global menu'
      	WRITE(*,*)'   9 = Go to plotter menu'
      	WRITE(*,*)'  10 = Gibbs/Newton switch'
      	WRITE(*,*)'  11 = Set Contour Limits'
      	WRITE(*,*)'  15 = Reaction facing (variance = 1 only)'
	write(*,*)'  20 = Begin/Save menu (Save Illustrator file)'
      	write(*,*)'---------------------------------'
      	write(*,*)'   -1 = Change mineral assemblage'
      	write(*,*)'   -2 = Set output length'
	write(*,*)'   -4 = Pick line color'
      	WRITE(*,*)' CHOOSE OPTION'
      	read(*,*)(jopt)

! 	if(jopt.lt.0)call global(jopt)
! ----------------------------------------------
	if(jopt.eq.-1)then
!      	if we are changing the assemblage, we must zero out any "new" variables
!      	because there is no way to determine whether they will still be in the 
!      	list after a mineral is removed.
      	NVAR=NVAR-numNew
      	NEQ=NEQ-numNew
      	numNew=0
       	call change
 !      	call rexn		now called inside change
	go to 1
	endif
! ----------------------------------------------
	if(jopt.eq.-2)then
		call AWE_setoutput		! Output length
	       	goto 1
	endif
! ---------------------------------------
	if(jopt.eq.-4)then
!		myColor = 1		! Keep the same myColor as previously
		call PickAWEColor(myColor)
		CurrentColor = myColor				! CurrentColor is in PlotStuff.inc common block
		goto 1
		endif
! ----------------------------------------------
      
	if(jopt.eq.0) return

! ------Begin/Save routine-------------
	if(jopt.eq.20)then
		call GibbsBegin
		go to 1
	endif
! ------Set contour limits routine-------------
	if(jopt.eq.11)then
		call SetContourLimit
		go to 1
	endif


! ------Reaction facing routine-------------
	if(jopt.eq.15)then
!		call RxnFacing
		go to 1
	endif

!  ----inewton------------------------------------------
	if(jopt.eq.10)then
		call GibbsorNewton
		go to 1
	endif
!  ----SELECT PEN-----------------------------------
	if(jopt.eq.1)then
	      GO TO 1
	endif
! -----CHOOSE VARIABLE ---------------------------
	if(jopt.eq.2)then
!      CHOOSE VARIABLE using absolute values
        !zinc4 = zinc

200	continue
      	if(nvar-neq.le.1)then
        	call FSS_alert('ALERT!!',' Variance is <= 1.  No parameters can be held constant')
        	go to 1
        	endif
	write(*,*)'----------------------'
	write(*,*)'Current values of variables'
	do 220 i = 1,nvar
	WRITE(*,221)I,VN1(I),VN2(I),ALLX(1,i)
221	FORMAT(I3,1X,A2,A16,5F12.4,I5)
220	continue

        NeedToPick = Nvar-Neq-1
      	do 205 i = 1,NeedToPick
      	Picked(i) = ncont(i+1)
205   	continue
	myColor = currentColor
      	Call PickAbsContour(nvar,needToPick,Picked,valuesPicked,myColor,iok)
!      	Call PickMonitors(nvar,needToPick,Picked,valuesPicked,iok)
	currentColor = myColor			! currentColor is in PlotStuff.inc common block
        if(iok.eq.1)go to 1
	write(*,*)' Picked and valuesPicked ',picked(1),valuesPicked(1)
!      	if(zinc4.le.0.01.or.zinc4.ge.20.)zinc4 = 1.
!      	zinc = zinc4
	zinc = 1.0d0
      	XINCST=(xmax-xmin)*zinc/100.D0
      	YINCST=(ymax-ymin)*zinc/100.D0
!      	PROGRAM USES BARS
!      	adjust because pressure axis IS  KBARS

      	IF(NYPLT.EQ.2)YINCST=YINCST*1000.D0
!     	set the variable convert, which is the tolerance for angles in routine compute
!     	convert = (5 degrees)*zinc expressed in radians
!        	(5 degrees = .087266 radians)
!      	this way if zinc = 2 the tolerance is 10 deg, if zinc = 3 tol = 15 etc. 
      	convert = .087266d0 * zinc

! 	First, we must use steps routine to get to the desired contour
! Question: should I code this for Newton only, or for either? Let's try for either
! 	Calculate the deltas, based on the absolutes in valuesPicked
        do 206 i = 1,NeedToPick
	! use i+1 because the first monitor parameter is the X or Y axis
        SMon(i+1) = Picked(i)				! used in steps
!        SDel(i+1) = ValuesPicked(i)/Float(SNstep)	! used in steps
        SDel(i+1) = ValuesPicked(i)-allX(1,picked(i))		! used in steps
       	ncont(i+1) = Picked(i)				! used in contour routine
206     continue

	call StepToAllX(nXPlt,iok)
	if(iok.eq.1)go to 1
	if(iNewton.eq.1)then
		call NewtContour(nXPlt,1)		! contour routine for Newton's method
	else
		call GibbsContour(nXPlt,1)
	endif
!       	CALL RESET(2,1)				! reset to starting values 11/3/00

      	GO TO 200
	endif

! -----CHOOSE VARIABLE ---------------------------
	if(jopt.eq.3)then
!      CHOOSE VARIABLE using deltas

      	if(nvar-neq.le.1)then
        	call FSS_alert('ALERT!!',' Variance is <= 1.  No parameters can be held constant')
        	go to 1
        	endif

        NeedToPick = Nvar-Neq-1
      	do 305 i = 1,NeedToPick
      	Picked(i) = ncont(i+1)
305   	continue
!        zinc4 = zinc
!      	Call PkCon(nvar,VN1,VN2,NeedToPick,Picked,zinc4,iok)
!     	Call PkCon(nvar,VN1,VN2,NeedToPick,Picked,iok)
	myColor = currentColor
      	Call PickContour(nvar,NeedToPick,Picked,myColor,iok)
 	currentColor = myColor
       	if(iok.eq.0)then
        	do 306 i = 1,NeedToPick
        	ncont(i+1) = Picked(i)
306     	continue

!      		if(zinc4.le.0.01.or.zinc4.ge.20.)zinc4 = 1.
!      		zinc = zinc4
		zinc = 1.
      		XINCST=(xmax-xmin)*zinc/100.D0
      		YINCST=(ymax-ymin)*zinc/100.D0
!       	PROGRAM USES BARS
!       	adjust because pressure axis IS  KBARS

      		IF(NYPLT.EQ.2)YINCST=YINCST*1000.D0
!      		set the variable convert, which is the tolerance for angles in routine compute
!      		convert = (5 degrees)*zinc expressed in radians
!            		(5 degrees = .087266 radians)
!      		this way if zinc = 2 the tolerance is 10 deg, if zinc = 3 tol = 15 etc. 
      		convert = .087266d0 * zinc

      	endif

      	GO TO 1
	endif
! -----PLOT CONTOUR-----------------------------
	if(jopt.eq.4)then
	if(iNewton.eq.1)then
		call NewtContour(nXPlt,1)		! contour routine for Newton's method
		go to 1
		else
		call GibbsContour(nXPlt,1)
		endif

	endif
! -----STEP---------------------------------------
	if(jopt.eq.5)then
      		idraw = 0
      		CALL STEPS
      		idraw = 1
      		GO TO 1
	endif
! ---REFERENCE POINTS--------------------------
	if(jopt.eq.6)then
       		WRITE(*,*)' REFERENCE POINT OPTIONS'
       		WRITE(*,*)'   0 = Return  '
       		WRITE(*,*)'   1 = Set reference point'
       		WRITE(*,*)'   2 = Return to reference'
       		WRITE(*,*)'   3 = Return to starting conditions'
       		READ(*,*)IREFER
       		IF(IREFER.EQ.1) THEN
        		DO 610 I=1,nvar
           		ALLX(4,I)=ALLX(1,I)
610      		CONTINUE
       		ENDIF
       		IF(IREFER.EQ.2) CALL RESET(4,1)
       		IF(IREFER.EQ.3) CALL RESET(2,1)
       		GO TO 1
	endif
!-----PRINT--------------------------------------
	if(jopt.eq.7)then
      		CALL PRINTT(1)
      		GO TO 1
	endif
! ----GLOBAL------------------------------------------
	if(jopt.eq.8)then
		CALL GLOBAL(jopt)
         	GO TO 1
	endif
!----------plotter stuff------------------------------
	if(jopt.eq.9)then
		CALL PLOTIN
        	go to 1
	endif
!------------------------------------------------
	go to 1
      END

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE GibbsContour(firstMon,iDrawLine)
    	use AWE_Interfaces
	use MyCanvas

! 	i is "firstmon", which is not needed in the GibbsContour routine (only NewtonContour)
! 	iDrawLine = 1 to draw line automatically; = 0 not to draw line
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
	include "Singular.inc"
! *****************************************
!      local variables
      INTEGER*4 i,j,L,iplot,izero,PlotCounterSave,iDrawLine,firstMon
      character*64 ppp
      real*4 xxx,yyy

! --------------------------------------------------------


!     	set some defaults for Contour routine
	idraw=1
       	ISYMB=0
!      	DEFAULT PEN # and symbol
!      	Default increment size
!      	ZINC=3		! increment set in calling routine; in common/singular/
      	XINCST=(xmax-xmin)*zinc/100.D0
      	YINCST=(ymax-ymin)*zinc/100.D0
!       PROGRAM USES BARS
!      	adjust because pressure axis IS  KBARS
      	IF(NYPLT.EQ.2)YINCST=YINCST*1000.D0
!      	set the variable convert, which is the tolerance for angles in routine compute
!      	convert = (5 degrees)*zinc expressed in radians
!            (5 degrees = .087266 radians)
!      	this way if zinc = 2 the tolerance is 10 deg, if zinc = 3 tol = 15 etc. 
      	convert = .087266d0 * zinc




! -----PLOT CONTOUR using Gibbs method (differential thermodynamics)-----------------------------
400   CONTINUE

!      set up some necessary information for automatic contour routine
!      define a slope of "one" as a diagonal line on the plot
      ONE=(ymax-ymin)/(xmax-xmin)
!       fix because pressure is different  KBARS
      IF(NYPLT.EQ.2)one=one*1000.D0

!       PLOT CONTOUR
!       IPLOT IS FLAG HERE, =1 FIRST PASS, =2 SECOND PASS
      iplot=1
!      NSTEP is number of finite difference steps
      nstep=1
!      THIS SETS REFERENCES FOR CONTOUR
      DO 403 I=1,nvar
!         3 IS WHERE CONTOUR STARTED
         ALLX(3,I)=ALLX(1,I)
403   CONTINUE
      mon(1)=NXPLT
      do 401 i=2,6
         deltax(I)=0.D0
401      mon(I)=ncont(I)
!      Set up array IPOINT to contain pointers to non-monitor parameters
      J=0
      Do 410 I=1,NVAR
      DO 411 L=1,NVAR-NEQ
      IF(MON(L).eq.I)go to 410
411   continue
      J=J+1
      IPOINT(J) = I
410   CONTINUE
	call zeroiLong
	ilong(1) = 1	! just to be sure we don't generate to much output

!      -----START MAIN LOOP-----------
	Plotcounter=0
450   continue


!      ICONT IS FLAG FOR COMPUTE
!          ICONT=0  COMPUTE CALLED FROM SOMEWHERE ELSE
!              1  FIRST COMPUTATION FROM HERE
!              2  SECOND AND LATER FROM HERE
      ICONT=1
!     MOVE PLOTTER PEN TO START
      XXX=ALLX(1,NXPLT)
      YYY=ALLX(1,NYPLT)
      IF(NYPLT.EQ.2)YYY=YYY/1000.
  !    call FSS_SetPort(3)
!	!call !scolor(icolor)
	if(isymb.eq.1)then
		call symb(XYPlot,XXX,YYY,1.,1.)      ! plot a symbol if needed
		endif
	!call !scolor(1)
!     Write X and Y to save arrays
	    Plotcounter=Plotcounter+1
	    if(Plotcounter.gt.plotCounterMax)then
	       call FSS_Alert('ALERT!!','Too many points: redimension SaveRxn')
	       Plotcounter=1
	       go to 999
	       endif
	    Plotsave(Plotcounter,1)=XXX
	    Plotsave(Plotcounter,2)=YYY
	    PenSave(PlotCounter) = 0

!      The X axis  (NXPLT) IS ALWAYS INDEPEND (INCREMENTED) VARIABLE to start
!      The Y axis  (NYPLT) IS DEPEND (CALCULATED) VARIABLE
      do 453 L=1,NEQ
      if(ipoint(L).eq.NYPLT)then
        ndep=L
        go to 452
      endif
453   continue    
      write(*,*)'ndep not found in subroutine Contour (loop 453)'
      pause 'hit return to continue'
452   continue

      mon(1)=NXPLT

!      set up for INITIAL CALCULATION WITH VERY SMALL STEP (1% of maximum step)
!      set initial values of xinc and yinc, depending on which half of the contour we are drawing
      if(iplot.eq.1)then
         XINC=XINCST/100.
         YINC=YINCST/100.
      endif
      if(iplot.eq.2)then
         XINC=-XINCST/100.
         YINC=-YINCST/100.
      endif
      izero=0
        call c2ompute(izero)
      if(izero.eq.1) go to 999

!     Go!

      ANGLE=0
!     WRITE(*,*)'  slope,   xinc,   yinc,   ANGLE'
! ----------------START MINOR LOOP

      
460   CONTINUE


      CDELX=ALLX(1,NXPLT)-ALLX(5,NXPLT)   !store the last value of X and Y
      CDELY=ALLX(1,NYPLT)-ALLX(5,NYPLT)

!      NXPLT IS ALWAYS INDEPEND (INCREMENTED) VARIABLE
!      NYPLT IS DEPEND (CALCULATED) VARIABLE
!      SIZE OF INCREMENT IS CALCULATED IN sub COMPUTE.  IT IS DEPENDANT
!         ON THE CURVATURE OF THE LINE

!     TEMP=ANGLE*180/3.1415
!      DISPLAY IN DEGREES
!     WRITE(*,462)slope,xinc,yinc,TEMP
462   format(4E15.5)
      ICONT=2
      IZERO=0
        call c2ompute(izero)
!       NOTE ALLX(5,L) IS SET IN sub COMPUTE
!      PLOT PREVIOUS POSITION
       xxx=ALLX(5,NXPLT)
       YYY=ALLX(5,NYPLT)
       IF (NYPLT.EQ.2)yyy=YYY/1000.
	if(isymb.eq.1)then
		call symb(XYPlot,XXX,YYY,1.,1.)	! Draw a symbol to show progress
		endif
!     Write X and Y to save arrays
	Plotcounter=Plotcounter+1
	if(Plotcounter.gt.plotCounterMax)then
	  call FSS_Alert('ALERT!!','Too many points: redimension SaveRxn')
	  Plotcounter=1
	  go to 999
	  endif
	Plotsave(Plotcounter,1)=XXX
	Plotsave(Plotcounter,2)=YYY
	PenSave(PlotCounter) = 1

      if(izero.eq.1)go to 470       ! was everything ok in compute?

!     check to see if user escaped with escape key
!      call trapit(izero)
      if(izero.eq.1)then
!            if(jps.eq.1.and.idraw.eq.1)CALL PPENUP
            write(*,*)' User terminated contour operation'
            pause 'Hit return to continue...'
              call reset(3,2)
            ICONT=0                 ! reset ICONT to zero because we are done contouring
            go to 999
            endif

!     check to see if pressure is too low
!       IF(PB.LT.250)go to 470

!      check to see if point is out of bounds
      XXX=ALLX(1,NXPLT)
      IF(NXPLT.EQ.2)XXX=XXX/1000.
      if(xmax.gt.xmin)then          
            if(XXX.gt.xMaxBound.or.XXX.lt.xMinBound)go to 470      ! then maximum value is on right
            else
            if(XXX.lt.xMaxBound.or.XXX.gt.xMinBound)go to 470      ! then maximum value is on left
            endif
      YYY=ALLX(1,NYPLT)
      IF(NYPLT.EQ.2)YYY=YYY/1000.
      if(ymax.gt.ymin)then
            if(YYY.GT.yMaxBound.OR.YYY.LT.yMinBound) go to 470     ! YMAX on top
            else
            if(YYY.LT.yMaxBound.OR.YYY.GT.yMinBound) go to 470     ! Ymax on bottom
            endif

! 	Check for mass balance error
	if(imass.eq.1)then
		call CheckMassBal
		endif


!      	--IN BOUNDS--
          GO TO 460
470   	CONTINUE
!      	-- out of bounds.--
!      	check to see if both halves of plot are finished
      	IF(IPLOT.EQ.1)then	! First half is plotted - go plot second half
!          	if(jps.eq.1.and.idraw.eq.1)CALL PPENUP
         	IPLOT=2
         	call reset(3,2)
! 		save value of PlotCounter so we can plot label at this point later
		PlotCounterSave = PlotCounter
         	go to 450
      	endif

! 	both halves calculated.  Draw plot in picture
        call reset(3,2)
	PenSave(PlotCounter+1) = 0	! Just to force a penup in postscript routine
! 	Now replot in the Picture
	if(iDrawLine.eq.0)go to 999
!	call FSS_OpenPicture(3)
!	!call !scolor(icolor)
	do 480 i = 1,PlotCounter
	XXX = PlotSave(i,1)
	YYY = PlotSave(i,2)
	ipen = PenSave(i)
	Call Plot(XYPlot,XXX,YYY,ipen)
	call pplot(XXX,YYY,ipen)
	if(PenSave(i+1).eq.0)Call Ppenup
480	continue

! 	label the part of the curve just drawn with the value of the 
! 	variable held constant
	XXX = PlotSave(PlotcounterSave,1)
	YYY = PlotSave(plotcounterSave,2)
!       IF (NYPLT.EQ.2)yyy=YYY/1000.	! Plotsave is already in bars
        call write_millimoles(ppp)
!        WRITE(ppp,476)(ALLX(1,NCONT(I)),I=2,NVAR-NEQ)
!476     FORMAT(6(f10.3))
	call TextOnPlot(XYPlot,xxx,yyy,ppp,12)
	call pslab(ppp,(xxx-.03*(xmax-xmin)),yyy,12,0)

	!call !scolor(1)
!	call FSS_ClosePicture(3)

      	ICONT=0
999	continue
	!call !scolor(1)
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      SUBROUTINE write_millimoles (ppp)
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
! *****************************************
!      local variables
      INTEGER*4 i
      real*4 temp(10)
      character*64 ppp

	do 10 i = 2,nvar-neq
	if(Vn1(ncont(i))(1:2).eq."M_")then
		temp(i) = ALLX(1,ncont(i))*1000
		else
		temp(i) = ALLX(1,ncont(i))
		endif
10	continue
        WRITE(ppp,476)(temp(i),I=2,NVAR-NEQ)
!        WRITE(ppp,476)(ALLX(1,NCONT(I)),I=2,NVAR-NEQ)
476     FORMAT(6(f10.4))

	return
	end



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE NewtContour(firstMon,iDrawLine)
    	use AWE_Interfaces
	use MatrixArrays
	use MyCanvas
! 	firstmon
! 	iDrawLine = 1 to draw line automatically; = 0 not to draw line

      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
	include "Singular.inc"
! 	include "Solute.inc"
! *****************************************
!     local variables
      INTEGER*4 i,j,iplot,izero,PlotCounterSave,firstMon,iDrawLine
      character*64 ppp
      real*4 xxx,yyy

! -------------------------------------------------------

!     	set some defaults for Contour routine
!      	DEFAULT PEN # and symbol
      	idraw=1
!      	ISYMB=0			! don't draw a symbol
      	ISYMB=1			! this will draw a symbol at every point for debugging

! 	zinc is in common /singular/ and is the % of X and Y axis for contour interval
      	XINCST=(xmax-xmin)*zinc/100.D0
	XINC = XINCST
      	YINCST=(ymax-ymin)*zinc/100.D0
      	IF(NYPLT.EQ.2)YINCST=YINCST*1000.D0		! convert to bars for calculation
	YINC = YINCST
	

! ----PLOT CONTOUR-----------------------------
! 	This routine plots a contour using Newton's method
! 	It is similar to the "Contour" routine algorithm, except there we go through an elaborate back-up and replot
! 	iteration to ensure that when the contour is tight (large curvature), we make very small steps.
! 	Here, we just pick a value of one axis (NXPLOT or NYPLOT) and calculate the other.
! 	The variable picked as "independent" depends on the slope of the contour.  If the slope is steep (large dY/dX), then the Y
! 	variable is independent.  If it is flat, the X is independent.
400   CONTINUE

! 	First make a call to compute to find the correct P-T-X values at the starting conditions
! 	Also find the slope of the contour (dY/dX)


!     set up some necessary information for automatic contour routine
!     define a slope of "one" as a diagonal line on the plot
      ONE=(ymax-ymin)/(xmax-xmin)
      IF(NYPLT.EQ.2)one=one*1000.D0		! pressure is calculated in bars but plotted in kbar


      NSTEP=1				! only do 1 step in this routine.  Bounds are checked below

!      IPLOT IS FLAG: =1 FIRST PASS (increasing X or Y from start), =2 SECOND PASS (decreasing X or Y from start)
      iplot=1
	Plotcounter=0
      ICONT=0				! in this routine we call COMPUTE as if this is a STEPS operation (no contouring done in Compute)

! 	fill up the monitor parameters that are NOT either X or Y axis
	

500	continue				! start calculating a contour

      do 401 i=2,nvar-neq
      deltax(i)=0.			! for this first increment, we set deltas to 0
401   mon(i)=ncont(i)
      mon(1) = firstMon			! set first monitor to value passed through (either X or Y)

!     Set up array IPOINT and ndep
	if(firstMon.eq.nXplt)then
	 	call SetMon(nYplt)
	else
		call SetMon(nXplt)
	endif

      IF(iLong(4).eq.1)then
		write(12,*)' Monitors are:', (mon(i),i=1,nvar-neq)
		write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)
		endif

!     THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
      DO 403 i=1,nvar
      ALLX(3,i)=ALLX(1,i) 		! ALLX(3,I) IS WHERE calculations STARTED-save
403   CONTINUE

           
        IZERO=0
        call c2ompute(izero)		! initial call to get value of slope of contour "slope"
        IF(IZERO.EQ.1)GO TO 470


!     THIS SETS REFERENCES FOR CONTOUR
      DO 404 I=1,nvar
         ALLX(3,I)=ALLX(1,I)		! where the contour started
404   CONTINUE

!     MOVE PLOTTER PEN TO START
      XXX=ALLX(1,NXPLT)
      YYY=ALLX(1,NYPLT)
	IF(NYPLT.EQ.2)YYY=YYY/1000.
	if(isymb.eq.1)then
		call symb(XYPlot,XXX,YYY,1.,1.)      ! plot a symbol if needed
		endif
	    Plotcounter=Plotcounter+1
	    if(Plotcounter.gt.plotCounterMax)then
	       call FSS_Alert('ALERT!!','Too many points: redimension SaveRxn')
	       Plotcounter=1
	       go to 999
	       endif
	    Plotsave(Plotcounter,1)=XXX
	    Plotsave(Plotcounter,2)=YYY
	    PenSave(PlotCounter) = 0

! 	Now we have the starting P-T-X consistent with thermo data
	slope = xx(ndep,1)	! ndep is either nXplt or nYplt
! 	and the slope at the starting point

! 	Based on value of slope, determine whether nXplt (shallow) or nYplt (steep) should be the independent variable
! 	Determine whether NXPLT or NYPLT is the independent variable, based on slope calculated in last increment
       IF (DABS(SLOPE).GT.ONE)THEN
!         STEEP SLOPE: USE (y axis=NYPLT) AS independent variable TO CALCULATE DELTA NXPLT
          DELTAX(1)=YINC
          mon(1)=NYPLT
 	    call SetMon(nXplt)
       ELSE
!         SHALLOW SLOPE; USE x axis=NXPLT AS independent variable TO CALCULATE DELTA NYPLT
          DELTAX(1)=XINC
          mon(1)=NXPLT
 	    call SetMon(nYplt)
       ENDIF
      IF(iLong(4).eq.1)then
		write(12,*)' Slope, Monitors:',Slope, (mon(i),i=1,nvar-neq)
		endif

	if(slope.lt.0.0d0)yinc = -yinc	! if negative slope, then yinc must be (-) to start
		
! 	Calculate the rest of the contour

460   CONTINUE	! loop to plot half contour

!        Store new values for each independent variable (monitor parameter)
! 	This is where Neton's method increments composition, T or P
         DO 420 J=1,NVAR-NEQ
         i=mon(j)
         ALLX(1,i)=ALLX(1,i) + DELTAX(J)
420      continue
         call SetTPX


      IZERO=0
      call c2ompute(izero)
      if(izero.gt.0)go to 470       ! was everything ok in compute?



!     check to see if pressure is too low
!      IF(PB.LT.250)go to 470

!     check to see if point is out of bounds
      XXX=ALLX(1,NXPLT)
      IF(NXPLT.EQ.2)XXX=XXX/1000.
      if(xmax.gt.xmin)then          
            if(XXX.gt.xMaxBound.or.XXX.lt.xMinBound)go to 470      ! then maximum value is on right
            else
            if(XXX.lt.xMaxBound.or.XXX.gt.xMinBound)go to 470      ! then maximum value is on left
            endif
      YYY=ALLX(1,NYPLT)
      IF(NYPLT.EQ.2)YYY=YYY/1000.
      if(ymax.gt.ymin)then
            if(YYY.GT.yMaxBound.OR.YYY.LT.yMinBound) go to 470     ! YMAX on top
            else
            if(YYY.LT.yMaxBound.OR.YYY.GT.yMinBound) go to 470     ! Ymax on bottom
            endif

! 	Check for mass balance error
	if(imass.eq.1)then
		call CheckMassBal
		endif


!     --IN BOUNDS--Plot this position
!      NOTE ALLX(5,L) IS SET IN sub COMPUTE
      xxx=ALLX(1,NXPLT)
      YYY=ALLX(1,NYPLT)
      IF (NYPLT.EQ.2)yyy=YYY/1000.

	if(isymb.eq.1)then
		call symb(XYPlot,XXX,YYY,1.,1.)	! Draw a symbol to show progress
		endif
!     Write X and Y to save arrays
	Plotcounter=Plotcounter+1
	if(Plotcounter.gt.plotCounterMax)then
	  call FSS_Alert('ALERT!!','Too many points: redimension SaveRxn')
	  Plotcounter=1
	  go to 999
	  endif
	Plotsave(Plotcounter,1)=XXX
	Plotsave(Plotcounter,2)=YYY
	PenSave(PlotCounter) = 1

	
! 	Determine whether NXPLT or NYPLT is the independent variable, based on slope calculated in last increment
	slope = xx(ndep,1)	! ndep is either nXplt or nYplt
	if(mon(1).eq.nYplt)slope = 1.d0/slope	! invert slope if Y is independent variable
       IF (DABS(SLOPE).GT.ONE)THEN
!         STEEP SLOPE: USE (y axis=NYPLT) AS independent variable TO CALCULATE DELTA NXPLT
!         CDELY is the last change in Y.  If this is negative then the next change must also be negative
          CDELY=ALLX(1,NYPLT)-ALLX(5,NYPLT)
          IF(CDELY.LT.0.D0)YINC=-Dabs(YINC)
          DELTAX(1)=YINC
          mon(1)=NYPLT
 	    call SetMon(nXplt)
       ELSE
!         SHALLOW SLOPE; USE x axis=NXPLT AS independent variable TO CALCULATE DELTA NYPLT
!         NOTE CDELX IS LAST INCREMENT OF x axis
          CDELX=ALLX(1,NXPLT)-ALLX(5,NXPLT)
          IF (CDELX.LT.0.0D0)XINC=-Dabs(XINC)
          DELTAX(1)=XINC
          mon(1)=NXPLT
 	    call SetMon(nYplt)
       ENDIF



      GO TO 460		! loop back for another increment of this contour


! ------------------------------------------------------
!        done with this half of contour
470   CONTINUE	! out of bounds here
!     check to see if both halves of plot are finished
      IF(IPLOT.EQ.1)then	! First half is plotted - go plot second half
!         if(jps.eq.1.and.idraw.eq.1)CALL PPENUP
         IPLOT=2
         call reset(3,2)
	   XINC = -XINCST
	   YINC = -YINCST
! 	save value of PlotCounter so we can plot label at this point later
	PlotCounterSave = PlotCounter
      go to 500
      endif

! 	both halves calculated.  Draw plot in picture
      call reset(3,2)
	if(iDrawLine.eq.0)go to 999
	if(plotCounter.eq.0)go to 999	! don't try to plot if there are no points
	PenSave(PlotCounter+1) = 0	! Just to force a penup in postscript routine
! 	Now replot in the Picture
	do 480 i = 1,PlotCounter
	XXX = PlotSave(i,1)
	YYY = PlotSave(i,2)
	ipen = PenSave(i)
	Call Plot (XYPlot,XXX,YYY,ipen)
	call pplot(XXX,YYY,ipen)
	if(PenSave(i+1).eq.0)Call Ppenup 
480	continue

! 	label the part of the curve just drawn with the value of the 
! 	variable held constant
	XXX = PlotSave(PlotcounterSave,1)
	YYY = PlotSave(PlotcounterSave,2)
!         call plot(XYPlot,XXX,YYY,0)		Don't need to move to the spot
!         call FSS_move(-42,0)
        call write_millimoles(ppp)
!         WRITE(ppp,476)(ALLX(1,NCONT(I)),I=2,NVAR-NEQ)
!476      FORMAT(6(f10.3))
	call TextOnPlot(XYPlot,xxx,yyy,ppp,12)
	call pslab(ppp,(xxx-.03*(xmax-xmin)),yyy,12,0)


999	continue
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	subroutine SetMon(jdep)
      implicit none

! *****************************************
 	include "Assemb.inc"
	include "Monit.inc"
	include "Singular.inc"
! *****************************************
!     local variables
	integer*4 i,j,l,jdep

!     Set up array IPOINT to contain pointers to non-monitor parameters
      j=0
      Do 410 i=1,NVAR
      DO 411 l=1,NVAR-NEQ
      IF(MON(l).eq.i)go to 410
411   continue
      j=j+1
      IPOINT(j) = I
410   CONTINUE

!     The X axis  (NXPLT) IS ALWAYS INDEPEND (INCREMENTED) VARIABLE to start
!     The Y axis  (NYPLT) IS DEPEND (CALCULATED) VARIABLE
! 	This loop sets initial ndep to the y-axis value
      do 453 l=1,NEQ
!       if(ipoint(L).eq.NYPLT)then
      if(ipoint(l).eq.jdep)then
        ndep=l
        go to 452
      endif
453   continue    
452	continue

	return
	end
	
	




! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE StepToAllX(firstMon,izero)
    	use AWE_Interfaces
	use MatrixArrays
      implicit none

! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
 	include "Singular.inc"
! 	include "Solute.inc"
	include "TPSave.inc"
! *****************************************
!      	local variables
      	integer*4 i,j,istep,izero,jsymb,firstMon
	real*8 newtonStepOld



! 	First determine whether the slope of the contour is steep or gentle

!     	define a slope of "one" as a diagonal line on the plot
      	ONE=0.5*(ymax-ymin)/(xmax-xmin)
      	IF(NYPLT.EQ.2)one=one*1000.D0	! pressure is calculated in bars but plotted in kbar


      	ICONT=3				! in this routine we call COMPUTE to get slopes only

! 	fill up the monitor parameters that are NOT either X or Y axis
      	do 401 i=2,nvar-neq
!      	deltax(i)=0.			! for this first increment, we set deltas to 0
401   	mon(i)=ncont(i)
      	mon(1) = firstMon		! set first monitor to value passed through (either X or Y)
!        deltaX(1)=0.
	
!     	Set up array IPOINT and ndep
	if(firstMon.eq.nXplt)then
	 	call SetMon(nYplt)
	else
		call SetMon(nXplt)
	endif

      	IF(iLong(4).eq.1) write(12,*)' Non-Monitors are:',(IPOINT(I),I=1,NEQ)


!     	THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
      	DO 403 i=1,nvar
      	ALLX(3,i)=ALLX(1,i) 		! ALLX(3,I) IS WHERE calculations STARTED-save
403   	CONTINUE
        
        IZERO=0
        call c2ompute(izero)		! initial call to get value of slope of contour "slope"
        IF(IZERO.EQ.1)return
!      	call trapit(izero)	! check to see if user escaped with escape key
      	if(izero.eq.1)then
            write(*,*)' User terminated operation'
            pause 'Hit return to continue...'
            call reset(3,1)
	    return
            endif

! 	Now look at slope and determine whether we should increment nXPlt or nYPlt to step to the initial conditions
! 	slope = xx(ndep,1)	! ndep is either nXplt or nYplt (calculated in compute and passed in common/singlular

       IF (DABS(SLOPE).GT.ONE)THEN
!         STEEP SLOPE: USE (y axis=NYPLT) AS independent variable TO CALCULATE DELTA NXPLT
          DELTAX(1)=0.0
          mon(1)=NYPLT
 	    call SetMon(nXplt)
       ELSE
!         SHALLOW SLOPE; USE x axis=NXPLT AS independent variable TO CALCULATE DELTA NYPLT
          DELTAX(1)=0.0
          mon(1)=NXPLT
 	    call SetMon(nYplt)
       ENDIF
		
! 	Now increment using these monitors (constant nXplt or nYPlt) and find the desired P-T-X point


400   CONTINUE
	jsymb=3
!      	set monitor parameters
	if(iNewton.eq.1)then
		nStep = 1	! set to a large value, just to be safe
		newtonStepOld = newtonStep
		newtonStep = .05
	else
		nStep = 50	! set to a large value, just to be safe
	endif

      	do 412 i=2,nvar-neq	! deltax(1)= 0  (held constant)
      	deltax(i)=sdel(i)/Dfloat(nStep)
412   	mon(i)=smon(i)


!      	first draw the starting point on the screen if idraw=1
      	Plotcounter=0
	call stepplot(1)

!      	THIS SETS REFERENCES FOR resetting (only if there is a user abort (esc key))
      	DO 413 I=1,nvar
!         ALLX(3,I) IS WHERE calculations STARTED-save
         ALLX(3,I)=ALLX(1,I)
413   	CONTINUE


	icont = 0		! call contour and calculate steps
      do 499 istep=1,nstep

      	if(inewton.eq.1)then
!       	Store new values for each independent variable (monitor parameter)
! 		This is where Neton's method increments composition, T or P
        	DO 420 J=1,NVAR-NEQ
         	i=mon(j)
         	ALLX(1,i)=ALLX(1,i) + DELTAX(J)
420      	continue
         	call SetTPX
      	endif
         


        IZERO=0
        call c2ompute(izero)
        IF(IZERO.EQ.1)go to 499

!      check to see if user escaped with escape key
!      call trapit(izero)
      if(izero.eq.1)then
            write(*,*)' User terminated operation'
            pause 'Hit return to continue...'
            call reset(3,1)
	    go to 499
            endif


!      plot this finite difference iteration
	call stepplot(2)

! 	Check for mass balance error
	if(imass.eq.1)then
		call CheckMassBal
		endif

499   continue

	if(inewton.eq.1)newtonStep = newtonStepOld
! 	we should be at the correct P-T-X spot

	return
	end







