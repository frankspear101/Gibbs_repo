! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PLOTIN
	USE AWE_INTERFACES
	use MyCanvas
	implicit none
!	TYPE(AWE_Canvas) :: XYPlot
! *****************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
! *****************************************
!      local variables
      INTEGER iplotin,i,j,jaxis,iok,Pick2(5,2),it,jcol,k,kCur,myColor,Picked(10)
      CHARACTER TempTitle*80,dummy*4
	Character*80 TheList(200),WhatToPick
	Integer*4 TheOne,ListSize,iselect


1     CONTINUE
      write(*,*)
      write(*,*)' Current axes are'
      write(*,6001)NXPlt,VN1(NxPLT),VN2(NxPLT)
      write(*,6002)NyPLT,VN1(NyPLT),VN2(NyPLT)
6001  FORMAT(' X-axis = ',I4,5x,A2,A16)
6002  FORMAT(' Y-axis = ',I4,5x,A2,A16)
      write(*,*)
      WRITE(*,*)' ****************************'
      WRITE(*,*)' Plotting options menu'
      WRITE(*,*)'  0 = Return'
      WRITE(*,*)'  1 = Read new plot definition from configuration file'
      WRITE(*,*)'  ------- X-Y diagrams---------------------'
      WRITE(*,*)'  2 = Draw axes for X-Y diagrams'
      WRITE(*,*)'  3 = Choose plotting variables for X-Y diagrams'
!       write(*,*)'      T-P diagrams'
!       write(*,*)'      T-X diagrams'
!       write(*,*)'      P-X diagrams'
      WRITE(*,*)'  ------- Ternary diagrams-----------------'
      write(*,*)'  4 = Draw axes for ternary diagrams'
      write(*,*)'  5 = Select plotting variables for ternary diagram'
      WRITE(*,*)'  ------- Miscellaneous---------------------'
!      write(*,*)'  6 = Get X-Y position of mouse'
 !     write(*,*)'  7 = Save graphic window as a PICT file'
      write(*,*)'  8 = Save graphic window as a Postscript file'
      write(*,*)'  9 = Pick pen color'
      WRITE(*,*)' CHOOSE'
      iplotin=0
      read(*,*)(iplotin)

      if(iplotin.eq.0)go to 999
      

      if(iplotin.eq.9)then
!	call PenColor
      goto 1
      endif

      if(iplotin.eq.6)then
! -------------------------------------------------------
!      Get X-Y position of mouse
! -------------------------------------------------------
!       call GetXY(Xloc,Yloc,iflag)
      goto 1
      endif
! -------------------------------------------------------
! ---------------option 7----------------------
      if(iplotin.eq.7)then
    !        call savepicture
            go to 1
            endif

! ---------------option 8----------------------
      if(iplotin.eq.8)then
!           PSFileName = Trim(filein)//'.ps'        	! default ps file name
           call psopcl(5)    				! save old PostScript scratch file
           GO TO 1
           ENDIF

! 
! -----------------------------------------------
!      Read new plot information from configuration file

      	if(iplotin.eq.1)then

!      	Build a list of analysis titles for listmanager
      	rewind(23)      			! rewind the Gibbs.PlotDefinitions file
      	READ(23,'(A)')dummy			! read the header

       DO 416 J=1,200
          READ(23,'(A)')dummy          	! read the row of stars
          READ(23,'(A)',END=4170)TempTitle
          WRITE(TheList(j),423)j
423      format(I4)
          TheList(j) = Trim(TheList(j))//' '//Trim(TempTitle)

          DO 415 I=1,5
            READ(23,'(A)',END=4170)dummy
415         CONTINUE
416      CONTINUE
4170     CONTINUE



!     Put the list up and let the user pick one plot type
      iSelect = 1
      ListSize = j-1
      WhatToPick = 'Select a plot definition '

      call AWE_Pick_One(TheOne,ListSize,TheList,WhatToPick,iok,iSelect)
      if(iok.ne.0)go to 1     			! user hit cancel

	jaxis = TheOne
        iBig = jaxis         			! iBig is in common block so we can set default on startup
      	REWIND(23)
      	READ(23,'(A)')dummy			! read the header

!         STARTING OVER, AGAIN IGNORE FILE TITLES
         DO 419 I=1,(jaxis-1)*7
            READ(23,'(A)')dummy     		! ignore the unwanted plot definitions
419      CONTINUE

         READ(23,'(A)')dummy     !Read line of stars
         READ(23,'(A)')PlTitle
         READ(23,*)NXPLT,NYPLT
         READ(23,'(A)') XLAB
         READ(23,*)XOR,XMIN,XMAX,XLEN,NXSTEP,NXDEC
         READ(23,'(A)') YLAB
         READ(23,*)YOR,YMIN,YMAX,YLEN,NYSTEP,NYDEC

	PlotIndex(1,1) = NXPLT		
	PlotIndex(1,2) = NYPLT
	NoTieMinerals = 1

	write(*,*)' Plot data '
	write(*,*)Pltitle
	write(*,*)nxplt,nyplt
	write(*,*)xlab
	write(*,*)xor,xmin,xmax,xlen,nxstep,nxdec
	write(*,*)ylab
	write(*,*)yor,ymin,ymax,ylen,nystep,nydec
	write(*,*)' '
	write(*,*)' '

      go to 1

      endif


! --PLOT AXES-------------------------------------------------------

      if(iplotin.eq.2)then
		call SetPlot(iok)
		if(iok.eq.1)go to 1
		myColor = currentColor		! save the previous color
		currentColor = 1		! Black for axes
		CALL PlotAxes(XYPlot)
		if(itrip.eq.1)then
			!call Al2SiO5(3,0)       		!3 plots Berman's triple point
			endif
		if(ialphabeta.eq.1)then
			call abqtz(0)
			endif
		currentColor = myColor		! set back to the previous color
		iternary = 0
		xMinBound = xmin
		xMaxBound = xmax
		yMinBound = ymin
		yMaxBound = ymax
		go to 1  
      		ENDIF
         
! ------------change plotting variables-------------------------------------
      if(iplotin.eq.3)then

! 	we have an XY plot and we want to pick monitor params to plot

	write(*,*)
	write(*,*)' How many minerals do you want to plot?'
	write(*,*)' choose 1 to plot a single mineral'
	write(*,*)' 2-5 to plot tie lines'
	read(*,*)NoTieMinerals

	if(NoTieMinerals.ge.1.and.NoTieMinerals.le.5)then

!	do 107 i=1,NoTieMinerals

!		Pick2(1,1) = 0	     		! force no default
!		Pick2(1,2) = 0	     		! force no default
! 		Picked(1) = NXPlt
! 		Picked(2) = NYPlt
		!call PkAxe(1,2,Picked,iok)
		call PkXYAxe(nvar,NoTieMinerals,Pick2,iok)
		if(iok.eq.1)go to 1

		do 107 i = 1,NoTieMinerals		
		PlotIndex(i,1) = Pick2(i,1)	! Here PlotIndex contain monitor params
		PlotIndex(i,2) = Pick2(i,2)
! 		Store the names of the plotting variables so they can be reset if a new problem is read in
		plotIndexName(i,1) = trim(VN1(plotIndex(i,1)))//VN2(plotIndex(i,1))
		plotIndexName(i,2) = trim(VN1(plotIndex(i,2)))//VN2(plotIndex(i,2))
	
		NXPLT = PlotIndex(i,1)		! to maintain compatibility with Contour routine
		NYPlt = PlotIndex(i,2)
		axisIsPhase(1,1) = 0		! default is plot an independent variable
		axisIsPhase(2,1) = 0
107	continue

108	continue
	if(imass.eq.1)then
		JCOL = TANDP + NX
		if(NXPLT.gt.jcol.or.NYPLT.gt.jcol)then
			write(*,*)' Are you plotting a mineral moles/volume/radius?'
			write(*,*)' 0 = no; 1 = yes'
			read(*,*)i
			if(i.ne.0)then
				write(*,*)' Please specify which phase to plot (sorry I am so stupid)'
			!	DO 390 k = 1,numPh
				Do 390 kCur = 1,numPh
				k = asmCurrent(kCur)
				write(*,200)K,PHNAME(K)
200     			FORMAT(' ',I4,4X,A32)
390				continue
				read(*,*)K
				if(K.le.0.or.K.gt.numPh)then
					write(*,*)' You goofed. Try again'
					go to 108
				endif
				write(*,*)' Which axis? Type 1 for X axis or 2 for Y axis'
				read(*,*)i
				axisIsPhase(i,2)=k		! contains the phase number to plot
				write(*,*)' Choose how to plot'
				write(*,*)' 0 = Oops, I goofed'
				write(*,*)' 1 = Moles of phase'
				write(*,*)' 2 = Volume of phase'
				write(*,*)' 3 = Radius of phase'
				read(*,*)axisIsPhase(i,1)	! switch as to how to plot the phase
			endif
		endif
	endif							! end if imass = 1
      endif             					! end choose tie minerals
	iternary = 0


      go to 1
      endif	! end option 3


!     --Plot TRIANGULAR diagrams ---------------------------------------------
      if(iplotin.eq.4)then
401	continue
	write(*,*)'0 = return'
	write(*,*)'1 = pick plotting coordinates'
	write(*,*)'2 = draw triangle'
	read(*,*)it
	if(it.eq.0)go to 1

	if(it.eq.1)then
!	Choose apex of triangle
	picked(1) = 0
	!call PkAxe(2,3,Picked,iok)	!2 here specifies pick from CompPlot list
	if(iok.eq.1)go to 1
	
	do 410 i = 1,3
410	Apx(i) = Picked(i)

	Llab = CompPlotName(Apx(1))
	Rlab = CompPlotName(Apx(2))
	Tlab = CompPlotName(Apx(3))
		
	go to 401
	endif
	
	if(it.eq.2)then


      !call settriangle(XorT,YorT,Llab,Rlab,Tlab,scale,itic,Pltitle,iok)

      if(iok.eq.1)go to 1  !cancel button was selected

!     do PostScript file stuff
      write(*,*)
      i=0
	!call YesNo(i,'Save PostScript File?')
	if (i.eq.2)go to 1		!cancel
      if(i.eq.1)then	!yes
!            PSFileName = Trim(filein)//'.ps'        !default ps file name
            call psopcl(5)    !save old PostScript scratch file
            endif
!      call volget(voltemp)
!      call volset(volrefnum)
!      close(PSUNIT)         !close old PostScript scratch file
      call psopcl(4)    !open a new scratch file
!      call volset(voltemp)


	iternary = 1
      CALL USER(xorT,0.,1.,10.,yorT,0.,1.,10.)
!	Call FSS_OpenPicture(3)	!sets up for a new plot, or continue with an old one

	!call triaxis(Llab,Rlab,Tlab,scale,itic,Pltitle)
!	if(jps.eq.1)call tripaxis(Llab,Rlab,Tlab,scale,itic,Pltitle)

!	Call FSS_ClosePicture(3)
      go to 1
	endif
	
	go to 1


      endif
      
!------------change plotting variables-------------------------------------
      if(iplotin.eq.5)then

!	we have an ternary plot and we want to pick composition params to plot

	write(*,*)' How many minerals do you want to plot?'
	write(*,*)' choose 1 to plot a single mineral'
	write(*,*)' 2-5 to plot tie lines'
	read(*,*)NoTieMinerals

	if(NoTieMinerals.ge.1.and.NoTieMinerals.le.5)then

		Picked(1) = 0	     !force no default
		!call PkAxe(3,NoTieMinerals,Picked,iok)	!3 here specifies pick from mineral list
		if(iok.eq.1)go to 1
		
		do 209 i = 1,NoTieMinerals
		PlotIndex(i,1) = Picked(i)	!Here PlotIndex contains mineral picked
209		continue


      endif             !end choose tie minerals
	iternary = 1
      go to 1
      endif

!-----End and return--------------------

999   continue


      return            !plot done, return from whence we came
! -----------------------------
      END

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Gplot(Xp,ipenG,jpsG)
	USE AWE_INTERFACES
	use MyCanvas

!	Routine to plot in either X-Y cartesian or A-B-! barycentric coords
!	jps is jps switch.  =1 if we are to make a ps call
!	The switch iternary is in common block Plotted
!	iternary = 0 for cartesian plot
!	iternary = 1 for ternary plot
      implicit none
!	TYPE(AWE_Canvas) :: XYPlot
!*****************************************
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
!*****************************************
!     local variables

	integer*4 ipenG,jpsG
	real*4 Xp(4),Left,Right,Top,sum


!****NOTE: here jps is used as a local variable

	select case (iternary)

	case(0)		!plot cartesian coordinates
	
	Call Plot(XYPlot,Xp(1),Xp(2),ipenG)
         if(jpsG.eq.1) then
            call pplot(Xp(1),Xp(2),ipenG)
	    endif	
	
	case(1)		!plot ternary diagram
!	Xp,Yp,Zp are to be normalized to 1.0 and plotted on a ternary
!	with A = left, B = up, ! = right
!	assumes xor,yor is lower left corner of triangle
!	assumes xmin = ymin = 0
!	assumes xmax=ymax = 1
	sum = Xp(1)+Xp(2)+Xp(3)
	left  = Xp(1)/sum
	right = Xp(2)/sum
	top   = Xp(3)/sum

	!call triplot(XYPlot,left,right,top,scale,ipenG)
         if(jpsG.eq.1) then
  	    !call tripplot(left,right,top,scale,ipenG)
	    endif	
	
	case default
	
	end select

	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Gsymb(Xp,isymbG,isizeG,jpsG)
	USE AWE_INTERFACES
	use MyCanvas

!	Routine to plot in either X-Y cartesian or A-B-! barycentric coords
!	jps is jps switch.  =1 if we are to make a ps call
!	The switch iternary is in common block Plotted
!	iternary = 0 for cartesian plot
!	iternary = 1 for ternary plot
      implicit none
!	TYPE(AWE_Canvas) :: XYPlot
!*****************************************
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
!*****************************************
!     local variables

	integer*4 isymbG,isizeG,jpsG
	real*4 Xp(4),Left,Right,Top,sum


!****NOTE: here jps is used as a local variable

	select case (iternary)

	case(0)		!plot cartesian coordinates
	
        call symb(XYPlot,Xp(1),Xp(2),1.,1.)      !plot a symbol to show where we are
!         if(jpsG.eq.1) then
!            call psymb(Xp(1),Xp(2),isymbG,isizeG)
!	    endif	
	
	case(1)		!plot ternary diagram
!	Xp,Yp,Zp are to be normalized to 1.0 and plotted on a ternary
!	assumes xor,yor is lower left corner of triangle
!	assumes xmin = ymin = 0
!	assumes xmax=ymax = 1
	sum = Xp(1)+Xp(2)+Xp(3)
	left  = Xp(1)/sum
	right = Xp(2)/sum
	top   = Xp(3)/sum

	!call trisymb(XYPlot,left,right,top,scale,isymbG,isizeG)
         if(jpsG.eq.1) then
  	    !call tripsymb(left,right,top,scale,isymbG,isizeG)
	    endif	
	
	case default
	
	end select

	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PRINTT(IPRINT)
      implicit none
!      ROUTINE TO PRINT CURRENT DATA ON TERMINAL AND DISK
! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
 	include "Output.inc"
! *****************************************
!      local variables
      integer*4 i,ii,j,k,L,LL,IZERO,Lsite,iprint,kCur
      real*8 DX,SDX
      real*4 sum

!       COMPUTE MOLAR VOLUME OF PHASE FOR NEW COMPOSITION, T, P AND PF
            izero=0
            CALL AllKduTPX(1,izero)
!      note: vp0, vp1 and vp2 are volumes of phases based on the moles present
!            these are not modes!!

!      IPRINT = 0 Do not print monitors (causes a bomb unless NVAR and NEQ are set
!      IPRINT = 1 print everything

1     continue

      IF(KFLU.EQ.0)then
         pfluid=pb
         pfstart=pstart
         endif

	write(12,*)
	write(12,*)
        write(12,*)'******************************************************************************'  ! these characters signal the start of the data set
	if(PFluidSwitch.gt.0)then
        	write(12,*)'          T         P        Pfluid'
        	write(12,33)TSTART,pstart,pfstart
33      	FORMAT(' START  ',3F8.1)
        	write(12,34)TC,PB,PFLUID
34      	FORMAT(' CURRENT',3F8.1)
		else			! Pfluid = Plithostatic, so don't print out Pfluid
        	write(12,*)'          T         P   '
        	write(12,33)TSTART,pstart
        	write(12,34)TC,PB
		endif

!	This code causes a bomb if NVAR and NEQ are not set
!	In that case, call PRINTT(0) and this code won't be executed
	if(iPrint.eq.1)then
      		write(12,*)
      		write(12,*)'--MONITOR PARAMETERS-- NSTEP= ',NSTEP,'     IMASS=',IMASS
      		write(12,9005)(MON(I),VN1(MON(I)),VN2(MON(I)),I=1,NVAR-NEQ)
9005  		FORMAT('    ',16(2X,I3,1X,A2,A4))
      		write(12,9006)(DELTAX(I)*NSTEP,I=1,NVAR-NEQ)
9006  		FORMAT(' ',16F12.3)
      		write(12,*)
		endif
!      write out system components
      write(12,*)' System components'
      write(12,4321)(CoName(i),i=1,NC)
4321  format(12A8)
	if(imass.eq.1)then
		write(12,*)' Mass balance = ON'
	else
		write(12,*)' Mass balance = OFF'
	endif
	
      write(12,*)

      write(12,*)' # OF PHASES, numPh= ',numPh
      write(12,*)'   Mineral compositions'
!      write(12,*)'            #phco Moles   dM(i)  sumdM(i) vol(i)                    FRACTL'
      write(12,*)'    PhCo          xPhCo(i)    dX(i)  sumdX(i) '

      L=TANDP+NX
      LL=TANDP
!	pause 'pause 1'
!     write out information on each phase
!      DO 390 k = 1,numPh
	Do 390 kCur = 1,numPh
	k = asmCurrent(kCur)
!         write(12,*)        ! put in a blank line
         L=L+1
         write(12,200)MINREC(K),PHNAME(K)
200      FORMAT(' ',I4,1x,A32)
         do 550 j=1,numPhCo(k)
         DX=0
         SDX=0
         IF(numPhCo(K).eq.1)go to 7004
            if(j.eq.1)THEN
             do 551 I=1,numPhCo(K)-1       ! this loop for dependent component
               DX = DX-(ALLX(1,LL+i)-ALLX(5,LL+i))
               SDX=SDX-(ALLX(1,LL+i)-ALLX(2,LL+i))
551          CONTINUE
            ELSE
             LL=LL+1
             DX = ALLX(1,LL)-ALLX(5,LL)
             SDX= ALLX(1,LL)-ALLX(2,LL)
            ENDIF
7004        CONTINUE
            write(12,501)phCoName(k,j),xPhCo(k,j),dX,SDX
501         FORMAT(' ',8x,A8,3x,3F12.6)
550      CONTINUE

      	if(numSiteAtom(K).gt.0)then
      		call XtoSite(K)
      		write(12,706)(SiteAtomName(k,Lsite),Lsite = 1,numSiteAtom(K))
706   		format(4x,20A12)
      		write(12,703)(SiteAtom(k,Lsite),Lsite = 1,numSiteAtom(K))
703   		format (20F12.6)
      		endif

390   CONTINUE
	!pause 'pause 2'

      if(numNew.gt.0)then
      		write(12,*)'  '
      		write(12,*)'  '
      		write(12,*)' Values of new variables'
!       	ii = TandP + NX
!       	if(imass.eq.1)ii = ii + numPh
!               DX = DX-(ALLX(1,LL+i)-ALLX(5,LL+i))
!               SDX=SDX-(ALLX(1,LL+i)-ALLX(2,LL+i))
!      	do 750 i = ii+1,numNew+ii
!      	write(12,751)VN1(I),VN2(I),ALLX(1,i ),(ALLX(1,i)-ALLX(5,i)),(ALLX(1,i)-ALLX(2,i))
! 751   	FORMAT(A2,A16,6F12.4)
! 750   	continue
      		do 750 i = 1,numNew
		ii = newVarPtr + i - 1
      		write(12,751)VN1(ii),VN2(ii),ALLX(1,ii),(ALLX(1,ii)-ALLX(5,ii)),(ALLX(1,ii)-ALLX(2,ii))
751   		FORMAT(A2,A16,6F12.4)
750   		continue
      endif

	!pause 'pause 3'
!      write out phase moles, volumes and bulk composition
      if(imass.eq.1)then

! 	calculate mode
	sum = 0
!	do 423 k = 1,numPh
	Do 423 kCur = 1,numPh
	k = asmCurrent(kCur)
	sum = sum + vp0(k)
423	continue
!	do 422 k = 1,numPh
	Do 422 kCur = 1,numPh
	k = asmCurrent(kCur)
	mode(k) = vp0(k)*100./sum
422	continue

	write(12,*)
	write(12,*)
	write(12,*)' Moles and volumes of phases'
!       	write(12,*)'   Mineral                  mMoles   dM(i)  sumdM(i)'
      	write(12,414)
414	format(T30,'switch',T60,'Volumes')
      	write(12,413)
413	format(T6,'Mineral',T28,'mMoles',T36,'fXliz',T45,'Mode',T55,'V0',T69,'V1',T83,'V2',T97,'mMoles change')
!	do 410 k = 1,numPh
	Do 410 kCur = 1,numPh
	k = asmCurrent(kCur)
	L = molesPhPtr + k - 1
! 	write(12,411)MINREC(K),PHNAME(K),mp0(k),ALLX(1,L)-ALLX(5,L),ALLX(1,L)-ALLX(2,L)
! 411      FORMAT(I4,1x,A16,5x,3F14.3)
	write(12,411)MINREC(K),PHNAME(K),mp0(k)*1000.,fractl(k),mode(k),vp0(k),vp1(k),vp2(k),1000*(ALLX(1,L)-ALLX(5,L))
411      FORMAT(I4,1x,A16,1x,F12.3,I4,F12.3,3E14.5,F12.3)
410	continue

! 	write(12,*)
! 	write(12,*)
! 	write(12,*)' Volumes of phases and modes'
!       	write(12,*)'   Mineral                     Mode        V0(k)     v1PhCo(k)       v2PhCo(k)    fractl'
! 	do 420 k = 1,numPh
! 	write(12,421)MINREC(K),PHNAME(K),mode(k),vp0(k),vp1(k),vp2(k),FRACTL(k)
! 421	FORMAT(I4,1x,A16,5x,F12.3,3E14.5,I5)
! 420	continue

	write(12,*)'  '
	write(12,*)'  '
	write(12,*)' Bulk composition - calculated from mineral composition and moles'
	write(12,600)(coname(i),i=1,nc)
600	format(10x,24a8)
!       compute and print delta moles of system components
	do 605 i=1,nc
605   	moles(i)=0
      	j=tandp + nx + numPh
      	do 606 i=1,iopen
      	ii=i+j
      	moles(iopena(i))=ALLX(1,ii)-ALLX(2,ii)
606   	continue
      	write(12,603)(1000.*moles(i),i=1,nc)
603   	format('dmMoles   ',24F8.3)
!      	compute and print moles and weight percent of system components
      	call bulkcomp()
      	write(12,601)(1000.*moles(i),i=1,nc)
601   	format(' mMoles   ',24F8.3)
      	write(12,602)(wtpct(i),i=1,nc)
602   	format(' Wt%     ',24F8.3)

	if(bulkCompSwitch)then
		write(12,*)' '
		write(12,*)'Bulk composition = ', bulkCompTitle
		write(12,*)' Bulk composition - input wt% and moles'
		write(12,600)(coname(i),i=1,nc)
	      	write(12,601)(1000.*bulkCompMoles(i),i=1,nc)
	      	write(12,602)(bulkCompWt(i),i=1,nc)
	endif
      endif	! end imass = 1

!      write out a few blank lines to improve readability
      write(12,*)' '
      write(12,*)'  '
      write(12,*)'  '

	if(iLong(8).eq.1)call PrintMinComp	! print mineral comps in wt%

	if(ilong(9).eq.1)call writeAsmToMineralFile

      return
      end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PRINTTMoles()
      implicit none
!      ROUTINE TO PRINT molesonly
! *****************************************
	include "Assemb.inc"
	include "Monit.inc"
 	include "Output.inc"
! *****************************************
!      local variables
      integer*4 k,L,IZERO,kCur
      real*4 sum
!       COMPUTE MOLAR VOLUME OF PHASE FOR NEW COMPOSITION, T, P AND PF
            izero=0
            CALL AllKduTPX(1,izero)
!      note: vp0, vp1 and vp2 are volumes of phases based on the moles present
!            these are not modes!!
! 	calculate mode
	sum = 0
!	do 423 k = 1,numPh
	Do 423 kCur = 1,numPh
	k = asmCurrent(kCur)
	sum = sum + vp0(k)
423	continue
!	do 422 k = 1,numPh
	Do 422 kCur = 1,numPh
	k = asmCurrent(kCur)
	mode(k) = vp0(k)*100./sum
422	continue
	write(12,*)
	write(12,*)
	write(12,*)' Moles and volumes of phases'
      	write(12,414)
414	format(T30,'switch',T60,'Volumes')
      	write(12,413)
413	format(T6,'Mineral',T28,'mMoles',T36,'fXliz',T45,'Mode',T55,'V0',T69,'V1',T83,'V2',T97,'mMoles change')
!	do 410 k = 1,numPh
	Do 410 kCur = 1,numPh
	k = asmCurrent(kCur)
	L = molesPhPtr + k - 1
	write(12,411)MINREC(K),PHNAME(K),mp0(k)*1000.,fractl(k),mode(k),vp0(k),vp1(k),vp2(k)
411      FORMAT(I4,1x,A16,1x,F12.3,I4,F12.3,3E14.5,F12.3)
410	continue
!      write out a few blank lines to improve readability
      write(12,*)' '
      write(12,*)'  '
      return
      end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE WriteAsmToMineralFile
      implicit none
!      ROUTINE TO to write all phases in this assemblage to a file
!	order of elements is the same as in the Thermodynamic data file being used
! *****************************************
	include "Assemb.inc"
! *****************************************
!      local variables
      integer*4 i,k,kCur
	real*4 compositionToWrite(40)


	call PhaseComp
	! fill in an array with the composition to write
	! First zero out all entries

!	do 30 K = 1,numPh
	Do 30 kCur = 1,numPh
	k = asmCurrent(kCur)
	
	do 10 i = 1,noSysCoInFile
	compositionToWrite(i) = 0.0
10	continue

	do 20 i = 1,NC
	compositionToWrite(NOEL(i)) = PhComp(K,i)
20	continue
	!now write record
	write(16,1010)phName(k),(compositionToWrite(i),i=1,noSysCoInFile)
1010	format(A16,5x,40F10.6)

30	continue



      return
      end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE PrintMinComp
      implicit none
!      ROUTINE TO PRINT Wt Pct mineral compositions
! *****************************************
	include "Assemb.inc"
! *****************************************
!      local variables
      integer*4 i,k

	call PhaseComp

        write(12,*)'***********************************************'  
	write(12,*)' Weight Percent Oxide'
	write(12,101)(PhName(asmCurrent(k)),k = 1,numPh)
101	format(15x,24A12)
	do 10 i = 1,nc
	write(12,102)CoName(i),(PhWtOx(asmCurrent(k),i),k = 1,numPh)
102	format(A10,24f12.4)
10	continue

        write(12,*)'***********************************************'  
	write(12,*)' Weight Percent Element'
	write(12,101)(PhName(asmCurrent(k)),k = 1,numPh)
	do 20 i = 1,nc
	write(12,102)CoName(i),(PhWtEl(asmCurrent(k),i),k = 1,numPh)
20	continue

      return
      end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine PrintMinAssemblage
	
! *****************************************
	include "Assemb.inc"
 	include "Gibbsfiles.inc"
! *****************************************

! 	write the specifics of this new problem to output window #12 
	write(12,*)
	write(12,*)
        write(12,*)'###########################################################'
        write(12,*)'###########################################################'
        Write(12,*)' New Input file'
!       call timdat(ppp)       			! puts time and date into array ppp
        write(12,*)' Input file name: ',FILEIN
        write(12,*)'Assemblage title: ',Trim(AssembTitle)
        write(12,2045)
2045    format(' ')

	call printt(0)        			! option 1 prints to the file and screen

	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE NewtonStepSave()
      implicit none
!      ROUTINE TO PRINT CURRENT DATA to a file - shorthand
! *****************************************
	include "Assemb.inc"
! *****************************************
!      local variables
      integer*4 k,j,kCur

!	do 10 k = 1,numPh
	Do 10 kCur = 1,numPh
	k = asmCurrent(kCur)

	IF (MINREC(k).EQ.30) THEN				! this only writes out for garnet = 30
		write(43,*)(xPhCo(k,j),j=1,numPhCo(K)),vp1(k)
		go to 99
		ENDIF
10	continue
99	continue
	return
	end
		
