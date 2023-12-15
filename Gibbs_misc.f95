! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine Variance
      implicit none
!     Subroutine to calculate the thermodynamic variance
! ****************************************
 	include "Assemb.inc"
! ****************************************

      	gibbsVariance = TANDP + NX - NRX
      	duhemsVariance = nvar - neq
		      write(*,*)' Gibbs variance  = ',gibbsVariance
      	if(imass.eq.1)Write(*,*)' Duhems variance = ',duhemsVariance

	return
	
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine ZeroiLong
      implicit none
!     Subroutine to zero the iLong array (output options
	include "Output.inc"
      	ilong(1)=0
      	ilong(2)=0
      	ilong(3)=0
      	ilong(4)=0
      	ilong(5)=0
      	ilong(6)=0
		iLong(7)=0
		iLong(8)=0
		iLong(9)=0
		iLong(10)=0
		ilong(11)=0
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine FluidPresent
      implicit none
!     Subroutine to determine if one of the minerals is a fluid
! 	returns value of KFlu (in Assemb.inc)
! 	KFlu = 0 if no fluid is present
! 	KFlu = mineral number in list if a  fluid is present
! ****************************************
	include "Assemb.inc"
! ****************************************
	integer*4 K,kCur
!--------------------------------------------------------
! 	check to see if a fluid is present from the fluid numbers list
! 	Fluid numbers in data base are
! 	2	H2O only (WHAAR equation)
! 	3	Mixed H2O-CO2 fluid (Kerrick and Jacobs w/WHAAR for H2O)
! 	4	Mixed H2O-CO2 fluid (not working yet)

!       Check to see if a fluid is really in data file
	Do 10 kCur = 1,numPh
	k = asmCurrent(kCur)
	IF(MINREC(K).EQ.2.or.Minrec(K).eq.3.or.Minrec(K).eq.4.or.Minrec(K).eq.2002)GO TO 20
10	continue
! 	A fluid is not present
! 	call FSS_alert('ALERT!!','A fluid is not present in the data file.')
	KFlu = 0
	return

! 	if here, then a fluid is present.  Go to appropriate subroutine to calculate volume
20	continue
	KFlu = K
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE hpload
      
!     subroutine to load Hydrostatic fluid pressure data

      implicit none
	include "Gibbsfiles.inc"
!
      real*4 h
      integer*4 j,i
! ***********/HPfluid/**************************************************
      REAL*4 HP(200,24)
      common /HPfluid/ HP
! ********************************************************************
! 	HPFile = 'HydrostaticP.table'
      open(31,file=HydroPFile,status='old')
	read(31,*)		! read title
	read(31,*)		! read header
	
      do 3 j=1,200
        read(31,*)h,(HP(j,i),i=1,24)
3	continue
	close(31)
      return
      end



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE HPinterpolate(h,tc,FluidPressure,izero)

!     subroutine to extract Hydrostatic fluid pressure from the table

      implicit none
!     local variables
      integer*4 i,ii,j,jj,izero
	real*4 extraH,extraT,h,tc,fluidpressure
	character*255 ppp
! ***********/HPfluid/**************************************************
      REAL*4 HP(200,24)
      common /HPfluid/ HP
! ********************************************************************
! ------- h between 500 and 100,000 meters -------------
! ------- TC between 50 and 1200. deg !
!
      if(h.GE.500..AND.h.LE.100000..AND.TC.GE.50..AND.TC.LE.1200.)then
!     get index for Tlow (the T just below the actual Tc of the point
      i= int((tc+.01)/50.)
      ii=i+1
!     protect against array out of bounds if TC=1000.
      if(ii.gt.24)ii=24
	extraT = TC - float(i*50)
!     get index for Hlow (the H just below the actual H of the point
      j= int((h+.01)/500.)
      jj=j+1
!     protect against array out of bounds if H=100000.
      if(jj.gt.200)jj=200
	extraH = H - float(j*500)


	FluidPressure = HP(j,i) + ((HP(jj,i)-HP(j,i))/500.)*extraH + ((HP(j,ii)-HP(j,i))/50.)*extraT
	izero = 0
	else
!     if here then TC and PB are out of bounds

!      WRITE(*,*)' Tc or H are out of range for Hydrostatic Fluid Pressure table'
!      WRITE(*,*)' Tc and H equal: ',tc,h
      WRITE(ppp,*)' Tc or H are out of range for Hydrostatic Fluid Pressure table.   Tc and H equal: ',tc,h
	call fss_alert('ALERT!!',ppp)
!     set flag
      izero=1
      endif

      return

      end
	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine GetFluidDensity(izero)
      implicit none
!     Subroutine to determine the density of the fluid
! 	Assumes the fluid is mineral number KFlu (from assemblage list)

! ****************************************
	include "Assemb.inc"
! ****************************************
	Integer*4 izero
	real*4 h4,tc4,fp4
	real*8 h
!--------------------------------------------------------

! 	check to see if a fluid is present from the fluid numbers list
! 	Fluid numbers in data base are
! 	2	H2O only (WHAAR equation)
! 	3	Mixed H2O-CO2 fluid (Kerrick and Jacobs w/WHAAR for H2O)
! 	4	Mixed H2O-CO2 fluid (not working yet)

	izero = 0


! 	calculate depth (in meters) from rock density and P
	h = 100.0d0*PB/(9.8d0*RockDensity)
	h4 = h
	Tc4 = TC
	call HPinterpolate(h4,tc4,FP4,izero)
	if(izero.eq.1)return
	PFluid = FP4
	FluidDensity = 100.0D0*PFluid/(9.8d0*h)		!units are g/cm^3
	
	return
	end



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine PhaseComp
      implicit none
!     Subroutine to determine the phase composition 
! ****************************************
	include "Assemb.inc"
! ****************************************

!     local variables
      integer*4 i,j,k,kCur
      real*4 sum
! *******************************************************************
!     Figure out the compositions of the minerals
	Do 431 kCur = 1,numPh
	k = asmCurrent(kCur)
      do 432 i = 1,NC
      PhComp(K,i) = 0.
      do 430 j=1,numPhCo(K)
      PhComp(K,i) = PhComp(K,i) + xPhCo(k,j)*comp(k,j,i)
430   continue
432   continue
431   continue

	Do 150 kCur = 1,numPh
	k = asmCurrent(kCur)
!     convert moles to wt % oxides
      sum=0.D0
      do 120 i=1,nc
      PhWtOx(K,i)=PhComp(K,i)*molwt(i)/NumCatInOxide(i)
120   sum=sum + PhWtOx(k,i)
      if(sum.le.0.)go to 140
       sum=100.D0/sum
      do 130 i=1,nc     
 	PhWtOx(k,i)=PhWtOx(k,i)*sum
	PhWtEl(k,i)=PhWtOx(k,i)*OxToElWt(i)
130	continue
140   continue
150	continue
	return
	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine PhaseCompALL
      implicit none
!     Subroutine to determine the phase composition 
! ****************************************
	include "Assemb.inc"
! ****************************************

!     local variables
      integer*4 i,j,k,kCur
      real*4 sum
! *******************************************************************
!     Figure out the compositions of the minerals
	Do 431 kCur = 1,numPhMIF
	k = asmCurrent(kCur)
      do 432 i = 1,NC
      PhComp(K,i) = 0.
      do 430 j=1,numPhCo(K)
      PhComp(K,i) = PhComp(K,i) + xPhCo(k,j)*comp(k,j,i)
430   continue
432   continue
431   continue

	Do 150 kCur = 1,numPhMIF
	k = asmCurrent(kCur)
!     convert moles to wt % oxides
      sum=0.D0
      do 120 i=1,nc
      PhWtOx(K,i)=PhComp(K,i)*molwt(i)/NumCatInOxide(i)
120   sum=sum + PhWtOx(k,i)
      if(sum.le.0.)go to 140
       sum=100.D0/sum
      do 130 i=1,nc     
 	PhWtOx(k,i)=PhWtOx(k,i)*sum
	PhWtEl(k,i)=PhWtOx(k,i)*OxToElWt(i)
130	continue
140   continue
150	continue
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine CheckMassBal
      implicit none
!     Subroutine to check that mass balance is still within tolerance
! ****************************************
	include "Assemb.inc"
	include "Output.inc"
! ****************************************

!     local variables
      integer*4 i,j,k,kCur
!      real*8 pcterror(30)
! *******************************************************************
!     Compute number of moles of each system component
	do 100 i=1,nc
	moles(i)=0.D0
	Do 110 kCur = 1,numPh
	k = asmCurrent(kCur)
	do 115 j=1,numPhCo(K)
	moles(i) = moles(i) + MP0(K)*xPhCo(k,j)*comp(k,j,i)
115	continue
110	continue
! 	Check difference from starting value
	go to 100
!	pctError(i) = Dabs(   (moles(i)-molesStart(i))/molesStart(i)  )
!	if(pctError(i).gt.0.1)then
!		write(12,*)' element ',coname(i),moles(i),pctError(i)
!		if(masserroralert.eq.0)then
!	call fss_alert('Mass balance error is greater than 10%.Try using smaller step sizes. (This message will not appear again.)')
!			masserroralert = 1
!			endif
!		endif
100	continue
	return
	END



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     
      Subroutine BULKCOMP ()
      implicit none
!     Subroutine to determine the bulk composition of a rock given the 
!     mineral compositions and modes
! ****************************************
	include "Assemb.inc"
! ****************************************

!     local variables
      integer*4 i,j,k,kCur
      real*8 sum
! *******************************************************************
!     Compute number of moles of each system component
      do 100 i=1,nc
      wtpct(i)=0.D0
      moles(i)=0.D0
	Do 100 kCur = 1,numPh
	k = asmCurrent(kCur)
      do 100 j=1,numPhCo(K)
100   moles(i) = moles(i) + MP0(K)*xPhCo(k,j)*comp(k,j,i)
!     convert moles to wt % oxides
      sum=0.D0
      do 120 i=1,nc
      wtpct(i)=moles(i)*molwt(i)/NumCatInOxide(i)
120   sum=sum + wtpct(i)
      if(sum.le.0.)go to 140
      sum=100.D0/sum
! 	sum = 1.D0		! don't normalize
      do 130 i=1,nc     
130   wtpct(i)=wtpct(i)*sum
140   continue
	return
      END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE names
      implicit none
!     THIS SUBROUTINE IS DESIGNED TO SET UP variable names

! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
! ****************************************

!     Local variables
      integer*4 i,j,k,L,kCur
      CHARACTER VNAME(15)*2,BLANK*4
      common /local3/vname,blank

      BLANK='    '
!      VNAME(1)='T'
!      vname(2)='P'
!      vname(3)='u'
!      vname(4)='X'
!      vname(5)='M'
!      vname(6)='m'
!      vname(7)='1'
!      vname(8)='2'
!      vname(9)='3'

      VNAME(1)='T_'
      vname(2)='P_'
      vname(3)='u_'
      vname(4)='X_'
      vname(5)='M_'
      vname(6)='m_'
      vname(7)='1_'
      vname(8)='2_'
      vname(9)='3_'

!     THIS SECTION SETS UP THE VARIABLE NAMES

      VN1(1)=VNAME(1)
      VN2(1)=BLANK
      VN1(2)=VNAME(2)
      VN2(2)=BLANK
      IF(KFLU.NE.0)THEN
         VN1(3)=VNAME(2)
         VN2(3)=PHNAME(KFLU)
      ENDIF
5006  L=TANDP
	Do 5004 kCur = 1,numPh
	k = asmCurrent(kCur)
            DO 5005 J=2,numPhCo(K)
             L=L+1
             VN1(L)=VNAME(4)
5005         VN2(L)=phCoName(k,J)
5004  CONTINUE
      IF (IMASS.EQ.1) then
	 molesPhPtr = L + 1
	Do 5001 kCur = 1,numPh
	k = asmCurrent(kCur)
         L=L+1
         VN1(L)=VNAME(5)
5001     VN2(L)=PHNAME(K)
         ENDIF
!     names of open system components
      if(iopen.ne.0) then
	    openPtr = L + 1
            do 5100 i=1,iopen
            L=L+1
            vn1(L)=vname(6)
            vn2(L)=coname(iOpena(i))
5100        continue
      endif
!     NAMES  of new variables
      newVarPtr = L + 1			! pointer to new variables in AllX array - KEEP!
	do 20 kCur = 1,numPh
	k = asmCurrent(kCur)
	if(includeNew(k).eq.1)then
		do 21 i = 1,numNewInPhK(k)
	!      do 20 I = 1,numNew
	        L=L+1
!       	 jj = Jnew(i)
!       	  VN1(L)=' '
        	VN1(L)='__'
        	VN2(L)=newVarName(k,i)
!       	 VN2(L)=NewName(Jnew(i))
21		continue
		endif
20    continue
! 	now place dashes in all unused names, so they don't show up by accident
	do 30 i = L+1,50
	VN1(i) = '--'
	vn2(i) = '----'
30	continue
      return
      end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE RESET(IRESET,ISUB)
    	use AWE_Interfaces
	use MyCanvas
      implicit none
!     SUBROUTINE TO RESET THE MINERAL COMPOSITIONS TO THEIR STARTING
!     VALUES
! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
! ****************************************

!     local variables
      integer*4 j,k,L,izero,ireset,isub,kCur
      real*4 XXX,YYY
      character*128 ppp
!
!
!
!     VALUES FOR IRESET
!        2 = STARTING CONDITIONS
!        3 = Reference point for Contour routine
!        4 = REFERENCE POINT-user selected
!        5 = PREVIOUS finite difference POINT
!        6 = IN sub CONTOUR ( go to 1/2 of last finite
!             difference step)

!                  allx(6,L) not used

!     Values of ISUB specify from where this subroutine was called
!     ISUB = 1  called from STEP
!     ISUB = 2  called from CONTOUR
!     ISUB = 3  called from TERNARY
!     ISUB = 4  called from CONTOUR
!
!
! --------------------------------------------------
!
	DO 1490 L=1,nvar
	IF(IRESET.EQ.6)THEN
		ALLX(1,L)=(ALLX(1,L)+ALLX(5,L))/2
		ELSE
		ALLX(1,L)=ALLX(IRESET,L)
		ALLX(5,L)=ALLX(1,L)     !set previous finite difference point equal to current point
		ENDIF
1490	CONTINUE
	L=TANDP
	Do 25 kCur = 1,numPh
	k = asmCurrent(kCur)
	IF(numPhCo(K).EQ.1)GO TO 424
	xPhCo(k,1)=1
	DO 24 J=2,numPhCo(K)
	L=L+1
	xPhCo(k,J)=ALLX(1,L)
	xPhCo(k,1)=xPhCo(k,1)-xPhCo(k,J)
24	CONTINUE
424	CONTINUE
25	CONTINUE
	!     RESET MINERAL ABUNDANCE MOLES
	IF(IMASS.EQ.0) GO TO 30
	L=TANDP+NX
	Do 35 kCur = 1,numPh
	k = asmCurrent(kCur)
	L=L+1
	MP0(K)=ALLX(1,L)
	MP1(K)=MP0(K)
	MP2(K)=MP0(K)
35	CONTINUE
30	CONTINUE
!
	TC=ALLX(1,1)
	PB=ALLX(1,2)
	IF(KFLU.NE.0)PFLUID=ALLX(1,3)

1000	continue
	if(isub.eq.3)go to 1020
	IF(IRESET.NE.6)THEN
	!     calls to plot only to move plotter pen back to starting conditions
		xxx=ALLX(1,NXPLT)
		YYY=ALLX(1,NYPLT)
		IF (NYPLT.EQ.2)yyy=YYY/1000.
		call symb(XYPlot,XXX,YYY,1.,1.)    !draw a symbol on the screen just to show we got there
		ENDIF
1020		continue

	if(imass.eq.1.and.ireset.ne.6)then
	!        PHASE VOLUMES RECOMPUTED HERE
		izero=0
		call AllKduTPX(1,izero)
	!        returns in J/bar; display in cc/mole
	!        moles per cm**3 (note VMOL is in joule/bar =CC/10)
	!        note: if this code is changed, it must also be changed in routine CHANGE
		Do 300 kCur = 1,numPh
		k = asmCurrent(kCur)
		vp0(k)=mp0(K)*vmol(k)*10.D0
		vp1(k)=vp0(k)
		vp2(k)=vp0(k)
		if(FRACTL(k).eq.1.OR.FRACTL(k).eq.2)then
			ppp = 'Warning: The phase '//Trim(phname(K))//' is fractionating. Volumes and moles may not be reset correctly.'
			endif
300		continue
		masserroralert = 0  	!switch for massbalance error alert (reset to start conditions)
		endif

	RETURN
	END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine SavRxn
      
!     Subroutine to save to a file a digitized reaction curve
      implicit none
! ****************************************
	include "TPSave.inc"
! ****************************************


      integer*4 iunit,iopt,i
      character*32 flname
      character*32 fspearhead
      character*64 rxnname
      character*32 prompt

      integer*4 iopen
      data iopen /0/
      data iunit /76/

      data fspearhead /'FSpear digitized data file      '/
! ************main menu options******************************

	save
      if(iopen.eq.0)then      !a data file is not yet opened
1           continue
            write(*,*)'Data file is not yet opened.  Make a choice'
            write(*,*)'0 = return'
            write(*,*)'1 = open new file'
            write(*,*)'2 = open old file (for append)'
	    read(*,*)(iopt)
            if(iopt.eq.0)return

            if(iopt.eq.1)then
                prompt = 'Input new file name'
                flname = 'T-P output.dig'
		open(iunit,File='',status='NEW')
                write(iunit,*)FSpearhead
                iopen=1
 
 	        elseif(iopt.eq.2)then
      	        OPEN(IUNIT,FILE='',STATUS='OLD')
                do
                read(iunit,*,end=9)
                repeat
 9              continue
                Backspace(iunit)
                iopen=1

                else
                go to 1  !only if user typed in bad option number
                endif

          endif

      write(*,*)'Input a title for this reaction'
      read(*,'(A64)')Rxnname
      write(iunit,'(A64)')Rxnname
      do 100 i = 1,TPcounter
      write(iunit,1000)i,TPSave(i,1),TPSave(i,2)
100   continue
1000  format(I5,2f12.2)
      write(iunit,1000)0,0,0

      return
      end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SetTPX
      implicit none
!
! 	This routine was part of old Subroutine Compute and was moved to
! 	here when Newton's method was introduced.  
! 	The purpose is to take the new values of variables as stored in ALLX
! 	and place them into the correct compositional arrays, P and T.
!
! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
 	include "Newton.inc"
! ****************************************
!     local variables
      integer*4 j,k,L,kCur
      real*8 DELX,SUM


!
!     Set NEW P AND T
      TC=ALLX(1,1)
      TK=TC+273.15D0
      PB=ALLX(1,2)
      IF(KFLU.NE.0)then
      	PFLUID=ALLX(1,3)
	else
	PFluid = PB
	endif
!     set new values of mineral composition in phase component array X
      L=TANDP                 ! Start counting ALLX array at first independent component
	Do 3600 kCur = 1,numPh
	k = asmCurrent(kCur)
      sum=0.0d0
      DO 3601 J = 2,numPhCo(K) ! loop from second to last phase component for this mineral
      L=L+1
      xPhCo(k,j)=ALLX(1,L)
      sum=sum+xPhCo(k,j)
3601  continue
      xPhCo(k,1)=1.0D0-sum
3600  continue


!     compute new mole fractions of phases
      IF(IMASS.EQ.0)GO TO 3620
!        COMPUTE NEW Mi FOR EACH PHASE
         L=TANDP+NX
	Do 3607 kCur = 1,numPh
	k = asmCurrent(kCur)
            L=L+1

	    DELX  =ALLX(1,L) - ALLX(5,L)
	    select case(FRACTL(k))
	    case (0) 	! no fractional crystallization 
            	MP0(K)=MP0(K)+DELX
            	mp1(k)=mp1(k)+delx
            	mp2(k)=mp2(k)+delx
	    case (1)	! Fraction of phase, but does not leave system
		if(bulkCompSwitch)then
!	      If(iNewton.eq.1)then
	      	write(*,*)'^^^^^^^^^^^^^^^^^^'
!	      	write(*,*)'Fractional crystallization does not work with Newton=on'
	      	write(*,*)'Fractional crystallization does not work with fixed bulk composition'
	      	pause 'hit return to continue'
	      	endif
	    	mp0(k) = 0.0d0
		ALLX(1,L) = 0.0d0
            	mp1(k)=mp1(k)+delx
            	mp2(k)=mp2(k)+delx
	    case(2)	! Fractionats and leaves the system
		if(bulkCompSwitch)then
!	      If(iNewton.eq.1)then
	      	write(*,*)'^^^^^^^^^^^^^^^^^^'
!	      	write(*,*)'Fractional crystallization does not work with Newton=on'
	      	write(*,*)'Fractional crystallization does not work with fixed bulk composition'
	      	pause 'hit return to continue'
	      	endif
		mp0(k)=0.D0
		mp1(k)=0.D0
		ALLX(1,L) = 0.0d0
            	mp2(k)=mp2(k)+delx
	   case default
	   end select
!	Code moved from PRINTT on 10/10/99. Look out for any problems this may cause.
            VP0(K)=MP0(K)*VMOL(K)*10.D0   !note that molar volume of phase (VMOL) is in J/mol
            vp1(k)=mp1(k)*vmol(k)*10.D0 ! and VP0, VP1, and VP2 are in cm^3
            vp2(k)=mp2(k)*vmol(k)*10.D0
	    
3607     CONTINUE
!
!        Phases that are undergoing fractonal crystallization have their
!        molar quantities set to 0
!        FRACTL(k)=0 phase does not fractionate
!        FRACTL(k)=1 phase fractionates, but does not leave system
!            (so that mp1 is computed including this phase)
!        FRACTL(k)=2 phase fractionates, and leaves the system
!           (so that mp0 and mp1 are computed without this phase)
!        mp0(k) stores moles of phases that do not fractionate (case 0)
!        mp1(k) stores moles of phases for case 0 and 1
!        mp2(k) stores moles of all phases, regardless of type

3620  CONTINUE


      RETURN
      END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE Solve_for_H(ThmoFile)
	use MatrixArrays
      implicit none
!
!     written July, 1994 by F. Spear
!     ROUTINE TO calculate enthalpies of phase components from an input dataset
!
! ****************************************
	include "Assemb.inc"
	include "Output.inc"
!	include "Solute.inc"
! ****************************************
!     local variables
      integer*4 i,ii,j,jj,k,ier,izero,ivol,numSolve,kCur,kSolve(20),jSolve(20),kk
      REAL*8 DETA,R,TEMPH,TEMPS,TEMPK
      Character*32 ThmoFile         !thermodynamic data file name
!
      DATA R/8.3144D0/
!

          call printt(0)        !option 1 prints to the file and screen

!         RECOMPUTE REACTION MATRIX

	  iLong(7) = 1
          CALL REXN
          iLong(7) = 0





!
      TK=TC+273.15D0
!
!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST
!
!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      izero=0
      ivol=0            !compute everything
      K=0

      CALL AllKduTPX(ivol,izero)
      if(izero.eq.1)return
!
!     OUTPUT THERMO DATA TO DISK
!


!      iout=9
350       CONTINUE
 	Do 390 kCur = 1,numPh
	k = asmCurrent(kCur)
          write(12,*)
          write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200       FORMAT(' ',I2,F8.1,4X,A32)
          write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
          IF(IMASS.EQ.1)then
          write(12,507)MP0(K),VP0(K)
          endif
507       FORMAT(' MOLES =',F10.5,' VOL =',F10.3)
          write(12,501)(phCoName(k,J),J=1,numPhCo(K))
501       FORMAT('          ',6(A4,6X))
          write(12,504)(xPhCo(k,J),J=1,numPhCo(K))
504       FORMAT('  COMP ',6(F10.4))
          write(12,505)(Dexp(lnAct(k,J)),J=1,numPhCo(K))
505       FORMAT(' Activ  ',6(E10.3))
          write(12,506)(lnAct(k,J),J=1,numPhCo(K))
506       FORMAT('  Ln(a) ',6(F10.5))
          write(12,514)(hPhCoZero(k,J),J=1,numPhCo(K))
514       Format(' hPhCoZero',6E15.7)
          write(12,515)(HATTP(k,J),J=1,numPhCo(K))
515       Format(' HATTP',6E15.7)
          write(12,509)(sPhCoZero(k,J),J=1,numPhCo(K))
509       Format(' SZERO',6F15.5)
          write(12,510)(SATTP(k,J),J=1,numPhCo(K))
510       Format(' SATTP',6F15.5)
          write(12,508)(vPhCoZero(k,J),J=1,numPhCo(K))
508       Format(' VZERO',6F15.5)
          write(12,513)(VATTP(k,J),J=1,numPhCo(K))
513       Format(' VATTP',6F15.5)
          write(12,516)(GATTP(k,J),J=1,numPhCo(K))
516       Format(' GATTP',6E15.7)


511       FORMAT(12E12.5)
512	continue
!         write(12,408)' ACP',(aPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' BCP',(bPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' CCP',(cPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' DCP',(dPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' ECP',(ePhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' FCP',(fPhCoCp(k,J),J=1,numPhCo(K))
!         write(12,408)' GCP',(gPhCoCp(k,J),J=1,numPhCo(K))
408       format(' ',a5,6E10.3)
551       CONTINUE
390       CONTINUE
!
!

!
!     THIS PART OF THE PROGRAM SETS UP THE MATRIX
!
!     ZERO MASTER ARRAY
!
      DO 1000 J=1,NVAR
      DO 1000 I=1,NEQ
1000  AA(I,J)=0.0D0
!
!
!     Finished with extra equations
!
!
!     Set up a data vector (Y) that contains H - TS + R T lnKeq for each reaction
!

      	ii = np-nrx      !note that reaction coefficients are stored in the end of array ARX
      	DO 1010 i=1,NRX
      	TEMPH=0.D0
      	TEMPS=0.D0
      	TEMPK=0.D0
      	jj = 0
	Do 1015 kCur = 1,numPh
	k = asmCurrent(kCur)
      	do 1015 J = 1,numPhCo(k)
!     	S reaction
 	jj = jj + 1
      	TEMPS =TEMPS + ARX(i+ii,jj)*SatTP(k,j)
!     	LnK reaction
      	TEMPK = TEMPK + ARX(i+ii,jj)*lnAct(k,j)
!     	H reaction
      	TEMPH =TEMPH + ARX(i+ii,jj)*HatTP(k,j)
1015  	continue      
      	YY(i)=(TEMPH - TK * TEMPS + R*TK*TEMPK)
1010	continue

!     Now choose the phase components to solve for H

      WRITE(*,*)'Phase components in this problem'
!	ier=0
	Do 105 kCur = 1,numPh
	k = asmCurrent(kCur)
      Do 105 j=1,numPhCo(k)      
!	ier = ier + 1
      WRITE(*,*)k,j,phCoName(k,j)
105   continue

	write(*,*)' Choose phase components to solve for H'
	write(*,*)' You must choose this many components: ',nrx
	write(*,*)'Phase-K, Component-J (2 integers)'
	do i=1,nrx
		write(*,*)' Number ',i
		read(*,*)kSolve(i),jSolve(i)
		write(*,*)'Phase, Component ',kSolve(i),jSolve(i)
		end do

!     Now construct a matrix (AA) with the appropriate coefficients for each phase component to solve.

!	pause 'this code will not work.....solve for H in Gibbs_misc.f'
      numSolve = nrx+1
      do 1105 i=1,nrx
      TempH = 0
 !	The problem is that this loop is for the pointer arrays, which no longer exist.
 !	I need to correctly index the right phase component.
	jj = 0
	do 1110 kCur = 1,numPh
	k = asmCurrent(kcur)
	do 1112 j = 1,numPhCo(k)
	jj = jj + 1
	do 1120 kk = 1,nrx
	if(kSolve(kk).eq.k.and.jSolve(kk).eq.j)then
	        A(i,kk) = -ARX(ii+i,jj)       !the A matrix coefficient is the reaction coefficient
	        TempH = TempH + ARX(ii+i,jj)*hPhCoZero(k,j)    !Subtract hPhCoZero from data vector because we need to solve for this
                                                !Note it was added in above already
		endif
!      do 1110 j = 1,np
!      do 1120 k = 1,nrx
!      if(nsolve(k).eq.j)then        !this is one of the components to solve for H
!        A(i,k) = -ARX(ii+i,j)       !the A matrix coefficient is the reaction coefficient
!        TempH = TempH + ARX(ii+i,j)*hPhCoZero(k,j)    !Subtract hPhCoZero from data vector because we need to solve for this
                                                !Note it was added in above already
!        endif
1120  	continue
1112	continue
1110  	continue
      	A(i,numSolve) = YY(i) - TempH      
1105  	continue
            
!
!     Output matrix and modified matrix, if desired
!
1320  FORMAT(' ',10E16.4)
!
      write(12,*)' '
      write(12,*)' '
      write(12,*)' MASTER MATRIX'
      write(12,1321)(phCoName(kSolve(i),jSolve(i)),i=1,nrx)
1321  format(6x,10A16)
      DO 1312 I=1,NRX
      write(12,*)'--------------------'
      write(12,1320)(A(I,J),J=1,numSolve),YY(i)
1312  continue
!
!     Find solution
!
      DETA=0.D0
      IER=0
      CALL REDUCE (NRX,numSolve,DETA,IER)
      IF (IER.EQ.1) then

      WRITE(*,*)' ************ ERROR **************************'
      WRITE(*,*)' Matrix failed to invert in SUBROUTINE REDUCE'
      WRITE(*,*)' The most likely cause is that the variables'
      WRITE(*,*)' specified as monitor parameters are not linearly'
      WRITE(*,*)' independent.  Note that only the number of'
      WRITE(*,*)' intensive variables equal to the phase rule'
      WRITE(*,*)' variance of the assemblage are independent.'
      WRITE(*,*)' Return to previous menu and try again.'
      write(*,*)' '
      izero=1
      pause 'Hit return to continue...'
      return
      endif
!
!
      write(12,*)'   '
      write(12,*)' RESULTS '
      write(12,*)'DETERMINANT=',DETA

!     Write out solution vector
      write(12,*)
      write(12,*)'Phase components        H'
      Do 1605 i=1,nrx
!      write(12,*)phCoName(k,Nsolve(i)),XX(i,1)
      write(12,*)phCoName(kSolve(i),jSolve(i)),XX(i,1)
1605  continue
      write(12,*)

      RETURN
      END
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!	subroutine timdat(iunit)
!	integer *4 MM,DD,YY,SEC,Hrs,Min
!        character*255 ppp


!	call date(MM,DD,YY)
!	Call TIME(sec)
!	Hrs=INT((sec+1)/3600)
!	sec=sec-3600*hrs
!	Min=int(SEC)/60
!	sec=sec-60*min
!	write(12,10)MM,DD,YY,Hrs,min,sec
!        call printf(ppp)
!10	format('Date: ',I2,'-',I2,'-',I2,'    Time: ',I2,':',I2,':',I2)
!	return
!	end
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE WHOLEROCK
    	use AWE_Interfaces
	use MatrixArrays
!     Subroutine WholeRock
!     Calcultes whole rock reactions based on the Jacobian
      implicit none

! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
!	include "Solute.inc"
! ****************************************
!     local variables
      integer*4 i,j,jj,k,L,LL,inorm,izero,iout,thermov,&
     & imassold,nvarold,neqold,kCur
      real*8 delT,delP
      character*32 FLNM
      REAL*8 Sphase(15),Vphase(15),SumS,SumV

! ****************************************



!     compute thermodynamic variance
!     If Thermov = 1 then we have a univariant curve and the reaction balance
!     is independent of the amounts of any phases present (only on their compositions)
!     So we will arbitrarily assign masses of 1 to each phase and proceed.
!     If Thermov>= 2 then the results depend on the amounts of each phase
!     present and we must force the user to go back and specify the amounts of
!     each phase.
!     If Thermov <1 then we just abort, because there is no code for this eventuality

      thermov = TANDP + NX - NRX

      if(thermov.le.0)then
            write(*,*)' The thermodynamic variance is ',thermov
            write(*,*)' This routine only works if the thermodynamic variance is >= 1'
            write(*,*)' Hit return to continue'
            pause
            return
            endif

! ----Compute whole wrock reaction for univariant reactions ------------------------------------------

      if(Thermov.eq.1)then          ! Reaction balancing for univariant reactions

1000   CONTINUE

      MON(1)=1
      Mon(2)=TANDP+NX+1       ! Pick first M variable as monitor

!     For univariant reactions, the whole rock reaction is independent of the
!     mass of phases, so assign an arbitrary mass of 1 to each and proceed
!     but only do this if Imass = 0 (that is, we have not already specified masses

      imassold=imass
      nvarold=nvar
      neqold=neq

      nvar=nvar+numPh
      neq=neq+nc

      if(imass.eq.0)then
	Do 1005 kCur = 1,numPh
	k = asmCurrent(kCur)
      MPSTART(K)=1
      MP0(K)=1
      MP1(k)=1
      MP2(K)=1
      VPSTART(k)=1
      VP0(k)=1
      VP1(k)=1
      VP2(k)=1
      VMOL(k)=1
1005  continue
      imass=1
      endif

      call names

      NSTEP=1
      DELTAX(1)=0
      DELTAX(2)=0


      WRITE(*,*)' *********************************'
      write(*,*)' Thermodynamic variance = ',Thermov
      WRITE(*,*)' CURRENT MONITOR PARAMETERS ARE'
      WRITE(*,5)(MON(I),VN1(MON(I)),VN2(MON(I)),I=1,NVAR-NEQ)
      WRITE(*,6)(ALLX(1,MON(I)),I=1,NVAR-NEQ)
      WRITE(*,6)(DELTAX(I)*NSTEP,I=1,NVAR-NEQ)
      WRITE(*,*)'NSTEP =',NSTEP,'  PEN IS',NOPEN,'   idraw =',idraw
      WRITE(*,*)' CURRENT AXES ARE  ',NXPLT,NYPLT
      write(*,7)VN1(NXPLT),VN2(NXPLT),VN1(NYPLT),VN2(NYPLT)
      WRITE(*,*)' *********************************'
      WRITE(*,*)' Routine to compute whole rock reactions for a univariant reaction'


      WRITE(*,*)'The procedure used is to compute the Jacobian for the assemblage at the current conditions'
      WRITE(*,*)'to get the change of moles of 1 phase relative to others.'
      WRITE(*,*)'The relative molar changes are then normalized to give a chemical reaction.'
      WRITE(*,*)'The chemical reaction is then multiplied by the S and V of each phase to give the facing.'

      write(*,*)' Hit return to continue'
      pause

      FLNM = 'WholeRockRxn.1'

!     Set up array IPOINT to contain pointers to non-monitor parameters
      J=0
      Do 1110 i=1,NVAR
      DO 1111 L=1,NVAR-NEQ
      IF(MON(L).eq.i)go to 1110
1111   continue
      J=J+1
      IPOINT(J) = I
1110   CONTINUE
      IF(iLong(4).eq.1) write(*,*)(IPOINT(I),I=1,NEQ)
!
        IZERO=0
        call c2ompute(izero)
!     now we have the Jacobian stored in array XX (from common block solute)
!     write results to console and data file
         iout=12
1510     CONTINUE
         WRITE(IOUT,*)
         WRITE(IOUT,1507)
1507     FORMAT('    Jacobian matrix:')
         WRITE(IOUT,1502)(VN1(MON(i)),VN2(MON(i)),i=1,nvar-neq)
1502     FORMAT(' Monitors:',5x,15(A2,A4,9x))
         write(iout,1506)(deltax(I),i=1,nvar-neq)
1506     format('      ',15(f15.4))
         JJ=TANDP + NX + 1      !start with the first molar variable
!     The XX array contains
!     XX(1,L)    = dP/d...
!     XX(2,L)    = dX1/d...
!     XX(NX+1,L) = dXnx/d...
!     XX(NX+2,L) = dM2/d...

         WRITE(IOUT,1501)2,VN1(2),VN2(2),(XX(1,LL),LL=1,nvar-neq)
1501        FORMAT(' ',I5,2x,A2,A4,' = ',15E15.5)

         DO 1520 i=2,NX+1
         JJ=i+1
         WRITE(IOUT,1501)jj,VN1(jj),VN2(jj),(XX(i,LL),LL=1,nvar-neq)
1520     CONTINUE

        do 1521 i = NX+2,NVAR-2
        jj=i+2
        WRITE(IOUT,1501)jj,VN1(jj),VN2(jj),(XX(i,LL),LL=1,nvar-neq)
1521  continue

!         IF(IOUT.EQ.30)GO TO 1500
!            IOUT=30
!            if(iok.eq.0)GO TO 1510
1500     CONTINUE
!
1150   continue

!      write(*,*)' Which phase or component should be used'
!      write(*,*)' for normalization (from Jacobian list)'
!      read(*,*)inorm
!1200   continue
!      write(*,*)' Specify dT and dP path for reaction (0,0 to exit)'
!      read(*,*)delT,delP
!      if(delT.eq.0.and.delp.eq.0)then
!            close(30)
!            return
!            endif

!     compute reaction normalized to phase 1 and store in AA(i,1)
      AA(1,1)= 1.0
      j=1
      do 1205 i=2+nx,neq
      j=j+1
      AA(j,1)= XX(i,2)
1205   continue

!
!     Compute Sphase and Vphase as Sphase = xPhCo(i)*iSattp(i)

      SumS=0.0D0
      SumV=0.0D0
	Do 1301 kCur = 1,numPh
	k = asmCurrent(kCur)
      Sphase(K)=0.0d0
      Vphase(K)=0.0D0
      do 1302 j = 1,numPhCo(k)
            Sphase(K)=Sphase(K)-xPhCo(k,j)*dudTPX(k,j,1)
            Vphase(K)=Vphase(K)+xPhCo(k,j)*dudTPX(k,j,2)
1302  continue

!     multiply each Mphase by Sphase or Vphase and store in AA(i,2 and AA(i,3)
      AA(K,2)=Sphase(K)
      AA(K,3)=AA(K,1)*Sphase(K)
      AA(K,4)=Vphase(K)
      AA(K,5)=AA(K,1)*Vphase(K)
!     sum S and V
      SumS = SumS+AA(K,3)
      SumV = SumV+AA(K,5)
1301  continue
   

!
!     output results
         iout=12
1610     CONTINUE
         WRITE(IOUT,*)
         WRITE(IOUT,1607)
1607     FORMAT('    Balanced reactions:')

!         WRITE(IOUT,1602)(VN1(MON(I)),VN2(MON(I)),i=1,nvar-neq)
!1602     FORMAT(' Monitors:',5x,15(A2,A4,9x))
!         write(iout,1606)delT,delP
!1606     format('      ',15(f15.4))

      write(iout,*)'    Normalized Reaction    Sphase          S          Vphase          V'

!     write out Xi terms
!         JJ=2
!         DO 1620 i=1,NP-numPh
!            JJ=JJ+1
!            WRITE(IOUT,1601)VN1(JJ),VN2(JJ),(AA(i,LL),LL=1,1)
!1620     CONTINUE

!     write out Mi terms
      jj=tandP+nx
	Do 1621 kCur = 1,numPh
	k = asmCurrent(kCur)
            JJ=JJ+1
            WRITE(IOUT,1601)VN1(JJ),VN2(JJ),(AA(k,LL),LL=1,5)
1621     CONTINUE
1601        FORMAT(' ',A2,A4,' : ',7E15.5)

      write(iout,1631)SumS, SumV
1631  format(' Sums     ',15x,15x,e15.5,15x,e15.5)



!         IF(IOUT.EQ.30)GO TO 1600
!            IOUT=30
!            if(iok.eq.0)GO TO 1610
1600     CONTINUE

      write(*,*)' Hit return to continue'
      pause


1999  continue
      close(30)
1998  continue
      imass=imassold
      nvar=nvarold
      neq=neqold
      return


      endif       ! end reaction balancing for divariant or greater reactions



! ----Compute whole wrock reaction------------------------------------------

      if(thermov.ge.2)then          ! reaction balancing for divariant or greater reactions

      MON(1)=1
      MON(2)=2

      NSTEP=1
      DELTAX(1)=0
      DELTAX(2)=0
1     CONTINUE
      WRITE(*,*)' *********************************'
      write(*,*)' Thermodynamic variance = ',Thermov
      WRITE(*,*)' CURRENT MONITOR PARAMETERS ARE'
      WRITE(*,5)(MON(I),VN1(MON(I)),VN2(MON(I)),I=1,NVAR-NEQ)
5     FORMAT(' ',3x,6(2X,I3,1X,A2,A4))
      WRITE(*,6)(ALLX(1,MON(I)),I=1,NVAR-NEQ)
      WRITE(*,6)(DELTAX(I)*NSTEP,I=1,NVAR-NEQ)
6     FORMAT(' ',6F13.4)
      WRITE(*,*)'NSTEP =',NSTEP,'  PEN IS',NOPEN,'   idraw =',idraw
      WRITE(*,*)' CURRENT AXES ARE  ',NXPLT,NYPLT
      write(*,7)VN1(NXPLT),VN2(NXPLT),VN1(NYPLT),VN2(NYPLT)
7     FORMAT(15x,2(5x,A2,A4))
      WRITE(*,*)' *********************************'
      WRITE(*,*)' Routine to compute whole rock reactions'

      if(NVAR-NEQ.NE.2.or.IMASS.NE.1.OR.TANDP.NE.2)then
      write(*,*)'Sorry, but this calculation is only possible in systems where the variance = 2 and mass balance is used (imass=1)'
           write(*,*) 'Hit return to continue'
            pause
            return
            endif

100   CONTINUE
      WRITE(*,*)'The procedure used is to compute the Jacobian for the assemblage at the current conditions.'
      WRITE(*,*)' A value of T and P is then specified and multiplied by the Jacobian to produce changes in '
      WRITE(*,*)'the composition and moles of the phases. ' 
      WRITE(*,*)' The partials of composition with T and P are multiplied by the number of'
      WRITE(*,*)'moles of each phase in the reference assemblage and the results normalized, if desired. ' 
      WRITE(*,*)'  A second option presents the reaction as a change from assemblage A1 to A2 :'
      write(*,*)'  A1 + B1 + C1 = A2 + B2 + C2 where A, B, and ! are phases'

      write(*,*)' Hit return to continue'
	pause
	
      FLNM = 'WholeRockRxn.1'

!	open(30,File='',status="NEW")

!     Set up array IPOINT to contain pointers to non-monitor parameters
      J=0
      Do 110 i=1,NVAR
      DO 111 L=1,NVAR-NEQ
      IF(MON(L).eq.i)go to 110
111   continue
      J=J+1
      IPOINT(J) = I
110   CONTINUE
      IF(iLong(4).eq.1) write(*,*)(IPOINT(I),I=1,NEQ)
!
        IZERO=0
        call c2ompute(izero)
!     now we have the Jacobian stored in array XX (from common block solute)
!     write results to console and data file
         iout=12
510      CONTINUE
         WRITE(IOUT,*)
         WRITE(IOUT,507)
507     FORMAT('    Jacobian matrix:')
         WRITE(IOUT,502)(VN1(MON(I)),VN2(MON(I)),i=1,nvar-neq)
502     FORMAT(' Monitors:',5x,15(A2,A4,9x))
         write(iout,506)(deltax(I),i=1,nvar-neq)
506     format('      ',15(f15.4))
         JJ=2
         DO 520 i=1,NEQ
            JJ=JJ+1
         WRITE(IOUT,501)i,VN1(JJ),VN2(JJ),(XX(i,LL),LL=1,nvar-neq)
501        FORMAT(' ',I5,2x,A2,A4,' = ',15E15.5)
520     CONTINUE
!         IF(IOUT.EQ.30)GO TO 500
!            IOUT=30
!            GO TO 510
500      CONTINUE
!
150   continue
      write(*,*)' Which phase or component should be used'
      write(*,*)' for normalization (from Jacobian list)'
      read(*,*)inorm
200   continue
      write(*,*)' Specify dT and dP path for reaction (0,0 to exit)'
      read(*,*)delT,delP
      if(delT.eq.0.and.delp.eq.0)then
!            close(30)
            return
            endif
!     compute reaction stoichiometry and store in AA(i,1)
      do 205 i=1,neq
      AA(i,1)=XX(i,1)*delT + XX(i,2)*delP
205   continue
!     normalize to 1.0 of phase (inorm) and store in AA(i,2)
      if(AA(inorm,1).eq.0)go to 211
      do 210 i=1,neq
      AA(i,2)=AA(i,1)/AA(inorm,1)
210   continue
211   continue
!     multiply exchange vectors by moles of the phase
      L=0
	Do 220 kCur = 1,numPh
	k = asmCurrent(kCur)
         IF(numPhCo(K).EQ.1)GO TO 220
         DO 225 J=2,numPhCo(K)
            L=L+1
            AA(L,1)=AA(L,1)*MP0(K)
225         AA(L,2)=AA(L,2)*MP0(K)
220   CONTINUE
!
!     Now compute before - after type of reaction
!     "before" is the number of moles of each phase at the starting conditions
!     "after" is the number of moles of each phase at T+T and P+P
!     store results in AA(3,i)-before and AA(4,i)-after
!     store normalized results in AA(5,i) and AA(6,i)
!     set arrays equal to 0
      do 250 i=1,neq
      AA(i,3)=0.
      AA(i,4)=0.
      AA(i,5)=0.
      AA(i,6)=0.
250   continue

!      do 230 i=1,numPh
!      K=i+np-numPh
!      AA(K,3)=MP0(i)
!      AA(K,4)=MP0(i) + AA(K,1)
!230   continue
!      if(AA(inorm,4).eq.0)go to 241
!      do 240 i=1,numPh
!      K=i+np-numPh
!      AA(k,5)=AA(k,3)/AA(inorm,4)
!      AA(k,6)=AA(k,4)/AA(inorm,4)
!240   continue
!241   continue

	Do 230 kCur = 1,numPh
	k = asmCurrent(kCur)
      i=k+np-numPh
      AA(i,3)=MP0(k)
      AA(i,4)=MP0(k) + AA(i,1)
230   continue
      if(AA(inorm,4).eq.0)go to 241
	Do 240 kCur = 1,numPh
	k = asmCurrent(kCur)
      i=k+np-numPh
      AA(i,5)=AA(i,3)/AA(inorm,4)
      AA(i,6)=AA(i,4)/AA(inorm,4)
240   continue
241   continue

!
!     output results
         iout=12
610     CONTINUE
         WRITE(IOUT,*)
         WRITE(IOUT,607)
607     FORMAT('    Balanced reactions:')
         WRITE(IOUT,602)(VN1(MON(I)),VN2(MON(I)),i=1,nvar-neq)
602     FORMAT(' Monitors:',5x,15(A2,A4,9x))
         write(iout,606)delT,delP
606     format('      ',15(f15.4))
      write(iout,*)'               EXCHANGE NOTATION           BEFORE-AFTER   NOTATION '
      write(iout,*)'            Absolute   Normalized         Absolute                         Normalized '
         JJ=2
         DO 620 i=1,NP-numPh
            JJ=JJ+1
            WRITE(IOUT,601)VN1(JJ),VN2(JJ),(AA(i,LL),LL=1,2)
620     CONTINUE
         DO 621 i=NP-numPh+1,NEQ
            JJ=JJ+1
            WRITE(IOUT,601)VN1(JJ),VN2(JJ),(AA(i,LL),LL=1,6)
621     CONTINUE
601        FORMAT(' ',A2,A4,' : ',E15.5,F10.3,2E15.5,2f10.3)
!         IF(IOUT.EQ.30)GO TO 600
!            IOUT=30
!            GO TO 610
600     CONTINUE
        go to 200

      endif       !end reaction balancing for divariant or greater reactions


      return


      END



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine XtoSite(K)
      implicit none
!     ROUTINE to compute site occupanicies of multisite minerals from mole fractions of components
! ****************************************
	include "Assemb.inc"
! ****************************************
!     local variables
      integer L,j,K

!     calculate atoms on site information
      if(numSiteAtom(K).eq.0)go to 699
      do 611 L = 1,numSiteAtom(K)
      SiteAtom(k,L) = 0
      do 610 j = 1,numPhCo(K)
      SiteAtom(k,L) = SiteAtom(k,L) + XToSiteAtom(k,L,j)*xPhCo(k,j)
!     Note that the difficult subscripting is due to:
!     (1)  X and phCoUsedinFile are linear arrays of only those phase components that are used
!            phCoUsedinFile contains integers refering to the original list of components from the thermodynamic data file
!     (2)  SiteAtom is a linear array of all the site fractions (used and unused) from the thermo file
!     (3)  XtoSiteAtom is a 2-D array
!           subscript 1 is the same as SiteAtom subscripts
!           subscript 2 refers to the original list of components in the thermo file, and hence needs
!                 the index from phCoUsedinFile
610   continue
611   continue

699   continue
      return
      end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine WriteOutJacobian()
	use MatrixArrays
	implicit none
	include "Assemb.inc"
 	include "Monit.inc"
 	include "Newton.inc"
! 	include "Solute.inc"
	integer*4 i,j,jj,ll,nsolve
	character soln1*2,soln2*4

	write(12,*)
	write(12,3507)
3507    FORMAT('    Jacobian matrix:')
	if(inewton.eq.0)then
		write(12,3502)(VN1(MON(I)),VN2(MON(I)),i=1,nvar-neq)
		nsolve = nvar
		else
		soln1 = 'So'
		soln2 = 'ln V'
		write(12,3502)(VN1(MON(I)),VN2(MON(I)),i=1,nvar-neq),soln1,soln2
		nsolve = nvar + 1
		endif		

3502    FORMAT(' Monitors:',5x,15(A2,A4,9x))
	write(12,3506)(deltax(I),i=1,nvar-neq)
3506    format('      ',15(f15.4))
	JJ=0
	DO 3520 I=1,NVAR
	DO 3523 J=1,NVAR-NEQ
3523    IF(I.EQ.MON(J))GO TO 3520
	JJ=JJ+1
	write(12,3501)VN1(I),VN2(I),(XX(JJ,LL),ll=1,nsolve-neq)	! writes out soln in last XX
3501    FORMAT(' ',A2,A4,' = ',15E15.5)
3520    CONTINUE
	return
	end

