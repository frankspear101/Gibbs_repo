	Module MyCanvas
	use AWE_Interfaces
	TYPE(AWE_Canvas) :: XYPlot		! The plotting canvas
	end Module
      PROGRAM GIBBS
	USE AWE_INTERFACES
	use MyCanvas
      implicit none
!	TYPE(AWE_Canvas) :: XYPlot

!           Coded by F.S. Spear 2/85
!           Macintosh version written by F. Spear 2/89
!           modified for contour routine by T. Menard 6-9/88
!           modified for open system behavior by F. Spear 3/89
!           modified for lnKeq equations by F. Spear 2/90
!           modified for linear arrays by F. Spear 9/90
!           modified for V=f(P,T) by M.J. Kohn 9/90
!           modified for Margules phases by M.J. Kohn 9/90
!               Kerrick and Jacobs mixed H2O-CO2 fluids
!           modified for Adobe Illustrator output by M.J. Kohn 1/91

!           Version 3.0 involves the major modifications necessary
!               to incorporate multisite mixing F. Spear (January, 1993)
!               v 3.0 also calculates Gibbs matrix by using partial
!                 molar quantities (, S and V) multiplied by
!                 coefficients matrix
!           Version 4.0 incorporates a Newton-Rapheson routine to
!               solve for P-T conditions using enthalpy data
!               plus code to solve for partial derivatives of complex
!               solutions using numerical partial derivatives
!	Version 5 is a complete rewriting of the storage of data.
!		Rather than using linear arrays and pointers (efficient storage but a coding nightmare)
!		Version 5 uses indicies of (phase,phaseComponent) - more storage but much easier to keep track of
!		This version may also eventually use pointers to the "currentAssemblage" rather than
!			loading and reloading every time a new assemblage is tried. This will be much more efficient
!			for things like pseudosections
!	November, 2016. Converted code toe FORTRAN 95 and implemented AWE interfaces
!		This allowed retiring the old Carbon interface routines
!     	Logical unit numbers used are:
!     	1 = thermodynamic data file 
!     	3 = configuration file
! 	23 = Gibbs.PlotDefinitions file (contains X-Y coordinate plot definitions)
!     	4 = RESTORE/SAVE file (*.sav)
!     	5 = input unit number in subroutine BEGIN
!     	26 = PT path input file
!     	9 = system console (on Macintosh)
!     	11 = Graphics window unit number
!     	12 = Output window unit number
!     	14 = New Input File window unit number
!     	15 = DiffGibbs output window unit number
!     	19 = DiffGibbs_D_Coef.dat file unit number
!     	22 = Output of *.in type of save file (SUBROUTINE CHANGE)
!     	51 = Postscript scratch file unit number
!     	21 = Postscript unit number (for saveas)
! c*******************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
	include "Singular.inc"
!	include "MakeMyGrid2.inc"
	include "Tangent.inc"
	include "MatrixSingular.inc"	
! c*******************************************
!     LOCAL VARIABLES
	logical exist
	integer*4 ioption,i,numBC
	character*32 ConfigTitle,rockfile(16)
	character*16 dummy
	character*64 thisThmoFile
	integer*4 skipKFlu
	common /skipswitch/skipKFlu
      	data iternary,NoTieMinerals/0,1/	! initialize tie line plotting to zero
!     	A LISTING OF THE VARIABLES USED IN THIS PROGRAM IS IN THE
!       GIBBS.TXT FILE INCLUDED WITH THE PROGRAM
!	 open configuration (preferences) file

	if(MaxPC.ne.MaxPCMonit)then		! just a check to make sure the 2 .inc files are in sync
		call FSS_Alert('Alert','MaxPC (Assemb.inc) =/ MaxPCMonit (Monit.inc) - Fix!!!')
		endif

      	CONFIG='./Gibbs_Essentials/Gibbs.FIG'
! 	Check to see that the configuration file is present.  If not, display a message
	inquire(file = config, EXIST = exist)
	if(.not.exist)then
		write(*,*)'The configuration file is not present.  Abort the mission'
		pause 'Hit return to continue'
		stop
		endif
      	open(3,file=Config,status='old')
      	read(3,'(A)') ConfigTitle
      	READ(3,'(A)')dummy          	! read the row of dashes
      	call ReadThomoFileList(thisThmoFile)
	thmoFile = thisThmoFile
	write(*,*)' ----------------'
      	READ(3,'(A)')dummy          	! read the row of stars
	read(3,*)numBC	! number of bulkrockcomposition files
	do 10 i = 1,numBC
	read(3,*)rockFile(i)			! temporary storage of file names
	write(*,*)i,rockFile(i)
10	continue
	write(*,*)'Pick the bulk composition file you wish to use'
	read(*,*)numBC
	bulkRockCompositionFile = rockFile(numBC)		! this is the file name used in Sub PickBulkComposition
	bulkRockCompositionFile = './Gibbs_Essentials/'//BulkRockCompositionFile		! this is the file name used in Sub PickBulkComposition

      	READ(3,'(A)')dummy          	! read the row of stars
      	READ(3,*)inewton,tolnewton,newtonstep
      	newtonStepMax = newtonStep
! 	inewton=1 uses Newton's method to solve integrated equations; 0 is Gibbs method                                  
! 	Tolnewton is default tolerance for convergence using Newton's method 
      	READ(3,'(A)')dummy          	! read the row of stars
!     	read plotting default values
      	Read(3,*)iBig,iTrip,iAlphaBeta     ! number for default plot in config file
      	READ(3,'(A)')dummy          	! read the row of stars
      	READ(3,'(A)')dummy          	! read A HEADER
	read(3,*)TcStart,TcEnd,Tinc
	read(3,*)PbStart,PbEnd,Pinc
! 	done reading Gibbs.fig file
      	open(23,file='./Gibbs_Essentials/Gibbs.PlotDefinitions',status='old')	! open file with plot definitions (unit 23)
      	READ(23,'(A)')dummy          	! read the header
      	do 215 i = 1,(iBig-1)*7     	! skip the ones we don't need
      	read(23,'(A)',end=218)dummy
215   	continue
      	go to 219         		! No end of file (yet)
218   	continue    			! we hit an end of file - problem!
  !    	call FSS_alert('The default plot number exceeds the number of plot def&
 !    	&initions in the configuration file.  Plot number 1 is being used')
      	rewind(23)			! start over and choose iBig = 1 as default plot
      	READ(23,'(A)')dummy          	! read the header
      	iBig = 1
219   	continue
!     	Now read the plot definition
      	READ(23,'(A)')dummy          	! read the row of stars
      	READ(23,'(a)',end=218)PlTitle
      	READ(23,*,end=218)NXPLT,NYPLT
      	READ(23,'(a)',end=218)XLAB
      	read(23,*,end=218)xor,xmin,xmax,xlen,nxstep,nxdec
      	READ(23,'(a)',end=218)YLAB
      	read(23,*,end=218)yor,ymin,ymax,ylen,nystep,nydec
!     	Open thermodynamic data file
	call OpenThermoFile
! 	Check to see that the file is present.  If not, display a message
	HydroPFile = './Gibbs_Essentials/HydrostaticP.table'
	inquire(file = HydroPFile, EXIST = exist)
	if(.not.exist)then
		call fss_alert('ALERT!!','The file (HydrostaticP.table) is not present in the Gibbs_Essentials folder.')
		stop
		else
! 		load hydrostatic fluid pressure table
		call hpload
		endif
! 	Check to see that the file is present.  If not, display a message
	DiffGibbsFile = './Gibbs_Essentials/DiffGibbs_D_Coeffs.dat'
	inquire(file = DiffGibbsFile , EXIST = exist)
	if(.not.exist)then
!		call fss_alert('The file (DiffGibbs_D_Coeffs.dat) is not present in the Gibbs_Essentials folder.&
 !    & If you try to run a DiffGibbs problem, the program will crash')
		endif

!     	initialize variables
      	icont=0
      	numNew=0
      	IOPEN=0
	iCareifMIS = 1		! i care if matrix is singular = yes (print error message)
	call ZeroiLong
      	ilong(1)=1		! error output
      	ilong(11)=1		! echo input file reading new assemblage
	newtonStepsOutput = 0			! flag to PRINTT every Newton Step (only for overstepping calculation)
	idoingpseudo = 0
	iDoingGrid = 0			! flag to alert activity calculations to stop if finite difference derivative fails at a subsystem
						! Only used in Subroutine dlnAdX in file Gibbs_Thermocalc.f95
! 	Kfluid variables
	skipKFlu = 0			! allows skipping KFluid questions from AFM routine only
	PFluidSwitch = 0
	KFlu = 0
	RockDensity = 2.7		! gm/cm^3
	bulkCompSwitch = .false.	! No bulk composition is open at start up 
!      	initialize postscript output
!      	SET DEFAULTS FOR MONITOR AND DELTAX
	ndep = 1
      	NSTEP=1
      	SNSTEP=1
      	DO 3001 I=1,16
         SMON(I)=I
         MON(I)=I
	 ncont(i) = i		! contour variable
         DELTAX(I)=0.d0
         SDEL(I)=0.d0
3001  	CONTINUE
!     	set defaults for plotting
      	idraw=1
      	ispace=' '
      	isymb=0
      	nopen=0
      	IPEN=3
! 	Open output window
	OPEN(12,FILE = 'Gibbs3 OutPut',ACCESS = 'window, 800, 800')
	call AWE_MoveWindow(12,750,0)
!      	draw initial plot with defaults
!      		CALL USER(xor,XMIN,XMAX,xlen,yor,YMIN,YMAX,ylen) -- now in Sub PlotAxes
!     		call axis(nxstep,NXDEC,XLAB,nystep,NYDEC,YLAB,PLTITL)-- now in Sub PlotAxes
	ColorFileName = './Gibbs_Essentials/GibbsColors.txt'
	call ReadColorFile		! set up AWE and PS colors
	CurrentColor = 1		! default = Black
	CurrentColor = 1		! Black for the axes

!	plotting axes on startup
!	call PlotAxes(XYPlot)


  	!if(itrip.eq.1) call Al2SiO5(3,0)       ! option 3 is Berman and Brown triple point (UBC)
  	!if(iAlphaBeta.eq.1) call abqtz(0)
! 	set default bounds for contour operations as plotting limits
	xMinBound = xmin
	xMaxBound = xmax
      	if(NYPlt.eq.2.and.ymin.eq.0.0)then
		yMinBound = 0.1			! If Y axis is P=0, set lower P bound to 100 bars because at P less than 100 bars, H2O model is no good
	else
		yMinBound = ymin
	endif
	yMaxBound = ymax
! 	set defaults for ternary plot
	Llab = 'Left'
	Rlab = 'Right'
	Tlab = 'Top'
	!Pltitle = 'Ternary plot'
	scale = 1.0
	itic = 1
	PlotIndex(1,1) = NXPLT
	PlotIndex(1,2) = NYPLT
	PlotIndex(1,3) = 3
	PlotIndex(1,4) = 4
      	isfileopen = 0
	! these are initial values for grain boundary energies.
	! They can be changed in subroutine ThermoData
	HseGB =  10000.  	! I'm not sure what these numbers should be (or even what the units are)
	HdeGB =1000000.

! ----- MAIN CONTROL MENU ----------------------------
	call AWE_activateWindow(6)
100   CONTINUE
	ioption=1         					! default
	write(*,*)' ***********************************'
	write(*,*)'Thermo file: ',THMOFILE         		! thermodynamic data file
     	If(inewton.eq.1)write(*,*)'Newton''s method'
     	If(inewton.eq.0)write(*,*)'Gibbs'' method'
	WRITE(*,*)' ***********************************'
	WRITE (*,*)' MAIN MENU OPTIONS:'
	WRITE (*,*)'   1 = Begin/save problem'
	WRITE (*,*)'  -----------------------------'
	WRITE (*,*)'   2 = Single steps'
!	WRITE (*,*)'  22 = Auto Single steps'
	WRITE (*,*)'   3 = Contour X-Y diagrams'
 	WRITE (*,*)'   4 = MAD modeling (Pseudosection)'
!	WRITE (*,*)'   4 = Make my grid ' ! - new version'
!	WRITE (*,*)'   5 = Grow Garnet '
! 	WRITE (*,*)'  55 = Grow garnet with diffusion'
!	WRITE (*,*)'   6 = DiffGibbs (Garnet growth with diffusion)'
!	WRITE (*,*)'  66 = DiffNewton (Garnet growth with diffusion)'
!	WRITE (*,*)'   7 = Whole rock reaction balancing'
	WRITE (*,*)'  -----------------------------'
	WRITE (*,*)'   8 = Go to global menu'
	WRITE (*,*)'   9 = Plotting menu'
!	WRITE (*,*)'  10 = Plot digitized reactions'
	WRITE (*,*)'  11 = Thermodynamic data menu'
	WRITE (*,*)'  -----------------------------'
 !	WRITE (*,*)'  12 = Pseudosections from Scratch'
! 	WRITE (*,*)'  13 = Forward modeling'
 !	WRITE (*,*)'  14 = Tweek calculations'
 !	WRITE (*,*)'  16 = Plot Keq curves'
 !	WRITE (*,*)'  17 = Forward modeling EBC = Nucleation density'
 !	WRITE (*,*)'  18 = Overstep garnet composition'
 !	WRITE (*,*)'  19 = Grid calculations'
	WRITE (*,*)'  -----------------------------'
	WRITE (*,*)' CHOOSE OPTION'
	read(*,*)ioption

! ----------------- 1 ---------------------
	if(ioption.eq.1)then
	call GibbsBegin
	endif
! ------------------ 2 --------------------
      IF(IOPTION.EQ.2)then
            if(isfileopen.eq.0)Then
                call fss_alert('ALERT!!','You must open an input file first')
                go to 100
                endif
       CALL STEPS
      endif
! ------------------ 22 --------------------
      IF(IOPTION.EQ.22)then
            if(isfileopen.eq.0)Then
                call fss_alert('ALERT!!','You must open an input file first')
                go to 100
                endif
!       CALL AutoSTEPS
      endif
! ------------------ 3 --------------------
      if(ioption.eq.3)then
            if(isfileopen.eq.0)Then
                call FSS_alert('ALERT!!','You must open an input file first')
                go to 100
                endif
      call contour
      endif
! 	-------------- 4 ------------------
!	if(ioption.eq.4)then
!	call MakeMyGrid2
!	go to 100
!	endif
! -------------------5---------------------
!      IF(IOPTION.EQ.5)CALL GroGrt

! -------------------55---------------------
!      IF(IOPTION.EQ.55)CALL GroGrt2

! ------------------ 6  -------------------
!      IF(IOPTION.EQ.6)call DiffGibbs
! ------------------ 66  -------------------
!      IF(IOPTION.EQ.66)call DiffNewton
! ------------------ 7 --------------------
      IF(IOPTION.EQ.7)then
            if(isfileopen.eq.0)Then
                call FSS_alert('ALERT!!','You must open an input file first')
                go to 100
                endif
	      CALL WHOLERock
	      endif
! ------------------ 8 --------------------
      IF(IOPTION.EQ.8)CALL GLOBAL(0)
! ------------------ 9 --------------------
      IF(IOPTION.EQ.9)CALL PLOTIN
! 	-------------- 10 ------------------
	if(ioption.eq.10)then
!	call digitplot
	go to 100
	endif
! 	-------------- 11 ------------------
	if(ioption.eq.11)then
	call thermodata
	go to 100
	endif
! 	-------------- 12 ------------------
 	if(ioption.eq.12)then
! 	call MakePseudoFromScratch
 	go to 100
 	endif
! 	-------------- 13 ------------------
	if(ioption.eq.13)then
!	call ForwardModel
	go to 100
	endif
! 	-------------- 14 ------------------
	if(ioption.eq.14)then
!	call Tweek2
	go to 100
	endif
! 	-------------- 4 ------------------
	if(ioption.eq.4)then
	call Tangent
	go to 100
	endif
! 	-------------- 44 ------------------
!	if(ioption.eq.44)then
!	call TangentPFluid
!	go to 100
!	endif
! 	-------------- 45 ------------------
!	if(ioption.eq.45)then
!	call TangentOSGrow
!	go to 100
!	endif
! 	-------------- 16 ------------------
	if(ioption.eq.16)then
!	call PlotKeq
	go to 100
	endif
! 	-------------- 17 ------------------
	if(ioption.eq.17)then
!	call ForwardNucDen
	go to 100
	endif
! 	-------------- 18 ------------------
	if(ioption.eq.18)then
!	call OverstepGarnetMain()
	go to 100
	endif
! 	-------------- 19 ------------------
	if(ioption.eq.19)then
!	call TangentGrid()
	go to 100
	endif

      go to 100
      END

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine GibbsBegin
    	use AWE_Interfaces
	implicit none
! c*******************************************
	include "Assemb.inc"
	include "Gibbsfiles.inc"
	include "Monit.inc"
	include "Newton.inc"
	include "Output.inc"
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
	include "SaveRxn.inc"
	include "Singular.inc"
!	include "MakeMyGrid2.inc"
! c*******************************************
      integer*4 ibegin,inew,i,istat

100	continue
        ibegin=0       					! default 
        WRITE(*,*)' ********************************'
	WRITE(*,*)' OPEN/SAVE menu'
        WRITE(*,*)' ********************************'
        WRITE(*,*)'   0 = Return'
        WRITE(*,*)'   1 = Select bulk composition from disk file BulkRockAnalyses'
        WRITE(*,*)'   2 = READ reference assemblage (MIF) from disk file'
	WRITE(*,*)'  -----------------------------'
        WRITE(*,*)'   5 = SAVE current problem as an input file'
	WRITE(*,*)'  -----------------------------'
        WRITE(*,*)'   6 = SAVE Graphics as an Adobe Illustrator (v3) file'
	WRITE(*,*)'  -----------------------------'
        WRITE(*,*)'   9 = SAVE Assemblage file and open a new one'
        WRITE(*,*)' CHOOSE OPTION'
      	read(*,*)ibegin

! ---------------------------
	if(ibegin.eq.0)return
! ---------------------------
       IF(IBEGIN.EQ.1)THEN				! Load a bulk composition
		call PickBulkComp
		go to 100
           	ENDIF
! ---------------------------
         IF(IBEGIN.EQ.2)THEN
            	INEW=0
            	numNew=0
            	i = 1			! this instructs subroutine to open a new file
            	CALL BEGIN(i,INEW)
            	IF(INEW.EQ.1)then       			! if true, we actually started a new problem
              		isfileopen = 1
			call PrintMinAssemblage
			call change				! to pick the assemblage
             		CALL REXN	
!			call PrintMinAssemblage
			call Printt(1)
              		go to 100
              		endif
         	endif
! ---------------------------
       IF(IBEGIN.EQ.5)THEN
           	if(isfileopen.eq.0)then
!       			write(*,*)'You must open a file before it can be saved.'
	                call fss_alert('ALERT!!','You must open an input file first')
             		go to 100
             		endif  
           	CALL SAVEIN(0,0)
           	GO TO 100
           	ENDIF
! ---------------------------
       IF(IBEGIN.EQ.6)THEN				! save postscript file
 !          	PSFileName = Trim(filein)//'.ai'        	! default ps file name
           	call psopcl(5)    				! save old PostScript scratch file
           	GO TO 100
           	ENDIF
! ---------------------------
	if(ibegin.eq.9)then
		close(16)
		open(unit=16,file="",status="NEW",iostat=istat)
			if(istat.ne.0)go to 100
		write(16,*)noSysCoInFile+1
      		WRITE(16,1010)(CoNameInFile(i),i=1,NoSysCoInFile)
1010  		format('Min                     ',40(4x,A4,2x))
         	go to 100
		endif
! ---------------------------
         go to 100
	end
	
	
	
