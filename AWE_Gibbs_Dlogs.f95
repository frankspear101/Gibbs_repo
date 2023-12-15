!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine PickMonitors(nvar,NeedToPick,Picked,ValuesPicked,iok)

!     Program to pick one or more monitors from the Gibbs variable list

!     NVAR =      Number of variables in list
!     VN1(50)     First part of variable name
!     VN2(50)     Second part of variable name
!     NeedToPick  Number of variables to pick
!     Picked(6)   variables that have been picked
 
!     iOk = 0     Normal return
!     iOk = 1     User hit cancel

!**************************************************************
    	use AWE_Interfaces
      	implicit none
!*****************************************
	include "Monit.inc"
!*****************************************
!      local variables
	integer*4 i,j,iok,nvar,NeedToPick,Picked(15),result
	Real*4 ValuesPicked(15)
	character*16 TheList(nvar)
	Character*16 editChars
!--------------------------------------------------------------
    	type(AWE_FormDialog)     :: PickMonDialog	        ! the dialog
    	type(AWE_FormLabel)      :: dialoglabel           ! static text
	type(AWE_FormLineEdit)   :: lineEditMon(15)       ! text edit box with default text
	type(AWE_FormLineEdit)   :: lineEditIterations         ! text box with default text
!    	type(AWE_FormComboBox)   :: comboBox(15)        ! drop down menu
    	type(AWE_FormComboBox),allocatable :: comboBox(:)        ! drop down menu



	allocate(comboBox(needToPick))
!       if(NeedToPick.gt.15)then
!       		call FSS_Alert('Gibbs Alert','Not enough dialog items for this many monitors.')
! 		iok = 1
!       		return
!       		endif

! 	write(*,*)nvar,NeedToPick,iok
! 	do i = 1,NeedToPick
! 		write(*,*)Picked(i),valuespicked(i)
! 		end do
! 	pause 'Input data -- hit return'

	PickMonDialog%title = "Pick monitor parameters"
    	call AWE_createDialog(PickMonDialog)

!     	Make the list of variables
      	do j = 1,nvar
		TheList(j) = '  '
	!      	WRITE(TheList(j),9015)j,VN1(j),VN2(j)
	!9015  	FORMAT(I3,1X,A2,A16)
		WRITE(TheList(j),9015)VN1(j),VN2(j)
	9015  	FORMAT(A2,A8)
	!      	TheList(j) = Trim(TheList(j))
		end do

	do i = 1,needToPick
		allocate(comboBox(i)%items(nvar))
		do j = 1,nvar
		!	comboBox(i)%items(j) = trim(TheList(j))
			comboBox(i)%items(j) = TheList(j)
			end do
	! 	this sets default monitors	
		comboBox(i)%selected = Picked(i)			! this should set the selection at the last one chosen
		call AWE_addToDialog(comboBox(i), PickMonDialog)
	!	Add the edit field for the delta
		lineEditMon(i)%title = "Delta value ="
		if(Vn1(Picked(i))(1:2).eq."M_")then
			write(editChars,'(F12.4)')ValuesPicked(i)*1000  ! write millimoles for phase amounts
			else
			write(editChars,'(F12.4)')ValuesPicked(i)
			endif
		lineEditMon(i)%text = editChars
		call AWE_addToDialog(lineEditMon(i), PickMonDialog)
		end do


!     write default value for Nstep
	lineEditIterations%title = "Iterations"
	write(editChars,'(I4)')SNstep
	lineEditIterations%text = editChars
	call AWE_addToDialog(lineEditIterations, PickMonDialog)
      

	result = AWE_showDialog(PickMonDialog)
! 	write(*,*)'result = ',result
! 	pause 'hit return to continue'
    	if (result == 1) then
		iok = 0
		do i = 1,NeedToPick
			read(lineEditMon(i)%text,*)ValuesPicked(i)
			Picked(i) = comboBox(i)%selected
			if(Vn1(Picked(i))(1:2).eq."M_")then
				ValuesPicked(i)=ValuesPicked(i)/1000.0
				endif
			end do
 		read(lineEditIterations%text,*)snStep
		else
		iok = 1
   		endif

! 	write(*,*)nvar,NeedToPick,iok
! 	do i = 1,NeedToPick
! 		write(*,*)Picked(i),valuespicked(i)
! 		end do
! 	pause 'output data -- hit return'

! 	write(*,*)'about to return'
! 	pause 'hit return to continue'
      return           
      end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine PickAbsContour(nvar,NeedToPick,Picked,ValuesPicked,myColor,iok)
!      Subroutine PickAbsContour(nvar,NeedToPick,Picked,ValuesPicked,zinc4,iok)

!     Program to pick one or more monitors from the Gibbs variable list

!     NVAR =      Number of variables in list
!     VN1(50)     First part of variable name
!     VN2(50)     Second part of variable name
!     NeedToPick  Number of variables to pick
!     Picked(6)   variables that have been picked
 
!     iOk = 0     Normal return
!     iOk = 1     User hit cancel

!**************************************************************
    	use AWE_Interfaces
      	implicit none
!*****************************************
	include "Monit.inc"
	INCLUDE "PlotStuff.inc"	
!*****************************************
!      local variables
	integer*4 i,j,iok,nvar,NeedToPick,Picked(10),result,myColor
	Real*4 ValuesPicked(10)
	character*16 TheList(nvar)
	Character*16 editChars
!--------------------------------------------------------------
    	type(AWE_FormDialog)     :: PickConDialog	        ! the dialog
    	type(AWE_FormLabel)      :: dialoglabel           ! static text
	type(AWE_FormLineEdit)   :: lineEditMon(10)       ! text edit box with default text
	type(AWE_FormLineEdit)   :: lineEditIterations         ! text box with default text
!    	type(AWE_FormComboBox)   :: comboBox(10)        ! drop down menu
    	type(AWE_FormComboBox),allocatable :: comboBox(:)        ! drop down menu
    	type(AWE_FormComboBox)   :: ColorComboBox        ! drop down menu


 	allocate(comboBox(needToPick))
!      if(NeedToPick.gt.10)then
!       		call FSS_Alert('Gibbs Alert','Not enough dialog items for this many monitors.')
! 		iok = 1
!       		return
!       		endif

	PickConDialog%title = "Absolute contour routine"
    	call AWE_createDialog(PickConDialog)
!	add static text
	dialoglabel%text = "Pick variable to contour:"
	call AWE_addToDialog(dialoglabel,PickConDialog)

!     	Make the list of variables
      	do j = 1,nvar
		TheList(j) = '  '
	!      	WRITE(TheList(j),9015)j,VN1(j),VN2(j)
	!9015  	FORMAT(I3,1X,A2,A16)
		WRITE(TheList(j),9015)VN1(j),VN2(j)
	9015  	FORMAT(A2,A8)
	!      	TheList(j) = Trim(TheList(j))
		end do

	do i = 1,needToPick
		allocate(comboBox(i)%items(nvar))
		comboBox%title = "Contour variable:"
		do j = 1,nvar
		!	comboBox(i)%items(j) = trim(TheList(j))
			comboBox(i)%items(j) = TheList(j)
			end do

	! 	this sets default monitors	
		comboBox(i)%selected = Picked(i)			! this should set the selection at the last one chosen
		call AWE_addToDialog(comboBox(i), PickConDialog)
	!	Add the edit field for the contour value
		lineEditMon(i)%title = "Value to contour ="
		if(Vn1(Picked(i))(1:2).eq."M_")then
			write(editChars,'(F12.4)')ValuesPicked(i)*1000  ! write millimoles for phase amounts
			else
			write(editChars,'(F12.4)')ValuesPicked(i)
			endif
		lineEditMon(i)%text = editChars
		call AWE_addToDialog(lineEditMon(i), PickConDialog)
		end do

! add comboBox drop down menu
	ColorComboBox%title = "Select your color"
    	allocate(ColorComboBox%items(AWEnumColors))
	ColorComboBox%selected = myColor			! this should set the selection at the last one chosen
	do i = 1,AWEnumColors
		ColorComboBox%items(i) = trim(AWEColorName(i))
		end do
    	call AWE_addToDialog(ColorComboBox, PickConDialog)

	result = AWE_showDialog(PickConDialog)
    	if (result == 1) then
		iok = 0
		do i = 1,NeedToPick
 		read(lineEditMon(i)%text,*)ValuesPicked(i)
		Picked(i) = comboBox(i)%selected
		if(Vn1(Picked(i))(1:2).eq."M_")then
			ValuesPicked(i) = ValuesPicked(i)/1000.0d0	! rescale moles from millimoles to moles
			endif
		end do
		myColor = ColorComboBox%selected
! 		read(lineEditIterations%text,*)snStep
		else
		iok = 1
   		endif

      return           
      end


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine PickContour(nvar,NeedToPick,Picked,myColor,iok)
!      Subroutine PickAbsContour(nvar,NeedToPick,Picked,ValuesPicked,zinc4,iok)

!     Program to pick one or more monitors from the Gibbs variable list

!     NVAR =      Number of variables in list
!     VN1(50)     First part of variable name
!     VN2(50)     Second part of variable name
!     NeedToPick  Number of variables to pick
!     Picked(6)   variables that have been picked
 
!     iOk = 0     Normal return
!     iOk = 1     User hit cancel

!**************************************************************
    	use AWE_Interfaces
      	implicit none
!*****************************************
	include "Monit.inc"
	INCLUDE "PlotStuff.inc"	
!*****************************************
!      local variables
	integer*4 i,j,iok,nvar,NeedToPick,Picked(10),result,myColor
!	Real*4 ValuesPicked(10)
	character*16 TheList(nvar)
!--------------------------------------------------------------
    	type(AWE_FormDialog)     :: PickConDialog	        ! the dialog
    	type(AWE_FormLabel)      :: dialoglabel           ! static text
	type(AWE_FormLineEdit)   :: lineEditMon(10)       ! text edit box with default text
	type(AWE_FormLineEdit)   :: lineEditIterations         ! text box with default text
    	type(AWE_FormComboBox),allocatable :: comboBox(:)        ! drop down menu
!    	type(AWE_FormComboBox)   :: comboBox(10)        ! drop down menu
    	type(AWE_FormComboBox)   :: ColorComboBox        ! drop down menu


	allocate(comboBox(needToPick))
!       if(NeedToPick.gt.10)then
!       		call FSS_Alert('Gibbs Alert','Not enough dialog items for this many monitors.')
! 		iok = 1
!       		return
!       		endif

	PickConDialog%title = "Contour routine"
    	call AWE_createDialog(PickConDialog)
!	add static text
	dialoglabel%text = "Pick variable to hold constant:"
	call AWE_addToDialog(dialoglabel,PickConDialog)

!     	Make the list of variables
      	do j = 1,nvar
		TheList(j) = '  '
	!      	WRITE(TheList(j),9015)j,VN1(j),VN2(j)
	!9015  	FORMAT(I3,1X,A2,A16)
		WRITE(TheList(j),9015)VN1(j),VN2(j)
	9015  	FORMAT(A2,A8)
	!      	TheList(j) = Trim(TheList(j))
		end do

	do i = 1,needToPick
		allocate(comboBox(i)%items(nvar))
		comboBox%title = "Contour variable:"
		do j = 1,nvar
			comboBox(i)%items(j) = TheList(j)
			end do

	! 	this sets default monitors	
		comboBox(i)%selected = Picked(i)			! this should set the selection at the last one chosen
		call AWE_addToDialog(comboBox(i), PickConDialog)
		end do

! add comboBox drop down menu
	ColorComboBox%title = "Select your color"
    	allocate(ColorComboBox%items(AWEnumColors))
	ColorComboBox%selected = myColor			! this should set the selection at the last one chosen
	do i = 1,AWEnumColors
		ColorComboBox%items(i) = trim(AWEColorName(i))
		end do
    	call AWE_addToDialog(ColorComboBox, PickConDialog)


	result = AWE_showDialog(PickConDialog)
    	if (result == 1) then
		iok = 0
		do i = 1,NeedToPick
	! 		read(lineEditMon(i)%text,*)ValuesPicked(i)
			Picked(i) = comboBox(i)%selected
			end do
		myColor = ColorComboBox%selected
! 		read(lineEditIterations%text,*)snStep
		else
		iok = 1
   		endif

      return           
      end


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine FSS_alert(title, text)
    	use AWE_Interfaces
!	interface Alert
!	subroutine AWE_alertBox(title, text)
	character(len=*) :: title, text
	call AWE_alertBox(title,text)
	end
!	end subroutine AWE_alertBox
!	end interface
!	title is used as the title of the alert box text is the text that will be displayed in it.
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine GibbsorNewton
    	use AWE_Interfaces
	implicit none

!     Program to set Gibbs/Newton switch
!****************************************
	include "Newton.inc"
!****************************************
!      	local variables
      	integer*4 result
    	character *32 editChars
!********/Newton/************************
!      real*8 tolnewton,newtonstep
!      integer*4 inewton
!      COMMON /Newton/tolnewton,newtonstep,inewton
!************************************

    	type(AWE_FormDialog)     :: GibbsOrNewtonDialog	        ! the dialog
    	type(AWE_FormLabel)      :: dialoglabel           ! static text
	type(AWE_FormLineEdit)   :: lineEditNewtTol       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditNewtStep         ! text box with default text
	TYPE(AWE_FormRadioButtons) :: NewtButtons
	! NewtButtons%title = "Model"
	allocate (NewtButtons%items(2))
	NewtButtons%items(1) = "Gibbs Method (differential thermodynamics)"
	NewtButtons%items(2) = "Newton Method (integral thermodynamics)"

! set the dialog title and create it
    	GibbsOrNewtonDialog%title = "Gibbs or Newton dialog"
    	call AWE_createDialog(GibbsOrNewtonDialog)


!	Set Gibbs/Newton switch
	if(inewton.eq.0)then	
       	        NewtButtons%selected = 1    !using Gibbs method
	        else        
       	        NewtButtons%selected = 2    !using Newton's method
	        endif
	CALL AWE_AddToDialog(NewtButtons, GibbsOrNewtonDialog)

!	Set existing tolnewton
	lineEditNewtTol%title = "Newton Tolerance (.1E-3 is default)"
	write(editChars,10)tolNewton
10	format(e8.3)
	lineEditNewtTol%text = editChars
	call AWE_addToDialog(lineEditNewtTol, GibbsOrNewtonDialog)

	lineEditNewtStep%title = "Newton Step Size (0.4 is default)"
	write(editChars,11)newtonStep
11	format(F8.3)
	lineEditNewtStep%text = editChars
	call AWE_addToDialog(lineEditNewtStep, GibbsOrNewtonDialog)

! show the dialog. 
!         If result==0, the Cancel button was clicked
!         If result==1, the OK button was clicked
	result = AWE_showDialog(GibbsOrNewtonDialog)
    	if (result == 1) then
!		iok = 0
 		read(lineEditNewtTol%text,*)tolNewton
 		read(lineEditNewtStep%text,*)newtonStep
		else
!		iok = 1
   		endif
	if(NewtButtons%selected == 1)then
		iNewton = 0		!using Gibbs method
		else
		iNewton = 1 
		endif  
 		
	return
	end

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine PickPT_grid(TCstart,PBstart,TCend,PBend,Tinc,Pinc,iok)
! 	Routine to set T and P for AFM diagrams in Gibbs
    	use AWE_Interfaces
	implicit none
! --------------------------------------------------------------------
!     local variables
      integer*4 iok,result
    	character *16 editChars
	real*4 TCstart,PBstart,TCend,PBend,Pinc,Tinc
    	type(AWE_FormDialog)     :: GridDialog	        ! the dialog
    	type(AWE_FormLabel)      :: dialoglabel           ! static text
	type(AWE_FormLineEdit)   :: lineEditTCstart       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditTCend         ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditPBstart       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditPBend         ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditTinc          ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditPinc          ! text box with default text


! set the dialog title and create it
    GridDialog%title = "PT grid dialog"
    call AWE_createDialog(GridDialog)

! add static text

	dialoglabel%text = "Set P and T endpoints and increments"
	call AWE_addToDialog(dialoglabel, GridDialog)

	lineEditTCstart%title = "TC start"
	write(editChars,*)TCstart
	lineEditTCstart%text = editChars
	call AWE_addToDialog(lineEditTCstart, GridDialog)

	lineEditTCend%title = "TC end"
	write(editChars,*)TCend
	lineEditTCend%text = editChars
	call AWE_addToDialog(lineEditTCend, GridDialog)

	lineEditTinc%title = "T  inc"
	write(editChars,*)Tinc
	lineEditTinc%text = editChars
	call AWE_addToDialog(lineEditTinc, GridDialog)

	lineEditPBstart%title = "PB start"
	write(editChars,*)PBstart
	lineEditPBstart%text = editChars
	call AWE_addToDialog(lineEditPBstart, GridDialog)

	lineEditPBend%title = "PB end"
	write(editChars,*)PBend
	lineEditPBend%text = editChars
	call AWE_addToDialog(lineEditPBend, GridDialog)

	lineEditPinc%title = "P  inc"
	write(editChars,*)Pinc
	lineEditPinc%text = editChars
	call AWE_addToDialog(lineEditPinc, GridDialog)

! show the dialog. 
!         If result==0, the Cancel button was clicked
!         If result==1, the OK button was clicked
	result = AWE_showDialog(GridDialog)
    	if (result == 1) then
		iok = 0
 		read(lineEditTCstart%text,*)TCstart
 		read(lineEditTCend%text,*)TCend
 		read(lineEditTinc%text,*)Tinc
 		read(lineEditPBstart%text,*)PBstart
 		read(lineEditPBend%text,*)PBend
 		read(lineEditPinc%text,*)Pinc
		else
		iok = 1
   		endif

      return           
      end
	

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine AWE_SetOutput
    	use AWE_Interfaces

!     Program to set Gibbs output options
!**************************************************************
      implicit none

! ***************************************************************
! *****************************************
	include "Output.inc"		
! !****************************************
!       local variables
      integer*4 i,result

! --------------------------------------------------------------
! 	Variable INTEGER iLong(10) in commonblock /OUTPUT/ contains switches for output as follows
! 	0 = off; 1 = on
! 	1 = Error output only
! 	2 = Jacobian Matrix
! 	3 = Thermodyanmic data
! 	4 = Running sum
! 	5 = Transformation matrix
! 	6 = Activity calculations
! 	7 = REXN independent reactions

! dialog and item types
    	type(AWE_FormDialog)     :: OutputDialog	        ! the dialog
    	type(AWE_FormLabel)      :: label           ! static text
!    	type(AWE_FormComboBox)   :: comboBox        ! drop down menu
    	type(AWE_FormCheckBox)   :: CheckBox(9)        ! Check box

! set the dialog title and create it

	OutPutDialog%title = "Select program output options"
    	call AWE_createDialog(OutPutDialog)

! add static text

	label%text = "Select the output options"
	call AWE_addToDialog(label, OutputDialog)

	CheckBox(1)%title = 'Error output only'
	CheckBox(1)%checked = .false.
	CheckBox(2)%title = 'Jacobian matrix'
	CheckBox(2)%checked = .false.
	CheckBox(3)%title = 'Thermodynamic data for every mineral'
	CheckBox(3)%checked = .false.
	CheckBox(4)%title = 'Running sum of each step'
	CheckBox(4)%checked = .false.
	CheckBox(5)%title = 'Transformation matrix'
	CheckBox(5)%checked = .false.
	CheckBox(6)%title = 'Activity calculations'
	CheckBox(6)%checked = .false.
	CheckBox(7)%title = 'REXN: independent reaction output on each step'
	CheckBox(7)%checked = .false.
	CheckBox(8)%title = 'Print phase compositions in Wt% on (7)'
	CheckBox(8)%checked = .false.
	CheckBox(9)%title = 'Write minerals to .asm file on (7)'
	CheckBox(9)%checked = .false.
	
!	This is where the existing output flags are set
	do 10 i = 1,9
	if(iLong(i).eq.1)then
		CheckBox(i)%checked = .true.
		endif
	call AWE_addToDialog(CheckBox(i),OutputDialog)
10	continue



! show the dialog. 
!         If result==0, the Cancel button was clicked
!         If result==1, the OK button was clicked
	result = AWE_showDialog(OutputDialog)
    	if (result == 1) then
		do 20 i = 1,9
		if(CheckBox(i)%checked)then
			iLong(i) = 1
			else
			iLong(i) = 0
			endif
20		continue
		endif


      return           
      end





!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine AWE_Pick_one(TheOne,ListLength,TheList,WhatToPick,iok,iselect)
! purpose:    to performance the pick a item from a list via. list utility
! ------------------------------------------------------
! name        type   I/O   description
! ------------------------------------------------------
! theOne:     int*4   O  the item number picked
! ListLength: int*4   I  the length of the list
! The List:   char*80 I  the List to be displayed and to be picked one item from
! WhattoPick:            prompt string
! IOK         int*4   O  flag to indicate the pick status
! iselect:    int*4  I/O default initial position 0 = use last selection; 1 = use selection #1
! ------------------------------------------------------

    	use AWE_Interfaces
    	implicit none
	character*80 TheList(*)
	character*80 WhatToPick
	Integer*4 TheOne,ListLength,iok,iselect,i,result

! 	dialog and item types
    	type(AWE_FormDialog)     :: PickList	        ! the dialog
    	type(AWE_FormLabel)      :: label           ! static text
    	type(AWE_FormComboBox)   :: comboBox        ! drop down menu
    	type(AWE_FormCheckBox)   :: CheckBox        ! Check box

! 	set the dialog title and create it
	PickList%title = "Pick one"
	call AWE_createDialog(PickList)

! 	add static text
	label%text = WhatToPick
	call AWE_addToDialog(label, PickList)

! 	add comboBox drop down menu

!    	comboBox%title = "Select your poison"
	allocate(comboBox%items(ListLength))

	if(iselect.eq.0)then
		comboBox%selected = TheOne			! this should set the selection at the last one chosen
		else
		comboBox%selected = 1
		endif
	do 10 i = 1,ListLength
	comboBox%items(i) = trim(TheList(i))
10	continue	

	call AWE_addToDialog(comboBox, PickList)

	! Add a checkbox (should fix the Cancel vs OK default)
!	CheckBox%title = 'Check box for a surprise!'
!	CheckBox%checked = .false.
!	call AWE_addToDialog(CheckBox,PickList)

! 	show the dialog. 
!         If result==0, the Cancel button was clicked
!         If result==1, the OK button was clicked
	result = AWE_showDialog(PickList)
    	if (result == 1) then
      		print *,"You picked...            ", comboBox%items(comboBox%selected)
		TheOne = comboBox%selected
		iok = 0
! 		if(CheckBox%checked)then
! 			call FSS_Alert('Good Job','Contratulations! You have won a trip to the Funny Farm') 
! 			endif
		else
		iok = 1
		endif

	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine SetContourLimit()
    	use AWE_Interfaces
	implicit none
! --------------------------------------------------------------------
	include "PlotStuff.inc"
	include "PlotGibbs.inc"
! --------------------------------------------------------------------
!     local variables
      integer*4 iok,result
    	character *16 editChars
    	type(AWE_FormDialog)     :: ContourLimitsDialog	        ! the dialog
    	type(AWE_FormLabel)      :: dialoglabel           ! static text
	type(AWE_FormLineEdit)   :: lineEditXMIN       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditXMAX         ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditYMIN       ! text box with default text
	type(AWE_FormLineEdit)   :: lineEditYMAX         ! text box with default text


! set the dialog title and create it
    ContourLimitsDialog%title = "Contour Limits"
    call AWE_createDialog(ContourLimitsDialog)

! add static text

	dialoglabel%text = "Set Contour Limits"
	call AWE_addToDialog(dialoglabel, ContourLimitsDialog)

	lineEditXMIN%title = "X minimum"
	write(editChars,*)xMinBound
	lineEditXMIN%text = editChars
	call AWE_addToDialog(lineEditXMIN, ContourLimitsDialog)

	lineEditXMAX%title = "X maximum"
	write(editChars,*)xMaxBound
	lineEditXMAX%text = editChars
	call AWE_addToDialog(lineEditXMAX, ContourLimitsDialog)

	lineEditYMIN%title = "Y minimum"
	write(editChars,*)YMinBound
	lineEditYMIN%text = editChars
	call AWE_addToDialog(lineEditYMIN, ContourLimitsDialog)

	lineEditYMAX%title = "Y maximum"
	write(editChars,*)YMaxBound
	lineEditYMAX%text = editChars
	call AWE_addToDialog(lineEditYMAX, ContourLimitsDialog)


! show the dialog. 
!         If result==0, the Cancel button was clicked
!         If result==1, the OK button was clicked
	result = AWE_showDialog(ContourLimitsDialog)
    	if (result == 1) then
		iok = 0
 		read(lineEditXMIN%text,*)XminBound
 		read(lineEditXMAX%text,*)XmaxBound
 		read(lineEditYMIN%text,*)YminBound
 		read(lineEditYMAX%text,*)YmaxBound
		else
		iok = 1
   		endif

      return           
      end
	
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine PkXYAxe(nvar,NeedToPick,Pick2,iok)

!c     Program to pick one or more things from a list

!c	itype = 1 pick monitor parameters
!c	itype = 2 pick CompPlot variables
!c	itype = 3 pick minerals from list

!c
!c     NeedToPick  Number of variables to pick
!c     Picked(10)   variables that have been picked
!c
!c     iOk = 0     Normal return
!c     iOk = 1     User hit cancel

!c**************************************************************
    	use AWE_Interfaces
	implicit none

!c***************************************************************
!c*****************************************
!	include "Assemb.inc"
	include "Monit.inc"
!	include "PlotStuff.inc"
!c*****************************************
!c      local variables
      integer*4 i,j,iok,NeedToPick,Pick2(5,2),ibox,result,nvar
	character*16 TheList(nvar)
!--------------------------------------------------------------
    	type(AWE_FormDialog)     :: PickAxeDialog	        ! the dialog
    	type(AWE_FormLabel)      :: dialoglabel           ! static text
	type(AWE_FormLineEdit)   :: lineEditMon(10)       ! text edit box with default text
	type(AWE_FormLineEdit)   :: lineEditIterations         ! text box with default text
!    	type(AWE_FormComboBox)   :: comboBox(10)        ! drop down menu
    	type(AWE_FormComboBox),allocatable :: comboBox(:)        ! drop down menu
!c--------------------------------------------------------------


	allocate(comboBox(needToPick))
!       if(NeedToPick.gt.10)then
!       		call FSS_Alert('ALERT!!','Not enough dialog items for this many monitors.')
! 		iok = 1
! 		return
! 		endif

	allocate(comboBox(needToPick+1))

	PickAxeDialog%title = "Pick stuff to plot"
    	call AWE_createDialog(PickAxeDialog)



!c     Make the list
!c	We make a list of ALLX variables for X-Y plots
!c	and a list of Calculated component names for ternary and tetrahedral plots
!c	and a list of minerals for other places

!	X-Y plots - list ALLX variables
	do j = 1,nvar
		TheList(j) = '  '
		WRITE(TheList(j),9015)j,VN1(j),VN2(j)
	9015	FORMAT(I3,1X,A2,A16)
		TheList(j) = Trim(TheList(j))
		end do

	ibox = 0
	do i = 1,needToPick
		ibox = ibox + 1
		allocate(comboBox(ibox)%items(nvar))
		allocate(comboBox(ibox+1)%items(nvar))
		do j = 1,nvar
			comboBox(ibox)%items(j) = TheList(j)
			comboBox(ibox+1)%items(j) = TheList(j)
			end do
		comboBox(ibox)%title = 'X axis to plot'			! this should set the selection at the last one chosen
	! 	this sets default monitors	
		comboBox(ibox)%selected = Pick2(i,1)			! this should set the selection at the last one chosen
		call AWE_addToDialog(comboBox(ibox), PickAxeDialog)

		ibox = ibox + 1
		comboBox(ibox)%title = 'Y axis to plot'			! this should set the selection at the last one chosen
	! 	this sets default monitors	
		comboBox(ibox)%selected = Pick2(i,2)			! this should set the selection at the last one chosen
		call AWE_addToDialog(comboBox(ibox), PickAxeDialog)
		dialoglabel%text = '--------------------'
		call AWE_AddToDialog(dialoglabel,PickAxeDialog)
		end do



	result = AWE_showDialog(PickAxeDialog)
    	if (result == 1) then
		iok = 0
		ibox = 0
		do i = 1,NeedToPick
			ibox = ibox + 1
			Pick2(i,1) = comboBox(ibox)%selected
			ibox = ibox + 1
			Pick2(i,2) = comboBox(ibox)%selected
			end do

		else
		iok = 1
   		endif

      return           
      end





