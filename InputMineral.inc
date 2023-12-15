!********************************************************************
!	common block to store the InputMineral information from a master input Mineral file
!c	the InputMineralFile is defined in include file "files.inc" as character*64
	integer*4 inputMineral(20),noMinInArray,noInputMineral
	character*255 inputMineralArray(200)	!internal array to store InputMineralFile
	common/InputMineralCom/inputMineral,inputMineralArray,noMinInArray,noInputMineral
!C********************************************************************
