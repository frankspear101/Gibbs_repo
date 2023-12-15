! **********/SOLUTE/************************************************
!     COMMON      BLOCK SOLUTE
	integer*4 maxSolute
	parameter (maxSolute = 500)
      REAL*8 A(maxSolute,maxSolute),XX(maxSolute,maxSolute),AA(maxSolute,maxSolute),YY(maxSolute),Asav(maxSolute,maxSolute)
      COMMON /SOLUTE/A,AA,XX,YY,Asav
! *******************************************************************
