!	100 = number of total phase components
!	100 = number of total system compontnts
!	integer maxPCCoTr; parameter(maxPCCoTr=150)		! maximum total phase components in a problem 
	integer maxPCCoTr; parameter(maxPCCoTr=100)		! maximum total phase components in a problem 
	! maxPCCoTr must be the same value as maxPC in Assemb.inc and maxPCMonit -- if you change it here, you must change it there 
      REAL*8 A(maxPCCoTr,maxPCCoTr),B(maxPCCoTr,maxPCCoTr)
      integer*4 pointer(maxPCCoTr)			! keeps track of the order of phase components when rows are swapped
	common /CoTransCommon/a,b,pointer
