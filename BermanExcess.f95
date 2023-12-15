	Subroutine BermanExcess(n,numWTerms,WG,WIndex,X,uexcess)
!	call BermanExcess(numMargulesSiteCats(k,jj),numMargulesWterms(k,jj),MargulesWG,WIndex,Xmarg,uexcessBerman(jj))
!	routine to calculate µexcess using Berman and Brown (1984) and Berman (1990) generalized Margules expansion
!	See also equation 7-152 of Spear (1993)
!	This routine is set up to calculate uexcess for each distinct site
!	the contribution to each phase component will then be calculate in the phase component loop in the calling routine (Subroutine Activity)

	implicit none
	integer*4 n,m,i,j,k,L,Q,index,kk,jj
	integer*4 numWTerms
!	integer*4 WIndex(*,*,*)
!	real*8 WG(*,*,*),X(*)
	integer*4 WIndex(20,20,20)	! dimension = numMargulesWtermsMax
	real*8 WG(20,20,20)		! dimension = numMargulesWtermsMax
	real*8 X(9)			! dimension = numMargulesSiteCatsMax
	real*8 uex,uexcess(9)		! dimension = numMargulesSiteCatsMax


!	write(*,*)'i   j   k   Q'
	do 5 m = 1,n		! loop for the µexcess for each individual atom on the site
!	write(*,*)
!	write(*,*)' Doing component m = ',m
!	write(*,*)'i   j   k   Q'
	uex = 0.0d0
	do 10 i = 1,n-1
	do 15 j = i,n
	do 20 k = j,n
	if(k.eq.i)then
		goto 20		! k is not equal to i
		endif
	Q = 0
	if(m.eq.i)Q=Q+1
	if(m.eq.j)Q=Q+1
	if(m.eq.k)Q=Q+1
	
!	write(*,*)i,j,k,Q
!	pause 'hit return 1'
	index = i*100 + 10*j + k		! e.g. 112,122 etc.
!	do 25 L = 1,numWterms(kk,jj) 		! loop through every W term to find the right one
!	if(index.eq.Windex(kk,jj,L))then	
	do 25 L = 1,numWterms 			! loop through every W term to find the right one
	if(index.eq.Windex(i,j,k))then	
!		sum = sum + WG(i,j,k)*(Q*Xi*Xj*Xk/Xm - 2*Xi*Xj*Xk)
		uex = uex + WG(i,j,k)*(float(Q)*X(i)*X(j)*X(k)/X(m) - 2*X(i)*X(j)*X(k))
		go to 20
		endif
25	continue
20	continue
15	continue
10	continue
!	pause 'hit return 2'
	uexcess(m) = uex
5	continue

	return
	end
