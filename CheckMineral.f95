! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine CheckStoichiometry(K,Xtemp,iflag)
!     Subroutine to check that the stochiometry of a mineral is not violated
!	There are several ways this must be checked, depending on the mineral
!	For a simple, 1-site phase (PhaseType = 1
!		0<X<1
!	For a site fraction phse (PhaseType = 3)
!		for each site fraction 0<X<1
!
!	This routine is called from block 5000 in subroutine Compute
!		and in block 5000 in routine ParallelToTangent
!
!	k = phase of interest (from current assemblage arrays)
!	xTemp = mole fractions of components (not necessarily sites, but components - some can be negative (e.g. reciprocal solutions)
!	returns:
!		iflag = 0	mineral is OK
!		iflag = 1	mineral is not OK
      implicit none
! ****************************************
	include "Assemb.inc"
! ****************************************
      	integer*4 K,iflag,j
      	real*8 Xtemp(20)
	iflag = 0
	call CheckSiteFractions(k,xTemp,iflag)
	return
!	if here, then we didn't recognize the mineral. Put out a warning
	write(*,*)'+++++++++++++++++++++++++++++++++'
	write(*,*)'Mineral not recognized in Subroutine CheckStoichiometry (file CheckMineral.f)'
	write(*,*)' Not sure what to do.....code needs fixing'
	write(*,*)'Phase = ',k,phName(k),minRec(k)
	write(*,*)'PhaseType(k) = ',PhaseType(k)
	do 50 j = 1,numPhCo(K)
	write(*,*)phCoName(k,j),xPhCo(k,j),xTemp(j)
50	continue
	iflag = 1
	pause 'Hit return to continue...Line 38 of Sub CheckStoichiometry'
	return
	end	
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine CheckSiteFractions(K,Xtemp,iflag)
!     Subroutine to check that the site fractions of minerals don't go negative
!	called from block 5000 in subroutine Compute
!	This is only for minerals special case = 9
      implicit none
! ****************************************
	include "Assemb.inc"
! ****************************************
      	integer*4 K,iflag,j,L
      	real*8 Xtemp(20)
	real*8 SiteFraction(20)

!      	Calculate mole fractions of cations on sites
!     	in terms of mole fractions of thermodynamic components
! 	This is the same code as in subroutine XtoSite except here I am
! 	using Xtemp for compositions instead of X because I need to calculate
! 	the finite difference derivative w/r/t X and Xtemp is passed through from calling routine
!      	LL = SiteAtomPointer(K)-1
      	do 611 L = 1,numSiteAtom(K)
      	SiteFraction(L) = 0.0d0
      	do 610 j = 1,numPhCo(K)
      	SiteFraction(L) = SiteFraction(L) + XToSiteAtom(k,L,j)*Xtemp(j)/SiteMultiplicity(k,L)
610   continue
611   continue
!	for testing only
!		do 9910 L = 1,numSiteAtom(K)
!		write(*,*)SiteAtomName(LL+L),SiteFraction(L)
!9910		continue
!	if(minrec(k).eq.41)then			! this is chloritoid
!		write(12,*)k,phName(k)
!		write(12,706)(SiteAtomName(k,L),L=1,numSiteAtom(k))
!		write(12,703)(SiteFraction(L),L=1,numSiteAtom(k))
!706   		format(4x,20A12)
!703   		format (20F12.6)
!		endif

! 	Check to see that no atoms in use are less than zero
	do 10 L = 1,numSiteAtom(K)
        SiteFraction(L) = SiteFraction(L)+1.d-60	! if zero, make very slightly larger
! 	if(SiteFraction(L).le.0.0D0)then
!		Try this code to see if we can avoid bombing because of 0 site fraction
!		Jan 22, 2019
!		Code did not work because newtonStep didn't decrease and component became negative.
!		SiteFraction(L) = 1.0d-20
!		endif
 	if(SiteFraction(L).le.0.0D0)then
!	for testing only
!		write(12,*)k,phName(k)
!		write(12,706)(SiteAtomName(k,LLL),LLL=1,numSiteAtom(k))
!		write(12,703)(SiteFraction(LLL),LLL=1,numSiteAtom(k))
!706   		format(4x,20A12)
!703   		format (20F12.6)
 		iflag = 1
 		return
 		endif
10	continue
	iflag = 0
	return
	end

