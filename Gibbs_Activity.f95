! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      Subroutine Activity(K,uMin,Xtemp,TKtemp,PBtemp,R,iflag)
!     Activity of a generic multi site mineral, as defined in data file
! 	Subroutine returns R*T*ln(act) of components in uMin
      implicit none
! ****************************************
	include "Assemb.inc"
	include "Output.inc"
! ****************************************
      	integer*4 K,iflag,j,L,n,i,i1,i2,j1,k1,jj,ii
      	real*8 uMin(PhCoMax),Xtemp(PhCoMax),TKtemp,Pbtemp,R,RTK,lnatemp,uexcess,sum,pkzero1,pkzero2,TR
	real*8 SiteFraction(numSiteAtomMax),PhiAFS(PhCoMax),tempsum,ASFtotal(PhCoMax),WGHP11(50),WGHP11temp(50)
	real*8 WGrecip,uRecip,uRecipTemp

	real*8 MargulesWG(numMargulesWTermsMax,numMargulesWTermsMax,numMargulesWTermsMax),XMarg(numMargulesSiteCatsMax)
	real*8 uexcessBerman(numSiteAtomMax),uexcessTemp(numMargulesSiteCatsMax),uexcessMargules
	integer*4 WIndex(numMargulesWTermsMax,numMargulesWTermsMax,numMargulesWTermsMax),index,index1
	real*8 Xfac,XFe,XpOl,XFeSp,XMgSp,XAlSp,XTiSp,XFe3Sp,AlTiFe3TotSp,xsp,ysp,zsp
 	real*8 pValue(PhCoMax),lnaIdeal(PhCoMax),RTlnGamma(PhCoMax)
 	common /gamma/pValue,lnaIdeal,RTlnGamma
 
      	data TR/298.15D0/

	if(iLong(6).eq.1)then
		write(12,*)
		write(12,*)'TK = ',TKtemp,'    Pbar = ',Pbtemp
		write(12,*)'Min ',Minrec(K),phname(K)
		endif

	RTK = R*TKtemp


!-------------------------
!	This section sets the appropriate composition terms into the array Xtem2 for use in calculating the excess energy
	select case (DatasetKey(k))
	case (1)		! Spac
!		Ideal activity calculations
	!      	Calculate mole fractions of cations on sites
	!     	in terms of mole fractions of thermodynamic components
	! 	This is the same code as in subroutine XtoSite except here I am
	! 	using Xtemp for compositions instead of X because I need to calculate
	! 	the finite difference derivative w/r/t X and Xtemp is passed through from calling routine
		do L = 1,numSiteAtom(K)
			SiteFraction(L) = 0
				do j = 1,numPhCo(K)
				SiteFraction(L) = SiteFraction(L) + XToSiteAtom(k,L,j)*Xtemp(j)/SiteMultiplicity(k,L)
				end do
			end do
! 		Check to see that no atoms in use are less than zero
		do L = 1,numSiteAtom(K)
			SiteFraction(L) = SiteFraction(L)+1.d-20	! if zero, make very slightly larger
			if(SiteFraction(L).le.1.0D-30)go to 999
			end do

		if(iLong(6).eq.1)then
			write(12,*)' Component X values'
			write(12,1016)(Xtemp(j),j=1,numPhCo(K))
1016			format(15E15.8)
			write(12,*)' Site Fractions '
			write(12,*)(SiteFraction(L),L=1,numSiteAtom(K))
			endif
! 		If here, then all site atoms are > 0.0

!		Calculate Margules excess values for every site
		if(numMargulesSites(k).gt.0)then
			do jj = 1,numSiteAtom(k)
				uexcessBerman(jj) = 0.0d0
				end do

			do 1210 jj = 1,numMargulesSites(k)
	
			do i1 = 1,numMargulesSiteCats(k,jj)
				uexcessBerman(i1) = 0.0d0
					do j1 = 1,numMargulesSiteCats(k,jj)
						do k1 = 1,numMargulesSiteCats(k,jj)
						MargulesWG(i1,j1,k1) = 0.0D0
						Windex(i1,j1,k1) = 0				! integer
						end do
					end do
				end do
			

			do L = 1,numMargulesSiteCats(k,jj)
				Xmarg(L) = SiteFraction(MargulesSiteCats(k,jj,L))
				uexcessTemp(L) = 0.0d0
				end do

			!write(*,*)' Margules calcs. Phase = ',k,phName(k)
			do L = 1,numMargulesWterms(k,jj)
				index = MargulesWindex(k,jj,L)
				index1 = index
				! write(*,*)L,index
				i1 = index1/100
				index1 = index1 - i1*100
				j1 = index1/10
				index1 = index1 - j1*10
				k1 = index1/1
				! MargulesWG(i1,j1,k1) = MargulesWH(k,jj,L) - MargulesWS(k,jj,L)*(TKtemp - TR) + MargulesWV(k,jj,L)*PbTemp
				MargulesWG(i1,j1,k1) = MargulesWH(k,jj,L) - MargulesWS(k,jj,L)*(TKtemp) + MargulesWV(k,jj,L)*PbTemp
				WIndex(i1,j1,k1) = index
				!write(12,*)index,MargulesWG(i1,j1,k1)
				end do

			!pause ' take a look and hit return'
			uexcessTemp = 0.0d0
			i1 = numMargulesSiteCats(k,jj)
			j1 = numMargulesWterms(k,jj)
			call BermanExcess(i1,j1,MargulesWG,WIndex,Xmarg,uexcessTemp)
			do i1 = 1,numMargulesSiteCats(k,jj)
				!MargulesSiteCats(phMax,numMargulesSitesMax,numMargulesSiteCatsMax)
				uexcessBerman(MargulesSiteCats(k,jj,i1)) = uexcessTemp(i1)
				end do
			!write(12,*)'Phase,site,uexcess = ',k,jj,uexcessBerman(jj)
1210			continue
			if(iLong(6).eq.1)then
				write(12,*)'uexcessBerman '
				write(12,1256)(SiteAtomName(k,jj),jj=1,numSiteAtom(k))
1256				Format(20(7x,A8))
				write(12,1257)(uexcessBerman(jj),jj=1,numSiteAtom(k))
1257				Format(20F15.3)
				endif

			endif		! end Margules calculations


!		****************************** phase component loop
! 		loop for every phase component
		Do 1020 j = 1,numPhCo(K)	! calculate  for each phase component
!		Calculate ideal activities
! 		If ActivityConstant = 0, then we are using a molecular model
		if(ActivityConstant(k,j).eq.0.0)then
			Lnatemp = Dlog(Xtemp(j))
			else
! 			Ideal activity part is made from expression
! 			R*TK*(ln(ActivityConstant) + alpha1*ln(X1) + alpha2*ln(X2) + alpha3*ln(X3) + .....)
			Lnatemp = Dlog(ActivityConstant(k,j))
			do 1030 L = 1,numSiteAtom(K)
			lnatemp = lnatemp + Alpha(k,j,L)*Dlog(SiteFraction(L))
1030			continue
			endif
	

!	Reciprocal free energy calculations
			uRecip = 0.0d0
			do  1110 ii = 1,numReciprocal(k)
			WGrecip = WHrecip(k,ii) - WSrecip(k,ii)*(TKtemp-TR) + WVrecip(k,ii)*Pbtemp
			do 1130 L = 1,numSiteAtom(K)
			if(RecipIndex(k,j,ii,L).eq.0)go to 1130
			if(RecipIndex(k,j,ii,L).eq.1)then		! either +1 or -1
				uRecipTemp = WGrecip * SiteFraction(L)
				else
				uRecipTemp = WGrecip * (1.0d0 - SiteFraction(L))
				endif				
1130			continue
			uRecipTemp = uRecipTemp * float(RecipConstant(k,j,ii))		! either +1 or -1
			uRecip = uRecip + uRecipTemp	
1110			continue


!		Add in Margules contribution
		uexcessMargules = 0.0d0
		do 1040 jj = 1,numMargulesSites(k)
		do 1046 i1 = 1,numMargulesSiteCats(k,jj)
		L = MargulesSiteCats(k,jj,i1)				! this is the cation number
		if(XtoSiteAtom(k,L,j).ne.0)then				! L is the site atom
			uexcessMargules = uexcessMargules + uexcessBerman(L)
			go to 1046					! loop out because we only add uexess once for each site
			endif
1046		continue
1045		continue
1040		continue



1250		continue		! jump here for melt models
		uMin(j) = lnatemp*RTK + uexcessMargules + uRecip

		if(iLong(6).eq.1)then
!			write(12,*)'Ph comp,lnatemp,lnatemp*RTK,uRecip,uexcess,utotal'
			write(12,1055)j,lnatemp,lnatemp*RTK,uRecip,uexcessMargules,uMin(j)
1055			Format('Component,lnatemp,lnatemp*RTK,uRecip,uexcess,utotal',I5,8E15.5)
			endif


1020		continue		! end of loop for this phase component

		iflag = 0
		return


	case (2)	! HP98 activity models

	!      	Calculate mole fractions of cations on sites
	!     	in terms of mole fractions of thermodynamic components
	! 	This is the same code as in subroutine XtoSite except here I am
	! 	using Xtemp for compositions instead of X because I need to calculate
	! 	the finite difference derivative w/r/t X and Xtemp is passed through from calling routine
		do 2011 L = 1,numSiteAtom(K)
		SiteFraction(L) = 0
		do 2010 j = 1,numPhCo(K)
		SiteFraction(L) = SiteFraction(L) + XToSiteAtom(k,L,j)*Xtemp(j)/SiteMultiplicity(k,L)
2010   		continue
2011   		continue
! 		Check to see that no atoms in use are less than zero
		do 2015 L = 1,numSiteAtom(K)
        	SiteFraction(L) = SiteFraction(L)+1.d-20	! if zero, make very slightly larger
 		if(SiteFraction(L).le.1.0D-30)go to 999
! 		if(SiteFraction(L).le.1.0D-30)then
!			SiteFraction(L) = 1.0d-30
!			endif
2015		continue
! 		If here, then all site atoms are > 0.0


		do 2030 n = 1,numSFterms(k)
!		WGHP11(n) = SFWH(k,n) - SFWS(k,n)*TKtemp + SFWV(k,n)*Pbtemp/1000.0d0	Wv is KJ/Kb or J/bar -- all my calcs are in J/bar
		WGHP11(n) = SFWH(k,n) - SFWS(k,n)*(TKtemp-TR) + SFWV(k,n)*Pbtemp
2030		continue

		if(iLong(6).eq.1)then
			write(12,*)'Min ',Minrec(K),' Nonideal terms = ',numSFterms(k)
			write(12,*)'  WH   WS    WV    ID1   ID2    WG'
			do 2045 n = 1,numSFterms(k)
			write(12,3043)SFWH(k,n),SFWS(k,n),SFWV(k,n),SFID(k,n,1),SFID(k,n,2),WGHP11(n)
2043			format(3F12.3,2I8,F15.3)
2045			continue
			endif

! 		loop for every phase component in the current problem
		Do 2020 j = 1,numPhCo(K)	! calculate  for each phase component
		select case (PhaseType(k))
		case(105)		! HP11 melt model ideal activities
!			This code REQUIRES that the first phase component in the melt is always H2O
!			ASF parameters for HP11 melt are all = 1, so Xtemp(j) = PhiAFS(j)
!			For the melt the site fraction is equal to the mole fraction of the component
			if(j.eq.1)then
				lnatemp = Xtemp(1)**2			! this one for H2OL
				else
				lnatemp = Xtemp(j)*(1.0d0-Xtemp(1))	! Xtemp2(1) is H2OL
				endif
				lnatemp = Dlog(lnatemp)
		case(103)		! HP11 ideal ionic activities
!			Calculate ideal activities
! 			If ActivityConstant = 0, then we are using a molecular model
			if(ActivityConstant(k,j).eq.0.0)then
				Lnatemp = Dlog(Xtemp(j))
				else
! 				Ideal activity part is made from expression
!	 			R*TK*(ln(ActivityConstant) + alpha1*ln(X1) + alpha2*ln(X2) + alpha3*ln(X3) + .....)
				Lnatemp = Dlog(ActivityConstant(k,j))
				do 2021 L = 1,numSiteAtom(K)
				lnatemp = lnatemp + Alpha(k,j,L)*Dlog(SiteFraction(L))
2021				continue
				endif
		case default
			write(*,*)' PhaseType error in thermodyanmic data file for phase:'
			write(*,*)k,phName(k)
			write(*,*)' PhaseType must be >200 for HP11 files. PhaseType for this phase = ',PhaseType(k)
			write(*,*)' This must be fixed and the program restarted'
			pause 'Hit return to continue'
		end select
	
	
		uexcess = 0.0d0
		! Equation 20 of Powell and Holland 1993 Am Mineral
		! and Holland and Powell 2003

!		----- We are calculating the uexcess for phase component j
!		Note that for the set of independent components, pkzero is either 0 or 1
!		However, for dependent components it may be a fraction 
!			(pkzero is the mole fraction of the dependent component in the component of interest)
!		But since Gibbs only deals with independent component sets, this coding should be fine
		tempsum = 0.0d0
		do n = 1,numSFterms(k)
			i1 = SFID(k,n,1)
			pkzero1 = 0.0d0
			if(j.eq.i1)pkzero1 = 1.0d0
			i2 = SFID(k,n,2)
			pkzero2 = 0.0d0
			if(j.eq.i2)pkzero2 = 1.0d0
			!tempsum = tempsum + WGHP11(n)*ASFtotal(j)*(pkzero1 - PhiAFS(i1))*(pkzero2 - PhiAFS(i2))
			tempsum = tempsum + WGHP11(n)*(pkzero1 - Xtemp(i1))*(pkzero2 - Xtemp(i2))
			end do
		uexcess = - tempsum   

2250		continue		! jump here for melt models
		uMin(j) = lnatemp*RTK + uexcess

		if(iLong(6).eq.1)then
			write(12,2055)j,exp(lnatemp),lnatemp,lnatemp*RTK,uexcess,uMin(j)
2055			Format('Component,act,lnatemp,lnatemp*RTK,uexcess,utotal',I5,5E15.5)
			endif

2020		continue	! end of loop for every phase component


		iflag = 0
		return



	case(3)		! for HP11  we need the size-paramater adjusted mole fractions of components
! 		This code is for excess models that use symmetric formalism (moles of components rather than cations on sites)

	!      	Calculate mole fractions of cations on sites
	!     	in terms of mole fractions of thermodynamic components
	! 	This is the same code as in subroutine XtoSite except here I am
	! 	using Xtemp for compositions instead of X because I need to calculate
	! 	the finite difference derivative w/r/t X and Xtemp is passed through from calling routine
		do L = 1,numSiteAtom(K)
			SiteFraction(L) = 0
			do j = 1,numPhCo(K)
				SiteFraction(L) = SiteFraction(L) + XToSiteAtom(k,L,j)*Xtemp(j)/SiteMultiplicity(k,L)
				end do
			end do

! 		Check to see that no atoms in use are less than zero
		do L = 1,numSiteAtom(K)
        		SiteFraction(L) = SiteFraction(L)+1.d-20	! if zero, make very slightly larger
 			if(SiteFraction(L).le.1.0D-30)go to 999
			end do
! 		If here, then all site atoms are > 0.0


! 		loop for every phase component in the current problem
		Do 3020 j = 1,numPhCo(K)	! calculate  for each phase component
		! I need to change this toXS model type so new melts can be added in the data file rather than the code here
		select case (PhaseType(k))
		case(205)		! HP11 melt model ideal activities
!			This code REQUIRES that the first phase component in the melt is always H2O
!			ASF parameters for HP11 melt are all = 1, so Xtemp(j) = PhiAFS(j)
!			For the melt the site fraction is equal to the mole fraction of the component
			select case (j)		! j is the component of the melt
				case(1)		! H2O
					lnatemp = Xtemp(1)**2			! this one for H2OL
				case(8)		! FaL
					Xfac = (1.0d0-Xtemp(1))
					XFe = Xtemp(8)/(Xtemp(8)+Xtemp(7))
					XpOl = Xtemp(8)+Xtemp(7)
					lnatemp = Xfac*XpOl*XFe**5
				case(7)		! FoL
					Xfac = (1.0d0-Xtemp(1))
					XFe = Xtemp(8)/(Xtemp(8)+Xtemp(7))
					XpOl = Xtemp(8)+Xtemp(7)
					lnatemp = Xfac*XpOl*(1.0d0 - XFe)**5
				case default
					Xfac = (1.0d0-Xtemp(1))
					lnatemp = Xtemp(j)*Xfac	! All other components
				end select					

!			lnatemp = Dlog(lnatemp)
			lnaIdeal(j) = Dlog(lnatemp)
!			pValue(j) = PhiAFS(j)
			pValue(j) = Xtemp(j)
		case(203)		! HP11 ideal ionic activities
!			Calculate ideal activities
! 			If ActivityConstant = 0, then we are using a molecular model
			if(ActivityConstant(k,j).eq.0.0)then
				Lnatemp = Dlog(Xtemp(j))
				else
! 				Ideal activity part is made from expression
!	 			R*TK*(ln(ActivityConstant) + alpha1*ln(X1) + alpha2*ln(X2) + alpha3*ln(X3) + .....)
				Lnatemp = Dlog(ActivityConstant(k,j))
				do L = 1,numSiteAtom(K)
					lnatemp = lnatemp + Alpha(k,j,L)*Dlog(SiteFraction(L))
					end do
				endif
			lnaIdeal(j) = lnatemp			! log is already taken
!			pValue(j) = PhiAFS(j)
			pValue(j) = Xtemp(j)

		case(206)		! HP11 Spinel activities
			! calculate ideal activities
! 			do l = 1,2
! 				SiteFractionSp(L) = SiteFraction(L)	! Site fractions are normalized to 1.0--These calculations require the number of cations on the site
! 				end do
! 			do l = 3,6
! !				SiteFractionSp(L) = 2.0D0*SiteFraction(L)	! Site fractions are normalized to 1.0--These calculations require the number of cations on the site
! 				SiteFractionSp(L) = SiteFraction(L)	! Site fractions are normalized to 1.0--These calculations require the number of cations on the site
! 				end do
			AlTiFe3TotSp = (SiteFraction(3) + SiteFraction(4) + 2.0d0*SiteFraction(5))	
!			xsp = (SiteFraction(2)+SiteFraction(6))/(SiteFraction(2) + SiteFraction(1) + SiteFraction(6))	! x - this includes Fe in the Al,Ti site -- WRONG!
			xsp = (SiteFraction(2))/(SiteFraction(2) + SiteFraction(1))		! x
			ysp = SiteFraction(3)/AlTiFe3TotSp					! y
			zsp = 2.0D0*SiteFraction(5)/AlTiFe3TotSp				! z
			XAlSp  = ysp								
			XFe3Sp = 1.0d0 - ysp - zsp
			XTiSp  = zsp				
			XMgSp = 1.0d0 - xsp
			XFeSp = xsp
! 			write(12,*)'Spinel'
! 			write(12,*)(SiteFraction(l),l=1,6)
! 			write(12,*)'xsp = ',xsp
! 			write(12,*)'ysp = ',ysp
! 			write(12,*)'zsp = ',zsp
! 			write(12,*)'XAlSp  = ',XAlSp			
! 			write(12,*)'XFe3Sp = ',XFe3Sp			
! 			write(12,*)'XTiSp  = ',XTiSp			
! 			write(12,*)'XMgSp  = ',XMgSp			
! 			write(12,*)'XFeSp  = ',XFeSp			
			select case(j)		! the phase component
				case(1)		! hercynite
				lnatemp = XFeSp * XAlSp				
				pValue(1) = ysp + (xsp - 1.0d0)*(1.0d0 + zsp)
				case(2)		! Spinel
				lnatemp = XMgSp * XAlSp				
				pValue(2) = (1.0d0 - xsp)*(1.0d0 + zsp)
				case(3)		! magnetite
				lnatemp = XFeSp * XFe3Sp
				pValue(3) = 1.0d0 - ysp - zsp
				case(4)		! ulvospinel
				lnatemp = XFeSp * XTiSp
				pValue(4) = zsp
				case default
					call fss_alert('Alert',' Case not found in spinel activity code')
				end select
			if(iLong(6).eq.1)then
				write(12,3056)j,lnatemp,pValue(j),(SiteFraction(L),L=1,6)
3056				Format('Component,IDact,pValue ',I5,12E16.7)
				endif

			lnaIdeal(j) = Dlog(lnatemp)
		case default
			write(*,*)' PhaseType error in thermodyanmic data file for phase:'
			write(*,*)k,phName(k)
			write(*,*)' PhaseType must be >200 for HP11 files. PhaseType for this phase = ',PhaseType(k)
			write(*,*)' This must be fixed and the program restarted'
			pause 'Hit return to continue'
		end select
	
3020	continue	! end loop on every phase component


		! Calculate the ASF scaling parameters and scale the compositions to these values
		sum = 0.0d0
		do j = 1,numPhCo(k)	! Loop through only the phase components being used
			ASFtotal(j) = (ASF(k,j,1) - ASF(k,j,2)*TKtemp + ASF(k,j,3)*Pbtemp)
!			sum = sum + Xtemp(j)*ASFtotal(j)		! Scale the composition to the ASF size
			sum = sum + pValue(j)*ASFtotal(j)		! Scale the composition to the ASF size
			end do
		! Now scale according to the size parameter - this is the term Phi in the equation at the bottom of the 
		! left column on page 493 of H&P 2003
		do j = 1,numPhCo(k)
!			PhiAFS(j) = Xtemp(j)*ASFtotal(j)/sum
			PhiAFS(j) = pValue(j)*ASFtotal(j)/sum
			end do
		if(iLong(6).eq.1)then
			Write(12,*)' ASF parameters (1,2,3, ASFtotal, Xtemp, PhiAFSScaled)'
			do j = 1,numPhCo(k)
				write(12,3024)ASF(k,j,1),ASF(k,j,2),ASF(k,j,3),ASFtotal(j),xtemp(j),PhiAFS(j)
3024				format(3F10.3,4F15.7)
				end do
			endif
		! At this point the array PhiAFS(j) contains the size-parameter (Phi) scaled compositions of each phase component
		! I should change the array name from PhiAFS to Phi
		! Calculate the WG terms at T and P
		do n = 1,numSFterms(k)
			WGHP11temp(n) = SFWH(k,n) - SFWS(k,n)*TKtemp + SFWV(k,n)*Pbtemp
			end do

!		Apply size parameter scaling to W terms Wij* = Wij/(ASF(i)+ASF(j))	= page 493 of HP 2003
! 			note that we will multiply by ASF(phase component of interest) below
		do n = 1,numSFterms(k)
			WGHP11(n) = 2.0d0*WGHP11temp(n) /(ASFtotal(SFID(k,n,1)) + ASFtotal(SFID(k,n,2)))
			end do

		if(iLong(6).eq.1)then
			write(12,*)'Min ',Minrec(K),' Nonideal terms = ',numSFterms(k)
			write(12,*)'  WH            WS           WV        ID1     ID2       WGraw        WGASFscaled'
                                   !25100.000     -10.800       0.338       1       2      42966.622
			do n = 1,numSFterms(k)
				write(12,3043)SFWH(k,n),SFWS(k,n),SFWV(k,n),SFID(k,n,1),SFID(k,n,2),WGHP11temp(n),WGHP11(n)
3043				format(3F12.3,2I8,2F15.3)
				end do
			endif

		! Calculate RTlnGamma	
		do j = 1,numPhCo(K)
			uexcess = 0.0d0			! this is RTlnGamma
			! This code is for the ASF model of
			! Equation 20 of Powell and Holland 1993 Am Mineral
			! and equation 2 Holland and Powell 2003

			!----- We are calculating the uexcess for phase component j
			! Note that for the set of independent components, pkzero is either 0 or 1
			! However, for dependent components it may be a fraction 
			! 	(pkzero is the mole fraction of the dependent component in the component of interest)
			! But since Gibbs only deals with independent component sets, this coding should be fine
			tempsum = 0.0d0			! temporary sum for RTlnGamma
			do n = 1,numSFterms(k)
				i1 = SFID(k,n,1)
				pkzero1 = 0.0d0
				if(j.eq.i1)pkzero1 = 1.0d0
				i2 = SFID(k,n,2)
				pkzero2 = 0.0d0
				if(j.eq.i2)pkzero2 = 1.0d0
				tempsum = tempsum + WGHP11(n)*ASFtotal(j)*(pkzero1 - PhiAFS(i1))*(pkzero2 - PhiAFS(i2))
!				tempsum = tempsum + WGHP11(n)*ASFtotal(j)*(pkzero1 - pValue(i1))*(pkzero2 - pValue(i2))
				if(iLong(6).eq.1)then
					write(12,*)' j,n,i1,pkzero1,i2,pkzero2,tempsum ',j,n,i1,pkzero1,i2,pkzero2,tempsum
					endif
				end do
			uexcess = - tempsum   
			RTlnGamma(j) = uexcess
	!		uMin(j) = lnatemp*RTK + uexcess
			uMin(j) = lnaIdeal(j)*RTK + uexcess
			if(iLong(6).eq.1)then
				write(12,3055) &
				     j,Dexp(lnaIdeal(j)),lnaIdeal(j),lnaIdeal(j)*RTK,RTlnGamma(j),exp(RTlnGamma(j)/RTK),uMin(j)
3055				Format('Component,IDact,lnIDact,lnIDact*RTK,RTlnGamma,Gamma,utotal',I5,10E15.5)
				endif
			end do


	iflag = 0
	return

	case default

	write(*,*)' Bad dataset type (DatasetKey). Please fix the thermodynamic data file'
	Pause 'Abort program'
	stop 

	end select


250	continue



999	continue	! a site atom has a zero or negative value
        iflag=1
	if(iLong(1).eq.1)then

	      	write(*,*)' In subroutine Activity'
		write(*,*)' A cation on a site has a zero or negative value'
		Write(*,*)PhName(K)
		write(*,*)'Atom      value'
		do L = 1,numSiteAtom(K)
			write(*,*)SiteAtomName(k,l),SiteFraction(L)
			end do
		write(*,*)
		write(*,*)'Component    value'
		do j = 1,numPhCo(K)
			write(*,*)phCoName(k,j),Xtemp(j),xPhCo(k,j)
			end do
	 	write(*,*)'Tc,Pb= ',TC,PB
	      	write(*,*) 'Hit return to continue'
	      	pause 'Line 458 in Sub Activity'
		endif
        return


      end
