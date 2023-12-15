! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE LeastSquares_H()
!	use MatrixArrays
      implicit none
!
!     written August 2017 by F. Spear
!     ROUTINE TO calculate enthalpies of phase components from a set of input datasets
!	using a weighted least squares routine
!
! ****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
!	include "Solute.inc"
	include 'LSQ.inc'
! ****************************************
!     local variables
	integer*4 i,ii,j,jj,k,ier,izero,ivol,numSolve,kCur,kk
	real*8 a(MaxEq,MaxEq),YY(MaxEq),Ysigma(MaxEq),soln(MaxEq),Yweight,AMasterRxn(MaxEq,MaxEq)
	integer*4 indexPhCoNames(MaxEq),MasterSolve(MaxEq)
	character*16 MasterPhCoNames(MaxEq),MIFPhCoNames(MaxEq)
	REAL*8 R,TEMPH,TEMPS,TEMPK
	character*128 Title,MasterMIFfile,MIFfile
	!character*128 filein
	integer*4 status,inew,numMasterPhCo,numEquationsThisMIF,numEquationsTotal,numMIFPhCo,ierr
!
	DATA R/8.3144D0/
!

!	First read in a MIF that contains all of the phase components that will be considered

	write(*,*)' Solve for H(Tref,Pref) using a least-squares approach'
	write(*,*)' '
	write(*,*)' There is a configuration file that contains the names of all of the files that are used in this routine'
	write(*,*)' That file is opened and the rest of the files should open automatically'
	write(*,*)' The configuration file contains the names of these files....'
	write(*,*)' 	(1) A  MIF that contains all of the phases and phase components'
	write(*,*)'    		under consideration for this problem. This is used to set up the arrays.'
	write(*,*)'    		The compositions of the phases in the MIF do not matter'
	write(*,*)'	(2) The name of the first input file. This (and subsequent) files contain '
	write(*,*)'		(a) the P and T of the experiment or natural assemblage'
	write(*,*)'		(b) the phases and compositions of phases in the assemblage'
	write(*,*)'		Note that every dataset used in the least squares calculations is defined in a separate input file'
	write(*,*)'	(3) The next input file...'
	write(*,*)'	(4) etc until EOF'
	write(*,*)' '
	Write(*,*)' IMPORTANT: all of the above files must be in the same folder'
	write(*,*)' '
	write(*,*)' '
	write(*,*)' Open the configuration file'
	open(31,file='',status='old',iostat=status)
      	if(status.ne.0)then
      		write(*,*)' Configuration file did not open'
		pause ' hit return to continue'
      		return
		endif
	read(31,*)Title
	read(31,*)MasterMIFfile
	write(*,*)' Contents of MasterMIFfile'
	write(*,*)title
	write(*,*)MasterMIFfile
      	OPEN(5,FILE=MasterMIFfile,STATUS='OLD',action='read',iostat=status)
 !       INQUIRE (5, NAME=filein)
!        write(*,*)filein
      	if(status.ne.0)then
      		write(*,*)' Master MIF file did not open'
		pause ' hit return to continue'
      		return
		endif
	inew = 0
	Call BEGIN (2,INEW)	! 2 says file is already open
	write(12,*)'Back from Sub Begin - line 75 of Sub LeastSquares_for_H'
	if (inew.eq.0)then
		write(*,*)' Problem reading Master MIF file'
		pause ' hit return to continue'
		return
		endif	

	do 9 i = 1,16		
	deltax(i) = 0.		! zero out values
9	continue
	nstep = 1
!	WRITE(*,*)NVAR,NEQ
!	PAUSE ' HIT RETURN'
!	NVAR = 1
!	NEQ = 0
	write(12,*) ' MasterMIF file'
	write(12,*)MasterMIFfile
          call printt(0)        !option 1 prints to the file and screen

!	make a list of all of the phase components in the MasterMIF file
	WRITE(12,*)'Phase components in this problem'
	WRITE(*,*)'Phase components in this problem'
	jj = 0
	Do 10 kCur = 1,numPh
	k = asmCurrent(kCur)
	Do 10 j=1,numPhCo(k)      
!	ier = ier + 1
	jj = jj + 1
	MasterPhCoNames(jj) = phCoName(k,j)
	write(12,*)jj,MasterPhCoNames(jj),phCoName(k,j),k,j
	write(*,*)jj,MasterPhCoNames(jj),phCoName(k,j),k,j
10   	continue
	numMasterPhCo = jj
	i = 0
	write(*,*)'Specify which components you wish to solve for'
!	write(*,*)'Phase-K, Component-J (2 integers)'
!	write(*,*)' End the list with two zeros (0,0)'
	write(*,*)' A single integer (left column) '
	write(*,*)' End the list with a zero (0)'
11	continue
	i = i + 1
!	read(*,*)kSolve(i),jSolve(i)
	read(*,*)MasterSolve(i)
!	read(*,*)j
!      	Nsolve(i)=j
!	if(kSolve(i).eq.0.and.jSolve(i).eq.0)then
	if(MasterSolve(i).eq.0)then
		numSolve = i - 1
		go to 12
		endif
!      	write(12,*)'Phase, Component ',kSolve(i),jSolve(i)
!      	write(*,*)'Phase, Component ',kSolve(i),jSolve(i)
      	write(12,*)'PhaseComponent ',MasterSolve(i),MasterPhCoNames(MasterSolve(i))
      	write(*,*)'PhaseComponent ',MasterSolve(i),MasterPhCoNames(MasterSolve(i))
	go to 11
12	continue	

!     ZERO MASTER ARRAY
!
      DO 15 J=1,MaxEq
      DO 15 I=1,MaxEq		! we don't know how many equations there will be until all files are read in
      				! so just zero them all out
15    A(I,J)=0.0D0


	numEquationsThisMIF = 0	! this will count the total number of equations in this MIF
	numEquationsTotal = 0	! count the total number of equations

!------------------------------------------------------------
!	Now read in the MIFs for each assemblage
!	We loop back here for every MIF
20	continue
	read(31,*,end = 99)MIFfile,Yweight
	write(*,*)MIFfile,Yweight
	write(12,*)'          '
	write(12,*)'          '
	write(12,*)'#####################################################'
	write(12,*)'#####################################################'
	write(12,*)'New MIF input file'
	write(12,*)MIFfile,Yweight
      	OPEN(5,FILE=MIFfile,STATUS='OLD',action='read',iostat=status)
      	if(status.ne.0)then
      		write(*,*)' MIF file did not open'
		write(*,*)' MIF file name = ',MIFfile
		pause ' hit return to continue'
      		return
		endif
	inew = 0
	Call BEGIN (2,INEW)	! 2 says file is already open
	if (inew.eq.0)then
		write(*,*)' Problem reading MIF file'
		write(*,*)' MIF file name = ',MIFfile
		pause ' hit return to continue'
		return
		endif	
          call printt(0)        !option 1 prints to the file and screen
!         COMPUTE REACTION MATRIX

	  iLong(7) = 0
          CALL REXN
          iLong(7) = 0

	WRITE(12,*)'Phase components in this problem'
	jj = 0
	Do 30 kCur = 1,numPh
	k = asmCurrent(kCur)
	Do 30 j=1,numPhCo(k)      
!	ier = ier + 1
	jj = jj + 1
	MIFPhCoNames(jj) = phCoName(k,j)
	write(12,*)k,j,MIFPhCoNames(jj),phCoName(k,j)
30   	continue
	numMIFPhCo = jj

	write(12,*)' '
	write(12,*)' Linearly independent reactions:'
	write(12,2110)(MIFPhCoNames(J),J=1,NP)
2110     FORMAT(' ',4X,30(A8))
	
      	ii = np-nrx      !note that reaction coefficients are stored in the end of array ARX
	DO 429 I=1,nrx
	write(12,2120)(ARX(i+ii,J),J=1,NP)
2120     FORMAT(' ',30F8.3)
429      CONTINUE

!	construct index array that relates MIF phase components to MasterMIF phase components
	do 50 i = 1,numMIFPhCo
	do 51 j = 1,numMasterPhCo
	if(MIFPhCoNames(i).eq.MasterPhCoNames(j))then
		indexPhCoNames(i) = j
		go to 50
		endif
51	continue
	write (*,*)' Problem setting up index array. Could not find component in MasterMIF'
	write(*,*)' Absent component = ',MIFPhCoNames(i)
	pause 'Hit return to continue'
50	continue


!
      TK=TC+273.15D0
      numEquationsThisMIF = NRX
!
!     THESE CALLS CALCULATE THE THERMODYNAMIC
!     PROPERTIES OF THE PHASES OF INTEREST
!
!     SUBROUTINE COMPUTES LnActivity,G,H,S,AND V OF ALL PHASE components
      izero=0
      ivol=0            !compute everything
      K=0
      CALL AllKduTPX(ivol,izero)
      if(izero.eq.1)then
		write(*,*)' Problem calculating thermodynamic properties in ALLKduTPX'
		pause 'Hit return to continue'
		return
		endif

350    CONTINUE
	Do 390 kCur = 1,numPh
	k = asmCurrent(kCur)
        write(12,*)
        write(12,200)MINREC(K),SITMUL(K),PHNAME(K)
200     FORMAT(' ',I10,F8.1,4X,A32)
        write(12,'(1X,A15,1F10.3)')'MOLAR VOLUME =',VMOL(K)
        write(12,501)(phCoName(k,J),J=1,numPhCo(K))
501     FORMAT('                   ',20(A12,3X))
        write(12,504)(xPhCo(k,J),J=1,numPhCo(K))
504     FORMAT('  COMP    ',20(F15.4))
        write(12,505)(Dexp(lnAct(k,J)),J=1,numPhCo(K))
505     FORMAT(' Activ    ',20(E15.3))
        write(12,506)(lnAct(k,J),J=1,numPhCo(K))
506     FORMAT('  Ln(a)   ',20(F15.5))
          write(12,514)(hPhCoZero(k,J),J=1,numPhCo(K))
514       Format(' HPhCozero',20E15.7)
          write(12,515)(HATTP(k,J),J=1,numPhCo(K))
515       Format(' HatTP    ',20E15.7)
          write(12,509)(sPhCoZero(k,J),J=1,numPhCo(K))
509       Format(' Szero    ',20F15.5)
          write(12,510)(SATTP(k,J),J=1,numPhCo(K))
510       Format(' SatTP    ',20F15.5)
          write(12,508)(vPhCoZero(k,J),J=1,numPhCo(K))
508       Format(' Vzero    ',20F15.5)
          write(12,513)(VATTP(k,J),J=1,numPhCo(K))
513       Format(' VatTP    ',20F15.5)
          write(12,516)(GATTP(k,J),J=1,numPhCo(K))
516       Format(' GatTP    ',20E15.7)
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
      	YY(numEquationsTotal + i)=(TEMPH - TK * TEMPS + R*TK*TEMPK)
1010	continue

!     Now choose the phase components to solve for H


!     Now construct a matrix (AA) with the appropriate coefficients for each phase component to solve.

!      numSolve (the number of H values to solve for) is input above
      do 1105 i=1,nrx		! loop through every equation in this assemblage
      TempH = 0
	jj = 0
	do 1110 kCur = 1,numPh
	k = asmCurrent(kcur)
	do 1112 j = 1,numPhCo(k)
	jj = jj + 1
	AMasterRxn(numEquationsTotal + i,indexPhCoNames(jj)) = ARX(ii+i,jj)
	do 1120 kk = 1,numSolve		! Set up A matrix with stoichiometric coefficients for each Href to solve for
	if(MasterSolve(kk).eq.indexPhCoNames(jj))then
!	if(kSolve(kk).eq.k.and.jSolve(kk).eq.j)then
!	        A(i,kk) = -ARX(ii+i,jj)       !the A matrix coefficient is the reaction coefficient
!	        A(numEquationsTotal + i,kk) = -ARX(ii+i,jj)       !the A matrix coefficient is the reaction coefficient
!	        TempH = TempH + ARX(ii+i,jj)*hPhCoZero(k,j)    !Subtract hPhCoZero from data vector because we need to solve for this
                                                !Note it was added in above already
!	        A(numEquationsTotal + i,indexPhCoNames(kk)) = -ARX(ii+i,jj)       !the A matrix coefficient is the reaction coefficient
	        A(numEquationsTotal + i,kk) = -ARX(ii+i,jj)       !the A matrix coefficient is the reaction coefficient
	        TempH = TempH + ARX(ii+i,jj)*hPhCoZero(k,j)    !Subtract hPhCoZero from data vector because we need to solve for this


		endif
1120  	continue
1112	continue
1110  	continue
!	Adjust the YY(solution) vector to remove the Href for the components to solve for
	YY(numEquationsTotal + i) = YY(numEquationsTotal + i) - TempH
!	Ysigma(numEquationsTotal + i) = Yweight		! the weighting factor for every equation in this MIF is the same
							!  The idea is that the weight goes to the uncertainty for a particular assemblage
!	Ysigma(numEquationsTotal + i) = 1.0d0		! the weighting factor for every equation in this MIF is the same
	Ysigma(numEquationsTotal + i) = Yweight*YY(numEquationsTotal+1)		!Ysigma is calculated as a fraction of YY. The fraction is Yweight.		
1105  	continue

	numEquationsTotal = numEquationsTotal + nrx
        go to 20	! loop back and read the next MIF

99	continue	! when we get here, we have hit EOF reading unit 31 (the config file) and all of the data have been read in    
	close (31)

	write(12,*)' '
	write(12,*)' '
	write(12,*)'******************************************'
	write(12,*)'Master Reaction coefficient matrix'
	write(12,56)(MasterPhCoNames(j),j=1,numMasterPhCo)
56	format(' ',4X,30(A8,2X))
	do 55 i = 1,numEquationsTotal
	write(12,57)(AMasterRxn(i,j),j=1,numMasterPhCo)
57	format(' ',30F10.3)
55	continue
	write(12,*)'******************************************'
	write(12,*)' '
	write(12,*)' '



	Call LinearLSQ(numEquationsTotal,numSolve,a,YY,Ysigma,soln,ierr)

	IF (IER.EQ.1) then
		WRITE(*,*)' ************ ERROR **************************'
		WRITE(*,*)' Matrix failed to invert in SUBROUTINE LinearLSQ'
		write(*,*)' '
		pause 'Hit return to continue...'
		return
		endif

!     Write out solution vector
	write(12,*)
	write(12,*)'Phase components        H'
	Do 1605 i=1,numSolve
!	write(12,*)phCoName(kSolve(i),jSolve(i)),soln(i)
	write(12,*)MasterPhCoNames(MasterSolve(i)),soln(i)
1605  continue
	write(12,*)
	write(*,*)' Results are in the output window'
	pause 'Take a look and hit return'

	RETURN
	END

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine LinearLSQ(numrows,numcolumns,aMatrix,Y,Ysigma,soln,ierr)
	implicit none
!	Subroutine to calculate linear least squares
!	A * x = Y
!	x = (AT*A)-1 * AT * Y
	include 'LSQ.inc'
	integer*4 numrows,numcolumns,i,j,k,l,ierr,ii
	real*8 aMatrix(MaxEq,MaxEq),aTrans(MaxEq,MaxEq),aTa(MaxEq,MaxEq),ataOriginal(MaxEq,MaxEq),aTaaT(MaxEq,MaxEq)
	real*8 cormat_data(MaxEq,MaxEq),cormat(MaxEq,MaxEq)
	real*8 Y(MaxEq),soln(MaxEq),model(MaxEq),residual(MaxEq),Ysigma(MaxEq),weight(MaxEq),solnsigma(MaxEq)
	real*8 Ywt(MaxEq),aWt(MaxEq,MaxEq)
	real *8 sum,deta,chi2,chi2reduced,avgwt,sigma2

!	Compute weighted Y and A matrix
	do 20 i = 1,numRows
	Ywt(i) = Y(i)/Ysigma(i)
	do 20 j = 1,numColumns
	Awt(i,j) = aMatrix(i,j)/Ysigma(i)
20	continue

	write(12,*)' '
	write(12,*)' '
	write(12,*)' '
	write(12,*)' '
	write(12,*)'***********************************************************************'
	write(12,*)'***********************************************************************'
	write(12,*)'Least squares routine'
	write(12,*)'Input Data'
	write(12,*)'            Y vector    Ysigma     A matrix.....'
	do 82 i = 1,numRows
!	write(12,185)(aMatrix(l,k),k=1,numColumns),y(l)
	write(12,185)Y(i),Ysigma(i),(aMatrix(i,j),j=1,numColumns)
!181	format(A8,30F12.5)
82	continue

	write(12,*)''
	write(12,*)'Weighted matrix and data vectors'
!	write(12,*)'Input Data = A matrix , Y vector'
	write(12,*)'            Y vector    Ysigma     A matrix.....'
	do 84 i = 1,numRows
!	write(12,185)(awt(i,j),j=1,numColumns),y(i)
	write(12,185)Ywt(i),Ysigma(i),(Awt(i,j),j=1,numColumns)
!181	format(A8,30F12.5)
84	continue

!	compute aTranspose
	do 100 i = 1,numRows
	do 100 j = 1,numColumns
!	aTrans(j,i) = aMatrix(i,j)		! unweighted code
	aTrans(j,i) = aWt(i,j)			! weighted code
100	continue

!	compute At*A
	do 110 i = 1,numColumns
	do 112 j = 1,numRows
	sum = 0.0d0
	do 115 k = 1,numRows
!	sum = sum + aTrans(i,k)*aMatrix(k,j)	! unweighted code
	sum = sum + aTrans(i,k)*aWt(k,j)	! weighted code
115	continue
	aTa(i,j) = sum			! this one gets modified in Minvers and is returned as the inverse
	aTaOriginal(i,j) = sum			
112	continue
110	continue

!	write(12,*)'A transpose matrix'
!	do 83 k = 1,numColumns
!	write(12,185)(aTrans(k,l),l=1,numRows)
185	format(T9,30E15.5)
!83	continue
	write(12,*)'ATA matrix'
	do 87 k = 1,numColumns
	write(12,185)(aTaOriginal(k,l),l=1,numColumns)
87	continue


!	CALCULATE CORRELATION MATRIX, CORMAT (equation 11-28 of Bevington)
	DO 2137 I=1,numColumns
	DO 2137 J=1,numColumns
	CORMAT_data(I,J)=aTa(I,J)/SQRT(aTa(I,I)*aTa(J,J))
2137	continue

	write(12,*)' Correlation matrix'
	DO 5055 i=1,numColumns
5055	write(12,5005)(Cormat(i,ii),ii=1,numColumns)


!	Compute Inverse of aTa
!	Subroutine Minverse(A,Ma,DetA,ier)

	if(numColumns.eq.1)then
		aTa(1,1) = 1.0d0/aTa(1,1)	! with only 1 unknown the inverse is just this
		else
		ierr = 0
		call Minverse(aTa,numColumns,DetA,ierr)
		if(ierr.eq.1)then
			write(*,*)'Matrix failed to invert in routine Minverse'
			write(*,*)' Called from LinearLSQ'
			pause 'Hit return to continue'
			ierr = 1
			return
			endif
		endif
		
!	Compute (ATA)-1 * AT
	do 130 i = 1,numColumns
	do 132 j = 1,numRows
	sum = 0.0D0
	do 134 k = 1,numColumns
	sum = sum + aTa(i,k)*aTrans(k,j)	
134	continue
	aTaaT(i,j) = sum
132	continue
130	continue
!	compute solution vector
	do 140 i = 1,numColumns
	sum = 0.0d0
	do 142 j = 1,numRows
	sum = sum + aTaaT(i,j)*yWt(j)
142	continue
	soln(i) = sum
140	continue
!	calculate model and residuals
	do 150 i = 1,numRows
	sum = 0.0d0
	do 152 j = 1,numColumns
	sum = sum + aMatrix(i,j)*soln(j)
152	continue
	model(i) = sum
	residual(i) = y(i) - model(i)
150	continue

!	calculate Chi2 for fit (see equation 11.3 of Bevington)
	sum=0
	do 2240 j=1,numRows
	sum=sum + ((1./Ysigma(j))*Residual(j))**2
2240	continue
	Chi2=sum
	Chi2reduced=Chi2/(numRows-numColumns)

!	calculate sigma squared (sigma2) for fit (see equn 11.1 of Bevington)
!	first calculate weighting factors for each datum (equn 11.2 of Bevington)
!	compute average weight
	sum=0
	do 2250 j=1,numRows
	sum=sum+(1./Ysigma(j))**2
2250	continue
	avgwt=sum/numRows
!	now compute weight for each datum
	do 2255 j=1,numRows
	weight(j)=((1./Ysigma(j))**2)/avgwt
2255	continue

!	Now we can compute sigma squared
	sum=0
	do 2260 j = 1,numRows
	sum=sum+weight(j)*Residual(J)**2	
2260	continue
	Sigma2=sum*(1./(numRows-numColumns))	!num-nvar = variance of system



!	CALCULATE ERRORS ON PARAMETER ESTIMATES
	DO 2330 i=1,numColumns
2330	SolnSigma(I)=DSQRT(aTa(i,i))	! note that aTa is now the inverse of aTa

!	CALCULATE CORRELATION MATRIX, CORMAT (equation 11-28 of Bevington)
	DO 2350 I=1,numColumns
	DO 2350 J=1,numColumns
	CORMAT(I,J)=aTa(I,J)/DSQRT(aTa(I,I)*aTa(J,J))	! note that aTa is now the inverse of aTa
2350	continue


	write(12,*)'(ATA)-1 (inverse) matrix'
	do 85 k = 1,numColumns
	write(12,185)(aTa(k,l),l=1,numColumns)
85	continue
	write(12,*)'(ATA)-1*AT matrix'
	do 86 k = 1,numColumns
	write(12,185)(aTaaT(k,l),l=1,numRows)
86	continue
	write(12,*)' '
	write(12,*)'---------------------'
	write(12,*)'Solution = (ATA)-1*AT*Y and parameter errors'
	write(12,185)(soln(l),l=1,numColumns)
	write(12,185)(solnSigma(l),l=1,numColumns)
	write(12,*)'---------------------'
	write(12,*)' '
!	write(12,*)'Model = A*soln'
!	write(12,185)(Model(l),l=1,numRows)
!	write(12,*)'Residual = Y - Model'
!	write(12,185)(Residual(l),l=1,numRows)

	write(12,*)
	write(12,*)' Covariance matrix'
	do 5050 i=1,numColumns
5050	write(12,5005)(aTa(i,ii),ii=1,numColumns)	! note that aTa is now the inverse of aTa

	write(12,*)
	write(12,*)'         Chi^2 for fit   = ',Chi2
	write(12,*)' Reduced Chi^2 for fit   = ',Chi2reduced
	write(12,*)'         Sigma^2 for fit = ',sigma2

	write(12,*)
	WRITE (12,*)'      Y          YMODEL        RESID        ERROR            RELATIVE(%)'
	DO 5060 J=1,NumRows
	WRITE (12,5005)Y(j),model(j),Residual(j),Ysigma(j),100*Dabs(Residual(j))/Ysigma(j)
5060	continue
5000	FORMAT (20F13.5)
5005	format(20e13.5)


	return
	end



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Subroutine Minverse(A,Ma,DetA,ier)
!	routine to invert the matrix A(Ma,Ma)
!	Ma is the part of A that is used
!	DetA is the determinant
!	ier is an error switch
!		ier=0 matrix inverted ok
!		ier=1 matrix is singular (failed to invert)
!	maximum size of matrix that can be inverted is 24x24 without changing
!		dimension statement
	implicit none
	include 'LSQ.inc'
	real*8 A(MaxEq,MaxEq),DetA,Piv,piv1,Temp,Xq,Test
	integer*4 i,j,k,L,m,ijkl,ir(MaxEq),ic(MaxEq),MA,ier
! ------------------------------------------------------
	IER=0		!initialize error switch to 0 (ok)
	DO 1 I=1,MA
	IR(I)=0
	IC(I)=0
1	continue
	DETA=1.0
	DO 123 IJKL=1,MA
!	this code finds the largest pivot element
	I=0
	J=0
	TEST=0.0
	 DO 105 K=1,MA
	IF (IR(K).NE.0) GO TO 105
	DO 104 L=1,MA
	IF(IC(L).NE.0) GO TO 104
	XQ=DABS(A(K,L))
	IF(XQ.LT.TEST) GO TO 104
	I=K
	J=L
 	TEST=XQ
104	CONTINUE
105 	CONTINUE
!	Now the pivot element is A(i,j)
	PIV=A(I,J)
	DETA=PIV*DETA
	IF (PIV.EQ.0.0) then	!if true then matrix is singular-abort
		ier=1
		return
		endif
	IR(I)=J
	IC(J)=I
	PIV=1.0D0/PIV
	DO 5 K=1,MA
 	A(I,K)=A(I,K)*PIV
5	continue
	A(I,J)=PIV
	DO 9 K=1,MA
	IF (K.EQ.I) GO  TO 9
	PIV1=A(K,J)
 	DO 8 L=1,MA
 	A(K,L)=A(K,L)-PIV1*A(I,L)
8	continue
	A(K,J)=PIV1
9	CONTINUE
	PIV1=A(I,J)
	DO 11 K=1,MA
	A(K,J)=-PIV*A(K,J)
11	continue
	A(I,J)=PIV1
123	CONTINUE
	DO 16 I=1,MA
	K=IC(I)
	M=IR(I)
	IF (K.EQ.I) GO TO 16
	DETA=-DETA
	DO 14 L=1,MA
	TEMP=A(K,L)
	A(K,L)=A(I,L)
	A(I,L)=TEMP
14	continue
	DO 15 L=1,MA
	TEMP=A(L,M)
	A(L,M)=A(L,I)
 	A(L,I)=TEMP
15	continue
	IC(M)=K
	IR(K)=M
16	CONTINUE
	RETURN
	END
	
