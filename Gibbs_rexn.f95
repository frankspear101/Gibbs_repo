! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE REXN
      implicit none
!     ROUTINE TO COMPUTE A LINEARLY INDEPENDENT SET OF STOICHIOMETRIC
!     RELATIONS AMONG PHASE COMPONENTS

! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
! c*****************************************
! --------------------------------------------------------------------
	include "Tangent.inc"
! --------------------------------------------------------------------

!     Local arrays
!     DIMENSIONED FOR 12 SYSTEM COMPONENTS AND 25 PHASE COMPONENTS
      INTEGER*4 L,j,jj,k,i,ll,kCur,jdep
      REAL*8 A(maxPC,maxPC),TEMP,FACT
      common /REXNcommon/A,Temp,Fact
      	Character VN3(maxPC)*4
	common/VN3names/VN3
! c*****************************************
!     Set up array of phase component names
	jdep = 1	! dependent phase component is always 1
	L=0
	jj = 0
!	do 5002 k = 1,numPh
	Do 5002 kCur = 1,numPh
	k = asmCurrent(kCur)
	DO 5002 J=1,numPhCo(k)
	if(OSPhase(k).eq.0)then
		jj = jj + 1
	        VN3(jj)=phCoName(k,J)
		else
		if(j.gt.1)then
			jj = jj + 1
		        VN3(jj)='*'//phCoName(k,J)
			endif
		endif		
5002  CONTINUE

!     THIS Routine IS DESIGNED TO COMPUTE A LINEARLY INDEPENDENT
!     SET OF REACTIONS AMONG PHASE COMPONENTS IN A SYSTEM OF
!     NC SYSTEM COMPONENTS AND NP PHASE COMPONENTS
!     USES METHOD OF GAUSS-JORDAN REDUCTION OPERATING ON AN IDENTITY
!     MATRIX

	IF (NP.GT.maxPC) THEN
		WRITE(*,*) 'Dimension of array ARX in subroutine REXN is too small -----stop'
		pause 'Hit return to continue'
		stop
		endif
	L = 0
!      DO 10 L=1,NP
!	do 10 k = 1,numPh
	Do 12 kCur = 1,numPh
	k = asmCurrent(kCur)
	do 13 j = 1,numPhCo(k)
	if(OSPhase(k).eq.0)then
!		EQ phase
		L = L+1
        	DO 10 ll=1,NC
        	A(L,ll)=COMP(k,j,ll)
10    		CONTINUE
		else
!		OS phase
		L = L+1
        	DO 11 ll=1,NC
        	A(L,ll)=COMP(k,j,ll)-COMP(k,jdep,ll)
11    		CONTINUE
		ENDIF		

13	continue
12	continue

				
!     	ECHO INPUT DATA
      	IF(iLong(7).eq.1)then
483      	CONTINUE
         	write(12,*)
         	write(12,*)
         	write(12,*)'INPUT DATA'
         	write(12,*)
         	write(12,*)
         	write(12,2000)(CONAME(I),I=1,NC)
2000     	FORMAT(' ',10X,20(4X,A4))
         	DO 14 I=1,NP
            	write(12,2001)VN3(I),(A(I,J),J=1,NC)
2001        	FORMAT(' ',A8,2X,20(F8.3))
14       	CONTINUE

		endif
		
		
      	DO 50 J=1,NP
        DO 50 I=1,NP
50      ARX(I,J)=0.
      	DO 51 J=1,NP
51      ARX(J,J)=1.
      	DO 392 J=1,NC
        JJ=J+1
        IF(DABS(A(J,J)).GT..00001)GO TO 2307
        DO 2006 K=J+1,NP
2006    IF(DABS(A(K,J)).GT.0.00001)GO TO 2106

! 	If here, then the assemblage is degenerate. Sometimes this does not cause a problem
! 	but sometimes it does.
! 	If we are doing a "PseudoFromScratch", then we want to avoid these degeneracies.
	If (ipseudo.eq.1)then
		idegen = 1
		return
	endif
         WRITE(*,*)' WARNING----'
         WRITE(*,*)' The number of system components specified is'
         WRITE(*,*)' greater than the number required to describe'
         WRITE(*,*)' the chemical variability of the phases'
         WRITE(*,*)' (the last column of array A in SUBROUTINE REXN'
         WRITE(*,*)' is all zero. COL=',J
         WRITE(*,*)' NOTE: usually, this does not cause a problem'
         write(*,*) 'Hit return to continue...'
         read(*,*)
         GO TO 392
2106     CONTINUE
!        SWITCH ROW K WITH ROW J
         DO 225 I=1,NC
            TEMP=A(J,I)
            A(J,I)=A(K,I)
225         A(K,I)=TEMP
         DO 2306 I=1,NP
            TEMP=ARX(J,I)
            ARX(J,I)=ARX(K,I)
2306        ARX(K,I)=TEMP
2307     CONTINUE
         DO 2308 I=JJ,NP
            FACT=A(I,J)/A(J,J)
            DO 2309 K=1,NP
2309         ARX(I,K)=ARX(I,K)-ARX(J,K)*FACT
            DO 130 K=1,NC
130          A(I,K)=A(I,K)-A(J,K)*FACT
2308     CONTINUE
392   CONTINUE
!     ALL DONE WITH ROW REDUCTION
!     CHECK TO SEE HOW MANY ROWS ARE ALL ZERO--START FROM BOTTOM UP
      DO 2311 I=NP,1,-1
         DO 2311 J=1,NC
            IF(DABS(A(I,J)).GT.0.0001)GO TO 410
2311  CONTINUE
410   CONTINUE

!     	THE ITH ROW HAS NON ZERO ELEMENTS
!     	CALCULATE # INDEPENDENT REACTIONS=NRX
      	NRX = NP - I
!     	CALCULATE TOTAL # EQUATIONS=NEQ
!      	NEQ = NRX
!      	IF(IMASS.EQ.1)NEQ = NEQ + NC
!      	if(numNew.gt.0)NEQ = NEQ + numNew
!      	if(PFluidSwitch.eq.2)NEQ = NEQ + 1





!     	OUTPUT SECTION

      IF(iLong(7).eq.1)then
420      CONTINUE
         write(12,*)
         write(12,*)'NC=',NC,'  NP=',NP,'  # REACTIONS=',NRX
         write(12,*)
         write(12,*)'REDUCED COMPOSITION MATRIX'
         write(12,2100)(CONAME(I),I=1,NC)
 2100    FORMAT(' ',10X,30(4X,A4))
         DO 449 I=1,NP
         write(12,2101)VN3(I),(A(I,J),J=1,NC)
2101     FORMAT(' ',A8,2X,30(F8.3))
449      CONTINUE
         write(12,*)
         write(12,*)
         write(12,*)'REDUCED PHASE MATRIX'
         write(12,2110)(VN3(J),J=1,NP)
2110     FORMAT(' ',4X,50(A4,4X))
         DO 429 I=1,NP
         write(12,2120)(ARX(I,J),J=1,NP)
         IF(I.eq.NP-NRX)then
           write(12,*)' '
           write(12,*)' Linearly independent reactions:'
           write(12,2110)(VN3(J),J=1,NP)
!           write(12,*)' '
           endif
2120     FORMAT(' ',50F8.3)
429      CONTINUE
	endif      

      return
      end


! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE SysCoToPhCoTrans()
      implicit none
!     ROUTINE TO COMPUTE transformation matrix between the phase components and the system components
!	Note that it uses the same code from REXN but then adds the last bit to do the back substitution
!	to get the transformation matrix.

! c*****************************************
	include "Assemb.inc"
	include "Monit.inc"
	include "Output.inc"
! c*****************************************
! --------------------------------------------------------------------
	include "Tangent.inc"
! --------------------------------------------------------------------

!     Local arrays
!     DIMENSIONED FOR 12 SYSTEM COMPONENTS AND 25 PHASE COMPONENTS
      INTEGER*4 L,j,jj,k,i,pointerTemp,nrxn,kCur
!      REAL*8 A(50,26),TEMP,FACT
      Character VN3(maxPC)*4
!      common /REXNcommon/A,Temp,Fact,VN3
      real*8 TEMP,FACT
!      REAL*8 A(50,26),B(50,26)
!	common /CoTransCommon/a,b
	include "CoTransCommon.inc"

! c*****************************************

!	np is the number of phase components for this problem - equal to the number of rows (equations)
!	nc is the number of system components - equal to the number of columns

!     Set up array of phase component names
	L=0
	jj = 0
!	do 5002 k = 1,numPh
	Do 5002 kCur = 1,numPh
	k = asmCurrent(kCur)
	DO 5002 J=1,numPhCo(k)
	jj = jj + 1
        VN3(jj)=phCoName(k,J)
!	pause 'this code wont work - line 237 of Gibbs_rexn.f'
        pointer(jj) = jj			! the initial order of phase components
5002  CONTINUE





!     THIS Routine IS DESIGNED TO COMPUTE A LINEARLY INDEPENDENT
!     SET OF REACTIONS AMONG PHASE COMPONENTS IN A SYSTEM OF
!     NC SYSTEM COMPONENTS AND NP (np) PHASE COMPONENTS
!     USES METHOD OF GAUSS-JORDAN REDUCTION OPERATING ON AN IDENTITY
!     MATRIX (this part is Subroutine REXN)
!	Then the back substitution on the ncxnc part of the A matrix is done to get the identity matrix
!		at which point we have equations that give the system components in terms of the phase components

      IF (np.GT.maxPC) THEN
         WRITE(*,*) 'Dimension of array ARX in subroutine REXN is too small -----stop'
         WRITE(*,*) 'In SysCoToPhCoTrans routine'
         pause 'Hit return to continue'
         stop
      endif

!     Figure out the compositions of the minerals
 	Do 1431 kCur = 1,numPh
	k = asmCurrent(kCur)
      do 1432 i = 1,NC
      PhComp(K,i) = 0.
      do 1430 j=1,numPhCo(K)
      PhComp(K,i) = PhComp(K,i) + xPhCo(k,j)*comp(k,j,i)
1430   continue
1432   continue
1431   continue


!	set up A matrix with each row defining the composition of a phase component
!	in terms of the system components
!      DO 10 j=1,np
	jj = 0
!	do 10 k = 1,numPh
	Do 10 kCur = 1,numPh
	k = asmCurrent(kCur)
	do 10 j = 1,numPhCo(k)
	jj = jj + 1
        DO 10 L=1,NC
            A(jj,L)=COMP(k,j,L)
10    CONTINUE
!     ECHO INPUT DATA
      IF(tanOut(5).eq.1)then

483      CONTINUE
         write(12,*)
         write(12,*)
	 write(12,*)' SysCoToPhCoTrans routine'
         write(12,*)'    INPUT DATA'
         write(12,*)
         write(12,*)
         write(12,2000)(CONAME(I),I=1,NC)
2000     FORMAT(' ',10X,20(4X,A4))
         DO 12 I=1,np
            write(12,2001)VN3(I),(A(I,J),J=1,NC)
2001        FORMAT(' ',A8,2X,20(F8.3))
12       CONTINUE
	endif
	! Make B the identity matrix to start
      DO 50 J=1,np
         DO 50 I=1,np
50          B(I,J)=0.
      DO 51 J=1,np
51       B(J,J)=1.


!	reduce every column
	DO 392 J=1,NC
        JJ=J+1
	IF(DABS(A(J,J)).GT..00001)GO TO 2307		! check to see if the pivot is non-zero
	DO 2006 K=J+1,np
2006	IF(DABS(A(K,J)).GT.0.00001)GO TO 2106		! find a non-zero element in the jth column
! 	If here, then the assemblage is degenerate. Sometimes this does not cause a problem
! 	but sometimes it does.
! 	If we are doing a "PseudoFromScratch", then we want to avoid these degeneracies.
	If (ipseudo.eq.1)then
		idegen = 1
		return
		endif
         WRITE(*,*)' WARNING----'
         WRITE(*,*)' The number of system components specified is'
         WRITE(*,*)' greater than the number required to describe'
         WRITE(*,*)' the chemical variability of the phases'
         WRITE(*,*)' (the last column of array A in SUBROUTINE REXN'
         WRITE(*,*)' is all zero. COL=',J
         WRITE(*,*)' NOTE: usually, this does not cause a problem'
         write(*,*) 'Hit return to continue...'
         pause
         GO TO 392

!	found a non-zero pivot element in row k
2106     CONTINUE
!        SWITCH ROW K WITH ROW J
         DO 225 I=1,NC			! switch the entire row of nc elements
            TEMP=A(J,I)
            A(J,I)=A(K,I)
225         A(K,I)=TEMP
         DO 2306 I=1,np			! switch the entire row of np phase components
            TEMP=B(J,I)
            B(J,I)=B(K,I)
2306        B(K,I)=TEMP
	pointerTemp = pointer(j)
	pointer(j) = pointer(k)
	pointer(k) = pointerTemp
!	reduce the rest of the pivot row and all columns below
2307     CONTINUE
         DO 2308 I=JJ,np
            FACT=A(I,J)/A(J,J)
            DO 2309 K=1,np
2309         B(I,K)=B(I,K)-B(J,K)*FACT
            DO 130 K=1,NC
130          A(I,K)=A(I,K)-A(J,K)*FACT
2308     CONTINUE
392   CONTINUE
!     ALL DONE WITH ROW REDUCTION
!     CHECK TO SEE HOW MANY ROWS ARE ALL ZERO--START FROM BOTTOM UP
      DO 2311 I=np,1,-1
         DO 2311 J=1,NC
            IF(DABS(A(I,J)).GT.0.0001)GO TO 410
2311  CONTINUE
410   CONTINUE

!     	THE ITH ROW HAS NON ZERO ELEMENTS
!     	CALCULATE # INDEPENDENT REACTIONS=NRX
      	NRXn=np-I
!     	CALCULATE TOTAL # EQUATIONS=NEQ
!      	NEQ=NRX
!      	IF(IMASS.EQ.1)NEQ=NEQ+NC
!      	if(numNew.gt.0)NEQ = NEQ + numNew
!      	if(PFluidSwitch.eq.2)NEQ = NEQ+1


!	output the results thus far for debugging
	if(tanOut(5).eq.1)then
		write(12,*)
		write(12,*)' Half way reduced - A matrix '
        	write(12,2000)(CONAME(I),I=1,NC)
         	DO 447 I=1,np
            	write(12,2001)VN3(I),(A(I,J),J=1,NC)
447       	CONTINUE
		write(12,*)
		write(12,*)' Half way reduced - B matrix '
         	write(12,2110)(VN3(J),J=1,np)
!         	write(12,2110)(VN3(pointer(J)),J=1,np)
         	DO 448 I=1,np
         	write(12,2120)(B(I,J),J=1,np)
448		continue
		endif
! 	Now continue reduction in the backwards manner to generate the transformation matrix

!      DO 492 J=np (Np),2,-1	! back reduce every column - starting with the bottom last element
!	I think we only back substitute the ncxnc matrix, so I changed "Np (np)" to "nc" 2/7/2010
      DO 492 J=nc,2,-1	! back reduce every column - starting with the bottom last element
         JJ=J-1

! 	Check to see if the pivot element is greater than zero
         IF(DABS(A(J,J)).GT..00001)GO TO 407

! 	check for a column of all zeros
         DO 416 K=J+1,nc
416      IF(DABS(A(K,J)).GT.0.00001)GO TO 418
         WRITE(*,*)' Back substitution --- WARNING----'
         WRITE(*,*)' The number of system components specified is'
         WRITE(*,*)' greater than the number required to describe'
         WRITE(*,*)' the chemical variability of the phases'
         WRITE(*,*)' (the last column of array A in SUBROUTINE REXN'
         WRITE(*,*)' is all zero. COL=',J
         WRITE(*,*)' NOTE: usually, this does not cause a problem'
         write(*,*) 'Hit return to continue...'
	pause
	GO TO 492

418      CONTINUE
!        SWITCH ROW K WITH ROW J
         DO 425 I=1,nc
            TEMP=A(J,I)
            A(J,I)=A(K,I)
425         A(K,I)=TEMP
         DO 426 I=1,np
            TEMP=B(J,I)
            B(J,I)=B(K,I)
426         B(K,I)=TEMP
	pointerTemp = pointer(j)
	pointer(j) = pointer(k)
	pointer(k) = pointerTemp

! 	normal reduction
407     CONTINUE
	! working on the jjth column (jj = j-1)
         DO 408 I=JJ,1,-1
            FACT=A(I,J)/A(J,J)
            DO 409 K=1,np  
409         B(I,K)=B(I,K)-B(J,K)*FACT
            DO 430 K=1,nc
430          A(I,K)=A(I,K)-A(J,K)*FACT
408     CONTINUE


492   CONTINUE

! 	Finally, we turn the A matrix into the identity matrix by dividing each equation by A(j,j)

	do 510 i = 1,nc   ! loop on each row (leave the null space alone)
	Fact = A(i,i)
	A(i,i) = A(i,i)/fact

	do 515 j = 1,np     ! loop on each column in the non-null space
	B(i,j) = B(i,j)/fact
515	continue

510	continue

! 	We should have the transformation matrix (inverse) in B(Np=np x nc)
! 	and the identity matrix in A





!     	OUTPUT SECTION

      IF(tanOut(5).eq.1)then
420      CONTINUE
         write(12,*)
         write(12,*)'NC=',NC,'  np=',np,'  # REACTIONS=',NRXn
         write(12,*)
         write(12,*)'REDUCED COMPOSITION MATRIX - should be the identity matrix'
         write(12,2100)(CONAME(I),I=1,NC)
 2100    FORMAT(' ',10X,30(4X,A4))
         DO 449 I=1,np
!         write(12,2101)VN3(pointer(I)),(A(I,J),J=1,NC)
         write(12,2101)VN3(I),(A(I,J),J=1,NC)
2101     FORMAT(' ',A8,2X,30(F8.3))
449      CONTINUE
         write(12,*)
         write(12,*)
         write(12,*)'REDUCED PHASE MATRIX'
!         write(12,2110)(VN3(pointer(J)),J=1,np)
         write(12,2110)(VN3(J),J=1,np)
2110     FORMAT(' ',4X,30(A4,4X))
         DO 429 I=1,np
         write(12,2120)(B(I,J),J=1,np)
         IF(I.eq.np-NRXn)then
           write(12,*)' '
           write(12,*)' Linearly independent reactions:'
           write(12,2110)(VN3(J),J=1,np)
!           write(12,*)' '
           endif
2120     FORMAT(' ',30F8.3)
429      CONTINUE
	endif      

      return
      end


! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!      SUBROUTINE CoTrans(A,B,Na,Np)
      SUBROUTINE CoTrans(Na,Np)
      implicit none
!     ROUTINE TO calculate component transformation given the XtoSite matrix

! 	A is the XtoSite matrix (na x np)
! 	B is the SitetoX matrix (np x na) (starts out as na x na identity matrix)

! 	Na is number of atoms in a single phase
! 	Np is number of phase components in a single phase


!     Local arrays
! 	dimensioned for 15 siteatoms and 10 phase components
!       Na <= 15, Np <= 10

      INTEGER*4 Na,Np,i,j,jj,k,Nco
!	  REAL*8 A(12,10),B(12,12),TEMP,FACT      
      real*8 TEMP,FACT
!      REAL*8 A(50,26),B(50,26)
!	common /CoTransCommon/a,b
	include "CoTransCommon.inc"
!*****************************************

!     USES METHOD OF GAUSS-JORDAN REDUCTION OPERATING ON AN IDENTITY
!     MATRIX

!      IF (NP.GT.12) THEN
       IF (NP.GT.100) THEN
 	write(*,*)'Dimension of array B in subroutine CoTrans is too small -----stop'
! 	call fss_alert('Dimension of array B in subroutine REXN is too small -----stop')
	pause
         stop
      endif


! 	Set up identity matrix
      DO 10 J=1,Na
         DO 10 I=1,Na
10          B(I,J)=0.
      DO 11 J=1,Na
11       B(J,J)=1.


      DO 392 J=1,Np	!reduce every column
         JJ=J+1

! 	Check to see if the pivot element is greater than zero
         IF(DABS(A(J,J)).GT..00001)GO TO 307

! 	check for a column of all zeros
         DO 316 K=J+1,Na
316      IF(DABS(A(K,J)).GT.0.00001)GO TO 318
         WRITE(*,*)' WARNING----'
         WRITE(*,*)' The number of system components specified is'
         WRITE(*,*)' greater than the number required to describe'
         WRITE(*,*)' the chemical variability of the phases'
         WRITE(*,*)' (the last column of array A in SUBROUTINE REXN'
         WRITE(*,*)' is all zero. COL=',J
         WRITE(*,*)' NOTE: usually, this does not cause a problem'
         write(*,*) 'Hit return to continue...'
         read(*,*)
         GO TO 392

318      CONTINUE
!        SWITCH ROW K WITH ROW J
         DO 325 I=1,Np			! nP is the number of phae components in this mineral
            TEMP=A(J,I)
            A(J,I)=A(K,I)
325         A(K,I)=TEMP
         DO 326 I=1,Na			! nA is the number of site cations
            TEMP=B(J,I)
            B(J,I)=B(K,I)
326         B(K,I)=TEMP

! 	normal reduction
307     CONTINUE
         DO 308 I=JJ,Na
            FACT=A(I,J)/A(J,J)
            DO 309 K=1,Na
309         B(I,K)=B(I,K)-B(J,K)*FACT
            DO 330 K=1,Np
330          A(I,K)=A(I,K)-A(J,K)*FACT
308     CONTINUE


392   CONTINUE


!     ALL DONE WITH ROW REDUCTION
!     CHECK TO SEE HOWMANY ROWS ARE ALL ZERO--START FROM BOTTOM UP
      DO 2311 I=Na,1,-1
         DO 2311 J=1,Np
            IF(DABS(A(I,J)).GT.0.0001)GO TO 410
2311  CONTINUE
410   CONTINUE

!     THE ITH ROW HAS NON ZERO ELEMENTS
!     CALCULATE # INDEPENDENT components = Nco
      Nco=i

	if(Nco.ne.Np)then
	write(*,*)'Calculated number of phase components not equal to input number.  Subroutine CoTran.'
	write(*,*)' Nco (calc), np(input) = ',nco,np
! 	call fss_alert('Calculated number of phase components not equal to input number.  Subroutine CoTran.)
	pause
	stop
	endif


! 	Now continue reduction in the backwards manner to generate the transformation matrix

      DO 492 J=Np,2,-1	! back reduce every column - starting with the bottom last element
         JJ=J-1

! 	Check to see if the pivot element is greater than zero
         IF(DABS(A(J,J)).GT..00001)GO TO 407

! 	check for a column of all zeros
         DO 416 K=J+1,Na
416      IF(DABS(A(K,J)).GT.0.00001)GO TO 418
         WRITE(*,*)' WARNING----'
         WRITE(*,*)' The number of system components specified is'
         WRITE(*,*)' greater than the number required to describe'
         WRITE(*,*)' the chemical variability of the phases'
         WRITE(*,*)' (the last column of array A in SUBROUTINE REXN'
         WRITE(*,*)' is all zero. COL=',J
         WRITE(*,*)' NOTE: usually, this does not cause a problem'
         write(*,*) 'Hit return to continue...'
	pause
	GO TO 492

418      CONTINUE
!        SWITCH ROW K WITH ROW J
         DO 425 I=1,Np
            TEMP=A(J,I)
            A(J,I)=A(K,I)
425         A(K,I)=TEMP
         DO 426 I=1,Na
            TEMP=B(J,I)
            B(J,I)=B(K,I)
426         B(K,I)=TEMP

! 	normal reduction
407     CONTINUE
         DO 408 I=JJ,1,-1
            FACT=A(I,J)/A(J,J)
            DO 409 K=1,Na  
409         B(I,K)=B(I,K)-B(J,K)*FACT
            DO 430 K=1,Np
430          A(I,K)=A(I,K)-A(J,K)*FACT
408     CONTINUE


492   CONTINUE

! 	Finally, we turn the A matrix into the identity matrix by dividing each equation by A(j,j)

	do 510 i = 1,np   ! loop on each row
	Fact = A(i,i)
	A(i,i) = A(i,i)/fact

	do 515 j = 1,na     ! loop on each column
	B(i,j) = B(i,j)/fact
515	continue

510	continue

! 	We should have the transformation matrix (inverse) in B(Np x Na)
! 	and the identity matrix in A



      return
      end
