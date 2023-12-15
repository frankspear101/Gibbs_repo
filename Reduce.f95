! !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      SUBROUTINE REDUCE(NEQ,NVAR,DETA,IER)
	use MatrixArrays
      implicit none

!	include "Solute.inc"
	include "MatrixSingular.inc"

!     ROUTINE TO SOLVE A SYSTEM OF LINEAR EQUATIONS BY GAUSS-JORDAN REDUCTION
!     A IS THE INPUT MATRIX
!     NEQ IS THE DIMENSIONS (NEQxNEQ)
!     THE NEQ TO NVAR COLUMNS OF A CONTAIN THE DATA VECTOR(S)
!     XX(I,J) IS THE SOLUTION VECTOR(S); THE J COLS CONTAIN SOLUTIONS FOR EACH
!     OF THE DATA VECTORS CONTAINED IN THE NEQ-NVAR COLUMNS OF A
!     DETA IS THE DETERMINANT
!     IER IS A CONDITION VARIABLE =0 IF MATRIX IS NONSINGULAR, =1 IF MATRIX IS
!          SINGULAR
      integer NCOND,N2,i,M,L,K,Nvar,Neq,ier,j,ii
      real*8 B,DETA,zero,C,D
!      IMPLICIT REAL*8 (A-H,O-Z)
!
      DATA ZERO/1.0D-250/
!      DATA ZERO/1.0D-150/
!
	matrixIsSingular = 0
!	write(12,*)'  '
!	write(12,*)'  '
!	write(12,*)'New reduction  '
 
 
  100 DETA=1.0D0
      NCOND=NVAR-NEQ
      N2=NEQ-1
      DO 109 I=1,N2
       M=I
       B=A(I,I)
       DO 102 L=I,NEQ
          C=DABS (B)
          IF(C-DABS(A(L,I)))101,102,102
  101     M=L
          B=A(L,I)
  102  CONTINUE
       IF(DABS(B).LT.ZERO)then		! Row M is pivot. B is value of pivot. Check to see that it is non-zero
       	!	write(*,*)'Line 39  ',B
       		GO TO 114
       		endif
       IF(I-M)104,105,104
  104  DETA=-DETA
!        PIVOT IS LOCATED AS Mth ROW
!        SWITCH Mth AND Ith ROWS
  105  DO 106 L=I,NVAR
          B=A(I,L)
          A(I,L)=A(M,L)
          A(M,L)=B
  106  CONTINUE
       B=A(I,I)
!        AFTER REDUCING A ROW, DIVIDE ENTIRE ROW BY PIVOT TO MAKE PIVOT 1.0
       DO 107 L=I,NVAR
          A(I,L)=A(I,L)/B
  107  CONTINUE
       DETA=DETA*B
       K=I+1
       DO 108 M=K,NEQ
          D=A(M,I)
          DO 108 L=I,NVAR
             A(M,L)=A(M,L)-D*A(I,L)
  108  CONTINUE
!         	write(12,*)'--------------------'
!		DO 2314 ii=1,NEQ
!		write(12,2320)(A(ii,J),J=1,nvar)
!2314  		continue
  109 CONTINUE
      DETA=DETA*A(NEQ,NEQ)

!     BACK SUBSTITUTE TO GET SOLUTION
      IF(DABS(A(NEQ,NEQ)).LT.ZERO)then
      	!	write(*,*)'Line 69  ',NEQ,A(NEQ,NEQ)
      		GO TO 114
      		endif
      DO 150 I=1,NVAR-NEQ
  150  XX(NEQ,I)=A(NEQ,I+NEQ)/A(NEQ,NEQ)
!     WORK ON ROW K
      K=NEQ-1
!      if(K) 114,113,111
	if(K.lt.0)then
		write(*,*)'line 77  K=  ',K
		go to 114
		endif
	if(K.eq.0)then
		return
		endif

!     FIND THE SOLUTION FOR THE Ith DATA VECTOR
  111 M=K+1
      DO 160 I=1,NVAR-NEQ
       B=0.0D0
       DO 112 L=M,NEQ
          B=B+A(K,L)*XX(L,I)
  112  CONTINUE
       XX(K,I)=A(K,I+NEQ)-B
160    CONTINUE
      K=K-1
	if(K.lt.0)then
	!	write(*,*)'line 95  K=  ',K
		go to 114
		endif
	if(K.eq.0)then

		return
		endif
	go to 111


!113   RETURN

114   CONTINUE
!     IF HERE, THEN MATRIX IS SINGULAR
      IER=1
      DETA=0.0D0
	if(iCareifMIS.eq.1)then
		WRITE(*,*)'A(NEQ,NEQ) IS 0) - matrix is singular'
		write(*,*)zero,A(NEQ,NEQ)
         	write(12,*)'--------------------'
		DO 2315 ii=1,NEQ
		write(12,2320)(A(ii,J),J=1,nvar)
2315  		continue
2320  		FORMAT(' ',30E12.4)
		endif!      IER=1
	matrixIsSingular = 1
      RETURN
      END

