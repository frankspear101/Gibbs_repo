!*******/OUTPUT/*****************************************************
      INTEGER iLong(20),masserroralert,newtonStepsOutput
      COMMON /OUTPUT/iLong,masserroralert,newtonStepsOutput
!********************************************************************
!	ilong(1) = error output
!	ilong(2) = Jacobian matrix
!	ilong(3) = Thermodynamic data for every mineral
!	ilong(4) = Running sum at each step
!	ilong(5) = Transformation matrix
!	ilong(6) = Activity calculations
!	ilong(7) = REXN" independent reaction output on each step
!	ilong(8) = Print phase compositions in Wt% on (7)
!	ilong(9) = Write minerals to .asm file on (7)
!	ilong(10) = 
!	ilong(11) = echo information in reading file input 

