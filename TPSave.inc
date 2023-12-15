!*********/TPSave/************************
!	Xsave and Ysave save the X-Y plotting positions
!	PenSave saves the penstate for each point (0=up; 1=down)
      integer*4 TPcounter	!Total points in this iteration
      real*4 TPSave(1000,2)
      COMMON /TPSavestuff/TPSave,TPcounter
!*************************************
