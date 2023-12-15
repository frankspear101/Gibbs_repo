!*********/SaveRxn/************************
!	Xsave and Ysave save the X-Y plotting positions
!	PenSave saves the penstate for each point (0=up; 1=down)

	integer*4 plotCounterMax
	Parameter(plotCounterMax = 4000)
      	integer*4 plotCounter	! Total points in this iteration
      	real*4 plotSave(plotCounterMax,4)
      	Integer*2 penSave(plotCounterMax)
      	COMMON /SaveRxn/plotSave,penSave,plotCounter
!*************************************
