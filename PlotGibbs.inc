!***********/PlotGibbs/***********************************************
!
!	TYPE(AWE_Canvas) :: XYPlot		! The plotting canvas
!	X-Y plotting variables for Gibbs only
	Integer*4 nopen,ipen,irot,ibig,ispace,isymb,itrip,iAlphaBeta
!	REAL*4 xor,yor,xmin,ymin,xmax,ymax,xlen,ylen,width,height
	real*4 yMinBound,yMaxBound,xMinBound,xMaxBound	! limits on X-Y contouring (user adjustable)
  !    CHARACTER*80 XLAB,YLAB,PLTITL,Pltitle
!	Ternary plotting variables
	INTEGER*4 itic
	REAL*4 xorT,yorT,scale
	Character*40 Llab,Rlab,Tlab
	integer*4 NoCompPlotVars,Apx(4)
	real*4 CompPlotCoeff(20,24),CompPlotConst(20)
	character*8 CompPlotName(20)
	character*20 plotIndexName(5,4)
!	general plot variables
	integer*4 iternary,NoTieMinerals,PlotIndex(5,4),axisIsPhase(4,2)
	
      common /PlotGibbs/nopen,ipen,irot,ibig,ispace,isymb,&
     &   scale,iternary,axisIsPhase,&
     &   xorT,yorT,itic,Llab,Rlab,Tlab,&
     &   CompPlotCoeff,CompPlotConst,CompPlotName,&
     &   NoCompPlotVars,Apx,NoTieMinerals,PlotIndex,iTrip,iAlphaBeta,&
     &   yMinBound,yMaxBound,xMinBound,xMaxBound,plotIndexName
!*******************************************************************
