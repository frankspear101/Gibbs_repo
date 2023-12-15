!***************************************************************
!     ListManager requires this common block
	integer*4 listDim
	parameter (listDim = 1000)
      character*80 theList(listDim)
      Character*32 whatToPick
      integer*4 lWindowSize,theOne,listSize
      common/liststuff/theList,lWindowSize,listSize,theOne,whatToPick
!***************************************************************
