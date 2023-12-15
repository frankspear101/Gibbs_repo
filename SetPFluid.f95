!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      Subroutine SetPFluid

!     Program to set P Fluid
!	variable PFluidSwitch
!		PFluidSwitch = 0	Lithostatic (the default)
!		PfluidSwitch = 1	PFluid is independent (adds to variance)
!		PFluidSwitch = 2	PFluid is hydrostatic - this adds an additional equation and variable
!	Routine also alows user to set rock density (RockDensity)
!**************************************************************
      	implicit none

!***************************************************************
!*****************************************
	include "Assemb.inc"
!*****************************************
!      local variables
      	integer*4 iok,iChoice,izero
	real*8 PFluidTemp

500   continue

	write(*,*)' Setup Pfluid routine'
	write(*,*)'PfluidSwitch = ',PfluidSwitch
	write(*,*)'RockDensity  = ',RockDensity	
	write(*,*)'TC           = ',TC
	write(*,*)'TK           = ',TK
      	write(*,*)'PB           = ',PB
!	Find initial Pfluid.  This is done:
!	Lithostatic - Pfluid = Prock
!	Independent - Whatever user has set
!	Hydrostatic - Pf = Prock*fluiddensity/rockdensity
!		but we don't know fluid density unless we know Pfluid - so we must iterate
!		We do this through a call to subroutine fluiddensity

!	Set existing P Fluid - this can be user changed if PFluidSwitch = 2
	FluidDensity = 1
	select case(PFluidSwitch)	
		case(0)
		PFluidTemp = PB
		FluidDensity = 1
		case(1)
		PFluidTemp = PFluid
		FluidDensity = 1
		case(2)
!		Set existing P density (this is calculated at P and T from equation of state
		! note that this sub resets Pfluid
		call GetFluidDensity(izero)
		if(izero.eq.1)then
			call fss_alert('ALERT!!','P or T is out of bounds for fluid equation of state. Setting fluid density = 1')
			FluidDensity = 1
			endif
		PFluidTemp = PB * FluidDensity/RockDensity
		case default
		end select
	write(*,*)'FluidDensity = ',FluidDensity
	write(*,*)'PFluid       = ',PFluid
	write(*,*)'PFluidTemp   = ',PFluidTemp

10	continue
	write(*,*)' '
	write(*,*)' Choose your option for PfluidSwitch'
!		PFluidSwitch = 0	Lithostatic (the default)
!		PfluidSwitch = 1	PFluid is independent (adds to variance)
!		PFluidSwitch = 2	PFluid is hydrostatic - this adds an additional equation and variable
	write(*,*)'-1 = exit routine (no changes)'
	write(*,*)' 0 = Lithostatic pressure (the default)'
	write(*,*)' 1 = PFluid is independent (adds to variance)'
	write(*,*)' 2 = PFluid is hydrostatic (adds an additional variable and equation)'
	write(*,*)' Your choice?'
	read(*,*)ichoice
	

	select case (ichoice)
	case(-1)  	!exit
		iok=1
		return
	case(0)
		PFluidSwitch = 0
		PFluidTemp = PB
		PFluid = PB
		FluidDensity = 1
	case(1)
		write(*,*)'Input initial value for Pfluid'
		read(*,*)Pfluid
		PFluidTemp = PFluid
		FluidDensity = 1
		PFluidSwitch = 1

	case(2)
		PFluidSwitch = 2
		! note that this sub resets Pfluid
		call GetFluidDensity(izero)
		if(izero.eq.1)then
			call fss_alert('ALERT!!','P or T is out of bounds for fluid equation of state.  Pick another fluid option')
			go to 500
			endif
		PFluidTemp = PB * FluidDensity/RockDensity
		PFluid = PFluidTemp

	case default
		write(*,*)'bad choice....try again'
		go to 10
	end select
	go to 500

      end
