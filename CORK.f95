!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine hprk(pbars,t,xco2,fh2o,fco2)
!!-----------------------------------------------------------------------
!! hprk routine to compute h2o-co2 fugacities from the holland and
!! powell (cmp 1991)
!!                                     j. connolly 1992.
 
!!  for the unitiated input and output is done through common blocks
!!  cst5, cst11, and cst26, the relevant variables are:
 
!!  for input:
!!               pbars= pressure in bars
!!               t    = temperature in kelvins
!!               r    = gas constant
!!               xco2 = mole fraction of co2 in fluid phase
!!  for output:
!!               v    = molar volume (cm3/mole) at p and t
!!               fco2 = the natural log of the co2 fugacity
!!-----------------------------------------------------------------------
      implicit double precision (a-g,o-y),integer (h-m,z)
 
!!      common/ cst5  /pbars,t,xco2,u1,u2,tr,pr,r,ps
!!     *      / cst11 /fh2o,fco2
 
	data r /8.31441d0/


      rt = r*t/1000.d0
 
      p = pbars/1000.d0
 
      if (xco2.eq.1d0) then
 
         call crkco2 (pbars,t,fco2) 

      else if (xco2.eq.0.d0) then
 
         call crkh2o (pbars,t,vol,fh2o) 

      else

         call crkco2 (pbars,t,fco2) 
         call crkh2o (pbars,t,vol,fh2o) 
  
         xh2o = 1.d0 - xco2
 
         wco2 = (13.2 - .290 * t**0.5d0) * p**0.25d0
         wh2o = (7.0  - 0.15 * t**0.5d0) * p**0.25d0
 
         gco2 = xh2o*xh2o*(wco2+2.d0*xco2*(wh2o-wco2))/rt
         gh2o = xco2*xco2*(wh2o+2.d0*xh2o*(wco2-wh2o))/rt
 
         fco2 = fco2 + gco2 + dlog(xco2)
         fh2o = fh2o + gh2o + dlog(xh2o)
 
      end if
      end



!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine crkh2o (pbar,t,vcrk,fh2o)
!!-----------------------------------------------------------------------
!! compute ln(f[h2o], bar) and volume h2o (kj/bar) from CORK EoS Holland 
!! & Powell CMP 109:265-273. Input pbar - pressure (bars); tk - temp (K).
!!                                     J.A.D. Connolly, 1992.
!! Modified by F. Spear July, 98 to fix errors and update to HoPo(98) virial coeffs
!!    This code reproduces HoPo(98) tables to 0.003.  
!!    Not sure why the discrepancy occurs - could be roundoff if HoPo(98) use single precision
!!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 

      dimension x(3) 

      data b,r,p0 /1.465d0,8.31441d-3,2.d0/

      p = pbar/1.d3
      rt = r*t
      rtp = rt/p
      t12 = dsqrt(t)
      t695 = 695.0d0
      t673 = 673.0d0
      t673_t = t673 - t

!!      if (t.lt.t695) then
!!         psat = -13.627d-3 + 7.29395d-7*t**2 - 2.34622d-9*t**3 + 4.83607d-15*t**5
!!         if (p.lt.psat.and.t.lt.t673) then
!!		write(*,*)'gas'
!!	    ! P < Psat and T < 673  -- gas use a0 a7 a8 a9
!!            a = 1113.4d0 + 5.8487d0*t673_t - 2.1370d-2*t673_t**2 + 6.8133d-5*t673_t**3
!!         else
!!            if (t.lt.t673) then 
!!	       ! P > Psat and T < 673  - liquid - use a0 a1 a2 a3
!!               a = 1113.4d0 - 0.88517*t673_t + 4.53d-3*t673_t**2 - 1.3183d-5*t673_t**3
!!            else
!!		! P = anything and 673 < T < 695  - a gas or liquid, we use liquid equation - use a0 a4 a5 a6
!!               a = 1113.4d0 - .22291d0*(-t673_t) - 3.8022d-4*(-t673_t)**2 + 1.7791d-7*(-t673_t)**3
!!            end if
!!         end if
!!      else 
!!	 ! T > 695  - above critical T - use liquid
!!         a = 1113.4d0 - .22291d0*(-t673_t) - 3.8022d-4*(-t673_t)**2 + 1.7791d-7*(-t673_t)**3
!!      end if


      if (t.lt.t695) then
         psat = -13.627d-3 + 7.29395d-7*t**2 - 2.34622d-9*t**3 + 4.83607d-15*t**5
      endif
      
      if (t.lt.t673) then
 	 if(p.lt.psat)then
	    !write(*,*)'gas'
	    ! P < Psat and T < 673  -- gas use a0 a7 a8 a9
            a = 1113.4d0 + 5.8487d0*t673_t - 2.1370d-2*t673_t**2 + 6.8133d-5*t673_t**3
          else
	    !liquid
 	    ! P > Psat and T < 673  - liquid - use a0 a1 a2 a3
            a = 1113.4d0 - 0.88517*t673_t + 4.53d-3*t673_t**2 - 1.3183d-5*t673_t**3
          endif
       else
        ! T > 673  -  use liquid equation for all
	! use a0 a4 a5 a6
        a = 1113.4d0 - .22291d0*(-t673_t) - 3.8022d-4*(-t673_t)**2 + 1.7791d-7*(-t673_t)**3
       end if



      a1 = -rtp
      a2 = a/t12/p - b*(rtp+b)
      a3 = -a*b/t12/p

      call roots3 (a1,a2,a3,x,xmin,xmax,iroots) 

      if (iroots.eq.1) then
         vol = x(1)
      else 
         if (p.lt.psat) then
            vol = xmax
         else 
            vol = xmin
         end if
      end if 



      cc = a/b/rt/t12

      gam = vol/rtp -1.d0 - dlog((vol-b)/rtp) - cc*dlog(1.d0+b/vol)

       vcrk = vol

!!	virial correction above 2 kbar
      if (p.gt.p0) then
!!         ! = -3.02565d-2 - 5.343144d-6*t            HP(91)
!!         d = -3.2297554d-3 + 2.2215221d-6*t         HP(91)
	 av = 1.9853d-3
	 bv = -8.9090d-2
	 cv = 8.0331d-2
         dp = p - p0
!!         vcrk = vol + !*dsqrt(dp) + d*dp            HP(91)
	 Vvir = av*dp + bv*dsqrt(dp) + cv*dp**0.25
	 vcrk = vol + Vvir
!!         gam = gam +  (dp**1.5d0*2.d0*!/3.d0 + dp**2*d/2.d0)/rt     HP(91)
         gam = gam +  ((av/2.0d0)*dp*dp + bv*(2.0d0/3.0d0)*dp**1.5d0 + cv*(4.0d0/5.0d0)*dp**1.25d0)/rt
      end if 

!!	Fugacity integration below Tc
      if (t.lt.695.d0.and.p.gt.psat) then
!!	error in logic here.  The integration must be done to Psat using liquid and gas equations
!!	calculated a P = Psat.  The original code overestimated lnf by approx 3 log units at TC=200, Pbar = 500
         p = psat
	 rtpsat = rt/psat
!!         a1 = -rtp
         a1 = -rtpsat
         a2 = -(b*(rt+b*p)-a/t12)/p
         a3 = -(a*b/t12)/p

         call roots3 (a1,a2,a3,x,vmin,vmax,iroots)
         if (iroots.eq.1) go to 99	!there should always be 3 roots at the saturation pressure
!!	 note that there is never any virial correction at Psat because Psat << P0 (=2 kbar)

	!integration for liquid - use liquid a value from above
         gam2 = vmin/rtpsat-1.d0-dlog((vmin-b)/rtpsat)-cc*dlog(1.d0+b/vmin)

	!integration for gas - use gas a value if T<673, otherwise use liquid value
	if(t.lt.t673)then
          a = 1113.4d0 + 5.8487d0*t673_t - 2.1370d-2*t673_t**2 + 6.8133d-5*t673_t**3
	endif
         cc = a/b/rt/t12
         gam1 = vmax/rtpsat-1.d0-dlog((vmax-b)/rtpsat)-cc*dlog(1.d0+b/vmax)

         gam = gam1 - gam2 + gam

      end if 

99    continue
	fh2o = gam + dlog(pbar)
	return
      end
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine crkco2 (pbar,t,fco2)
!!-----------------------------------------------------------------------
!! compute ln(f[co2], bar) from Eq 8 of Holland 
!! & Powell CMP 109:265-273. Input pbar - pressure (bars); tk - temp (K).
!!                                     J.A.D. Connolly, 1992.
!!-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 

      data b,r,tc,pc /3.7852d0,8.31441d-3,304.2d0,0.0738d0/

      p = pbar/1.d3
      rt = r*t
      bp = b*p

      a = (5.45963d-5*tc - 8.6392d-6*t)*tc**(1.5d0)/pc
      ! = (-3.30558d-5*tc + 2.30524d-6*t)/pc**(1.5d0)
      d = (6.93054d-7*tc - 8.38293d-8*t)/pc/pc

      fco2 = dlog(pbar) + (bp + a/b/dsqrt(t)*(dlog(rt+bp)&
            - dlog(rt+2.d0*bp)) + 2.d0 & !*p**(1.5d0)/3.d0
            + (d/2.d0)*p*p)/rt 


!!	virial correction above 5 kbar - makes lnfco2 worse than without for 
!!	probably should only use if MRK is used
!!      if (p.gt.5.0d0) then
!!	 av = 5.40776d-3 - 1.59046d-6*t
!!	 bv = -1.78198d-1 + 2.45317d-5*t
!!	 cv = 0.0d0
!!         dp = p - 5.0d0
!!	 Vvir = av*dp + bv*dsqrt(dp) + cv*dp**0.25
!!	 vcrk = vol + Vvir
!!         fco2 = fco2 +  ((av/2.0d0)*dp*dp + bv*(2.0d0/3.0d0)*dp**1.5d0 + cv*(4.0d0/5.0d0)*dp**1.25d0)/rt
!!      end if 



      end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      subroutine roots3 (a1,a2,a3,x,vmin,vmax,iroots)
!---------------------------------------------------------
! returns real roots (in x) of: x**3 + a1*x**2 + a2*x + a3
!---------------------------------------------------------

      implicit double precision (a-h,o-z)  

      dimension x(3)

      qq = (a1**2-3.d0*a2)/9.d0
      rr = (2.d0*a1**3-9.d0*a1*a2+27.d0*a3)/54.d0

      a5  = a1/3.d0

      dif = qq**3-rr**2

      if (dif.ge.0.d0) then
         phi = dacos(rr/dsqrt(qq**3))
         a4  = -2.d0*dsqrt(qq)
         a6  = phi/3.d0
         dphi = 0.d0
         vmin = 1.d9
         vmax = -1.d9
         do 10 i = 1, 3
            v = a4*dcos(a6+dphi) - a5
            if (v.gt.vmax) vmax = v
            if (v.lt.vmin) vmin = v
            x(i) = v
10          dphi = dphi + 2.094395102497915d0
         iroots = 3
      else
         a7 = (dsqrt(-dif) + dabs(rr))**(1.d0/3.d0)
         x(1) = -rr/dabs(rr)*(a7 + qq/a7) - a5
         iroots = 1
      end if 

      end

