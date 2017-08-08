c Weizsaecker-Williams function, for various parametrizations:
c To initialize, initialize the real * 8 common
c      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
c xmuww^2 = cut in Q2
c thcww   = angular cut
c elphen  = electron energy (in the frame where the angular cut is applied)
c zminww  = minimum photon energy fraction (with respect to incoming electron
c zmaxww  = maximum photon energy fraction
c itypeww:
c   itypeww=1 -> cut on Q2
c   itypeww=2 -> cut on angle
c   itypeww=3 -> cut on Q2, log term only
c   itypeww=4 -> cut on angle, log term only
c   itypeww=5,6 -> cut on Q2, used only by jetelgen; neglects O(m2/Q2) terms
c
c
c Initialization can also be performed by
c call fww_ww0(iinput,ioutput)
c The initialization subroutine prompts into unit ioutput
c for input from unit iinput.
c
c
c
      subroutine fww_ww0(iinput,ioutput)
c wrapper to an entry
c (it avoid direct call to an entry, because
c it is more convenient to use the subroutine calling
c form, and some compilers do not tolerate calling a function
c as a subroutine)
      integer iinput, ioutput
      real * 8 fww_ww00,tmp
      tmp=fww_ww00(iinput,ioutput)
      end

      function fww_ww(x)
      implicit none
      real * 8 fww_ww,fww_ww00,fww_ww1,x,xx
      real * 8 pi,alfaem
      parameter(pi=3.14159265358979312D0,alfaem=1/137.d0)
      real * 8 aemo2pi,xme,xme2
      parameter(aemo2pi=alfaem/(2*pi),xme=0.511d-3,xme2=xme**2)
      real * 8 xmuww,xmuww2,q2eff,tmp
      integer iinput,ioutput,icoeffs

      integer itypeww
      real * 8 elphen,thcww
c ww common
      real * 8 zminww,zmaxww
      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
      
      goto 1
      entry fww_ww1(icoeffs,x)
c return the coefficient of log(q2eff) (icoeffs=2), return the remaining 
c term (icoeff=1); neglects O(m2/Q2) terms
      if(icoeffs.eq.1) then
         itypeww=5
      elseif(icoeffs.eq.2) then
         itypeww=6
      endif
      goto 1
c
c iinput is not used in what follows
      entry fww_ww00(iinput,ioutput)
c      
      write(*,*)'In the Weizsaecker-Williams function:'
      write(*,*)'enter 1 for cut in photon virtuality'
      write(*,*)'      2 to cut on angle'
      write(*,*)'      3 as in 1, only log term'
      write(*,*)'      4 as in 2, only log term'
      read(55,*)itypeww
      write(ioutput,661) itypeww
 661  format
     #  (i1,12x,'! 1=Q2, 2=angle, 3=Q2, log only, 4=angle, log only') 
c Necessary in order to use the same fww_ww function used by jetelgen.f
      if(itypeww.eq.1.or.itypeww.eq.3)then
         write(*,*)'enter the effective WW scale in GeV (upper limit'
         write(*,*)'of the absolute value of the photon virtuality)'
         read(55,*)xmuww
         write(ioutput,662) xmuww
 662     format(d10.4,3x,'! xmu_WW') 
         xmuww2=xmuww**2
      elseif(itypeww.eq.2.or.itypeww.eq.4)then
         write(*,*)'enter electron energy, theta_cut'
         read(55,*)elphen,thcww
         write(ioutput,663) elphen,thcww
 663     format(d10.4,1x,d10.4,1x,'! elphen,thcww') 
      else
         write(*,*)'error: undefined itypeww in fww'
         stop
      endif
 2    continue
      return
c end initialization
c
 1    continue
      xmuww2=xmuww**2
      if(itypeww.eq.1)then
c q2 cut
         q2eff = xmuww2
         tmp = (1+(1-x)**2) / x * log( q2eff*(1-x)/(xme*x)**2 )
     #        +2*xme2*x*( 1/q2eff-(1-x)/(xme*x)**2 )
         fww_ww = aemo2pi * tmp
      elseif(itypeww.eq.2)then
c angular cut
         q2eff = xme2*x**2+elphen**2*(1-x)**2*thcww**2
         tmp = (1+(1-x)**2) / x * log( q2eff/(xme*x)**2 )
     #        +2*(1-x) * ( xme2*x/q2eff - 1/x )
         fww_ww = aemo2pi * tmp
      elseif(itypeww.eq.3)then
c as in 1 with log term only
         q2eff = xmuww2
         tmp = (1+(1-x)**2) / x * log( q2eff*(1-x)/(xme*x)**2 )
         fww_ww = aemo2pi * tmp
      elseif(itypeww.eq.4)then
c as in 2 with log term only
         q2eff = xme2*x**2+elphen**2*(1-x)**2*thcww**2
         tmp = (1+(1-x)**2) / x * log( q2eff*(1-x)/(xme*x)**2 )
         fww_ww = aemo2pi * tmp
      elseif(itypeww.eq.5)then
c everything but the coefficient of log(q2eff) and O(m2/Q2) terms
         tmp = aemo2pi*(1+(1-x)**2)/x*log( (1-x)/(x**2*xme2) )
     #        -2*(1-x)/x 
         fww_ww = aemo2pi * tmp
      elseif(itypeww.eq.6)then
c coefficient of log(q2eff)
         fww_ww = aemo2pi*(1+(1-x)**2)/x
      else
         write(6,*)'error in function fww_ww'
         stop
      endif
      return
      end      
