c....calculates and returns x*F_b(x,qq) or x*F_c(x,qq)

      function phfbm(x,qq,xmh)
      implicit real*8(a-h,o-z)
      integer ifixed
      character*1 hvqs
      common/fbmine/xx,sc
      COMMON/W50512/QCDL4,QCDL5
c      common/hvqmass/xmb,xmc
      common/asscale/assc
      common/hvqtype/hvqs
      common/fixedorder/ifixed
      external convfbm
      data pi /3.1415927d0/, nf /5/

      if(qq.le.xmh) then
           phfbm = 0.
           return
      endif
    
      if(hvqs.eq.'c') then
           nf = 4
      endif

      xm = xmh
      xx = x
      sc = qq
            
cc      as = ALPHAS2(sc)/10.
cc      asxm = alphas2(xm)/10.
      as = alfas_p(sc**2,QCDL5,nf)
      asxm = alfas_p(xm**2,QCDL5,nf)

      b0 = (33.-2.*nf)/12./pi
c definition of t rescaled (only works at order as)
      t = 0.5/pi/b0*log(asxm/as) /assc
      
      if(ifixed.eq.1) then
c.....with fixed alphas       
        phfbm = x * as/2./pi*log(sc**2/xm**2)*
     #             dgauss(convfbm,x,1.d0,0.001d0)
      elseif(ifixed.eq.0) then
c......"partly resummed", with running alphas
         phfbm = x * t * dgauss(convfbm,x,1.d0,0.001d0)         
      endif
       
      end
      
            
      function convfbm(z)
      implicit real*8(a-h,o-z)
      real*8 xpdf(-6:6)
      common/fbmine/xx,sc
      write(*,*) ' fbm-gamma: option no longer maintained'
      stop
c      y = xx/z
c      call pftopdg(y,sc,xpdf)
c      gg = xpdf(0)/y
c      convfbm = 0.5*(z**2 + (1.-z)**2)*gg/z
      
      end
