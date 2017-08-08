c This subroutine returns x times u,ubarr,d,dbarr,s,c,b,t,g
c distribution functions at a scale q2 for proton or antiproton
      subroutine fonstru0(x,ih,q2,xuh,xubh,xdh,xdbh,xsh,xch,xbh,xth,
     # xgpro)
      implicit none
      real * 8 pi
      parameter (pi=3.141592653589793d0)
      real * 8 x,q2,xuh,xubh,xdh,xdbh,xsh,xch,xbh,xth,xgpro
      integer ih
      integer ih1,ih2
      common/hadr/ih1,ih2
      real * 8 h1pdf,h2pdf
      common/mypdfsets/h1pdf(20),h2pdf(20)
      integer nfl
      common/flnumb/nfl
      real * 8 assc
      common/asscale/assc
      include 'xsect.h'
      real * 8 hc2,zrs,zal,cm,cmu,cmp
      common/choi/hc2,zrs,zal,cm,cmu,cmp
c local variables
      integer ndns,iht
      real * 8 xmb,xmc
      common/hvqmass/xmb,xmc
      real * 8 xq2,xx
      real * 8 fx(-5:5)
c
      if(ih.eq.ih1) then
         if(ih1.eq.100) then
            xuh = 0
            xubh= 0
            xdh = 0
            xdbh= 0
            xsh = 0
            xch = 0
            xbh = 0
            xth = 0
            xgpro= 0
            xch = xx/137.d0/2/pi*4.d0/9.d0*3*(xx**2 + (1-xx)**2)
     #           *log(xq2/xmc**2)
            xbh = 0
         endif
         ndns=nint(h1pdf(3))
      elseif(ih.eq.ih2) then
         ndns=nint(h2pdf(3))
      endif
c more common convention: 1 for proton, -1 for antiproton
      iht=ih
      xq2 = q2
      xx = x
      if(istrsc.eq.1) then
         xq2=xq2/cm**2
      endif
      call fonllmlmpdf(ndns,iht,xq2,xx,fx,5)
      xuh=x*fx(1)
      xubh=x*fx(-1)
      xdh=x*fx(2)
      xdbh=x*fx(-2)
      xsh=x*fx(3)
      xch=x*fx(4)
      xbh=x*fx(5)
      xth=0
      xgpro=x*fx(0)
      xch=xch/assc
c      xbh=bot/10.
c      xbh = fb(alfai(q2),qq,x)
c      xbh = fbm(x,qq)
      xbh=xbh/assc
      if(abs(nfl).eq.4) then
	xbh=0d0
      endif
c      write(*,'(3(e12.5,1x))') x, bot/10., xbh
      end
 
c     alfa strong at 2 loops at scale q2
      double precision function alfai(q2)
      implicit none
      real * 8 q2
      integer nfl
      common/flnumb/nfl
      real * 8 lam5qcd
      common/lamqcd/lam5qcd
      real * 8 alfas_p
      alfai = alfas_p(q2,lam5qcd,nfl)
      end

      subroutine setlambda(lam5qcd)
      implicit none
      integer j,ih,iret
      character * 2 sche
      real * 8 lam5qcd,lam4qcd,xlam
      integer ih1,ih2
      common/hadr/ih1,ih2
      real * 8 h1pdf,h2pdf
      common/mypdfsets/h1pdf(20),h2pdf(20)
      real * 8 masch2,masbo2,masto2,lambda2
      common/alfacm/masch2,masbo2,masto2,lambda2
      integer ioutput
      common/out/ioutput 
      if(lam5qcd.le.0) then
         if(ih1.eq.0.or.abs(ih1).eq.1) then
            j=nint(h1pdf(3))
         elseif(ih2.eq.0.or.abs(ih2).eq.1) then
            j=nint(h2pdf(3))
         else
            write(*,*) ' don''t know how to guess Lambda'
            stop
         endif
         ih=1
         call pdfpar(j,ih,xlam,sche,iret)
         lam5qcd=xlam
      endif
      lam4qcd = lam5qcd*(5./lam5qcd)**(2./25.)*
     #     (2*log(5./lam5qcd))**(963./14375.)
      lambda2 = lam4qcd**2
 155  format(f5.4,16x,'! lambda5 (GeV) with the set of hadron 1')
 156  format(f5.4,16x,'! lambda5 (GeV) NB. My input!!!')
      end

c      subroutine pftopdg(xd,qd,fxp)
c      real * 8 xd,qd,fxp(*)
c      write(*,*) ' not using pdflib!'
c      stop
c      end
