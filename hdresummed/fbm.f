      subroutine fonstru(x,ih,q2,xuh,xubh,xdh,xdbh,xsh,xch,xbh,xth,
     # xgpro)
      implicit none
      include 'xsect.h'
      character*1 hvqs
      common/hvqtype/hvqs
      real * 8 xmb,xmc
      common/hvqmass/xmb,xmc
      real * 8 x,q2,xuh,xubh,xdh,xdbh,xsh,xch,xbh,xth,xgpro
      integer ih
      real * 8 qrkm
      real * 8 fbm
      call fonstru0(x,ih,q2,xuh,xubh,xdh,xdbh,xsh,xch,xbh,xth,xgpro)
      if(hvqs.eq.'b') then
         qrkm = xmb
      elseif(hvqs.eq.'c') then
         qrkm = xmc
      else
         write(*,*) ' fbm not setup to do quark ', hvqs
         stop
      endif
      if(iwhichsfh.ne.0.or.sqrt(q2).lt.qrkm) then
         if(hvqs.eq.'b') then
            xbh=fbm(x,q2,qrkm)
         elseif(hvqs.eq.'c') then
            xch=fbm(x,q2,qrkm)
         else
            write(*,*) ' fbm not setup to do quark ', hvqs
            stop
         endif
      endif
      end

c....calculates and returns x*F_h(x,q2) for a heavy quark of mass qrkm
      function fbm(x,q2,qrkm)
      implicit none
      real * 8 fbm,x,q2,qrkm
c parameters
      integer npoints
      parameter (npoints=200)
c static local
      real * 8 zarr(npoints),farr(npoints)
      real * 8 oldq2,oldqrkm,step,xmin,xmax
      integer na(1)
c local
      real * 8 ztmp
c functions
      real * 8 dfint
      data na/npoints/
      data oldq2,oldqrkm/-1.d0,-1.d0/
      data step/0.2d0/
      save oldq2,oldqrkm,step,xmin,xmax,zarr,farr
      if(q2.ne.oldq2.or.qrkm.ne.oldqrkm) then
         call setupfbm(q2,qrkm,zarr,farr,step,npoints)
         oldq2=q2
         oldqrkm=qrkm
         xmax=exp( npoints/2*step)/(exp( npoints/2*step)+1)
         xmin=exp(-npoints/2*step)/(exp(-npoints/2*step)+1)
      endif
      if(x.lt.xmin.or.x.gt.xmax)
     #    write(*,*) ' warning: x out of range in fbm ,=',x
      ztmp=log(x/(1-x))
      fbm=dfint(1,ztmp,na,zarr,farr)
      end

c....calculates and returns x*F_h(x,q2) for a heavy quark of mass qrkm
      subroutine setupfbm(q2,qrkm,zarr,farr,step,npoints)
      implicit none
      integer npoints
      real * 8 q2,qrkm,zarr(npoints),farr(npoints),step     
      real * 8 lam5qcd
      common/lamqcd/lam5qcd
      real * 8 xx,sc,asav,xm
      common/fbmine/xx,sc,asav,xm
      integer nfl
      common/flnumb/nfl
c static local
      real * 8 pi
c local
      integer j
      real * 8 t,asxm,as
c functions
      real * 8 alfas_p,dgauss,fbmconv
      external fbmconv
      data pi /3.1415927d0/
      save pi
      sc=sqrt(q2)
      xm = qrkm
      asxm = alfas_p(xm**2,lam5qcd,nfl)
      do j=1,npoints
         zarr(j) = (j-npoints/2)*step
         xx = exp(zarr(j))
         xx = xx/(xx + 1)
         as = alfas_p(q2,lam5qcd,nfl)
         asav=(as+asxm)/2
         t=asav/(2*pi)*log(q2/xm**2)
c......"partly resummed", with running alphas
c      fbm = x * t *         
         farr(j) = xx * t * dgauss(fbmconv,xx,1.d0,0.001d0)
c dgquad(fbmconv,xx,1.d0,4)
      enddo
      end
      
      function fbmconv(z)
      implicit real*8(a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      common/fbmine/xx,sc,asav,xm
      real * 8 x,q2,xuh,xubh,xdh,xdbh,xsh,xch,xbh,xth,xgpro
      integer ih
      include 'xsect.h'
      y = xx/z
      x=xx
      q2=sc**2
      call fonstru0(y,ih,q2,xuh,xubh,xdh,xdbh,xsh,xch,xbh,xth,
     # xgpro)
      gg = xgpro/y
      fbmconv = 0.5d0*(z**2 + (1-z)**2)*gg/z
c Add NL bit, if required (iwhichsfh=2) or if we are simply supplying
c for the lack of mu<m values in the pdf. In fact, the last step
c should only hold for NL pdf's, to be fixed in the future.
      if(iwhichsfh.eq.2.or.(iwhichsfh.eq.0.and.sc.lt.xm)) then
         fbmconv=fbmconv+asav/(2*pi)*p1sgf(z)*gg/z
      endif
      end

      function p1sgf(x)
      implicit none
      real * 8 p1sgf,x
      real * 8 tf,cf,cg,nf,omx,lgomx,lgx,pgfx,pgfmx,s2,pi
      parameter (pi=3.14159265358979312D0)
      real * 8 ddilog
      tf = 1.0d0/2.0d0
      cf = 4.0d0/3.0d0
      cg = 3
      nf = 1
      omx = 1-x
      lgomx = log(omx)
      lgx = log(x)
      pgfx = x**2+omx**2
      pgfmx = (x+1)**2+x**2
      s2 = (log(x/(x+1))**2-log(x+1)**2)/2.0d0-log(x+1)*log(x/(x+1))+2*d
     1   dilog(x/(x+1))-pi**2/6.0d0
      p1sgf = cg*nf*tf*(lgx*(1.36d2*x/3.0d0+(-3.8d1)/3.0d0)+1.4d1*x/9.0d
     1   0+lgx**2*(-8*x-2)+4.0d1/(9.0d0*x)+2*pgfmx*s2+(-lgx**2+4.4d1*lgx
     2   /3.0d0-2*lgomx**2+4*lgomx+pi**2/3.0d0+(-2.18d2)/9.0d0)*pgfx-4*
     3   lgomx+1.82d2/9.0d0)+cf*nf*tf*(lgx*(4*x-1)+lgx**2*(2*x-1)-9*x+(2
     4   *lgx**2-4*lgomx*lgx+4*lgx+2*lgomx**2-4*lgomx+(-2.0d0)*pi**2/3.
     5   0d0+10)*pgfx+4*lgomx+4)
      end
