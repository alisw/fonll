      subroutine fonllpdfset(parm,hpdf)
      implicit none
      character * 2 sche
      real * 8 xlam
      real * 8 qcdl4,qcdl5
      integer ihphres
      common/hadr/ihphres
      common/w50512/qcdl4,qcdl5
      character * 20 parm(20)
      real * 8 hpdf(20)
      integer ih,ndns
      common /pdflib2mlm/ih,ndns
      integer iret
      ih=ihphres
      ndns=hpdf(3)+hpdf(2)*100
      call pdfpar(ndns,ih,xlam,sche,iret)
      if(iret.ne.0) stop
      qcdl5=xlam
c In photon_wps8 the result does not depend upon qcdl4;
c however, qcdl4 cannot be set to 0
      qcdl4=qcdl5
      end



      subroutine fonllstructm(x,q,upv,dnv,usea,dsea,str,
     #             chm,bot,top,gl)
      implicit none
      real * 8 x,q,upv,dnv,usea,dsea,str,
     #             chm,bot,top,gl
      integer ih,ndns
      common /pdflib2mlm/ih,ndns
      real * 8 fx(-5:5),q2,xs
      q2=q**2
      xs=x
      call fonllmlmpdf(ndns,ih,q2,xs,fx,5)
      usea=x*fx(-1)
      upv=x*fx(1)-usea
      dsea=x*fx(-2)
      dnv=x*fx(2)-dsea
      str=x*fx(3)
      chm=x*fx(4)
      bot=x*fx(5)
      top=0
      gl=x*fx(0)
      end

      subroutine fonllmlmpdf(ndns,ih,q2,xs,fx,nf)
      implicit none
      integer ndns,ih,nf
      real * 8 q2,xs,fx(-nf:nf),xxs
      real * 8 xpdf(-6:6)
      real * 8 q
      integer k
      integer atomicnumber,pset,AA
      common/shadowing/atomicnumber
      real * 8 A,ruv,rdv,ru,rd,rs,rc,rb,rt,rg
      q=sqrt(q2)
      if(ih.eq.5) then
         call ppgrv_int(xs,q,xpdf,nf,ndns)
         do k=-nf,nf
            fx(k)=xpdf(k)/xs
         enddo
      else
c.....include nuclear effects
      if(ih.eq.0.and.atomicnumber.gt.2) then
         call mlmpdf8(ndns,1,q2,xs,fx,nf)
	 if (atomicnumber .lt. 1000) then
c........fonll crashes when used with eks98 in some situations,
c        apparently because xs gets too cloose to 1 and an index
c        in eks98 goes out of boundary. This should fix it.
c        MC, 8/5/2007
              xxs = xs
              if ( xs .ge. 1d0-1d-8) then
	          xxs = xs - 1d-8
	      endif
	      A = atomicnumber
              call eks98(xxs,q,A,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
	 else
c... this calls the new EPS09. Atomic number and pset are to be extracted
c... from an input of the form AA*1000+pset
c... (note that A was a real in eks98, but an integer in eps09, hence AA)
              AA = atomicnumber/1000
              pset = atomicnumber - 1000*AA              
              call eps09(2,pset,AA,xs,q,ruv,rdv,ru,rd,rs,rc,rb,rg)
	 endif
	     	     
	 fx(0) = fx(0)*rg                           ! gluon
	 fx(1) = (fx(1) - fx(-1))*ruv + fx(-1)*ru   ! uval + usea
         fx(2) = (fx(2) - fx(-2))*rdv + fx(-2)*rd   ! dval + dsea
	 fx(-2) = fx(-2)*rd                         ! dsea
	 fx(-1) = fx(-1)*ru                         ! usea
	 fx(3) = fx(3)*rs                           ! strange
	 fx(-3) = fx(3)
	 fx(4) = fx(4)*rc                           ! charm
	 fx(-4) = fx(4)
	 fx(5) = fx(5)*rb                           ! bottom
	 fx(-5) = fx(5)
      else
         call mlmpdf8(ndns,ih,q2,xs,fx,nf)
      endif
      endif
c      write(*,*) ndns,ih,xs,fx(0),fx(1),fx(2),fx(3),fx(4)
      end

* This program convolutes the parton distribution functions of the
* photon with a photon distribution function, say a Weizsaecker-Williams,
* thereby producing a set of electron parton distrib. functions.
*
* MC, Nov. 1996
*

      subroutine ppgrv_int(x,qq,xpdf,nfl,ndns)
      implicit none
      integer nfl,ndns
      real * 8 x,qq,xpdf(-6:6)
c nag parameters
*      integer n,lck,lwrk
*      parameter (n=300)
*      parameter (LCK=N+4,LWRK=6*N+16)
*      real*8 coeff(lck,-6:6),lamda(lck)
*      common/splines/coeff,lamda
c
      integer ndnsc
      common/fonllelpdfc/ndnsc
c ww common
      integer itypeww
      real * 8 xmuww,thcww,zminww,zmaxww,elphen
      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
c
      integer j,ifail
c
      integer ini,ndnsold
      real * 8 qqold
      data ini/0/
      save ini, ndnsold, qqold
      if ((ini.eq.0).or.(qq.ne.qqold).or.ndns.ne.ndnsold)then
         ndnsc=ndns
         call pdg_spline(qq,nfl)
         ini=1
         qqold = qq
         ndnsold = ndns
      endif

      do j= -nfl,nfl
        ifail=1
        call getspline(x,j,xpdf(j),ifail)
        if(ifail.ne.0.and.x.gt..95d0) then
              xpdf(j) = 0d0
        elseif(ifail.ne.0.and.x.le..95d0) then
              print*,'**** spline fails ****',ifail,x,j
              stop
        endif
c.....back to the struct. funct. from its log
c        xpdf(j) = exp(xpdf(j))

c.....this to take care of fit oscillations beyond the endpoint.
      if(x.ge.zmaxww) xpdf(j) = 0.d0

c.....anti-flavour = flavour (and gluon=gluon, of course.....)
c.....occhio, mica vero per i valence quark....., cambia i do loop su
c.....j e togli questa riga, nel caso.
cc        xpdf(-j) = xpdf(j)
      enddo

c      print*,x
c      print*,xpdf

      end

      subroutine pdg_spline(qq,nfl)
      implicit none
      real * 8 qq
      integer nfl
      integer jfl,iret,iarg,jarg
      real * 8 xz,val,vv,ee
      integer n
      parameter        (n=300)
      real * 8 xx,scale
      common/xscale/xx,scale
      integer jj
      common/flv/jj
c ww common
      integer itypeww
      real * 8 xmuww,thcww,zminww,zmaxww,elphen
      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
c local static
      real * 8 abscis
      dimension abscis(0:n)
      real * 8 elpdf
      dimension elpdf(0:n,-6:6)
      real*8 a(n,-6:6),b(n,-6:6),c(n,-6:6),d(n,-6:6)
      integer interptp(-6:6)
c local
      real * 8 x,z,ylow,xmin
      integer j,i
c functions
      real * 8 dgauss,convol
      integer locatd
      external convol
      save abscis,elpdf,interptp
      save a,b,c,d
      scale=qq
      print*,' In pdg_spline, scale = ',scale

c.....log range for x: x from 10^-4 to 1 => log(x) from -9 to 0
      ylow = -9d0

      do j = -nfl,nfl
c log interpolation
         interptp(j)=1
         jj=j
 11      continue
         do i=0,n
           abscis(i) = log(zmaxww)+(ylow - ylow*float(i)/float(N+1))
           x = exp(abscis(i))
           xx=x
           xmin = max(x,zminww)
           elpdf(i,j) = x*dgauss(convol,xmin,zmaxww,1d-3)
           if(interptp(j).eq.1) then
              if(elpdf(i,j).le.0) then
                 interptp(j)=0
                 goto 11
              endif
              elpdf(i,j) = log(elpdf(i,j))
           endif
         enddo
         call dcspln(n,abscis,1,elpdf(0,j),n,0,
     #    a(1,j),b(1,j),c(1,j),d(1,j))
         if(interptp(j).eq.1) then
            print*,' Now fitting flavour ',j,' log interp.'
         else
            print*,' Now fitting flavour ',j,' lin. interp.'
         endif
      enddo
      return
      entry getspline(xz,jfl,val,iret)
c locatd returns
c i if z=absciss(i-1)
c 0 if z<absciss(0)
c -i if abscis(i-1)<z<abscis(i)
c -n+1 if z>ascis(n)
      z=log(xz)
      i=locatd(abscis(0),n+1,z)
      if(i.eq.0) then
         iret=-1
         return
      elseif(i.eq.-(n+1)) then
c above zmaxww
         val = 0
         iret=0
         return
      endif
c if z=abscis(n), i=(n+1); set it to -n (as if abscis(n-1)<z<abscis(n))
      if(i.eq.n+1) i=-n
      i=abs(i)-1
      val=a(i+1,jfl)+b(i+1,jfl)*(z-abscis(i))
     #  +c(i+1,jfl)*(z-abscis(i))**2
     #  +d(i+1,jfl)*(z-abscis(i))**3
      if(interptp(jfl).eq.1) val=exp(val)
      iret=0
      return
      entry getarrs(iarg,jarg,vv,ee)
      ee=elpdf(iarg,jarg)
      vv=exp(abscis(iarg))
      if(interptp(jarg).eq.1) ee=exp(ee) 
      end

      function convol(z)
      implicit none
      real * 8 convol,z
      real * 8 x,scale
      common/xscale/x,scale
      integer j
      common/flv/j
      real * 8 pdg(-6:6)
      integer ndnsc
      common/fonllelpdfc/ndnsc
      real * 8 y,q2
      real * 8 fww_ww
      y = x/z
      q2=scale**2
      call mlmpdf8(ndnsc,4,q2,y,pdg,6)
      convol = fww_ww(z)/z*pdg(j)
      end


