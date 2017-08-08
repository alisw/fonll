      implicit none
      real * 8 xmq,xmq2
      common/cmass/xmq,xmq2
      real * 8 xmh,xmh2
      integer mheqfl
      common/chmass/xmh,xmh2,mheqfl
      real * 8 ebeam2,ebeam1
      common/cenergies/ebeam1,ebeam2
      real * 8 sh,ycm
      common/clab/sh,ycm
      real * 8 ep
      common/cep/ep
      real * 8 xnorm
      common/cxnorm/xnorm
      integer imode
      common/cimode/imode
      integer ifrag,ifrframe
      common/cifrag/ifrag,ifrframe
      real * 8 eymin,eymax
      common/ceylim/eymin,eymax
      real * 8 ptmin,ptmax
      common/cptlim/ptmin,ptmax
      integer ncall0,nitn0
      common/cvegas/ncall0,nitn0
      integer kmode
      common/ckmode/kmode
      integer ipt,iimode
      real * 8 tmp,xl,xu,err,pt,pt2,ey,av,sd
      real * 8 fragfun,dgauss
      external fragfun
      open(unit=21,file='fragfonll.log')
      call loadresult
c
      write(6,*)
     # 'enter -1 for lower band, 0 for central, 1 for upper band'
      read(5,*) kmode
      write(21,*) kmode,'  ! -1 for low, 0 for central, 1 for high'
      write(6,*)'enter 0 for no fragmentation,'
      write(6,*)'enter 1 for Peterson fragmentation,'
      write(6,*)'enter 2 for (1-x) x^ep form'
      read(5,*) ifrag
      write(21,*) ifrag,'  ! 0:no fr., 1: Pet., 2: (1-x)x^ep'
      if(ifrag.ne.0) then
         write(6,*)'enter 0 to perform the fragmentation in the lab'
         write(6,*)
     # '      1 to perform the fragmentation at y=0 (default)'
         read(5,*)ifrframe
         write(21,*)ifrframe,'  ! 0 for fragm. in lab, 1 in y=0 frame'
         write(6,*)'enter heavy meson mass, <0 to set it = quark mass'
         read(5,*)xmh
         write(21,*)xmh,'    ! meson mass (<0 -> =quark mass'
         if(xmh.lt.0.d0) then
            mheqfl=1
         else
            mheqfl=0
         endif
         write(6,*)'enter Peterson parameter'
         read(5,*)ep
         write(21,*)ep,'     ! ep param. for fragm.'
c Compute the normalization of the fragmentation function
         xnorm=1.d0
         xl=0.d0
         xu=1.d0
         err=1.d-8
         tmp=dgauss(fragfun,xl,xu,err)
         xnorm=1/tmp
         write(6,*)'integral of the fragmentation function=',tmp
         write(6,*)' '
      else
         mheqfl=1
      endif
c
      ptmin=0.d0
      ptmax=0.d0
      eymin=0.d0
      eymax=0.d0
      write(6,*)'enter 0 to compute d sigma/d pt^2 d eta'
      write(6,*)'      1 to compute d sigma/d pt^2 d y'
      write(6,*)'      2 to compute d sigma/d pt^2, etamin<eta<etamax'
      write(6,*)'      3 to compute d sigma/d pt^2, ymin<y<ymax'
      write(6,*)'      4 to compute d sigma/d eta, ptmin<pt<ptmax'
      write(6,*)'      5 to compute d sigma/d y, ptmin<pt<ptmax'
      write(6,*)'      6 to compute sigma, cuts on pt and eta'
      write(6,*)'      7 to compute sigma, cuts on pt and y'
      write(6,*)'      8 to compute d sigma/d pt d eta'
      write(6,*)'      9 to compute d sigma/d pt d y'
      write(6,*)'     10 to compute d sigma/d pt, etamin<eta<etamax'
      write(6,*)'     11 to compute d sigma/d pt, ymin<y<ymax'
      read(5,*)iimode
      write(21,*) iimode,'   ! 0 for dsig/ddpt^2 d eta ...'
      open(unit=31,file='fragmfonll.dat',status='unknown')
      if(iimode.eq.0) then
      write(31,7)'pt2,       eta,       dsig/dpt^2deta,error'
      elseif(iimode.eq.1) then
      write(31,7)'pt2,       y,         dsig/dpt^2dy,error'
      elseif(iimode.eq.2) then
      write(31,7)'pt2,       etamin,    etamax,    dsig/d pt^2,error'
      elseif(iimode.eq.3) then
      write(31,7)'pt2,       ymin,      ymax,      dsig/d pt^2,error'
      elseif(iimode.eq.4) then
      write(31,7)'eta,       ptmin,     ptmax,     dsig/deta, error'
      elseif(iimode.eq.5) then
      write(31,7)'y,         ptmin,     ptmax,     dsig/d y,  error'
      elseif(iimode.eq.6) then
      write(31,7)
     # 'etamin,    etamax,    ptmin,     ptmax,     sigma,     error'
      elseif(iimode.eq.7) then
      write(31,7)
     # 'ymin,      ymax,      ptmin,     ptmax,     sigma,     error'
      elseif(iimode.eq.8) then
      write(31,7)'pt,        eta,       dsig/dptdeta,error'
      elseif(iimode.eq.9) then
      write(31,7)'pt,        y,         dsig/dptdy,error'
      elseif(iimode.eq.10) then
      write(31,7)'pt,        etamin,    etamax,    dsig/d pt, error'
      elseif(iimode.eq.11) then
      write(31,7)'pt,        ymin,      ymax,      dsig/d pt, error'
      endif
 7    format(a)
 33   continue
      imode=iimode
      if(imode.eq.0)then
        write(6,*)             'enter pt2, eta'
        read(5,*,err=999)pt2,ey
        write(21,*)      pt2,ey,'   ! pt2, eta'
        if(pt2.le.0) stop
        pt=sqrt(pt2)
      elseif(imode.eq.1)then
        write(6,*)             'enter pt2, y'
        read(5,*,err=999)pt2,ey
        write(21,*)      pt2,ey,'   ! pt2, y'
        if(pt2.le.0) stop
        pt=sqrt(pt2)
      elseif(imode.eq.2)then
        write(6,*)                      'enter pt2, etamin, etamax'
        read(5,*,err=999)pt2,eymin,eymax
        write(21,*)      pt2,eymin,eymax,'   ! pt2, etamin, etamax'
        if(pt2.le.0) stop
        pt=sqrt(pt2)
      elseif(imode.eq.3)then
        write(6,*)                      'enter pt2, ymin, ymax'
        read(5,*,err=999)pt2,eymin,eymax
         write(21,*)     pt2,eymin,eymax,'   ! pt2, ymin, ymax'
        if(pt2.le.0) stop
        pt=sqrt(pt2)
      elseif(imode.eq.4)then
        write(6,*)                     'enter eta, ptmin, ptmax'
        read(5,*,err=999)ey,ptmin,ptmax
        write(21,*)      ey,ptmin,ptmax,'   ! eta, ptmin, ptmax'
        if(ptmax.le.0) stop
      elseif(imode.eq.5)then
        write(6,*)                     'enter y, ptmin, ptmax'
        read(5,*,err=999)ey,ptmin,ptmax
        write(21,*)      ey,ptmin,ptmax,'   ! y, ptmin, ptmax'
        if(ptmax.le.0) stop
      elseif(imode.eq.6)then
        write(6,*)'enter etamin, etamax, ptmin, ptmax'
        read(5,*,err=999)eymin,eymax,ptmin,ptmax
        write(21,*)      eymin,eymax,ptmin,ptmax,
     #             '   ! etamin, etamax, ptmin, ptmax'
        if(ptmax.le.0) stop
      elseif(imode.eq.7)then
        write(6,*)'enter ymin, ymax, ptmin, ptmax'
        read(5,*,err=999)eymin,eymax,ptmin,ptmax
        write(21,*)      eymin,eymax,ptmin,ptmax,
     #             '   ! ymin, ymax, ptmin, ptmax'
        if(ptmax.le.0) stop
      elseif(imode.eq.8)then
        write(6,*)'enter pt, eta'
        read(5,*,err=999)pt,ey
        write(21,*) pt,ey,'   ! pt, eta'
        if(pt.le.0) stop
      elseif(imode.eq.9)then
        write(6,*)            'enter pt, y'
        read(5,*,err=999)pt,ey
        write(21,*)      pt,ey,'   ! pt, y'
        if(pt.le.0) stop
      elseif(imode.eq.10)then
        write(6,*)                     'enter pt, etamin, etamax'
        read(5,*,err=999)pt,eymin,eymax
        write(21,*)      pt,eymin,eymax,'   ! pt, etamin, etamax'
        if(pt.le.0) stop
      elseif(imode.eq.11)then
        write(6,*)                     'enter pt, ymin, ymax'
        read(5,*,err=999)pt,eymin,eymax
        write(21,*)     pt,eymin,eymax,'   ! pt, ymin, ymax'
        if(pt.le.0) stop
      else
        write(6,*)'ERROR: this mode is not implemented'
        stop
      endif
      if(ptmin.eq.0.d0)ptmin=1.d-6
c imode=8-11 similar to 0-3: treat them together
      ipt=2
      if(imode.ge.8)then
        ipt=1
        imode=imode-8
      endif
c Vegas parameters
      nitn0=10
      if(imode.lt.2)then
        ncall0=4000
      elseif(imode.ge.2.and.imode.lt.6)then
        ncall0=20000
      else
        ncall0=100000
      endif
c
      call fragdiff(pt,ey,av,sd)
      write(6,*)'   '
      write(6,*)'result=',av,' +- ',sd,' (only integration error)'
c Save results to a file
      if(ptmin.eq.1.d-6)ptmin=0.d0
      if(imode.lt.2)then
        if(ipt.eq.2)then
          write(31,'(4(d10.4,1x))') pt2,ey,av,sd
        else
          write(31,'(4(d10.4,1x))') pt,ey,2*pt*av,2*pt*sd
        endif
      elseif(imode.ge.2.and.imode.lt.4)then
        if(ipt.eq.2)then
          write(31,'(5(d10.4,1x))') pt2,eymin,eymax,av,sd
        else
          write(31,'(5(d10.4,1x))') pt,eymin,eymax,2*pt*av,2*pt*sd
        endif
      elseif(imode.ge.4.and.imode.lt.6)then
        write(31,'(5(d10.4,1x))') ey,ptmin,ptmax,av,sd
      else
        write(31,'(6(d10.4,1x))') eymin,eymax,ptmin,ptmax,av,sd
      endif
      goto 33
 999  continue
      end


      subroutine fragdiff(pt,ey,av0,sd0)
c wrapper to fragdiff0 to eventually find maximum and minimum of the cross section
      implicit none
      real * 8 pt,ey,av0,sd0
      integer kmode
      common/ckmode/kmode
c
      integer nymx,nptmx,maxfiles
      parameter(nymx=50,nptmx=250,maxfiles=10)
      real * 8 yval(nymx,maxfiles),zseq(nptmx,maxfiles),
     #         ypxsec(nymx*nptmx,maxfiles),
     #         xm(maxfiles),ffact(maxfiles),ebeam1(maxfiles),
     #         ebeam2(maxfiles),zmax1(maxfiles),zmax2(maxfiles),
     #         ptmax(maxfiles),bspli(nymx*nptmx,maxfiles)
	integer ny(maxfiles),npt(maxfiles),mode(maxfiles)
      integer lfiles
      common/cdsig/yval,zseq,ypxsec,xm,ffact,ebeam1,ebeam2,
     #             zmax1,zmax2,ptmax,bspli,ny,npt,lfiles,mode
c
      real * 8 xebeam2,xebeam1
      common/cenergies/xebeam1,xebeam2
      real * 8 xmq,xmq2
      common/cmass/xmq,xmq2
      real * 8 xmh,xmh2
      integer mheqfl
      common/chmass/xmh,xmh2,mheqfl
      real * 8 sh,ycm
      common/clab/sh,ycm
      integer jwhichf
      common/cjwhichf/jwhichf
c local variables
      real * 8 tmpav0,tmpsd0
      do jwhichf=1,lfiles
         xmq=xm(jwhichf)
         xmq2=xmq**2
         xebeam1=ebeam1(jwhichf)
         xebeam2=ebeam2(jwhichf)
         sh=4*xebeam1*xebeam2
         ycm=0.5d0*log(xebeam1/xebeam2)
         if(mheqfl.eq.1) then
            xmh=xmq
            xmh2=xmh**2
         endif
         call fragdiff0(pt,ey,tmpav0,tmpsd0)
         if(jwhichf.gt.1) then
            if(kmode*tmpav0.gt.kmode*av0) then
               av0=tmpav0
               sd0=tmpsd0
            endif
         else
            av0=tmpav0
            sd0=tmpsd0
         endif
c only one iteration in this case; first file has central values
         if(kmode.eq.0) goto 999
      enddo
 999  end


      subroutine fragdiff0(pt,ey,av0,sd0)
c Returns the cross section after fragmentation and integration(s):
c  imode=0 -> d sigma/d pt^2 d eta
c  imode=1 -> d sigma/d pt^2 d y
c  imode=2 -> d sigma/d pt^2 (etamin<eta<etamax)
c  imode=3 -> d sigma/d pt^2 (ymin<y<ymax)
c  imode=4 -> d sigma/d eta (ptmin<pt<ptmax)
c  imode=5 -> d sigma/d y (ptmin<pt<ptmax)
c  imode=6 -> sigma (ptmin<pt<ptmax, etamin<eta<etamax)
c  imode=7 -> sigma (ptmin<pt<ptmax, ymin<y<ymax)
c ey is either eta (imode=0,4) or y (imode=1,5), or is unused
      implicit none
      real * 8 pt,ey,av0,sd0
      integer imode
      common/cimode/imode
      real * 8 pt0,ey0
      common/kine/pt0,ey0
      integer ncall0,nitn0
      common/cvegas/ncall0,nitn0
      real*8 xl,xu,acc
      integer ndim,ncall,itmx,nprn
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      integer ifrag,ifrframe
      common/cifrag/ifrag,ifrframe
      integer j
      real*8 av,sd,chi2,tmp
      real*8 xintegrand
      external xintegrand
c
      pt0=pt
      ey0=ey
      acc=1.d-4
      do j=1,10
        xl(j)=0.d0
        xu(j)=1.d0
      enddo
      nprn=1
      itmx=nitn0
      ncall=ncall0
      if(imode.lt.2)then
        ndim=1
      elseif(imode.ge.2.and.imode.lt.6)then
        ndim=2
      else
        ndim=3
      endif
      if(ifrag.eq.0) ndim=ndim-1
      if(ndim.eq.0) then
         av=xintegrand(tmp,1.d0)
         sd=0
      else
         call vegas(xintegrand,av,sd,chi2)
      endif
      av0=av
      sd0=sd
      return
      end


      function xintegrand(xx,weight)
c This is the function that Vegas integrates
      implicit none
      real*8 xintegrand,weight
      real*8 xx(*)
      real * 8 pt0,ey0
      common/kine/pt0,ey0
      integer imode
      common/cimode/imode
      real * 8 eymin,eymax
      common/ceylim/eymin,eymax
      real * 8 ptmin,ptmax
      common/cptlim/ptmin,ptmax
      real * 8 xmq,xmq2
      common/cmass/xmq,xmq2
      real * 8 xmh,xmh2
      integer mheqfl
      common/chmass/xmh,xmh2,mheqfl
      integer ifrag,ifrframe
      common/cifrag/ifrag,ifrframe
      integer iret,izind
      real * 8 xjac,zfragmin,zfrag,y,ey,eta,xmt,pt
      real * 8 eta_to_y,xsecfrag
c
      xjac=1.d0
c index of fragmentation variable: defaults to xx(2)
      izind=2
      if(imode.eq.0)then
c ey0 is eta
        pt=pt0
        eta=ey0
        y=eta_to_y(ey0,pt,xmh2)
c index of fragmentation variable: defaults to xx(2)
        izind=1
      elseif(imode.eq.1)then
c ey0 is y
        pt=pt0
        y=ey0
        izind=1
      elseif(imode.eq.2.or.imode.eq.3)then
c ey is either eta (imode=2) or y (imode=3); ey0 is unused
        pt=pt0
        ey=eymin+(eymax-eymin)*xx(1)
        xjac=xjac*(eymax-eymin)
        if(imode.eq.2)then
          eta=ey
          y=eta_to_y(ey,pt,xmh2)
        else
          y=ey
        endif
        izind=2
      elseif(imode.eq.4.or.imode.eq.5)then
c ey0 is either eta (imode=4) or y (imode=5); pt0 is unused
        pt=ptmin+(ptmax-ptmin)*xx(1)**2
        xjac=xjac*(ptmax-ptmin)*2*xx(1)
        if(imode.eq.4)then
          eta=ey0
          y=eta_to_y(ey0,pt,xmh2)
        else
          y=ey0
        endif
      elseif(imode.eq.6.or.imode.eq.7)then
c ey is either eta (imode=6) or y (imode=7); pt0 is unused
        pt=ptmin+(ptmax-ptmin)*xx(1)**2
        xjac=xjac*(ptmax-ptmin)*2*xx(1)
        ey=eymin+(eymax-eymin)*xx(2)
        xjac=xjac*(eymax-eymin)
        if(imode.eq.6)then
          eta=ey
          y=eta_to_y(ey,pt,xmh2)
        else
          y=ey
        endif
        izind=3
      endif
      call szmin(pt,y,zfragmin,iret)
      if(iret.eq.0)then
        if(zfragmin.lt.0.d0.or.zfragmin.gt.1.d0)then
          write(6,*)'ERROR: fatal error in xintegrand: zmin=',
     #              zfragmin
          stop
        endif
        if(ifrag.eq.0) then
c no fragmentation
           zfrag=1
        else
           zfrag=1.d0-(1.d0-zfragmin)*xx(izind)**2
           xjac=xjac*(1.d0-zfragmin)*2*xx(izind)
        endif
c The cross section is given as d/d pt^2 dy; multiply for the relevant jacobian
        if(imode.eq.0.or.imode.eq.2.or.imode.eq.4.or.imode.eq.6)then
c This is dy/deta
          xmt=sqrt(pt**2+xmh2)
          xjac=xjac*pt*cosh(eta)/(xmt*cosh(y))
        endif
        if(imode.ge.4)then
c This is dpt2/dpt
          xjac=xjac*2*pt
        endif
        xintegrand=xjac*xsecfrag(pt,y,zfrag)
      else
        xintegrand=0.d0
      endif
      return
      end


      function xsecfrag(pt,y,zfrag)
c Master function: xsecfrag == dsigma/dpt^2 dy, with pt and y relevant
c to the fragmented quark, to be integrated over the relevant range. It 
c is given by the product of the fragmentation function (fragfun), times
c the cross section for bare quarks (dsigdpt2dy), times a suitable
c jacobian, required to pass from the variables relevant to the
c fragmented quark (pt, y) to those relevant to the bare 
c quark (pthat, yhat). The starting point is the formula
c   sigma(p) = \int dz frag(z) sigma_bare(p^),
c and d/d^2pt dy = E d^3/d^3 p. p^ is defined as follows:
c   p^ = B^-1( B(p)/z ).
c In B(p)/z it is understood that only the 3-vector is divided by z,
c while the energy is obtained through a mass-shell condition (the
c mass is always that of the quark). B is a suitable boost. For ifrframe=0 
c then B==1, and for ifrframe=1 then B is such that B(p) has rapidity set 
c equal to zero
      implicit none
      real * 8 xsecfrag,pt,y,zfrag
      real * 8 fragfun,dsigdpt2dy
      real * 8 sh,ycm
      common/clab/sh,ycm
      real * 8 xmq,xmq2
      common/cmass/xmq,xmq2
      real * 8 xmh,xmh2
      integer mheqfl
      common/chmass/xmh,xmh2,mheqfl
      real * 8 ebeam2,ebeam1
      common/cenergies/ebeam1,ebeam2
      integer ifrag,ifrframe
      common/cifrag/ifrag,ifrframe
      real * 8 eh,plh,pthat,ehat,plhat,yhat,xmt,xmth,xjac,epart
      integer kmode
      common/ckmode/kmode
c
      xmt=sqrt(pt**2+xmq2) 
      xmth=sqrt(pt**2+xmh2) 
      eh  =xmth*cosh(y)
      plh =xmth*sinh(y)
      pthat=pt/zfrag
      if(ifrframe.eq.0)then
        plhat=plh/zfrag
        ehat = sqrt(pthat**2+plhat**2+xmq2)
        yhat=0.5d0*log( (ehat+plhat)/(ehat-plhat) )
      elseif(ifrframe.eq.1)then
        yhat=y
        ehat=sqrt(pthat**2+xmq2)*cosh(yhat)
      else
        write(6,*)'ERROR: option not implemented: ifrframe=',ifrframe
        stop
      endif
c Check if zfrag is not out of range
      epart=sqrt(pthat**2+xmq2)*cosh(yhat-ycm)
      if(epart.le.sqrt(sh)/2.d0)then
c xjac=d^3 p^/d^3p * E/E^
        if(ifrframe.eq.0)then
          xjac=eh/(zfrag**3*ehat)
        else
          xjac=1/zfrag**2
        endif
        xsecfrag=xjac*fragfun(zfrag)*
     #           dsigdpt2dy(pthat,yhat)
      else
        xsecfrag=0.d0
      endif
      return
      end


      function eta_to_y(eta,pt,xmq2)
      implicit none
      real * 8 eta_to_y,eta,pt,xmq2,chy,acosh,tmp
c
      chy=sqrt( (pt**2*cosh(eta)**2+xmq2)/(pt**2+xmq2) )
c What follows holds since sinh(eta)=mt/pt*sinh(y)
      if(eta.gt.0)then
        tmp=acosh(chy)
      else
        tmp=-acosh(chy)
      endif
      eta_to_y=tmp
      return
      end


      function y_to_eta(y,pt,xmq2)
      implicit none
      real * 8 y_to_eta,y,pt,xmq2,cheta,acosh,tmp
c
      cheta=sqrt( (pt**2+xmq2)*cosh(y)**2-xmq2)/pt
c What follows holds since sinh(eta)=mt/pt*sinh(y)
      if(y.gt.0)then
        tmp=acosh(cheta)
      else
        tmp=-acosh(cheta)
      endif
      y_to_eta=tmp
      tmp = log(sqrt(cosh(y)**2+xmq2/pt**2*sinh(y)**2)
     #   +sqrt(1+xmq2/pt**2)*sinh(y))
      if(abs(tmp-y_to_eta).gt.1.d-7) then
         write(*,*) ' y_to_eta: troubles...'
         stop
      endif
      return
      end


      function acosh(chx)
      implicit none
      real * 8 acosh,chx
c
      acosh=log(chx+sqrt(chx**2-1))
      return
      end


      function fragfun(z)
      implicit none
      real * 8 fragfun,z
      real * 8 ep
      common/cep/ep
      real * 8 xnorm
      common/cxnorm/xnorm
      integer ifrag,ifrframe
      common/cifrag/ifrag,ifrframe
      real * 8 tmp
      if(ifrag.eq.0) then
         fragfun=1
      elseif(ifrag.eq.1) then
c
c Use this form to avoid numerical problems at z=0 and z=1
         tmp=z*(1-z)**2/((1-z)**2+ep*z)**2
         fragfun=xnorm*tmp
      elseif(ifrag.eq.2) then
         fragfun=xnorm*(1-z)*z**ep
      endif
      end


      subroutine szmin(pt,y,zmin,iret)
c Returns the lower bound in the z integration range, which is obtained
c by imposing E_quark<sqrt{S}/2 in CM frame of the the collinding particles,
c when ifrframe=1, and a loose (i.e. less than the real one) lower
c bound by imposing E_quark<sqrt{S} in the lab frame when ifrframe=0.
c The routine checks if pt and y are inside (iret=0) or outside (iret=1) 
c the physical region. Remember that pt and y are meson variables; their
c relationship to the quark variables is given in function xsecfrag
      implicit none
      integer iret
      real * 8 pt,y,zmin
      real * 8 sh,ycm
      integer ifrag,ifrframe
      common/cifrag/ifrag,ifrframe
      common/clab/sh,ycm
      real * 8 xmq,xmq2
      common/cmass/xmq,xmq2
      real * 8 xmh,xmh2
      integer mheqfl
      common/chmass/xmh,xmh2,mheqfl
      real * 8 ypart
c
      ypart=y-ycm
      iret=0
      if((pt**2+xmh2)*cosh(ypart)**2.gt.sh/4.d0)then
        iret=1
        zmin=-1.d0
      else
        iret=0
        if(ifrframe.eq.0)then
          zmin=sqrt((pt**2*cosh(y)**2+xmh2*sinh(y)**2)/(sh-xmq2))
        elseif(ifrframe.eq.1)then
          zmin=pt/sqrt(sh/(4*cosh(ypart)**2)-xmq2)
        else
          write(6,*)'ERROR in szmin: ifrframe=',ifrframe
          stop
        endif
      endif
      return
      end


      subroutine loadresult
      implicit none
      integer itype,iproc,iret
      real * 8 cparam
c
      integer nymx,nptmx,maxfiles
      parameter(nymx=50,nptmx=250,maxfiles=10)
      real * 8 yval(nymx,maxfiles),zseq(nptmx,maxfiles),
     #         ypxsec(nymx*nptmx,maxfiles),
     #         xm(maxfiles),ffact(maxfiles),ebeam1(maxfiles),
     #         ebeam2(maxfiles),zmax1(maxfiles),zmax2(maxfiles),
     #         ptmax(maxfiles),bspli(nymx*nptmx,maxfiles)
	integer ny(maxfiles),npt(maxfiles),mode(maxfiles)
      integer lfiles
      common/cdsig/yval,zseq,ypxsec,xm,ffact,ebeam1,ebeam2,
     #             zmax1,zmax2,ptmax,bspli,ny,npt,lfiles,mode
c
      write(*,*)'enter 0 for hadroproduction, 1 for photoproduction'
      read(*,*) iproc
      write(21,*) iproc,' ! 0 for H-H, 1 for ph/el-H'
      if(iproc.eq.1) then
         write(*,*)'enter 0 for matched cross section'
         write(*,*)'      1 for pointlike massive'
         write(*,*)'      2 for hadronic massive'
         write(*,*)'      3 for pointlike massless'
         write(*,*)'      4 for hadronic massless'
         write(*,*)'      5 for pointlike resummed'
         write(*,*)'      6 for hadronic resummed'
         write(*,*)'      7 for pointlike+hadronic massive'
         write(*,*)'      8 for pointlike+hadronic massless'
         write(*,*)'      9 for pointlike+hadronic resummed'
      else
         write(*,*)'enter 0 for matched cross section'
         write(*,*)'      1 for hadronic massive'
         write(*,*)'      2 for hadronic massless'
         write(*,*)'      3 for hadronic resummed'
      endif         
      read(*,*)itype
      write(21,*)itype,' ! 0 for matched, etc.'
      if(itype.eq.0)then
        write(*,*)'enter c, parameter for matching'
        read(*,*)cparam
        write(21,*)cparam,'  ! c param. for matching'
      endif
      write(6,*)'enter names of file with data (one per line)'
      write(6,*)'an empty line to terminate.'
      write(6,*)'the first file must correspond to the central value'
      write(6,*)'the remaining files are used to compute the upper'
      write(6,*)'and lower range of the cross section'
      do lfiles=1,maxfiles
         call loadfile(itype,iproc,cparam,iret,yval(1,lfiles),
     #        zseq(1,lfiles),ypxsec(1,lfiles),xm(lfiles),
     #        ffact(lfiles),ebeam1(lfiles),ebeam2(lfiles),zmax1(lfiles),
     #        zmax2(lfiles),ptmax(lfiles),ny(lfiles),npt(lfiles),
     #        mode(lfiles))
         if(iret.ne.0) then
            goto 10
         endif
      enddo
      write(*,*) ' too many files'
      stop
 10   continue
      lfiles=lfiles-1
      end

      subroutine loadfile(itype,iproc,cparam,iret,yval,zseq,
     #     ypxsec,xm,ffact,ebeam1,ebeam2,zmax1,zmax2,ptmax,ny,npt,
     #     mode)
      implicit none
c parameters
      integer narr
      integer nline
      parameter (nline=200)
      parameter (narr=10)
      integer nymx,nptmx,nmx
      parameter (nymx=50,nptmx=250)
      parameter (nmx=nymx*nptmx)
c arguments
      integer itype,iproc,iret,mode
      real * 8 cparam
      real * 8 yval(nymx),zseq(nptmx),ypxsec(nymx*nptmx),
     #         xm,ffact,ebeam1,ebeam2,zmax1,zmax2,ptmax
      integer npt,ny
c common 
      integer jptseqfl
      common/cjptseqfl/jptseqfl
c local
      character * 70 file
      character * 200 line
      real * 8 rarr(narr)
      integer karr,icount,ii,jj,kk,i,j,iy,ipt
      real * 8 xentries(10,nmx)
      real * 8 errxsec(nymx,nptmx)
      real * 8 pt,ptn,y,yn,xtmp,xm2,tmp1,tmp2,
     #  ymax,ymin,z
c statement functions
      integer ind1,ind2
c function to get index in 1-dim array from 2-dim indices
c      ind2(i,j)=(i-1)*ny+j
      ind1(i,j)=(i-1)*npt+j
      ind2(i,j)=i+(j-1)*ny
c Set this to 1 (interpolation parameters not yet computed)
      mode=1
      write(*,*) ' enter file (empty line to terminate)'
      read(5,'(a)') file
      write(21,'(a)') file
      if(file.eq.'') then
         iret=-1
         return
      endif
c Read the results produced by fonllgrid.f.
      open(unit=11,file=file,status='old')
      read(11,*) ebeam1,zmax1,ebeam2,zmax2,ffact,xm,ptmax,npt
      xm2=xm**2
      icount=1
      do while(icount.le.nmx)
         read(11,'(a)',err=199,end=200) line
         call reads(line,nline,rarr,narr,karr)
         if(karr.gt.narr) then
            write(*,*) ' error: more than narr elements in each row'
            stop
         else
            do ii=1,karr
               xentries(ii,icount)=rarr(ii)
            enddo
            do ii=karr+1,narr
c     set to zero (photoproduction stuff)
               xentries(ii,icount)=0
            enddo
         endif
         icount=icount+1
      enddo
c     if file is finished, jump
      read(11,*,end=200)
      write(*,*) 'ERROR: found more than ',nmx,' values'
      stop
 199  continue
      write(*,*) 'ERROR: cannot read line ',icount
      stop
 200  continue
      close(11)
      icount=icount-1
c     Now order the results in y and pt
 201  continue
      write(6,*)'icount=',icount
      do ii=1,icount-1
         do jj=1,icount-ii
            pt=xentries(1,jj)
            y =xentries(2,jj)
            ptn=xentries(1,jj+1)
            yn =xentries(2,jj+1)
            if(y.eq.yn.and.pt.eq.ptn) then
               write(6,*)
     #              'WARNING: at least two pairs (y,pt) are equal'
               xtmp=0.d0
               do kk=1,karr
                  xtmp=xtmp+abs(xentries(kk,jj)-xentries(kk,jj+1))
               enddo
               if(xtmp.gt.0.d0)then
                  write(6,*)
     #                 'ERROR: The computed cross sections disagree'
                  stop
               else
                  write(6,*)'dropping one entry'
               endif
               call delete(xentries,jj+1,icount)
               goto 201
            endif
            if(y.gt.yn.or.(y.eq.yn.and.pt.gt.ptn))then
               call swap(xentries,jj,jj+1)
            endif
         enddo
      enddo
c     Checks of consistency;
c     first count pt points
      npt=2
      do while(npt.le.nmx.and.xentries(2,npt).eq.xentries(2,1))
         npt=npt+1
      enddo
      npt=npt-1
      ny=icount/npt
c     icount must be a multiple of npt
      if(ny*npt.ne.icount)then
         write(6,*)'ERROR: Mismatch # 1 in input'
         stop
      endif
      write(6,*)'npt,ny=',npt,ny
c     check that subsequent npt entries have the same rapidity
      do iy=1,ny
         do ipt=2,npt
            if(xentries(2,ind1(iy,ipt))
     #           .ne.xentries(2,ind1(iy,1))) then
               write(6,*)'ERROR: Mismatch # 1 in input'
               stop
            endif
         enddo
      enddo
c fill rapidity array
      ymax=yval(1)
      ymin=ymax
      do iy=1,ny
         yval(iy)=xentries(2,ind1(iy,1))
         ymax=max(ymax,yval(iy))
         ymin=min(ymin,yval(iy))
      enddo
c check for consistent pt points
      jptseqfl=0
      do iy=1,ny
         do ipt=1,npt
            z=dble(ipt-1)/(npt-1)
            call ptpoints(ebeam1,zmax1,ebeam2,zmax2,
     #           ffact,xentries(2,ind1(iy,ipt)),xm,ptmax,0,z,pt)
            if(abs(pt-xentries(1,ind1(iy,ipt)))/pt.gt.1.e-4)then
               jptseqfl=1
               goto 11
            endif
         enddo
      enddo
 11   continue
      if(jptseqfl.ne.0) then
         do iy=1,ny
            do ipt=1,npt
               if(xentries(1,ind1(1,ipt)).ne.xentries(1,ind1(iy,ipt)))
     #              then
                  write(*,*) ' Error in loadfile: mismatch in input'
                  stop
               endif
            enddo
         enddo
      endif
c set pt parameters
      do ipt=1,npt
         pt=xentries(1,ind1(1,ipt))
         call ptpoints(ebeam1,zmax1,ebeam2,zmax2,
     #         ffact,xentries(2,ind1(1,ipt)),xm,ptmax,1,zseq(ipt),pt)
      enddo
c     Define the matched cross section
      do iy=1,ny
         do ipt=1,npt
            kk=ind1(iy,ipt)
            call xsecmtch(itype,iproc,cparam,
     #           xentries(1,kk),xm2,
     #           xentries(3,kk),xentries(4,kk),
     #           xentries(5,kk),xentries(6,kk),
     #           xentries(7,kk),xentries(8,kk),
     #           xentries(9,kk),xentries(10,kk),
     #           tmp1,tmp2)
            ypxsec(ind2(iy,ipt))=tmp1
            errxsec(iy,ipt)=tmp2
         enddo
      enddo
c     check for too large errors.
c Look for >5% errors, relative to the maximum value of the
c cross section near each point
c     do iy=1,ny
c     do ipt=1,npt
      end

      subroutine delete(xentries,index,icount)
      real * 8 xentries(10,*)
      integer index,icount
      integer l,k
      do k=index,icount-1
         do l=1,10
            xentries(l,k)=xentries(l,k+1)
         enddo
      enddo
      icount=icount-1
      end

      subroutine swap(xentries,jj,kk)
      implicit none
      integer jj,kk
      real * 8 xentries(10,*)
      real * 8 tmp
      integer l
      do l=1,10
         tmp=xentries(l,kk)
         xentries(l,kk)=xentries(l,jj)
         xentries(l,jj)=tmp
      enddo
      end

      subroutine xsecmtch(itype,iproc,cparam,
     #      pt,xm2,hdmv,hdml,hdrs,ehrs,phmv,phml,phrs,eprs,
     #                    xsec,err)
c Eq.(5.2) when itype=0
      implicit none
      integer itype,iproc
      real * 8 cparam,pt,xm2,phmv,hdmv,phml,hdml,phrs,hdrs
      real * 8 eprs,ehrs,xsec,err
      real * 8 xsmear
c
      err=0.d0
      if(iproc.eq.1) then
         if(itype.eq.0)then
            xsec=phmv+hdmv+xsmear(pt,xm2,cparam)*(phrs+hdrs-phml-hdml)
            err=xsmear(pt,xm2,cparam)*sqrt(eprs**2+ehrs**2)
         elseif(itype.eq.1)then
            xsec=phmv
         elseif(itype.eq.2)then
            xsec=hdmv
         elseif(itype.eq.3)then
            xsec=phml
         elseif(itype.eq.4)then
            xsec=hdml
         elseif(itype.eq.5)then
            xsec=phrs
            err=eprs
         elseif(itype.eq.6)then
            xsec=hdrs
            err=ehrs
         elseif(itype.eq.7)then
            xsec=phmv+hdmv
         elseif(itype.eq.8)then
            xsec=phml+hdml
         elseif(itype.eq.9)then
            xsec=phrs+hdrs
            err=sqrt(eprs**2+ehrs**2)
         else
            write(6,*)'wrong itype:',itype
            stop
         endif
      elseif(iproc.eq.0) then
         if(itype.eq.0)then
            xsec=hdmv+xsmear(pt,xm2,cparam)*(hdrs-hdml)
            err=xsmear(pt,xm2,cparam)*ehrs
         elseif(itype.eq.1)then
            xsec=hdmv
         elseif(itype.eq.2)then
            xsec=hdml
         elseif(itype.eq.3)then
            xsec=hdrs
         else
            write(6,*)'wrong itype:',itype
            stop
         endif
      else
         write(*,*) ' wrong iproc'
         stop
      endif
      end

      function xsmear(pt,xm2,cparam)
c Eq.(5.1), with c==cparam
      implicit none
      real * 8 xsmear,pt,xm2,cparam
c
      xsmear=pt**2/(pt**2+(cparam)**2*xm2)
      return
      end

      subroutine strnum(string,num)
c- writes the number num on the string string starting at the blank
c- following the last non-blank character
      character * (*) string
      character * 20 tmp
      l = len(string)
      write(tmp,'(i15)')num
      j=1
      dowhile(tmp(j:j).eq.' ')
        j=j+1
      enddo
      ipos = istrl(string)
      ito = ipos+1+(15-j)
      if(ito.gt.l) then
         write(*,*)'error, string too short'
         write(*,*) string
         stop
      endif
      string(ipos+1:ito)=tmp(j:)
      end

      function istrl(string)
c returns the position of the last non-blank character in string
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      istrl = i
      end

      subroutine strcat(str1,str2,str)
c concatenates str1 and str2 into str. Ignores trailing blanks of str1,str2
      character *(*) str1,str2,str
      l1=istrl(str1)
      l2=istrl(str2)
      l =len(str)
      if(l.lt.l1+l2) then
          write(*,*) 'error: l1+l2>l in strcat'
          write(*,*) 'l1=',l1,' str1=',str1
          write(*,*) 'l2=',l2,' str2=',str2
          write(*,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end

      function dsigdpt2dy(pt,y)
      implicit none
      real * 8 dsigdpt2dy,pt,y
c
      integer nymx,nptmx,maxfiles
      parameter(nymx=50,nptmx=250,maxfiles=10)
      real * 8 yval(nymx,maxfiles),zseq(nptmx,maxfiles),
     #         ypxsec(nymx*nptmx,maxfiles),
     #         xm(maxfiles),ffact(maxfiles),ebeam1(maxfiles),
     #         ebeam2(maxfiles),zmax1(maxfiles),zmax2(maxfiles),
     #         ptmax(maxfiles),bspli(nymx*nptmx,maxfiles)
	integer ny(maxfiles),npt(maxfiles),mode(maxfiles)
      integer lfiles
      common/cdsig/yval,zseq,ypxsec,xm,ffact,ebeam1,ebeam2,
     #             zmax1,zmax2,ptmax,bspli,ny,npt,lfiles,mode
c
      integer jwhichf
      common/cjwhichf/jwhichf
      integer ier,j
      real * 8 z,val
      if(jwhichf.gt.lfiles.or.jwhichf.lt.1) then
         write(*,*) ' dsigdpt2dy: jwhichf=',jwhichf
         stop
      endif
      call ptpoints(ebeam1(jwhichf),zmax1(jwhichf),
     # ebeam2(jwhichf),zmax2(jwhichf),
     # ffact(jwhichf),y,xm(jwhichf),ptmax(jwhichf),1,z,pt)
      if(y.lt.yval(1,jwhichf)) then
c         write(*,*) ' requested point below grid in y'
c         stop
         goto 998
      elseif(y.gt.yval(ny(jwhichf),jwhichf)) then
c         write(*,*) ' requested point above grid in y'
c         stop
         goto 998
      endif
      if(z.lt.zseq(1,jwhichf)) then
c         write(*,*) ' requested point below grid in pt'
c         stop
         goto 998
      elseif(z.gt.zseq(npt(jwhichf),jwhichf)) then
c         write(*,*) ' requested point above grid in pt'
c         stop
         goto 998
      endif
      if(mode(jwhichf).eq.1) then
         do j=1,ny(jwhichf)*npt(jwhichf)
            if(ypxsec(j,jwhichf).le.0) ypxsec(j,jwhichf)=1.d-25
            ypxsec(j,jwhichf)=log(ypxsec(j,jwhichf))
         enddo
      endif
      call rgbi3p(mode(jwhichf),ny(jwhichf),npt(jwhichf),
     #  yval(1,jwhichf),zseq(1,jwhichf),ypxsec(1,jwhichf),
     #  1,y,z,val,ier,bspli(1,jwhichf))
      if(ier.ne.0) then
         write(*,*) ' bispline failure'
         stop
      endif
      mode(jwhichf)=2
      dsigdpt2dy=exp(val)
      goto 999
 998  dsigdpt2dy=0
 999  end

         

c Program to read numbers from strings
      subroutine reads(string,nstr,rarr,narr,karr)
c string: input string of nstr characters
c it outputs in rarr the real numbers in string,
c karr is the number of real numbers found
      implicit none
      integer nstr,narr,karr
      real * 8 rarr(narr)
      character *(*) string
c
      character * 15 numstr
      integer istr,istart,ios
      karr=0
      istr=1
 1    continue
c no more numbers:
      if(istr.gt.nstr) goto 999
c skip blanks
 10   if(string(istr:istr).eq.' ') then
         if(istr.eq.nstr) goto 999
         istr=istr+1
         goto 10
      endif
c first non-blank
      istart=istr
c find next non-blank
 20   if(string(istr:istr).ne.' ') then
         if(istr.eq.nstr) goto 22
         istr=istr+1
         goto 20
      endif
      istr=istr-1
 22   continue
      if(istr-istart+1.gt.15) then
         write(*,*) ' error: number too long'
         stop
      endif
      numstr=string(istart:istr)
      karr=karr+1
      read(unit=numstr,iostat=ios,fmt='(bn,f15.0)') rarr(karr)
c if not a valid number ignore silently
      if(ios.ne.0) then
         karr=karr-1
         goto 999
      endif
      istr=istr+1
      goto 1
 999  end

c Program to read integers from strings
      subroutine ireads(string,nstr,iarr,narr,karr)
      implicit none
      integer narr
      integer iarr(narr)
      character * (*) string
      integer nstr,karr,j
c local variables
      real * 8 rarr(20)
      call reads(string,nstr,rarr,20,karr)
      if(karr.gt.narr.or.(narr.gt.20.and.karr.eq.20)) then
         write(*,*) ' more than 20 numbers in input?'
         stop
      endif
      do j=1,karr
         iarr(j)=nint(rarr(j))
         if(rarr(j)-iarr(j).ne.0) then
            write(*,*) ' not an integer in input!'
            stop
         endif
      enddo
      end

c program to read numbers from input line
      subroutine iread(iun,iarr,nel,nread)
      dimension iarr(*)
      character * 300 str
      read(iun,'(a)') str
      call ireads(str,300,iarr,nel,nread)
      end

c program to read numbers from input line
      subroutine rread(iun,rarr,nel,nread)
      real * 8  rarr(*)
      character * 300 str
      read(iun,'(a)') str
      call reads(str,300,rarr,nel,nread)
      end
      
