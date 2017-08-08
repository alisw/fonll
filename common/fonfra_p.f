      subroutine fonfra(ifragmode,xc,ipi,qp2
     #     ,xdup,xdubp,xddp,xddbp,xdsp,xdcp
     #     ,xdcbp,xdbp,xdbbp,xdtp,xdtbp,xdgp)
      implicit none
c     This subroutine returns x*fragmentation functions at the scale
c     qp2, x=xc. The variable ipi IS NOT USED
c     xdup is for u, xdubp is for ubarr, xddp is
c     for d, xddbp is for dbarr, xdsp is for s, xdcp is for c, xdbp
c     is for b, xdbbp is for bbarr, xdtp is for t, xdtbp is for tbarr,
c     and xdgp is for g.
c The flag ifragmode controls the behaviour of the function:
c ifragmode=0: Normal behaviour
c ifragmode=1: only the density of the selected heavy flavour (b or c) is
c              non-zero, and is the same as in ifragmode=0
c ifragmode=2: only the density of the selected heavy flavour (b or c) is
c              non-zero, and is equal to 1 (delta function in x space)
c
c The selected heavy quark is the value of common/hvqtype/hvqs,
c equal to 'b' or 'c'.
c
c Common blocks to set before calling this function
c    common/alfabeta/alfa,beta,npsm
c         Parameters for non pert. smearing
c         npsm=1: do pert. smearing x^beta*(1-x)^alfa
c    common/asscale/assc
c         Factor to reduce as (for testing only)
c    common/chat/verbose
c         verbosity for testing
c    common/cmu0/cmu0
c         initial evolution scale factor
c         maybe should be set to 1 (it is never varied)
c    common/flnumb/nfl
c         number of flavours
c    common/frschemec/frscheme
c         Fragmentation scheme
c    common/hvqmass/xmb,xmc
c         masses of heavy flavour
c    common/hvqtype/hvqs
c         'b' or 'c'
c    common/lamqcd/alam5qcd
c          Lambda_MSBAR_5 value
c    common/sudakov/isuda
c          include or not sudakov resummation
c
c         
c    in xsect.h:
c      iloopas,iloopfr,iwhichsfh
c      ialtev
c
c    common/bndary/xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2,nnf
c    common/integr/step,ermx,nmax
c    common/spline/cheavy,cgluon,clight,cantih,lamda
c    common/xvalue/x
c    common/xvalue/xx                                       
c
c
c
      integer ifragmode,ipi
      real * 8 xc,qp2,xdup,xdubp,xddp,xddbp,xdsp,xdcp,
     #     xdcbp,xdbp,xdbbp,xdtp,xdtbp,xdgp
      integer m1,m2,m3,m,lck,lwrk
      parameter (m1=0,m2=50,m3=50)
      parameter (m=m1+m2+m3,lck=m+4,lwrk=6*m+16)      
      real * 8  x(1)
      character*1 hvqs
      common/hvqtype/hvqs
c save values in case of repeated calls with same arguments
      real * 8 xcold,xdupold,xdbpold,xdgpold,xdbbpold,xdcpold,xdcbpold
      save xcold,xdupold,xdbpold,xdgpold,xdbbpold,xdcpold,xdcbpold
c Communication with interpolating functions.
      real * 8 cheavy(lck),cgluon(lck),clight(lck),cantih(lck),
     #     lamda(lck)
      common/spline/cheavy,cgluon,clight,cantih,lamda
      integer ini
      real * 8 qp2old,xdhq,xdha
      integer ifailh,ifaila,ifailg,ifaill,ifail
      external e02bbf
      data ini/0/
      save ini,qp2old
      if(ifragmode.eq.2) then
c return 1 for the selected heavy quark, 0 for all others
         if(hvqs.eq.'b') then
            xdbp = 1d0
            xdcp = 0d0
         else
            xdbp = 0d0
            xdcp = 1d0
         endif
         xdgp = 0d0
         xdup = 0d0
         xdubp = xdup
         xddp = xdup
         xddbp = xdup
         xdsp = xdup
         xdcbp = xdup
         xdbbp = xdup
         xdtp = 0d0
         xdtbp = 0d0
         return
      endif
c use previous values in this cases
      if(ini.ne.0.and.xc.eq.xcold.and.qp2.eq.qp2old) then
         xdbp = xdbpold
         xdgp = xdgpold
         xdup = xdupold
         xdubp = xdup
         xddp = xdup
         xddbp = xdup
         xdsp = xdup
         xdcp = xdcpold
         xdcbp = xdcbpold
         xdbbp = xdbbpold
         xdtp = 0d0
         xdtbp = 0d0
         goto 998
      endif
c At the first call, or if qp2 has changed, call the routine that prepares
c the x grid, setting the common/spline/ arrays
      if ((ini.eq.0).or.(qp2.ne.qp2old))then
         call fonfra_spline(qp2)
         ini=1
         qp2old=qp2
      endif
      x(1) = xc
      ifailh=1
      ifaila=1
      ifailg=1
      ifaill=1
c get interpolated values for heavy quark, heavy antiquark, gluon
c and sea
      call e02bbf(m+4,lamda,cheavy,x(1),xdhq,ifailh)
      call e02bbf(m+4,lamda,cantih,x(1),xdha,ifaila)
      call e02bbf(m+4,lamda,cgluon,x(1),xdgp,ifailg)
      call e02bbf(m+4,lamda,clight,x(1),xdup,ifaill)
      ifail = ifailh+ifailg+ifaill+ifaila
      if(hvqs.eq.'b') then
         xdbp = xdhq
         xdbbp = xdha
         xdcp = xdup
         xdcbp = xdup
      elseif(hvqs.eq.'c') then
c         xdbp = xdup
c         xdbbp = xdup
c the above seems inconsistent; bug corrected by PN 6/2001
         xdbp = 0
         xdbbp = 0
         xdcp = xdhq
         xdcbp = xdha
      else
         write(*,*) ' error: heavy quark type=',hvqs,' in fonfra'
         stop
      endif
      xdubp = xdup
      xddp = xdup
      xddbp = xdup
      xdsp = xdup
      xdtp = 0d0
      xdtbp = 0d0
c
      xcold = xc
      xdupold = xdup
      xdbpold = xdbp
      xdbbpold = xdbbp
      xdcpold = xdcp
      xdcbpold = xdcbp
      xdgpold = xdgp
 998  continue
      if(ifragmode.eq.1) then
         xdup=0
         xdubp=0
         xddp=0
         xddbp=0
         xdsp=0
         xdgp=0
         xdcbp=0
         xdbbp=0
         xdtp=0
         xdtbp=0
         if(hvqs.eq.'b') then
            xdcp=0
         elseif(hvqs.eq.'c') then
            xdbp=0
         else
            write(*,*) ' unimplemented flavour ',hvqs
            stop
         endif
      endif
      end

      subroutine fonfra_spline(qp2)
      implicit none
      real * 8 qp2
      integer m1,m2,m3,m,lck,lwrk
      parameter (m1=0,m2=50,m3=50)
      parameter (m=m1+m2+m3,lck=m+4,lwrk=6*m+16)
      real*8 yheavy(m), ygluon(m),ylight(m),ybott(m),ycharm(m)
      real*8 yantibot(m),yantich(m),yantih(m)
      character*1 hvqs
      common/hvqtype/hvqs
      real * 8 cheavy(lck),cgluon(lck),clight(lck),cantih(lck),
     #     lamda(lck)
      common/spline/cheavy,cgluon,clight,cantih,lamda
      real * 8 wrk(lwrk), x(m)
      logical verbose
      common/chat/verbose
      integer ifailh,ifaila,ifailg,ifaill,ifail
      real * 8 xdup,xdubp,xddp,xddbp,xdsp
     #        ,xdcp,xdcbp,xdbp,xdbbp,xdtp,xdtbp,xdgp
      integer i
      real * 8 xl,xu,yl,yu,y
      external e01baf
      data xl/.01d0/,xu/.9999d0/
c     data xl/.2d0/,xm/.97d0/,xu/.9999d0/
c     data xl/.05d0/,xm/.97d0/,xu/.99d0/
      save xl,xu
      yl=log((1-xl)/xl)
      yu=log((1-xu)/xu)
      do i = 1,m
         y=yl+(i-1)*(yu-yl)/(m-1)
         x(i) = 1/(1+exp(y))
         call fonfra_true(x(i),qp2,xdup,xdubp,xddp,xddbp,xdsp
     #        ,xdcp,xdcbp,xdbp,xdbbp,xdtp,xdtbp,xdgp)
         ybott(i) = xdbp
         yantibot(i) = xdbbp
         ycharm(i) = xdcp
         yantich(i) = xdcbp
         ygluon(i) = xdgp
         ylight(i) = xdup
      enddo
      if(hvqs.eq.'b') then
         do i =1,m
            yheavy(i) = ybott(i)
            yantih(i) = yantibot(i)
         enddo
      else
         do i =1,m
            yheavy(i) = ycharm(i)
            yantih(i) = yantich(i)
         enddo
      endif
      ifailh=-1
      ifailg=-1
      ifaill=-1
      ifaila=-1
      call e01baf(m,x,yheavy,lamda,cheavy,lck,wrk,lwrk,ifailh)
      call e01baf(m,x,yantih,lamda,cantih,lck,wrk,lwrk,ifaila)
      call e01baf(m,x,ygluon,lamda,cgluon,lck,wrk,lwrk,ifailg)
      call e01baf(m,x,ylight,lamda,clight,lck,wrk,lwrk,ifaill)
      ifail =  ifailh+ifailg+ifaill+ifaila
      if(verbose) then
         write(*,*)' End of E01BAF, ifail = ',ifail,', Q = ',sqrt(qp2)
      endif
c plot the fragmentation functions
c      call plotfra(qp2,m,lamda,cheavy,cantih,cgluon,clight)
      end
      subroutine plotfra(qp2,m,lamda,cheavy,cantih,cgluon,clight)
      implicit none
      integer m,ifail,j,k
      real * 8 lamda(*),cheavy(*),cantih(*),cgluon(*),clight(*)
      character * 8 cqp2
      real * 8 xx,xval,qp2
      write(cqp2,'(i8)')nint(sqrt(qp2))
      do j=1,8
         if(cqp2(j:j).ne.' ') then
            k=j
            goto 10
         endif
      enddo
 10   continue
      open(unit=11,file='fragpl'//cqp2(k:)//'.top',status='unknown')
      write(11,*)' set order x dummy y'
      write(11,*)' ( q= ',sqrt(qp2)
      write(11,*)' ( D_Q '
      do xx=0.001d0,0.999d0,0.001d0
         call e02bbf(m+4,lamda,cheavy,xx,xval,ifail)
         write(11,*) xx,1-xx,xval
      enddo
      write(11,*) ' join'
      write(11,*)' ( D_Q_bar '
      do xx=0.001d0,0.999d0,0.001d0
         call e02bbf(m+4,lamda,cantih,xx,xval,ifail)
         write(11,*) xx,1-xx,xval
      enddo
      write(11,*) ' join'
      write(11,*)' ( D_g '
      do xx=0.001d0,0.999d0,0.001d0
         call e02bbf(m+4,lamda,cgluon,xx,xval,ifail)
         write(11,*) xx,1-xx,xval
      enddo
      write(11,*) ' join'
      write(11,*)' ( D_light '
      do xx=0.001d0,0.999d0,0.001d0
         call e02bbf(m+4,lamda,clight,xx,xval,ifail)
         write(11,*) xx,1-xx,xval
      enddo
      write(11,*) ' join'
      close(11)
      end

      subroutine fonfra_true(xc,qp2,xdup,xdubp,xddp,xddbp,xdsp
     #     ,xdcp,xdcbp,xdbp,xdbbp,xdtp,xdtbp,xdgp)
      implicit none
      real * 8 xc,qp2,xdup,xdubp,xddp,xddbp,xdsp,xdsbp
     #     ,xdcp,xdcbp,xdbp,xdbbp,xdtp,xdtbp,xdgp
c     This subroutine returns x*fragmentation functions at the scale
c     qp2, x=xc.
c     xdup is for u, xdubp is for ubarr, xddp is
c     for d, xddbp is for dbarr, xdsp is for s, xdcp is for c, xdbp
c     is for b, xdbbp is for bbarr, xdtp is for t, xdtbp is for tbarr,
c     and xdgp is for g.
      integer n,maxfn,maxnp,maxnf
      parameter (n=100,maxnf=10,maxnp=2*maxnf+1)
      real * 8 res(maxnp)
      character*1 hvqs
      common/hvqtype/hvqs
      real * 8 step,ermx
      integer nmax
      common/integr/step,ermx,nmax
      real * 8 alam5qcd,xmb,xmc
      common/lamqcd/alam5qcd
      common/hvqmass/xmb,xmc
      include 'xsect.h'
      real * 8 cmu0
      common/cmu0/cmu0
      integer nfl
      common/flnumb/nfl
      integer leadfl,nf,lcoef
      real * 8 q0,qf,sclf,xlam,xmq
c     Integration parameters
      nmax = 7                  ! enter initial # of step of integration
      step = 1.d0/nmax
      ermx = 1d-3               ! maximum error
      leadfl = iloopfr		! 1 for leading, 2 for next-to-leading
******************number of flavours ********************************
      nf=abs(nfl)
***********************************************************************
      if(hvqs.eq.'b') then
         xmq = xmb
      else
         if(hvqs.eq.'c') then
            xmq = xmc
         else
            print*,'Non-implemented kind of heavy flavour'
            stop
         endif
      endif
c.... 'enter Q initial, Q final, and mu/Q'
      q0 = cmu0*xmq
      qf = sqrt(qp2)
      sclf = 1d0
      xlam = alam5qcd           ! 'enter lambda_5'
c     c      xlam = .2d0		! 'enter lambda_5'
c.....computes evolution
      call frag(xc,q0,qf,xlam,res,sclf,nf,leadfl,lcoef)
      xdup = xc*res(6)
c     xdup = 0d0
      xdubp = xc*res(7)
c     xdubp = 0d0
      xddp = xdup
      xddbp = xdubp
      xdsp = xdup
      xdsbp = xdubp
      xdcp = xc*res(4)
c     xdcp = 0d0
      xdcbp = xc*res(5)
c     xdcbp = 0d0
      xdbp = xc*res(2)
      xdbbp = xc*res(3)
      if(abs(nfl).eq.4) then
         xdbp = 0d0
         xdbbp = 0d0
      endif
      xdtp = 0d0
      xdtbp = 0d0
      xdgp = xc*res(1)
c     xdgp = 0d0
c     write(*,'(5(1x,e14.6))') xc,res(4),res(5),res(1),res(2)
      end




c Notation used in evolution program:
c Splitting functions:
c
c | Pgg         Pgq/(2 n_f) | |g|
c |                         | | |
c | 2 n_f Pqg   Pqq         | |s|
c
c s = q_1+qbar_1+q_2+qbar_2+...
c
c Thus Pgq has to be interpreted as a gluon splitting into a quark
c (not like: finding a gluon in the quark, which corresponds to a quark
c splitting into a gluon (DIS)
c
c In leading order:
c
c Pgq = 2 nf Tf (x^2+(1-x)^2)
c Pqg =      Cf (1+(1-x)^2)/x
c In next-to-leading order, the splitting kernels correspond to the paper
c by Furmanski and Petronzio. The above notation coorespond to the
c notation of that paper, except that they use F as upper component
c and g as lower component.
c The corresponding notation for structure functions would be
c |Pqq  Pgq||s|
c |Pqg  Pgg||g|
c where s has the same meaning as above
c
      subroutine frag(x,q0,qf,xlam,res,sclf,nf,leadfl,lcoef)
c      implicit real * 8 (a-h,o-z)
      implicit none
c arguments
      real * 8 x,q0,qf,xlam,res(*),sclf
      integer nf,leadfl,lcoef
c
      real * 8 xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2
      integer nnf
      common/bndary/xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2,nnf
      real * 8 step,ermx
      integer nmax
      common/integr/step,ermx,nmax
c static local
      real * 8 x0,x1
c local
      real * 8 err
      integer np,nit
c functions
      real * 8 alfas_p
      external evol1
      data x0,x1/0.d0,1.d0/
      nnf = nf
      xxsclf2 = sclf**2
      xxqsq0 = q0*q0
      xxqsqf = qf*qf*xxsclf2
      xxalf0 = alfas_p(xxqsq0,xlam,nnf)
      xxalff = alfas_p(xxqsqf,xlam,nnf)
c-------------np is the number of parton types
      np = 2*nf+2
c.....air gives the evolved fragm. funct. by performing numerically
c.....an inverse mellin transform
      call air(evol1,x0,x1,step,ermx,err,nmax,nit,res,np,
     #     x,leadfl,lcoef)
c     write(*,*)'error=',err,'for ',nit,'steps'
      call parton(res,nf)
      end

      subroutine parton(fpart,nf)
c------------------------------------------------------------------------------
c     Converts the parton densities given as:
c                                      (q = quark, a = antiquark, s = singlet
c     fpart(1) = g
c     (2) = s                          ( singlet = q_1 + a_1 + q_2 + a_2 + ...
c     (3) = (q_1 + a_1)/2 - s/(2*nf)
c     (4) = (q_1 - a_1)/2
c     ...
c     (2*nf+1) = (q_nf + a_nf)/2 - s/(2*nf)
c     (2*nf+2) = (q_nf - a_nf)/2
c
      implicit none
c argumens
      real * 8 fpart(*)
      integer nf
c
      integer j,last
      real * 8 s,s1,q,a
      last = 2 * nf + 1
      s    = fpart(2)
      s1   = s/(2*nf)
      do j=3,2*nf+1,2
         q = fpart(j) + fpart(j+1) + s1
         a = q - 2*fpart(j+1)
         fpart(j-1) = q
         fpart( j ) = a
      enddo
      return
      end

c.....the subroutine strfnc contains the initial conditions for the
c.....fragmentation functions of quarks and gluon, given in the
c.....moment space (refer to Mele-Nason, Nucl. Phys. B361 (1991) 626 )
      subroutine strfnc(cn,csn,cfl)
c arguments
      implicit none
      complex * 16 cn,csn(2),cfl(10)
c
      real * 8 xcf,xtf,eulergamma,logscale
      complex*16 constants_ini,coverlap_ini,delta,
     #           delta_ini
      real*8 b0ini,xnmaxi
      parameter (xcf=4.d0/3.d0,xtf=0.5d0)
      real * 8 xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2
      integer nnf
      common/bndary/xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2,nnf
      real*8 alfa,beta
      integer npsm
      common/alfabeta/alfa,beta,npsm
      real * 8 alam5qcd
      common/lamqcd/alam5qcd
      include 'xsect.h'
      integer isuda
      common/sudakov/isuda
c----------------------------------------------------
c     csn(1) is the gluon,
c     csn(2) is the singlet quark, = sum(q+q_bar)
c     cfl(1) (q1+q1_bar)/2-csn(2)/(2*nf)
c     cfl(2) (q1-q1_bar)/2
c     cfl(3) (q2+q2_bar)/2-csn(2)/(2*nf)
c     cfl(4) (q2-q2_bar)/2
c     ......
c     cfl(2nf-1)= (qnf+qnf_bar)/2-csn(2)/(2*nf)
c     cfl(2nf)  = (qnf-qnf_bar)/2
c
      real * 8 xmb,xmc
      common/hvqmass/xmb,xmc
      character*1 hvqs
      common/hvqtype/hvqs
      character * 2 frscheme
      common/frschemec/frscheme
c local
      real*8 qq,b0,xpi,xmq,xlam,xl2,xlgq0,xlgm0,xc,xm2
      integer nf,j
      complex * 16 one,cpsi0u,cpsi1u,s1,s2,cpqq0,cdin1,dnp,clgn,clgq1,
     # clgm1,clgm2,quarkb,aqb,quarkc,aqc,quarkl
c functions
      real*8 alfas_p
      complex * 16 cpsi0,cpsi1,cdnp,cpeterson,csudev,csudin
      data one/(1.d0,0.d0)/
      xpi = 4.d0 * datan(1.d0)
      nf = nnf
      if(hvqs.eq.'b') then
         xmq = xmb
      else
         xmq = xmc
      endif
      
c-----------------------------------
      cpsi0u = cpsi0(one)
      cpsi1u = cpsi1(one)
      s1 =  cpsi0(cn+1) - cpsi0u
      s2 = -cpsi1(cn+1) + cpsi1u
      cpqq0 = 1.5d0 + 1/cn/(cn+1) - 2.d0*s1
      cdin1 =  -2*(s1**2-1/(cn*(cn+1))*s1
     #     +1/(cn+1)**2+s2)
     #     +2-1/(cn*(cn+1))+2*s1
c.....non perturbative smearing
      if(npsm.eq.1) then
c     dnp = anorm*cdnp(alfa,beta,cn)
         dnp = cdnp(alfa,beta,cn)
c      elseif(npsm.eq.2) then
c     dnp = anorm*cpeterson(alfa,cn)
c         dnp = cpeterson(alfa,cn)
      else
         dnp = 1
      endif
c     c      print*,dnp


c.....Sudakov terms for the initial condition
      if (isuda.eq.1) then
c---------------------------------------------------
c     lambda 1 loop e' circa lambda_MS_bar/2
         xlam = alam5qcd
         xl2 = (xlam/2)**2
         qq = xxqsqf
         b0 = (33.-2.*nf)/xpi/12.d0
         xc = xcf/(xpi*b0)
         xm2 = xmq**2
         clgn  = log(cn)
         xlgq0 = log(qq/xl2)
         clgq1 = xlgq0-clgn
         xlgm0 = log(xm2/xl2)
         clgm1 = xlgm0-clgn
         clgm2 = clgm1-clgn
c-------------------------------------------------------
c     Sudakov form factor nell'evoluzione da m^2 a q^2.
c     I termini di ordine alfa^0 e alfa^1 vengono sottratti
c     dall'esponente, perche' appaiono di gia' nei
c     kernel di evoluzione leading e NL.
c
c...  ordine alpha. log(N) gia' sottratto
         csudev = exp(xc*(
     #        clgq1 * log(clgq1/xlgq0)
     #        - clgm1 * log(clgm1/xlgm0)
     #        ))
         if(iloopfr.eq.2) then
            csudev = csudev*exp(xc*( -.5d0*clgn**2/xlgq0
     #           +.5d0*clgn**2/xlgm0))
         endif
         csudin = exp(xc*(clgm1*log(clgm1)
     #        -.5d0*clgm2*log(clgm2) - .5d0*xlgm0*log(xlgm0)
     #        ))
         if(iloopfr.eq.2) then
            csudin = csudin + xc*.5d0*clgn**2/xlgm0
         endif
	 
c....new sudakov implementation. Modern version	 
      elseif(isuda.ge.9) then
        eulergamma=0.5772156649d0
	if(isuda.eq.12) then
c	  clgn = -log((1/cn + alam5qcd/xmq)/(1 + alam5qcd/xmq))
         b0ini = (33.-2.*(nf-1))/xpi/12.d0
         xnmaxi = 1.25*exp(-1/(2*xxalf0*b0ini))
         clgn = -log((1/cn + xnmaxi)/(1 + xnmaxi))
        else
	  clgn = log(cn)
	endif
c        logscale = log(xxqsq0/xmq/xmq)
        logscale = 0	! fix q0 = m
c...expansion of delta	                                                           
        coverlap_ini = 1.+ xcf/xpi*xxalf0*(-clgn**2 - clgn*(logscale
     #             +2*eulergamma - 1.))
c...constants from the initial condition
        constants_ini =  xcf/xpi*xxalf0*(-eulergamma**2 
     #         - xpi**2/6 + 7./4. - (eulergamma - 3./4.)*
     #           (logscale - 1.) )
c...resummation factor
c...we have decided to call the sudakov always with nlf flavours,
c...i.e. 3 for charm and 4 for bottom. Hence, we use here nf-1
c...This is now consistent with HVQF
        delta_ini = delta(clgn,xxalf0,nf-1)
      else
         csudev = 1.
         csudin = 1.
      endif
c-----------------------------------
c     gluon
      if(iloopfr.eq.2) then
         csn(1) = xxalf0/2d0/xpi*xtf*(1/cn - 2/(cn+1.d0) + 2/(cn+2.d0))*
     #        log(xxqsq0/xmq**2)
         csn(1) = csn(1)*dnp
      else
         csn(1) = 0
      endif
c     b
      if(hvqs.eq.'b') then
         if(iloopfr.eq.2) then
            if(frscheme.eq.'DL') then
               quarkb = 1.d0 + xxalf0/2d0/xpi*xcf*
     #              (log(xxqsq0/xmb**2)*cpqq0 )
            else
               quarkb = 1.d0 + xxalf0/2d0/xpi*xcf*
     #              (log(xxqsq0/xmb**2)*cpqq0 + cdin1 )
            endif
         elseif(iloopfr.eq.1) then
            quarkb=1
         else
            write(*,*) ' iloopfr=',iloopfr,' in strfnc'
            stop
         endif
c.....Sudakov
         if(isuda.eq.1) then
            quarkb = quarkb*csudev*csudin
         elseif(isuda.ge.9) then
            quarkb = delta_ini*exp(quarkb - coverlap_ini)
         endif	 
         quarkb = quarkb*dnp
         aqb = 0d0
      else
         quarkb = 0d0
         aqb = 0d0
      endif
c     c
      if(hvqs.eq.'c') then
         if(iloopfr.eq.2) then
            if(frscheme.eq.'DL') then
               quarkc = 1.d0+ xxalf0/2d0/xpi*xcf*
     #              (log(xxqsq0/xmc**2)*cpqq0 ) 
            else
               quarkc = 1.d0 + xxalf0/2d0/xpi*xcf*
     #              (log(xxqsq0/xmc**2)*cpqq0 + cdin1 )
            endif
         elseif(iloopfr.eq.1) then
            quarkc=1
         else
            write(*,*) ' iloopfr=',iloopfr,' in strfnc'
            stop
         endif
c     xlam = alam5qcd
c     xxal = alfas_p(xxqsqf/4,xlam,nnf)
c     print*,sqrt(xxqsqf/4), xlam, nnf,xxal
c     quarkc = 1.d0 + xxal/2d0/xpi*xcf*
c     #         (log(xxqsq0/xmc**2)*cpqq0 + cdin1 )
c.....Sudakov
         if(isuda.eq.1) then
            quarkc = quarkc*csudev*csudin
         elseif(isuda.ge.9) then
            quarkc = delta_ini*exp(quarkc - coverlap_ini)
         endif	 
         quarkc = quarkc*dnp
*********************
c     quarkc = dnp
*********************
c     aqc = quarkc
         aqc = 0d0
      else
         quarkc = 0d0
         aqc = 0d0
      endif
c     light
c     c      quarkl = xal*exp(-clgam(xbl+cn)+clgam(cn-1)+dlgl)
      quarkl = 0d0
c
c     singlet
      csn(2) = quarkb + aqb + quarkc + aqc + 6*quarkl
c-----
      s1 = csn(2)/(2*nf)
c--   cfl(1) = ( q_1 + a_1 )/2 - s1
      cfl(1) = (quarkb + aqb)/2d0 - s1
c--   cfl(2) = ( q_1 - a_1 )/2
      cfl(2) = (quarkb - aqb)/2d0
c--   cfl(3)
      cfl(3) = (quarkc + aqc)/2d0-s1
      cfl(4) = (quarkc - aqc)/2d0
      do j = 5, 2*nf - 1,2
         cfl( j ) = quarkl - s1
         cfl(j+1) = 0
      enddo
      return
      end


c...the Delta resummation factor
      complex*16 function delta(clogn,as,nf)
      implicit none
      real * 8 ca,cf,pi,tr
      parameter (ca=3d0,cf=4d0/3d0,pi=3.141592653589793d0,tr=0.5d0)
      real*8 b0
      complex * 16 clogn
      real * 8 as
      integer nf
      complex*16 lambda,g1,g2
      real*8 logscale,logmurmuf2
      integer isuda
      common/sudakov/isuda
      integer inf
      data inf/0/
      save inf,b0
      if(inf.eq.0.or.inf.ne.nf) then
         b0 = (11*ca-4*tr*nf)/12/pi
         inf=nf
      endif      

*********************************************************

c      logmurmuf2 = log(mur**2/muf**2)
c      logscale = log(muf**2/xm2)
      logmurmuf2 = 0
      logscale = 0

********************************************************
      
      lambda = b0*as*clogn
      
       if(isuda.eq.3) then
c...DLA
         delta = exp(-cf/pi*as*clogn**2)
       elseif(isuda.eq.4.or.isuda.eq.7) then
c....Leading terms
         delta = exp(clogn*g1(lambda,nf))
       else
c....full NLL expression
         delta = exp(clogn*g1(lambda,nf)
     #        + g2(lambda,logmurmuf2,logscale,nf))
       endif
       
c...no exp
c      delta = clogn*g1(lambda) + g2(lambda)

c....leading order only
c      delta = exp(clogn*g1(lambda)) 
 
c      print*,cn,delta,lambda,b0,as,clogn
      
c      delta = log(cn)*ca*as/pi
      
      
      end     

c.....the g1 function
      complex*16 function g1(lambda,nf)
      implicit none
      integer nf
      complex*16 lambda
      real * 8 ca,cf,pi,tr
      parameter (ca=3d0,cf=4d0/3d0,pi=3.141592653589793d0,tr=0.5d0)
      real*8 b0
      character*2 channel
      complex * 16 one
      integer inf
      data one/(1.d0,0.d0)/
      data inf/0/
      save inf,one,b0
      if(inf.eq.0.or.inf.ne.nf) then
         b0 = (11*ca-4*tr*nf)/12/pi
         inf=nf
      endif      
      
      g1 = - cf/2/pi/b0/lambda*(2*lambda+(1.d0-2*lambda)*
     #       log(1.d0-2*lambda))

      end
      
c.....the g2 function
      complex*16 function g2(lambda,logmurmuf2,logscale,nf)
      implicit none
      complex * 16 lambda
      real * 8 logmurmuf2,logscale
      integer nf
      real * 8 ca,cf,pi,tr,eulergamma
      parameter (ca=3d0,cf=4d0/3d0,pi=3.141592653589793d0,tr=0.5d0)
      parameter (eulergamma=0.5772156649d0)
      complex*16 o2l,log12l,log12l2l
      real*8 k,b0,b1
      integer inf
      data inf/0/
      save inf,b0,b1,k
      if(inf.eq.0.or.inf.ne.nf) then
         b0 = (11*ca-4*tr*nf)/12/pi
         b1 = (17*ca**2 - 10*ca*tr*nf - 6*cf*tr*nf)/24/pi**2
         k  = (67./18. - pi**2/6.)*ca - 10./9.*tr*nf
         inf=nf
      endif      
c      g2 = -eulergamma*2*log(1.-2*lambda) + b1/b0**2*(2*lambda + 
c     #        log(1.-2*lambda) + 0.5*log(1.-2*lambda)**2)
c     #	   - k/2/pi/b0*(2*lambda + log(1.-2*lambda))
c     #	   - log(1.-2*lambda)*logscale
c     #     - (log(1.-2*lambda) + 2*lambda)*log(mur**2/muf**2)

      o2l = 1.d0 - 2*lambda
      log12l = log(o2l)
      log12l2l = log12l + 2*lambda
c now given in argument list
c      logmurmuf2 = log(mur**2/muf**2)
      
      g2 = +eulergamma*log12l 
     #     - b1/2/b0**2*(log12l2l + 0.5*log12l**2)
     #	   + k/4/pi/b0*log12l2l
     #	   + 0.5*log12l*logscale 
     #     + 0.5*log12l2l*logmurmuf2

      g2 = cf/pi/b0*g2
      
      g2 = g2 - cf/2/pi/b0*log(1.-2*lambda)      
      
      end
      
      







      function cdnp(alfa,beta,cn)
      implicit none
      real * 8 alfa,beta,alfanp,betanp,dlogam
      complex*16 cn,cdnp,wlogam,logcdnp,ccn
      ccn = cn
      alfanp = alfa
      betanp = beta
      logcdnp =  dlogam(alfanp+betanp+2)+wlogam(betanp+ccn)
     #     -wlogam(alfanp+betanp+1+ccn)-dlogam(betanp+1)
c.....D_NP non normalizzata a 1
c     logcdnp = wlogam(betanp+ccn)+dlogam(alfanp+1)-
c     #          wlogam(betanp+ccn+alfanp+1)
      cdnp = exp(logcdnp)
      end

c     complex*16 function cdnp(alfa,beta,cn)
c     implicit real*8(a-h,o-z)
c     complex*16 cn,ca,cx,hypergeom16
c     data xpi/3.141592653589793/
c     ca = (6d0,0d0)
c     cx = (.5d0,0d0)
c     eps = 1d-4
c     cdnp = 1d0/64d0*(
c     #       16d0/(cn+1)*hypergeom16(ca,cn+1,cn+2,cx,eps,1000)
c     #      -64d0/(cn+2)*hypergeom16(ca,cn+2,cn+3,cx,eps,1000)
c     #      +152d0/(cn+3)*hypergeom16(ca,cn+3,cn+4,cx,eps,1000)
c     #      -208d0/(cn+4)*hypergeom16(ca,cn+4,cn+5,cx,eps,1000)
c     #      +141d0/(cn+5)*hypergeom16(ca,cn+5,cn+6,cx,eps,1000)
c     #      -42d0/(cn+6)*hypergeom16(ca,cn+6,cn+7,cx,eps,1000)
c     #      +5d0/(cn+7)*hypergeom16(ca,cn+7,cn+8,cx,eps,1000))
c     as = 0.23d0
c     r02 = 0.82d0
c     ampsi = 3.097d0
c     cdnp = 64d0/27d0/xpi*as**2*r02/ampsi**3*cdnp
c     end

*       complex*16 function cpeterson(eps,cn)
*       implicit complex*16 (c)
*       implicit real*8 (a-b,d-h,o-z)
*       complex*16 hypergeom
*       common/xvalue/xx
*       data ci/(0.d0,1.d0)/
*       data ini /0/
*       save ini,anorm,epsold
*       cone = (1d0,0d0)
*       ctwo = (2d0,0d0)
*       cthree = (3d0,0d0)
*       cnp4 = cn + 4.
*       if(ini.eq.0.or.eps.ne.epsold) then
*          root = sqrt(4*eps-eps**2)
*          anorm = (eps**2 - 6*eps + 4.d0)/(4.d0-eps)/root*
*      #        (atan(eps/root) + atan((2.-eps)/root)) + 0.5d0*log(eps) +
*      #        1.d0/(4.d0-eps)
*          ini = 1
*          epsold = eps
*       endif
*       areal = eps/2.
*       aimag = eps/2.*sqrt(4./eps-1.)
*       cx1 = areal + ci*aimag
*       cocx1 = 1./cx1
*       cx2 = areal - ci*aimag
*       cocx2 = 1./cx2
*       acc = 1d-4
* c     nit = 160
*       nit = 150
*       ilim = ilimit(eps)
*       recn = dble(cn)
*       aimcn = imag(cn)
*       xhigh = 0.97
* c     if(abs(cn).le.ilimit(eps)) then
*       if(recn.le.ilim.and.aimcn.le.ilim.and.xx.lt.xhigh) then
* c     if(abs(cn).le.60) then
*          cpeterson = 2/anorm/(cn+1)/(cn+2)/(cn+3)*(
*      #     hypergeom(ctwo,cthree,cnp4,cocx1,acc,nit)/cx1**2/(cx1-cx2)**2
*      #   + hypergeom(ctwo,cthree,cnp4,cocx2,acc,nit)/cx2**2/(cx1-cx2)**2
*      #     +2*hypergeom(cone,cthree,cnp4,cocx1,acc,nit)/cx1/(cx1-cx2)**3
*      #     -2*hypergeom(cone,cthree,cnp4,cocx2,acc,nit)/cx2/(cx1-cx2)**3
*      #        )
*       else
* c..   up to ninth order
*          cpeterson =   2/anorm*(-3628800. - 80640*cn + 10886400*eps +
*      #        954720*cn*eps +
*      -        48600*cn**2*eps + 1080*cn**3*eps - 7620480*eps**2 -
*      -    1623456*cn*eps**2 - 232200*cn**2*eps**2 - 21000*cn**3*eps**2 -
*      -        1080*cn**4*eps**2 - 24*cn**5*eps**2 + 1451520*eps**3 +
*      -    787824*cn*eps**3 + 278820*cn**2*eps**3 + 64144*cn**3*eps**3 +
*      -        9495*cn**4*eps**3 + 871*cn**5*eps**3 + 45*cn**6*eps**3 +
*      -        cn**7*eps**3)/
*      -        (cn*(1 + cn)*(2 + cn)*(3 + cn)*(4 + cn)*(5 + cn)*(6 + cn)*
*      -        (7 + cn)*(8 + cn)*(9 + cn)*eps**5)
*       endif
*       end
* 
*       function ilimit(eps)
*       parameter (npoints=14)
*       implicit real (a-h,o-z)
*       real*8 eps
* c.....the limit of validity in N of the hypergeometric evaluation of the mellin
* c.....transform of the Peterson happens to be function of epsilon. After this
* c.....value, the transform of the series expansion around 1 should be used.
* c.....Here follows the limiting value, for each value of eps. The limiting
* c.....values are given as (n,n), as function of n. The modulus is then taken.
*       dimension veceps(npoints), vecn(npoints)
*       data veceps / 0.005, 0.01, 0.015, 0.02, 0.03, 0.04,
*      #     0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
*      #     0.12, 0.15 /
*       data vecn  /  157., 110., 89., 78., 63., 55., 46.,
*      #     42., 37., 42., 36., 33., 30., 29. /
*       save epsold, ilimitold
*       epssin = sngl(eps)
*       if(epssin.eq.epsold) then
*          ilimit = ilimitold
*          return
*       endif
*       epsold = epssin
* c.... for a given value of eps, we now perform an interpolation:
*       ideg = 2
*       an = divdif(vecn,veceps,npoints,epssin,ideg)
*       ilimit = an
* c.... then we calculate the modulus of the complex number (an,an)
* c     ilimit = sqrt(2.*an**2)
*       ilimitold = ilimit
*       end
* 
*       complex*16 FUNCTION HYPERGEOM(A,B,C,X,EPS,NMAX)
* c.....evaluates the hypergeomemtric function 2F1 for |x|>1
*       implicit none
*       REAL*8 EPS,rerel,aimrel
*       integer i, nmax
*       complex*16 a,b,c,x,m,aa,nn,poch,di,cpart1,cpart2,csum1,csum2,clgam
*       complex*16 cpsi0,cterm,lgpoch
*       poch(aa,nn) = exp(clgam(aa+nn)-clgam(aa))
*       lgpoch(aa,nn) = clgam(aa+nn)-clgam(aa)
*       m = b-a
*       cpart1 = exp(clgam(c))*(-x)**(-b)/exp(clgam(b))/exp(clgam(c-a))
*       csum1 = 0.
*       do i=0,nmax
*          di = dcmplx(i)
*          cterm =  exp(lgpoch(a,m+di) + lgpoch(1-c+a,m+di)-clgam(di+1)-
*      #        clgam(m+di+1))*x**(-di)*(log(-x) + cpsi0(1+m+di)
*      #        + cpsi0(1+di) - cpsi0(b+di) - cpsi0(c-b-di))
*          csum1 = csum1 + cterm
*          rerel = dble(cterm)/dble(csum1)
*          aimrel = imag(cterm)/imag(csum1)
*          if(abs(rerel).lt.eps.and.abs(aimrel).lt.eps) then
*             goto 10
*          endif
*       enddo
*       WRITE(*,'(1x,"NMAX=",I10,1X,"reached in HYPERGEOM")')NMAX
*       print*,c
*       print*,cterm,csum1
*       print*,rerel,aimrel
*       print*,' '
*  10   continue
*       cpart2 = (-x)**(-a)*exp(clgam(c)-clgam(b))
*       csum2 = 0.
*       do i=0,aint(dble(m-1))
*          di = dcmplx(i)
*          csum2 = csum2 +  exp(clgam(m-di))*poch(a,di)/exp(clgam(di+1))/
*      #        exp(clgam(c-a-di))*x**(-di)
*       enddo
*       hypergeom = cpart1*csum1 + cpart2*csum2
*       RETURN
*       END
* 
*       complex*16 FUNCTION HYPERGEOM16(A,B,C,X,EPS,NMAX)
*       REAL*8 EPS
*       complex*16 a,b,c,x
* C--
* C     1992 MAR 25
* C-----------------------------------------------------------------------
*       complex*16 DI,ADD,cone
* c     IF (dreal(C).LE.0.D0) THEN
* c     WRITE(*,*)'non positive Re(C) in HYPERGEOM16'
* c     STOP
* c     ENDIF
*       cone = (1d0,0d0)
*       ADD = cone
*       HYPERGEOM16=ADD
*       DO 10 I=0,NMAX
*          DI=dcmplx(I)
*          ADD=ADD*(A+DI)*(B+DI)*X/((C+DI)*(cone+DI))
*          HYPERGEOM16=HYPERGEOM16+ADD
* c     print*,di,abs(add),hypergeom16
*  10      IF (abs(ADD).LE.EPS) RETURN
* c     WRITE(*,'(5HNMAX=,I10,1X,21Hreached in HYPERGEOM16)')NMAX
*          RETURN
*       END

c-------------------------------------------------
c     The function evol1 integrated between zero and one
c     returns the parton densities at x=xx in the array xxf,
c     with xxf(1) the gluon, xxf(2) the singlet quark content,
c     and xxf(3:2 nf-1) the nonsinglet densities. They are organized as
c     xxf(3) -> q-q_bar, xxf(4) = (q+q_bar)_i - (q+q_bar)_j
c     It calls:
c     1) evmat(cn, xxqsq, xxqsqi, xxlam2, nf, cpsn, cpfm,cpfp)
c     to get the evolution matrix for complex n (cn), from xxqsqi to
c     xxqsq
c     with given Lambda=xxlam2,nf as number of flavours, and returns
c     in cpsn(2,2) the singlet matrix evolution, in cpfm, cpfp the
c     evolution factors for non-singlet densities:
c     cpfm for q - q_bar
c     cpfp for (q+q_bar)_i - (q+q_bar)_j
c     2) strfn(cn,csn,cfl) to get the structure functions, where csn(2)
c     returns the gluon and singlet quark, cfl the nonsinglet:
c     cfl(1) = q + a ( non-singlet )
c     cfl(2) = q - a
c
      subroutine evol1(xxn,xxf,xx,leadfl,lcoef)
c      implicit complex * 16 (a-h,o-w), real * 8 (x-z)
      implicit none
c      arguments
      real * 8 xxn,xxf(*),xx
      integer leadfl,lcoef
c
      real * 8 xpi
      parameter (xpi=3.141592653589793d0)
c
      real * 8 xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2
      integer nnf
      common/bndary/xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2,nnf
      real * 8 x
      common/xvalue/x
c
      complex * 16  cpsn(2,2),csn(2),cfl(20),coefq,coefg,cpfm,cpfp
c for testing
      complex * 16  cpsnxx(2,2),cpfmxx,cpfpxx
c
      complex * 16 ci,cn,cjac,cttt
      integer k1,k2,init,nf,j
      real * 8 xn,xomx,xlgx,z,xjac
      data ci/(0.d0,1.d0)/
      data k1,k2/2,1/,init/0/
      xn = xxn
      x  = xx
      nf = nnf
      if(nf.gt.10) then
         write(*,*)'too many light flavours in evol1, nf=', nf
         stop
      endif
      xomx = 1-x
      xlgx = - log(x)
c-------------------------------------------------------------------------
c     cn = 1+1/log(x) - z(1 - %i), z=(1-t)/t/log(x)
c
c     - first set of variable change
      z = -2.0d0*log(xn)/xlgx
      xjac = 2.0d0/xn/xlgx
c     - further variable change
c     - first choice ( im part goes to - infty as fast as real part)
c     cn = 3/(3+xlgx) + 1/xomx - z*(1-ci)
c...  for the Peterson
c     cn = 3/(3+xlgx) + 1/xomx**.25 - z*(1-ci)
c     cjac = - xjac * (1-ci)
c     - second choice ( im part goes to - infty slower than real part)
      cn = 3/(3+xlgx) + 1/xomx - z*(1-ci) - z**2 * xomx
      cjac = - xjac * ((1-ci) + 2*z*xomx)
c     - divide by %pi %i for inverse laplace
      cttt = cjac * exp(cn*xlgx)/xpi/ci
c     - get structure functions.
      call strfnc(cn, csn, cfl)
c     - singlet evolution matrix cpsn(2,2)
      call evmat(cn, leadfl, cpsn, cpfm,cpfp)
c Begin Debug
c      call oevmat(cn, leadfl, cpsnxx, cpfmxx, cpfpxx)
c      write(*,*) ' evmat test, n=',cn
c      write(*,*) cpfmxx/cpfm
c      write(*,*) cpfpxx/cpfp
c      write(*,*) cpsnxx(1,1)/cpsn(1,1)
c      if(cpsn(1,2).ne.0)
c     # write(*,*) cpsnxx(1,2)/cpsn(1,2)
c      if(cpsn(1,2).ne.0)
c     # write(*,*) cpsnxx(1,2)/cpsn(1,2)
c      write(*,*) cpsnxx(2,2)/cpsn(2,2)
c      read(*,*)
c End Debug
c     coefficient functions
      call coeffun(cn,nf,xxalff,coefq,coefg,xxsclf2,leadfl,lcoef)
      xxf(1) = coefg*( cttt *
     2     ( cpsn(1, 1) * csn(1) +        cpsn(1, 2) * csn(2) ) )
      xxf(2) = coefq*( cttt *
     2     ( cpsn(2, 1) * csn(1) +        cpsn(2, 2) * csn(2) ) )
c     - non-singlet evolution factor cpfl
      do j = 1,2*nf,2
c     - q + q_bar
         xxf(j+2) = coefq*(cttt*cpfp*cfl( j ))
c     - q - q_bar
         xxf(j+3) = coefq*(cttt*cpfm*cfl(j+1))
      enddo
      return
      end

      subroutine fragmell(cnn,leadfl,xxf)
c returns the fragmentation function cnn moment
      implicit complex * 16 (a-h,o-w), real * 8 (x-z)
      dimension cpsn(2,2), csn(2), cfl(20)
      dimension xxf(*)
      common/bndary/xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2,nnf
      data ci/(0.d0,1.d0)/
      data k1,k2/2,1/,init/0/
      data xpi/3.141592653589793d0/
      nf = nnf
      cn = cnn
      lcoef=0
c     print*,'cn = ',cn
c     - get structure functions.
      call strfnc(cn, csn, cfl)
c     - singlet evolution matrix cpsn(2,2)
      call evmat(cn, leadfl, cpsn, cpfm,cpfp)
c     coefficient functions
      call coeffun(cn,nf,xxalff,coefq,coefg,xxsclf2,leadfl,lcoef)
      xxf(1) = coefg*
     2     ( cpsn(1, 1) * csn(1) +        cpsn(1, 2) * csn(2) )
      xxf(2) = coefq*
     2     ( cpsn(2, 1) * csn(1) +        cpsn(2, 2) * csn(2) )
c     - non-singlet evolution factor cpfl
      do j = 1,2*nf,2
c     - q + q_bar
         xxf(j+2) = cpfp*cfl( j )
c     - q - q_bar
         xxf(j+3) = cpfm*cfl(j+1)
      enddo
      return
      end


      subroutine evmat(ccn,leadfl,cmsn,cfm,cfp)
c      implicit complex * 16 (a-h,o-w), real * 8 (x-z)
      implicit none
c arguments
      complex * 16 ccn,cmsn(2,2),cfm,cfp
      integer leadfl
c
      real * 8 xpf,xpi,xcf
      parameter (xpf=0.5d0,xcf=4.d0/3.d0,xpi=3.141592653589793d0)
      complex * 16 one
      parameter (one=(1.d0,0.d0))
      real * 8 xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2
      integer nnf
      common/bndary/xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2,nnf
      real*8 assc
      common/asscale/assc
      character * 2 frscheme
      common/frschemec/frscheme
      include 'xsect.h'
c local
      complex * 16 cn,s,s1,s2,d,cnorm,cnormt,exxp,exxm,exfl,
     #        cdin1,cspfm1,cspfp1,cspfl,clamp,clamm,
     #        pp(2,2),pm(2,2),cspsn(2,2),cspsn1(2,2),can(2,2),
     #        cp1p(2,2),cp1m(2,2),cm1p(2,2),cm1m(2,2)
      integer j,k,l,m
      real * 8 xao2p0,xao2p1,xb0,xb1,xt
      integer nf
c functions
      complex * 16 cpsi0,cpsi1
      cn = ccn
      xao2p0 = xxalf0 /2/xpi
      xao2p1 = xxalff /2/xpi
      nf = nnf
c--   b0 (Furmansky e Petronzio notation)
      xb0 = 11 -2*nf/3.d0
      xb1 = 102 - 38*nf/3.d0
c--   find xt (evolution parameter)
      xt = 2/xb0*log(xao2p0/xao2p1)
c.....!!!!!
      xt = xt/assc
c.... !!!!!
c--   Call the splitting function (lowest order)
      call spmat (cn,nf,cspsn ,cspfl )
******************MIO *********************
c     cspsn(1,1) = 0d0
c     cspsn(2,1) = 0d0
c     cspsn(1,2) = 0d0
ccc   cspsn(2,2) = 0d0
***************************************
c--   In the timelike case, correct for flavour factors
      cspsn(1,2)=cspsn(1,2)/(2*nf)
      cspsn(2,1)=cspsn(2,1)*(2*nf)
      if(leadfl.eq.1) then
c--   Find trace part of cspsn
         s = ( cspsn(1,1) + cspsn(2,2) )/2
         d = ( cspsn(1,1) - cspsn(2,2) )/2
c--   Find norm of traceless part
         cnorm = sqrt( d**2 + cspsn(1,2)*cspsn(2,1) )
         cnormt = cnorm * 2
c--   Find projection operators on positive and negative cnorm
         d = d/cnormt
         pp(1,1) = xpf + d
         pp(2,2) = xpf - d
         pp(1,2) = cspsn(1,2)/cnormt
         pp(2,1) = cspsn(2,1)/cnormt
         pm(1,1) = pp(2,2)
         pm(2,2) = pp(1,1)
         pm(1,2) = - pp(1,2)
         pm(2,1) = - pp(2,1)
c--   Find evolution factors
         clamp = s+cnorm
         clamm = s-cnorm
         exxp = exp( xt*clamp )
         exxm = exp( xt*clamm )
         exfl = exp( xt*cspfl )
c--   leading terms
         do j=1,2
            do k=1,2
               cmsn(j,k) = pp(j,k)*exxp + pm(j,k)*exxm
            enddo
         enddo
         cfp = exfl
         cfm = exfl
         return
      endif
c--   find subleading terms
      call spmat1(cn,nf,cspsn1,cspfm1,cspfp1)
******************MIO *********************
c     cspsn1(1,1) = 0d0
c     cspsn1(2,1) = 0d0
c     cspsn1(1,2) = 0d0
ccc   cspsn1(2,2) = 0d0
***************************************
c--   In the timelike case, correct for flavour factors
      cspsn1(1,2)=cspsn1(1,2)/(2*nf)
      cspsn1(2,1)=cspsn1(2,1)*(2*nf)
c-- scheme change
      if(frscheme.eq.'DL') then
c correction to initial condition (Mele, Nason)
         s1 =  cpsi0(cn+1) - cpsi0(one)
         s2 = -cpsi1(cn+1) + cpsi1(one)
         cdin1 =  -2*(s1**2-1/(cn*(cn+1))*s1
     #     +1/(cn+1)**2+s2)
     #     +2-1/(cn*(cn+1))+2*s1
c correction matrix: 1+as/(2 pi)*can
c which sign? According to Nason-Webber we need a - here
c (1+as/(2 pi) can) X old D's = new D's
         can(1,1)=0
         can(1,2)=0
         can(2,1)=0
         can(2,2)=-cdin1*xcf
c modify next to leading splitting kernel accordingly
         cspfm1=cspfm1-xb0/2*can(2,2)
         cspfp1=cspfp1-xb0/2*can(2,2)
         do j=1,2
         do k=1,2
            do l=1,2
c commutator piece
               cspsn1(j,k)=cspsn1(j,k)
     #            + can(j,l)*cspsn(l,k) - cspsn(j,l)*can(l,k)
            enddo
c term from derivative of as
            cspsn1(j,k)=cspsn1(j,k)-xb0/2*can(j,k)
         enddo
         enddo
      endif
c -- end scheme change
      do j=1,2
         do k=1,2
            cspsn1(j,k)=cspsn1(j,k)-xb1/(2*xb0)*cspsn(j,k)
         enddo
      enddo
c Expand about different point when ialtev=1
      if(ialtev.eq.1) then
         do j=1,2
            do k=1,2
               cspsn(j,k)=cspsn(j,k)+xao2p0*cspsn1(j,k)
            enddo
         enddo
      endif
c--   Find trace part of cspsn
      s = ( cspsn(1,1) + cspsn(2,2) )/2
      d = ( cspsn(1,1) - cspsn(2,2) )/2
c--   Find norm of traceless part
      cnorm = sqrt( d**2 + cspsn(1,2)*cspsn(2,1) )
      cnormt = cnorm * 2
c--   Find projection operators on positive and negative cnorm
      d = d/cnormt
      pp(1,1) = xpf + d
      pp(2,2) = xpf - d
      pp(1,2) = cspsn(1,2)/cnormt
      pp(2,1) = cspsn(2,1)/cnormt
      pm(1,1) = pp(2,2)
      pm(2,2) = pp(1,1)
      pm(1,2) = - pp(1,2)
      pm(2,1) = - pp(2,1)
c--   Find evolution factors
      clamp = s+cnorm
      clamm = s-cnorm
      exxp = exp( xt*clamp )
      exxm = exp( xt*clamm )
      exfl = exp( xt*cspfl )
c--   leading terms
      do j=1,2
         do k=1,2
            cmsn(j,k) = pp(j,k)*exxp + pm(j,k)*exxm
         enddo
      enddo
c-- compute cp1p,cm1m,cp1m,cm1p
      do j=1,2
         do k=1,2
            cp1p(j,k)=0
            cp1m(j,k)=0
            cm1p(j,k)=0
            cm1m(j,k)=0
            do l=1,2
               do m=1,2
                  cp1p(j,k)=pp(j,l)*cspsn1(l,m)*pp(m,k)+cp1p(j,k)
                  cp1m(j,k)=pp(j,l)*cspsn1(l,m)*pm(m,k)+cp1m(j,k)
                  cm1p(j,k)=pm(j,l)*cspsn1(l,m)*pp(m,k)+cm1p(j,k)
                  cm1m(j,k)=pm(j,l)*cspsn1(l,m)*pm(m,k)+cm1m(j,k)
               enddo
            enddo
         enddo
      enddo
      do j=1,2
         do k=1,2
            cmsn(j,k) = cmsn(j,k) + 
     #    (2/xb0*(xao2p0-xao2p1)-ialtev*xt*xao2p0)
     #           *( cp1p(j,k)*exxp + cm1m(j,k)*exxm ) +
     #    cp1m(j,k)*(
     #      (xao2p1*exxm-xao2p0*exxp)/(clamm-clamp-xb0/2) -
     #           ialtev*xao2p0*(exxm-exxp)/(clamm-clamp) ) +
     #    cm1p(j,k)*(
     #      (xao2p1*exxp-xao2p0*exxm)/(clamp-clamm-xb0/2) -
     #           ialtev*xao2p0*(exxp-exxm)/(clamp-clamm) )
         enddo
      enddo
      if(ialtev.eq.1) then
         cfp = exfl*exp(2/xb0*(xao2p0-xao2p1)*(cspfp1-xb1/xb0/2*cspfl))
         cfm = exfl*exp(2/xb0*(xao2p0-xao2p1)*(cspfm1-xb1/xb0/2*cspfl))
      else
         cfp = exfl*(1+(2/xb0*(xao2p0-xao2p1)*(cspfp1-xb1/xb0/2*cspfl)))
         cfm = exfl*(1+(2/xb0*(xao2p0-xao2p1)*(cspfm1-xb1/xb0/2*cspfl)))
      endif
      return
      end


      subroutine oevmat(ccn,leadfl,cmsn,cfm,cfp)
c      implicit complex * 16 (a-h,o-w), real * 8 (x-z)
      implicit none
c arguments
      complex * 16 ccn,cmsn(2,2),cfm,cfp
      integer leadfl
c
      real * 8 xpf,xpi,xcf
      parameter (xpf=0.5d0,xcf=4.d0/3.d0,xpi=3.141592653589793d0)
      complex * 16 one
      parameter (one=(1.d0,0.d0))
      complex * 16 pp(2,2),pm(2,2),cspsn(2,2),cspsn1(2,2)
      real * 8 xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2
      integer nnf
      common/bndary/xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2,nnf
      real*8 assc
      common/asscale/assc
      character * 2 frscheme
      common/frschemec/frscheme
c local
      complex * 16 cn,s,s1,s2,d,cnorm,cnormt,cfl,exxp,exxm,exfl,
     #        can(2,2),cdin1,cspfm1,cspfp1,cspfl,clamp,clamm
      integer j,k,l,m
      real * 8 xao2p0,xao2p1,xb0,xb1,xt
      integer nf
c functions
      complex * 16 cpsi0,cpsi1
      cn = ccn
      xao2p0 = xxalf0 /2/xpi
      xao2p1 = xxalff /2/xpi
      nf = nnf
c--   b0 (Furmansky e Petronzio notation)
      xb0 = 11 -2*nf/3.d0
      xb1 = 102 - 38*nf/3.d0
c--   find xt (evolution parameter)
      xt = 2/xb0*log(xao2p0/xao2p1)
c.....!!!!!
      xt = xt/assc
c.... !!!!!
c--   Call the splitting function (lowest order)
      call spmat (cn,nf,cspsn ,cspfl )
******************MIO *********************
c     cspsn(1,1) = 0d0
c     cspsn(2,1) = 0d0
c     cspsn(1,2) = 0d0
ccc   cspsn(2,2) = 0d0
***************************************
c--   In the timelike case, correct for flavour factors
      cspsn(1,2)=cspsn(1,2)/(2*nf)
      cspsn(2,1)=cspsn(2,1)*(2*nf)
c--   Find trace part of cspsn
      s = ( cspsn(1,1) + cspsn(2,2) )/2
      d = ( cspsn(1,1) - cspsn(2,2) )/2
c--   Find norm of traceless part
      cnorm = sqrt( d**2 + cspsn(1,2)*cspsn(2,1) )
      cnormt = cnorm * 2
c--   Find projection operators on positive and negative cnorm
      d = d/cnormt
      pp(1,1) = xpf + d
      pp(2,2) = xpf - d
      pp(1,2) = cspsn(1,2)/cnormt
      pp(2,1) = cspsn(2,1)/cnormt
      pm(1,1) = pp(2,2)
      pm(2,2) = pp(1,1)
      pm(1,2) = - pp(1,2)
      pm(2,1) = - pp(2,1)
c--   Find evolution factors
      clamp = s+cnorm
      clamm = s-cnorm
      exxp = exp( xt*clamp )
      exxm = exp( xt*clamm )
      exfl = exp( xt*cspfl )
c--   leading terms
      do j=1,2
         do k=1,2
            cmsn(j,k) = pp(j,k)*exxp + pm(j,k)*exxm
         enddo
      enddo
      cfp = exfl
      cfm = exfl
      if(leadfl.eq.1) return
c--   find subleading terms
      call spmat1(cn,nf,cspsn1,cspfm1,cspfp1)
******************MIO *********************
c     cspsn1(1,1) = 0d0
c     cspsn1(2,1) = 0d0
c     cspsn1(1,2) = 0d0
ccc   cspsn1(2,2) = 0d0
***************************************
c--   In the timelike case, correct for flavour factors
      cspsn1(1,2)=cspsn1(1,2)/(2*nf)
      cspsn1(2,1)=cspsn1(2,1)*(2*nf)
c-- scheme change
      if(frscheme.eq.'DL') then
c correction to initial condition (Mele, Nason)
         s1 =  cpsi0(cn+1) - cpsi0(one)
         s2 = -cpsi1(cn+1) + cpsi1(one)
         cdin1 =  -2*(s1**2-1/(cn*(cn+1))*s1
     #     +1/(cn+1)**2+s2)
     #     +2-1/(cn*(cn+1))+2*s1
c correction matrix: 1+as/(2 pi)*can
c which sign? According to Nason-Webber we need a - here
c (1+as/(2 pi) can) X old D's = new D's
         can(1,1)=0
         can(1,2)=0
         can(2,1)=0
         can(2,2)=-cdin1*xcf
c modify next to leading splitting kernel accordingly
         cspfm1=cspfm1-xb0/2*can(2,2)
         cspfp1=cspfp1-xb0/2*can(2,2)
         do j=1,2
         do k=1,2
            do l=1,2
c commutator piece
               cspsn1(j,k)=cspsn1(j,k)
     #            + can(j,l)*cspsn(l,k) - cspsn(j,l)*can(l,k)
            enddo
c term from derivative of as
            cspsn1(j,k)=cspsn1(j,k)-xb0/2*can(j,k)
         enddo
         enddo
      endif
c -- end scheme change
c--
      do j=1,2
         do k=1,2
c--   inner sum indeces
            do l=1,2
               do m=1,2
                  cmsn(j,k) = cmsn(j,k) +
     #             2/xb0*(xao2p0-xao2p1)*
     #             (pp(j,l)*cspsn1(l,m)*pp(m,k)*exxp +
     #             pm(j,l)*cspsn1(l,m)*pm(m,k)*exxm ) +
     #             pp(j,l)*cspsn1(l,m)*pm(m,k)*(xao2p1*exxm-xao2p0*exxp)
     #             /(clamm-clamp-xb0/2)   +
     #             pm(j,l)*cspsn1(l,m)*pp(m,k)*(xao2p1*exxp-xao2p0*exxm)
     #                 /(clamp-clamm-xb0/2)
               enddo
            enddo
            cmsn(j,k) = cmsn(j,k) +
     #           2/xb0*(xao2p0-xao2p1)*(-xb1/xb0/2)*
     #           (pp(j,k)*clamp*exxp+pm(j,k)*clamm*exxm)
         enddo
      enddo
      cfp = cfp + 2/xb0*(xao2p0-xao2p1)*exfl*
     #     (cspfp1-xb1/xb0/2*cspfl)
      cfm = cfm + 2/xb0*(xao2p0-xao2p1)*exfl*
     #     (cspfm1-xb1/xb0/2*cspfl)
      return
      end

c----------------------------------------------------------------
c     Splitting function in lowest order.
c     Return value: cspsn(2,2), cspfl
c
      subroutine spmat(ccn,nnf,cspsn,cspfl)
      implicit complex * 16 (a-h,o-w), real * 8 (x-z)
      dimension cspsn(2,2)
      data xca,xcf,xtf/3.d0,1.333333333333d0,.5d0/
      data xgame/5.772156649015329D-1/
      cn = ccn
      nf = nnf
      xb02pi = (11*xca-4*xtf*nf)/6.d0
      cnp1 = cn+1
      cnp2 = cnp1 + 1
      cnm1 = cn-1
      cs1n = cpsi0(cn) + xgame
      cs1np1= cs1n + 1/cn
      cs1np2= cs1np1 + 1/cnp1
      cncnp1 = cn*cnp1
      cttt  = (cncnp1 + 2)/cncnp1
      cspsn(1,1) = 2*xca*(-cs1np1+11.d0/12+1/cn/cnm1+1/cnp1/cnp2)
     #     -2*nf*xtf/3
      cspsn(1,2) = 2*nf*xtf*cttt/cnp2
      cspsn(2,1) = xcf*cttt/cnm1
      cspfl      = xcf *( 1.5d0 +1/cncnp1-2*cs1np1 )
      cspsn(2,2) = cspfl
      return
      end

c     Adaptive integration for a vector function
      subroutine air(f,xx0,xx1,sstep,eermax,eerr,nmax,nit,res,nd,
     #     arg3,iarg4,iarg5)
      parameter (nn=50)
      implicit real * 8 (a-h,o-z)
      dimension xg6(6), wg6(6), xg10(10), wg10(10)
      dimension res(nd), tmp6(nn),tmp10(nn),xv(nn)
      data init/0/, zero,one/0.d0,1.d0/,
     x     xredd,xredu,ss/.1d0,.9d0,.9d0/
      save init, zero,one,xredd,xredu,ss,xg6,wg6,xg10,wg10
      if(init.eq.0) then
         call gauleg(zero,one,xg6,wg6,6)
         call gauleg(zero,one,xg10,wg10,10)
         init = 1
      endif
      if(nd.gt.nn) then
         write(*,*)'air: too large a vector with ,'
     #        ,nd,'components (',nn,'max)'
         stop
      endif
      x0   = xx0
      x1   = xx1
      step  = sstep
      x01  = x1-x0
      if(x01.lt.0) then
         step   = -abs(step)
         isign  = -1
      else
         step = abs(step)
         isign  = 1
      endif
      step0 = step
      ermax = eermax
      nit = 0
c     - Start.
      call tozero(res,nd)
      err  = 0
 1    continue
      if(isign*x1.le.isign*x0) then
         eerr = err
         return
      endif
      if(isign*x01.le.isign*step) then
         step  = x01
      endif
      call tozero(tmp6,nd)
      call tozero(tmp10,nd)
      do j=1,6
         x = xg6(j)*step+x0
         call f(x,xv,arg3,iarg4,iarg5)
         call addmul(tmp6,xv,wg6(j)*step,nd)
      enddo
      do j=1,10
         x = xg10(j)*step + x0
         call f(x,xv,arg3,iarg4,iarg5)
         call addmul(tmp10,xv,wg10(j)*step,nd)
c     debug
c     write(15,'(10(1x,d10.4))') x,(xv(jjj),jjj=1,nd)
c     endebug
      enddo
      tmper = absdf(tmp10,tmp6,nd)
      if(tmper.lt.ermax) then
         call addvec(res,tmp10,nd)
         nit = nit+1
c     c         if(nit.gt.nmax) write(*,*) nit,'step'
         err = err+tmper
         if(step.eq.x01) then
            eerr = err
            return
         endif
         ermax = ermax - tmper
         x0   = x0 + step
         x01  = x1 - x0
         if(tmper.eq.0) then
            step = step0
         else
            step =
     #        isign*min(isign*step0,isign*step*(ermax/tmper)**xredu*ss)
         endif
         goto 1
      else
         step = step*(ermax/tmper)**xredd*ss
         goto 1
      endif
 131  write(*,*)'air: failed to converge at point',x1
      eerr = 1.0d30
      return
      end

      subroutine tozero(vec,n)
      implicit double precision (a-h,o-z)
      dimension vec(n)
      do j=1,n
         vec(j)=0
      enddo
      return
      end

      function absdf(vec1,vec2,n)
      implicit double precision (a-h,o-z)
      dimension vec1(n), vec2(n)
      ttt=0
      do j=1,n
         ttt = ttt + (vec1(j)-vec2(j))**2
      enddo
      absdf = sqrt(ttt)
      return
      end

      subroutine addmul(vec1,vec2,val,n)
      implicit double precision (a-h,o-z)
      dimension vec1(n),vec2(n)
      do j=1,n
         vec1(j) = vec1(j) + vec2(j)*val
      enddo
      return
      end

      subroutine addvec(vec1,vec2,n)
      implicit double precision (a-h,o-z)
      dimension vec1(n),vec2(n)
      do j=1,n
         vec1(j) = vec1(j) + vec2(j)
      enddo
      return
      end

c     gaussian integration for a vector function
      subroutine gau(f,xx0,xx1,sstep,eermax,eerr,res,nd,
     #     arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
      parameter (nn=100)
      implicit real * 8 (a-h,o-z)
      dimension xg(nn),wg(nn)
      dimension res(nd), tmp(nn)
      data init/0/, zero,one/0.d0,1.d0/,
     x     xredd,xredu,ss/.1d0,.9d0,.9d0/
      save init, zero,one,xredd,xredu,ss
      if(init.eq.0) then
         call gauleg(zero,one,xg,wg,nn)
         init = 1
      endif
      if(nd.gt.nn) then
         write(*,*)'gau: too large a vector with ,'
     #        ,nd,'components (',nn,'max)'
         stop
      endif
      x0   = xx0
      x1   = xx1
      x01  = x1 - x0
      call tozero(res,nd)
      do j = 1,nn
         x = xg(j)*x01 + x0
         call f(x,tmp,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)
c     debug
c     write(12,'(10(1x,d10.4))') x,(tmp(jjj),jjj=1,nd)
c     endebug
         call addmul(res,tmp,wg(j)*x01,nd)
      enddo
      return
      end

c-------------------------------------------------------------
c     Program to calculate alfa strong with nf flavours,
c     as a function of lambda with 5 flavors.
c     The value of alfa is matched at the thresholds q = mq.
c     When invoked with nf < 0 it chooses nf as the number of
c     flavors with mass less then q.
c
      function alfas_p(q2,xlam,inf)
      implicit real * 8 (a-h,o-z)
      common/hvqmass/xmb,xmc
      include 'xsect.h'
      common/asscale/assc
      data iflag/0/,pi/3.1415927d0/
c     data xmb/5.d0/,xmc/1.5d0/
      save iflag,b3,bp3,b4,bp4,b5,bp5,xlc,xlb,xllc,xllb,c45,c35
      if(iflag.eq.0) then
         b5  = (33-2*5)/pi/12
         bp5 = (153 - 19*5) / pi / 2 / (33 - 2*5)
         b4  = (33-2*4)/pi/12
         bp4 = (153 - 19*4) / pi / 2 / (33 - 2*4)
         b3  = (33-2*3)/pi/12
         bp3 = (153 - 19*3) / pi / 2 / (33 - 2*3)
c.....for leading order alfa_s
         if(iloopas.eq.1) then
            bp3 = 0d0
            bp4 = 0d0
            bp5 = 0d0
         endif
         xlc = 2 * log(xmc/xlam)
         xlb = 2 * log(xmb/xlam)
         xllc = log(xlc)
         xllb = log(xlb)
         c45  =  1/( 1/(b5 * xlb) - xllb*bp5/(b5 * xlb)**2 )
     #        - 1/( 1/(b4 * xlb) - xllb*bp4/(b4 * xlb)**2 )
         c35  =  1/( 1/(b4 * xlc) - xllc*bp4/(b4 * xlc)**2 )
     #        - 1/( 1/(b3 * xlc) - xllc*bp3/(b3 * xlc)**2 ) + c45
         iflag = 1
      endif
      q   = sqrt(q2)
      xlq = 2 * log( q/xlam )
      xllq = log( xlq )
      nf = inf
      if( nf .lt. 0) then
         if( q .gt. xmb ) then
            nf = 5
         elseif( q .gt. xmc ) then
            nf = 4
         else
            nf = 3
         endif
      endif
      if    ( nf .eq. 5 ) then
         alfas_p = 1/(b5 * xlq) -  bp5/(b5 * xlq)**2 * xllq
      elseif( nf .eq. 4 ) then
         alfas_p = 1/( 1/(1/(b4 * xlq) - bp4/(b4 * xlq)**2 *xllq)+c45 )
      elseif( nf .eq. 3 ) then
         alfas_p = 1/( 1/(1/(b3 * xlq) - bp3/(b3 * xlq)**2 * xllq)+c35)
      else
         print *,'error in alfa: unimplemented # of light flavours',nf
         stop
      endif
      alfas_p = alfas_p/assc
      return
      end

      function dlogam(x)
      implicit real * 8 (a-z)
      dlogam = dlgam(x)
      end

      function wlogam(x)
      implicit complex * 16 (a-h,o-z)
      wlogam = clgam(x)
      end

      function dlgam(xx)
c     real logarithm of gamma function
      implicit real * 8 (a-h,o-z)
      real * 8 cof(6),stp,gt,g,cst
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     #     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data gt,g/5.5d0,5.0d0/
      data cst/4.081061466d0/
      x = xx - 1
      xpgt = x + gt
      xmp5  = x + .5d0
      s0 = 1
      do 1 j=1,6
         x = x + 1
         tmp = cof(j)/x
         s0  = s0 + tmp
 1    continue
      r10 = log(s0)
      dlgam = xmp5*(log(xpgt)-1) + r10 - cst
      return
      end

      function clgam(xx)
c     complex logarithm of gamma function, cut along the negative
c     real axis. Good approximation for positive real part.
      implicit complex * 16 (a-h,o-z)
      real * 8 cof(6),stp,gt,g,cst
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     #     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data gt,g/5.5d0,5.0d0/
      data cst/4.081061466d0/
      x = xx - 1
      xpgt = x + gt
      xmp5  = x + .5d0
      s0 = 1
      do 1 j=1,6
         x = x + 1
         tmp = cof(j)/x
         s0  = s0 + tmp
 1    continue
      r10 = log(s0)
      clgam = xmp5*(log(xpgt)-1) + r10 - cst
      return
      end

c---------------------------------------------------------------
c     Poligamma functions
c     cpsi0(xx) = diff(log(gamma(x),x),
c     cpsi1(xx) = diff(log(gamma(x),x,2), etc.
c
      function cpsi0(xx)
      implicit complex * 16 (a-h,o-z)
      real * 8 cof(6),stp,gt,g
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     #     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data gt,g/5.5d0,5.0d0/
      x = xx - 1
      xpgt = x + gt
      s0 = 1
      s1 = 0
      do 1 j=1,6
         x = x + 1
         tmp = cof(j)/x
         s0  = s0 + tmp
         tmp = -tmp/x
         s1  = s1 + tmp
 1    continue
      r10 = s1/s0
      goxpgt = g/xpgt
      cpsi0 = log(xpgt) - goxpgt + r10
      return
      end

      function cpsi1(xx)
      implicit complex * 16 (a-h,o-z)
      real * 8 cof(6),stp,gt,g
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     #     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data gt,g/5.5d0,5.0d0/
      x = xx - 1
      xpgt = x + gt
      s0 = 1
      s1 = 0
      s2 = 0
      do 1 j=1,6
         x = x + 1
         tmp = cof(j)/x
         s0  = s0 + tmp
         tmp = -tmp/x
         s1  = s1 + tmp
         tmp = -2*tmp/x
         s2  = s2 + tmp
 1    continue
      r10 = s1/s0
      goxpgt = g/xpgt
      r20 = s2/s0
      r102 = r10*r10
      cpsi1 = (1+goxpgt)/xpgt - r102 + r20
      return
      end

      function cpsi2(xx)
      implicit complex * 16 (a-h,o-z)
      real * 8 cof(6),stp,gt,g
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     #     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data gt,g/5.5d0,5.0d0/
      x = xx - 1
      xpgt = x + gt
      s0 = 1
      s1 = 0
      s2 = 0
      s3 = 0
      do 1 j=1,6
         x = x + 1
         tmp = cof(j)/x
         s0  = s0 + tmp
         tmp = -tmp/x
         s1  = s1 + tmp
         tmp = -2*tmp/x
         s2  = s2 + tmp
         tmp = -3*tmp/x
         s3  = s3 + tmp
 1    continue
      r10 = s1/s0
      goxpgt = g/xpgt
c--   cpsi0 = log(xpgt) - goxpgt + r10
      r20 = s2/s0
      r102 = r10*r10
c--   cpsi1 = (1+goxpgt)/xpgt - r102 + r20
      r30 = s3/s0
      cpsi2 = -(1+2*goxpgt)/(xpgt*xpgt) + (2*r102-3*r20)*r10 + r30
      return
      end

      subroutine coeffun(cn,nf,xxalf,coefq,coefg,xxsclf2,leadfl,lcoef)
c     Dummy Coefficient functions for e+e- --> hadrons
c      implicit complex * 16 (a-h,o-z)
      implicit none
      real * 8 xxalf,xxsclf2
      complex * 16 cn,coefq,coefg
      integer nf,leadfl,lcoef
      coefq=1
      coefg=1
c For e^+e^-, in the Delta scheme (DL) the correction for the annihilation
c scheme is
c      character * 2 frscheme
c      common/frschemec/frscheme
c      complex * 16 s1,s2,cdin1,cpsi0,cpsi1
c      real * 8 xcf, xpi
c      parameter (xcf=4.d0/3.d0,xpi=3.141592653589793d0)
c      if(frscheme.eq.'DL') then
c         s1 =  cpsi0(cn+1) - cpsi0(one)
c         s2 = -cpsi1(cn+1) + cpsi1(one)
c         cdin1 =  -2*(s1**2-1/(cn*(cn+1))*s1
c     #     +1/(cn+1)**2+s2)
c     #     +2-1/(cn*(cn+1))+2*s1
c         coefq=coefq+xxalf/(2*pi)*xcf*cdin1
c      endif
      end

      subroutine ecoeffun(cn,nf,xxalf,coefq,coefg,xxsclf2,leadfl,lcoef)
c      Coefficient functions for e+e- --> hadrons
c      implicit complex * 16 (a-h,o-z)
      implicit none
      real * 8 xxalf,xxsclf2
      complex * 16 cn,coefq,coefg
      integer nf,leadfl,lcoef
      character * 2 frscheme
      common/frschemec/frscheme
      complex * 16 s1,s2,cdin1,cpsi0,cpsi1,anq,ang,one
      real * 8 xcf, xpi
      parameter (xcf=4.d0/3.d0,xpi=3.141592653589793d0,one=(1.d0,0.d0))
      coefq=1
      coefg=0
c For e^+e^-, in the Delta scheme (DL) the correction for the annihilation
c scheme is
      if(frscheme.eq.'DL') then
         s1 =  cpsi0(cn+1) - cpsi0(one)
         s2 = -cpsi1(cn+1) + cpsi1(one)
         cdin1 =  -2*(s1**2-1/(cn*(cn+1))*s1
     #     +1/(cn+1)**2+s2)
     #     +2-1/(cn*(cn+1))+2*s1
         coefq=coefq+xxalf/(2*xpi)*xcf*cdin1
      endif
      call acfun(cn,anq,ang)
      coefg=coefg+xxalf/(2*xpi)*ang
      coefq=coefq+xxalf/(2*xpi)*anq
      end


      subroutine acfun(cn,anq,ang)
c Coefficient functions for e+e- --> hadrons
c They are given by the MS_BAR result for the processes
c
c gamma* --> q+(g+q_bar) =
c [1+as/(2 pi)*antq]*3/8*{1+costh^2} + as/(2 pi)[anlq/2]*3/4*{sinth^2}
c 
c gamma* --> g+(q+q_bar) =
c [as/(2 pi)*antg]*3/8*{1+costh^2} + as/(2 pi)[anlg/2]*3/4*{sinth^2}
c
c     anq = antq+anlq/2
c     ang = antg+anlg/2
c
c      implicit complex * 16 (a-h,o-z)
      implicit none
c arguments
      complex * 16 cn,anq,ang
c parameters
      complex * 16 cone
      real * 8 xcf,pi,pi2
      parameter (pi=3.141592653589793d0,pi2=pi*pi,
     #     xcf=4.d0/3.d0,cone=(1.d0,0.d0)) 
c static local
      integer ini
c local
      complex * 16 psi01,psi11,psi21,psi0np1,psi1np1,s1n,s2n,
     # antq,antg,anlq,anlg
c functions
      complex * 16 cpsi0,cpsi1,cpsi2
      data ini/0/
      save ini,psi01,psi11,psi21
      if(ini.eq.0) then
         psi01=cpsi0(cone)
         psi11=cpsi1(cone)
         psi21=cpsi2(cone)
         ini=1
      endif
c- transverse
      psi0np1=cpsi0(cn+1)
      psi1np1=cpsi1(cn+1)
      s1n= psi0np1-psi01
      s2n=-psi1np1+psi11
c from Mele+Nason, NP B361(91)626, 3.5 for al, a.12 for antq+anlq/2
      antq = xcf*(
     #       -2*(2*cn+1)/(cn*(cn+1))**2-4*psi1np1+s1n**2
     #       -1/(cn*(cn+1))*s1n+1/(cn+1)**2+s2n
     #        +3.d0/2*s1n-3.d0/2/(cn+1)+2*pi2/3-9.d0/2)
c from Furmanski+Petronzio, Z. Phys. C11(1982)293
      antg = xcf*(
     #       (1/(cn-1)+2/((cn-1)*cn*(cn+1)) )*
     #       (psi01-psi0np1)-4/(cn-1)**2+4/cn**2-3/(cn+1)**2 )
c  longitudinal
c Mellin transform of a constant: int( x^(n-1),x,0,1)=1/n
      anlq = xcf*2 / cn
      anlg = xcf * 4/(cn*(cn-1))
c -2 to normalize to the total (quark) cross section
      anq = antq+anlq/2-2
      ang = antg+anlg/2
      end
