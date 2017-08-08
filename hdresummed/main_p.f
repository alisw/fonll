CCCCCCCCCCCCCCCCCCCCCCsection efficace CCCCCCCCCCCCCCCCCCCCCCCCCC
C     *===================================================================
C     
C     
C     PHENOMENOLOGY PROGRAM TO CALCULATE HADRON - HADRON CROSS-SECTIONS
C     NUMERICALLY, STARTING FROM MATRIX ELEMENTS AT O(ALFAS**3).
C     FOR INCLUSIVE PARTICLE PRODUCTION.
C     ALL ENERGIES IN GEV
C==================================================================*/
      subroutine hdrs(result,error)
      implicit none
      real * 8 result,error
*     .. Parameters ..
      integer indim,imaxpts,ilenwrk
      parameter        (indim=3,imaxpts=3000*indim,ilenwrk=(indim+2)
     +     *(1+imaxpts/(2**indim+2*indim*indim+2*indim+1)))
      integer indim2,imaxpts2,ilenwrk2
      parameter        (indim2=2,imaxpts2=2000*indim2,
     +     ilenwrk2=(indim2+2)*(1+imaxpts2/(2**indim2
     +     +2*indim2*indim2+2*indim2+1)))
      integer ilw,iliw
      parameter        (ilw=1200,iliw=ilw/4)
C*******COMMUNICATION AVEC LES FONCTIONS DE STRUCTURE**************
      integer ipion
      common/pioncm/ipion
C*******************************************************************
      real * 8 pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      real * 8 m,mp,mu
      common/fonllscale/m,mp,mu
      integer iflag,ichoi
      real * 8 al,cq,v1,v2,v3,v4
      common/orde/iflag,ichoi,al,cq,v1,v2,v3,v4
      real * 8 zq
      common/fact/zq
C*******COMMUNICATION AVEC  INPUTP ********************************
      integer ipt,irap
      real * 8 ptevo,rap,enh1,enh2
      common/evopt/ipt,irap,ptevo(50),rap(50),enh1,enh2
      integer jmar
      common/valu/jmar
      integer isigm
      common/hdinput/isigm
      include 'xsect.h'
      real * 8 hc2,zrs,zal,cm,cmu,cmp
      common/choi/hc2,zrs,zal,cm,cmu,cmp
      integer isubtr
      common/isubtr0/isubtr
      real * 8 xmb,xmc
      common/hvqmass/xmb,xmc
      real * 8 alfa,beta
      integer npsm
      common/alfabeta/alfa,beta,npsm
      logical verbose
      common/chat/verbose  
      real * 8 sech1_delta
      common/sech/sech1_delta
      real * 8 xes1_delta
      common/xes/xes1_delta
      character*1 hvqs
      common/hvqtype/hvqs
      integer nfl
      common/flnumb/nfl
      integer imsvknfl,imsvscfl
      common/masskin/imsvknfl,imsvscfl
      integer ioutput
      common/out/ioutput
      real * 8 precint
      common/precintc/precint
      real * 8 lam5qcd
      COMMON/lamqcd/lam5qcd
      integer ih1,ih2
      common/hadr/ih1,ih2
      real * 8 h1pdf,h2pdf
      common/mypdfsets/h1pdf(20),h2pdf(20)
      integer iasyflag
      common/ciasyflag/iasyflag
      real * 8 qrkm,qrkm2,qm,pt,amt,rp,rpcm,teta,eta,raplim,
     #     csect,cerror,ebeam,pq,eq,bigt,bigu,bigs,bigv,bigw,nf,
     #     snf,coeff,sech1,err_sech1,
     #     err_sech1_delta,xes1,
     #     err_xes1,err_xes1_delta,sigto,facin,
     #     rs,alfas_p,accum,accerr,epsint,epsga,
     #     zmin,zmax
c function called
      character * 3 pkgname
c cernlib functions
c dgauss1 is a copy of dgauss, to avoid recursive calls
      real * 8 dasinh,dgauss1
c
      integer ik,ij
      real * 8 dplusp,dpluspv,ddel1p,ddel1pv,dplusp0,ddel1p0
      external dplusp,dpluspv,ddel1p,ddel1pv,dplusp0,ddel1p0
      pi=3.141592654d0
c************ No asymmetry: (q+qbar)/2 *********
      iasyflag=0
c*********this calls the input file ******************************
      call hdrsparam
c***************************************
      call fict
      ipion=1
c rs=root(S)
      rs=zrs
      cq=zq
      al=zal
      gs=rs**2
      n=3.
      vc=(n)**2-1.
      cf=4./3.d0
      v1=vc**2/n
      v2=vc/n
      v3=(n**4-1.)/2./n**2
      v4=vc**2/2./n**2
C******************************************************************
c******mass of the produced heavy quark
      if(hvqs.eq.'b') then
         qrkm = xmb
      else
         qrkm = xmc
      endif
      qrkm2 = qrkm**2
      if(imsvknfl.eq.1) then
         qm = qrkm
      else
         qm = 0d0
      endif
c****************************************
c        Loop on p_t
c...Dec 17 2007: commented out output to unit 32, hdrsout.tmp
c...Essentially useless in FONLL, and leaves the file behind
cc32      open(unit=32,file='hdrsout.tmp',status='UNKNOWN')
cc32      call toend(32)
cc32      write(32,'(a)')' en1      en2      cm  cmu cmp '//
cc32     # 'pt       rp       tot        err'
c reopen for each pt point, so we can read it while it runs'
cc32      close(32)
      do ik=1,ipt
         pt=ptevo(ik)
         amt = sqrt(qm**2 + pt**2)
c.... loop on rapidity
         do ij = 1,irap
cc32            open(unit=32,file='hdrsout.tmp',status='UNKNOWN')
cc32            call toend(32)
            rp = rap(ij)
c.......boost the rapidity to the cm frame
            if(enh1.eq.enh2) then
               rpcm = rp
            else  
               rpcm = rp - log(enh1/enh2)/2.
            endif

c....check if out of *massive* kinematical limit
c....Exit in case
            raplim = dasinh(sqrt(gs/4/(pt**2+qrkm**2) - 1d0))
            if(abs(rpcm).gt.raplim) then
               print*,' Out of kinematical limits'
               csect = 0
               cerror = 0
               goto 5556    
            endif

            if(rpcm.ne.0) then
cccmc               teta = atan(pt/amt/sinh(rpcm))
cc NB. teta evaluated equating pseudorapidity to rapidity (line below) 
cc is equivalent to the expression above when qm = 0
            teta = 2.*atan(exp(-rpcm))
            else
               teta = pi/2
            endif
            eta=2*pt/rs
            ebeam = rs/2
            pq = pt/sin(teta)
            eq = sqrt(pq**2 + qm**2)
            bigt = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
            bigu = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
            bigs = 4*ebeam**2
            bigv = 1d0 + bigt/bigs
            bigw = -bigu/(bigs + bigt)
c electron is coming from the right (see input.f)
c            if(pkgname().eq.'mlm'.and.ih2.eq.3)then
c               call elpdf_userz(nint(h2pdf(3)),zmin,zmax)
c               if(-bigt/(bigs + bigu).gt.zmax)then
c                  print*,' z range is empty'
c                  csect = 0
c                  cerror = 0
c                  goto 5556    
c               endif
c            endif
**********************
c.... these are equivalent to bigv,bigw when qm=0 and teta
c.....is calculated from rpcm by equating the latter to pseudorapidity
c            GV=1.-PT/RS/SIN(TETA)*(1.-COS(TETA))
c            GW=PT**2/GS/GV/(1.-GV)
            gv = bigv
            gw = bigw
c            print*, 'gv, gw = ', gv,gw
c            print*, 'bigv, bigw = ', bigv,bigw
**********************
            if(verbose) write(*,11) pt,rp,1-gv + gv*gw 
 11         format(1x,'pt = ',f6.2,' y = ',f5.2,'   x3min = ',d12.6)  
            if(verbose) then
               print*,'lam5qcd = ',lam5qcd
               print*,'alpha_s(p_t) = ',alfas_p(pt**2,lam5qcd,-5)
            endif
            pt2=pt**2
            if(imsvscfl.eq.1) then        
               m=dsqrt(pt2 + qrkm2)*cm
               mu=dsqrt(pt2 + qrkm2)*cmu
               mp=dsqrt(pt2 + qrkm2)*cmp
            else    
               m=dsqrt(pt2)*cm
               mu=dsqrt(pt2)*cmu
               mp=dsqrt(pt2)*cmp
            endif
            if(nfl.lt.0) then
               nf = snf(m**2)
            else
               nf = abs(nfl)
            endif
            gtr=nf/2.d0

C     IFLAG=1 (RESP 2) CORRESPOND A ALPHAS 1 BOUCLE (2 BOUCLES)
            iflag=iloopas
c
            if(isubtr.ne.1) then
               coeff=0
            else
c.....here you use the second moment of the fragmentation function
c.....to regularize the pole in x=1
               call subtract(mp**2,coeff)	
            endif
c
            accum=0
            accerr=0
            epsga=precint
            epsint = epsga*10.
C******************************************************************
            if(isubtr.eq.1) then
               sech1_delta =  dgauss1(ddel1p0,0.d0,1.d0,epsga)  * coeff
               err_sech1_delta=abs(sech1_delta*epsga)
               accum=accum+sech1_delta
               accerr=sqrt(accerr**2+err_sech1_delta**2)
            else
               sech1_delta = 0d0
               err_sech1_delta = 0d0
            endif
            if(verbose)
     #       write(*,*) ' sech1_delta=',sech1_delta,'+-',err_sech1_delta
C******************************************************************
            if(ias3term.ne.0.and.isubtr.eq.1) then
               call dovegas(2,dplusp0,imaxpts2,accum,accerr,epsga,
     #              coeff,xes1_delta,err_xes1_delta)
               accum=accum+xes1_delta
               accerr=sqrt(accerr**2+err_xes1_delta**2)
            else
               xes1_delta = 0d0
               err_xes1_delta = 0d0
            endif
            if(verbose)
     #       write(*,*) ' xes1_delta=',xes1_delta,'+-',err_xes1_delta
c**************     
            call dovegas(2,ddel1pv,imaxpts2,accum,accerr,epsga,
     #           1.d0,sech1,err_sech1)
            if(verbose)
     #           write(*,*) ' sech1=',sech1,'+-', err_sech1
            accum=accum+sech1
            accerr=sqrt(accerr**2+err_sech1**2)
c******************************************************************
            if(ias3term.ne.0) then
               call dovegas(3,dpluspv,imaxpts,accum,accerr,epsint,
     #           1.d0,xes1,err_xes1)
               accum = accum+xes1
               accerr = sqrt(accerr**2+err_xes1**2)
            else
               xes1 = 0
               err_xes1 = 0
            endif
            if(verbose)
     #           write (*,*) ' xes1',xes1,'+-',err_xes1
            if(verbose) then
               write(*,*) ' tot del =',sech1+sech1_delta,' +- ',
     #              sqrt(err_sech1**2+err_sech1_delta**2)
               write(*,*) ' tot plus=',xes1+xes1_delta,' +- ',
     #              sqrt(err_xes1**2+err_xes1_delta**2)
            endif
            if(verbose)
     #           write (*,*) ' unnormalized total',accum,'+-',accerr
            sigto=accum*hc2/pi/gs
            cerror=accerr*hc2/pi/gs
C     SIGTO EST E*DSIGMA/D3P POUR LE TERME  BORN+CORRECTIONS RADIATIVES
            if (isigm.eq.1) then
C     POUR OBTENIR DSIGMA/DY/DPT2
               facin=pi
            else if (isigm.eq.2) then
               facin=1.d0
            else if (isigm.eq.3) then
C     POUR OBTENIR DSIGMA/DPT/DY
               facin=pi*2.d0*pt
            endif
c     
            csect=facin*sigto
            cerror = cerror*facin
            if(verbose) then
               write (*,*) 'normalized_total ',csect,' +- ',cerror
               write (*,211) csect,cerror
            endif
 5556       if(ioutput.gt.0)
     #      write(ioutput,213) pt,rp,csect,cerror
cc32            write(32,132) enh1,enh2,cm,cmu,cmp,pt,rp,csect,cerror
cc32 132        format(1x,2(d8.2,1x),3(f3.1,1x),2(d8.2,1x),d10.4,1x,d7.1)
cc32            close(32)
C     /* END LOOP ON rapidity */
         enddo
C     /* END LOOP ON P_T */
      enddo
      result=csect
      error=cerror
 211  format(2(3x,d12.6))
 213  format(2(3x,f8.4),2(3x,d12.6))
 1    format(3x,'DS/DY/DPT2',//)
 2    format(3x,'E DS/D3P=',//)
 3    format(3x,'DS/DY/DPT=',//)
c.....check for pdflib errors
c     call pdfsta
      end

      subroutine subtract(q2,coeff)
      implicit none
c      implicit real*8 (a-h,o-z)
      integer nfl
      common/flnumb/nfl
      real * 8 q2,coeff
      real * 8 xxf(50)
      character*1 hvqs
      include 'xsect.h'
      complex * 16 ctmp
      integer ioutput
      common/out/ioutput
      logical verbose
      common/chat/verbose        
      real * 8 qp2
      common/myscale/qp2
      integer ipolemom
      common/ipole/ipolemom
      real * 8 xfragmom
      external xfragmom
      real * 8 pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/hvqtype/hvqs
      integer iasyflag
      common/ciasyflag/iasyflag
      real * 8 xmax,a,b,epsrel,sumr,xhint
c function called
      real * 8 dgauss
      qp2 = q2
      xmax = 1.-gv+gv*gw
      epsrel = 1.0d-5
      a = 0
      b = xmax
      xhint = dgauss(xfragmom,a,b,epsrel)
c
      ctmp = ipolemom
      call fragmell(ctmp,iloopfr,xxf)
      call parton(xxf,nfl)
      if(hvqs.eq.'b') then
         if(iasyflag.eq.0) then
            sumr = (xxf(2)+xxf(3))/2
         else
            sumr = xxf(2)
         endif
      elseif(hvqs.eq.'c') then
         if(iasyflag.eq.0) then
            sumr = (xxf(4)+xxf(5))/2
         else
            sumr = xxf(4)
         endif
      elseif(hvqs.eq.'g') then
         sumr = xxf(1)
      endif
c
      if(iasyflag.eq.0) then
         coeff = (sumr - xhint)*2
      else
         coeff = sumr - xhint
      endif
      if(verbose) write(*,12) xmax,xhint,sumr,coeff  
 12   format('xmax = ',f7.5,' xhint = ',f6.4,' mom2 = ',f6.4,
     &     ' coeff = ',f8.5)     
      end

      function xfragmom(x)
      implicit real*8(a-h,o-z)
      character*1 hvqs
      common/myscale/qp2
      common/hvqtype/hvqs
      common/ipole/ipolemom
      integer ifragmode
      integer iasyflag
      common/ciasyflag/iasyflag
      xc = x
c normal behaviour of fonfra
      ifragmode=0
      call fonfra(ifragmode,xc,ipi,qp2,xdup,xdubp,xddp,xddbp,xdsp,xdcp
     #     ,xdcbp,xdbp,xdbbp,xdtp,xdtbp,xdgp)
      if(hvqs.eq.'b') then
         if(iasyflag.eq.0) then
            xfragmom = (xdbp+xdbbp)/2
         else
            xfragmom = xdbp
         endif
      elseif(hvqs.eq.'c') then
         if(iasyflag.eq.0) then
            xfragmom = (xdcp+xdcbp)/2
         else
            xfragmom = xdcp
         endif
      elseif(hvqs.eq.'g') then
         xfragmom = xdgp
      else
         print*,'Non-implemented kind of heavy flavour'
         stop
      endif
      xfragmom = xfragmom*x**(ipolemom-2)
      end

      function ddel1pv(xx,wgt)
      implicit none
      real * 8 ddel1pv,wgt,ddel1p,xx(*),xxx(2)
      xxx(1)=xx(1)
      xxx(2)=1-xx(2)**2
      ddel1pv = 2*xx(2)*ddel1p(2,xxx)
      end

      function ddel1p(indim,xx)
      implicit none
      integer indim
      real * 8 ddel1p,ddel1,xx(indim)
       integer isubtr
      common/isubtr0/isubtr
      if(isubtr.eq.1) then
         ddel1p = ddel1(indim,xx,1)
         ddel1p = ddel1p - ddel1(indim,xx,2)
      else
         ddel1p = ddel1(indim,xx,1)
      endif
      end

      function ddel1p0(x)
      implicit none
      real * 8 ddel1p0,x,ddel1,xx(1)
      xx(1)=x
      ddel1p0 = ddel1(1,xx,3)
      end

      function ddel1(indim,xx,isubflag)
c      implicit real*8 (a-h,l-z)
      implicit none
      integer indim,isubflag
      real * 8 xx(indim),ddel1
      real * 8 pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      real * 8 m,mp,mu
      common/fonllscale/m,mp,mu
      integer iflag,ichoi
      real * 8 al,cq,v1,v2,v3,v4
      common/orde/iflag,ichoi,al,cq,v1,v2,v3,v4
      integer j0
      common/edf/j0
      integer ih1,ih2
      common/hadr/ih1,ih2
      real * 8 grrt,grrc
      common/fctd/grrt,grrc
      integer ichan,ichnum
      common/channels/ichan(16),ichnum
      integer ipolemom
      common/ipole/ipolemom
      include 'xsect.h'
      real * 8 mu2,m2,mp2,x3min,x3max,x3frag,x3jac,x3,vmin,vmax,
     #     v,a,un,ghd,bx1,bx2,cc,tmp,xjac,xlv,xlvmin,xlvmax

      character * 2 frscheme
      common/frschemec/frscheme
c functions invoked
      real * 8 f0,fdel,fvwpl,fvlo,alfai
c local
      real * 8 grrt1,grrt2,grrt3,grrt4
     #     ,grrt5,grrt6,grrt7,grrt8,grrt9,grrt10,grrt11,grrt12,grrt13
     #     ,grrt14,grrt15,grrt16
     #     ,grrc1,grrc2,grrc3,grrc4
     #     ,grrc5,grrc6,grrc7,grrc8,grrc9,grrc10,grrc11,grrc12,grrc13
     #     ,grrc14,grrc15,grrc16
      integer ifragmode,i
      xjac=1
      mu2=mu**2
      m2=m**2
      mp2=mp**2
      x3min=1.-gv+gv*gw
      x3max = 1d0
      if(isubflag.eq.1) then
c normal fragm. mode
         ifragmode=0
         x3frag = x3min+(x3max-x3min)*xx(2)
         x3 = x3frag
         x3jac = x3max-x3min
      elseif(isubflag.eq.2) then 
c only heavy quark->heavy quark from fragm. function, for pole subtraction
         ifragmode=1
         x3frag = x3min+(x3max-x3min)*xx(2)
         x3 = 1d0
         x3jac = x3max-x3min
c struct divides already by x3^3
         x3jac = x3jac*x3frag**3
c weight to get ipolemom moment
         x3jac = x3jac*x3frag**(ipolemom-2)
      elseif(isubflag.eq.3) then
c only heavy quark=1 in fragm. function to compensate for pole subtraction
         ifragmode=2
         x3frag = 1d0
         x3    = 1d0
         x3jac = x3frag**3
      else
         write(*,*) ' error: improper flag in ddel1'
         stop
      endif
      vmin=gv*gw/x3
      vmax=1.d0-(1.d0-gv)/x3
c some importance sampling: instead of
c      v=vmin+(vmax-vmin)*xx(1)
c      xjac = vmax-vmin
c use: log(v/(1-v)) 
      xlvmin = log(vmin/(1.d0-vmin))
      xlvmax = log(vmax/(1.d0-vmax))
      xlv = xlvmin + (xlvmax-xlvmin)*xx(1)
      xjac = xjac*(xlvmax-xlvmin)
      v = 1.d0-1.d0/(1.d0+exp(xlv))
      xjac = xjac*exp(xlv)/(1.d0+exp(xlv))**2
      a=gv*gw/v/x3
      un=1.d0
      ghd=0.d0
      if(a.ge.un) then
         goto 11
      endif
      bx1=gv*gw/v/x3
      bx2=(1.d0-gv)/(1.d0-v)/x3
      call stru(ifragmode,bx1,bx2,x3frag,ih1,ih2,m2,mp2
     #     ,grrt1,grrt2,grrt3,grrt4
     #     ,grrt5,grrt6,grrt7,grrt8,grrt9,grrt10,grrt11,grrt12,grrt13
     #     ,grrt14,grrt15,grrt16)
      call stru(ifragmode,bx2,bx1,x3frag,ih2,ih1,m2,mp2
     #     ,grrc1,grrc2,grrc3,grrc4
     #     ,grrc5,grrc6,grrc7,grrc8,grrc9,grrc10,grrc11,grrc12,grrc13
     #     ,grrc14,grrc15,grrc16)
      do i = 1,ichnum
         j0 = ichan(i)
         if (j0.eq.16.or.j0.eq.15) then
            cc=vc**2
         else if (j0.eq.14.or.j0.eq.13.or.j0.eq.10.or.j0.eq.9.or.j0
     #           .eq.8) then
            cc=vc*n
         else
            cc=n**2
         endif
         if (j0.eq.16) then
            grrt=grrt16
            grrc=grrc16
         else if (j0.eq.15) then
            grrt=grrt15
            grrc=grrc15
         else if (j0.eq.14) then
            grrt=grrt14
            grrc=grrc14
         else if (j0.eq.13) then
            grrt=grrt13
            grrc=grrc13
         else if (j0.eq.12) then
            grrt=grrt12
            grrc=grrc12
         else if (j0.eq.11) then
            grrt=grrt11
            grrc=grrc11
         else if (j0.eq.10) then
            grrt=grrt10
            grrc=grrc10
         else if (j0.eq.9) then
            grrt=grrt9
            grrc=grrc9
         else if (j0.eq.8) then
            grrt=grrt8
            grrc=grrc8
         else if (j0.eq.7) then
            grrt=grrt7
            grrc=grrc7
         else if (j0.eq.6) then
            grrt=grrt6
            grrc=grrc6
         else if (j0.eq.5) then
            grrt=grrt5
            grrc=grrc5
         else if (j0.eq.4) then
            grrt=grrt4
            grrc=grrc4
         else if (j0.eq.3) then
            grrt=grrt3
            grrc=grrc3
         else if (j0.eq.2) then
            grrt=grrt2
            grrc=grrc2
         else if (j0.eq.1) then
            grrt=grrt1
            grrc=grrc1
         endif
         if(grrt.ne.0.or.grrc.ne.0) then
            if (ias2term.ne.0) ghd=f0(v,x3)+ghd
            if (ias3term.eq.1) then
              ghd=(fdel(v,x3)+fvwpl(un,v,x3)*dlog(1.d0-a)+fvlo(un,v,x3)*
     #           (dlog(1.d0-a))**2/2.d0)/(8*cc*(1.d0-v))*alfai(mu2)
     #            + ghd
C  COMMENT *******************************************************
C        J0=16==>  PROCESSUS : G  G   ---> QJ
C        J0=15==>  PROCESSUS : G  G   ---> G
C        J0=14==>  PROCESSUS : QI G   ---> G
C        J0=13==>  PROCESSUS : QI G   ---> QI
C        J0=12==>  PROCESSUS : QI QBI ---> G
C        J0=11==>  PROCESSUS : QI QBI ---> QI
C        J0=10==>  PROCESSUS : QI G   ---> QBI
C        J0=9 ==>  PROCESSUS : QI G   ---> QBK
C        J0=8 ==>  PROCESSUS : QI G   ---> QK
C        J0=7 ==>  PROCESSUS : QI QI  ---> G
C        J0=6 ==>  PROCESSUS : QI QI  ---> QI
C        J0=5 ==>  PROCESSUS : QI QBI ---> QK
C        J0=4 ==>  PROCESSUS : QI QBK ---> G
C        J0=3 ==>  PROCESSUS : QI QBK ---> QI
C        J0=2 ==>  PROCESSUS : QI QK  ---> G
C        J0=1 ==>  PROCESSUS : QI QK  ---> QI
C*******************************************************************
            endif
c leave this feature for debugging: ias3term=2 does only the
c change of scheme contribution
            if (frscheme.eq.'DL'
     #           .and.(ias3term.eq.1.or.ias3term.eq.2).and.
     #           j0.ne.15.and.j0.ne.14.and.j0.ne.12.and.j0.ne.7.and.
     #           j0.ne.4.and.j0.ne.2) then
               tmp=1-v+gv*gw/x3
               ghd=ghd+v*f0(v,x3)*(-1.d0/v*alfai(mu2)*cf/(2*pi)*
     #              (log(1.d0-tmp)*(tmp**2+2*tmp-1.d0)-2*tmp
     #              +2*log(1.d0-tmp)**2))
            endif
         endif
      enddo
 11   continue
      ddel1=ghd*alfai(mu2)**2*x3jac*xjac
      return
      end

      function dpluspv(xx,wgth)
      implicit none
      real * 8  xx(*),xxx(3),dpluspv,wgth,dplusp
c 1 corresponds to v, 2 to w, 3 to x3
      xxx(1)=xx(1)
      xxx(2)=1-xx(2)**2
      xxx(3)=1-xx(3)**2
      dpluspv=4*xx(2)*xx(3)*dplusp(3,xxx)
      end

      function dplusp(indim,xx)
      implicit none
      integer indim
      real * 8 dplusp,xx(indim)
      integer isubtr
      common/isubtr0/isubtr
      real * 8 dplus,dplusp1,dplusp2
      if(isubtr.eq.1) then
         dplusp1 = dplus(indim,xx,1)
         dplusp2 = dplus(indim,xx,2)
	 dplusp = dplusp1 - dplusp2
      else
         dplusp = dplus(indim,xx,1)
      endif
      end

      function dplusp0(xx,wgth)
      implicit none
      real * 8 xx(*),xxx(2),wgth,dplusp0
      real * 8 dplus
      xxx(1)=xx(1)
      xxx(2)=1-xx(2)**2
      dplusp0 = 2*xx(2)*dplus(2,xxx,3)
      end

      function dplus(indim,xx,isubflag)
      implicit none
      integer indim,isubflag
      real * 8 xx(indim),dplus
      real * 8 pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      real * 8 m,mp,mu
      common/fonllscale/m,mp,mu
      integer iflag,ichoi
      real * 8 al,cq,v1,v2,v3,v4
      common/orde/iflag,ichoi,al,cq,v1,v2,v3,v4
      integer j0
      common/edf/j0
      integer ih1,ih2
      common/hadr/ih1,ih2
      real * 8 gppt,gppc
      common/fctc/gppt,gppc
      real * 8 grrt,grrc
      common/fctd/grrt,grrc
      integer ichan,ichnum
      common/channels/ichan(16),ichnum
      integer ipolemom
      common/ipole/ipolemom
      character * 2 frscheme
      common/frschemec/frscheme
      include 'xsect.h'
      real * 8 mu2,m2,mp2,x3min,x3max,x3frag,x3jac,x1,x2,x3,vmin,vmax,
     #     v,wmin,wmax,w,un,ghe,bx1,bx2,cc,cco,cpl,tmp,tmp1,tmp2,y,vv,
     # xjac,xlv,xlvmin,xlvmax
      integer i,ifragmode
c functions invoked
      real * 8 fresc,fvwpl,fvlo,alfai,f0
      real * 8 grrt1,grrt2,grrt3,grrt4
     #     ,grrt5,grrt6,grrt7,grrt8,grrt9,grrt10,grrt11,grrt12,grrt13
     #     ,grrt14,grrt15,grrt16
     #     ,grrc1,grrc2,grrc3,grrc4
     #     ,grrc5,grrc6,grrc7,grrc8,grrc9,grrc10,grrc11,grrc12,grrc13
     #     ,grrc14,grrc15,grrc16
     #     ,gppt1,gppt2,gppt3,gppt4
     #     ,gppt5,gppt6,gppt7,gppt8,gppt9,gppt10,gppt11,gppt12,gppt13
     #     ,gppt14,gppt15,gppt16
     #     ,gppc1,gppc2,gppc3,gppc4
     #     ,gppc5,gppc6,gppc7,gppc8,gppc9,gppc10,gppc11,gppc12,gppc13
     #     ,gppc14,gppc15,gppc16
      xjac = 1.d0
      mu2=mu**2
      m2=m**2
      mp2=mp**2
      x3min=1.d0-gv+gv*gw
      x3max= 1d0 
      if(isubflag.eq.1) then
c normal fragm. mode
         ifragmode=0
         x3frag = x3min+(x3max-x3min)*xx(3)
         x3 = x3frag
         x3jac = x3max-x3min
      elseif(isubflag.eq.2) then 
c only heavy quark->heavy quark from fragm. function, for pole subtraction
c the difference between isubflag=1 and isubflag=2 is integrated
         ifragmode=1
         x3frag = x3min+(x3max-x3min)*xx(3)
         x3 = 1.d0
         x3jac = x3max-x3min
c struct divides already by x3^3
         x3jac = x3jac*x3frag**3
c weight to get appropriate moment
         x3jac = x3jac*x3frag**(ipolemom-2)
      elseif(isubflag.eq.3) then
c used to add back the contribution of isubflag=2. It is afterwards multiplied by the
c appropriate moment of the quark component of the fragmentation function.
c Thus heavy quark=1 in fragm. function, all others are set to zero
         ifragmode=2
         x3frag = 1.d0
         x3    = 1.d0
         x3jac = x3frag**3
      else
         write(*,*) ' error: improper flag in dplusp'
         stop
      endif
      vmin=gv*gw/x3
      vmax=1.d0-(1.d0-gv)/x3
c some importance sampling: instead of
c      v=vmin+(vmax-vmin)*xx(1)
c      xjac = vmax-vmin
c use: log(v/(1-v)) as integration variable 
      xlvmin = log(vmin/(1.d0-vmin))
      xlvmax = log(vmax/(1.d0-vmax))
      xlv = xlvmin + (xlvmax-xlvmin)*xx(1)
      xjac = xjac*(xlvmax-xlvmin)
      v = 1.d0-1.d0/(1.d0+exp(xlv))
      xjac = xjac*exp(xlv)/(1.d0+exp(xlv))**2
c
      wmin=gv*gw/x3/v
      wmax=1.d0
      w=wmin+(wmax-wmin)*xx(2)
      xjac = xjac*(wmax-wmin)

      un=1.d0
      x1=gv*gw/v/w/x3
      x2=(1.d0-gv)/(1.d0-v)/x3
      call stru(ifragmode,x1,x2,x3frag,ih1,ih2,m2,mp2
     #     ,gppt1,gppt2,gppt3,gppt4
     #     ,gppt5,gppt6,gppt7,gppt8,gppt9,gppt10,gppt11,gppt12,gppt13
     #     ,gppt14,gppt15,gppt16)
      call stru(ifragmode,x2,x1,x3frag,ih2,ih1,m2,mp2
     #     ,gppc1,gppc2,gppc3,gppc4
     #     ,gppc5,gppc6,gppc7,gppc8,gppc9,gppc10,gppc11,gppc12,gppc13
     #     ,gppc14,gppc15,gppc16)
      bx1=gv*gw/v/x3
      bx2=(1.d0-gv)/(1.d0-v)/x3
      call stru(ifragmode,bx1,bx2,x3frag,ih1,ih2,m2,mp2
     #     ,grrt1,grrt2,grrt3,grrt4
     #     ,grrt5,grrt6,grrt7,grrt8,grrt9,grrt10,grrt11,grrt12,grrt13
     #     ,grrt14,grrt15,grrt16)
      call stru(ifragmode,bx2,bx1,x3frag,ih2,ih1,m2,mp2
     #     ,grrc1,grrc2,grrc3,grrc4
     #     ,grrc5,grrc6,grrc7,grrc8,grrc9,grrc10,grrc11,grrc12,grrc13
     #     ,grrc14,grrc15,grrc16)
      ghe=0.d0
      do i = 1,ichnum
         j0 = ichan(i)
         if (j0.eq.16.or.j0.eq.15) then
            cc=vc**2
         else if (j0.eq.14.or.j0.eq.13.or.j0.eq.10.or.j0.eq.9.or.j0
     #           .eq.8) then
            cc=vc*n
         else
            cc=n**2
         endif
         if (j0.eq.16) then
            gppt=gppt16
            gppc=gppc16
            grrt=grrt16
            grrc=grrc16
         else if (j0.eq.15) then
            gppt=gppt15
            gppc=gppc15
            grrt=grrt15
            grrc=grrc15
         else if (j0.eq.14) then
            gppt=gppt14
            gppc=gppc14
            grrt=grrt14
            grrc=grrc14
         else if (j0.eq.13) then
            gppt=gppt13
            gppc=gppc13
            grrt=grrt13
            grrc=grrc13
         else if (j0.eq.12) then
            gppt=gppt12
            gppc=gppc12
            grrt=grrt12
            grrc=grrc12
         else if (j0.eq.11) then
            gppt=gppt11
            gppc=gppc11
            grrt=grrt11
            grrc=grrc11
         else if (j0.eq.10) then
            gppt=gppt10
            gppc=gppc10
            grrt=grrt10
            grrc=grrc10
         else if (j0.eq.9) then
            gppt=gppt9
            gppc=gppc9
            grrt=grrt9
            grrc=grrc9
         else if (j0.eq.8) then
            gppt=gppt8
            gppc=gppc8
            grrt=grrt8
            grrc=grrc8
         else if (j0.eq.7) then
            gppt=gppt7
            gppc=gppc7
            grrt=grrt7
            grrc=grrc7
         else if (j0.eq.6) then
            gppt=gppt6
            gppc=gppc6
            grrt=grrt6
            grrc=grrc6
         else if (j0.eq.5) then
            gppt=gppt5
            gppc=gppc5
            grrt=grrt5
            grrc=grrc5
         else if (j0.eq.4) then
            gppt=gppt4
            gppc=gppc4
            grrt=grrt4
            grrc=grrc4
         else if (j0.eq.3) then
            gppt=gppt3
            gppc=gppc3
            grrt=grrt3
            grrc=grrc3
         else if (j0.eq.2) then
            gppt=gppt2
            gppc=gppc2
            grrt=grrt2
            grrc=grrc2
         else if (j0.eq.1) then
            gppt=gppt1
            gppc=gppc1
            grrt=grrt1
            grrc=grrc1
         endif
         if(gppt.ne.0.or.gppc.ne.0.or.grrt.ne.0.or.grrc.ne.0) then
            if(ias3term.eq.1) then
c... MC Dec 17, 2007: shift w to avoid 1/(1-w) division by zero
c....error below when w approaches 1 at very large rapidity (like 7-8 or so)
	       if (w .eq. 1d0) w = 1d0-1d-16
               cpl=((fvwpl(w,v,x3)/w-fvwpl(un,v,x3))/(1.d0-w)+(
     #           fvlo(w,v,x3)/w- fvlo(un,v,x3))*dlog(1.d0-w)/(1.d0-w))/
     #           (8*cc*(1.d0-v) )
               cco=fresc(w,v,x3)/w/(8*cc*(1.d0-v))
               ghe=cpl+cco + ghe
            endif
            if (frscheme.eq.'DL'
     #           .and.(ias3term.eq.1.or.ias3term.eq.2).and.
     #           j0.ne.15.and.j0.ne.14.and.j0.ne.12.and.j0.ne.7.and.
     #           j0.ne.4.and.j0.ne.2) then
               y=v*w+(1.d0-v)
               vv=v*w/y
c     in f0: common/fctd/gppt,gppc
c     here : common/fctd/grrt,grrc
c     First term must have structure functions at 3-body kinematics
               tmp1=grrt
               tmp2=grrc
               grrt=gppt
               grrc=gppc
c v/y is the jacobian, 1/y^2 is a remnant of 1/x3^3 *x3*D(x3)
c in the structure functions
               tmp=v/y**3*f0(vv,y*x3)
               grrt=tmp1
               grrc=tmp2
               tmp=tmp-v*f0(v,x3)
           tmp=tmp*(-2*log(1.d0-y)-1.d0)*(1.d0+y**2)/(1.d0-y)*cf/(2*pi)
               ghe=ghe+tmp
            endif
         endif
      enddo
      dplus=ghe*alfai(mu2)**3*x3jac*xjac
      return
      end

      real*8 function fdel(v,x3)
      implicit real*8 (a-h,l-z)
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/fonllscale/m,mp,mu
      common/orde/iflag,ichoi,al,cq,v1,v2,v3,v4
      common/fctd/grrt,grrc
      bx1=gv*gw/v/x3
      bx2=(1.d0-gv)/(1.d0-v)/x3
      shd=bx1*bx2*gs
      mu2=mu**2
      un=1.d0
      fkel=grrt*avdel(v,shd)
      fkelc=grrc*(avdel(1.d0-v,shd)+dlog(v/(1d0-v))*avwpl(un,1d0-v,shd)
     #     +.5d0*avlo(un,1.d0-v,shd)*(dlog((1.d0-v)/v))**2)*(1.d0-v)/v
      fdel=(fkel+fkelc)/shd
      return
      end

      real*8 function fvwpl(w,v,x3)
      implicit real*8 (a-h,l-z)
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/fonllscale/m,mp,mu
      common/orde/iflag,ichoi,al,cq,v1,v2,v3,v4
      common/fctc/gppt,gppc
      common/fctd/grrt,grrc
      x1=gv*gw/v/w/x3
      x2=(1.d0-gv)/(1.d0-v)/x3
      vx=1.d0-v*w
      wx=(1.d0-v)/(1.d0-v*w)
      sh=x1*x2*gs
      if (w.eq.1.d0) then
         fppt=grrt
         fppc=grrc
      else
         fppt=gppt
         fppc=gppc
      endif
      rvwpl=avwpl(w,v,sh)*fppt
      rvwplc=(avwpl(wx,vx,sh)+avlo(wx,vx,sh)*dlog(v/vx))*vx/v*fppc
      fvwpl=(rvwpl+rvwplc)/sh
      return
      end

      real*8 function fvlo(w,v,x3)
      implicit real*8 (a-h,l-z)
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/fonllscale/m,mp,mu
      common/orde/iflag,ichoi,al,cq,v1,v2,v3,v4
      common/fctc/gppt,gppc
      common/fctd/grrt,grrc
      x1=gv*gw/v/w/x3
      x2=(1.d0-gv)/(1.d0-v)/x3
      vx=1.d0-v*w
      wx=(1.d0-v)/(1.d0-v*w)
      sh=x1*x2*gs
      if (w.eq.1.d0) then
         fppt=grrt
         fppc=grrc
      else
         fppt=gppt
         fppc=gppc
      endif
      rvwlo=avlo(w,v,sh)*fppt
      rvwloc=avlo(wx,vx,sh)*fppc*vx/v
      fvlo=(rvwlo+rvwloc)/sh
      return
      end

      real*8 function fresc(ww,vv,x3)
      implicit real*8 (a-h,l-z)
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/fonllscale/m,mp,mu
      common/orde/iflag,ichoi,al,cq,v1,v2,v3,v4
      common/fctc/gppt,gppc
cmc
cmc Artifical cutoff, meant to tame problems arising at large
cmc rapidities, large energy and small pt (e.g. LHC). See ChangeLog.
cmc July 2007     
cmc
cmc Actually perhaps not needed after implementing struv_plus_avgo below?!?
cmc Septemner 2007
c      if ( ww .gt. 0.999 .and. vv .gt. 0.999) then
c         w = 0.999
c	 v = 0.999
c      else
         v = vv
	 w = ww
c      endif
cmc
      x1=gv*gw/v/w/x3
      x2=(1.d0-gv)/(1.d0-v)/x3
      vx=1.d0-v*w
      wx=(1.d0-v)/(1.d0-v*w)
      sh=x1*x2*gs
c... old implementation. It fails in some limits, like when v 
c... and/or w gets very close to 1
c      rresc=(struv(w,v,x3,sh)+avgo(w,v))*gppt
c      rrescc=(struv(wx,vx,x3,sh)+avgo(wx,vx))*gppc
c      write(9,*) "AA ", vv,ww,struv(w,v,x3,sh)+avgo(w,v)
c      write(9,*) "BB ", vx,wx,struv(wx,vx,x3,sh)+avgo(wx,vx)
c.... new implementation (september 2007). Seems to behave better 
c.... where the old implementation fails
      rresc=struv_plus_avgo(w,v,x3,sh)*gppt
      rrescc=struv_plus_avgo(wx,vx,x3,sh)*gppc

      fresc=(rresc+rrescc)/sh
      return
      end

C     THIS PROGRAM ALLOWS TO MOVE FROM A FACTORIZATION SCHEME TO ANOTHER
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     FRAGMENTATION   US(CQ=0) - MSBAR
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function fqqd(x)
      implicit real*8 (a-h,l-z)
      pi=4.*datan(1.d0)
      pi2=pi**2
      fqqd=4./3.*(-9./2.+2.*pi2/3.)
      return
      end

      real*8 function fqqw(x)
      implicit real*8 (a-h,l-z)
      fqqw=4./3.*(-3./2.+2.*(1.+x**2)*log(x)+(1.-x)**2*3./2.)
      return
      end

      real*8 function fqql(x)
      implicit real*8 (a-h,l-z)
      fqql=4./3.*(1.+x**2)
      return
      end

      real*8 function fqgl(x)
      implicit real*8 (a-h,l-z)
      fqgl=0.
      return
      end

      real*8 function fqgd(x)
      implicit real*8 (a-h,l-z)
      fqgd=0.
      return
      end

      real*8 function fqgw(x)
      implicit real*8 (a-h,l-z)
      fqgw=0.
      return
      end

      real*8 function fggl(x)
      implicit real*8 (a-h,l-z)
      fggl=0.
      return
      end

      real*8 function fggd(x)
      implicit real*8 (a-h,l-z)
      fggd=0.
      return
      end

      real*8 function fggw(x)
      implicit real*8 (a-h,l-z)
      fggw=0.
      return
      end

      real*8 function fgql(x)
      implicit real*8 (a-h,l-z)
      fgql=0.
      return
      end

      real*8 function fgqd(x)
      implicit real*8 (a-h,l-z)
      fgqd=0.
      return
      end

      real*8 function fgqw(x)
      implicit real*8 (a-h,l-z)
      fgqw=0.
      return
      end

C     THIS PROGRAM ALLOWS TO MOVE FROM A FACTORIZATION SCHEME TO ANOTHER
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     JMAR=0 OU 2    DIFFERENCE BETWEEN US(CQ=0) -  MSBARRE SCHEME
C     JMAR=3         DIFFERENCE BETWEEN US(CQ=0) -  MARTINELLI
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function cgqd(x)
      implicit real*8 (a-h,l-z)
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/valu/jmar
      pi2=pi**2
      if (jmar.eq.0.or.jmar.eq.2) then
         cgqd=0.d0
      else if (jmar.eq.1) then
         cgqd=0.d0
      else if (jmar.eq.3) then
         cgqd=-4./3.d0*(9./2.d0+pi2/3.d0)
      endif
      return
      end

      real*8 function cgqw(x)
      implicit real*8 (a-h,l-z)
      common/valu/jmar
      if (jmar.eq.0.or.jmar.eq.2) then
         cgqw=0.d0
      else if (jmar.eq.1) then
         cgqw=0.d0
      else if (jmar.eq.3) then
         cgqw=4./3.d0*(-3./2.d0-(1.+x**2)*dlog(x)+(1.-x)*(3.+2.*x))
      endif
      return
      end

      real*8 function cgql(x)
      implicit real*8 (a-h,l-z)
      common/valu/jmar
      if (jmar.eq.0.or.jmar.eq.2) then
         cgql=0.d0
      else if (jmar.eq.1) then
         cgql=0.d0
      else if (jmar.eq.3) then
         cgql=4./3.d0*(1.+x**2)
      endif
      return
      end

      real*8 function cqqd(x)
      implicit real*8 (a-h,l-z)
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/valu/jmar
      pi2=pi**2
      if (jmar.eq.0.or.jmar.eq.2) then
         cqqd=-(9./2.d0+pi2/3.d0)*cf
      else if (jmar.eq.1) then
         cqqd=0.d0
      else if (jmar.eq.3) then
         cqqd=0.d0
      endif
      return
      end

      real*8 function cqqw(x)
      implicit real*8 (a-h,l-z)
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/valu/jmar
      if (jmar.eq.0.or.jmar.eq.2) then
         cqqw=cf*(-3./2.d0+(3.+2.*x)*(1.-x)-(1.+x**2)*dlog(x))
      else if (jmar.eq.1) then
         cqqw=0.d0
      else if (jmar.eq.3) then
         cqqw=0.d0
      endif
      return
      end

      real*8 function cqql(x)
      implicit real*8 (a-h,l-z)
      common/valu/jmar
      if (jmar.eq.0.or.jmar.eq.2) then
         cqql=4./3.d0*(1.+x**2)
      else if (jmar.eq.1) then
         cqql=0.d0
      else if (jmar.eq.3) then
         cqql=0.d0
      endif
      return
      end

      real*8 function cqgd(x)
      implicit real*8 (a-h,l-z)
      common/valu/jmar
      if (jmar.eq.0.or.jmar.eq.2) then
         cqgd=0.d0
      else if (jmar.eq.1) then
         cqgd=0.d0
      else if (jmar.eq.3) then
         cqgd=0.d0
      endif
      return
      end

      real*8 function cqgw(x)
      implicit real*8 (a-h,l-z)
      common/valu/jmar
      if (jmar.eq.0.or.jmar.eq.2) then
         cqgw=0.d0
      else if (jmar.eq.1) then
         cqgw=0.d0
      else if (jmar.eq.3) then
         cqgw=-(1.-x)/2.d0*(-(x**2+(1.-x)**2)*dlog(x)+8.*x*(1.-x)-1.)
      endif
      return
      end

      real*8 function cqgl(x)
      implicit real*8 (a-h,l-z)
      common/valu/jmar
      if (jmar.eq.0.or.jmar.eq.2) then
         cqgl=0.d0
      else if (jmar.eq.1) then
         cqgl=0.d0
      else if (jmar.eq.3) then
         cqgl=-(1.-x)/2.d0*(x**2+(1.-x)**2)
      endif
      return
      end

      real*8 function cggd(x)
      implicit real*8 (a-h,l-z)
      common/valu/jmar
      if (jmar.eq.0.or.jmar.eq.2) then
         cggd=0.d0
      else if (jmar.eq.1) then
         cggd=0.d0
      else if (jmar.eq.3) then
         cggd=0.d0
      endif
      return
      end

      real*8 function cggw(x)
      implicit real*8 (a-h,l-z)
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/valu/jmar
      nf=gtr*2.d0
      if (jmar.eq.0.or.jmar.eq.2) then
         cggw=0.d0
      else if (jmar.eq.1) then
         cggw=0.d0
      else if (jmar.eq.3) then
         cggw=2.*nf*(1.-x)/2.*(-(x**2+(1.-x)**2)*dlog(x)+8.*x*(1.-x)-1.)
      endif
      return
      end

      real*8 function cggl(x)
      implicit real*8 (a-h,l-z)
      common/cons/pi,gs,gv,gw,n,gtr,cf,pt2,vc
      common/valu/jmar
      nf=gtr*2.d0
      if (jmar.eq.0.or.jmar.eq.2) then
         cggl=0.d0
      else if (jmar.eq.1) then
         cggl=0.d0
      else if (jmar.eq.3) then
         cggl=2.*nf*(1.-x)/2.d0*(x**2+(1.-x)**2)
      endif
      return
      end

      real*8 function snf(scale)
      implicit real*8 (a-h,l-z)
      common/alfacm/masch2,masbo2,masto2,lambda2
      if(scale.gt.masto2) then
         snf=6.
      else if(scale.gt.masbo2) then
         snf=5.
      else if(scale.gt.masch2) then
         snf=4.
      endif
      return
      end

      subroutine fict
      implicit real*8 (a-z)
      common/fact/zq
      zq=0.
      end

