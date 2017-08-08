c      subroutine getmemberlha(mem)
c returns the member number after the last call to PDFSET (LHAGLUE)
c WARNING! works with 5.3.0; may break with different LHA versions!!!
c      integer mem
c      character*272 LHANAME
c      integer LHASET, LHAMEMB
c      common/LHAPDF/LHANAME, LHASET, LHAMEMB
c      mem=lhamemb
c      end

c Front end to mlmpdf0: it stores the arguments and return values of
c the nrec most recent calls to mlmpdf0. When invoked it looks in the
c stored calls; if a match its found, its return value is used.
c In this framework it is found that nrec=8 would be enough.
      subroutine mlmpdf8(ns,ih,xmu,x,fx,nf)
      implicit none
      integer ns,ih,nf
      real * 8 xmu,x,fx(-nf:nf)
      integer nrec
      parameter (nrec=10)
      real * 8 oxmu(nrec),ox(nrec),ofx(-6:6,nrec)
      integer ons(nrec),oih(nrec)
      integer irec
      save oxmu,ox,ofx,ons,oih,irec
c set to impossible values to begin with
      data ox/nrec*-1d0/
      data irec/0/
      integer j,k
      do j=1,nrec
         if(x.eq.ox(j)) then
            if(xmu.eq.oxmu(j)) then
               if(ns.eq.ons(j).and.ih.eq.oih(j)) then
                  do k=-nf,nf
                     fx(k)=ofx(k,j)
                  enddo
                  return
               endif
            endif
         endif
      enddo
      irec=irec+1
      if(irec.gt.nrec) irec=1
      ons(irec)=ns
      oih(irec)=ih
      oxmu(irec)=xmu
      ox(irec)=x
      call mlmpdf0(ns,ih,xmu,x,ofx(-6,irec),6)
      do k=-nf,nf
         fx(k)=ofx(k,irec)
      enddo
      end



      


      subroutine mlmpdf0(ndns,ih,q2,xs,fx,nf)
      implicit none
      integer ndns,ih,nf
      character * 20 parm(20)
      real * 8 val(20),tmp
      real * 8 q2,xs, fx(-nf:nf),afx(-6:6)
      integer nset,ns1,ns2
      save nset,ns1,ns2
      integer j
      character*10 lhaversion
      character*3 lhavs
      integer lhavi
      integer ini
      logical afgbug
      data ini/0/, afgbug/.false./
      save ini,afgbug
      if(ini.eq.0) then
         nset=0
c        check for AFG bug (only present in LHAPDF < 5.8.5)	 
         call getlhapdfversion(lhaversion)
c        extract the numbers  
         lhavs = lhaversion(1:1) // lhaversion(3:3) // lhaversion(5:5)
c        convert to integer
         read(lhavs,*) lhavi
c        check if bug is present	 
         if ( lhavi .lt. 585 ) then
            afgbug = .true.
            write(*,*) "LHAPDF version = ", lhaversion, " < 5.8.5, fixin
     #g the AFG bug"
         endif
         ini=1
      endif
      if(nset.eq.0) then
         nset=1
         parm(1)='DEFAULT'
         val(1)=ndns
         ns1=ndns
         call PDFSET(parm,val)
         if ( ih .eq. 4 .or. ih .eq. 5 ) then   ! photon or electron
           call evolvepdfpm(1,xs,sqrt(q2),0d0,0,afx)        
         else
          call evolvepdf(xs,sqrt(q2),afx)
         endif
         goto 998
      elseif(nset.eq.1) then
         if(ndns.eq.ns1) then
         if ( ih .eq. 4 .or. ih .eq. 5 ) then   ! photon or electron
           call evolvepdfpm(1,xs,sqrt(q2),0d0,0,afx)        
         else
           call evolvepdf(xs,sqrt(q2),afx)
         endif
             goto 998
         else
            nset=2
            parm(1)='DEFAULT'
            val(1)=ndns
            ns2=ndns
            call PDFSET(parm,val)
            if ( ih .eq. 4 .or. ih .eq. 5 ) then   ! photon or electron
              call evolvepdfpm(2,xs,sqrt(q2),0d0,0,afx)        
            else
              call evolvepdfm(2,xs,sqrt(q2),afx)
            endif
            write(*,*) ' mlmpdf0: PDF 1:'
            call GetDescM(1)
            write(*,*) ' mlmpdf0: PDF 2:'
            call GetDescM(2)
            write(*,*) ' ***************** '
            goto 998
         endif
      elseif(nset.eq.2) then
         if(ndns.eq.ns1) then
            if ( ih .eq. 4 .or. ih .eq. 5 ) then   ! photon or electron
              call evolvepdfpm(1,xs,sqrt(q2),0d0,0,afx)        
            else
              call evolvepdfm(1,xs,sqrt(q2),afx)
            endif
         elseif(ndns.eq.ns2) then
            if ( ih .eq. 4 .or. ih .eq. 5 ) then   ! photon or electron
              call evolvepdfpm(2,xs,sqrt(q2),0d0,0,afx)        
            else
              call evolvepdfm(2,xs,sqrt(q2),afx)
            endif
         else
            write(*,*) ' mlmpdf0: ndns is neither ns1 nor ns2: ',
     #            ndns,ns1,ns2
            stop
         endif
      endif
 998  continue
c our convention: 1/xs, since lhapdf returns parton momentum densities
      do j=-nf,nf
         fx(j)=afx(j)/xs
      enddo

c Switch to mlmpdf convention: 1=up, 2=down; LHAPDF: 1=down, 2=up.
c Same thing as in jetpdflib.f
      if (abs(ih).eq.1) then        ! proton
c switch to mlm convention: exch. up with down
          tmp=fx(1)
          fx(1)=fx(2)
          fx(2)=tmp
          tmp=fx(-1)
          fx(-1)=fx(-2)
          fx(-2)=tmp
      elseif (abs(ih).eq.2) then    ! neutron (actually not used)
          continue
      elseif (abs(ih).eq.3) then    ! pion
          tmp = fx(-2)
          fx(-2) = fx(2)
          fx(2)  = tmp
      elseif (abs(ih).eq.4) then     ! photon
         tmp     = fx(-1)
         fx(-1) = fx(-2)
         fx(-2) = tmp
         fx(1)=fx(-1)
         fx(2)=fx(-2)
      elseif (ih.eq.0) then
c this is (n+p)/2 in our convention
         fx(1)  = 0.5 * ( fx(1)+fx(2) )
         fx(-1) = 0.5 * ( fx(-1)+fx(-2) )
         fx(2)  = fx(1)
         fx(-2) = fx(-1)
      else
         write(*,*) ' ih was', ih,' instead of 0, +-1, +-2, +-3 or 4'
         stop
      endif

c correct for a bug in LHAPDF specific to AFG densities in the photon:
c the quarks must be divided by a factor of two.
c The same bug was present in PDFLIB (carried over with same code)
c NOTE: the bug has finally been fixed in LHAPDF v. 5.8.5 on 2/2/2011
c FONLL >= v1.3.3 comments out the lines below, and is expected to be 
c used with LHAPDF >= v5.8.5
      if((ih.eq.4).and.ndns.eq.363.and.afgbug.eqv..true.) then         
         do j=-nf,nf
            if(j.ne.0) fx(j)=fx(j)/2.d0
          enddo
      endif

c Finally, exchange particles with antiparticles if needed
      if (ih.lt.0) then
         do j=1,nf
            tmp=fx(j)
            fx(j)=fx(-j)
            fx(-j)=tmp
         enddo
      endif 
      end

      subroutine pdfpar(ns,ih,xlam,scheme,iret)
      implicit none
      integer ns,ih
      real * 8 xlam
      real * 8 fx(-6:6)
      character * 2 scheme
      integer iret
      character * 20 parm(2)
      real * 8 val(20)
      real * 8 QCDL4,QCDL5
      real*8 mz, asmz,newlam5,alphasPDF,genericxlambd
      parameter (mz=91.1876d0)
      COMMON/W50512/QCDL4,QCDL5
      integer iwritten,iwritten1,ihtmp
      save iwritten,iwritten1
      data iwritten /0/, iwritten1/0/
c dummy call to setup pdf if not yet done
      if(ih.eq.5) then    ! mlmpdf0 doesn't like ih=5
          ihtmp = 4
      else
          ihtmp=ih
      endif
      call mlmpdf0(ns,ihtmp,0.5d0,100d0,fx,6)
      parm(1)='DEFAULT'
      val(1)=ns
      call pdfset(parm,val)
c....LHAPDF v5.3.0 and 5.3.1 introduce the CTEQ65 and CTEQ65c sets,
c    but fail to set lambda_QCD values in the parameter files, returning 0
c    This fixes it and returns "reasonable values"
c      if((ns.ge.10350.and.ns.le.10390).or.
c     #   (ns.ge.10450.and.ns.le.10456)) then
c         QCDL5 = 0.226d0
c	 QCDL4 = 0.326d0
c        if (iwritten .le. 5 ) then
c        write(*,*) ' pdfpar: ********** WARNING !!!! ******************'
c	write(*,*) ' lamda^5 set to 0.226 for CTEQ65 set'
c        endif
c      endif	 
c      call getlam5(0,xlam) ! alternative way of setting xlam from LHAPDF

      asmz=alphasPDF(mz)  ! ges alphas from the lhapdf package
      if ( iwritten1 .lt. 10 ) then   ! avoid writing this out a thousand times when calling photons
         write(*,*) ' check: alpha_s(Mz)=',asmz
         iwritten1 = iwritten1 + 1
      endif
      newlam5=genericxlambd(asmz,mz,5) ! extract lambda5
      if (dabs((newlam5-QCDL5)/newlam5).gt.0.01) then
        write(*,*) ' pdfpar: ********** WARNING !!!! ******************'
        write(*,*) 
     #  'Calculated lam5 differs from common block one.'
        write(*,*) 'Calculated = ', newlam5
        write(*,*) 'Common block = ', QCDL5
        write(*,*) 'Using calculated one.'
      endif
      xlam=newlam5
      scheme='MS'
      if (iwritten .le. 5 ) then
       write(*,*) ' pdfpar: ********** WARNING !!!! ******************'
       write(*,*) ' PDF Scheme set to MSbar; if you are using a DIS pdf'
       write(*,*) ' the result will be wrong'
       if (iwritten .eq. 5 ) then
         write(*,*) ' *** Further warnings will be suppressed ***'
       endif
      iwritten = iwritten + 1
      endif
      iret=0 ! apparently gfortran doesn't return this (needed by setlam52 in hdms)
      end

      subroutine prntsf
      write(*,*) ' wrong pdf '
      stop
      end


C------- ALPHA QCD -------------------------------------
c Program to calculate alfa strong with nf flavours,
c as a function of lambda with 5 flavors.
c The value of alfa is matched at the thresholds q = mq.
c When invoked with nf < 0 it chooses nf as the number of
c flavors with mass less then q.
c
      function alfas(q2,xlam,inf)
      implicit real * 8 (a-h,o-z)
      data olam/0.d0/,pi/3.14159d0/
      data xmb/5.d0/,xmc/1.5d0/
      save
      if(xlam.ne.olam) then
        olam = xlam
        b5  = (33-2*5)/pi/12
        bp5 = (153 - 19*5) / pi / 2 / (33 - 2*5)
        b4  = (33-2*4)/pi/12
        bp4 = (153 - 19*4) / pi / 2 / (33 - 2*4)
        b3  = (33-2*3)/pi/12
        bp3 = (153 - 19*3) / pi / 2 / (33 - 2*3)
        xlc = 2 * log(xmc/xlam)
        xlb = 2 * log(xmb/xlam)
        xllc = log(xlc)
        xllb = log(xlb)
        c45  =  1/( 1/(b5 * xlb) - xllb*bp5/(b5 * xlb)**2 )
     #        - 1/( 1/(b4 * xlb) - xllb*bp4/(b4 * xlb)**2 )
        c35  =  1/( 1/(b4 * xlc) - xllc*bp4/(b4 * xlc)**2 )
     #        - 1/( 1/(b3 * xlc) - xllc*bp3/(b3 * xlc)**2 ) + c45
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
        alfas = 1/(b5 * xlq) -  bp5/(b5 * xlq)**2 * xllq
      elseif( nf .eq. 4 ) then
        alfas = 1/( 1/(1/(b4 * xlq) - bp4/(b4 * xlq)**2 * xllq) + c45 )
      elseif( nf .eq. 3 ) then
        alfas = 1/( 1/(1/(b3 * xlq) - bp3/(b3 * xlq)**2 * xllq) + c35 )
      else
        print *,'error in alfa: unimplemented # of light flavours',nf
        stop
      endif
      end


c.....extract lambda from alpha_s(q) value
      function genericxlambd(as,q,nf)
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)
      b  = (33-2*nf)/pi/12
      bp = (153 - 19*nf) / pi / 2 / (33 - 2*nf)
      t  = 1/b/as
    1 xlt = log(t)
      ot = t
c-----------------------------------------------------------
c Value and Derivative of alfa with respect to t
      as0  = 1/b/t - bp*xlt/(b*t)**2
      as1  = - 1/b/t**2 -bp/b**2*(1-2*xlt)/t**3
      t  = (as-as0)/as1 + t
      if(abs(ot-t)/ot.gt.0.00000001d0)goto 1
      genericxlambd = q/exp(t/2)
      return
      end

