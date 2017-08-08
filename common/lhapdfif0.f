c Front end to mlmpdf0 it stores the arguments and return values of
c the nrec most recent calls to mlmpdf0. When invoked it looks in the
c stored calls; if a match its found, its return value is used.
c In this framework it is found that nrec=8 would be enough.
c This provides a remarkable increase in spead (better than a factor of 3)
c when cteq6 pdf are used.
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
      save parm
      real * 8 val(20),tmp
      real * 8 q2,xs, fx(-nf:nf),afx(-6:6)
      integer ns
      save ns
      integer j
      integer ini
      data ini/0/
      save ini
      if(ini.eq.0) then
         parm(1)='DEFAULT'
         ini=1
         ns=0
      endif
      if(ns.ne.ndns) then
         val(1)=ndns
         ns=ndns
         call PDFSET(parm,val)
      endif
      if ( ih .eq. 4 .or. ih .eq. 5 ) then   ! photon or electron
   	   call evolvepdfp(xs,sqrt(q2),0d0,0,afx)        
      else
           call evolvepdf(xs,sqrt(q2),afx)
      endif
      
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
      if((ih.eq.4).and.ndns.eq.363)then
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
      integer iret,mem
      character * 20 parm(2)
      real * 8 val(20)
      real * 8 QCDL4,QCDL5
      COMMON/W50512/QCDL4,QCDL5
      integer iwritten
      save iwritten
      data iwritten /0/
c dummy call to setup pdf if not yet done
      call mlmpdf0(ns,ih,0.5d0,100d0,fx,6)
      parm(1)='DEFAULT'
      val(1)=ns
cc      call pdfset(parm,val)
c....LHAPDF v5.3.0 and 5.3.1 introduce the CTEQ65 and CTEQ65c sets,
c    but fail to set lambda_QCD values in the parameter files, returning 0
c    This fixes it and returns "reasonable values"
      if((ns.ge.10350.and.ns.le.10390).or.
     #   (ns.ge.10450.and.ns.le.10456)) then
         QCDL5 = 0.226d0
	 QCDL4 = 0.326d0
        if (iwritten .le. 5 ) then
        write(*,*) ' pdfpar: ********** WARNING !!!! ******************'
	write(*,*) ' lamda^5 set to 0.226 for CTEQ65 set'
        endif
      endif	 
c      call getlam5(0,xlam) ! alternative way of setting xlam from LHAPDF
      xlam=QCDL5
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
