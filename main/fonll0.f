      subroutine fonll0(icalctype,prefix,jprefix,
     #  ibeam1,ener1, nptype1, nngroup1, nnset1,
     #  ibeam2,ener2, nptype2, nngroup2, nnset2,
     #  xmq,xxlam,ffact,fren,pt,ylab,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
      implicit none
c prefix for generated files
      character * 30 prefix
c jprefix is last non-blank in prefix
      integer jprefix
      integer ibeam1,ibeam2,nptype1,nptype2
     #     ,nngroup1,nngroup2,nnset1,nnset2,icalctype
      real * 8 ener1,ener2,xmq,xxlam,ffact,fren,pt,ylab
      real * 8 hdmassive,hdmassless,hdresummed,hdreserr,sigmafonll,
     #         phmassive,phmassless,phresummed,phreserr
c ww common
      integer itypeww
      real * 8 xmuww,thcww,zminww,zmaxww,elphen
      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
c local variables
      character *1 hf
      character * 3 cbeam1,cbeam2
      character * 6 cbeams
      character * 2 sche1,sche2,scheph
      integer ib1rs,ib2rs,nlf,nf
      real * 8 xmprot,xlam,xlam1,xlam2,ebeamcm,yofcm,ycm,csmear
c
      integer iret,im0,irapph,ibeam,nptype,
     #        ngroup,nset,ngroup1,nset1,ngroup2,nset2
c logfile unit
      integer itm
      data itm/57/
      data xmprot/0.938d0/
      data csmear/5.d0/

c common for nuclear structure functions (see fonllmlmpdf in fonllpdf.f)
      integer atomicnumber
      common/shadowing/atomicnumber
      if((ibeam1.le.-2.and.ibeam2.le.-2).and.(ibeam1.ne.ibeam2)) then
            write(*,*) 'Different nuclei: not implemented'
	    stop
      endif
      atomicnumber = 0   
      if(ibeam1.le.-2) then 
         atomicnumber = - ibeam1
	 ibeam1 = 0
      endif	 
      if(ibeam2.le.-2) then 
         atomicnumber = - ibeam2
	 ibeam2 = 0
      endif 
c
c check that they are not both photons or electrons
c (no gamma-gamma process available here)
      if( (ibeam1.eq.4.or.ibeam1.eq.5) .and.
     #    (ibeam2.eq.4.or.ibeam2.eq.5)  ) then
         write(*,*) ' two incoming gamma/electron not available'
         stop
      endif
c convert ngroup, nset to our convention for pdflib numbers;
c this is necessary to link pdflibs through the jetpdflib interface
      nset1 = 100*nngroup1+nnset1
      ngroup1 = 0
      nset2 = 100*nngroup2+nnset2
      ngroup2 = 0
c 
c Find value of lambdaQCD if default has been requested
      xlam=xxlam
      call pdfpar(nset1,ibeam1,xlam1,sche1,iret)
      call pdfpar(nset2,ibeam2,xlam2,sche2,iret)
      if(xlam.le.0) then
         if(xlam1.eq.xlam2) then
            xlam=xlam1
         elseif(abs(ibeam1).le.2.and.abs(ibeam2).gt.2) then
            write(*,*)
     #        'Picking Lambda according to nucleon strf. (beam1)'
            xlam=xlam1
         elseif(abs(ibeam2).le.2.and.abs(ibeam1).gt.2) then
            write(*,*)
     #        'Picking Lambda according to nucleon strf. (beam2)'
            xlam=xlam2
         else
            write(*,*)
     # 'you have two incoming nucleon beams with different pdfs'
            write(*,*) ' Unable to set default value for Lambda'
            write(*,*) ' Choose a value.'
            stop
         endif
      endif
      if(   (abs(ibeam1).lt.4.and.sche1.ne.'MS')
     # .or. (abs(ibeam2).lt.4.and.sche2.ne.'MS')) then
         write(*,*)
     # 'Only MSbar scheme is implemented here!'
         write(*,*) ' Found ',sche1,' and ',sche2
         stop
      elseif(abs(ibeam1).ge.4.or.abs(ibeam2).ge.4) then
         if(abs(ibeam1).ge.4) then
            scheph=sche1
         else
            scheph=sche2
         endif
         if(scheph.eq.'DG') then
c it is called DG by jetpdflib
            scheph='DI'
         endif
         if(scheph.ne.'MS'.and.icalctype.ne.1) then
            write(*,*) 
     # 'This photon scheme is not implemented in the resummed code'
            write(*,*) ' Found: ',scheph
            stop
         endif
         if(scheph.ne.'DI'.and.scheph.ne.'MS') then
            write(*,*)
     # 'This photon scheme is not implemented in the resummed code'
            stop
         endif
      endif
c find # of flavours
      if(xmq.lt.3) then
         write(*,*) ' assuming charm quark'
         nlf=3
         nf=4
         hf='c'
      elseif(xmq.lt.7) then
         write(*,*) ' assuming bottom quark'
         nlf=4
         nf=5
         hf='b'
      else
         write(*,*) ' top not implemented'
         stop
      endif
c no bottom for resummed photon!
      if( icalctype.ne.2.and.(ibeam1.eq.4.or.ibeam1.eq.5 .or.
     #     ibeam2.eq.4.or.ibeam2.eq.5 ) .and. hf.eq.'b') then
         write(*,*) ' Resummed Photoproduction program'
         write(*,*) ' does not work for bottom'
         stop
      endif
c
c Find CM beam energies and quark rapidity
      if(ener1.gt.0.and.ener2.gt.0) then
         ebeamcm=0.5d0*sqrt((ener1+ener2)**2-(ener1-ener2)**2)
         write(*,*) 'CM energy=',2*ebeamcm
         yofcm=0.5d0*log( ((ener1+ener2)+(ener1-ener2))/
     #                    ((ener1+ener2)-(ener1-ener2))  )
      elseif(ener2.eq.0) then
         if(ibeam2.ne.1.and.ibeam2.ne.2.and.ibeam2.ne.0) then
            write(*,*)' beam 2 at fixed target not possible'
            stop
         endif
c neglect proton-neutron mass differences
         ebeamcm=0.5d0*sqrt(2*xmprot*ener1)
         yofcm=0.5d0*log((2*ener1+xmprot)/xmprot)
      elseif(ener1.eq.0) then
         write(*,*)' beam 1 must have positive rapidity'
         write(*,*)' cannot be the target'
         stop
      endif
c Find cm rapidity of quark
      ycm = ylab - yofcm

c Check the phase space limits. If out of them, immediately return zeros
c Calculate ymax from E < sqrt(s)/2
      call checkphsp(ibeam1,ibeam2,ebeamcm,xmq,pt,ycm,iret)
      if(iret.eq.1) then
         hdmassive=0
         hdmassless=0
         hdresummed=0
         hdreserr=0
         phmassive=0
         phmassless=0
         phresummed=0
         phreserr=0
	 sigmafonll=0
	 return
      endif

c Symbolic names for beam type
      call sybname(ibeam1,cbeam1)
      call sybname(ibeam2,cbeam2)
c ************************************************************
c *********   HADRONIC MASSIVE *******************************
c ************************************************************

c prepare input file
      do im0=0,1
c im0=0: massive, 1: massless
      if(im0.eq.0) then
         open(unit=itm, file=prefix(1:jprefix)//'hdmv.tmp')
      else
         open(unit=itm, file=prefix(1:jprefix)//'hdml.tmp')
      endif
      write(itm,'(a)') '80          ! # of gaussian points'
      write(itm,*) nptype1,ngroup1,nset1,
     #                      ' ! nptype, ngroup, nset (beam1)'
      write(itm,*) nptype2,ngroup2,nset2,
     #                      ' ! nptype, ngroup, nset (beam2)'
      write(itm,*) xlam, ' ! Lambda_5'
      write(itm,'(a)') 'dfpy       ! task'
      write(itm,*) '0          ! 1 to convert dsig/dpt^2 --> dsig/dpt'    
      write(itm,*) nlf,'       !# of light flavours'
      write(itm,*) '1    ! 1 for the nlf+1  scheme, 0 for nlf scheme'
      write(itm,*) '1    ! alpha_S --> alpha_S/rscalpha'
      write(itm,'(a,d14.8,a)') cbeam1//' '//cbeam2//'   ',2*ebeamcm,
     # '     ! process and CM energy'
      write(itm,*) '1      ! 1 for CM, 2 for FT, 3 for lab, y and eta'
      write(itm,'(a)') 'al          ! gg, qq, qg, any other for all'
      if(im0.eq.0) then
         write(itm,*) '1      '//
     #   '! 1 mu0=sqrt(m^2+pt^2), 2 mu0=pt, 3 sequence (in mu0)'
      elseif(im0.eq.1) then
         write(itm,*) '2      '//
     #   '! 1 mu0=sqrt(m^2+pt^2), 2 mu0=pt, 3 sequence (in mu0)'
      endif
      write(itm,*) ffact,'       ! xmu_fact(hadron)=ffact*mu0'
      write(itm,*) '0            ! 0 normal, 1 for tests'
      write(itm,*) fren,'     ! xmu_ren=fren*mu0'
      write(itm,*) '-1  ! mass value, <0 to set rapidity and pt values'
      write(itm,*) ycm
      write(itm,*) '*'
      write(itm,*) pt
      write(itm,*) '*'
      if(im0.eq.0) then
         write(itm,*) xmq, ' 0  ! 0 for massive, 1 for massless'
      else
         write(itm,*) xmq, ' 1  ! 0 for massive, 1 for massless'
         write(itm,*) '1    ! convert pt to mt'
      endif
      write(itm,'(a)') '.'
      close(itm)
      enddo
c ************************************************************
c *********   HADRONIC RESUMMED *******************************
c ************************************************************
      open(unit=itm,file=prefix(1:jprefix)//'hdrs.tmp')
c type convention different in this program
c      if(ibeam1.eq.1) then
c         ib1rs=0
c      elseif(ibeam1.eq.-1) then
c         ib1rs=1
c      elseif(ibeam1.eq.4) then
c         ib1rs=2
c      elseif(ibeam1.eq.5) then
c         ib1rs=3
c      else
c         write(*,*) 'beam 1 type not implemented in hdrs'
c         stop
c      endif
c      if(ibeam2.eq.1) then
c         ib2rs=0
c      elseif(ibeam2.eq.-1) then
c         ib2rs=1
c      elseif(ibeam2.eq.4) then
c         ib2rs=2
c      elseif(ibeam2.eq.5) then
c         ib2rs=3
c      else
c         write(*,*) 'beam 2 type not implemented in hdrs'
c         stop
c      endif
      ib1rs=ibeam1
      ib2rs=ibeam2
      write(itm,*) '''hdres.tmp'''
      write(itm,*) ' 4                   ! pole moment'
      write(itm,*) '1                    ! isigm'
      write(itm,*) '2                    ! iloopfr'
      write(itm,*) '0                    ! iwhichsfh'
      write(itm,*) '2                    ! iloop in alfas'
      write(itm,*) '1                    ! ias2term'
      write(itm,*) '1                    ! ias3term'
      write(itm,*) '2                    ! jmar'
      write(itm,*) ' ''MS''                ! fragmentation scheme'
      write(itm,*) '0                    ! alternative evmat'
      write(itm,*) pt
      write(itm,*) -1
      write(itm,*) '1             ! to convert pt->sqrt(pt^2+m^2)'
      write(itm,*) ycm
      write(itm,*) 2000
      write(itm,*) ib1rs, ebeamcm, nptype1, ngroup1, nset1
c     #, '    had. 1, ener., pdf set'
      write(itm,*) ib2rs, ebeamcm, nptype2, ngroup2, nset2
c     #, '    had. 2, ener., pdf set'
      write(itm,*) xlam,'      ! Lambda_5 (GeV)'
      write(itm,*) nf,'         ! number of flavours'
      write(itm,*) ffact,fren,1
c     #, '   ! fact. and ren. scale factors'
      write(itm,*)'0        ! 0 for mu0*cm, 1 for mu0 in str. funct.'
      if(hf.eq.'b') then
         write(itm,*)1.5,xmq,'      ! charm mass, bottom mass'
      elseif(hf.eq.'c') then
         write(itm,*) xmq, 5, '      ! charm mass, bottom mass'
      else
         write(*,*) ' flavour not implemented in hdrs:',hf
      endif
      write(itm,'(a)')hf//'        ! flavour with singular behaviour'
      write(itm,*)'1       ! isubtr for pole subtraction 1=yes, 0=no'
      write(itm,*)' 1.0    ! rescaling factor for fragm. initial scale'
      write(itm,*)'12     ! isuda: Sudakov form factors (1) or not (0)'
      write(itm,*)'0     0.000     0.000 ! npsm, alfa, beta'
      write(itm,*)'16                   ! # of channel to include'
      write(itm,*)'  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16'
      write(itm,*)' 0.1000E-02          ! precision of the integrals'
      write(itm,*)'   1.000             ! alpha_s scale factor'
      close(itm)
c ************************************************************
c The following for one incoming photon or electron
c ************************************************************
      if(    ibeam1.eq.4.or.ibeam1.eq.5
     #  .or. ibeam2.eq.4.or.ibeam2.eq.5 ) then
         if(ibeam1.eq.4.or.ibeam1.eq.5) then
c positive photon rapidity
            irapph=1
c photon or electron energy in the lab
            elphen=ener1
c hadron beam:
            ibeam=ibeam2
            nptype=nptype2
            ngroup=ngroup2
            nset=nset2
            cbeams=cbeam1//cbeam2
         else
c negative photon rapidity
            irapph=-1
c photon or electron energy in the lab
            elphen=ener2
c hadron beam:
            ibeam=ibeam1
            nptype=nptype1
            ngroup=ngroup1
            nset=nset1
            cbeams=cbeam2//cbeam1
         endif
c ************************************************************
c *********   Pointlike massive and massless *****************
c ************************************************************
         do im0=0,1
         if(im0.eq.0) then
               open(unit=itm,file=prefix(1:jprefix)//'phmv.tmp')
         else
               open(unit=itm,file=prefix(1:jprefix)//'phml.tmp')
         endif
         write(itm,*) '80       ! # of gaussian points'
         write(itm,*) nset,'   ! nptype, ngroup, nset'
         write(itm,*) xlam,'    ! Lambda_5'
         write(itm,*) '''',scheph,'''   ! scheme for photon'
         write(itm,'(a)') 'dfpy         ! task'
         write(itm,*)'0       ! 1 to convert dsig/dpt^2 --> dsig/dpt'
         write(itm,*) nlf,'        ! # of light flavours'
         write(itm,*)'1    ! 1 for the nlf+1  scheme, 0 for nlf scheme'
         write(itm,*)'0.1000E+01   ! alpha_S --> alpha_S/rscalpha'
         write(itm,'(a,4x,d14.8,a)') cbeams, 2*ebeamcm,
     #     '   ! beams and cm energy'
         write(itm,*) '1   ! 1 for CM, 2 for FT, 3 for lab, y and eta'
c Weizsaeker-Williams parameters in case of electron beam
         if(cbeams(1:2).eq.'el') then
            write(itm,*) zminww,zmaxww,'    ! zminww,zmaxww'
         endif
         write(itm,'(a)')'+al   ! pg, pq, any other for all'
         write(itm,*)'1      ! 1 mu0=sqrt(m^2+pt^2)'
         write(itm,*) ffact, '  ! xmu_fact(hadron)=ffact*mu0'
         write(itm,*)'0            ! 0 normal, 1 for tests'
         write(itm,*) ffact, '   ! xmu_fact(photon)=ffactph*mu0'
         write(itm,*) fren,  '   ! xmu_ren=fren*mu0'
         write(itm,*) '-1   !mass value, to set rapidity and pt values'
         write(itm,*) ycm*irapph
         write(itm,'(a)')'*'
         write(itm,*) pt
         write(itm,'(a)')'*'
         if(im0.eq.0) then
            write(itm,*) xmq, ' 0  ! 0 for massive, 1 for massless'
         else
            write(itm,*) xmq, ' 1  ! 0 for massive, 1 for massless'
            write(itm,*) '1    ! convert pt to mt'
         endif
         write(itm,*) '-100.             ! 3*charge, -100 for default'
         write(itm,'(a)')'.'
         close(itm)
         enddo
c ************************************************************
c *********   Pointlike resummed *****************************
c ************************************************************
         open(unit=itm,file=prefix(1:jprefix)//'phrs.tmp')
         write(itm,*)'''phres.tmp'''
         write(itm,*)'0 0 0   ! 0 for no smearing'
c in fact, only charm is implemented; just in case
c we implement bottom in the future, leave this in
         if(hf.eq.'b') then
            write(itm,*)1.5,xmq,'      ! charm mass, bottom mass'
         elseif(hf.eq.'c') then
            write(itm,*) xmq, 5, '      ! charm mass, bottom mass'
         else
            write(*,*) ' flavour not implemented in phrs:',hf
         endif
         write(itm,*)'1    ! rescaling for fragm. initial scale'
         write(itm,*) hf
         write(itm,*)'12    ! no Sudakov (1=yes)'
         write(itm,*)'2    ! loops in alphas'
         write(itm,*)'2    ! order of computation'
         write(itm,*) fren, ffact
         write(itm,*) 
     #        ibeam, ebeamcm, nptype,ngroup,nset,ebeamcm
c     #  ,' !ih.,E_H,ntype,ngroup,nset,E_el/gam'
         write(itm,*) xlam,'         ! lambdaQCD, <=0 for default'
         write(itm,*)'1       ! # of pt points'
         write(itm,*) pt
         write(itm,*) '1        ! 1 for pt->sqrt(pt^2+m^2'
         write(itm,*) '1        ! number of rapidity points'
         write(itm,*) -ycm*irapph
         if(cbeams(1:2).eq.'el') then
c an electron
            write(itm,*) '1           ! use Weizsaeker-WIlliams'
            write(itm,*) zminww,zmaxww,'    ! zminww,zmaxww'
         else
            write(itm,*) '0         ! a photon'
         endif
         write(itm,*) '''MS''      ! FF factorization scheme'
         write(itm,*) '0             ! real or fake hvq pdf'
         write(itm,*) '1           ! alpha_s scale factor'
         write(itm,*) '0           ! resummed pdf and ff'
         write(itm,*) '3000        ! number of vegas calls'
         write(itm,*) '10          ! # of Vegas iterations'
         close(itm)
      endif

c Run the program
      if(icalctype.eq.1.or.icalctype.eq.2.or.icalctype.eq.4) then
         open(unit=55,file=prefix(1:jprefix)//'hdmv.tmp')
         call hdms(hdmassive)
         hdmassive=hdmassive*1.d6
         close(55)
      else
         hdmassive=0
      endif
      if(icalctype.eq.1.or.icalctype.eq.4) then
         open(unit=55,file=prefix(1:jprefix)//'hdml.tmp')
         call hdms(hdmassless)
         hdmassless=hdmassless*1.d6
         close(55)
      else
         hdmassless=0   
      endif
      if(icalctype.eq.1.or.icalctype.eq.3) then
         open(unit=55,file=prefix(1:jprefix)//'hdrs.tmp')
         call hdrs(hdresummed,hdreserr)
         close(55)
      else
         hdresummed=0
         hdreserr=0
      endif
      if(    ibeam1.eq.4.or.ibeam1.eq.5
     #  .or. ibeam2.eq.4.or.ibeam2.eq.5 ) then
         if(icalctype.eq.1.or.icalctype.eq.2) then
            open(unit=55,file=prefix(1:jprefix)//'phmv.tmp')
            call phms(phmassive)
            phmassive=phmassive*1.d6
            close(55)
         else
            phmassive=0
         endif
         if(icalctype.eq.1) then
            open(unit=55,file=prefix(1:jprefix)//'phml.tmp')
            call phms(phmassless)
            phmassless=phmassless*1.d6
            close(55)
         else
            phmassless=0
         endif
         if(icalctype.eq.1.or.icalctype.eq.3) then
            open(unit=55,file=prefix(1:jprefix)//'phrs.tmp')
            call phrs(phresummed,phreserr)
            close(55)
         else
            phresummed=0
            phreserr=0
         endif
      else
         phmassive=0
         phmassless=0
         phresummed=0
         phreserr=0
      endif
c Cross section in pb/GeV^2
      sigmafonll=hdmassive+phmassive+
     #  (hdresummed-hdmassless+phresummed-phmassless)*pt**2/
     #             (pt**2+(csmear*xmq)**2)
      end


      subroutine sybname(ibeam,cbeam)
      implicit none
      character * 3 cbeam
      integer ibeam
      if(ibeam.eq.1) then
         cbeam='pr+'
      elseif(ibeam.eq.-1) then
         cbeam='pr-'
      elseif(ibeam.eq.0) then
         cbeam='nu+'
      elseif(ibeam.eq.2) then
         cbeam='n +'
      elseif(ibeam.eq.-2) then
         cbeam='n -'
      elseif(ibeam.eq.3) then
         cbeam='pi+'
      elseif(ibeam.eq.-3) then
         cbeam='pi-'
      elseif(ibeam.eq.4) then
         cbeam='ph '
      elseif(ibeam.eq.5) then
         cbeam='el '
      else
         write(*,*) ' invalid beam',ibeam
         stop
      endif
      end

      subroutine sybnumb(cbeam,ibeam)
      implicit none
      character * 3 cbeam
      integer ibeam
      if(cbeam.eq.'pr+') then
         ibeam=1
      elseif(cbeam.eq.'pr-') then
         ibeam=-1
      elseif(cbeam.eq.'nu+') then
         ibeam=0
      elseif(cbeam.eq.'n +') then
         ibeam=2
      elseif(cbeam.eq.'n -') then
         ibeam=-2
      elseif(cbeam.eq.'pi+') then
         ibeam=3
      elseif(cbeam.eq.'pi-') then
         ibeam=-3
      elseif(cbeam.eq.'ph') then
         ibeam=4
      elseif(cbeam.eq.'el') then
         ibeam=5
      else
         write(*,*) ' invalid beam',cbeam
         stop
      endif
      end

c Evaluate phase space limits.
c Also taking care of boost to photon rest frame in e-p case
c Calculate ymax from E < sqrt(s)/2
      subroutine checkphsp(ibeam1,ibeam2,ebeamcm,xmq,pt,ycm,iret)
      implicit none
      integer ibeam1,ibeam2,iret
      real * 8 ebeamcm,xmq,pt,ycm
c ww common
      integer itypeww
      real * 8 xmuww,thcww,zminww,zmaxww,elphen
      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
c local
      real * 8 e1,e2,yofcm,ycmp,ymax,ebeamcmp
c functions
      real * 8 dasinh
      e1=ebeamcm
      e2=e1
      if(ibeam1.eq.5) e1=e1*zmaxww
      if(ibeam2.eq.5) e2=e2*zmaxww
      ebeamcmp=0.5d0*sqrt((e1+e2)**2-(e1-e2)**2)
      yofcm=0.5d0*log( ((e1+e2)+(e1-e2))/
     #                    ((e1+e2)-(e1-e2))  )
c Find cm rapidity of quark
      ycmp = ycm - yofcm
      ymax = dasinh(sqrt(ebeamcmp**2/(pt**2+xmq**2) - 1d0))
      if(abs(ycmp).gt.ymax) then
         iret=1
      else
         iret=0
      endif
      end

