************************************************************************
*                                                                      *
*                    FONLL            version 1.1 (Nov 2003)           *
*                                                                      *
*      Program to calculate heavy quark transverse momentum and        *
*      rapidity distributions in hadron-hadron and photon-hadron       *
*      collisions, matching Fixed Order next-to-leading order terms    *
*      and Next-to-Leading-Log large-p_T resummation.                  *
*                                                                      *
*      by Matteo Cacciari, Stefano Frixione, Paolo Nason               *
*                                                                      *
*      Built upon the following papers/codes:                          *
*                                                                      *
*      S. Frixione and P. Nason,                                       *
*           "Phenomenological study of charm photoproduction at HERA", *
*           JHEP 0203 (2002) 053 [hep-ph/0201281]                      *
*      M. Cacciari, S. Frixione and P. Nason,                          *
*           "The p_T Spectrum in Heavy Flavor Photoproduction",        *
*           JHEP 0103 (2001) 006 [hep-ph/0102134]                      *
*      M. Cacciari, M. Greco and P. Nason,                             *
*           "The p_T Spectrum in Heavy Flavor Hadroproduction",        *
*           JHEP 9805 (1998) 007 [hep-ph/9803400]                      *
*      M. Cacciari and M. Greco,                                       *
*           "Charm Photoproduction via Fragmentation",                 *
*           Z.Phys. C69 (1996) 459 [hep-ph/9505419]                    *
*      M. Cacciari and M. Greco,                                       *
*           "Large p_T Hadroproduction of Heavy Quarks",               *
*           Nucl.Phys. B421 (1994) 530 [hep-ph/9311260]                *
*      B. Mele and P. Nason,                                           *
*           "The Fragmentation Function for Heavy Quarks in QCD",      *
*           Nucl.Phys. B361 (1991) 626                                 *
*      P. Nason, S. Dawson and R.K. Ellis,                             *
*           "The One-Particle Inclusive Differential Cross-Section     *
*            for Heavy Quark Production in Hadronic Collisions",       *
*           Nucl.Phys. B327 (1989) 49                                  *
*      F. Aversa, P. Chiappetta, M. Greco and J.-Ph. Guillet,          *
*           "QCD Corrections to Parton-Parton Scattering Processes",   *
*           Nucl.Phys. B327 (1989) 105                                 *
*      R.K. Ellis and P. Nason,                                        *
*           "QCD Radiative Corrections to the Photoproduction          *
*            of Heavy Quarks",                                         *
*           Nucl.Phys. B312 (1989) 551                                 *
*      P. Aurenche, R. Baier, A. Douiri, M. Fontannaz and D. Schiff,   *
*           "Scheme Invariant Higher Order QCD Predictions for         *
*            Large p_T Photoproduction Reactions",                     *
*           Nucl.Phys. B286 (1987) 553                                 *
*                                                                      *
************************************************************************
c Driving program for FONLL calculations
      program fonll
      implicit none
c parameters
      integer nymax,nptmax
      parameter (nymax=100,nptmax=100)
c common
c ww
      integer itypeww
      real * 8 xmuww,thcww,zminww,zmaxww,elphen
      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
      integer jptseqfl
c which pt sequence
      common/cjptseqfl/jptseqfl
c local
c prefix for generated files
      character * 30 prefix
c jprefix is last non-blank in prefix
      integer jprefix
      integer ibeam1,ibeam2,nptype1,nptype2,ngroup1,ngroup2,nset1,nset2
     # ,icalctype
      real * 8 ptmax,z,zmax1,zmax2
      real * 8 yseq(nymax),zseq(nptmax)
      real * 8 ener1,ener2,xmq,xlam,ffact,fren,pt,ylab
      real * 8 hdmassive,hdmassless,hdresummed,hdreserr,sigmafonll,
     #         phmassive,phmassless,phresummed,phreserr
      integer ny,npt,jy,jpt
c Unit used for temporary files. Keep closed
c before calling fonll0
      integer itm
      data itm/57/
      write(*,*) ' enter prefix for this run'
      read(*,*) prefix
      do jprefix=len(prefix),1,-1
         if(prefix(jprefix:jprefix).ne.' ') goto 11
      enddo
 11   continue
      open(unit=itm,file=prefix(1:jprefix)//'fonll.log')
      write(*,*) ' pdf type is described by 3 numbers'
      write(*,*) ' pdflib uses all 3; other may use less'
      write(*,*) ' All energy values are in GeV'
c
      write(*,*) ' Beam type are:'
      write(*,*) ' -1 for p_bar, 1 for p, 0 for (n+p)/2, 2 for n,'
      write(*,*) ' 3 for pi^+, -3 for pi^-'
      write(*,*) ' 4 for photon, 5 for electron'
c
      write(*,*) ' beam1 has positive rapidity'
c
      write(*,*)
      write(*,*) 'enter beam 1 type, energy, 3 pdf numbers'
      read(*,*)  ibeam1,ener1, nptype1, ngroup1, nset1
      write(itm,*) ibeam1,ener1, nptype1, ngroup1, nset1,
     # ' ! beam1: type, ener.,  nptype, ngroup, nset'
      write(*,*)
      write(*,*) 'enter beam 2 type, energy, 3 pdf numbers'
      read(*,*)  ibeam2,ener2, nptype2, ngroup2, nset2
      write(itm,*) ibeam2,ener2, nptype2, ngroup2, nset2,
     # ' ! beam2: type, ener.,  nptype, ngroup, nset'
c
      write(*,*)
      write(*,*) ' heavy quark mass'
      read(*,*) xmq
      write(itm,*)xmq,' ! heavy quark mass'
      write(*,*)
      write(*,*) 'LambdaQCD_5 in GeV, (<=0 for default'
      read(*,*) xlam
      write(itm,*) xlam, ' ! Lambda_5, <=0 for default'
c
      write(*,*)
      write(*,*) 'fact. and ren. scale factors'
      read(*,*) ffact,fren
      write(itm,*) ffact,fren,' !  fact. and ren. scale factors'
c
      write(*,*) 'enter rapidity sequence, 10000 to end'
      do ny=1,nymax
         read(*,*) yseq(ny)
         if(yseq(ny).eq.10000) goto 10
         write(itm,*) yseq(ny),'     ! y value'
      enddo
      write(*,*) ' too many rapidity points; increase nymax'
      stop
 10   continue
      write(itm,*) yseq(ny), '     ! y value to end sequence'
      ny=ny-1
      write(*,*) 'enter 0 if you want the pt sequence such that'
      write(*,*) '  mt=mtmin*exp((j-1)*log(mtmax/mtmin)/(npt-1))'
      write(*,*) 'anything else to enter the sequence youself'
      read(*,*) jptseqfl
      write(itm,*) jptseqfl, ' ! 0 for auto pt seq.'
      if(jptseqfl.eq.0) then
         write(*,*) 'enter maximum pt and # of pt points'
         read(*,*) ptmax,npt
         write(itm,*) ptmax,npt, ' ! pt,y_lab'
         if(npt.gt.nptmax) then
            write(*,*) 'too many pt points'
            stop
         endif
         do jpt=1,npt
            zseq(jpt)=dble(jpt-1)/(npt-1)
         enddo
      else
         write(*,*)'enter pt sequence, <0 to terminate'
         do npt=1,nptmax
            read(*,*) zseq(npt)
            write(itm,*) zseq(npt),'   ! pt points'
            if(zseq(npt).lt.0) goto 21
         enddo
         write(*,*) ' too many pt points'
         stop
 21      npt=npt-1
         ptmax=zseq(npt)
      endif
c     
      write(*,*)
      write(*,*) ' enter 1 for FONLL, 2 for Fixed Order, 3 for '//
     # 'resummed in the massless limit (no power suppressed mass terms)'
      read(*,*) icalctype
      write(itm,*) icalctype, '   ! icalctype'
      if(ibeam1.eq.5.or.ibeam2.eq.5) then
         write(*,*) ' Weiszaeker-Williams parameters'
         write(*,*)'enter 1 for cut in photon virtuality'
         write(*,*)'      2 to cut on angle'
         write(*,*)'      3 as in 1, only log term'
         write(*,*)'      4 as in 2, only log term'
         read(*,*) itypeww
         write(itm,*) itypeww,'    ! itype ww (1-4)'
         if(itypeww.eq.1.or.itypeww.eq.3)then
            write(*,*)'enter the effective WW scale in GeV (upper limit'
            write(*,*)'of the absolute value of the photon virtuality)'
            read(*,*) xmuww
            write(itm,*) xmuww, '   ! effective ww scale'
         elseif(itypeww.eq.2.or.itypeww.eq.4)then
            write(*,*)'theta_cut'
            read(*,*) thcww
            write(itm,*) thcww, '   ! angular cut for ww electron'
         endif
         write(*,*)
     #  ' enter minimum and maximum energy(photon)/energy(electron)'
         read(*,*) zminww,zmaxww
         write(itm,*) zminww,zmaxww,'! zminww,zmaxww'
      endif
      close(itm)
c End of input
      zmax1=1
      zmax2=1
      if(ibeam1.eq.5) zmax1=zmaxww
      if(ibeam2.eq.5) zmax2=zmaxww
c rapidity loop
      open(unit=12,file=prefix(1:jprefix)//'.out')
c on initial run
      write(12,'(7(d12.6,1x),i3)') ener1,zmax1,ener2,zmax2,
     #  ffact,xmq,ptmax,npt
      close(12)
      do jy=1,ny
         ylab=yseq(jy)
c pt loop
         do jpt=1,npt
            z=zseq(jpt)
            call ptpoints(ener1,zmax1,ener2,zmax2,
     #           ffact,ylab,xmq,ptmax,0,z,pt)
            call fonll0(icalctype,prefix,jprefix,
     #           ibeam1,ener1, nptype1, ngroup1, nset1,
     #           ibeam2,ener2, nptype2, ngroup2, nset2,
     #           xmq,xlam,ffact,fren,pt,ylab,
     #           hdmassive,hdmassless,phmassive,phmassless,
     #           hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
c write output
            open(unit=12,file=prefix(1:jprefix)//'.out')
            call toend(12)
            write(12,'(10(d12.6,1x),i3)')
     #           pt,ylab,hdmassive,hdmassless,hdresummed,hdreserr,
     #           phmassive,phmassless,phresummed,phreserr
            close(12)
            open(11,file=prefix(1:jprefix)//'.outlog')
c wind file to end
            call toend(11)
            write(11,111)ibeam1,ener1, nptype1, ngroup1, nset1,
     #           ibeam2,ener2, nptype2, ngroup2, nset2,
     #           xmq,xlam,ffact,fren,pt,ylab,sigmafonll,
     #           hdmassive,hdmassless,phmassive,phmassless,
     #           hdresummed,hdreserr,phresummed,phreserr
            write(11,*)
            close(11)
         enddo
      enddo
 111  format('beam1=',i2,',e1=',d12.6,',pdf1=',i1,',',i1,',',i3,
     #      ',beam2=',i2,',e2=',d12.6,',pdf2=',i1,',',i1,',',i3,/,
     #       'mq=',d12.6,',l5=',d12.6,',ff=',d12.6,',fr=',d12.6,/,
     #       'pt,y,fonll',3(1x,d12.6),/,
     #       'hdmv=',d12.6,',hdml=',d12.6,',phmv=',d12.6,',phml=',d12.6,
     # / ,'hdrs=',d12.4,'+-',d8.2,',phrs=',d12.6,'+-',d8.2)
      end
