************************************************************************
*                                                                      *
*                    FONLL             version 1.3.3 (Apr 2014)        *
*                                                                      *
*      Program to calculate heavy quark transverse momentum and        *
*      rapidity distributions in hadron-hadron and photon-hadron       *
*      collisions, matching Fixed Order next-to-leading order terms    *
*      and Next-to-Leading-Log large-p_T resummation.                  *
*                                                                      *
*      by Matteo Cacciari, Stefano Frixione, Paolo Nason               *
*                                                                      *
*      Citation for this work:                                         *
*                                                                      *
*      M. Cacciari, S. Frixione and P. Nason,                          *
*           "The p_T Spectrum in Heavy Flavor Photoproduction",        *
*           JHEP 0103 (2001) 006 [hep-ph/0102134]                      *
*      M. Cacciari, M. Greco and P. Nason,                             *
*           "The p_T Spectrum in Heavy Flavor Hadroproduction",        *
*           JHEP 9805 (1998) 007 [hep-ph/9803400]                      *
*                                                                      *
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
c verion number string
      include 'version.h'
c prefix for generated files
      character * 30 prefix
c jprefix is last non-blank in prefix
      integer jprefix
      integer ibeam1,ibeam2,nptype1,nptype2,ngroup1,ngroup2,nset1,nset2
     # ,icalctype
      real * 8 ener1,ener2,xmq,xlam,ffact,fren,pt,ylab
      real * 8 hdmassive,hdmassless,hdresummed,hdreserr,sigmafonll,
     #         phmassive,phmassless,phresummed,phreserr
      real*8 ratio,reserr
c ww common
      integer itypeww
      real * 8 xmuww,thcww,zminww,zmaxww,elphen
      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
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
      write(*,*)
      write(*,*) 'enter pt,y_lab'
      read(*,*) pt,ylab
      write(itm,*) pt,ylab, ' ! pt,y_lab'
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
      call fonll0(icalctype,prefix,jprefix,
     #  ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
c....used to check higher orders below
      if ( hdmassive+phmassive .gt. 0d0 ) then
          ratio = (hdresummed+phresummed-(hdmassless+phmassless))/
     #   (hdmassive+phmassive)
         reserr = hdreserr + phreserr
      else
          ratio = 0d0
          reserr = 0d0
      endif
      open(11,file=prefix(1:jprefix)//'.outlog',access='append')
c... commented out because a bug (#474382) in gfortran 4.4
c... produces an infinite loop
c wind file to end
c      call toend(11)
c... now using the line above with access='append', which seems to
c... work both in g77 and gfortran 4.4. Note however that the official
c... gfortran syntax is supposed to be
c...      open(11,file=prefix(1:jprefix)//'.outlog',position='append')
c...
      write(11,'(a,a,a)') versionstring
      write(11,111)ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,sigmafonll,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,ratio
      write(11,*)
      close(11)
 111  format('beam1=',i2,',e1=',d12.6,',pdf1=',i1,',',i1,',',i5,
     #      ',beam2=',i2,',e2=',d12.6,',pdf2=',i1,',',i1,',',i5,/,
     #       'mq=',d12.6,',l5=',d12.6,',ff=',d12.6,',fr=',d12.6,/,
     #       'pt,y,fonll',3(1x,d12.6),/,
     #       'hdmv=',d12.6,',hdml=',d12.6,',phmv=',d12.6,',phml=',d12.6,
     #       / ,'hdrs=',d12.4,'+-',d8.2,',phrs=',d12.6,'+-',d8.2,/,
     #       '(rs-ml)/mv=',d8.2)
      open(11,file=prefix(1:jprefix)//'.out')
c...check if higher orders are very large
c      if(dabs(ratio).gt.5d0)  then
c...check if numerical error on resummed is large
      if(dabs(reserr).gt.0.5*dabs(hdmassive+phmassive))  then
        write(11,'(1x,2(f14.6,1x),e18.9,2x,e15.6,2x,a)') 
     #                 pt,ylab,sigmafonll,reserr,'*'  
c        write(11,'(1x,2(f14.6,1x),e18.9,2x,e15.6,2x,a)') 
c     #                 pt,ylab,hdmassive+phmassive,reserr,'*'  
      else
        write(11,'(1x,2(f14.6,1x),e18.9,2x,e15.6)') 
     #                 pt,ylab,sigmafonll,reserr
      endif      
      close(11)
      write(*,*) ' matched result:', sigmafonll
      write(*,*) ' partial results:'
      write(*,*) ' hadronic mv, ml, rs:',
     # hdmassive,hdmassless,hdresummed
      if(    ibeam1.eq.4.or.ibeam1.eq.5
     #     .or. ibeam2.eq.4.or.ibeam2.eq.5 ) then
         write(*,*) ' pointlike photon mv, ml, rs:',
     #        phmassive, phmassless,phresummed
      endif
      end
