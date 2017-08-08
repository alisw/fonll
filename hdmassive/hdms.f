c Program to calculate heavy quarks hadroproduction.
c Written by Paolo Nason
c Modified by MC on 15-5-97 to evaluate also massless limit
c
c
c
      subroutine hdms(result)
c-------------------------------------------------------
c ********************** MAIN PROGRAM ******************
c-------------------------------------------------------
c
      implicit none
      real * 8 result
      integer nymx,nptmx
      parameter(nymx=50,nptmx=250)
      character * 70 file
      character * 6 task
      character * 9 vareq, tmpstr*30

      real * 8  xy(nymx),xpt(nptmx),a(nptmx,nymx),b(nptmx,nymx),
     # c(nptmx,nymx)
c-----------------------------------------------
c Common for
c            vlf=number of light flavors
c            vca=casimir adjoint
c            vcf=casimir fermion
c            vda=dimension adjoint
c            vtf=Tf (1/2 for fundamental)
c            vbf=(n^2-4)*(n^2-1)/16/n (for fundamental)
c            zg=coupling constant
c            xm=heavy quark mass
c            xmu=subtraction scale
c
      character * 2 beam*6,proc,asy
      integer nfl, lead
      real * 8 xmu2,xmuph2,xcsi,xcsiph,as
      common/hvqtrs/xmu2,xmuph2,xcsi,xcsiph,as,nfl,lead,proc,beam,asy
      real * 8 xlramu,ffact
      integer irunsc,istrsc
      common/rensca/xlramu,ffact,irunsc,istrsc
      integer ipt2topt
c------------------------------------------------------
c  common for massless limit
      real * 8 xm2true,rotrue,arotrue
      integer im0
      common/true/xm2true,rotrue,arotrue,im0
c------------------------------------------------------
c boost from CM to lab frame
      real * 8 ycm
      common/hera/ycm
c
      real * 8 xlam
      common/lambda/xlam
c common for flag nlfp1sch = 1 if we use the nlf+1 scheme
c for alfa, 0 otherwise
      integer nlfp1sch
      common/asnf/nlfp1sch
c set number for parton densities
      integer nset1,nset2
      common/pdfs/nset1,nset2
c local variables
      integer maxfcn,nptype1,ngroup1,nptype2,ngroup2,ny,j,npt,ill,k,
     # mufl,icmflag,lfile
      real * 8 e,s,fren,ptfix,xm,xm2,ro,ymax,xmtrue,pp2,xmu,xmure2,
     # y,tmp1,tmp2,tmp3,tmp4,etamin,etamax,xmhad1,xmhad2,phad1,phad2,
     # xmuph,ffactph,rscalpha,xlam1,xlam2
      character * 2 sche
      logical logfile
c functions
      real * 8 dfpy,alfas,dfpxf,dfp,dfy,ctpy,dfxf,ctxdfp
c-----------------------------------------------
c  Conversion factor from Gev**(-2) to microbarns.
c
      real * 8 xhc,pi
      integer iun7,iun8,iun9,iun10
      data xhc/389.3857/
      data pi/3.141 592 653 589 793/
      data iun7/7/,iun8/8/,iun9/9/,iun10/10/
c-----------------------------------------------
c On ibm only, inhibits underflow exceptions.
c
c      call xuflow(0)
c--------------------------------------------------
      logfile = .false.

c # of gaussian points
c
      if (logfile) open(unit=66,file='hdms-log.tmp',status='unknown')
      write(*,'(1x,a)')'enter maxfcn = # of gaussian points'
      read(55,*) maxfcn
      if (logfile) write(66,'(i3,10x,a)')maxfcn,'! # of gaussian points'
c In what follows, nptype is not used, except for getting Lambda_QCD and 
c scheme associated to the pdf set given by nset. Entering nptype=0 does
c not give problem when linking to mlmpdf, while a correct value should
c be used when linking to pdflib. In either cases, MLM labeling conventions
c for particles should be used. Ngroup is unused even when linking to
c pdflib; in that case, nset will be set equal to 100*#group+#set.
c See jetpdflib.f for details.
      write(*,*)
     # 'In what follows: nptype=1 -> proton, 4-> photon, 5-> electron'
      write(*,*)
     # 'ngroup is not used.'
      write(*,'(1x,a)')'enter nptype, ngroup, nset (beam1)'
      read(55,*) nptype1, ngroup1, nset1
      if (logfile) write(66,'(3(1x,i3),1x,a)') nptype1, ngroup1, nset1,
     #    '! nptype, ngroup, nset (beam1)'
      call setlam52(nset1,nptype1,xlam1,sche)
c In this code, only MSbar is implemented; '**' only comes when linking
c to PDFLIBs; in such a case, if a set given in DIS scheme is chosen, the
c result is erroneous, and no warning is given
      if(sche.eq.'**')then
        sche='MS'
c This value of xlam forces the code to ask for xlam to be used (not
c trusting the values given by PDFLIB)
        xlam1=-1.234d0
      endif
      if(sche.ne.'MS'.and.sche.ne.'DG')then
         write(*,*)'scheme not implemented'
         stop
      endif
      write(*,'(1x,a)')'enter nptype, ngroup, nset (beam2)'
      read(55,*) nptype2, ngroup2, nset2
      if (logfile) write(66,'(3(1x,i3),1x,a)') nptype2, ngroup2, nset2,
     #    '! nptype, ngroup, nset (beam2)'
      call setlam52(nset2,nptype2,xlam2,sche)
      if(sche.eq.'**')sche='MS'
      if(sche.ne.'MS'.and.sche.ne.'DG')then
         write(*,*)'scheme not implemented'
         stop
      endif
c One might want to set Lambda to a value different from that of the pdfs
      write(*,*)'enter Lambda_5 (0 for default)'
      read(55,*) tmp1
      if(tmp1.gt.0.d0) then
         xlam=tmp1
      elseif(xlam1.eq.xlam2)then
         xlam=xlam1
      else
         write(*,*)'cannot extablish default value for Lambda_5'
         write(*,*)'pdf set 1 wants',xlam1
         write(*,*)'pdf set 2 wants',xlam2
         write(*,*)'enter a specific value.'
         stop
      endif
c---------------------------------------------------
c Calculate values tables.
c
      write(*,'(1x,a,/,a)')
     # 'enter dfpy,dfpxf,ctpy,dfp,dfy,dfxf,ctxdfp for task selection,',
     # ' precede by - if you want rapidity by row'
c- read and edit (convert to lowercase and untab)
      call getstr(task)
      if (logfile) write(66,'(a,7x,a)') task,'! task'
      if(task(1:1).eq.'-') then
         vareq='y   pt'
         task = task(2:)
      else
         vareq='pt  y'
      endif
      if(task.eq.'dfpy') then
         write(*,*)' d sigma /( dy d (pT**2) )'
      elseif(task.eq.'dfpxf') then
         write(*,*)' d sigma /( dxf d (pT**2) )'
      elseif(task.eq.'ctpy') then
         write(*,*)' sigma( pt>pt_cut, y<y_cut)'
      elseif(task.eq.'dfp') then
         write(*,*)' d sigma/ d(pT**2) for y<y_cut'
      elseif(task.eq.'dfy') then
         write(*,*)' d sigma/ dy for pt>pt_cut'
      elseif(task.eq.'dfxf') then
         write(*,*)' d sigma/ d xf for pt>pt_cut'
         write(*,*)
     #   ' y value interpreted as xf values, xf = p_l/p_l_max in CM'
      elseif(task.eq.'ctxdfp') then
         write(*,*)' d sigma/ d(pt**2) for xf>xf_cut'
         write(*,*)
     #   ' y value interpreted as xf values, xf = p_l/p_l_max in CM'
      else
         write(*,*) 'undefined task'
         stop
      endif
      ipt2topt=0
      if(task.eq.'dfpy' .or. task.eq.'dfpxf' .or. 
     #   task.eq.'dfp' .or. task.eq.'dfp2')then
         write(*,*)' enter 1 to convert dsig/dpt^2 --> dsig/dpt'
         read(55,*)ipt2topt
         if (logfile) write(66,'(i1,12x,a)') ipt2topt,
     #        '! 1 to convert dsig/dpt^2 --> dsig/dpt' 
         if(ipt2topt.eq.1)
     #     write(*,*)' results will be expressed as dsig/dpt'
      endif
      write(*,'(1x,a)')'enter # of light flavours'
      read(55,*) nfl
      if (logfile) write(66,'(i1,12x,a)') nfl,'! # of light flavours'
      write(*,*)' enter 1 for the nlf+1  scheme, 0 for nlf scheme)'
      read(55,*) nlfp1sch
      if (logfile) write(66,'(i1,12x,a)') nlfp1sch,
     # '! 1 for the nlf+1  scheme, 0 for nlf scheme'
      write(*,*)' enter the rescaling factor for alpha_S'
      write(*,*)' alpha_S --> alpha_S/rscalpha'
      read(55,*)rscalpha
      if (logfile) write(66,'(d10.4,3x,a)') rscalpha,
     #   '! alpha_S --> alpha_S/rscalpha'
      write(*,'(1x,a)')
     # 'enter process and CM energy; instead of CM energy, one can',
     # 'also enter hadron momenta and hadron masses (enter m_had<0',
     # 'for default values). Valid entries are:',
     # 'pr+pr- 1800  (proton antiproton at E_cm=1800)',
     # 'pr+nu+ 100 0 -1 -1 (100 GeV protons on nucleons (i.e. (p+n)/2)',
     # '                    masses are set to 0.938 GeV)',
     # 'pr+nu+ 100 0 0 -1 (100 GeV protons on nucleons (i.e. (p+n)/2 )',
     # '     kinematics of the incoming hadron in the massless limit)',
     # 'ph pr+ 200   (photon-proton at E_cm=200 GeV)',
     # 'ph pr+ 20 820 0 (photon-proton with k_gamma=20 GeV)',
     # '                 k_pr=820 GeV, and neglecting proton mass)',
     # 'ph nu+ 30    (photon-nucleon ( i.e. (p+n)/2 ) at E_cm=30 GeV)',
     # 'ph nu+ 100 0 -1 -1 (100 GeV photon on nucleon at fixed target,',
     # '                    with m_ph=0, m_nu=0.938 GeV)',
     # 'el pr+ 300   ( electron proton at E_cm=300 GeV)'
      call getstr(tmpstr)
      if(tmpstr(4:4).eq.' ') then
c if one enters pr+ pr-, get rid of blank in 4rth position
         tmpstr=tmpstr(1:3)//tmpstr(5:30)
      endif
      if (logfile) write(66,'(a,a)') tmpstr(1:30),
     #                               '! process and CM energy'
      beam = tmpstr(1:6)
      read(tmpstr(7:),fmt=*,err=512,end=512) phad1,phad2,xmhad1,xmhad2
      if(xmhad1.lt.0) then
         if(beam(1:2).eq.'pr'.or.beam(1:2).eq.'nu') then
            xmhad1=0.938
         elseif(beam(1:2).eq.'pi') then
            xmhad1=0.140
         elseif(beam(1:2).eq.'ph'.or.beam(1:2).eq.'el') then
            xmhad1=0
         else
            write(*,*) ' error, beam',beam,'invalid'
            stop
         endif
      endif
      if(xmhad2.lt.0) then
         if(beam(1:2).eq.'pr'.or.beam(1:2).eq.'nu') then
            xmhad2=0.938
         elseif(beam(1:2).eq.'pi') then
            xmhad2=0.140
         elseif(beam(1:2).eq.'ph'.or.beam(1:2).eq.'el') then
            xmhad2=0
         else
            write(*,*) ' error, beam',beam,'invalid'
            stop
         endif
      endif
      write(*,*) ' hadron masses:',xmhad1,xmhad2
      e=sqrt( 2*sqrt(xmhad1**2+phad1**2)*sqrt(xmhad2**2+phad2**2)
     #       +xmhad1**2+xmhad2**2+2*phad1*phad2 )
      goto 513
 512  read(tmpstr(7:),fmt=*) e
 513  continue
      write(*,*) beam, e
      s = e**2
c Find now if rapidity is relative to CM or lab
      write(*,*) ' enter 1,2 or 3 to specify frame:'
      write(*,*) ' 1 for rapidity relative to CM of incoming hadrons'
      write(*,*) ' 2 to CM of hadron 2 (i.e., fixed target)'
      write(*,*) ' 3 in lab. (according to momenta entered earlier)'
      read(55,*) icmflag
      if (logfile) write(66,'(i1,12x,a)') icmflag,
     # '! 1 for CM, 2 for FT, 3 for lab, y and eta'
      if(icmflag.eq.1) then
         ycm=0
      elseif(icmflag.eq.2) then
c beam energy (in the lab frame):
         tmp1=(s-xmhad1**2-xmhad2**2)/(2*xmhad2)
c beam momentum (in the lab frame):
         tmp2=sqrt(tmp1**2-xmhad1**2)
         ycm=0.5d0*
     #     log( ((tmp1+xmhad2)+tmp2)/((tmp1+xmhad2)-tmp2) )
      elseif(
c (beam(1:2).eq.'ph'.or.beam(1:2).eq.'el').and.
     #              icmflag.eq.3) then
         tmp1=sqrt(phad1**2+xmhad1**2)+sqrt(phad2**2+xmhad2**2)
         tmp2=phad1-phad2
         ycm=0.5d0*log( (tmp1+tmp2) / (tmp1-tmp2) )
         write(*,*) '!!!!!!!!!!***** WARNING! *****!!!!!!!!!!!'
         write(*,*)
     #   ' in our convention hadron 1 has positive rapidity'
      else
         write(*,*) ' This option is not implemented'
         stop
      endif
      write(*,'(1x,a,2(/,a))')
     # 'enter gg, qq, qg, any other for all',
     # 'to compute (sigma(Q)+sigma(Qbar))/2 from the given channel',
     # '(-,-qq,-qg for total and partial asymmetry='
     # //'(sigma(Q)-sigma(Qbar))/2'
      call getstr(tmpstr)
      if (logfile) write(66,'(a,a)') tmpstr(1:13),
     #                               '! gg, qq, qg, any other for all'
      if(tmpstr(1:1).eq.'-') then
         proc = tmpstr(2:)
         asy  = 'as'
         write(*,*)'will return (sigma(Q)-sigma(Q_bar))/2'
      else
         write(*,*)'will return (sigma(Q)+sigma(Q_bar))/2'
         proc = tmpstr
      endif
      write(*,*)' enter 1 for reference scale mu0=sqrt(m^2+pt^2)'
      write(*,*)'       2 for reference scale mu0=pt'
      write(*,*)'       3 to be prompted for sequence in mu0'
      read(55,*) mufl
      if(mufl.ne.1.and.mufl.ne.2.and.mufl.ne.3) then
         write(*,*) ' must between 1 and 3!'
         stop
      endif
      if (logfile) write(66,'(i1,12x,a)') mufl,
     #   '! 1 mu0=sqrt(m^2+pt^2), 2 mu0=pt, 3 sequence (in mu0)'
c xmu_fact is the factorization scale of leg # 2 (proton at HERA in
c our conventions); see below for scale relevant to leg # 1
      write(*,*)' enter ffact, xmu_fact=ffact*mu0'
      read(55,*) ffact
      if (logfile) write(66,'(d10.4,3x,a)') ffact,'! xmu_fact=ffact*mu0'
      write(*,*) ' ffact =', ffact
      write(*,*) ' enter 1 to keep ffact=1 only in structure functions'
      write(*,*) ' (for debugging purposes only), 0 otherwise'
      read(55,*) istrsc
      if (logfile) write(66,'(i1,12x,a)') istrsc, 
     #                                    '! 0 normal, 1 for tests'
      ffactph=ffact
c Set the factorization scale on leg # 1 equal to that of leg # 2;
c in a more refined version, treat them independently, and uncomment
c the following lines
c      write(*,*)' enter ffactph, xmu_fact(photon)=ffactph*mu0'
c      read(55,*) ffactph
c      if (logfile) write(66,'(d10.4,3x,a)') ffactph,
c     #   '! xmu_fact(photon)=ffactph*mu0'
c      write(*,*) ' ffactph =', ffactph
c      if(ffactph.ne.ffact)then
c         write(*,*)'ffactph # ffact: option not implemented'
c         stop
c      endif
      write(*,*) 'enter fren, xmu_ren=fren*mu0'
      read(55,*)  fren
      if (logfile) write(66,'(d10.4,3x,a)') fren, '! xmu_ren=fren*mu0'
      write(*,*) ' fren =', fren
c      write(*,*)
c     #' enter 1 if you want mu to be recomputed for each pt,y point'
c      read(55,*) irunsc
c------------------------------------
c xlramu is always log(xmu_ren^2/xmu_fact^2)
c
      xlramu = 2*log(fren/ffact)
      if(mufl.eq.3) then
        write(*,*)' enter pt value'
        read(55,*)  ptfix
        if (logfile) write(66,'(d10.4,3x,a)') ptfix, '! pt value'
      endif
      write(*,'(1x,a)')
     #'enter a mass value, to set rapidity and pt values to',
     #'0,.2,.4,.6,.8*ymax, pt=.02,.2,.3,.6,1,1.2,1.5,2.2.5,3 * mass',
     #'enter a negative value to be prompted for y and pt sequence'
      read(55,*) xm
      if (logfile) write(66,'(d10.4,3x,a)') xm,
     # '! mass value, <0 to set rapidity and pt values'
      if(xm.ge.0) then
        xm2 = xm**2
        ro = 4*xm**2/s
c-----------------------------------------------
c Defaults.
c
        ymax = dlog( 1/dsqrt(ro/4)/2+dsqrt(1/(ro/4)/4-1) )
        ny = 5
c These are laboratory rapidities 
        xy ( 1 ) = 0.d0         + ycm
        xy ( 2 ) = 0.2d0 * ymax + ycm
        xy ( 3 ) = 0.4d0 * ymax + ycm
        xy ( 4 ) = 0.6d0 * ymax + ycm
        xy ( 5 ) = 0.8d0 * ymax + ycm
        write(*,'(1x,a,d10.4,/,7(1x,d10.4))')'ymax=',ymax,(xy(j),j=1,ny)
        npt = 11
        xpt ( 1 ) = (.02 * xm)
        xpt ( 2 ) = (0.2 * xm)
        xpt ( 3 ) = (0.4 * xm)
        xpt ( 4 ) = (0.6 * xm)
        xpt ( 5 ) = (0.8 * xm)
        xpt ( 6 ) = (1.0 * xm)
        xpt ( 7 ) = (1.2 * xm)
        xpt ( 8 ) = (1.5 * xm)
        xpt ( 9 ) = (2.0 * xm)
        xpt ( 10) = (2.5 * xm)
        xpt ( 11) = (3.0 * xm)
c      xpt ( 1 +11) = (3.5 * xm)
c      xpt ( 2 +11) = (4   * xm)
c      xpt ( 3 +11) = (4.5 * xm)
c      xpt ( 4 +11) = (5   * xm)
c      xpt ( 5 +11) = (5.5 * xm)
c      xpt ( 6 +11) = (6   * xm)
c      xpt ( 7 +11) = (6.5 * xm)
c      xpt ( 8 +11) = (7.0 * xm)
c      xpt ( 9 +11) = (7.5 * xm)
c      xpt ( 10+11) = (8.0 * xm)
c      xpt ( 11+11) = (8.5 * xm)
c      xpt ( 12+11) = (9.0 * xm)
c      xpt ( 13+11) = (9.5 * xm)
      else
        write(*,*)'enter y value(s), * to terminate the sequence'
        write(*,*)'NB: rapidities must be given in the LAB FRAME'
        do j=1,nymx
           read(55,'(a)',err=133) tmpstr
           if (logfile) write(66,'(a)') tmpstr
           if(index(tmpstr,'*').ne.0) goto 133
           read(tmpstr,*) xy(j)
        enddo
        write(*,*)'no more values.'
 133    ny = j-1
        if(mufl.ne.3)then
           write(*,*)'enter pt value(s), * to terminate the sequence'
        else
           write(*,*)'enter mu0 value(s), * to terminate the sequence'
        endif
        do j=1,nptmx
        read(55,'(a)',err=134) tmpstr
        if (logfile) write(66,'(a)') tmpstr
        if(index(tmpstr,'*').ne.0) goto 134
        read(tmpstr,*) xpt(j)
        enddo
        write(*,*)'no more values.'
 134    npt = j-1
      endif
c
c....modifica mia qui, MC
c
 211  write(*,'(1x,a)')
     #' enter mass of heavy quark followed by an integer=',
     #' 0 for massive, 1 for massless, 2 for massless leading log'
      read(55,*) xmtrue, im0
      if(im0.eq.0) then
         ill = 0
         xm = xmtrue
      elseif(im0.eq.1) then
         ill = 0
         xm = 0
      elseif(im0.eq.2) then
         im0 = 1
         xm = 0
         ill = 1
      else
         write(*,*) ' enter 0,1 or 2'
         goto 211
      endif
      if (logfile) write(66,'(d10.4,1x,i1,1x,a)') xmtrue, im0,
     #'! mass + flag, 0 for massive, 1 for massless, 2 for massless LL'
      if(im0.ne.0) then
         write(*,*) ' enter 1 if you want pt -> sqrt(pt^2+m^2)'
         read(55,*) j
         if (logfile) write(66,'(i2,11x,a)')j,
     #                             '!1 if you want pt -> sqrt(pt^2+m^2)'
         if(j.eq.1) then
            do j=1,npt
               xpt(j)=sqrt(xpt(j)**2+xmtrue**2)
            enddo
         endif
      endif
      xm2true = xmtrue**2
      rotrue = 4*xm2true/s
      xm2 = xm**2
c
      write(*,'(1x,a)')'enter filename to go, . to ouput to screen'
      call getstr(file)
      if (logfile) write(66,'(a)') file
      file = file(1:index(file,' '))
      write(*,'(1x,a)')file
      if(file.eq.'.') then
        iun7 = 6
        iun8 = 6
        iun9 = 6
        lfile = 0
      else
        do k=len(file),1,-1
           if(file(k:k).ne.' ')goto 235
        enddo
 235    continue
        open (iun7,file = file(1:k),status='UNKNOWN')
c,recl=132)
        open (iun8,file = file(1:k)//'_lo',status='UNKNOWN')
c,recl=132)
        open (iun9,file = file(1:k)//'_ra',status='UNKNOWN')
c,recl=132)
        lfile=k
      endif
c-------------------------------------
c       Max tp2mx is (1-ro)/4
c
c------------------------------------------------------------
c  Y is the heavy quark rapidity (CM frame), tp2 is pt^2/s;
c  we have:
c  t1*t2=tp2+ro/4
c  y=1/2*log(t2/t1)
c
c  therefore: tp2 < (1-ro)/4, equivalent to pt^2+m^2 < e^2
c
c  Define t12=t1*t2=tp2+ro/4.
c  From the condition t1+t2 < 1 we get:
c
c             Y > log(1/2/sqrt(t12)-sqrt(1/4/t12-1))
c             Y < log(1/2/sqrt(t12)+sqrt(1/4/t12-1))
c
c-------------------------------------
c Nested do loop; p_perp loop: ( or mu loop, if mufl = 3 )
c
      write(*,'(1x,a,d10.4)') 'xlam = ',xlam,' cme =',e,' m=',xm
c Close dialog log
      if (logfile) close(66)
c
      do 2 j = 1,npt
        pp2  = (xpt ( j ))**2
c When pt=0, the cross section diverges due to the 1/v singularity
        if(pp2.eq.0.d0)then
          write(*,*)'At p_T=0 the cross section is divergent at NLO'
          stop
        endif
c----------------------------------------------------
c       Set value of alpha for scale xmu;
c
        if(mufl.eq.3) then
          xmu = xpt(j) * ffact
          xmuph = xpt(j) * ffactph
          pp2 = ptfix**2
        elseif(mufl.eq.1)then
          xmu = sqrt(xm2+pp2) * ffact
          xmuph = sqrt(xm2+pp2) * ffactph
        elseif(mufl.eq.2)then
          xmu = sqrt(pp2) * ffact
          xmuph = sqrt(pp2) * ffactph
        endif
        xmu2 = xmu*xmu
        xmuph2 = xmuph*xmuph
        xcsi = log(xmu2/xm2true)
        xcsiph = log(xmuph**2/xm2true)
        xmure2 = xmu2*exp(xlramu)
c
        if(nlfp1sch.eq.1) then
           as   = alfas(xmure2,xlam,nfl+1)
           as   = as/rscalpha
        elseif(nlfp1sch.eq.0) then
           as   = alfas(xmure2,xlam,nfl)
           as   = as/rscalpha
        else
           write(*,*) ' error: flag for nlf+1 scheme should be 1 or 0'
           stop
        endif
c
        write(*,*) ' alfa=', as, 'at scale ',sqrt(xmure2)
c-------------------------------------------------------------------
c    Rapidity loop;
c
        do 1 k=1,ny
           y = xy(k)-ycm
           if(task.eq.'dfpy') then
              lead = 1
              b(j,k) = dfpy(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = dfpy(pp2,y,xm2,s,maxfcn)
              if(ill.eq.1) then
                 tmp1 = xm2true
                 tmp2 = rotrue
                 tmp3 = xcsi
                 tmp4 = xcsiph
                 xcsi = 0
                 xcsiph = 0
                 xm2true = xmu2
                 rotrue = 4*xm2true/s
                 a(j,k)=a(j,k)-dfpy(pp2,y,xm2,s,maxfcn)+b(j,k)
                 xcsiph = tmp4
                 xcsi = tmp3
                 rotrue = tmp2
                 xm2true = tmp1
              endif
           elseif(task.eq.'dfpxf') then
              lead = 1
              b(j,k) = dfpxf(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = dfpxf(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'dfp') then
              lead = 1
              b(j,k) = dfp(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = dfp(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'dfy') then
              lead = 1
              b(j,k) = dfy(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = dfy(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'ctpy') then
              lead = 1
              b(j,k) = ctpy(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = ctpy(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'dfxf') then
              lead = 1
              b(j,k) = dfxf(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = dfxf(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'ctxdfp') then
              lead = 1
              b(j,k) = ctxdfp(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = ctxdfp(pp2,y,xm2,s,maxfcn)
           endif
c------------------------------------------------------
c Normalize; xhc is the plank constant times the speed
c of light in Gev^2*microbarns; we have to divide
c by s**2, because in the functions ssp,ssd this
c factor is omitted, and multiply by pi for d^2 pt -> d (pt^2):
c
           a(j,k) = a(j,k)*pi*xhc/s**2
           b(j,k) = b(j,k)*pi*xhc/s**2
           if(b(j,k).ne.0) then
              c(j,k) = a(j,k)/b(j,k)
           else
              c(j,k) = 1
           endif
 1         continue
c---- Print the values.
        if(ipt2topt.ne.1)then
           write( 6,'(10(1x,a,1x,d14.8,/))')
     #          vareq ,xpt(j),(' tot=',a(j,k),' lo=',b(j,k),
     #          ' rat=',c(j,k),k=1,ny)
        else
           write( 6,'(10(1x,a,1x,d14.8,/))')
     #          vareq ,xpt(j),(' tot=',2*xpt(j)*a(j,k),
     #          ' lo=',2*xpt(j)*b(j,k),' rat=',c(j,k),k=1,ny)
        endif
    2 continue
      if (vareq(1:2).eq.'y ') then
        if(index(task,'x').ne.0) vareq(1:3)='xf'
        if(mufl.eq.3) vareq(5:)='xmu'
        write( iun7,'(1x,1h!,a,10(1x,e14.8))') vareq,(xpt(k),k=1,npt)
        write( iun8,'(1x,1h!,a,10(1x,e14.8))') vareq,(xpt(k),k=1,npt)
        write( iun9,'(1x,1h!,a,10(1x,e14.8))') vareq,(xpt(k),k=1,npt)
        do k=1,ny
           if(ipt2topt.ne.1)then
              write( iun7,'(10(1x,e14.8))') xy(k),(a(j,k),j=1,npt)
              write( iun8,'(10(1x,e14.8))') xy(k),(b(j,k),j=1,npt)
           else
              write( iun7,'(10(1x,e14.8))') xy(k),
     #                                     (2*xpt(j)*a(j,k),j=1,npt)
              write( iun8,'(10(1x,e14.8))') xy(k),
     #                                     (2*xpt(j)*b(j,k),j=1,npt)
           endif
           write( iun9,'(10(1x,e14.8))') xy(k),(c(j,k),j=1,npt)
        enddo
      else
        if(index(task,'x').ne.0) vareq(5:)='xf'
        if(mufl.eq.3) vareq(1:3)='xmu'
        write( iun7,'(1x,1h!,a,10(1x,e14.8))') vareq,(xy(k),k=1,ny)
        write( iun8,'(1x,1h!,a,10(1x,e14.8))') vareq,(xy(k),k=1,ny)
        write( iun9,'(1x,1h!,a,10(1x,e14.8))') vareq,(xy(k),k=1,ny)
        do j=1,npt
           if(ipt2topt.ne.1)then
              write( iun7,'(10(1x,e14.8))') xpt(j),(a(j,k),k=1,ny)
              write( iun8,'(10(1x,e14.8))') xpt(j),(b(j,k),k=1,ny)
           else
              write( iun7,'(10(1x,e14.8))') xpt(j),
     #                                     (2*xpt(j)*a(j,k),k=1,ny)
              write( iun8,'(10(1x,e14.8))') xpt(j),
     #                                     (2*xpt(j)*b(j,k),k=1,ny)
           endif
           write( iun9,'(10(1x,e14.8))') xpt(j),(c(j,k),k=1,ny)
        enddo
      endif
      if(lfile.gt.0) then
         open(iun10, form='unformatted',
     #        file=file(1:lfile)//'_raw',status='UNKNOWN')
c unformatted write of all data
         write(iun10) npt,ny
         write(iun10) (xpt(k),k=1,npt)
         write(iun10) (xy(k),k=1,ny)
         write(iun10) ((a(j,k),j=1,npt),k=1,ny)
         write(iun10) ((b(j,k),j=1,npt),k=1,ny)
      endif
      result=a(1,1)
      end

c--------------------------------------------------
c returns d sigma/dy (*s*4*2pi^2)
c with pt cut: pt2>pp2
c
      function dfy(pp2ct,y,xm2,s,maxfcn)
      implicit none
      real * 8 dfy,pp2ct,y,xm2,s
      integer maxfcn
      integer maxgau
      real * 8 par
      parameter (maxgau=160,par=.5d0)
      real * 8 xg(maxgau),wg(maxgau)
c local
      real * 8 t2ot1,p2mx,xint,v,pp2,xjac
      integer j
c functions
      real * 8 dfpy
c
      real * 8 vmn,vmx,zero,one
      integer ngau
      data vmn/0.d0/,vmx/1.d0/,ngau/0/
      data zero,one/0.d0,1.d0/
      t2ot1 = exp(2*y)
      p2mx = s*t2ot1/(1+t2ot1)**2 - xm2
c also if p2mx<0:
      if(pp2ct.ge.p2mx) then
        dfy = 0
        return
      endif
      if(ngau.eq.0) then
         if(maxfcn.gt.maxgau) then
            write(*,*)
     #      'maxfcn too large: increase maxgau and recompile'
            stop
         endif
         ngau = maxfcn
         call gauleg(zero,one,xg,wg,ngau)
      endif
c--- integration loop
      xint = 0
      vmn = ((pp2ct+xm2)/(xm2+p2mx))**(1/par)
      do 1 j=1,ngau
         v    = vmn + (vmx-vmn)*xg(j)
         pp2  = (1/v**par - 1)*(xm2+pp2ct) + pp2ct
         xjac = (vmx-vmn)*(xm2+pp2ct) * par/ v**(par+1)
         xint = xint+xjac*wg(j)*dfpy(pp2,y,xm2,s,maxfcn)
 1       continue
      dfy = xint
      return
      end

c--------------------------------------------------
c returns d sigma/d p2 (*s*4*2pi^2)
c with y cut: y<yct
c
      function dfp(pp2,yct,xm2,s,maxfcn)
      implicit none
      real * 8 dfp,pp2,yct,xm2,s
      integer maxfcn
      integer maxgau
      parameter (maxgau=160)
      real * 8 xg(maxgau),wg(maxgau)
c local
      real * 8 t12,ymax,yl,xint,y,xjac
      integer j
c function
      real * 8 dfpy
      real * 8 zero,one
      integer ngau
      data zero,one/0.d0,1.d0/,ngau/0/
      if(pp2.ge.s/4-xm2) then
        dfp = 0
        return
      endif
      if(ngau.eq.0) then
         if(maxfcn.gt.maxgau) then
            write(*,*)
     #      'maxfcn too large: increase maxgau and recompile'
            stop
         endif
         ngau = maxfcn
         call gauleg(zero,one,xg,wg,ngau)
      endif
c--- integration loop
c
c.....questo ro non e' mai usato ne' passato in common,
c.....quindi si puo' magari anche togliere
c
c      ro   = 4*xm2/s
c
      t12  = (pp2+xm2)/s
      ymax = dlog( 1/dsqrt(t12)/2+dsqrt(1/t12/4-1) )
      yl   = min(ymax,yct)
      xint = 0
      do 1 j=1,ngau
         y    = yl*(1-2*xg(j))
         xjac = 2*yl
         xint = xint+xjac*dfpy(pp2,y,xm2,s,maxfcn)*wg(j)
 1       continue
      dfp = xint
      return
      end

c--------------------------------------------------
c returns sigma (*s*4*2pi^2)
c with pt cut: pt2>pp2ct
c and  y  cut: y<yct
c
      function ctpy(pp2ct,yct,xm2,s,maxfcn)
      implicit none
      real * 8 ctpy,pp2ct,yct,xm2,s
      integer maxfcn
      integer maxgau
      real * 8 par
      parameter (maxgau=160,par=.5d0)
      real * 8 xg(maxgau),wg(maxgau)
c local
      real * 8 ro,t12,ymax,yl,xint,y,weighty,t2ot1,p2mx,v,
     #     pp2,weightp
      integer k,j
c functions
      real * 8 dfpy
      real * 8 vmn,vmx,zero,one
      integer ngau
      data vmn/0.d0/,vmx/1.d0/,ngau/0/
      data zero,one/0.d0,1.d0/
      if(pp2ct.ge.s/4-xm2) then
        ctpy = 0
        return
      endif
      if(ngau.eq.0) then
         if(maxfcn.gt.maxgau) then
            write(*,*)
     #      'maxfcn too large: increase maxgau and recompile'
            stop
         endif
         ngau = maxfcn
         call gauleg(zero,one,xg,wg,ngau)
      endif
c--- integration loops
      ro   = 4*xm2/s
c in old version was (pp2+xm2)/s; bug?
      t12  = (pp2ct+xm2)/s
      ymax = dlog( 1/dsqrt(t12)/2+dsqrt(1/t12/4-1) )
      yl   = min(ymax,yct)
c--- y loop
      xint = 0
      do 1 k=1,ngau
         y  = yl*(1-2*xg(k))
         weighty = yl * 2 * wg(k)
         t2ot1 = exp(2*y)
         p2mx = s*t2ot1/(1+t2ot1)**2 - xm2
         vmn = ((pp2ct+xm2)/(xm2+p2mx))**(1/par)
         do 1 j=1,ngau
            v    = vmn + (vmx-vmn)*xg(j)
            pp2  = (1/v**par - 1)*(xm2+pp2ct) + pp2ct
            weightp = wg(j)*(vmx-vmn)*(xm2+pp2ct)*par/v**(par+1)
            weightp = weightp*dfpy(pp2,y,xm2,s,maxfcn)
c            write(*,*)'weightp=',weightp
            xint = xint+weighty*weightp
 1          continue
      ctpy = xint
      return
      end

      function dfxf(pp2ct,xf,xm2,s,maxfcn)
      implicit none
      real * 8 dfxf,pp2ct,xf,xm2,s
      integer maxfcn
      integer maxgau
      real * 8 par
      parameter (maxgau=160,par=.5d0)
      real * 8 xg(maxgau),wg(maxgau)
c local
      real * 8 ro,b,xfb,pp2mx,xint,v,pp2,t1pt2,y,xjac
      integer j
c functions
      real * 8 dfpy
      real * 8 vmn,vmx,zero,one
      integer ngau
      data vmn/0.d0/,vmx/1.d0/,ngau/0/
      data zero,one/0.d0,1.d0/
      ro = 4*xm2/s
      b  = sqrt(1-ro)
      xfb = xf*b
      pp2mx = s*(1+xfb)*(1-xfb)/4-xm2
      if(pp2ct.gt. pp2mx) then
         dfxf = 0
         return
      endif
      if(ngau.eq.0) then
         if(maxfcn.gt.maxgau) then
            write(*,*)
     #      'maxfcn too large: increase maxgau and recompile'
            stop
         endif
         ngau = maxfcn
         call gauleg(zero,one,xg,wg,ngau)
      endif
c--- integration loop
      xint = 0
      vmn = ((pp2ct+xm2)/(xm2+pp2mx))**(1/par)
      do 1 j=1,ngau
         v    = vmn + (vmx-vmn)*xg(j)
         pp2  = (1/v**par - 1)*(xm2+pp2ct) + pp2ct
         t1pt2 = sqrt(xfb**2+4*(pp2+xm2)/s)
         y     = log( (t1pt2+xfb)/(t1pt2-xfb) )/2
c------------------------------------------------------------------
c The extra factor b/(t1+t2) is to go from d sigma/dy/dpp2 to
c d sigma /dxf/dpp2
c
         xjac = b/t1pt2
         xjac = xjac*(vmx-vmn)*(xm2+pp2ct) * par/ v**(par+1)
         xint = xint+xjac*wg(j)*dfpy(pp2,y,xm2,s,maxfcn)
 1       continue
      dfxf = xint
      return
      end

c--------------------------------------------------
c returns d sigma/d p2 (*s*4*2pi^2)
c with xf cut: xf<xfct
c
      function ctxdfp(pp2,xfct,xm2,s,maxfcn)
      implicit none
      real * 8 ctxdfp,pp2,xfct,xm2,s
      integer maxfcn
      integer maxgau
      parameter (maxgau=160)
      real * 8 xg(maxgau),wg(maxgau)
c local
      real * 8 ro,b,t12,xfbmax,xfbct,xint,xfb,t1pt2,y,xjac
      integer j
c functions
      real * 8 dfpy
      real * 8 zero,one
      integer ngau
      data zero,one/0.d0,1.d0/,ngau/0/
      if(pp2.ge.s/4-xm2) then
        ctxdfp = 0
        return
      endif
      if(ngau.eq.0) then
         if(maxfcn.gt.maxgau) then
            write(*,*)
     #      'maxfcn too large: increase maxgau and recompile'
            stop
         endif
         ngau = maxfcn
         call gauleg(zero,one,xg,wg,ngau)
      endif
c--- integration loop
      ro   = 4*xm2/s
      b    = sqrt(1-ro)
      t12  = (pp2+xm2)/s
      xfbmax = sqrt(1-4*t12)
      xfbct  = xfct*b
      if(xfbmax.le.xfbct) then
         ctxdfp = 0
         return
      endif
      xint = 0
      do 1 j=1,ngau
         xfb    = (xfbmax-xfbct)*xg(j) + xfbct
         t1pt2 = sqrt(xfb**2+4*t12)
         y     = log( (t1pt2+xfb)/(t1pt2-xfb) )/2
c------------------------------------------------------------------
c The extra factor 1/(t1+t2) is to go from d sigma/dy/dpp2 to
c d sigma /dxfb/dpp2
c
         xjac = (xfbmax-xfbct) * wg(j) / t1pt2
         xint = xint+xjac*dfpy(pp2,y,xm2,s,maxfcn)
 1       continue
      ctxdfp = xint
      return
      end

*       function dfpy(pp2,y,xm2,s,maxfcn)
*       implicit none
*       real * 8 dfpy,pp2,y,xm2,s
*       integer maxfcn
*       integer maxgau
*       parameter (maxgau=160)
*       real * 8 xg(maxgau),wg(maxgau)
*       character * 2 beam*6,proc,asy
*       real * 8 zmin,zmax
*       common/zlim/zmin,zmax
*       real * 8 xmu2,xmuph2,xcsi,xcsiph,as
*       integer nfl,lead
*       common/hvqtrs/xmu2,xmuph2,xcsi,xcsiph,as,nfl,lead,proc,beam,asy
*       real * 8 xlramu,ffact
*       integer irunsc,istrsc
*       common/rensca/xlramu,ffact,irunsc,istrsc
*       real * 8 xlam
*       common/lambda/xlam
*       real * 8 xm2true,rotrue,arotrue
*       integer im0
*       common/true/xm2true,rotrue,arotrue,im0
* c local
*       real * 8 ro,t12,t2,t1,xint,xintk
*       integer j,k
* c functions
*       real * 8 ssd,ssp
*       real * 8 zero,one
*       integer ngau
*       data zero,one/0.d0,1.d0/,ngau/0/
* c
* c....queste scale, e soprattutto xcsi, sono da lasciare massive
* c....modifica mia, MC
* c
* c      if(irunsc.eq.1) then
* c         xmu2    = (pp2+xm2true)*ffact**2
* c         xcsi    = log(xmu2/xm2true)
* c         xmure2 = xmu2*exp(xlramu)
* c         as      = alfas(xmure2,xlam,nfl+1)
* c      endif
* c
* c....questi ro,t12,t2 danno i limiti di sp. fasi, saranno da mettere massless
* c
*       ro   = 4*xm2/s
*       t12  = (pp2+xm2)/s
*       t2  = sqrt(exp(2*y)*t12)
*       t1  = t12/t2
*       if(t1+t2.gt.1) then
*          dfpy = 0
*          return
*       endif
*       if(ngau.eq.0) then
*          if(maxfcn.gt.maxgau) then
*             write(*,*)
*      #      'maxfcn too large: increase maxgau and recompile'
*             stop
*          endif
*          ngau = maxfcn
*          call gauleg(zero,one,xg,wg,ngau)
*       endif
*       xint = 0
*       do 1 j=1,ngau
*         xint = xint + ssd(xg(j),ro,t1,t2)*wg(j)
*    1  continue
*       if(lead.ne.1) then
*         do 2 j=1,ngau
*            xintk = 0
*            do 3 k=1,ngau
*                xintk = xintk + ssp(xg(j),xg(k),ro,t1,t2)*wg(k)
*    3           continue
*            xint = xint + xintk*wg(j)
*    2       continue
*       endif
*       dfpy = xint
*       return
*       end

      function dfpy(pp2,y,xm2,s,maxfcn)
      implicit none
      real * 8 dfpy,pp2,y,xm2,s
      integer maxfcn
      integer maxgau
      parameter (maxgau=5000)
      real * 8 xg(maxgau),wg(maxgau)
      character * 2 beam*6,proc,asy
      real * 8 zmin,zmax
      common/zlim/zmin,zmax
      real * 8 xmu2,xmuph2,xcsi,xcsiph,as
      integer nfl,lead
      common/hvqtrs/xmu2,xmuph2,xcsi,xcsiph,as,nfl,lead,proc,beam,asy
      real * 8 xlramu,ffact
      integer irunsc,istrsc
      common/rensca/xlramu,ffact,irunsc,istrsc
      real * 8 xlam
      common/lambda/xlam
      real * 8 xm2true,rotrue,arotrue
      integer im0
      common/true/xm2true,rotrue,arotrue,im0
c local
      real * 8 ro,t12,t2,t1,xint,xintk
      integer j,k
c functions
      real * 8 ssd_dgauss,ssp_vegas,ssp,ssd
      external ssd_dgauss,ssp_vegas
      real * 8 zero,one
      integer ngau
      data zero,one/0.d0,1.d0/,ngau/0/
      common/tovegas/ro,t1,t2
      real*8 xint2,xint2_err,accum,accerr,eps,dgauss1
      integer imaxpts
      real * 8 xl,xu,acc,chi2a
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      integer ndim,ncall,itmx,nprn
      save xg,wg,ngau
c
c....questi ro,t12,t2 danno i limiti di sp. fasi, saranno da mettere massless
c
      ro   = 4*xm2/s
      t12  = (pp2+xm2)/s
c      t2  = sqrt(exp(2*y)*t12)
      t2  = exp(y)*sqrt(t12)
      t1  = t12/t2
cc      t1 = sqrt(t12)*exp(-y)
c      print*,t1,t2,sqrt(xm2),sqrt(pp2)
      if(t1+t2.gt.1) then
         dfpy = 0
         return
      endif
      
      if ( im0 .eq. 0 .or. im0 .eq. 1) then    ! massive & massless
c      if ( im0 .eq. -1 .or. im0 .eq. -1) then ! this would always do vegas below
       if(ngau.eq.0) then
         if(maxfcn.gt.maxgau) then
            write(*,*)
     #      'maxfcn too large: increase maxgau and recompile'
            stop
         endif
         ngau = maxfcn
         call gauleg(zero,one,xg,wg,ngau)
       endif
       xint = 0
       do 1 j=1,ngau
        xint = xint + ssd(xg(j),ro,t1,t2)*wg(j)
   1   continue
       if(lead.ne.1) then
        do 2 j=1,ngau
           xintk = 0
           do 3 k=1,ngau
               xintk = xintk + ssp(xg(j),xg(k),ro,t1,t2)*wg(k)
   3           continue
           xint = xint + xintk*wg(j)
   2       continue
       endif
      
       dfpy = xint
      else        ! integration with vegas (never used)
                  ! vegas tends to give a NAN with the massive
		  ! contribution, because it probes too extremely
		  ! near the edges      
       eps = 1e-3
       xint = dgauss1(ssd_dgauss,0d0,1d0,eps)
      
       if (lead .ne. 1) then
        ndim = 2
        do k=1,ndim
          xl(k)=0
          xu(k)=1
        enddo
        nprn = 1
        accum=0
        accerr=0
	acc = eps
	ncall = 5000
	itmx = 2
         call vegas(ssp_vegas,xint2,xint2_err,chi2a)
	ncall = 10000
	itmx = 10
         call vegas1(ssp_vegas,xint2,xint2_err,chi2a) ! inizialize cumulants, not grid
       else
         xint2 = 0.d0
	 xint2_err = 0.d0
       endif
      endif
      dfpy = xint + xint2
      end 
      
      real*8 function ssd_dgauss(x)
      implicit none
      real*8 x,ssd,ro,t1,t2
      common/tovegas/ro,t1,t2
      ssd_dgauss = ssd(x,ro,t1,t2)
      end

      real*8 function ssp_vegas(x,wgt)
      implicit none
      real*8 x(2),wgt,ssp,ro,t1,t2
      common/tovegas/ro,t1,t2
      ssp_vegas = ssp(x(1),x(2),ro,t1,t2)
      end




      function dfpxf(pp2,xf,xm2,s,maxfcn)
      implicit none
      real * 8 dfpxf,pp2,xf,xm2,s
      integer maxfcn
c local
      real * 8 ro,b,t12,xfb,t1pt2,y,xjac
c functions
      real * 8 dfpy
      ro   = 4*xm2/s
      b    = sqrt(1-ro)
      t12  = (pp2+xm2)/s
      xfb  = xf*b
      t1pt2 = sqrt(xfb**2+4*t12)
      y     = log( (t1pt2+xfb)/(t1pt2-xfb) )/2
c------------------------------------------------------------------
c The extra factor b/(t1+t2) is to go from d sigma/dy/dpp2 to
c d sigma /dxf/dpp2
c
      xjac = b/t1pt2
      dfpxf = dfpy(pp2,y,xm2,s,maxfcn)*xjac
      return
      end

      function ssp(y1,y2,ro,t1,t2)
c------------------------------------------------------
c When integrated in the unit square in y(1:2), gives
c the contribution to the (inclusive cross section times s)
c for the production of a heavy quark coming from terms
c of the form pp(1/tx) and pl(log(tx)/tx)
c with distributions pp and pl defined so that their
c integral in tx in the unit interval vanishes.
c It calls a function rab(ap,apl,ro,t1,tx), which
c returns in ap the coefficient of p(1/tx), and in apl
c the coefficient of p(log(tx)/tx) in the invariant
c amplitude, and f1(x1),f2(x2) the parton distribution
c of the incoming particles.
c
      implicit none
      real * 8  ssp,y1,y2,ro,t1,t2
c minimum and maximum of Egamma/Eelectron in electroproduction
      real * 8 zmin,zmax
      common/zlim/zmin,zmax
      real * 8 xm2true,rotrue,arotrue
      integer im0
      common/true/xm2true,rotrue,arotrue,im0
c local
      real * 8 t3,tx,at3min,xlg3,at3,atx,xlgr,xjac1,xjac2,
     # r,at2,at1,x1,x2,aro,ap,apl,xlgtx,gf1,gf2
c functions
      real * 8 par
      data par/4.0d0/
c-------------------------------------------------------------------
c At the endpoints return zero.
c
      if(y2.ne.0.and.y2.ne.1.and.y1.ne.0.and.y1.ne.1) go to 1
      ssp=0
      return
c-------------------------------------------------------------------
c  twilda quantities are defined as:
c  ro_twilda  :aro
c  taux_twilda:atx
c  tau1_twilda:at1
c  ...
c
   1  continue
c-------------------------------------------------
c The integration:
c dx1 dx2 = d at1/at1 dat2/at2 *x1*x2 = d at1/at1 dat2/at2 * x1*x2
c d at1 d at2 = d at3 d r *at1^2/at3, (r=at2/at1) thus
c dx1 dx2 = d at3 dr /(r*at3) * x1*x2
      t3=t1+t2
      tx=1-t3
      xlg3=dlog(t3)*y2**par
      at3=dexp(xlg3)
      atx=1-at3
      if(atx.eq.0) then
        ssp=0
        return
      endif
      xlgr=dlog((at3-t1)/t1)*y1+dlog(t2/(at3-t2))*(1-y1)
      xjac1=dlog((at3-t1)*(at3-t2)/t1/t2)/at3
      xjac2=-dlog(t3)*at3*y2**(par-1)*par
      r=dexp(xlgr)
c the extra factor x1*x2 is included in the luminosities.
c the normalization factor 1/(x1*x2)^2 appearing in rab and vab
c accounts for the 1/shat^2 to be supplied to the standard formulae.
c The factor 1/s^2 is supplied in the main program in the final normalization.
c------------------------------------
c Find all values.
c
      at2=at3*r/(1+r)
      at1=at3/(1+r)
      x1=t2/at2
      x2=t1/at1
c....mia modifica, MC
c      if(im0.eq.1) then  
          arotrue=rotrue/x1/x2
c      else
          aro=ro/x1/x2
c      endif
c------------------------------------------------------------
c First call to amplitude values:
c
      call rab(ap,apl,aro,at1,at2,atx,x1,x2)
      ssp=xjac1*(ap+apl*dlog(atx))/atx
c------------------------------------------------------------
c Now redefine all quantities at tx=0.
c
      xlgr=dlog((1-t1)/t1)*y1+dlog(t2/(1-t2))*(1-y1)
      r=dexp(xlgr)
      at2=r/(1+r)
      at1=1/(1+r)
      x1=t2/at2
      x2=t1/at1
c....mia modifica, MC
c      if(im0.eq.1) then  
          arotrue=rotrue/x1/x2
c      else
          aro=ro/x1/x2
c      endif
      xjac1=dlog((1-t1)*(1-t2)/t1/t2)
c------------------------------------------------------------
c Call to amplitude values in the soft limit:
c int( (gf*atx+1)/atx, atx, 0, tx)  gives zero if 1/atx is
c regulated by +.
c same for int( (gf2*atx+1)*log(atx)/atx ,atx,0,tx)
c
      xlgtx=dlog(tx)
      gf1 = -(6*xlgtx-5)/tx/2.0+2*atx*(3*xlgtx-4)/tx**2+(-3.0)*atx**2*(2
     1   *xlgtx-3)/(2.0*tx**3)
      gf2 = -(18*xlgtx**2-30*xlgtx+19)/(tx*(6*xlgtx-11))/2.0+2*atx*(9*xl
     1   gtx**2-24*xlgtx+26)/(tx**2*(6*xlgtx-11))+(-9.0)*atx**2*(2*xlgtx
     2   **2-6*xlgtx+7)/(2.0*tx**3*(6*xlgtx-11))
      call rab(ap,apl,aro,at1,at2,0.d0,x1,x2)
      if(ap.eq.0.and.apl.eq.0)go to 2
      ssp=ssp+xjac1*
     # ( -ap/atx*(gf1*atx+1) - apl*dlog(atx)/atx*(gf2*atx+1) )
c---------------------------------------
c Total
c
  2   continue
      ssp=xjac2*ssp
      return
      end


      function ssd(y1,ro,t1,t2)
c------------------------------------------------------
c When integrated between zero and one in y1 returns
c the contribution to the inclusive cross section from the
c deltoid piece. It calls the value of the deltoid
c term in the invariant amplitude vab(at1,aro) and
c the structure functions f1,f2.
c
      implicit double precision (a-h,o-z)
      common/true/xm2true,rotrue,arotrue,im0
      ssd=0
c------------------------------------------------
c when at1~td watch for 1/v singulatity.
c
      td=dsqrt(1-ro/4/t1/t2)/2
c     td=0.5d0
      if(t1.lt.0.5d0-td.and.t2.lt.0.5d0-td) go to 10
c-------------------------------------------------
c Return zero at the boundaries.
c
      if(y1.eq.0.or.y1.eq.1)return
c-------------------------------------------------
c Go from at1 (tau_tilda_1) to y1.
c
      xlgr=dlog((1-t1)/t1)*y1+dlog(t2/(1-t2))*(1-y1)
      xjac=dlog((1-t1)*(1-t2)/t1/t2)
      r=dexp(xlgr)
c------------------------------------
c Find all values.
c
      at2=r/(1+r)
      at1=1-at2
      x1=t2/at2
      x2=t1/at1
c....mia modifica, MC
c      if(im0.eq.1) then  
          arotrue=rotrue/x1/x2
c      else
          aro=ro/x1/x2
c      endif
c------------------------------------------------
c integrand value:
c
       ssd=xjac*vab(aro,at1,at2,x1,x2)
      return
c------------------------------------------------
c Alternative variables of integration for 1/v
c singularity.
c Return zero at low boundary.
c
 10   continue
      if(y1.eq.0)return
c-------------------------------------------------
c Contribution for at1<.5
c
      z = dexp(dlog(t1/(td+.5d0-t1))*(1-y1)+dlog(1/td/2)*y1)
      at1 = (td+.5d0)*z/(1+z)
      at2 = 1-at1
      xjac = 1/at2/(1+z)*dlog((td+.5d0-t1)/2/td/t1)
      x1 = t2/at2
      x2 = t1/at1
c....mia modifica, MC
c      if(im0.eq.1) then  
          arotrue=rotrue/x1/x2
c      else
          aro=ro/x1/x2
c      endif
      ssd = xjac*vab(aro,at1,at2,x1,x2)
c----------------------------------------------------------
c Add contribution for at1>.5
c
      z = dexp(dlog(t2/(td+.5d0-t2))*(1-y1)+dlog(1/td/2)*y1)
      at2 = (td+.5d0)*z/(1+z)
      at1 = 1-at2
      xjac = 1/at1/(1+z)*dlog((td+.5d0-t2)/2/td/t2)
      x1 = t2/at2
      x2 = t1/at1
c....mia modifica, MC
c      if(im0.eq.1) then  
          arotrue=rotrue/x1/x2
c      else
          aro=ro/x1/x2
c      endif
      ssd = ssd+xjac*vab(aro,at1,at2,x1,x2)
      return
      end

      subroutine rab(ap,apl,ro,at1,at2,atx,x1,x2)
      implicit double precision (a-h,o-z)
      character * 2 beam*6,proc,asy
      common/hvqtrs/xmu2,xmuph2,xcsi,xcsiph,as,nfl,lead,proc,beam,asy
c.....common for massless limit
      common/true/xm2true,rotrue,arotrue,im0
c
      data pi/3.141 592 653 589 793/
      as3o2p  = as**3/(2*pi)
      call xlum(beam,proc,nfl,xmu2,xmuph2,x1,x2,
     #          xlgg,xlqa,xlaq,xlqg,xlgq,xlag,xlga)
      gg   = 0
      ggl  = 0
      qa   = 0
      qal  = 0
      qg   = 0
      qgl  = 0
      gq   = 0
      gql  = 0
      aq   = 0
      aql  = 0
      ag   = 0
      agl  = 0
      ga   = 0
      gal  = 0
c
c....mia modifica, inserimento funzioni nel limite massless
c
      if(im0.eq.1) then
c.....massless limit here
c
      if(asy.ne.'as') then
         if(proc.ne.'qq'.and.proc.ne.'qg') then
            gg  =  hqhpggm0(atx,at1,arotrue) +
     #             hqbpgg(atx,at1,ro)*xcsi
     #             + cthpgg(atx,at1,ro,nfl)
            ggl =  hqhlggm0(atx,at1,arotrue) 
     #             + cthlgg(atx,at1,ro,nfl)
         endif
         if(proc.ne.'gg'.and.proc.ne.'qg') then
            qa  =  hqhpqam0(atx,at1,arotrue) +
     #             hqbpqa(atx,at1,ro)*xcsi
     #             + cthpqa(atx,at1,ro,nfl)
            qal =  hqhlqam0(atx,at1,arotrue) 
     #             + cthlqa(atx,at1,ro,nfl)
         endif
         if(proc.ne.'gg'.and.proc.ne.'qq') then
            qg  =  hqhpqgm0(atx,at1,arotrue) +
     #             hqbpqg(atx,at1,aro)*xcsi
     #             + cthpqg(atx,at1,ro,nfl)
            qgl =  hqhlqgm0(atx,at1,arotrue) 
     #             + cthlqg(atx,at1,ro,nfl)
            gq  =  hqhpqgm0(atx,at2,arotrue) +
     #             hqbpqg(atx,at2,ro)*xcsi
     #             + cthpqg(atx,at2,ro,nfl)
            gql =  hqhlqgm0(atx,at2,arotrue) 
     #             + cthlqg(atx,at2,ro,nfl)
         endif
         aq = qa
         aql = qal
         ag = qg
         agl = qgl
         ga = gq
         gal = gql
      else
         if(proc.ne.'gg'.and.proc.ne.'qg') then
            qa   =  ashpqam0(atx,at1,arotrue)
            aq   =  - qa
         endif
         if(proc.ne.'gg'.and.proc.ne.'qq') then
            qg   =  ashpqgm0(atx,at1,arotrue)
            ag   =  - qg
            gq   =  ashpqgm0(atx,at2,arotrue)
            ga   =  - gq
         endif
      endif
c
      else
c
c....massive calculation here
c
      if(asy.ne.'as') then
         if(proc.ne.'qq'.and.proc.ne.'qg') then
            gg  =  hqhpgg(atx,at1,ro) +
     #             hqbpgg(atx,at1,ro)*xcsi+ cthpgg(atx,at1,ro,nfl)
            ggl =  hqhlgg(atx,at1,ro) + cthlgg(atx,at1,ro,nfl)
         endif
         if(proc.ne.'gg'.and.proc.ne.'qg') then
            qa  =  hqhpqa(atx,at1,ro) +
     #             hqbpqa(atx,at1,ro)*xcsi+ cthpqa(atx,at1,ro,nfl)
            qal =  hqhlqa(atx,at1,ro) + cthlqa(atx,at1,ro,nfl)
         endif
         if(proc.ne.'gg'.and.proc.ne.'qq') then
            qg  =  hqhpqg(atx,at1,ro) +
     #             hqbpqg(atx,at1,ro)*xcsi+ cthpqg(atx,at1,ro,nfl)
            qgl =  hqhlqg(atx,at1,ro) + cthlqg(atx,at1,ro,nfl)
            gq  =  hqhpqg(atx,at2,ro) +
     #             hqbpqg(atx,at2,ro)*xcsi+ cthpqg(atx,at2,ro,nfl)
            gql =  hqhlqg(atx,at2,ro) + cthlqg(atx,at2,ro,nfl)
         endif
         aq = qa
         aql = qal
         ag = qg
         agl = qgl
         ga = gq
         gal = gql
      else
         if(proc.ne.'gg'.and.proc.ne.'qg') then
            qa   =  ashpqa(atx,at1,ro)
            aq   =  - qa
         endif
         if(proc.ne.'gg'.and.proc.ne.'qq') then
            qg   =  ashpqg(atx,at1,ro)
            ag   =  - qg
            gq   =  ashpqg(atx,at2,ro)
            ga   =  - gq
         endif
      endif
c
      endif
c
      ap =   xlgg * gg   + xlqa * qa +
     #       xlaq * aq   + xlqg * qg  + xlgq * gq
     #     + xlag * ag   + xlga * ga
      apl =  xlgg * ggl  + xlqa * qal +
     #       xlaq * aql  + xlqg * qgl + xlgq * gql
     #     + xlag * agl  + xlga * gal
      ap  = as3o2p * ap  /  (x1 * x2)**2
      apl = as3o2p * apl /  (x1 * x2)**2
      return
      end

      function vab(ro,at1,at2,x1,x2)
      implicit double precision (a-h,o-z)
      character * 2 beam*6,proc,asy
c.....common for massless limit
      common/true/xm2true,rotrue,arotrue,im0
c
      common/hvqtrs/xmu2,xmuph2,xcsi,xcsiph,as,nfl,lead,proc,beam,asy
c------------------------------------------------------
c xlramu = log(mu_ren^2/mu_fact^2)
c
      common/rensca/xlramu,ffact,irunsc,istrsc
      integer nlfp1sch
      common/asnf/nlfp1sch
      data pi/3.141 592 653 589 793/
      as3o2p  = as**3/(2*pi)
      tpoas   = 2*pi/as
c---------------------------------------------
c 4 pi b0 * log(mu_ren^2/mu_fact^2)
c
      fpb0    = (33-2*nfl)*xlramu/3.d0
c correction to be included in the gg term to allow for the use
c of as(nf+1), and of F_g^{nf+1}
      if(nlfp1sch.eq.1) then
         cnfgg = 2.d0/3.d0*(-xlramu)
         cnfqa = 2.d0/3.d0*(-xcsi)
      elseif(nlfp1sch.eq.0) then
         cnfgg = 0
         cnfqa = 0
      else
         write(*,*) ' error: flag for nlf+1 scheme should be 1 or 0'
         stop
      endif
c      cnfgg = 0.
c      cnfqa = 0.
      call xlum(beam,proc,nfl,xmu2,xmuph2,x1,x2,
     #          xlgg,xlqa,xlaq,xlqg,xlgq,xlag,xlga)
      gg   = 0
      qa   = 0
      qg   = 0
      gq   = 0
      aq   = 0
      ag   = 0
      ga   = 0
c
c....mia modifica, inserimento funzioni nel limite massless
c
      if(im0.eq.1) then
c.....massless limit here
c
      if(lead.ne.1) then
         if(asy.ne.'as') then
            if(proc.ne.'qg'.and.proc.ne.'qq')
     #         gg  = (   (tpoas+fpb0+cnfgg)*hqh0gg(at1,ro)
     #                  + hqhdggm0(at1,arotrue,nfl) +
     #                    hqbdgg(at1,ro)*xcsi
     #                    + cthdgg(at1,ro,nfl) )
            if(proc.ne.'gg'.and.proc.ne.'qg')
     #         qa  = ( (tpoas+fpb0+cnfqa)*hqh0qa(at1,ro)
     #                + hqhdqam0(at1,arotrue,nfl) +
     #                  hqbdqa(at1,ro,nfl)*xcsi 
     #                  + cthdqa(at1,ro,nfl) )
            if(proc.ne.'gg'.and.proc.ne.'qq') then
               qg  =  cthdqg(at1,ro,nfl)
               gq  =  cthdqg(at2,ro,nfl)
            endif
            aq = qa
            ag = qg
            ga = gq
         else
            if(proc.ne.'gg'.and.proc.ne.'qg') qa = ashdqam0(at1,arotrue)
            aq = -qa
         endif
         vab = xlgg * gg  + xlqa * qa +
     #         xlaq * aq  + xlqg * qg  + xlgq * gq
     #       + xlag * ag  + xlga * ga
      else
         if(proc.ne.'qq'.and.proc.ne.'qg') 
     #           gg  = tpoas*hqh0gg(at1,ro)
         if(proc.ne.'gg'.and.proc.ne.'qg') 
     #           qa  = tpoas*hqh0qa(at1,ro)
         vab = xlgg * gg + (xlqa + xlaq)*qa
      endif
c
      else
c
c....massive calculation here
c

      if(lead.ne.1) then
         if(asy.ne.'as') then
            if(proc.ne.'qg'.and.proc.ne.'qq')
     #         gg  = (   (tpoas+fpb0+cnfgg)*hqh0gg(at1,ro)
     #                  + hqhdgg(at1,ro,nfl) +
     #                    hqbdgg(at1,ro)*xcsi+ cthdgg(at1,ro,nfl) )
            if(proc.ne.'gg'.and.proc.ne.'qg')
     #         qa  = ( (tpoas+fpb0+cnfqa)*hqh0qa(at1,ro)
     #                + hqhdqa(at1,ro,nfl) +
     #                  hqbdqa(at1,ro,nfl)*xcsi + cthdqa(at1,ro,nfl) )
            if(proc.ne.'gg'.and.proc.ne.'qq') then
               qg  =  cthdqg(at1,ro,nfl)
               gq  =  cthdqg(at2,ro,nfl)
            endif
            aq = qa
            ag = qg
            ga = gq
         else
            if(proc.ne.'gg'.and.proc.ne.'qg') qa  =  ashdqa(at1,ro)
            aq = -qa
         endif
         vab = xlgg * gg  + xlqa * qa +
     #         xlaq * aq  + xlqg * qg  + xlgq * gq
     #       + xlag * ag  + xlga * ga
      else
         if(proc.ne.'qq'.and.proc.ne.'qg') gg  = tpoas*hqh0gg(at1,ro)
         if(proc.ne.'gg'.and.proc.ne.'qg') qa  = tpoas*hqh0qa(at1,ro)
         vab = xlgg * gg + (xlqa + xlaq)*qa
      endif
c
      endif
c
      vab = as3o2p * vab /  (x1 * x2)**2
      return
      end

c------------------------------------------------------------------------
c Routines to find the luminosity in a hadronic collision.
c beam = '+pr+pr'  proton-proton
c        '+pr-pr'  proton-antiproton
c        '+pi+pr'  pion-proton
c etc.
c proc = 'qq','qg', 'gg' or anything else to indicate all processes.
c nf   = number of light flavours
c xmu2 = mu^2 (scale of structure functions)
c xlgg = G_1(x1)*G_2(x2)
c etc.  (q stands for quark, a for antiquark)
c
      subroutine xlum(beam,proc,nf,xmu2,xmuph2,x1,x2,
     #                xlgg,xlqa,xlaq,xlqg,xlgq,xlag,xlga)
      implicit double precision (a-h,o-z)
      character * 6 beam, proc*2
      integer nset1,nset2
      common/pdfs/nset1,nset2
c electron or photon are always beam 1
      call fxab2
     # (beam(1:3),nset1,proc,nf,xmuph2,x1,
     #  xfg1,xfu1,xub1,xfd1,xdb1,xfs1,xfc1,xfb1)
      call fxab2
     # (beam(4:6),nset2,proc,nf,xmu2,x2,
     #  xfg2,xfu2,xub2,xfd2,xdb2,xfs2,xfc2,xfb2)
      xlgg = 0
      xlqg  = 0
      xlgq  = 0
      xlag  = 0
      xlga  = 0
      xlqa  = 0
      xlaq  = 0
      if(proc.ne.'qg'.and.proc.ne.'qq') xlgg = xfg1 * xfg2
      if(proc.ne.'qq'.and.proc.ne.'gg') then
         xlqg  = ( xfu1 + xfd1 + xfs1 + xfc1 + xfb1 ) * xfg2
         xlgq  = ( xfu2 + xfd2 + xfs2 + xfc2 + xfb2 ) * xfg1
         xlag  = ( xub1 + xdb1 + xfs1 + xfc1 + xfb1 ) * xfg2
         xlga  = ( xub2 + xdb2 + xfs2 + xfc2 + xfb2 ) * xfg1
      endif
      if(proc.ne.'qg'.and.proc.ne.'gg') then
         tmp   =   xfs1 * xfs2 + xfc1 * xfc2 + xfb1 * xfb2
         xlqa  =  xfu1 * xub2 + xfd1 * xdb2 + tmp
         xlaq  =  xub1 * xfu2 + xdb1 * xfd2 + tmp
      endif
      return
      end

      function dmlin(fun,xl,xh,n,maxfcn,aerr,rerr,ier)
      implicit double precision (a-h,o-z)
      parameter (ngm=50,maxd = 4)
      dimension xl(1),xh(1),x(maxd),ik(maxd),ij(maxd),yy(maxd)
      dimension xg(ngm),wg(ngm)
      data ng/0/,xm/-.5d0/,xp/.5d0/
      idim = n
      ifc  = maxfcn
      if(ifc.gt.ngm) then
        m = (ifc-1)/ngm + 1
        ifc = ngm
      else
        m = 1
      endif
      if(ifc.ne.ng) then
        call gauleg(xm,xp,xg,wg,ifc)
        ng = ifc
      endif
    1 continue
      if( idim .eq. 1 ) then
        y = 0
        do 2 k=1,m
        do 2 j=1,ng
          x(1) = xl(1) + (xh(1)-xl(1))*(k+xg(j)-.5d0)/m
          y = y + (xh(1)-xl(1))/m * fun(1,x) * wg(j)
    2   continue
        go to 4
      endif
      j = 1
      k = 1
      y = 0
    3 continue
      x(idim) = xl(idim) +
     # (xh(idim)-xl(idim))*(k+xg(j)-.5d0)/m
c-------------------------------------------
c Save i,k,y before recursive call
c
      ij(idim) = j
      ik(idim) = k
      yy(idim) = y
      idim = idim - 1
c------------------------------------------ call
      go to 1
c------------------------------------------ return
    4 continue
c--------------------------------- Final exit condition.
      if( idim .ge. n) go to 10
c--------------------------------------
c Restore idim,j,k after recursive call.
c
      idim = idim + 1
      k = ik(idim)
      j = ij(idim)
c------------------------------- restore y + calculated value.
      y = yy(idim) + (xh(idim)-xl(idim))/m * y * wg(j)
c------------------------------- loop back.
      j = j + 1
      if ( j .le. ng ) go to 3
      j = 1
      k = k + 1
      if ( k .le. m ) go to 3
      go to 4
  10  continue
      dmlin = y
      return
      end
