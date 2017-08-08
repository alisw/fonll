c Program to calculate heavy quarks photoproduction.
c From hvqms, by Paolo Nason and M. Cacciari, 15-5-97
c Turned into a photoproduction program by S. Frixione and P. Nason, 20-7-99
c Modified on 29-7-99 by S. Frixione to take into account to boost peculiar
c of HERA physics. A new function (dfp2, dsig/dpt^2 with eta cut) has
c been introduced
c
c
c
      subroutine phms(result)
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
      integer nfl,lead
      real * 8 xmu2,xcsi,xcsiph,as,zeh
      common/hvqphtrs/xmu2,xcsi,xcsiph,as,zeh,nfl,lead,proc,beam,asy
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
c schemes for hadron and photon legs
      character * 2 schhad,schpho
      common/schemes/schhad,schpho
c minimum and maximum of Egamma/Eelectron in electroproduction
      real * 8 zmin,zmax
      common/zlim/zmin,zmax
c local variables
      integer maxfcn, nptype, ngroup, nset,ny,j,npt,ill,k,mufl,icmflag,
     #        lfile
      real * 8 e,s,fren,ptfix,xm,xm2,ro,ymax,xmtrue,pp2,xmu,xmure2,
     # y,tmp1,tmp2,tmp3,tmp4,etamin,etamax,pphoton,phadron,xmtarget,
     # xmuph,ffactph,rscalpha
      logical logfile
c functions
      real * 8 phdfpy,alfas,phdfpxf,phdfp,dfp2,
     #         phdfy,phctpy,phdfxf,phctxdfp
      character * 3 pkgname
c-----------------------------------------------
c  Conversion factor from Gev**(-2) to microbarns.
c
      real * 8 xhc,pi
      integer iun7,iun8,iun9, iun10, iun17,iun18,iun19
      data xhc/389.3857/
      data pi/3.141 592 653 589 793/
      data iun7/7/,iun8/8/,iun9/9/, iun10/10/,
     # iun17/17/,iun18/18/,iun19/19/
c-----------------------------------------------
c On ibm only, inhibits underflow exceptions.
c
c      call xuflow(0)
c--------------------------------------------------
      logfile = .false.

c # of gaussian points
c
      if (logfile) open(unit=66,file='phms-log.tmp',status='unknown')
      write(*,'(1x,a)')'enter maxfcn = # of gaussian points'
      read(55,*) maxfcn
      if (logfile) write(66,'(i3,10x,a)') maxfcn, 
     #                                   '! # of gaussian points'
      if(pkgname().eq.'mlm') then
         write(*,'(1x,a)')'enter pdf set'
         read(55,*) nset
c nptype and ngroup not used by mlmpdf;
         nptype=0
         ngroup=0
         if (logfile) write(66,'((1x,i3),9x,a)') nset,
     #    '! nptype, ngroup, nset'      
      else
         write(*,'(1x,a)')'enter nptype, ngroup, nset'
         read(55,*) nptype, ngroup, nset
c nptype and ngroup not used by pdflib too! make this clearer some day;
         if (logfile) write(66,'(3(1x,i3),1x,a)') nptype, ngroup, nset,
     #    '! nptype, ngroup, nset'      
      endif
      call selectpdf(nptype, ngroup, nset)
c     set lambda5
      call setlam5(xlam,schhad)
c The following lines are needed if the program is linked to PDFLIB,
c which does not return the scheme (and Lambda is possibly wrong)
      if(schhad.eq.'**')then
         write(*,*)' enter Lambda_5'
         read(55,*) xlam
         if (logfile) write(66,'(d10.4,3x,a)') xlam, '! Lambda_5'
         write(*,*)'enter scheme for hadron: ''DI'' or ''MS'''
         read(55,*) schhad
         if (logfile) write(66,'(1x,a4,8x,a)')
     #   ''''//schhad//'''','! scheme for hadron'
      else
c One might want to set Lambda to a value different from that of the pdfs
         write(*,*)' enter Lambda_5 (0 for default)'
         read(55,*) tmp1
         if(tmp1.gt.0.d0)xlam=tmp1
         if (logfile) write(66,'(d10.4,3x,a)') xlam, '! Lambda_5'
      endif
c assign photon scheme
      schpho='**'
      dowhile(schpho.ne.'MS'.and.schpho.ne.'DI')
         write(*,*)'enter now scheme for photon: ''DI'' or ''MS'''
         write(*,*)
     #   '(use DI when hadr. component is GRV-like, MS otherwise)'
         write(*,*)'( in practice DI for GRV and GRS )'
         read(55,*) schpho
      enddo
      if (logfile) write(66,'(1x,a4,8x,a)')
     #   ''''//schpho//'''','! scheme for photon'
      write(*,*)
     # ' enter dfpy,dfpxf,ctpy,dfp,dfp2,dfy,dfxf,ctxdfp',
     # ' for task selection;',
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
      elseif(task.eq.'dfp2') then
         write(*,*)' d sigma/ d(pT**2) for etamin<eta<etamax'
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
     # 'also enter photon momentum, hadron momentum and hadron mass',
     # '(enter m_had<0 for m_had=0.938 GeV). Valid entries are:',
     # 'ph pr+ 200 (photon-proton at E_cm=200 GeV)',
     # 'ph pr+ 20 820 0 (photon-proton with k_gamma=20 GeV',
     # '                 k_pr=820 GeV, and neglecting proton mass)',
     # 'ph nu+ 30 (photon-nucleon ( i.e. (p+n)/2 ) at E_cm=30 GeV',
     # 'ph nu+ 100 0 -1 (100 GeV photon on nucleon at fixed target)',
     # 'el pr+ 300   (electron-proton at E_cm=300 GeV)'
      call getstr(tmpstr)
      if (logfile) write(66,'(a,a)') tmpstr(1:30),
     #                               '! process and CM energy'
      beam = tmpstr(1:6)
      read(tmpstr(7:),fmt=*,err=512,end=512) pphoton,phadron,xmtarget
      if(xmtarget.lt.0) xmtarget=0.938
      e=sqrt( 2*pphoton*sqrt(xmtarget**2+phadron**2)
     #       +xmtarget**2+2*pphoton*phadron )
      goto 513
 512  read(tmpstr(7:),fmt=*) e
c assume it is a proton (none was specified)
      xmtarget=0.938
      pphoton=0.d0
      phadron=0.d0
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
c lab. beam energy:
         tmp1=(s-xmtarget**2)/(2*xmtarget)
         ycm=0.5d0*
     #     log( ((tmp1+xmtarget)+tmp1)/((tmp1+xmtarget)-tmp1) )
      elseif(icmflag.eq.3) then
         if(pphoton.eq.0.d0.and.phadron.eq.0.d0)then
           write(*,*)' do not know how to compute the boost: stop'
           stop
         endif
         tmp1=pphoton+sqrt(phadron**2+xmtarget**2)
         tmp2=pphoton-phadron
         ycm=0.5d0*log( (tmp1+tmp2) / (tmp1-tmp2) )
         write(*,*) '!!!!!!!!!!***** WARNING! *****!!!!!!!!!!!'
         write(*,*)
     #   ' in our convention the electron has positive rapidity'
      else
         write(*,*) ' should have entered 1,2 or 3'
         stop
      endif
c if electron beam, get z limits
      if(beam(1:2).eq.'el') then
c call initalization of Weiszaeker-Williams function
c input from standard input, output to 66
c         call fww_ww0(6,66)
 223     write(*,*)
     # ' enter zmin, zmax (minimum and maximum of Egamma/Eelectron)'
         read(55,*) zmin, zmax
         if (logfile) write(66,'(2(f5.3,1x),1x,a)') zmin,zmax,
     #        '! minimum and maximum of Egamma/Eelectron'
         if(zmin.ge.zmax.or.zmin.lt.0.or.zmax.gt.1) then
            write(*,*) ' must have 0<=zmin<zmax<=1'
            goto 223
         endif
      else
         zmax=1
         zmin=1
      endif
      write(*,'(1x,a,/,a)')
     # 'enter pg, pq, any other for all',
     # '(- for asymmetries, + without asymmetries)'
      call getstr(tmpstr)
      if (logfile) write(66,'(a,a)') tmpstr(1:13),
     #                               '! pg, pq, any other for all'
      if(tmpstr(1:1).eq.'-') then
         proc = tmpstr(2:)
         asy  = 'as'
         write(*,*)'will return (sigma(Q)-sigma(Q_bar))/2'
      elseif(tmpstr(1:1).eq.'+') then
         proc = tmpstr(2:)
         asy  = 'no'
         write(*,*)'will return (sigma(Q)+sigma(Q_bar))/2'
      else
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
      write(*,*)' enter ffact, xmu_fact(hadron)=ffact*mu0'
      read(55,*) ffact
      if (logfile) write(66,'(d10.4,3x,a)') ffact, 
     #   '! xmu_fact(hadron)=ffact*mu0'
      write(*,*) ' ffact =', ffact
      write(*,*) ' enter 1 to keep ffact=1 only in structure functions'
      write(*,*) ' (for debugging purposes only), 0 otherwise'
      read(55,*) istrsc
      if (logfile) write(66,'(i1,12x,a)') istrsc, 
     #     '! 0 normal, 1 for tests'
      write(*,*)' enter ffactph, xmu_fact(photon)=ffactph*mu0'
      read(55,*) ffactph
      if (logfile) write(66,'(d10.4,3x,a)') ffactph,
     #   '! xmu_fact(photon)=ffactph*mu0'
      write(*,*) ' ffactph =', ffactph
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
     # '! mass value, to set rapidity and pt values'
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
        write(*,*)'if task=dfp2, enter two values which will be'
        write(*,*)'interpreted as etamin and etamax in the lab'
        write(*,*)'WARNING: at HERA, out conventions are such that'
        write(*,*)'the photon is coming from the left'
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
         write(*,*)' This option is not implemented: enter the mass,'
         write(*,*)' and then 0,1 or 2'
         goto 211
      endif
      if (logfile) write(66,'(d10.4,1x,i1,1x,a)') xmtrue, im0,
     #'! mass + flag, 0 for massive, 1 for massless, 2 for massless LL'
      if(im0.ne.0) then
         write(*,*) ' enter 1 if you want pt -> sqrt(pt^2+m^2)'
         read(55,*) j
         if (logfile) write(66,'(i2,11x,a)')j,
     #      '!1 if you want pt -> sqrt(pt^2+m^2)'
         if(j.eq.1) then
            do j=1,npt
               xpt(j)=sqrt(xpt(j)**2+xmtrue**2)
            enddo
         endif
      endif
      write(*,*) ' enter 3*charge, -100 for default'
      read(55,*) zeh
      if (logfile) write(66,'(f6.0, 7x,a)') zeh, 
     #     '! 3*charge, -100 for default'
      if(zeh.eq.-100) then
         if(xmtrue.lt.3) then
            write(*,*) ' assuming a charm quark, ch=2/3'
            zeh=2.d0/3.d0
         elseif(xmtrue.lt.7) then
            write(*,*) ' assuming a bottom quark, ch=-1/3'
            zeh=-1.d0/3.d0
         elseif(xmtrue.lt.200) then
            write(*,*) ' assuming a top quark, ch=2/3'
            zeh=2.d0/3.d0
         endif
      else
         zeh=zeh/3
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
        open (7,file = file(1:k),status='UNKNOWN')
c,recl=132)
        open (8,file = file(1:k)//'_lo',status='UNKNOWN')
c,recl=132)
        open (9,file = file(1:k)//'_ra',status='UNKNOWN')
c,recl=132)
        if(vareq(1:1).eq.'y') then
           open (17,file = file(1:k)//'-y',status='UNKNOWN')
c,recl=132)
           open (18,file = file(1:k)//'-y'//'_lo',status='UNKNOWN')
c,recl=132)
           open (19,file = file(1:k)//'-y'//'_ra',status='UNKNOWN')
c,recl=132)
        endif
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

        if(task.eq.'dfp2')then
           if(j.eq.1) then
              if(ny.ne.2) then
                 write(*,*) ' need two pseudorapidity bounds'
                 stop
              endif
              ny=1
           endif
           etamin=xy(1)
           etamax=xy(2)
           lead = 1
           b(j,1) = dfp2(pp2,etamin,etamax,xm2,s,maxfcn)
           lead = 0
           a(j,1) = dfp2(pp2,etamin,etamax,xm2,s,maxfcn)
        endif
c-------------------------------------------------------------------
c    Rapidity/xf loop; for dfp2 it does only the normalization
c
        do 1 k=1,ny
           y = xy(k)-ycm
           if(task.eq.'dfpy') then
              lead = 1
              b(j,k) = phdfpy(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = phdfpy(pp2,y,xm2,s,maxfcn)
              if(ill.eq.1) then
                 tmp1 = xm2true
                 tmp2 = rotrue
                 tmp3 = xcsi
                 tmp4 = xcsiph
                 xcsi = 0
                 xcsiph = 0
                 xm2true = xmu2
                 rotrue = 4*xm2true/s
                 a(j,k)=a(j,k)-phdfpy(pp2,y,xm2,s,maxfcn)+b(j,k)
                 xcsiph = tmp4
                 xcsi = tmp3
                 rotrue = tmp2
                 xm2true = tmp1
              endif
           elseif(task.eq.'dfpxf') then
              lead = 1
              b(j,k) = phdfpxf(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = phdfpxf(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'dfp') then
              lead = 1
              b(j,k) = phdfp(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = phdfp(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'dfy') then
              lead = 1
              b(j,k) = phdfy(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = phdfy(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'ctpy') then
              lead = 1
              b(j,k) = phctpy(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = phctpy(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'dfxf') then
              lead = 1
              b(j,k) = phdfxf(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = phdfxf(pp2,y,xm2,s,maxfcn)
           elseif(task.eq.'ctxdfp') then
              lead = 1
              b(j,k) = phctxdfp(pp2,y,xm2,s,maxfcn)
              lead = 0
              a(j,k) = phctxdfp(pp2,y,xm2,s,maxfcn)
           endif
c------------------------------------------------------
c Normalize; xhc is the plank constant times the speed
c of light in Gev^2*microbarns; we have to divide
c by s**2, because in the functions phssp,phssd this
c factor is omitted, and multiply by pi for d^2 pt -> d (pt^2):
c
           a(j,k) = a(j,k)*pi*xhc/s**2
           b(j,k) = b(j,k)*pi*xhc/s**2
           if(b(j,k).ne.0) then
              c(j,k) = a(j,k)/b(j,k)
           else
              c(j,k) = 1
           endif
 1      continue
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
c unformatted write of all data
      if(lfile.gt.0) then
         open(iun10, form='unformatted',
     #        file=file(1:lfile)//'_raw',status='UNKNOWN')
         write(iun10) npt,ny
         write(iun10) (xpt(k),k=1,npt)
         write(iun10) (xy(k),k=1,ny)
         write(iun10) ((a(j,k),j=1,npt),k=1,ny)
         write(iun10) ((b(j,k),j=1,npt),k=1,ny)
      endif
c end unformatted write
      if (vareq(1:2).eq.'y ') then
        if(index(task,'x').ne.0) vareq(1:3)='xf'
        if(mufl.eq.3) vareq(5:)='xmu'
        write( iun17,'(1x,1h!,a,10(1x,e14.8))') vareq,(xpt(k),k=1,npt)
        write( iun18,'(1x,1h!,a,10(1x,e14.8))') vareq,(xpt(k),k=1,npt)
        write( iun19,'(1x,1h!,a,10(1x,e14.8))') vareq,(xpt(k),k=1,npt)
        do k=1,ny
           if(ipt2topt.ne.1)then
              write( iun17,'(10(1x,e14.8))') xy(k),(a(j,k),j=1,npt)
              write( iun18,'(10(1x,e14.8))') xy(k),(b(j,k),j=1,npt)
           else
              write( iun17,'(10(1x,e14.8))') xy(k),
     #                                     (2*xpt(j)*a(j,k),j=1,npt)
              write( iun18,'(10(1x,e14.8))') xy(k),
     #                                     (2*xpt(j)*b(j,k),j=1,npt)
           endif
           write( iun19,'(10(1x,e14.8))') xy(k),(c(j,k),j=1,npt)
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
     #              (2*xpt(j)*a(j,k),k=1,ny)
               write( iun8,'(10(1x,e14.8))') xpt(j),
     #              (2*xpt(j)*b(j,k),k=1,ny)
            endif
            write( iun9,'(10(1x,e14.8))') xpt(j),(c(j,k),k=1,ny)
         enddo
      endif
      result=a(1,1)
      end

c--------------------------------------------------
c returns d sigma/dy (*s*4*2pi^2)
c with pt cut: pt2>pp2
c
      function phdfy(pp2ct,y,xm2,s,maxfcn)
      implicit none
      real * 8 phdfy,pp2ct,y,xm2,s
      integer maxfcn
      integer maxgau
      real * 8 par
      parameter (maxgau=160,par=.5d0)
      real * 8 xg(maxgau),wg(maxgau)
c local
      real * 8 t2ot1,p2mx,xint,v,pp2,xjac
      integer j
c functions
      real * 8 phdfpy
c
      real * 8 vmn,vmx,zero,one
      integer ngau
      data vmn/0.d0/,vmx/1.d0/,ngau/0/
      data zero,one/0.d0,1.d0/
      t2ot1 = exp(2*y)
      p2mx = s*t2ot1/(1+t2ot1)**2 - xm2
c also if p2mx<0:
      if(pp2ct.ge.p2mx) then
        phdfy = 0
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
         xint = xint+xjac*wg(j)*phdfpy(pp2,y,xm2,s,maxfcn)
 1       continue
      phdfy = xint
      return
      end

c--------------------------------------------------
c returns d sigma/d p2 (*s*4*2pi^2)
c with y cut: y<yct
c
      function phdfp(pp2,yct,xm2,s,maxfcn)
      implicit none
      real * 8 phdfp,pp2,yct,xm2,s
      integer maxfcn
      integer maxgau
      parameter (maxgau=160)
      real * 8 xg(maxgau),wg(maxgau)
c local
      real * 8 t12,ymax,yl,xint,y,xjac
      integer j
c function
      real * 8 phdfpy
      real * 8 zero,one
      integer ngau
      data zero,one/0.d0,1.d0/,ngau/0/
      if(pp2.ge.s/4-xm2) then
        phdfp = 0
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
         xint = xint+xjac*phdfpy(pp2,y,xm2,s,maxfcn)*wg(j)
 1       continue
      phdfp = xint
      return
      end

c--------------------------------------------------
c returns d sigma/d p2 (*s*4*2pi^2)
c with eta cut: etamin<eta<etamax
c
      function dfp2(pp2,etamin,etamax,xm2,s,maxfcn)
      implicit none
      real * 8 dfp2,pp2,etamin,etamax,xm2,s
      integer maxfcn
      integer maxgau
      parameter (maxgau=160)
      real * 8 xg(maxgau),wg(maxgau)
      real * 8 ycm
      common/hera/ycm
c local
      real * 8 t12,ymax0,ymax,ymin,y,xint,xjac
      integer j
c function
      real * 8 phdfpy
      real * 8 zero,one
      integer ngau
      data zero,one/0.d0,1.d0/,ngau/0/
      if(pp2.ge.s/4-xm2) then
        dfp2 = 0
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
      t12  = (pp2+xm2)/s
c find y cuts corresponding to given eta cuts
      ymax0 = dlog( 1/dsqrt(t12)/2+dsqrt(1/t12/4-1) )
      call eta2y(etamin,pp2,xm2,ymin)
      call eta2y(etamax,pp2,xm2,ymax)
c now ymin/ymax are in the lab; go to CM
      ymin=max(ymin-ycm,-ymax0)
      ymax=min(ymax-ycm,ymax0)
c
      xint = 0
      do 1 j=1,ngau
         y    = ymin+(ymax-ymin)*xg(j)
         xjac = ymax-ymin
         xint = xint+xjac*phdfpy(pp2,y,xm2,s,maxfcn)*wg(j)
 1    continue
      dfp2 = xint
      return
      end

      function dfpeta(pp2,eta,xm2,s,maxfcn)
      implicit none
      real * 8 dfpeta,pp2,eta,xm2,s
      integer maxfcn
      real * 8 y
      real * 8 phdfpy
      call eta2y(eta,pp2,xm2,y)
      dfpeta=phdfpy(pp2,y,xm2,s,maxfcn)*pp2/(pp2+xm2)
      end

      subroutine eta2y(eta,p2,xm2,y)
      implicit none
      real * 8 eta,p2,xm2,y
      real * 8 p,p0,pperp,ppar
      pperp=sqrt(p2)
      p=pperp*cosh(eta)
      ppar=pperp*sinh(eta)
      p0=sqrt(p**2+xm2)
      y=0.5d0*log((p0+ppar)/(p0-ppar))
      end

c--------------------------------------------------
c returns sigma (*s*4*2pi^2)
c with pt cut: pt2>pp2ct
c and  y  cut: y<yct
c
      function phctpy(pp2ct,yct,xm2,s,maxfcn)
      implicit none
      real * 8 phctpy,pp2ct,yct,xm2,s
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
      real * 8 phdfpy
      real * 8 vmn,vmx,zero,one
      integer ngau
      data vmn/0.d0/,vmx/1.d0/,ngau/0/
      data zero,one/0.d0,1.d0/
      if(pp2ct.ge.s/4-xm2) then
        phctpy = 0
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
            weightp = weightp*phdfpy(pp2,y,xm2,s,maxfcn)
c            write(*,*)'weightp=',weightp
            xint = xint+weighty*weightp
 1          continue
      phctpy = xint
      return
      end

      function phdfxf(pp2ct,xf,xm2,s,maxfcn)
      implicit none
      real * 8 phdfxf,pp2ct,xf,xm2,s
      integer maxfcn
      integer maxgau
      real * 8 par
      parameter (maxgau=160,par=.5d0)
      real * 8 xg(maxgau),wg(maxgau)
c local
      real * 8 ro,b,xfb,pp2mx,xint,v,pp2,t1pt2,y,xjac
      integer j
c functions
      real * 8 phdfpy
      real * 8 vmn,vmx,zero,one
      integer ngau
      data vmn/0.d0/,vmx/1.d0/,ngau/0/
      data zero,one/0.d0,1.d0/
      ro = 4*xm2/s
      b  = sqrt(1-ro)
      xfb = xf*b
      pp2mx = s*(1+xfb)*(1-xfb)/4-xm2
      if(pp2ct.gt. pp2mx) then
         phdfxf = 0
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
         xint = xint+xjac*wg(j)*phdfpy(pp2,y,xm2,s,maxfcn)
 1       continue
      phdfxf = xint
      return
      end

c--------------------------------------------------
c returns d sigma/d p2 (*s*4*2pi^2)
c with xf cut: xf<xfct
c
      function phctxdfp(pp2,xfct,xm2,s,maxfcn)
      implicit none
      real * 8 phctxdfp,pp2,xfct,xm2,s
      integer maxfcn
      integer maxgau
      parameter (maxgau=160)
      real * 8 xg(maxgau),wg(maxgau)
c local
      real * 8 ro,b,t12,xfbmax,xfbct,xint,xfb,t1pt2,y,xjac
      integer j
c functions
      real * 8 phdfpy
      real * 8 zero,one
      integer ngau
      data zero,one/0.d0,1.d0/,ngau/0/
      if(pp2.ge.s/4-xm2) then
        phctxdfp = 0
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
         phctxdfp = 0
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
         xint = xint+xjac*phdfpy(pp2,y,xm2,s,maxfcn)
 1       continue
      phctxdfp = xint
      return
      end

      function phdfpy(pp2,y,xm2,s,maxfcn)
      implicit none
      real * 8 phdfpy,pp2,y,xm2,s
      integer maxfcn
      integer maxgau
      parameter (maxgau=160)
      real * 8 xg(maxgau),wg(maxgau)
      character * 2 beam*6,proc,asy
      real * 8 zmin,zmax
      common/zlim/zmin,zmax
      real * 8 xmu2,xcsi,xcsiph,as,zeh
      integer nfl,lead
      common/hvqphtrs/xmu2,xcsi,xcsiph,as,zeh,nfl,lead,proc,beam,asy
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
      real * 8 phssd,phssp,ssdmon,sspmon
      real * 8 zero,one
      integer ngau
      data zero,one/0.d0,1.d0/,ngau/0/
      save xg,wg,ngau
c
c....queste scale, e soprattutto xcsi, sono da lasciare massive
c....modifica mia, MC
c
c      if(irunsc.eq.1) then
c         xmu2    = (pp2+xm2true)*ffact**2
c         xcsi    = log(xmu2/xm2true)
c         xmure2 = xmu2*exp(xlramu)
c         as      = alfas(xmure2,xlam,nfl+1)
c      endif
c
c....questi ro,t12,t2 danno i limiti di sp. fasi, saranno da mettere massless
c
      ro   = 4*xm2/s
      t12  = (pp2+xm2)/s
      t2  = sqrt(exp(2*y)*t12)
      t1  = t12/t2
      if(t1+t2/zmax.gt.1) then
         phdfpy = 0
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
      if(beam(1:2).eq.'el') then
         xint = 0
         do j=1,ngau
            xint = xint + phssd(xg(j),ro,t1,t2)*wg(j)
         enddo
         if(lead.ne.1) then
            do j=1,ngau
               xintk = 0
               do k=1,ngau
                  xintk = xintk + phssp(xg(j),xg(k),ro,t1,t2)*wg(k)
               enddo
               xint = xint + xintk*wg(j)
            enddo
         endif
      elseif(beam(1:2).eq.'ph') then
         xint = ssdmon(ro,t1,t2)
         if(lead.ne.1) then
            do j=1,ngau
               xint  =  xint + sspmon(xg(j),ro,t1,t2)*wg(j)
            enddo
         endif
      endif
      phdfpy = xint
      return
      end

      function phdfpxf(pp2,xf,xm2,s,maxfcn)
      implicit none
      real * 8 phdfpxf,pp2,xf,xm2,s
      integer maxfcn
c local
      real * 8 ro,b,t12,xfb,t1pt2,y,xjac
c functions
      real * 8 phdfpy
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
      phdfpxf = phdfpy(pp2,y,xm2,s,maxfcn)*xjac
      return
      end

      function phssp(y1,y2,ro,t1,t2)
c------------------------------------------------------
c When integrated in the unit square in y(1:2), gives
c the contribution to the (inclusive cross section times s)
c for the production of a heavy quark coming from terms
c of the form pp(1/tx) and pl(log(tx)/tx)
c with distributions pp and pl defined so that their
c integral in tx in the unit interval vanishes.
c It calls a function phrab(ap,apl,ro,t1,tx), which
c returns in ap the coefficient of p(1/tx), and in apl
c the coefficient of p(log(tx)/tx) in the invariant
c amplitude, and f1(x1),f2(x2) the parton distribution
c of the incoming particles.
c
      implicit none
      real * 8  phssp,y1,y2,ro,t1,t2
c minimum and maximum of Egamma/Eelectron in electroproduction
      real * 8 zmin,zmax
      common/zlim/zmin,zmax
      real * 8 xm2true,rotrue,arotrue
      integer im0
      common/true/xm2true,rotrue,arotrue,im0
c local
      real * 8 t3,tx,at3min,at3,atx,rmin,rmax,xlgr,xjac1,xjac2,
     # r,at2,at1,x1,x2,aro,ap,apl,xlgtx,gf1,gf2,txmax
c functions
      real * 8 par
      data par/4.0d0/
c-------------------------------------------------------------------
c At the endpoints return zero.
c
      if(y2.ne.0.and.y2.ne.1.and.y1.ne.0.and.y1.ne.1) go to 1
      phssp=0
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
      at3min=t1+t2/zmax
      txmax=1-at3min
      at3=at3min**(y2**par)
      atx=1-at3
      if(atx.lt.1.d-14) then
        phssp=0
        return
      endif
      rmin=t2/(at3*zmax-t2)
      if(zmin.lt.t2/(at3-t1)) then
         rmax=(at3-t1)/t1
      else
         rmax=t2/(at3*zmin-t2)
      endif
      xlgr=dlog(rmax)*y1+dlog(rmin)*(1-y1)
c below the jacobian d r / dy1 /(r*at3)
      xjac1=dlog(rmax/rmin)/at3
c below the jacobian d at3/ d y2
      xjac2=-dlog(at3min)*at3*y2**(par-1)*par
      r=dexp(xlgr)
c the extra factor x1*x2 is included in the luminosities.
c the normalization factor 1/(x1*x2)^2 appearing in phrab and phvab
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
      call phrab(ap,apl,aro,at1,at2,atx,x1,x2)
      phssp=xjac1*(ap+apl*dlog(atx))/atx
c------------------------------------------------------------
c Now redefine all quantities at tx=0.
c
      rmin=t2/(zmax-t2)
      if(zmin.lt.t2/(1-t1)) then
         rmax=(1-t1)/t1
      else
         rmax=t2/(zmin-t2)
      endif
      xlgr=dlog(rmax)*y1+dlog(rmin)*(1-y1)
c below the jacobian d r / dy1 /(r*at3)
      xjac1=dlog(rmax/rmin)
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
c------------------------------------------------------------
c Call to amplitude values in the soft limit:
c int( (gf*atx+1)/atx, atx, 0, tx)  gives zero if 1/atx is
c regulated by +.
c same for int( (gf2*atx+1)*log(atx)/atx ,atx,0,tx)
c
      xlgtx=dlog(txmax)
      gf1 = -(6*xlgtx-5)/txmax/2.0+2*atx*(3*xlgtx-4)/txmax**2+
     #    (-3.0)*atx**2*(2*xlgtx-3)/(2.0*txmax**3)
      gf2 = -(18*xlgtx**2-30*xlgtx+19)/(txmax*(6*xlgtx-11))/2.0+
     #    2*atx*(9*xlgtx**2-24*xlgtx+26)/(txmax**2*(6*xlgtx-11))+
     #    (-9.0)*atx**2*(2*xlgtx**2-6*xlgtx+7)/
     #    (2.0*txmax**3*(6*xlgtx-11))
c
      call phrab(ap,apl,aro,at1,at2,0.d0,x1,x2)
      if(ap.eq.0.and.apl.eq.0)go to 2
      phssp=phssp+xjac1*
     # ( -ap/atx*(gf1*atx+1) - apl*dlog(atx)/atx*(gf2*atx+1) )
c---------------------------------------
c Total
c
  2   continue
      phssp=xjac2*phssp
      return
      end

      function sspmon(y2,ro,t1,t2)
c------------------------------------------------------
c When integrated in the unit segment in y1, gives
c the contribution to the (inclusive cross section times s^2)
c for the production of a heavy quark coming from terms
c of the form pp(1/tx) and pl(log(tx)/tx)
c with distributions pp and pl defined so that their
c integral in tx in the unit interval vanishes.
c It calls a function phrab(ap,apl,ro,t1,tx), which
c returns in ap the coefficient of p(1/tx), and in apl
c the coefficient of p(log(tx)/tx) in the invariant
c amplitude, and f1(x1),f2(x2) the parton distribution
c of the incoming particles.
c
      implicit none
      real * 8 sspmon,y2,ro,t1,t2
      real * 8 xm2true,rotrue,arotrue
      integer im0
      common/true/xm2true,rotrue,arotrue,im0
c local
      real * 8 t3,tx,at3,atx,at2,at1,xjac1,xjac2,x1,x2,aro,
     # ap,apl,xlgtx,gf1,gf2
c functions
      real * 8 par
      data par/4.0d0/
c-------------------------------------------------------------------
c At the endpoints return zero.
c
      if(y2.ne.0.and.y2.ne.1) go to 1
      sspmon=0
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
c dx2 = d at1/at1 *x2 = d at3/at1 *x2
c with at3=at1+t2 (at2=t2)
      t3=t1+t2
      tx=1-t3
      at3=t3**(y2**par)
      atx=1-at3
      if(atx.lt.1.d-14) then
        sspmon=0
        return
      endif
      at2=t2
      at1=at3-at2
      xjac1=1/at1
c below the jacobian d at3/ d y2
      xjac2=-dlog(t3)*at3*y2**(par-1)*par
      x1=1
      x2=t1/at1
c the extra factor x2 is included in the luminosisites.
c the normalization factor 1/(x1*x2)^2 appearing in phrab and phvab
c accounts for the 1/shat^2 to be supplied to the standard formulae.
c The factor 1/s^2 is supplied in the main program in the final normalization.
c....mia modifica, MC
c      if(im0.eq.1) then  
          arotrue=rotrue/x1/x2
c      else
          aro=ro/x1/x2
c      endif
c------------------------------------------------------------
c First call to amplitude values:
c
      call phrab(ap,apl,aro,at1,at2,atx,x1,x2)
      sspmon=xjac1*(ap+apl*dlog(atx))/atx
c------------------------------------------------------------
c Now redefine all quantities at tx=0.
c
      at2=t2
      at1=1-at2
      x1=1
      x2=t1/at1
c....mia modifica, MC
c      if(im0.eq.1) then  
          arotrue=rotrue/x1/x2
c      else
          aro=ro/x1/x2
c      endif
      xjac1=1/at1
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
      call phrab(ap,apl,aro,at1,at2,0.d0,x1,x2)
      if(ap.eq.0.and.apl.eq.0)go to 2
      sspmon=sspmon+xjac1*
     # ( -ap/atx*(gf1*atx+1) - apl*dlog(atx)/atx*(gf2*atx+1) )
c---------------------------------------
c Total
c
  2   continue
      sspmon=xjac2*sspmon
      return
      end

      function phssd(y1,ro,t1,t2)
c------------------------------------------------------
c When integrated between zero and one in y1 returns
c the contribution to the inclusive cross section from the
c deltoid piece. It calls the value of the deltoid
c term in the invariant amplitude phvab(at1,aro) and
c the structure functions f1,f2.
c
      implicit none
      real * 8 phssd,y1,ro,t1,t2
c minimum and maximum of Egamma/Eelectron in electroproduction
      real * 8 zmin,zmax
      common/zlim/zmin,zmax
      real * 8 xm2true,rotrue,arotrue
      integer im0
      common/true/xm2true,rotrue,arotrue,im0
c local
      real * 8 td,rmin,rmax,xlgr,xjac,r,at2,at1,x1,x2,aro,z
c functions
      real * 8 phvab
      phssd=0
c------------------------------------------------
c when at1~td watch for 1/v singulatity.
c
      td=dsqrt(1-ro/4/t1/t2)/2
c     td=0.5d0
c      if(t1.lt.0.5d0-td.and.t2.lt.0.5d0-td) go to 10
c-------------------------------------------------
c Return zero at the boundaries.
c
      if(y1.eq.0.or.y1.eq.1)return
c-------------------------------------------------
c Go from at1 (tau_tilda_1) to y1.
c
      rmin=t2/(zmax-t2)
      if(zmin.lt.t2/(1-t1)) then
         rmax=(1-t1)/t1
      else
         rmax=t2/(zmin-t2)
      endif
      xlgr=dlog(rmax)*y1+dlog(rmin)*(1-y1)
      xjac=dlog(rmax/rmin)
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
       phssd=xjac*phvab(aro,at1,at2,x1,x2)
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
      phssd = xjac*phvab(aro,at1,at2,x1,x2)
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
      phssd = phssd+xjac*phvab(aro,at1,at2,x1,x2)
      return
      end

      function ssdmon(ro,t1,t2)
c------------------------------------------------------
c Returns the contribution to the inclusive cross section from the
c deltoid piece. It calls the value of the deltoid
c term in the invariant amplitude phvab(at1,aro) and
c the structure functions f1,f2.
c
      implicit none
      real * 8 ssdmon,ro,t1,t2
      real * 8 xm2true,rotrue,arotrue
      integer im0
      common/true/xm2true,rotrue,arotrue,im0
c local
      real * 8 at2,at1,xjac,x1,x2,aro
c functions
      real * 8 phvab
      ssdmon=0
      at2=t2
      at1=1-at2
      xjac=1/at1
      x1=1
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
      ssdmon=xjac*phvab(aro,at1,at2,x1,x2)
      end

      subroutine phrab(ap,apl,ro,at1,at2,atx,x1,x2)
      implicit none
      real * 8 ap,apl,ro,at1,at2,atx,x1,x2
      real * 8 alfaem
      parameter (alfaem=1.d0/137)
      real * 8 pi
      parameter ( pi=3.141 592 653 589 793d0)
      real * 8 xmu2,xcsi,xcsiph,as,zeh
      integer nfl,lead
      character * 2 beam*6,proc,asy
      common/hvqphtrs/xmu2,xcsi,xcsiph,as,zeh,nfl,lead,proc,beam,asy
c.....common for massless limit
      real * 8 xm2true,rotrue,arotrue
      integer im0
      common/true/xm2true,rotrue,arotrue,im0
      character * 2 schhad,schpho
      common/schemes/schhad,schpho
c local
      real * 8 as2aemo2p,xlpg,xlpn,xlpc,xlpa,pg,pgl,pn,pnl,pc,pcl,
     # pa,pal
c functions
      real * 8 hqhppg,hqhppgm0,hqbppg,cthppg,hqhlpg,hqhlpgm0,
     # cthlpg,hqhppc,hqhppcm0,hqbppc,cthppc,hqhlpc,hqhlpcm0,
     # cthlpc,hqhppn,hqhppnm0,hqbppn,cthppn,hqhlpn,hqhlpnm0,
     # cthlpn,hqhppa,hqhppam0,hqhlpam0
c
      as2aemo2p  = as**2*alfaem/(2*pi)
      call xlumpho
     # (beam,proc,nfl,xmu2,x1,x2,xlpg,xlpn,xlpc,xlpa)
      pg   = 0
      pgl  = 0
      pn   = 0
      pnl  = 0
      pc   = 0
      pcl  = 0
      pa   = 0
      pal  = 0
c
c....mia modifica, inserimento funzioni nel limite massless
c
      if(im0.eq.1) then
c.....massless limit here
c
         if(asy.ne.'as') then
            if(proc.ne.'pq') then
               pg  =  hqhppgm0(atx,at1,arotrue) +
     #              hqbppg(atx,at1,ro)*xcsi
               pgl =  hqhlpgm0(atx,at1,arotrue) 
               if(schhad.eq.'DI')then
                  pg  = pg + cthppg(atx,at1,ro,nfl)
                  pgl = pgl + cthlpg(atx,at1,ro,nfl)
               endif
            endif
            if(proc.ne.'pg') then
               pc  =  hqhppcm0(atx,at1,arotrue) +
     #              hqbppc(atx,at1,ro)*xcsiph
               pcl =  hqhlpcm0(atx,at1,arotrue) 
               pn  =  hqhppnm0(atx,at1,arotrue) +
     #              hqbppn(atx,at1,ro)*xcsi
               pnl =  hqhlpnm0(atx,at1,arotrue) 
               if(asy.eq.'no')then
                  pa = 0.d0
                  pal = 0.d0
               else
                  pa  =  hqhppam0(atx,at1,arotrue)
                  pal =  hqhlpam0(atx,at1,arotrue)
               endif
               if(schhad.eq.'DI')then
                  pn  = pn + cthppn(atx,at1,ro,nfl)
                  pnl = pnl + cthlpn(atx,at1,ro,nfl)
               endif
               if(schpho.eq.'DI')then
                  pc  = pc + cthppc(atx,at1,ro,nfl)
                  pcl = pcl + cthlpc(atx,at1,ro,nfl)
               endif
            endif
         else
            pa  =  hqhppam0(atx,at1,arotrue)
            pal =  hqhlpam0(atx,at1,arotrue)
         endif
c     
      else
c     
c.... massive calculation here
c 
         if(asy.ne.'as') then
            if(proc.ne.'pq') then
               pg  =  hqhppg(atx,at1,ro) +
     #              hqbppg(atx,at1,ro)*xcsi
               pgl =  hqhlpg(atx,at1,ro) 
               if(schhad.eq.'DI')then
                  pg  = pg + cthppg(atx,at1,ro,nfl)
                  pgl = pgl + cthlpg(atx,at1,ro,nfl)
               endif
            endif
            if(proc.ne.'pg') then
               pc  =  hqhppc(atx,at1,ro) +
     #              hqbppc(atx,at1,ro)*xcsiph
               pcl =  hqhlpc(atx,at1,ro) 
               pn  =  hqhppn(atx,at1,ro) +
     #              hqbppn(atx,at1,ro)*xcsi
               pnl =  hqhlpn(atx,at1,ro) 
               if(asy.eq.'no')then
                  pa = 0.d0
                  pal = 0.d0
               else
                  pa  =  hqhppa(atx,at1,ro)
                  pal =  0.d0
               endif
               if(schhad.eq.'DI')then
                  pn  = pn + cthppn(atx,at1,ro,nfl)
                  pnl = pnl + cthlpn(atx,at1,ro,nfl)
               endif
               if(schpho.eq.'DI')then
                  pc  = pc + cthppc(atx,at1,ro,nfl)
                  pcl = pcl + cthlpc(atx,at1,ro,nfl)
               endif
            endif
         else
            pa  =  hqhppa(atx,at1,ro)
            pal =  0
         endif
      endif
c
      ap =   (xlpg * pg + xlpn * pn) *zeh**2 +
     #        xlpc * pc + xlpa * pa * zeh
      apl =   (xlpg * pgl + xlpn * pnl) *zeh**2 +
     #        xlpc * pcl + xlpa * pal * zeh
c the normalization factor 1/(x1*x2)^2 accounts for the 1/shat^2
c to be supplied to the standard formulae.
c The remaining factor 1/s^2 should be supplied in the main program.
      ap  = as2aemo2p * ap  /  (x1 * x2)**2
      apl = as2aemo2p * apl /  (x1 * x2)**2
      end

      function phvab(ro,at1,at2,x1,x2)
      implicit none
      real * 8 phvab,ro,at1,at2,x1,x2
c.....common for massless limit
      real * 8 xm2true,rotrue,arotrue
      integer im0
      common/true/xm2true,rotrue,arotrue,im0
      character * 2 schhad,schpho
      common/schemes/schhad,schpho
c
      real * 8 xmu2,xcsi,xcsiph,as,zeh
      integer nfl,lead
      character * 2 beam*6,proc,asy
      common/hvqphtrs/xmu2,xcsi,xcsiph,as,zeh,nfl,lead,proc,beam,asy
c------------------------------------------------------
c xlramu = log(mu_ren^2/mu_fact^2)
c
      real * 8 xlramu,ffact
      integer irunsc,istrsc
      common/rensca/xlramu,ffact,irunsc,istrsc
      real * 8 alfaem
      parameter (alfaem=1.d0/137)
      real * 8 pi
      parameter (pi=3.141 592 653 589 793d0)
      integer nlfp1sch
      common/asnf/nlfp1sch
c local
      real * 8 tpoas,fpb0,cnfpg,as2aemo2p,xlpg,xlpn,xlpc,xlpa,
     # pg,pn,pc
c functions
      real * 8 hqh0pg,hqhdpg,hqhdpgm0,hqbdpg,cthdpg,cthdpn,cthdpc
c
      tpoas   = 2*pi/as
c---------------------------------------------
c 4 pi b0 * log(mu_ren^2/mu_fact^2)
c
      fpb0    = (33-2*nfl)*xlramu/3.d0
c correction to be included in the gamma-g term to allow for the use
c of as(nf+1), and of F_g^{nf+1}. It is 1/2 of hadroproduction case
      if(nlfp1sch.eq.1) then
         cnfpg = 1.d0/3.d0*(-xlramu)
      elseif(nlfp1sch.eq.0) then
         cnfpg = 0
      else
         write(*,*) ' error: flag for nlf+1 scheme should be 1 or 0'
         stop
      endif
      as2aemo2p  = as**2*alfaem/(2*pi)
      call xlumpho
     # (beam,proc,nfl,xmu2,x1,x2,xlpg,xlpn,xlpc,xlpa)
      pg   = 0
      pn   = 0
      pc   = 0
c
c....mia modifica, inserimento funzioni nel limite massless
c
      if(im0.eq.1) then
c.....massless limit here
c
         if(lead.ne.1) then
            if(asy.ne.'as') then
               if(proc.ne.'pq') then
                  pg  = (tpoas+fpb0/2+cnfpg)*hqh0pg(at1,ro)
     #                 + hqhdpgm0(at1,arotrue,nfl) +
     #                 hqbdpg(at1,ro)*xcsi
                  if(schhad.eq.'DI')pg=pg+cthdpg(at1,ro,nfl) 
               endif
               if(proc.ne.'pg') then
                  pn = 0
                  pc = 0
                  if(schhad.eq.'DI')pn=pn+cthdpn(at1,ro,nfl)
                  if(schpho.eq.'DI')pc=pc+cthdpc(at1,ro,nfl)
               endif
            endif
            phvab = (xlpg * pg  + xlpn * pn)*zeh**2 + xlpc*pc
         else
            if(proc.ne.'pq') 
     #           pg  = tpoas*hqh0pg(at1,ro)
            phvab = xlpg * pg *zeh**2
         endif
c
      else
c
c....massive calculation here
c
         if(lead.ne.1) then
            if(asy.ne.'as') then
               if(proc.ne.'pq') then
                  pg  = (tpoas+fpb0/2+cnfpg)*hqh0pg(at1,ro)
     #                 + hqhdpg(at1,ro,nfl) +
     #                 hqbdpg(at1,ro)*xcsi
                  if(schhad.eq.'DI')pg=pg+cthdpg(at1,ro,nfl) 
               endif
               if(proc.ne.'pg') then
                  pn = 0
                  pc = 0
                  if(schhad.eq.'DI')pn=pn+cthdpn(at1,ro,nfl)
                  if(schpho.eq.'DI')pc=pc+cthdpc(at1,ro,nfl)
               endif
            endif
            phvab = (xlpg * pg  + xlpn * pn)*zeh**2 + xlpc*pc
         else
            if(proc.ne.'pq') 
     #           pg  = tpoas*hqh0pg(at1,ro)
            phvab = xlpg * pg *zeh**2
         endif
      endif
c
c the normalization factor 1/(x1*x2)^2 accounts for the 1/shat^2
c to be supplied to the standard formulae.
c The remaining factor 1/s^2 should be supplied in the main program.
      phvab = as2aemo2p * phvab /  (x1 * x2)**2
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
      subroutine xlumpho
     # (beam,proc,nf,xmu2,x1,x2,xlpg,xlpn,xlpc,xlpa)
c returns the luminosities times x1x2
      implicit none
      character * 6 beam, proc*2
      integer nf
      real * 8 xmu2,x1,x2,xlpg,xlpn,xlpc,xlpa
      real * 8 chu,chd,fww_ww
      parameter (chu=2.d0/3.d0)
      parameter (chd=-1.d0/3.d0)
c local
      real * 8 xfp,xfg2,xfu2,xub2,xfd2,xdb2,xfs2,xfc2,xfb2
      if(beam(1:2).eq.'ph') then
         xfp=1
      elseif(beam(1:2).eq.'el') then
         xfp = x1 * fww_ww(x1)
      endif
      call fxab
     # (beam(4:6),proc,nf,xmu2,x2,
     #  xfg2,xfu2,xub2,xfd2,xdb2,xfs2,xfc2,xfb2)
      xlpg  = 0
      xlpn  = 0
      xlpc  = 0
      xlpa  = 0
      if(proc.ne.'pq') xlpg = xfp * xfg2
      if(proc.ne.'pg') then
         xlpn  = (xfu2+xub2+xfd2+xdb2+2*(xfs2+xfc2+xfb2))*xfp
         xlpc  = ( (xfu2+xub2+2*xfc2)*chu**2
     #        +(xfd2+xdb2+2*(xfs2+xfb2))*chd**2 ) * xfp
         xlpa  = ((xfu2-xub2)*chu+(xfd2-xdb2)*chd) * xfp
      endif
      end
