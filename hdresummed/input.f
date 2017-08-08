cccccccccc    input file    ccccccccccccccccccccccc
      subroutine hdrsparam
      implicit none
c      implicit real*8 (a-h,l-z)
      logical verbose
      common/chat/verbose
      integer jmar
      common/valu/jmar
      integer ipt,irap
      real * 8 ptevo,rap,enh1,enh2
      common/evopt/ipt,irap,ptevo(50),rap(50),enh1,enh2
      integer ih1,ih2
      common/hadr/ih1,ih2
      integer isigm
      common/hdinput/isigm
      include 'xsect.h'
      real * 8 hc2,zrs,zal,cm,cmu,cmp
      common/choi/hc2,zrs,zal,cm,cmu,cmp
      real * 8 masch2,masbo2,masto2,lambda2
      common/alfacm/masch2,masbo2,masto2,lambda2
      integer ichan,ichnum
      common/channels/ichan(16),ichnum
      real * 8 xmb,xmc
      common/hvqmass/xmb,xmc
      real * 8 cmu0
      common/cmu0/cmu0
      character*1 hvqs
      common/hvqtype/hvqs
      real * 8 alfa,beta
      integer npsm
      common/alfabeta/alfa,beta,npsm
      integer isubtr
      common/isubtr0/isubtr
      real * 8 lam5qcd
      common/lamqcd/lam5qcd
      integer nfl
      common/flnumb/nfl
      integer imsvknfl,imsvscfl
      common/masskin/imsvknfl,imsvscfl
      integer ioutput
      common/out/ioutput
      integer isuda
      common/sudakov/isuda
      real * 8 h1pdf,h2pdf
      common/mypdfsets/h1pdf(20),h2pdf(20)
      real * 8 q2max,zmin,zmax
      integer iww
      common/wwpar/q2max,zmin,zmax,iww
      real * 8 assc
      common/asscale/assc
      integer ipolemom
      common/ipole/ipolemom
      real * 8 precint
      common/precintc/precint
      character * 2 frscheme
      common/frschemec/frscheme
c function called
      character * 3 pkgname
      integer istrl
c local variables
      integer i,iinput,iptmt
      character * 70 outfile
      real * 8 masch,masbo,masto,tmp
*************************************************
c   zal parameter (not allowed to change)
      zal=1.
*************************************************
      verbose = .true.
*************************************************
      iinput  = 55
c MC Dec 17 2007: set ioutput=15 to output a file hdres.tmp,
c set ioutput=0 to suppress output of file
c      ioutput = 15
      ioutput = 0
c      open(iinput,file='input',status='old')
c.....output filename
      write(*,*) ' enter name of output file'
      read(iinput,*) outfile
      if(ioutput.gt.0)
     #open(ioutput,file=outfile,status='unknown')
      if(ioutput.gt.0)
     #write(ioutput,*) ''''//outfile(1:istrl(outfile))//''''
c Value of the conversion factor: (hbarr*c)**2
      hc2=.389d+9
      write(*,*) ' enter moment of subtraction (integer>=2, default: 4)'
      read(iinput,*) ipolemom
      if(ioutput.gt.0)
     #write(ioutput,10) ipolemom
 10   format(i2,19x,'! pole moment')
C Choice of cross-section
c     isigm=1 ==> dsigma/dy/dpt2
c     isigm=2 ==> e*dsigma/d3p
c     isigm=3 ==> dsigma/dy/dpt
      write(*,*)' enter 1 for dsigma/dy/dpt2,'//
     #     ' 2 for e*dsigma/d3p, 3 for dsigma/dy/dpt'
      read(iinput,*) isigm
      if(ioutput.gt.0)
     #write(ioutput,11) isigm
 11   format(i1,20x,'! isigm')
      write(*,*)
     #  ' enter 1 (1loop) or 2 (2loop) for fragm. func.'
      read(iinput,*) iloopfr
      if(ioutput.gt.0)
     #write(ioutput,22) iloopfr
 22   format(i1,20x,'! iloopfr')
      write(*,*)
     #  ' enter 0 for normal F_b/c(x), 1 for 1 loop approximation,'
     # //'2 for 2 loop approx.'
      read(iinput,*) iwhichsfh
      if(ioutput.gt.0)
     #write(ioutput,23) iwhichsfh
 23   format(i1,20x,'! iwhichsfh')
      write(*,*) ' enter 1 for alfa-1loop, 2 for 2-loops expression'
      write(*,*) ' (2 is normal)'
      read(iinput,*) iloopas
      if(ioutput.gt.0)
     #write(ioutput,24) iloopas
 24   format(i1,20x,'! iloop in alfas')
c
      write(*,*)
     #' enter 1 to include O(as^2) terms in cross section, 0 to exclude'
      write(*,*) ' ( 1 is normal)'
      read(iinput,*) ias2term
      if(ioutput.gt.0)
     #write(ioutput,33) ias2term
 33   format(i1,20x,'! ias2term')
      write(*,*)
     #' enter 1 to include O(as^3) terms in cross section, 0 to exclude'
      write(*,*) ' ( 1 is normal)'
      read(iinput,*) ias3term
      if(ioutput.gt.0)
     #write(ioutput,34) ias3term
 34   format(i1,20x,'! ias3term')
c jmar factorization scheme dependence:
c 1=our factorization scheme cq=1 (c.f. nucl phys. b327,105,
c 2=msbarr factorization scheme (to use with hmrs, mt.....),
c 3=dis factorization scheme (to use with dflm, mt, ....)
cc      jmar=2
      write(*,*)' enter 1 for NUCL PHYS. B327,105 fact. scheme,'
     # //' 2 for MSbar, 3 for DIS ala DFLM'
      read(iinput,*) jmar
      if(ioutput.gt.0)
     #write(ioutput,44) jmar
 44   format(i1,20x,'! jmar')
c  Fragmentation scheme
 45   write(*,*)' enter DL for Delta type fragmentation scheme, '
     # //' MS for MSbar'
      read(iinput,*) frscheme
      if(frscheme.ne.'DL'.and.frscheme.ne.'MS') then
         write(*,*) ' MS or DL!'
         goto 45
      endif
      if(ioutput.gt.0)
     #write(ioutput,'(1x,a4,16x,a)') ''''//frscheme//'''',
     # '! fragmentation scheme'
 46   write(*,*)' for alternative evolution solution (1 on, 0 off)'
      write(*,*)' 0 is normal'
      read(iinput,*) ialtev
      if(ialtev.ne.0.and.ialtev.ne.1) then
         write(*,*) ' 1 or 0!'
         goto 46
      endif
      if(ioutput.gt.0)
     #write(ioutput,'(i1,20x,a)') ialtev,'! alternative evmat'
c   ipt number of points in pt
      ipt=0
      write(*,*)' enter pt values, < 0 to terminate'
 53   read(iinput,*) tmp
      if(tmp.gt.0) then
         ipt=ipt+1
         if(ipt.le.50) then
            ptevo(ipt)=tmp
            goto 53
         else
            write(*,*) ' too many pt points!'
            stop
         endif
      endif
      if(ioutput.gt.0)
     #write(ioutput,'(d16.10,5x,a)') (ptevo(i),'! pt values',i=1,ipt)
      if(ioutput.gt.0)
     #write(ioutput,'(d16.10,5x,a)') -1.d0,'! end pt values'
      write(*,'(a)') ' enter 1 for pt->sqrt(pt^2+m^2)'
      read(iinput,*) iptmt
      if(ioutput.gt.0)
     #write(ioutput,'(i2,19x,a)') iptmt,'! 1 for pt->sqrt(pt^2+m^2)'
c   irap number of points in rapidity
      irap=0
      write(*,*)' enter rapidity values, >=1000 to terminate'
 56   read(iinput,*) tmp
      if(tmp.lt.1000) then
         irap=irap+1
         if(irap.le.50) then
            rap(irap)=tmp
            goto 56
         else
            write(*,*) ' too many y points!'
            stop
         endif
      endif
      if(ioutput.gt.0)
     #write(ioutput,'(d16.10,5x,a)') (rap(i),'! y values',i=1,irap)
      if(ioutput.gt.0)
     #write(ioutput,'(d16.10,5x,a)') 2000.d0,'! end y values'
c The original package treats simmetrically hadrons 1 and 2; we
c decided that electrons and photons can only come from the right, 
c following HERA conventions and the pointlike resummed code
 71   continue
      write(*,*) 
     # 'enter beam 1: -1 pbar,0 (p+n)/2, 1 p, 2 n, 3 pi+, 4 gamma, 5 el'
      write(*,*)
     # ' 100 for photon with only gamma->c component (for testing)),'
      write(*,*) ' its energy in GeV,'
      if(pkgname().eq.'pdl')then
         write(*,*) ' its pdf (Nptype Ngroup Nset (PDFLIB)).'
      elseif(pkgname().eq.'mlm')then
         write(*,*) ' its pdf (0 0 ndns (two zeros followed by MLM #));'
         write(*,*) ' enter ndns<0 to get the list of the PDFS.'
      endif
      read(iinput,*) ih1,enh1,(h1pdf(i),i=1,3)
      if(pkgname().eq.'mlm'.and.h1pdf(3).lt.0) then
         call pdftype
         goto 71
      endif
      write(*,*) ' enter hadron 2 type'
      write(*,*) ' its energy in GeV,'
      if(pkgname().eq.'pdl')then
         write(*,*) ' its pdf (Nptype Ngroup Nset (PDFLIB)).'
      elseif(pkgname().eq.'mlm')then
         write(*,*) ' its pdf (0 0 ndns (two zeros followed by MLM #));'
         write(*,*) ' enter ndns<0 to get the list of the PDFS.'
      endif
      read(iinput,*) ih2,enh2,(h2pdf(i),i=1,3)
      if(ioutput.gt.0)
     #write(ioutput,77)
     # ih1,enh1,(h1pdf(i),i=1,3),ih2,enh2,(h2pdf(i),i=1,3)
 77   format(i2,2x,f9.4,1x,f2.0,1x,f2.0,1x,f6.0,1x,
     #       '! hadron 1, its energy, its PDF set',/,
     #i2,2x,f9.4,1x,f2.0,1x,f2.0,1x,f6.0,1x,
     #'! hadron 2, its energy, its PDF set')
c Check if the PDF set numbers don't contradict the packege we link to
c      if(    (h1pdf(1).ne.0 .and. pkgname().ne.'pdl')
c     #  .or. (h1pdf(1).eq.0 .and. pkgname().ne.'mlm')   ) then
c         write(*,*) ' PDF set does not match the package being used!'
c         stop
c      endif         
c WARNING: the following might give rise to inconsistencies. The file
c conv.f, which performs the convolution of the photon densities with
c the WW function, has ITS OWN VERSION of WW. In principle, one 
c should keep only one version, preferably that in ../common/fww_ww.
c.....whether to have the photon PDF convoluted with a Weizsaecker-Williams
c      if(pkgname().eq.'pdl')then
c         if(ih2.eq.2) then
c            write(*,*)'iww,q2max,zmin,zmax For photon in electron'
c            write(*,*)'iww=0: incoming photon, iww#0 incoming electron'
c            read(iinput,*) iww,q2max,zmin,zmax
c            write(ioutput,79) iww,q2max,zmin,zmax
c         endif
c      else
c the following should prevent undesired use of the code
      iww=-1
      q2max=0.d0
      zmin=0.d0
      zmax=0.d0
c      endif
 79   format(i1,2x,f6.2,2x,f8.4,2x,f8.4,3x,'! iww, q2max, zmin, zmax')
c.....QCD Lambda_5 in MSbar (if zero the h1pdf set value is used)
      write(*,*) ' enter lambda5 (0 for default)'
      read(iinput,*) lam5qcd
      if(ioutput.gt.0)
     #write(ioutput,155) lam5qcd
 155  format(f5.4,16x,'! lambda5 (GeV)')
c choice of distribution functions following pdflib
      call setlambda(lam5qcd)
c....number of flavours in alpha and struct. and fragm. funct.
      write(*,*)' enter number of light flavour in alpha, struct.'
     # // ' func., and fragm. func.'
      read(iinput,*) nfl
      if(ioutput.gt.0)
     #write(ioutput,35) nfl
 35   format(i1,20x,'! number of flavours')
c   zrs square root of s
      zrs= dsqrt((enh1+enh2)**2 - (enh1 - enh2)**2)
c      write(ioutput,88) zrs
c 88   format(f8.2,13x,'! zrs, center of mass energy')
c   cm factorization scale for distribution functions: m=cm*pt
cc      cm=1.
c   cmu renormalization scale(alphas): mu=cmu*pt
cc      cmu=1.
c   cmp factorization scale for fragmentation: mp=pt*zyy
cc      cmp=1.
      write(*,*)' enter cm, cmu, cmp'
      write(*,*)' Fact. scale for distribution functions: m=cm*pt'
      write(*,*)' Renormalization scale(alphas): mu=cmu*pt'
      write(*,*)' Factorization scale for fragmentation: mp=cmp*pt'
      read(iinput,*) cm, cmu, cmp
      if(ioutput.gt.0)
     #write(ioutput,99) cm, cmu, cmp
 99   format(3(f5.2,1x),3x,
     &'! cm, cmu, cmp, rescaling factors for scales')
c
      write(*,*)
     #' enter 0 for normal, 1 to keep pdf scale fixed in struct. func.'
      write(*,*) ' used only for testing purposes'
      read(iinput,*) istrsc
      if(ioutput.gt.0)
     #write(ioutput,100) istrsc
 100  format(i1,20x,'! 0 for mu0*cm, 1 for mu0 in str. funct.')
c mass values for c,b and t quarks (depend on the choice of
c distribution functions):
c for mrs,hmrs,kmrs: masto=1.d+10 (5 flavours), for dflm
c (resp. mt) masto=40 (resp. 90)
cc      masch=1.5d0
cc      masbo=4.75d0
      write(*,*) ' enter mass for charm and bottom'
      read(iinput,*) masch, masbo
	xmc = masch
	xmb = masbo
      if(ioutput.gt.0)
     #write(ioutput,111) masch, masbo
 111  format(2(f4.2,3x),7x,'! charm mass, bottom mass')
      masto=1.d+10
c thresholds for active flavour (depend on the choice
c of distribution functions)
      masch2=masch**2
      masbo2=masbo**2
      masto2=masto**2
c....choice of produced heavy quark (c or b (or g!!)) (Or, rather,
c....selects the singular fragmentation function which needs
c....pole subtraction)
      write(*,*) ' enter b or c for heavy flavour type'
      read(iinput,*) hvqs
      if(ioutput.gt.0)
     #write(ioutput,'(1x,a3,17x,a)') ''''//hvqs//'''',
     # '! flavour with singular behaviour'
c Keep it zero always, P.N. 10/11/2000
       imsvknfl=0
cc....massive (1) or non massive (0) kinematics
c      write(*,*) ' massive (1) or non massive (0) kinematics'
c      read(iinput,*) imsvknfl
c      write(ioutput,37) imsvknfl
c 37   format(i1,20x,'! massive (1) or non massive (0) kinematics')
c Keep it zero always, P.N. 10/11/2000
       imsvscfl=0
cc....massive (1) or non massive (0) ren/fact. scale
c      write(*,*) ' massive (1) or non massive (0) ren/fact scales'
c      read(iinput,*) imsvscfl
c      write(ioutput,39) imsvscfl
c 39   format(i1,20x,'! massive (1) or non massive (0) ren/fact. scale')
c....flag for pole subtraction
      write(*,*) ' 1 for pole subtr. (normal), 0 otherwise'
      read(iinput,*) isubtr
      if(ioutput.gt.0)
     #write(ioutput,166) isubtr
 166  format(i1,20x,'! isubtr for pole subtraction 1=yes, 0=no')
c....scale factor for the scale of the initial state fragm. funct.
c....mu0 = cmu0*quark mass
cc      cmu0 = 1d0
      write(*,*)' scale factor for initial scale of fragm. func.'
      read(iinput,*) cmu0
      if(ioutput.gt.0)
     #write(ioutput,133) cmu0
 133  format(f4.1,17x,'! rescaling factor for fragm. initial scale')
c.....switch for Sudakov form factors (isuda=0: no, isuda=1: yes)
      write(*,*) '0 to exclude, 1 to include sudakov resummation'
      read(iinput,*) isuda
      if(ioutput.gt.0)
     #write(ioutput,173) isuda
 173  format(i2,19x,'! isuda: Sudakov form factors (1) or not (0)')
c.....parameters for the non-perturbative smearing
c******* npsm: 0 -> no non pert. smearing; 1 -> non pert. smear. with
c******* alfa and beta given
cc	npsm = 0
cc	alfa = 1d0
cc	beta = 1000d0
 176  continue
      write(*,*)' parameters for non-perturbative smearing '
      write(*,*)' [0-> no smearing | 1->with smearing], alfa, beta'
      read(iinput,*) npsm, alfa, beta
      if(npsm.eq.0)then
         alfa=0.d0
         beta=0.d0
      elseif(npsm.ne.1)then
         write(*,*)' non-implemented option'
         goto 176
      endif
      if(ioutput.gt.0)
     #write(ioutput,177) npsm, alfa, beta
 177  format(1x,i1,1x,2(f9.3,1x),'! npsm, alfa, beta')
c....selection of active parton channels
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
c....everything
c	ichnum = 16
c	data ichan /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/
c....gluon final state only
c	ichnum = 6
c	data ichan /2,4,7,12,14,15/
c....gg -> Q and qqbar -> Q only
c	ichnum = 2
c	data ichan /5,16/
c	ichnum = 1
c	data ichan /12/
c....gg initial state only
c	ichnum = 2
c	data ichan /15,16/
c....qq initial state only
c	ichnum = 9
c	data ichan /1,2,3,4,5,6,7,11,12/
c....qg initial state only
c	ichnum = 5
c	data ichan /8,9,10,13,14/
      write(*,'(a)') ' Channel numbering:',
     #   ' 16==>  G  G   ---> QJ',
     #   ' 15==>  G  G   ---> G',
     #   ' 14==>  QI G   ---> G',
     #   ' 13==>  QI G   ---> QI',
     #   ' 12==>  QI QBI ---> G',
     #   ' 11==>  QI QBI ---> QI',
     #   ' 10==>  QI G   ---> QBI',
     #   ' 9 ==>  QI G   ---> QBK',
     #   ' 8 ==>  QI G   ---> QK',
     #   ' 7 ==>  QI QI  ---> G',
     #   ' 6 ==>  QI QI  ---> QI',
     #   ' 5 ==>  QI QBI ---> QK',
     #   ' 4 ==>  QI QBK ---> G',
     #   ' 3 ==>  QI QBK ---> QI',
     #   ' 2 ==>  QI QK  ---> G',
     #   ' 1 ==>  QI QK  ---> QI'
      write(*,*) ' enter the number of channels to include,'
      write(*,*) ' followed on the next line by the list of channels'
      read(iinput,*) ichnum
      if(ioutput.gt.0)
     #write(ioutput,'(i2,19x,a)') ichnum,
     # '! # of channel to include'
      read(iinput,*) (ichan(i),i=1,ichnum)
      if(ioutput.gt.0)
     #write(ioutput,144)  (ichan(i), i=1,ichnum)
 144  format(1x,16(i2,1x))
      write(*,*) ' relative precision of the integrals (try 0.001)'
      read(iinput,*) precint
      if(ioutput.gt.0)
     #write(ioutput,145)  precint
 145  format(1x,d10.4,10x,'! precision of the integrals')
      write(*,*)' alfas rescale factor: alfas->alfas times scalef.'
      write(*,*)' only for testing, always choose 1 for normal use'
      read(iinput,*) assc
      if(ioutput.gt.0)
     #write(ioutput,154)  assc
 154  format(f8.3,13x,'! alpha_s scale factor')
      if(ioutput.gt.0)
     #write(ioutput,188)
 188  format(/,'   **** Cross Section (pb) ****',/,
     #'    pt       eta        sigma          error',/)
      if(hvqs.eq.'c') then
         tmp=xmc
      else
         tmp=xmb
      endif
      if(iptmt.eq.1) then
         do i=1,ipt
            ptevo(i)=sqrt(ptevo(i)**2+tmp**2)
         enddo
      endif
      end
