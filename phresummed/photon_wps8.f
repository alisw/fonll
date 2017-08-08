c This file was originally named photon.f; now implements pole subtraction
c in the heavy quark fragmentation function (wps=WithPoleSubtraction)
      subroutine phrs(result,error)
c      implicit real*8(a-h,o-z)
      implicit none
      real * 8 result,error
      integer maxy,maxpt
      parameter (maxy=100,maxpt=100)
      REAL*8 LAM,IS,NC,NF,LLD,LLQ,LLG,HO,LLR,LQ
      REAL*8 LAMP
c      real*8 csivec(5)
      real*8 ptvec(maxpt),ylabvec(maxy)
      character*20 parm(20)
      character*70 outfile
      real * 8 alam5
      common/lambda5/alam5
      real * 8 rusc,usc,pusc
      COMMON/CUSC/rusc,USC,PUSC
      real * 8 pi,cas
      COMMON/CST/LAM,PI,CAS
      COMMON/CLA/LAMP
      real * 8 FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,DIF,PIP,PTL,CF
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF 
      real * 8 pka
      COMMON/KAON/PKA
      real * 8 fgg
      COMMON/FSM/FGG
      real * 8 CS1,CS2,CS3
      COMMON/CS/CS1,CS2,CS3
      real * 8 dgg
      COMMON/GLUN/DGG
      real * 8 FLCA
      integer ipi
      COMMON/CHARM/FLCA,IPI
      EXTERNAL SCOD,SLD,SCQ,SPDI,SPIA,SPIB,SCOG,SFU,SLLQ,SLLG,SCG
      EXTERNAL SFQ,SFG,SPDF,SLR,ALPS,FQ1,SLQ
c      EXTERNAL FOPTT,FOPTA
      EXTERNAL FOPTT
      real * 8 RAE,RBE,RCE,RTE
      COMMON/FOPT/RAE,RBE,RCE,RTE
      real * 8 par,capr
      COMMON/GPRO/PAR(30),CAPR(8,20,32)
      real * 8 par1,par2,capol,cavdm
      COMMON/GPHO/PAR1(30),PAR2(30),CAPOL(8,20,32),CAVDM(8,20,32)
      REAL*8 AZ,BZ,EPS,ETA,FOPTT,SOR
      real * 8 SCM,SSCM,PTM,YCM
      COMMON/NEWP/SCM,SSCM,PTM,YCM
      real * 8 PLLA
      COMMON/FLAG/HO,PLLA
      real * 8 qm,xmq
      integer ims
      common/quarkmass/qm,xmq,ims
c
      real*8 alfa,beta,xmb,xmc,cmu0,alam5qcd,qcdl4,qcdl5
      real*8 h1pdf(20),h2pdf(20)
      real*8 assc
      integer npsm,ifbm,ifixed
      logical verbose
      COMMON/W50512/QCDL4,QCDL5
      common/hvqmass/XMB,XMC
      integer nfl
      common/flnumb/nfl
      common/alfabeta/alfa,beta,npsm
      integer iloopas,iloopfr
      common/fonfrloop/iloopas,iloopfr
      common/cmu0/cmu0
      COMMON/lamqcd/alam5qcd
      common/mypdfsets/h1pdf,h2pdf
      integer ih
      COMMON/HADR/IH
      integer isuda
      common/sudakov/isuda
      common/chat/verbose
      character*2 scheme
      common/schemeflag/scheme
      character*1 hvqs
      common/hvqtype/hvqs
      common/fakepdf/ifbm
      common/asscale/assc
      common/fixedorder/ifixed
      real * 8 zmin1,zmax
      common/wwrange/zmin1,zmax
      real*8 res,err,cqcd,rcsi,csi,enh1,enhp,ee,ep,dy,zmin,ylab
      real * 8 tcm,ucm,zmi
      integer iinput,ioutput,iftout,massive,it,i,ipt,imtflag,irap
      integer iww,iplla,ncall0,itmp,nitn0,j
c function calls
      integer istrl
      real * 8 alps
      verbose = .true.
      pi = 4.*atan(1.d0)
      iinput  = 55
c MC Oct 1 2010: set ioutput=15 to output a file phres.tmp,
c set ioutput=0 to suppress output of file
c      ioutput = 15
      ioutput = 0
c      open(iinput,file='input',status='old')
      write(*,*)'enter name of outfile'
      read(iinput,*) outfile
      iftout = ioutput
      if(ioutput.gt.0)
     #   open(ioutput,file=outfile,status='unknown')
      if(ioutput.gt.0)
     #   write(ioutput,*) ''''//outfile(1:istrl(outfile))//''''
      write(*,*)
     # ' Parameters for non-perturbative fragmentation function'
      write(*,*)' npsm: 0 to set NPFF==1'
      write(*,*)'       1 for A(1-z)^alfa*z^beta, '//
     # 'int_0^1 A(1-z)^alfa*z^beta=1'
c      write(*,*)'       2 for anorm*peterson, epsilon=alfa'
      write(*,*)'enter npsm, alfa, beta'
      read(iinput,*) npsm, alfa, beta
      if(npsm.lt.0.or.npsm.gt.1)then
         write(*,*)' non-implemented option'
         stop
      endif
      if(npsm.eq.0)then
         alfa=0.d0
         beta=0.d0
      endif
      if(ioutput.gt.0)
     #write(ioutput,177) npsm, alfa, beta
 177  format(i1,1x,2(f9.3,1x),3x,'! npsm, alfa, beta')
      write(*,*)'enter charm and bottom masses (1.5 and 5 are default)'
      write(*,*)'WARNING: this program only works for charm'
      read(iinput,*) xmc,xmb
      if(ioutput.gt.0)
     #write(ioutput,111) xmc, xmb
 111  format(1x,2(f4.2,3x),6x,'! charm mass, bottom mass')
c....scale factor for the scale of the initial state fragm. funct.
c....mu0 = cmu0*quark mass
      write(*,*)'enter scale factor for the initial scale of FF, '//
     #          '1 is default'
      read(iinput,*) cmu0
      if(ioutput.gt.0)
     #write(ioutput,133) cmu0
 133  format(f4.1,17x,'! rescaling factor for fragm. initial scale')
      write(*,*) ' enter heavy flavour (c or b)'
      read(iinput,*) hvqs
      if(ioutput.gt.0)
     #write(ioutput,122) hvqs
 122  format(a1,20x,'! heavy flavour (c or b)')
      if(hvqs.eq.'b') then
        xmq = xmb
      else
        xmq = xmc
      endif
c Always use massless kinematics.
c Uncomment below if you want to experiment ...
      massive=0
c      write(*,*)'enter 1 for massive, 0 for massless kinematics'
c      read(iinput,*) massive
c      write(ioutput,123) massive
c 123  format(i3,20x,'! massive (1) or not (0) kinematics')
      if(massive.eq.1) then
        qm = xmq
        print*,' Massive kinematics'
      else
        qm = 0.
      endif 
c Always massless scales
      ims=0
c      write(*,*)'enter 1 for massive, 0 for massless scales'
c      read(iinput,*) ims
c      write(ioutput,125) ims
c 125  format(i3,20x,'! massive (1) or not (0) ren/fact. scales')

c.....switch for Sudakov form factors (isuda=0: no, isuda=1: yes)
      write(*,*)'enter 1 for Sudakov form factors, 0 otherwise'
      read(iinput,*) isuda
      if(ioutput.gt.0)
     #write(ioutput,173) isuda
 173  format(i1,20x,'! isuda: Sudakov form factors (1) or not (0)')
C NUMBER OF LOOPS IN ALPHAS AND IN THE EVOLUTION EQUATIONS (DEPENDS
C ON THE CHOICE OF THE DISTRIBUTION FUNCTIONS)
      write(*,*)'enter number of loops in alpha_S (1 or 2)'
      read(iinput,*) iloopas
      if(ioutput.gt.0)
     #write(ioutput,22) iloopas
 22   format(i1,20x,'! iloopas')

C ILOOPFR=1 LOWEST ORDER ONLY, ILOOPFR=2 HIGHER ORDER
      write(*,*)'enter 1 for LO, 2 for NLO'
      read(iinput,*) iloopfr
      if(ioutput.gt.0)
     #write(ioutput,33) iloopfr
 33   format(i1,20x,'! iloopfr')
      HO= iloopfr-1
c      iloopas = ho + 1
c      iloopfr = ho + 1

C************************************************************************
C     FLCA=1.   avec charme
C     IPI=2     set 2 des fonctions de fragmentation de JPHGUILLET
C***********************************************************************
      FLCA=1.
      IPI=2

c      IFTOUT = 42
c      OPEN( IFTOUT,FILE='SORTIE',STATUS='OLD' )
cc      PRINT *,LAM,LAMP,alam5
c      CALL WATE8
c      CALL WATE32
c      OPEN(UNIT=11,FILE='INPRO',STATUS='OLD') 
c      OPEN(UNIT=12,FILE='INPL',STATUS='OLD') 
c      OPEN(UNIT=13,FILE='INVDM',STATUS='OLD') 
c      READ(11,*) PAR
c      READ(12,*) PAR1
c      READ(13,*) PAR2
c      READ(11,222) CAPR
c      READ(12,222) CAPOL
c      READ(13,222) CAVDM
c 222  FORMAT(8E15.4)
C************************************************************************
C Toujours mettre FGG=FGQ=DGG=DQG=0.0 POUR SUPPRIMER LES TERMES RAJOUTES
C AUX ORDRES SUPERIEURS POUR DES DEFINITIONS NON UNIVERSELLES DE G/P ET H/G
C************************************************************************
      FQQ=0.
      FQG=0.
      FGG=0.
      FGQ=0.
C************************************************************************
C Ne jamais changer la valeur de FQP. Apres les modifications AEM -> MSbar,
C le programme n'est pas coherent avec FQP=1.
C************************************************************************
      FQP=0.
      DQQ=0.
      DGG=0.
      DQG=0.
      DGQ=0.
      CS2=1.
      PUSC=2.
      PKA=0.
      DIF=0.
      PIP=0.
      CF=4./3.
      NC=3.
      NF=4.
      PI=3.141592653589793d0
      CQCD=1.54
      IT=123456789
      CAS=12.*PI/25.

c hardwired:
      nfl=4

      write(*,*)'renorm. scale = pt*r, factorization scale = pt*f'
      write(*,*)'enter r and f, 1 and 1 are default'
      read(iinput,*) rcsi, csi
      if(ioutput.gt.0)
     #write(ioutput,99) rcsi, csi
 99   format(f5.2,1x,f5.2,10x,
     &'! rcsi, csi, rescaling factors for scales')

cc      do 55 icsi = 1,1
cc      csi = csivec(icsi)
cc      csi = 1.
cc      WRITE( IFTOUT, * ) ' '
cc      WRITE( IFTOUT, 299 ) csi
cc 299  FORMAT(' csi = ',F10.5)

      rusc = rcsi**2
      USC  = csi**2
c Stuff below commented away; ntype decides hadron type...
c ih is only used
c.....hadron type (1 - proton, 2 - photon)
c      write(*,*)'enter 1 for proton, 2 for photon (hadronic)'
c      read(iinput,*) ih
c      write(ioutput,553) ih
c 553   format(i3,18x,'! ih, proton (1) or photon (2)')
c
      write(*,*)' ntype, ngroup, nset are pdflib parameters'
     # //' 3 for pion, 4 for photon, 5 for electron'
      write(*,*) 
     # 'ih=hadron type=1 for proton, 2 for neutron, 0 for (n+p)/2,'
     # //' 3 for pion, 4 for photon, 5 for electron'
      write(*,*) 'enter ih, energy, ntype, ngroup, nset for hadron,'
     # //' and energy for photon/electron'
      read(iinput,*) ih,enh1,(h1pdf(i),i=1,3),enhp
      if(ioutput.gt.0)
     #write(ioutput,'(i1,1x,d11.5,1x,3f5.0,1x,d11.5,a)')
     #  ih,enh1,(h1pdf(i), i=1,3),enhp,
     # ' ! E_H,ntype,ngroup,nset,E_gamma/el'
c Photon pdf not used here
c      read(iinput,*) enhp,(h2pdf(i), i=1,3)
c      write(ioutput,77) enh1,(h1pdf(i), i=1,3),enhp,(h2pdf(i), i=1,3)
c 77   format(f6.2,1x,3f4.0,3x,
c     #       '! hadron energy, its PDF set',/,
c     #f6.2,1x,3f4.0,3x,'! photon energy, its PDF set')


      EE = enhp
      EP = enh1



c.....QCD Lambda_5 in MSbar (if zero the h1pdf set value is used)
      write(*,*)'enter Lambda_5 (0 for the default)'
      read(iinput,*) alam5

C CHOICE OF DISTRIBUTION FUNCTIONS FOLLOWING PDFLIB
      parm(1) = 'Nptype'
      parm(2) = 'Ngroup'
      parm(3) = 'Nset'
      CALL FONLLPDFSET(PARM,h1pdf)         

      if(alam5.eq.0) then
c...mistake here, corrected on 4/12/1996
cc         lam = QCDL4**2
         lam = QCDL4
         alam5 = qcdl5
         if(ioutput.gt.0)
     #    write(ioutput,155) alam5
      else
c final result seems (tested numerically...) to be independent of lam; 
c compensation occurs between the explict expression given in sigmaphoton0 
c and those hidden in the rest of the code
        lam = alam5*(5./alam5)**(2./25.)*
     #             (2*log(5./alam5))**(963./14375.)
        if(ioutput.gt.0)
     #   write(ioutput,156) alam5
      endif

      alam5qcd = alam5
      lamp = lam
         
 155  format(f5.4,16x,'! lambda5 (GeV) got from the set of hadron 1')
 156  format(f5.4,16x,'! lambda5 (GeV) given in input')


c      SCM=.94**2+4.*EP*EE
c Treat the proton as massless
      SCM=4.*EP*EE
      SSCM=SQRT(SCM)
      DY=LOG(EP/EE)/2.
c      WRITE( IFTOUT, 303 ) SSCM
c 303  FORMAT(' SQRT(S) = ',F10.5)
	
C   IPT NUMBER OF POINTS IN PT
cc      IPT=10
      write(*,*)'enter the number of pt points <=',maxpt
      read(iinput,*) ipt
      if(ipt.eq.0)stop
      if(ipt.gt.maxpt) then
         write(*,*) ipt,'>',maxpt,'!'
         stop
      endif
      write(*,*)'enter the pt values'
      read(iinput,*) (ptvec(i),i=1,ipt)
      if(ioutput.gt.0)
     #write(ioutput,551) ipt
 551  format(i3,18x,'! ipt, number of pt points')
      do i=1,ipt
         if(ioutput.gt.0)
     #    write(ioutput,*) ptvec(i)
      enddo
      write(*,*)
     #  ' enter 1 if you want to convert pt->sqrt(pt^2+m^2)'//
     #  ' (anything else for no)'
      read(iinput,*) imtflag
      if(ioutput.gt.0)
     #write(ioutput,'(i1,20x,a)') imtflag,'! 1 for pt->sqrt(pt^2+m^2)'
      if(imtflag.eq.1) then
         do i=1,ipt
            ptvec(i)=sqrt(ptvec(i)**2+xmq**2)
         enddo
      endif
C   irap NUMBER OF POINTS IN rapidity
cc      irap=1
      write(*,*)'enter the number of rapidity points <=',maxy
      read(iinput,*) irap
      if(irap.eq.0)stop
      if(irap.gt.maxy) then
         write(*,*) irap,'>',maxy,'!'
         stop
      endif
      write(*,*)'enter the rapidity values (lab frame)'
      read(iinput,*) (ylabvec(i),i=1,irap)
      if(ioutput.gt.0)
     #write(ioutput,66) irap
 66   format(i3,18x,'! irap, number of rapidity points')
      do i=1,irap
         if(ioutput.gt.0)
     #    write(ioutput,*) ylabvec(i)
      enddo      
c....flag for WW convolution, qmax^2, ymin, ymax (here called zmin and zmax)
      write(*,*)'enter 0 for a monochromatic photon,'
      write(*,*)'1 for Weizsaecker-Williams photon'
      read(iinput,*) iww
      if(ioutput.gt.0)
     #write(ioutput,'(i1,20x,a)' ) iww,
     # '! 0 for a monochromatic photon, 1 for WW'
      if(iww.eq.1) then
c         call fww_ww0(iinput,ioutput)
         write(*,*)'enter zmin,zmax'
         read(iinput,*) zmin,zmax
         if(ioutput.gt.0)
     #    write(ioutput,664) zmin,zmax
 664     format(2(1x,f7.3),5x,'! zmin, zmax in WW') 
         if(zmin.ge.zmax.or.zmin.lt.0.or.zmax.gt.1) then
            write(*,*) ' must have 0<=zmin<zmax<=1'
            stop
         endif
      endif
      iplla=0
c Do not use: what it does is not known
c      write(*,*)'enter 1 to add the resolved part, 0 otherwise'
c      read(iinput,*) iplla
c      write(ioutput,661) iplla
c 661  format(i3,18x,'! plla, resolved part added (1) or not (0)')
      plla = iplla
c....choice of the fragmentation function factorization scheme
      write(*,*)
     # 'enter the factorization scheme for fragmentation function'//
     # ' (MS, PO, EX)'
      read(iinput,*) scheme
      if(ioutput.gt.0)
     #write(ioutput,'(1x,a4,16x,a)') ''''//scheme//'''',
     # '! FF factorization scheme'
c....real hvq pdf (0), or fake one (phfbm) (1)
      write(*,*)
     # 'enter 1 to use fake pdfs, 0 otherwise (0 for real work)'
      read(iinput,*) ifbm
      if(ioutput.gt.0)
     #write(ioutput,159) ifbm
 159  format(i1,20x,'! real or fake hvq pdf')
      write(*,*)'enter alpha scale factor (1 for real work)'
      read(iinput,*) assc
      if(ioutput.gt.0)
     #write(ioutput,154)  assc
 154  format(f8.3,13x,'! alpha_s scale factor')
      write(*,*)'enter 0 to use resummed, 1 for fixed order pdf and FF'
      write(*,*)'(0 for real work)'
      read(iinput,*) ifixed
      if(ioutput.gt.0)
     #write(ioutput,161)  ifixed
 161  format(i1,20x,
     #  '! using resummed (0) or fixed order (1) pdf and ff')
c
      if(iww.eq.0)then
         ncall0=3000
      elseif(iww.eq.1)then
         ncall0=20000
      else
         write(*,*)'wrong iww value: ',iww
         stop
      endif
      write(*,*)'enter number of Vegas calls, -1 for default'
      read(iinput,*) itmp
      if(itmp.ne.-1)ncall0=itmp
      if(ioutput.gt.0)
     #write(ioutput,162) ncall0
 162  format(i5,16x,'! number of Vegas calls')
c
      write(*,*)'enter number of Vegas iterations'
      read(iinput,*) nitn0
      if(ioutput.gt.0)
     #write(ioutput,163) nitn0
 163  format(i3,18x,'! number of Vegas iterations')
c
c Initialization of the routine that saves pdf values
      call cacheinit()
c
      if ( iftout .gt. 0 ) then
         WRITE( IFTOUT, * ) ' '
         WRITE( IFTOUT, * ) '  pt           y_lab             x-sect 
     #(ds/dy/dpt2, pb)'
         WRITE( IFTOUT, * ) ' '
      endif
      
      DO 2 I=1,ipt
      ptm = ptvec(i)

        DO 3 J=1,irap
        ylab = ylabvec(j)
        YCM=-YLAB+DY
        
C*******************************************************************
C     YLAB CORRESPOND A LA CONVENTION DES EXPERIMENTATEURS (<O POUR LE PHOTON)
C     YCM CORRESPOND A MA CONVENTION (<O POUR LE PROTON)
C************************************************************************

      PRINT*,' ptm, ylab = ',PTM,YLAB
      TCM=-PTM*SSCM*EXP(-YCM)
      UCM=-PTM*SSCM*EXP(YCM)
c this neglects m_p^2 in the denominator, but keeps it in the numerator
c      ZMI=(.9-UCM)/(SCM+TCM)      
c this is consistent: neglects m_p^2 everywhere
      zmi=-ucm/(scm+tcm)
      zmin1 = MAX(zmin,zmi)
      print*,' zmi,zmin,zmin1 = ',zmi,zmin,zmin1
c
      if(iww.eq.1) then
c.....call for electroproduction
         if(zmin1.ge.zmax)then
            res=0.d0
            err=0.d0
         else
            call sigmaphotonwb(res,err,ncall0,nitn0)
         endif
      else
c.....call for photoproduction
         call sigmaphotonm(res,err,ncall0,nitn0)
      endif
      if(ioutput.gt.0)
     #write(ioutput,213) ptm,ylab,res*pi,err*pi
 213  FORMAT(2(3X,f8.4),3X,D12.6,' +-',d8.2)

      PRINT*,' res = ',RES,' err=',err
c      WRITE (IFTOUT,307)YLAB,RES,PLLA
 307  FORMAT( ' YLAB ',f12.6,/,' RES ',e12.6,/,' LEADING LOG ',f5.1,/)
      write(*,*)' alpha_s = ',alps(usc*ptm**2)
c      WRITE(IFTOUT,*)YLAB,RES,PLLA
  1   CONTINUE
  3   CONTINUE
  2   CONTINUE
 55   continue
        if ( iftout .gt. 0 ) CLOSE( IFTOUT )
        CLOSE( iinput )
   7  FORMAT ( )
      result=res*pi
      error=err*pi
      END




      FUNCTION ALPS(QE2)
      implicit real*8(a-h,o-z)
      common/lambda5/alam5
      real*8 alfas_p
c      COMMON/CST/LAM,PI,CAS
c      COMMON/CLA/LAMP
c      REAL*8 LAM,LAM2,LAMP
c      LAM2=LAMP**2
c      ALPSA=LOG(QE2/LAM2)
c      ALPS=CAS/ALPSA
c      ALPS=ALPS*(1.-ALPS*LOG(ALPSA)*154./(100.*PI))


      alps = alfas_p(qe2,alam5,4)



      RETURN
      END
      FUNCTION CA1(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      COMMON/GLUN/DGG
c      EXTERNAL FQ1
      EXTERNAL FQ1
      REAL*8 LSQ,LQC,LV,LV1,NF,NC
      COMMON/EC/FQA,FQB,FQC,FQE
      REAL*8 NF4
      NF4=4.
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q1=FQ1(1.d0,X22)
      Q2=FQB
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LQC=LOG(SH/Q1)
      LV=LOG(V)
      LV1=LOG(1.-V)
C
      CFA1=
     1 +CF*(7.-7.*V+7.*V**2-7.*V**(-1))
     1 +CF*LSQ*(-3.+3.*V-3.*V**2+3.*V**(-1))
     1 +CF*LV*(-7.+V+3.*V**2+3.*V**(-1))
     1 +CF*LV*LSQ*(-4.+4.*V-4.*V**2+4.*V**(-1))
     1 +CF*LV*LV1*(12.-16.*V+8.*V**2-4.*V**(-1))
     1 +CF*LV*LV1*FQQ*(-4.+4.*V-4.*V**2+4.*V**(-1))
     1 +CF*LV*FQQ*(-3.+3.*V-3.*V**2+3.*V**(-1))
     1 +CF*LV**2*(-8.+10.*V-6.*V**2+4.*V**(-1))
     1 +CF*LV**2*FQQ*(2.-2.*V+2.*V**2-2.*V**(-1))
     1 +CF*LV1*(8.-8.*V)
     1 +CF*LV1*LSQ*(4.-4.*V+4.*V**2-4.*V**(-1))
     1 +CF*LV1*FQQ*(3.-3.*V+3.*V**2-3.*V**(-1))
     1 +CF*LV1**2*(-12.+12.*V-4.*V**2+4.*V**(-1))
     1 +CF*LV1**2*FQQ*(2.-2.*V+2.*V**2-2.*V**(-1))
     1 +CF*FQQ*(-9.+9.*V-9.*V**2+9.*V**(-1))
C
       CFAP1=
     1 +CF*PI2*(-16./3.+22./3.*V-10./3.*V**2+4./3.*V**(-1))
     1 +CF*PI2*FQQ*(-2./3.+2./3.*V-2./3.*V**2+2./3.*V**(-1))
C
C
      CNA1=
     1 +LSQ*(2./3.*NF-2./3.*V*NF+2./3.*V**2*NF-2./3.*V**(-1)*NF)
     1 +LQC*(-2./3.*NF4+2./3.*V*NF4-2./3.*V**2*NF4+2./3.*V**(-1)*NF4)
     1 +NC*LSQ*(-11./3.+11./3.*V-11./3.*V**2+11./3.*V**(-1))
     1 +NC*LQC*(11./3.-11./3.*V+11./3.*V**2-11./3.*V**(-1))
     1 +NC*LV*(2.-2.*V)
     1 +NC*LV*LSQ*(-4.+4.*V-4.*V**2+4.*V**(-1))
     1 +NC*LV*LV1*(2.*V+2.*V**2-4.*V**(-1))
     1 +NC*LV**2*(-6.+5.*V-7.*V**2+8.*V**(-1))
     1 +NC*LV1*(-4.+4.*V)
     1 +NC*LV1**2*(5.-5.*V+V**2-V**(-1))
C
       CNAP1=
     1 +NC*PI2*(2.-3.*V+V**2)
      CA1=CFA1+CFAP1+CNA1+CNAP1
      CODEL=61./16.-PI**2/3.
      CADEL=DGG*NC*4.*(LV**2/2.-CODEL)*(V**2-V+1.-1./V)
      CA1=CA1+CADEL
      RETURN
      END
      FUNCTION CA2(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      COMMON/GLUN/DGG
c      EXTERNAL FQ1
      EXTERNAL FQ1
      REAL*8 LSQ,LV,LV1,NF,NC
      COMMON/EC/FQA,FQB,FQC,FQE
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q2=FQB
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LV=LOG(V)
      LV1=LOG(1.-V)
       CFA2=
     1 +CF*(3.-3.*V+3.*V**2-3.*V**(-1))
     1 +CF*LSQ*(-4.+4.*V-4.*V**2+4.*V**(-1))
     1 +CF*LV*(-4.+4.*V-4.*V**2+4.*V**(-1))
     1 +CF*LV*FQQ*(4.-4.*V+4.*V**2-4.*V**(-1))
     1 +CF*LV1*FQQ*(-4.+4.*V-4.*V**2+4.*V**(-1))
     1 +CF*FQQ*(-3.+3.*V-3.*V**2+3.*V**(-1))
C
       CNA2=
     1 +NC*LSQ*(-4.+4.*V-4.*V**2+4.*V**(-1))
     1 +NC*LV*(-12.+12.*V-12.*V**2+12.*V**(-1))
     1 +NC*LV1*(4.-4.*V+4.*V**2-4.*V**(-1))
      CA2=CFA2+CNA2
      CAW=DGG*NC*(4.*LV*V**2-4.*LV*V-4.*LV*V**(-1)+4.*LV)
      CA2=CA2+CAW
      RETURN
      END
      FUNCTION CA3(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      COMMON/GLUN/DGG
      REAL*8 NC
      PI2=PI**2
      W=WC
      V=VC
       CFA3=
     1 +CF*(-4.+4.*V-4.*V**2+4.*V**(-1))
     1 +CF*FQQ*(4.-4.*V+4.*V**2-4.*V**(-1))
C
       CNA3=
     1 +NC*(-8.+8.*V-8.*V**2+8.*V**(-1))
      CA3=CFA3+CNA3
      CALW=DGG*NC*(4.*V**2-4.*V-4.*V**(-1)+4.)
      CA3=CA3+CALW
      RETURN
      END
      FUNCTION CA4(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      COMMON/GLUN/DGG
c      EXTERNAL FQ1
      EXTERNAL FQ1
      REAL*8 LSQ,LQC,NF,NC,LV,LV1,LX1,LX2,LW,LW1,LR
      COMMON/EC/FQA,FQB,FQC,FQE
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q1=FQ1(1.d0,X22)
      Q2=FQB
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LQC=LOG(SH/Q1)
      LV=LOG(V)
      LV1=LOG(1.-V)
      V2=2.-V
      V1=1.-V
      X1=1.-V*W
      X2=1.-V+V*W
      W1=1.-W
      LX1=LOG(X1)
      LX2=LOG(X2)
      LW=LOG(W)
      LW1=LOG(1.-W)
      LR=LOG(X1/V1)
       CFA5=
     1 +CF*LV*(6.-6.*V+4.*V*X1**(-1)+12.*V*X2**(-3)-20.*V*X2**(-2)+18.*
     1 V*X2**(-1)-12.*V*W+8.*V*W**2+6.*V**2-12.*V**2*X2**(-3)+16.*V**2*
     1 X2**(-2)-10.*V**2*X2**(-1)+4.*V**2*W**2+4.*V**3*X2**(-3)-4.*V**3
     1 *X2**(-2)+2.*V**3*X2**(-1)+2.*V**(-1)-4.*V**(-1)*W+4.*V**(-1)*W*
     1 *2-4.*X1**(-1)-4.*X2**(-3)+8.*X2**(-2)-10.*X2**(-1)+8.*W-8.*W**2
     1 )
     1 +CF*LV*FQQ*(-4.+6.*V-4.*V*X1**(-1)+2.*V*W-4.*V**2-2.*V**2*W-2.*V
     1 **2*W**2+4.*X1**(-1))
     1 +CF*LV*DGQ*(-6.+4.*V-12.*V*X2**(-3)+20.*V*X2**(-2)-18.*V*X2**(-1
     1 )+2.*V*W-2.*V**2+12.*V**2*X2**(-3)-16.*V**2*X2**(-2)+10.*V**2*X2
     1 **(-1)+2.*V**2*W-2.*V**2*W**2-4.*V**3*X2**(-3)+4.*V**3*X2**(-2)-
     1 2.*V**3*X2**(-1)+4.*X2**(-3)-8.*X2**(-2)+10.*X2**(-1))
C
       CFA6=
     1 +CF*LX1*(-12.+20.*V+4.*V*W-12.*V**2-4.*V**2*W-4.*V**2*W**2)
C
       CFA7=
     1 +CF*LX2*(8.*V-12.*V*W1**(-1)+24.*V*X2**(-3)-40.*V*X2**(-2)+36.*V
     1 *X2**(-1)-4.*V*W-8.*V**2+4.*V**2*W1**(-1)-24.*V**2*X2**(-3)+32.*
     1 V**2*X2**(-2)-20.*V**2*X2**(-1)+4.*V**2*W+8.*V**3*X2**(-3)-8.*V*
     1 *3*X2**(-2)+4.*V**3*X2**(-1)-8.*V**(-1)*W1**(-1)+16.*W1**(-1)-8.
     1 *X2**(-3)+16.*X2**(-2)-20.*X2**(-1))
     1 +CF*LX2*DGQ*(-12.+8.*V-24.*V*X2**(-3)+40.*V*X2**(-2)-36.*V*X2**(
     1 -1)+4.*V*W-4.*V**2+24.*V**2*X2**(-3)-32.*V**2*X2**(-2)+20.*V**2*
     1 X2**(-1)+4.*V**2*W-4.*V**2*W**2-8.*V**3*X2**(-3)+8.*V**3*X2**(-2
     1 )-4.*V**3*X2**(-1)+8.*X2**(-3)-16.*X2**(-2)+20.*X2**(-1))
C
       CFA8=
     1 +CF*LV1*(16.-28.*V-4.*V*W+20.*V**2-4.*V**2*W+8.*V**2*W**2)
     1 +CF*LV1*FQQ*(4.-6.*V+4.*V*X1**(-1)-2.*V*W+4.*V**2+2.*V**2*W+2.*V
     1 **2*W**2-4.*X1**(-1))
C
       CFA9=
     1 +CF*LW*(4.-12.*V+16.*V*W1**(-1)-4.*V*W+8.*V**2-8.*V**2*W1**(-1)+
     1 4.*V**2*W+4.*V**2*W**2+4.*V**(-1)*W1**(-1)-12.*W1**(-1))
     1 +CF*LW*FQP*(-4.+4.*V-8.*V*W+8.*V*W**2+2.*V**(-1)-4.*V**(-1)*W+4.
     1 *V**(-1)*W**2+8.*W-8.*W**2)
C
       CFA10=
     1 +CF*LW1*(10.-14.*V+4.*V*X1**(-1)+12.*V*X2**(-3)-20.*V*X2**(-2)+1
     1 8.*V*X2**(-1)-12.*V*W+8.*V*W**2+14.*V**2-12.*V**2*X2**(-3)+16.*V
     1 **2*X2**(-2)-10.*V**2*X2**(-1)-8.*V**2*W+8.*V**2*W**2+4.*V**3*X2
     1 **(-3)-4.*V**3*X2**(-2)+2.*V**3*X2**(-1)+2.*V**(-1)-4.*V**(-1)*W
     1 +4.*V**(-1)*W**2-4.*X1**(-1)-4.*X2**(-3)+8.*X2**(-2)-10.*X2**(-1
     1 )+8.*W-8.*W**2)
     1 +CF*LW1*FQP*(4.-4.*V+8.*V*W-8.*V*W**2-2.*V**(-1)+4.*V**(-1)*W-4.
     1 *V**(-1)*W**2-8.*W+8.*W**2)
     1 +CF*LW1*FQQ*(-4.+6.*V-4.*V*X1**(-1)+2.*V*W-4.*V**2-2.*V**2*W-2.*
     1 V**2*W**2+4.*X1**(-1))
     1 +CF*LW1*DGQ*(-6.+4.*V-12.*V*X2**(-3)+20.*V*X2**(-2)-18.*V*X2**(-
     1 1)+2.*V*W-2.*V**2+12.*V**2*X2**(-3)-16.*V**2*X2**(-2)+10.*V**2*X
     1 2**(-1)+2.*V**2*W-2.*V**2*W**2-4.*V**3*X2**(-3)+4.*V**3*X2**(-2)
     1 -2.*V**3*X2**(-1)+4.*X2**(-3)-8.*X2**(-2)+10.*X2**(-1))
C
       CFA10=CFA10
     1 +CF*(4.-4.*V+8.*V*W-8.*V*W**2-2.*V**(-1)+4.*V**(-1)*W-4.
     1 *V**(-1)*W**2-8.*W+8.*W**2)
       CFA11=
     1 +CF*(-9.+4.*V-4.*V*X1**(-1)-28.*V*W1**(-1)*LR-36.*V*X2**(-3)+54.
     1 *V*X2**(-2)-8.*V*X2**(-1)+12.*V*W-14.*V*W**2+3.*V**2+12.*V**2*W1
     1 **(-1)*LR+36.*V**2*X2**(-3)-42.*V**2*X2**(-2)-6.*V**2*X2**(-1)-4
     1 .*V**2*W+6.*V**2*W**2-12.*V**3*X2**(-3)+10.*V**3*X2**(-2)+6.*V**
     1 3*X2**(-1)+6.*V**(-1)-12.*V**(-1)*W1**(-1)*LR+10.*V**(-1)*W-12.*
     1 V**(-1)*W**2+4.*X1**(-1)+28.*W1**(-1)*LR+12.*X2**(-3)-22.*X2**(-
     1 2)+8.*X2**(-1)-18.*W+24.*W**2)
     1 +CF*LSQ*(6.-6.*V+4.*V*X1**(-1)+12.*V*X2**(-3)-20.*V*X2**(-2)+18.
     1 *V*X2**(-1)-12.*V*W+8.*V*W**2+6.*V**2-12.*V**2*X2**(-3)+16.*V**2
     1 *X2**(-2)-10.*V**2*X2**(-1)+4.*V**2*W**2+4.*V**3*X2**(-3)-4.*V**
     1 3*X2**(-2)+2.*V**3*X2**(-1)+2.*V**(-1)-4.*V**(-1)*W+4.*V**(-1)*W
     1 **2-4.*X1**(-1)-4.*X2**(-3)+8.*X2**(-2)-10.*X2**(-1)+8.*W-8.*W**
     1 2)
       AUCFA11=
     1 +CF*FQP*(-24.*V*W+24.*V*W**2-12.*V**(-1)*W+12.*V**(-1)*W**2+24.*
     1 W-24.*W**2)
     1 +CF*FQQ*(1.-7.*V+8.*V*X1**(-1)+V*W+3.*V**2-V**2*W-3.*V**2*W**2-8
     1 .*X1**(-1))
     1 +CF*DGQ*(4.+12.*V*X2**(-3)-20.*V*X2**(-2)+12.*V*X2**(-1)-4.*V*W-
     1 12.*V**2*X2**(-3)+16.*V**2*X2**(-2)-4.*V**2*X2**(-1)+4.*V**3*X2*
     1 *(-3)-4.*V**3*X2**(-2)-4.*X2**(-3)+8.*X2**(-2)-8.*X2**(-1))+0.
       CFA11=CFA11+AUCFA11
       CNA5=
     1 +NC*LV*(-34.+46.*V-12.*V*X2**(-3)+28.*V*X2**(-2)-64.*V*X2**(-1)-
     1 20.*V*W+12.*V*W**2-20.*V**2+12.*V**2*X2**(-3)-20.*V**2*X2**(-2)+
     1 44.*V**2*X2**(-1)+32.*V**2*W-24.*V**2*W**2+8.*V**3-4.*V**3*X2**(
     1 -3)+4.*V**3*X2**(-2)-12.*V**3*X2**(-1)-16.*V**3*W+16.*V**3*W**2+
     1 4.*X2**(-3)-12.*X2**(-2)+32.*X2**(-1)+4.*W-4.*W**2)
C
       CNA6=
     1 +NC*LX1*(6.-8.*V+2.*V*W+2.*V**2-2.*V**2*W)
C
       CNA7=
     1 +NC*LX2*(-24.+30.*V+8.*V*W1**(-1)-24.*V*X2**(-3)+56.*V*X2**(-2)-
     1 96.*V*X2**(-1)+14.*V*W-6.*V**2-4.*V**2*W1**(-1)+24.*V**2*X2**(-3
     1 )-40.*V**2*X2**(-2)+64.*V**2*X2**(-1)-14.*V**2*W-8.*V**3*X2**(-3
     1 )+8.*V**3*X2**(-2)-16.*V**3*X2**(-1)+6.*V**(-1)*W1**(-1)-10.*W1*
     1 *(-1)+8.*X2**(-3)-24.*X2**(-2)+48.*X2**(-1))
C
       CNA8=
     1 +NC*LV1*(2.+2.*V+16.*V*X2**(-1)-4.*V**2-12.*V**2*X2**(-1)+4.*V**
     1 3*X2**(-1)-8.*X2**(-1))
C
       CNA9=
     1 +NC*LW*(-18.+28.*V-2.*V*W1**(-1)-16.*V*X2**(-1)-10.*V*W-18.*V**2
     1 -2.*V**2*W1**(-1)+12.*V**2*X2**(-1)+26.*V**2*W-16.*V**2*W**2+8.*
     1 V**3-4.*V**3*X2**(-1)-16.*V**3*W+16.*V**3*W**2+4.*V**(-1)*W1**(-
     1 1)+8.*X2**(-1))
     1 +NC*LW*FQP*(-2.+6.*V-12.*V*W+12.*V*W**2-8.*V**2+16.*V**2*W-16.*V
     1 **2*W**2+4.*V**3-8.*V**3*W+8.*V**3*W**2+4.*W-4.*W**2)
C
       CNA10=
     1 +NC*LW1*(-26.+40.*V-12.*V*X2**(-3)+28.*V*X2**(-2)-48.*V*X2**(-1)
     1 -18.*V*W+12.*V*W**2-22.*V**2+12.*V**2*X2**(-3)-20.*V**2*X2**(-2)
     1 +32.*V**2*X2**(-1)+30.*V**2*W-24.*V**2*W**2+8.*V**3-4.*V**3*X2**
     1 (-3)+4.*V**3*X2**(-2)-8.*V**3*X2**(-1)-16.*V**3*W+16.*V**3*W**2+
     1 4.*X2**(-3)-12.*X2**(-2)+24.*X2**(-1)+4.*W-4.*W**2)
     1 +NC*LW1*FQP*(2.-6.*V+12.*V*W-12.*V*W**2+8.*V**2-16.*V**2*W+16.*V
     1 **2*W**2-4.*V**3+8.*V**3*W-8.*V**3*W**2-4.*W+4.*W**2)
C
       CNA10=CNA10
     1 +NC*(2.-6.*V+12.*V*W-12.*V*W**2+8.*V**2-16.*V**2*W+16.*V
     1 **2*W**2-4.*V**3+8.*V**3*W-8.*V**3*W**2-4.*W+4.*W**2)
       CNA11=
     1 +NC*(8.-6.*V+4.*V*X1**(-1)+10.*V*W1**(-1)*LR+36.*V*X2**(-3)-62.*
     1 V*X2**(-2)+38.*V*X2**(-1)+4.*V*W-12.*V*W**2-6.*V**2-2.*V**2*W1**
     1 (-1)*LR-36.*V**2*X2**(-3)+46.*V**2*X2**(-2)-20.*V**2*X2**(-1)+4.
     1 *V**2*W**2+4.*V**3+12.*V**3*X2**(-3)-10.*V**3*X2**(-2)+2.*V**3*X
     1 2**(-1)+2.*V**(-1)*W1**(-1)*LR-4.*X1**(-1)-10.*W1**(-1)*LR-12.*X
     1 2**(-3)+26.*X2**(-2)-20.*X2**(-1)-4.*W+8.*W**2)
     1 +NC*LSQ*(-22.+34.*V-12.*V*X2**(-3)+28.*V*X2**(-2)-48.*V*X2**(-1)
     1 -12.*V*W+12.*V*W**2-20.*V**2+12.*V**2*X2**(-3)-20.*V**2*X2**(-2)
     1 +32.*V**2*X2**(-1)+24.*V**2*W-24.*V**2*W**2+8.*V**3-4.*V**3*X2**
     1 (-3)+4.*V**3*X2**(-2)-8.*V**3*X2**(-1)-16.*V**3*W+16.*V**3*W**2+
     1 4.*X2**(-3)-12.*X2**(-2)+24.*X2**(-1)+4.*W-4.*W**2)
     1 +NC*FQP*(-36.*V*W+36.*V*W**2+48.*V**2*W-48.*V**2*W**2-24.*V**3*W
     1 +24.*V**3*W**2+12.*W-12.*W**2)+0.
      CA4=CFA5+CFA6+CFA7+CFA8+CFA9+CFA10+CFA11
      CA4=CA4+CNA5+CNA6+CNA7+CNA8+CNA9+CNA10+CNA11
      CAC=DGG*NC*X2**(-1)*(4.*LV*V**3-12.*LV*V**2+12.*LV*V-
     . 4.*LV+4.*LW1*V**3-12.*LW1*V**2+12.*LW1*V-4.*LW1)+DGG*
     . NC*(-8.*LV*V+8.*LV-8.*LW1*V+8.*LW1)
      CAC=CAC-DGG*8.*NC*W*V1*LX2*(V*W/X2+X2/(V*W))/W1
      CA4=CA4+CAC
      RETURN
      END
      FUNCTION CB1(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
c      EXTERNAL FQ1
      EXTERNAL FQ1
      REAL*8 LSQ,LQC,LV,LV1,NF,NC
      COMMON/EC/FQA,FQB,FQC,FQE
      REAL*8 NF4
      NF4=4.
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q1=FQ1(1.d0,X22)
      Q2=FQA
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LQC=LOG(SH/Q1)
      LV=LOG(V)
      LV1=LOG(1.-V)
      CB1=-4./3.*NF+4./3.*V*NF-4./3.*V**2*NF
     1 +LQC*(4./3.*NF4-4./3.*V*NF4+2./3.*V**2*NF4)
     1 +NC*(134./9.-134./9.*V+67./9.*V**2)
     1 +NC*LQC*(-22./3.+22./3.*V-11./3.*V**2)
     1 +NC*LV*(-34./3.+34./3.*V-11./3.*V**2)
     1 +NC*LV*LV1*(14.-16.*V+10.*V**2)
     1 +NC*LV**2*(-14.+14.*V-9.*V**2)
     1 +NC*LV1*(2.-2.*V)
     1 +NC*LV1**2*(-3.+4.*V-3.*V**2)
     1 +CF*(-28.+28.*V-14.*V**2)
     1 +CF*LSQ*(12.-12.*V+6.*V**2)
     1 +CF*LV*(8.-8.*V)
     1 +CF*LV*LSQ*(16.-16.*V+8.*V**2)
      CB1=CB1
     1 +CF*LV*LV1*(-28.+32.*V-20.*V**2)
     1 +CF*LV*LV1*FQQ*(8.-8.*V+4.*V**2)
     1 +CF*LV*FQQ*(6.-6.*V+3.*V**2)
     1 +CF*LV*DQQ*(6.-6.*V+3.*V**2)
     1 +CF*LV**2*(36.-36.*V+22.*V**2)
     1 +CF*LV**2*FQQ*(-4.+4.*V-2.*V**2)
     1 +CF*LV**2*DQQ*(-4.+4.*V-2.*V**2)
     1 +CF*LV1*(2.+4.*V)
     1 +CF*LV1*LSQ*(-8.+8.*V-4.*V**2)
     1 +CF*LV1*FQQ*(-6.+6.*V-3.*V**2)
     1 +CF*LV1**2*(6.-8.*V+6.*V**2)
     1 +CF*LV1**2*FQQ*(-4.+4.*V-2.*V**2)
     1 +CF*FQQ*(18.-18.*V+9.*V**2)
     1 +CF*DQQ*(18.-18.*V+9.*V**2)
       CBP1=
     1 +NC*PI2*(-1./3.+4./3.*V-5./3.*V**2)
     1 +CF*PI2*(2.-4.*V+4.*V**2)
     1 +CF*PI2*FQQ*(4./3.-4./3.*V+2./3.*V**2)
     1 +CF*PI2*DQQ*(-8./3.+8./3.*V-4./3.*V**2)
      CB1=CB1+CBP1
      RETURN
      END
      FUNCTION CB2(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
c      EXTERNAL FQ1
      EXTERNAL FQ1
      REAL*8 LSQ,LV,LV1,NF,NC
      COMMON/EC/FQA,FQB,FQC,FQE
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q2=FQA
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LV=LOG(V)
      LV1=LOG(1.-V)
       CB2=
     1 +NC*(-22./3.+22./3.*V-11./3.*V**2)
     1 +NC*LV*(-16.+16.*V-8.*V**2)
     1 +NC*LV1*(8.-8.*V+4.*V**2)
     1 +CF*LSQ*(16.-16.*V+8.*V**2)
     1 +CF*LV*(48.-48.*V+24.*V**2)
     1 +CF*LV*FQQ*(-8.+8.*V-4.*V**2)
     1 +CF*LV*DQQ*(-8.+8.*V-4.*V**2)
     1 +CF*LV1*(-16.+16.*V-8.*V**2)
     1 +CF*LV1*FQQ*(8.-8.*V+4.*V**2)
     1 +CF*FQQ*(6.-6.*V+3.*V**2)
     1 +CF*DQQ*(6.-6.*V+3.*V**2)
      RETURN
      END
      FUNCTION CB3(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      REAL*8 NC
      PI2=PI**2
      W=WC
      V=VC
       CB3=
     1 +NC*(-8.+8.*V-4.*V**2)
     1 +CF*(32.-32.*V+16.*V**2)
     1 +CF*FQQ*(-8.+8.*V-4.*V**2)
     1 +CF*DQQ*(-8.+8.*V-4.*V**2)
      RETURN
      END
      FUNCTION CB4(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
c      EXTERNAL FQ1
      EXTERNAL FQ1
      REAL*8 LSQ,LQC,NF,NC,LV,LV1,LX1,LX2,LW,LW1,LR
      COMMON/EC/FQA,FQB,FQC,FQE
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q1=FQ1(1.d0,X22)
      Q2=FQA
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LQC=LOG(SH/Q1)
      LV=LOG(V)
      LV1=LOG(1.-V)
      V2=2.-V
      V1=1.-V
      X1=1.-V*W
      X2=1.-V+V*W
      W1=1.-W
      LX1=LOG(X1)
      LX2=LOG(X2)
      LW=LOG(W)
      LW1=LOG(1.-W)
      LR=LOG(X1/V1)
       CB5=
     1 +NC*LV*(8.-8.*V-4.*V*X1**(-1)-12.*V*X2**(-1)+4.*V**2+8.*V**2*X2*
     1 *(-1)+8.*V**2*W**2-4.*V**3*X2**(-1)+8.*X1**(-1)*V2**(-1)+8.*X2**
     1 (-1)*V2**(-1))
     1 +CF*LV*(-28.+28.*V+2.*V*X1**(-2)+6.*V*X1**(-1)+6.*V*X2**(-2)+26.
     1 *V*X2**(-1)-14.*V**2-6.*V**2*X2**(-2)-20.*V**2*X2**(-1)-20.*V**2
     1 *W**2+2.*V**3*X2**(-2)+10.*V**3*X2**(-1)-2.*X1**(-2)-16.*X1**(-1
     1 )*V2**(-1)-2.*X2**(-2)-16.*X2**(-1)*V2**(-1))
     1 +CF*LV*FQQ*(6.-8.*V-2.*V*X1**(-2)+2.*V*X1**(-1)-4.*V*W+4.*V**2+2
     1 .*V**2*W+2.*V**2*W**2+2.*X1**(-2))
     1 +CF*LV*DQQ*(6.-4.*V-6.*V*X2**(-2)-2.*V*X2**(-1)+4.*V*W+2.*V**2+6
     1 .*V**2*X2**(-2)+4.*V**2*X2**(-1)-2.*V**2*W+2.*V**2*W**2-2.*V**3*
     1 X2**(-2)-2.*V**3*X2**(-1)+2.*X2**(-2))
C
       CB6=
     1 +NC*LX1*(-4.+4.*V-2.*V**2-4.*V**2*W**2)
     1 +CF*LX1*(16.-16.*V+8.*V**2+8.*V**2*W**2)
C
       CB7=
     1 +NC*LX2*(-4.+4.*V-4.*V*W1**(-1)-2.*V**2+2.*V**2*W1**(-1)-4.*V**2
     1 *W**2+4.*W1**(-1))
     1 +CF*LX2*(4.-8.*V+12.*V*X2**(-2)+4.*V*X2**(-1)-8.*V*W+4.*V**2-12.
     1 *V**2*X2**(-2)-8.*V**2*X2**(-1)+4.*V**2*W+4.*V**2*W**2+4.*V**3*X
     1 2**(-2)+4.*V**3*X2**(-1)-4.*X2**(-2))
     1 +CF*LX2*DQQ*(12.-8.*V+16.*V*W1**(-1)-12.*V*X2**(-2)-4.*V*X2**(-1
     1 )+8.*V*W+4.*V**2-8.*V**2*W1**(-1)+12.*V**2*X2**(-2)+8.*V**2*X2**
     1 (-1)-4.*V**2*W+4.*V**2*W**2-4.*V**3*X2**(-2)-4.*V**3*X2**(-1)-16
     1 .*W1**(-1)+4.*X2**(-2))
C
       CB8=
     1 +NC*LV1*(2.*V*X1**(-1)+6.*V*X2**(-1)-4.*V**2*X2**(-1)+2.*V**3*X2
     1 **(-1)-4.*X1**(-1)*V2**(-1)-4.*X2**(-1)*V2**(-1))
     1 +CF*LV1*(-8.+8.*V-4.*V*X1**(-1)-12.*V*X2**(-1)-4.*V**2+8.*V**2*X
     1 2**(-1)-4.*V**3*X2**(-1)+8.*X1**(-1)*V2**(-1)+8.*X2**(-1)*V2**(-
     1 1))
     1 +CF*LV1*FQQ*(-6.+8.*V+2.*V*X1**(-2)-2.*V*X1**(-1)+4.*V*W-4.*V**2
     1 -2.*V**2*W-2.*V**2*W**2-2.*X1**(-2))
C
       CB9=
     1 +NC*LW*(4.-4.*V-2.*V*X1**(-1)+8.*V*W1**(-1)-6.*V*X2**(-1)+2.*V**
     1 2-4.*V**2*W1**(-1)+4.*V**2*X2**(-1)+4.*V**2*W**2-2.*V**3*X2**(-1
     1 )+4.*X1**(-1)*V2**(-1)-8.*W1**(-1)+4.*X2**(-1)*V2**(-1))
     1 +CF*LW*(-8.+8.*V+4.*V*X1**(-1)-16.*V*W1**(-1)+12.*V*X2**(-1)-4.*
     1 V**2+8.*V**2*W1**(-1)-8.*V**2*X2**(-1)-8.*V**2*W**2+4.*V**3*X2**
     1 (-1)-8.*X1**(-1)*V2**(-1)+16.*W1**(-1)-8.*X2**(-1)*V2**(-1))
C
       CB10=
     1 +NC*LW1*(4.-4.*V-2.*V*X1**(-1)-6.*V*X2**(-1)+2.*V**2+4.*V**2*X2*
     1 *(-1)+4.*V**2*W**2-2.*V**3*X2**(-1)+4.*X1**(-1)*V2**(-1)+4.*X2**
     1 (-1)*V2**(-1))
     1 +CF*LW1*(-20.+20.*V+2.*V*X1**(-2)+2.*V*X1**(-1)+6.*V*X2**(-2)+14
     1 .*V*X2**(-1)-10.*V**2-6.*V**2*X2**(-2)-12.*V**2*X2**(-1)-12.*V**
     1 2*W**2+2.*V**3*X2**(-2)+6.*V**3*X2**(-1)-2.*X1**(-2)-8.*X1**(-1)
     1 *V2**(-1)-2.*X2**(-2)-8.*X2**(-1)*V2**(-1))
     1 +CF*LW1*FQQ*(6.-8.*V-2.*V*X1**(-2)+2.*V*X1**(-1)-4.*V*W+4.*V**2+
     1 2.*V**2*W+2.*V**2*W**2+2.*X1**(-2))
     1 +CF*LW1*DQQ*(6.-4.*V-6.*V*X2**(-2)-2.*V*X2**(-1)+4.*V*W+2.*V**2+
     1 6.*V**2*X2**(-2)+4.*V**2*X2**(-1)-2.*V**2*W+2.*V**2*W**2-2.*V**3
     1 *X2**(-2)-2.*V**3*X2**(-1)+2.*X2**(-2))
C
       CB11=
     1 +NC*(-8./3.+8./3.*V-V*X1**(-2)-V*X1**(-1)-4.*V*W1**(-1)*LR-3.*V*
     1 X2**(-2)-11.*V*X2**(-1)+2./3.*V**2+2.*V**2*W1**(-1)*LR+3.*V**2*X
     1 2**(-2)+10.*V**2*X2**(-1)+4./3.*V**2*W**2-V**3*X2**(-2)-3.*V**3*
     1 X2**(-1)+X1**(-2)+4.*X1**(-1)+4.*W1**(-1)*LR+X2**(-2)+4.*X2**(-1
     1 ))
     1 +CF*(12.-12.*V-2.*V*X1**(-2)+8.*V*X1**(-1)+16.*V*W1**(-1)*LR-6.*
     1 V*X2**(-2)+16.*V*X2**(-1)-8.*V**2*W1**(-1)*LR+6.*V**2*X2**(-2)-8
     1 .*V**2*X2**(-1)+4.*V**2*W-2.*V**3*X2**(-2)+2.*X1**(-2)-8.*X1**(-
     1 1)-16.*W1**(-1)*LR+2.*X2**(-2)-8.*X2**(-1))
     1 +CF*LSQ*(-12.+12.*V+2.*V*X1**(-2)-2.*V*X1**(-1)+6.*V*X2**(-2)+2.
     1 *V*X2**(-1)-6.*V**2-6.*V**2*X2**(-2)-4.*V**2*X2**(-1)-4.*V**2*W*
     1 *2+2.*V**3*X2**(-2)+2.*V**3*X2**(-1)-2.*X1**(-2)-2.*X2**(-2))
     1 +CF*FQQ*(6.*V+4.*V*X1**(-2)-4.*V*X1**(-1)-4.*V*W-3.*V**2+V**2*W+
     1 3.*V**2*W**2-4.*X1**(-2)-2.*X1**(-1))
     1 +CF*DQQ*(-6.*V+18.*V*X2**(-1)+3.*V**2-18.*V**2*X2**(-1)-3.*V**2*
     1 W+3.*V**2*W**2+6.*V**3*X2**(-1)-6.*X2**(-1))+0.
      CB4=CB5+CB6+CB7+CB8+CB9+CB10+CB11
      RETURN
      END
      FUNCTION CC1(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/FSM/FGG
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
c      EXTERNAL FQ1
      EXTERNAL FQ1
      REAL*8 LSQ,LQC,LV,LV1,NF,NC
      COMMON/EC/FQA,FQB,FQC,FQE
      REAL*8 NF4
      CGG=1.5
      NF4=4.
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q1=FQ1(1.d0,X22)
      Q2=FQC
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LQC=LOG(SH/Q1)
      LV=LOG(V)
      LV1=LOG(1.-V)
       CC1=
     1 +NF4*(-2./3.*LQC+2./3.*V*LQC+1./3.*V**(-1)*LQC)
     1 +NF*LSQ*(2./3.-2./3.*V-1./3.*V**(-1))
     1 +NC*(11./3.*LQC-11./3.*V*LQC-11./6.*V**(-1)*LQC)
     1 +NC*LV*(-1.+V)
     1 +NC*LV*LV1*(8.-8.*V-4.*V**(-1))
     1 +NC*LV*LSQ*(-4.+4.*V+2.*V**(-1))
     1 +NC*LV**2*(-9.+17./2.*V+4.*V**(-1))
     1 +NC*LV1*(-1.+V)+NC*LV1*LSQ*(4.-4.*V-2.*V**(-1))
     1 +NC*LV1**2*(1./2.*V-1./2.*V**(-1))
     1 +NC*LSQ*(-11./3.+11./3.*V+11./6.*V**(-1))
     1 +CF*(7.-7.*V+9.*V*DQQ-4./3.*V*DQQ*PI2+4./3.*V*PI2-7./2.*V**(-1)+
     1 9./2.*V**(-1)*DQQ-2./3.*V**(-1)*DQQ*PI2+2./3.*V**(-1)*PI2-9.*DQQ
     1 +4./3.*DQQ*PI2-4./3.*PI2)
     1 +CF*LV*(-1.-2.*V+3.*V*DQQ+3./2.*V**(-1)+3./2.*V**(-1)*DQQ-3.*DQQ
     1 )+CF*LV*LSQ*(-4.+4.*V+2.*V**(-1))
     1 +CF*LV**2*(-2.+3.*V-2.*V*DQQ+2.*V**(-1)-V**(-1)*DQQ+2.*DQQ)
     1 +CF*LV1*(2.+V)+CF*LV1**2*(-2.+V+2.*V**(-1))
     1 +CF*LSQ*(-3.+3.*V+3./2.*V**(-1))
      R3D=LV1**2*V+1./2.*LV1**2*V**(-1)-LV1**2-2.*LV1*LV*V-
     . LV1*LV*V**(-1)+2.*LV1*LV+LV**2*V+1./2.*LV**2*V**(-1)
     . -LV**2-2.*CGG*V-CGG*V**(-1)+2.*CGG
      CC1=CC1+FGG*R3D*(-2*NC)
      RETURN
      END
      FUNCTION CC2(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/FSM/FGG
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
c      EXTERNAL FQ1
      EXTERNAL FQ1
      REAL*8 LSQ,LV,LV1,NF,NC
      COMMON/EC/FQA,FQB,FQC,FQE
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q2=FQC
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LV=LOG(V)
      LV1=LOG(1.-V)
       CC2=
     1 +NC*LV*(-12.+12.*V+6.*V**(-1))
     1 +NC*LV1*(4.-4.*V-2.*V**(-1))
     1 +NC*LSQ*(-4.+4.*V+2.*V**(-1))
     1 +CF*(3.-3.*V+3.*V*DQQ-3./2.*V**(-1)+3./2.*V**(-1)*DQQ-3.*DQQ)
     1 +CF*LV*(-4.+4.*V-4.*V*DQQ+2.*V**(-1)-2.*V**(-1)*DQQ+4.*DQQ)
     1 +CF*LSQ*(-4.+4.*V+2.*V**(-1))
      R3W=-2.*LV1*V-LV1*V**(-1)+2.*LV1+2.*LV*V+LV*V**(-1)-
     . 2.*LV
      CC2=CC2+FGG*R3W*(-2*NC)
      RETURN
      END
      FUNCTION CC3(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/FSM/FGG
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      REAL*8 NC
      PI2=PI**2
      W=WC
      V=VC
       CC3=
     1 +NC*(-8.+8.*V+4.*V**(-1))
     1 +CF*(-4.+4.*V-4.*V*DQQ+2.*V**(-1)-2.*V**(-1)*DQQ+4.*DQQ)
      R3L=2.*V+V**(-1)-2.
      CC3=CC3+FGG*R3L*(-2*NC)
      RETURN
      END
      FUNCTION CC4(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/FSM/FGG
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
c      EXTERNAL FQ1
      EXTERNAL FQ1
      REAL*8 LSQ,LQC,NF,NC,LV,LV1,LX1,LX2,LW,LW1,LR
      COMMON/EC/FQA,FQB,FQC,FQE
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q1=FQ1(1.d0,X22)
      Q2=FQC
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LQC=LOG(SH/Q1)
      LV=LOG(V)
      LV1=LOG(1.-V)
      V2=2.-V
      V1=1.-V
      X1=1.-V*W
      X2=1.-V+V*W
      W1=1.-W
      LX1=LOG(X1)
      LX2=LOG(X2)
      LW=LOG(W)
      LW1=LOG(1.-W)
      LR=LOG(X1/V1)
      CC11=
     1 +NC*(-6.+3.*V-12.*V*X1**(-3)+21.*V*X1**(-2)-11.*V*X1**(-1)-5.*V*
     1 W1**(-1)*LR-4.*V*X2**(-1)-2.*V*W+2.*V*W**2+6.*V**2*X1**(-3)-8.*V
     1 **2*X1**(-2)+2.*V**2*X1**(-1)+2.*V**2*X2**(-1)-V**(-1)*W1**(-1)*
     1 LR+2.*V1**(-1)+6.*X1**(-3)-13.*X1**(-2)+10.*X1**(-1)+4.*W1**(-1)
     1 *LR+2.*X2**(-1)+2.*W-4.*W**2)
     1 +NC*LSQ*(7.-9.*V+4.*V*X1**(-3)-10.*V*X1**(-2)+12.*V*X1**(-1)+6.*
     1 V*W-6.*V*W**2-2.*V**2*X1**(-3)+4.*V**2*X1**(-2)-4.*V**2*X1**(-1)
     1 +4.*V1**(-1)-8.*V1**(-1)*W+8.*V1**(-1)*W**2-2.*X1**(-3)+6.*X1**(
     1 -2)-12.*X1**(-1)+6.*W-6.*W**2)
     1 +NC*FQP*(6.*V*W-6.*V*W**2-12.*V1**(-1)*W+12.*V1**(-1)*W**2+6.*W-
     1 6.*W**2)
     1 +CF*(-9./2.+2.*V+12.*V*X1**(-3)-17.*V*X1**(-2)+8.*V*X1**(-1)-9./
     1 2.*V*DQQ-3./2.*V*DQQ*W-2.*V*W1**(-1)*LR+4.*V*X2**(-1)+3.*V*W-V*W
     1 **2-2.*V**2-6.*V**2*X1**(-3)+6.*V**2*X1**(-2)-V**2*X1**(-1)+3./2
     1 .*V**2*DQQ+3./2.*V**2*DQQ*W**2-2.*V**2*X2**(-1)-2.*V**2*W**2+3.*
     1 V**(-1)-4.*V**(-1)*W1**(-1)*LR+5.*V**(-1)*W-6.*V**(-1)*W**2-6.*X
     1 1**(-3)+11.*X1**(-2)-4.*X1**(-1)+3./2.*DQQ+4.*W1**(-1)*LR-2.*X2*
     1 *(-1)-6.*W+6.*W**2)
      CC11=CC11
     1 +CF*LSQ*(-6.+6.*V-4.*V*X1**(-3)+6.*V*X1**(-2)-6.*V*X1**(-1)-4.*V
     1 *X2**(-1)-4.*V*W+2.*V*W**2-4.*V**2+2.*V**2*X1**(-3)-2.*V**2*X1**
     1 (-2)+2.*V**2*X1**(-1)+2.*V**2*X2**(-1)+4.*V**2*W-4.*V**2*W**2+V*
     1 *(-1)-2.*V**(-1)*W+2.*V**(-1)*W**2+2.*X1**(-3)-4.*X1**(-2)+5.*X1
     1 **(-1)+2.*X2**(-1)+2.*W-2.*W**2)
     1 +CF*FQP*(-6.*V*W+6.*V*W**2+6.*V**2*W-6.*V**2*W**2-6.*V**(-1)*W+6
     1 .*V**(-1)*W**2+6.*W-6.*W**2)
     1 +CF*FQG*(-6.+12.*V-12.*V*X1**(-3)+18.*V*X1**(-2)-18.*V*X1**(-1)-
     1 6.*V*W-6.*V**2+6.*V**2*X1**(-3)-6.*V**2*X1**(-2)+6.*V**2*X1**(-1
     1 )+6.*V**2*W+6.*X1**(-3)-12.*X1**(-2)+12.*X1**(-1))
       CC5=
     1 +NC*LV*(13.-15.*V+4.*V*X1**(-3)-10.*V*X1**(-2)+16.*V*X1**(-1)+2.
     1 *V*W-6.*V*W**2-2.*V**2*X1**(-3)+4.*V**2*X1**(-2)-6.*V**2*X1**(-1
     1 )+4.*V1**(-1)-8.*V1**(-1)*W+8.*V1**(-1)*W**2-2.*X1**(-3)+6.*X1**
     1 (-2)-16.*X1**(-1)+6.*W-6.*W**2)
     1 +CF*LV*(-6.+6.*V-4.*V*X1**(-3)+6.*V*X1**(-2)-6.*V*X1**(-1)-V*DQQ
     1 +4.*V*DQQ*X2**(-1)+V*DQQ*W-4.*V*X2**(-1)-4.*V*W+2.*V*W**2-4.*V**
     1 2+2.*V**2*X1**(-3)-2.*V**2*X1**(-2)+2.*V**2*X1**(-1)+V**2*DQQ-2.
     1 *V**2*DQQ*X2**(-1)+V**2*DQQ*W**2+2.*V**2*X2**(-1)+4.*V**2*W-4.*V
     1 **2*W**2+V**(-1)-2.*V**(-1)*W+2.*V**(-1)*W**2+2.*X1**(-3)-4.*X1*
     1 *(-2)+5.*X1**(-1)+2.*DQQ-2.*DQQ*X2**(-1)+2.*X2**(-1)+2.*W-2.*W**
     1 2)
     1 +CF*LV*FQG*(3.-4.*V+4.*V*X1**(-3)-6.*V*X1**(-2)+6.*V*X1**(-1)+V*
     1 W+2.*V**2-2.*V**2*X1**(-3)+2.*V**2*X1**(-2)-2.*V**2*X1**(-1)-2.*
     1 V**2*W+V**2*W**2-2.*X1**(-3)+4.*X1**(-2)-5.*X1**(-1))
C
       CC8=
     1 +NC*LV1*(2.-6.*V+4.*V*W-8.*V*W**2+4.*V1**(-1)-8.*V1**(-1)*W+8.*V
     1 1**(-1)*W**2+8.*W-8.*W**2)
     1 +CF*LV1*(-4.+4.*V-2.*V**2)
     1 +CF*LV1*FQG*(-3.+4.*V-4.*V*X1**(-3)+6.*V*X1**(-2)-6.*V*X1**(-1)-
     1 V*W-2.*V**2+2.*V**2*X1**(-3)-2.*V**2*X1**(-2)+2.*V**2*X1**(-1)+2
     1 .*V**2*W-V**2*W**2+2.*X1**(-3)-4.*X1**(-2)+5.*X1**(-1))
C
       CC10=
     1 +NC*LW1*(9.-10.*V+4.*V*X1**(-3)-10.*V*X1**(-2)+12.*V*X1**(-1)+3.
     1 *V*W-6.*V*W**2-2.*V**2*X1**(-3)+4.*V**2*X1**(-2)-4.*V**2*X1**(-1
     1 )+4.*V1**(-1)-8.*V1**(-1)*W+8.*V1**(-1)*W**2-2.*X1**(-3)+6.*X1**
     1 (-2)-12.*X1**(-1)+6.*W-6.*W**2)
     1 +NC*LW1*FQP*(1.+V-2.*V*W+2.*V*W**2-2.*V1**(-1)+4.*V1**(-1)*W-4.*
     1 V1**(-1)*W**2-2.*W+2.*W**2)
     1 +CF*LW1*(-8.+6.*V-4.*V*X1**(-3)+6.*V*X1**(-2)-6.*V*X1**(-1)-V*DQ
     1 Q+4.*V*DQQ*X2**(-1)+V*DQQ*W-4.*V*X2**(-1)-4.*V*W+2.*V*W**2-6.*V*
     1 *2+2.*V**2*X1**(-3)-2.*V**2*X1**(-2)+2.*V**2*X1**(-1)+V**2*DQQ-2
     1 .*V**2*DQQ*X2**(-1)+V**2*DQQ*W**2+2.*V**2*X2**(-1)+8.*V**2*W-6.*
     1 V**2*W**2+V**(-1)-2.*V**(-1)*W+2.*V**(-1)*W**2+2.*X1**(-3)-4.*X1
     1 **(-2)+5.*X1**(-1)+2.*DQQ-2.*DQQ*X2**(-1)+2.*X2**(-1)+2.*W-2.*W*
     1 *2)
     1 +CF*LW1*FQP*(1.-V+2.*V*W-2.*V*W**2+V**2-2.*V**2*W+2.*V**2*W**2-V
     1 **(-1)+2.*V**(-1)*W-2.*V**(-1)*W**2-2.*W+2.*W**2)
     1 +CF*LW1*FQG*(3.-4.*V+4.*V*X1**(-3)-6.*V*X1**(-2)+6.*V*X1**(-1)+V
     1 *W+2.*V**2-2.*V**2*X1**(-3)+2.*V**2*X1**(-2)-2.*V**2*X1**(-1)-2.
     1 *V**2*W+V**2*W**2-2.*X1**(-3)+4.*X1**(-2)-5.*X1**(-1))
C
       CC10=CC10
     1 +CF*(3.-4.*V+4.*V*X1**(-3)-6.*V*X1**(-2)+6.*V*X1**(-1)+V
     1 *W+2.*V**2-2.*V**2*X1**(-3)+2.*V**2*X1**(-2)-2.*V**2*X1**(-1)-2.
     1 *V**2*W+V**2*W**2-2.*X1**(-3)+4.*X1**(-2)-5.*X1**(-1))
C
       CC10=CC10
     1 +NC*(1.+V-2.*V*W+2.*V*W**2-2.*V1**(-1)+4.*V1**(-1)*W-4.*
     1 V1**(-1)*W**2-2.*W+2.*W**2)
     1 +CF*(1.-V+2.*V*W-2.*V*W**2+V**2-2.*V**2*W+2.*V**2*W**2-V
     1 **(-1)+2.*V**(-1)*W-2.*V**(-1)*W**2-2.*W+2.*W**2)
C
       CC6=
     1 +NC*LX1*(-4.+7.*V-V*W+8.*V*W**2-4.*V1**(-1)+8.*V1**(-1)*W-8.*V1*
     1 *(-1)*W**2-8.*W+8.*W**2)
     1 +CF*LX1*(6.-4.*V+4.*V**2-4.*V**2*W+2.*V**2*W**2)
C
       CC7=
     1 +NC*LX2*(-3.+2.*V-2.*V*W1**(-1)+V*W+V**(-1)*W1**(-1)+2.*W1**(-1)
     1 )
     1 +CF*LX2*(2.-2.*V*DQQ-8.*V*DQQ*W1**(-1)+8.*V*DQQ*X2**(-1)+2.*V*DQ
     1 Q*W+4.*V*W1**(-1)-8.*V*X2**(-1)+2.*V**2*DQQ-4.*V**2*DQQ*X2**(-1)
     1 +2.*V**2*DQQ*W**2+4.*V**2*X2**(-1)-4.*V**(-1)*DQQ*W1**(-1)-2.*V*
     1 *(-1)*W1**(-1)+4.*DQQ+8.*DQQ*W1**(-1)-4.*DQQ*X2**(-1)-4.*W1**(-1
     1 )+4.*X2**(-1))
C
       CC9=
     1 +NC*LW*(5.-8.*V+4.*V*X1**(-1)+5.*V*W1**(-1)+3.*V*W-8.*V*W**2-2.*
     1 V**2*X1**(-1)+2.*V**(-1)*W1**(-1)+4.*V1**(-1)-8.*V1**(-1)*W+8.*V
     1 1**(-1)*W**2-4.*X1**(-1)-6.*W1**(-1)+8.*W-8.*W**2)
     1 +NC*LW*FQP*(-1.-V+2.*V*W-2.*V*W**2+2.*V1**(-1)-4.*V1**(-1)*W+4.*
     1 V1**(-1)*W**2+2.*W-2.*W**2)
     1 +CF*LW*(-2.-2.*V+2.*V*W1**(-1)-2.*V*W-2.*V**2*W**2+2.*V**(-1)*W1
     1 **(-1))
     1 +CF*LW*FQP*(-1.+V-2.*V*W+2.*V*W**2-V**2+2.*V**2*W-2.*V**2*W**2+V
     1 **(-1)-2.*V**(-1)*W+2.*V**(-1)*W**2+2.*W-2.*W**2)+0.
      CC4=CC5+CC6+CC7+CC8+CC9+CC10+CC11
      R3C=-LW1*X1**(-1)-2.*LW1*V+2.*LW1-LV*X1**(-1)-2.*LV*V
     . +2.*LV-2.*LR*(-W+1.)**(-1)*V-LR*(-W+1.)**(-1)*V**(-1
     . )+2.*LR*(-W+1.)**(-1)+LX1*X1**(-1)+2.*LX1*V-2.*LX1
      R3C=R3C+(LX1-LV1)*(X1+(V*W)**2/X1)*V1/V/W1
      CC4=CC4+FGG*R3C*(-2*NC)
      RETURN
      END
      FUNCTION CD4(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      EXTERNAL FQ1
      REAL*8 LSQ,LQC,NF,NC,LV,LV1,LX1,LX2,LW,LW1,LR
      COMMON/EC/FQA,FQB,FQC,FQE
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q1=FQ1(1.d0,X22)
      Q2=FQE
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LQC=LOG(SH/Q1)
      LV=LOG(V)
      LV1=LOG(1.-V)
      V2=2.-V
      V1=1.-V
      X1=1.-V*W
      X2=1.-V+V*W
      W1=1.-W
      LX1=LOG(X1)
      LX2=LOG(X2)
      LW=LOG(W)
      LW1=LOG(1.-W)
      LR=LOG(X1/V1)
      CD5=
     1 +NC*LV*(14.-10.*V+4.*V*W-4.*V*W**2+4.*V*X1**(-1)+4.*V*X2**(-1)-4
     1 .*V**2*X1**(-1)-4.*V**2*X2**(-1)+4.*V**(-2)-8.*V**(-2)*W+8.*V**(
     1 -2)*W**2-8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W**2+8.*V2**(-1)*X
     1 1**(-1)+8.*V2**(-1)*X2**(-1)-12.*W+12.*W**2-8.*X1**(-1)-8.*X2**(
     1 -1))
     1 +CF*LV*(8.-8.*V+8.*V*W+4.*V*W*DGQ-8.*V*W**2-16.*V*X1**(-2)+16.*V
     1 *X1**(-1)-16.*V*X2**(-2)+16.*V*X2**(-2)*DGQ+16.*V*X2**(-1)-24.*V
     1 *X2**(-1)*DGQ+8.*V*DGQ+8.*V**2-8.*V**2*W+8.*V**2*W**2-2.*V**2*W*
     1 *2*DGQ+8.*V**2*X1**(-2)+8.*V**2*X2**(-2)-8.*V**2*X2**(-2)*DGQ+8.
     1 *V**2*X2**(-1)*DGQ-2.*V**2*DGQ-16.*V2**(-1)*X1**(-1)-16.*V2**(-1
     1 )*X2**(-1)-8.*W+8.*W**2+8.*X1**(-2)+8.*X2**(-2)-8.*X2**(-2)*DGQ+
     1 16.*X2**(-1)*DGQ-10.*DGQ)
     1 +CF*LV*FQG*(-10.+12.*V-4.*V*W+16.*V*X1**(-2)-24.*V*X1**(-1)-4.*V
     1 **2+4.*V**2*W-2.*V**2*W**2-8.*V**2*X1**(-2)+8.*V**2*X1**(-1)-8.*
     1 X1**(-2)+16.*X1**(-1))
C
       CD6=
     1 +CF*LX1*(-4.+8.*V-8.*V*W-8.*V**2+8.*V**2*W-4.*V**2*W**2)
C
       CD7=
     1 +CF*LX2*(16.-16.*V+8.*V*W*DGQ-32.*V*X2**(-2)+32.*V*X2**(-2)*DGQ+
     1 48.*V*X2**(-1)-48.*V*X2**(-1)*DGQ+16.*V*DGQ-4.*V**2*W**2*DGQ+16.
     1 *V**2*X2**(-2)-16.*V**2*X2**(-2)*DGQ-16.*V**2*X2**(-1)+16.*V**2*
     1 X2**(-1)*DGQ-4.*V**2*DGQ+16.*X2**(-2)-16.*X2**(-2)*DGQ-32.*X2**(
     1 -1)+32.*X2**(-1)*DGQ-20.*DGQ)
C
       CD8=
     1 +NC*LV1*(2.+2.*V+4.*V*X1**(-1)+4.*V*X2**(-1)-4.*V**(-1)+4.*V**(-
     1 1)*X1**(-1)+4.*V**(-1)*X2**(-1)-4.*V2**(-1)*X1**(-1)-4.*V2**(-1)
     1 *X2**(-1)-4.*X1**(-1)-4.*X2**(-1))
     1 +CF*LV1*(-12.-8.*V*W+8.*V*W**2-8.*V*X1**(-1)-8.*V*X2**(-1)+4.*V*
     1 *2+16.*V**(-1)-8.*V**(-1)*X1**(-1)-8.*V**(-1)*X2**(-1)+8.*V2**(-
     1 1)*X1**(-1)+8.*V2**(-1)*X2**(-1)+8.*X1**(-1)+8.*X2**(-1))
     1 +CF*LV1*FQG*(10.-12.*V+4.*V*W-16.*V*X1**(-2)+24.*V*X1**(-1)+4.*V
     1 **2-4.*V**2*W+2.*V**2*W**2+8.*V**2*X1**(-2)-8.*V**2*X1**(-1)+8.*
     1 X1**(-2)-16.*X1**(-1))
C
       CD9=
     1 +NC*LW*(2.-2.*V+4.*V*W+8.*V*X1**(-1)-4.*V*X2**(-1)-4.*V**2*X1**(
     1 -1)+4.*V**(-1)*X1**(-1)-4.*V**(-1)*X2**(-1)+4.*V2**(-1)*X1**(-1)
     1 +4.*V2**(-1)*X2**(-1)-4.*W-12.*X1**(-1)+4.*X2**(-1))
     1 +NC*LW*FQP*(6.-2.*V+4.*V*W-4.*V*W**2+4.*V**(-2)-8.*V**(-2)*W+8.*
     1 V**(-2)*W**2-8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W**2-12.*W+12.
     1 *W**2)
     1 +CF*LW*(-8.+8.*V-8.*V*W-16.*V*X1**(-1)+8.*V*X2**(-1)+4.*V**2*W**
     1 2+8.*V**2*X1**(-1)-8.*V**(-1)*X1**(-1)+8.*V**(-1)*X2**(-1)-8.*V2
     1 **(-1)*X1**(-1)-8.*V2**(-1)*X2**(-1)+8.*W+24.*X1**(-1)-8.*X2**(-
     1 1))
     1 +CF*LW*FQP*(4.-4.*V+8.*V*W-8.*V*W**2+2.*V**2-4.*V**2*W+4.*V**2*W
     1 **2-8.*W+8.*W**2)
C
       CD10=
     1 +NC*LW1*(4.-4.*V*W**2-4.*V*X1**(-1)+8.*V*X2**(-1)-4.*V**2*X2**(-
     1 1)+4.*V**(-2)-8.*V**(-2)*W+8.*V**(-2)*W**2-8.*V**(-1)+16.*V**(-1
     1 )*W-16.*V**(-1)*W**2-4.*V**(-1)*X1**(-1)+4.*V**(-1)*X2**(-1)+4.*
     1 V2**(-1)*X1**(-1)+4.*V2**(-1)*X2**(-1)-8.*W+12.*W**2+4.*X1**(-1)
     1 -12.*X2**(-1))
     1 +NC*LW1*FQP*(-6.+2.*V-4.*V*W+4.*V*W**2-4.*V**(-2)+8.*V**(-2)*W-8
     1 .*V**(-2)*W**2+8.*V**(-1)-16.*V**(-1)*W+16.*V**(-1)*W**2+12.*W-1
     1 2.*W**2)
     1 +CF*LW1*(24.-24.*V+16.*V*W+4.*V*W*DGQ-8.*V*W**2-16.*V*X1**(-2)+3
     1 2.*V*X1**(-1)-16.*V*X2**(-2)+16.*V*X2**(-2)*DGQ+8.*V*X2**(-1)-24
     1 .*V*X2**(-1)*DGQ+8.*V*DGQ+12.*V**2-16.*V**2*W+12.*V**2*W**2-2.*V
     1 **2*W**2*DGQ+8.*V**2*X1**(-2)-8.*V**2*X1**(-1)+8.*V**2*X2**(-2)-
     1 8.*V**2*X2**(-2)*DGQ+8.*V**2*X2**(-1)*DGQ-2.*V**2*DGQ+8.*V**(-1)
     1 *X1**(-1)-8.*V**(-1)*X2**(-1)-8.*V2**(-1)*X1**(-1)-8.*V2**(-1)*X
     1 2**(-1)-16.*W+8.*W**2+8.*X1**(-2)-24.*X1**(-1)+8.*X2**(-2)-8.*X2
     1 **(-2)*DGQ+8.*X2**(-1)+16.*X2**(-1)*DGQ-10.*DGQ)
       AUCD10=
     1 +CF*LW1*FQP*(-4.+4.*V-8.*V*W+8.*V*W**2-2.*V**2+4.*V**2*W-4.*V**2
     1 *W**2+8.*W-8.*W**2)
     1 +CF*LW1*FQG*(-10.+12.*V-4.*V*W+16.*V*X1**(-2)-24.*V*X1**(-1)-4.*
     1 V**2+4.*V**2*W-2.*V**2*W**2-8.*V**2*X1**(-2)+8.*V**2*X1**(-1)-8.
     1 *X1**(-2)+16.*X1**(-1))
       CD10=CD10+AUCD10
C
       CD10=CD10
     1 +NC*(-6.+2.*V-4.*V*W+4.*V*W**2-4.*V**(-2)+8.*V**(-2)*W-8
     1 .*V**(-2)*W**2+8.*V**(-1)-16.*V**(-1)*W+16.*V**(-1)*W**2+12.*W-1
     1 2.*W**2)
     1 +CF*(-4.+4.*V-8.*V*W+8.*V*W**2-2.*V**2+4.*V**2*W-4.*V**2
     1 *W**2+8.*W-8.*W**2)
C
       CD10=CD10
     1 +CF*(-10.+12.*V-4.*V*W+16.*V*X1**(-2)-24.*V*X1**(-1)-4.*
     1 V**2+4.*V**2*W-2.*V**2*W**2-8.*V**2*X1**(-2)+8.*V**2*X1**(-1)-8.
     1 *X1**(-2)+16.*X1**(-1))
C
       CD11=
     1 +NC*(-6.+6.*V-4.*V*W+4.*V*W**2+16.*V*X1**(-2)-12.*V*X1**(-1)+16.
     1 *V*X2**(-2)-12.*V*X2**(-1)-8.*V**2*X1**(-2)+4.*V**2*X1**(-1)-8.*
     1 V**2*X2**(-2)+4.*V**2*X2**(-1)+24.*V**(-2)*W-24.*V**(-2)*W**2-48
     1 .*V**(-1)*W+48.*V**(-1)*W**2+28.*W-28.*W**2-8.*X1**(-2)+8.*X1**(
     1 -1)-8.*X2**(-2)+8.*X2**(-1))
     1 +NC*LSQ*(6.-2.*V+4.*V*W-4.*V*W**2+4.*V**(-2)-8.*V**(-2)*W+8.*V**
     1 (-2)*W**2-8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W**2-12.*W+12.*W*
     1 *2)
     1 +NC*FQP*(12.*V*W-12.*V*W**2-24.*V**(-2)*W+24.*V**(-2)*W**2+48.*V
     1 **(-1)*W-48.*V**(-1)*W**2-36.*W+36.*W**2)
     1 +CF*(8.-8.*V-4.*V*W*DGQ-16.*V*X2**(-2)*DGQ+24.*V*X2**(-1)*DGQ-4.
     1 *V*DGQ+4.*V**2+4.*V**2*W**2+8.*V**2*X2**(-2)*DGQ-8.*V**2*X2**(-1
     1 )*DGQ+8.*X2**(-2)*DGQ-16.*X2**(-1)*DGQ+8.*DGQ)
     1 +CF*LSQ*(24.-24.*V+8.*V*W-8.*V*W**2-16.*V*X1**(-2)+24.*V*X1**(-1
     1 )-16.*V*X2**(-2)+24.*V*X2**(-1)+8.*V**2-8.*V**2*W+8.*V**2*W**2+8
     1 .*V**2*X1**(-2)-8.*V**2*X1**(-1)+8.*V**2*X2**(-2)-8.*V**2*X2**(-
     1 1)-8.*W+8.*W**2+8.*X1**(-2)-16.*X1**(-1)+8.*X2**(-2)-16.*X2**(-1
     1 ))
       AUCD11=
     1 +CF*FQP*(24.*V*W-24.*V*W**2-12.*V**2*W+12.*V**2*W**2-24.*W+24.*W
     1 **2)
     1 +CF*FQG*(24.-36.*V+12.*V*W-48.*V*X1**(-2)+72.*V*X1**(-1)+12.*V**
     1 2-12.*V**2*W+24.*V**2*X1**(-2)-24.*V**2*X1**(-1)+24.*X1**(-2)-48
     1 .*X1**(-1))+0.
      CD11=CD11+AUCD11
      CD4=CD5+CD6+CD7+CD8+CD9+CD10+CD11
      RETURN
      END
      FUNCTION CE1(VC,WC,E1,E2)
      implicit real*8(a-h,o-z)
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      REAL*8 LV
      W=WC
      V=VC
      LV=LOG(V)
      QI=E1
      QF=E2
       CE1=
     1 +CF*QI**2*(-8./9.+8./9.*V+2./9.*V**2)
     1 +CF*QI**2*LV*(4./3.-4./3.*V+2./3.*V**2)
      RETURN
      END
      FUNCTION CE2(VC,WC,E1,E2)
      implicit real*8(a-h,o-z)
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      W=WC
      V=VC
      QI=E1
      QF=E2
       CE2=
     1 +CF*QI**2*(4./3.-4./3.*V+2./3.*V**2)
      RETURN
      END
      FUNCTION CE4(VC,WC,E1,E2)
      implicit real*8(a-h,o-z)
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      REAL*8 LSQ,LV,LW,LW1
      COMMON/EC/FQA,FQB,FQC,FQE
      QI=E1
      QF=E2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q2=FQA
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LV=LOG(V)
      LW=LOG(W)
      LW1=LOG(1.-W)
      CE11=
     1 +CF*QI**2*(-4./3.+4./3.*V-2./3.*V**2-4./3.*V**2*W**2)
     1 +CF*QF**2*(-4.*V*W+4.*V*W**2+24.*V**(-2)*W-24.*V**(-2)*W**2-48.*
     1 V**(-1)*W+48.*V**(-1)*W**2+28.*W-28.*W**2)
     1 +CF*QF**2*(6.-2.*V+4.*V*W-4.*V*W**2+4.*V**(-2)-8.*V**(-2)*W+
     1 8.*V**(-2)*W**2-8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W**2-12.*W+
     1 12.*W**2)*(LSQ-1.)
     1 +CF*QF**2*FQP*(12.*V*W-12.*V*W**2-24.*V**(-2)*W+24.*V**(-2)*W**2
     1 +48.*V**(-1)*W-48.*V**(-1)*W**2-36.*W+36.*W**2)
C
       CE5=
     1 +CF*QF**2*LV*(6.-2.*V+4.*V*W-4.*V*W**2+4.*V**(-2)-8.*V**(-2)*W+8
     1 .*V**(-2)*W**2-8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W**2-12.*W+1
     1 2.*W**2)
C
       CE9=
     1 +CF*QF**2*LW*FQP*(6.-2.*V+4.*V*W-4.*V*W**2+4.*V**(-2)-8.*V**(-2)
     1 *W+8.*V**(-2)*W**2-8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W**2-12.
     1 *W+12.*W**2)
C
       CE10=
     1 +CF*QF**2*LW1*(6.-2.*V+4.*V*W-4.*V*W**2+4.*V**(-2)-8.*V**(-2)*W+
     1 8.*V**(-2)*W**2-8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W**2-12.*W+
     1 12.*W**2)
     1 +CF*QF**2*LW1*FQP*(-6.+2.*V-4.*V*W+4.*V*W**2-4.*V**(-2)+8.*V**(-
     1 2)*W-8.*V**(-2)*W**2+8.*V**(-1)-16.*V**(-1)*W+16.*V**(-1)*W**2+1
     1 2.*W-12.*W**2)
      CE4=CE5+CE9+CE10+CE11
      RETURN
      END
      FUNCTION CF4(VC,WC,E1,E2)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      EXTERNAL FQ1
      REAL*8 LSQ,LQC,NF,NC,LV,LV1,LX1,LX2,LW,LW1,LR
      COMMON/EC/FQA,FQB,FQC,FQE
      PI2=PI**2
      W=WC
      V=VC
      QI=E1
      QF=E2
      FGQZQ=0.
      DQGY2=0.
      X22=(1.-GV)/(X3*(1.-V))
      Q1=FQ1(1.d0,X22)
      Q2=FQA
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LQC=LOG(SH/Q1)
      LV=LOG(V)
      LV1=LOG(1.-V)
      V2=2.-V
      V1=1.-V
      X1=1.-V*W
      X2=1.-V+V*W
      W1=1.-W
      LX1=LOG(X1)
      LX2=LOG(X2)
      LW=LOG(W)
      LW1=LOG(1.-W)
      LR=LOG(X1/V1)
       CF11=
     1 +CF*QI*QF*(-4.+4.*V+4.*V*W-8.*V*W**2-2.*V*X1**(-1)-4.*V*X2**(-1)
     1 +2.*V**2*X2**(-1)-4.*W+8.*W**2+2.*X1**(-1)+2.*X2**(-1))
     1 +CF*QI**2*(-3.+6.*V-2.*V*W+6.*V*W**2+6.*V*X2**(-2)-15.*V*X2**(-1
     1 )-3.*V**2+4.*V**2*W-6.*V**2*W**2-6.*V**2*X2**(-2)+12.*V**2*X2**(
     1 -1)-4.*V**3*W+4.*V**3*W**2+2.*V**3*X2**(-2)-3.*V**3*X2**(-1)+2.*
     1 W-4.*W**2-2.*X2**(-2)+6.*X2**(-1))
     1 +CF*QI**2*LSQ*(5.-11.*V+6.*V*W-6.*V*W**2-3.*V*X2**(-2)+10.*V*X2*
     1 *(-1)+10.*V**2-12.*V**2*W+12.*V**2*W**2+3.*V**2*X2**(-2)-8.*V**2
     1 *X2**(-1)-4.*V**3+8.*V**3*W-8.*V**3*W**2-V**3*X2**(-2)+2.*V**3*X
     1 2**(-1)-2.*W+2.*W**2+X2**(-2)-4.*X2**(-1))
     1 +CF*QI**2*FQP*(18.*V*W-18.*V*W**2-24.*V**2*W+24.*V**2*W**2+12.*V
     1 **3*W-12.*V**3*W**2-6.*W+6.*W**2)
      CF11=CF11
     1 +CF*QI**2*DQG*(3.*V*X2**(-2)*DQGY2-4.*V*X2**(-1)*DQGY2+2.*V*DQGY
     1 2-3.*V**2*X2**(-2)*DQGY2+2.*V**2*X2**(-1)*DQGY2+V**3*X2**(-2)*DQ
     1 GY2-X2**(-2)*DQGY2+2.*X2**(-1)*DQGY2-2.*DQGY2)
     1 +CF*QF**2*(-3.-6.*V*W+6.*V*W**2+2.*V*X1**(-2)-3.*V*X1**(-1)+4.*V
     1 1**(-1)*W-4.*V1**(-1)*W**2-2.*W-2.*X1**(-2)+6.*X1**(-1))
     1 +CF*QF**2*LSQ*(1.-3.*V+6.*V*W-6.*V*W**2-V*X1**(-2)+2.*V*X1**(-1)
     1 +4.*V1**(-1)-8.*V1**(-1)*W+8.*V1**(-1)*W**2+6.*W-6.*W**2+X1**(-2
     1 )-4.*X1**(-1))
     1 +CF*QF**2*FQP*(6.*V*W-6.*V*W**2-12.*V1**(-1)*W+12.*V1**(-1)*W**2
     1 +6.*W-6.*W**2)
     1 +CF*QF**2*FGQ*(2.*V*W*FGQZQ-X1**(-1)*FGQZQ)
C
       CF5=
     1 +CF*QI**2*LV*(5.-11.*V+6.*V*W-6.*V*W**2-3.*V*X2**(-2)+10.*V*X2**
     1 (-1)+10.*V**2-12.*V**2*W+12.*V**2*W**2+3.*V**2*X2**(-2)-8.*V**2*
     1 X2**(-1)-4.*V**3+8.*V**3*W-8.*V**3*W**2-V**3*X2**(-2)+2.*V**3*X2
     1 **(-1)-2.*W+2.*W**2+X2**(-2)-4.*X2**(-1))
     1 +CF*QF**2*LV*(1.-3.*V+6.*V*W-6.*V*W**2-V*X1**(-2)+2.*V*X1**(-1)+
     1 4.*V1**(-1)-8.*V1**(-1)*W+8.*V1**(-1)*W**2+6.*W-6.*W**2+X1**(-2)
     1 -4.*X1**(-1))
C
       CF8=
     1 +CF*QI*QF*LV1*(-4.*V*W+8.*V*W**2+4.*V**2*W-8.*V**2*W**2)
     1 +CF*QF**2*LV1*(-2.-2.*V+4.*V*W-8.*V*W**2+4.*V1**(-1)-8.*V1**(-1)
     1 *W+8.*V1**(-1)*W**2+8.*W-8.*W**2)
C
       CF9=
     1 +CF*QI*QF*LW*(4.*V*W-6.*V**2*W+12.*V**2*W**2-16.*V**2*W**3-4.*W)
     1 +CF*QI**2*LW*(2.-6.*V+4.*V*W+8.*V**2-12.*V**2*W+8.*V**2*W**2-4.*
     1 V**3+8.*V**3*W-8.*V**3*W**2)
     1 +CF*QI**2*LW*FQP*(1.-3.*V+6.*V*W-6.*V*W**2+4.*V**2-8.*V**2*W+8.*
     1 V**2*W**2-2.*V**3+4.*V**3*W-4.*V**3*W**2-2.*W+2.*W**2)
     1 +CF*QF**2*LW*(-2.-2.*V+4.*V*W-8.*V*W**2+4.*V1**(-1)-8.*V1**(-1)*
     1 W+8.*V1**(-1)*W**2+8.*W-8.*W**2)
     1 +CF*QF**2*LW*FQP*(-1.-V+2.*V*W-2.*V*W**2+2.*V1**(-1)-4.*V1**(-1)
     1 *W+4.*V1**(-1)*W**2+2.*W-2.*W**2)
C
       CF10=
     1 +CF*QI*QF*LW1*(-4.*V*W+10.*V**2*W-20.*V**2*W**2+16.*V**2*W**3+4.
     1 *W)
     1 +CF*QI**2*LW1*(5.-11.*V+6.*V*W-6.*V*W**2-3.*V*X2**(-2)+10.*V*X2*
     1 *(-1)+10.*V**2-12.*V**2*W+12.*V**2*W**2+3.*V**2*X2**(-2)-8.*V**2
     1 *X2**(-1)-4.*V**3+8.*V**3*W-8.*V**3*W**2-V**3*X2**(-2)+2.*V**3*X
     1 2**(-1)-2.*W+2.*W**2+X2**(-2)-4.*X2**(-1))
     1 +CF*QI**2*(-1.+3.*V-6.*V*W+6.*V*W**2-4.*V**2+8.*V**2*W-8
     1 .*V**2*W**2+2.*V**3-4.*V**3*W+4.*V**3*W**2+2.*W-2.*W**2)
     1 +CF*QF**2*LW1*(1.-3.*V+6.*V*W-6.*V*W**2-V*X1**(-2)+2.*V*X1**(-1)
     1 +4.*V1**(-1)-8.*V1**(-1)*W+8.*V1**(-1)*W**2+6.*W-6.*W**2+X1**(-2
     1 )-4.*X1**(-1))
     1 +CF*QF**2*(1.+V-2.*V*W+2.*V*W**2-2.*V1**(-1)+4.*V1**(-1)
     1 *W-4.*V1**(-1)*W**2-2.*W+2.*W**2)
C
       CF6=
     1 +CF*QI*QF*LX1*(-2.*V**2*W+4.*V**2*W**2)
     1 +CF*QF**2*LX1*(2.+2.*V-4.*V*W+8.*V*W**2-4.*V1**(-1)+8.*V1**(-1)*
     1 W-8.*V1**(-1)*W**2-8.*W+8.*W**2)
C
       CF7=
     1 +CF*QI*QF*(-2.*V**2*W*LX2+4.*V**2*W**2*LX2)
     1 +CF*QI**2*(-4.*V*W*LX2-6.*V*X2**(-2)*LX2+20.*V*X2**(-1)*LX2-10.*
     1 V*LX2+4.*V**2*W*LX2+6.*V**2*X2**(-2)*LX2-16.*V**2*X2**(-1)*LX2+4
     1 .*V**2*LX2-2.*V**3*X2**(-2)*LX2+4.*V**3*X2**(-1)*LX2+2.*X2**(-2)
     1 *LX2-8.*X2**(-1)*LX2+6.*LX2)+0.
      CF4=CF5+CF6+CF7+CF8+CF9+CF10+CF11
      LV=LV+2.*LX2
      CFC=DQG*CF*X2**(-1)*(-2.*LV*V**3+8.*LV*V**2-10.*LV*V+
     . 4.*LV-2.*LW1*V**3+8.*LW1*V**2-10.*LW1*V+4.*LW1)+DQG*
     . CF*X2**(-2)*(LV*V**3-3.*LV*V**2+3.*LV*V-LV+LW1*V**3-
     . 3.*LW1*V**2+3.*LW1*V-LW1)+DQG*CF*(4.*LV*W**2*V**3-4.*
     . LV*W**2*V**2-4.*LV*W*V**3+4.*LV*W*V**2+2.*LV*V**3-6.
     . *LV*V**2+8.*LV*V-4.*LV+4.*LW1*W**2*V**3-4.*LW1*W**2*
     . V**2-4.*LW1*W*V**3+4.*LW1*W*V**2+2.*LW1*V**3-6.*LW1*
     . V**2+8.*LW1*V-4.*LW1)
      CF4=CF4+QI**2*CFC
      LV=LV-2.*LX2
      RR3=2.*LV*X1**(-1)*V-4.*LV*X1**(-1)-LV*X1**(-2)*V+LV*
     . X1**(-2)+4.*LV*(-V+1.)**(-1)*W**2-4.*LV*(-V+1.)**(-1
     . )*W+2.*LV*(-V+1.)**(-1)-4.*LV*W**2*V-4.*LV*W**2+4.*
     . LV*W*V+4.*LV*W-2.*LV*V+2.*LV-2.*LX1*X1**(-1)*V+4.*
     . LX1*X1**(-1)+LX1*X1**(-2)*V-LX1*X1**(-2)-4.*LX1*(-V+
     . 1.)**(-1)*W**2+4.*LX1*(-V+1.)**(-1)*W-2.*LX1*(-V+1.)
     . **(-1)+4.*LX1*W**2*V+4.*LX1*W**2-4.*LX1*W*V-4.*LX1*W
     . +2.*LX1*V-2.*LX1+2.*LW1*X1**(-1)*V-4.*LW1*X1**(-1)-
     . LW1*X1**(-2)*V+LW1*X1**(-2)+4.*LW1*(-V+1.)**(-1)*W**
     . 2-4.*LW1*(-V+1.)**(-1)*W+2.*LW1*(-V+1.)**(-1)-4.*LW1
     . *W**2*V-4.*LW1*W**2+4.*LW1*W*V+4.*LW1*W-2.*LW1*V+2.*
     . LW1-68./3.*X1**(-1)*V**2+68.*X1**(-1)*V-136./3.*X1**
     . (-1)+68./3.*X1**(-2)*V**2-170./3.*X1**(-2)*V+34.*X1
     . **(-2)-34./3.*X1**(-3)*V**2+68./3.*X1**(-3)*V-34./3.
     . *X1**(-3)-68./3.*V+68./3.
      CF4=CF4+FGQ*QF**2*RR3*(-CF)
C
      LV=1.
      RR4=2.*LV*X1**(-1)*V-4.*LV*X1**(-1)-LV*X1**(-2)*V+LV*
     . X1**(-2)+4.*LV*(-V+1.)**(-1)*W**2-4.*LV*(-V+1.)**(-1
     . )*W+2.*LV*(-V+1.)**(-1)-4.*LV*W**2*V-4.*LV*W**2+4.*
     . LV*W*V+4.*LV*W-2.*LV*V+2.*LV
      CF4=CF4+QF**2*CF*RR4
      RETURN
      END
      FUNCTION CJ1(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      REAL*8 LV
      W=WC
      V=VC
      LV=LOG(V)
       CJ1=
     1 +CF*(-8./9.+8./9.*V+2./9.*V**2)
     1 +CF*LV*(4./3.-4./3.*V+2./3.*V**2)
      RETURN
      END
      FUNCTION CJ2(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      W=WC
      V=VC
       CJ2=+CF*(4./3.-4./3.*V+2./3.*V**2)
      RETURN
      END

      FUNCTION CJ4(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      EXTERNAL FQ1
      REAL*8 LSQ,LQC,NF,NC,LV,LV1,LX1,LX2,LW,LW1,LR
      COMMON/EC/FQA,FQB,FQC,FQE
      PI2=PI**2
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q1=FQ1(1.d0,X22)
      Q2=FQA
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LQC=LOG(SH/Q1)
      LV=LOG(V)
      LV1=LOG(1.-V)
      V2=2.-V
      V1=1.-V
      X1=1.-V*W
      X2=1.-V+V*W
      W1=1.-W
      LX1=LOG(X1)
      LX2=LOG(X2)
      LW=LOG(W)
      LW1=LOG(1.-W)
      LR=LOG(X1/V1)
      CJ11=
     1 +CF*(-34./3.+34./3.*V-8.*V*W+36.*V*W*FQP+8.*V*W**2-36.*V*W**2*FQ
     1 P+2.*V*X1**(-2)-5.*V*X1**(-1)+6.*V*X2**(-2)-19.*V*X2**(-1)-11./3
     1 .*V**2+4.*V**2*W-24.*V**2*W*FQP-22./3.*V**2*W**2+24.*V**2*W**2*F
     1 QP-6.*V**2*X2**(-2)+14.*V**2*X2**(-1)-4.*V**3*W+12.*V**3*W*FQP+4
     1 .*V**3*W**2-12.*V**3*W**2*FQP+2.*V**3*X2**(-2)-3.*V**3*X2**(-1)+
     1 24.*V**(-2)*W-24.*V**(-2)*W*FQP-24.*V**(-2)*W**2+24.*V**(-2)*W**
     1 2*FQP-48.*V**(-1)*W+48.*V**(-1)*W*FQP+48.*V**(-1)*W**2-48.*V**(-
     1 1)*W**2*FQP+4.*V1**(-1)*W-12.*V1**(-1)*W*FQP-4.*V1**(-1)*W**2+12
     1 .*V1**(-1)*W**2*FQP+24.*W-36.*W*FQP-24.*W**2+36.*W**2*FQP-2.*X1*
     1 *(-2)+8.*X1**(-1)-2.*X2**(-2)+8.*X2**(-1))
      CJ11=CJ11
     1 +CF*NC**(-1)*(4.*V*W-36.*V*W*FQP-8.*V*W**2+36.*V*W**2*FQP-4.*V*W
     1 1**(-1)*LR-V*X1**(-2)+V*X1**(-1)-3.*V*X2**(-2)+11.*V*X2**(-1)+V*
     1 *2+12.*V**2*W*FQP+6.*V**2*W**2-12.*V**2*W**2*FQP+4.*V**2*W1**(-1
     1 )*LR+3.*V**2*X2**(-2)-10.*V**2*X2**(-1)-V**3*X2**(-2)+3.*V**3*X2
     1 **(-1)-4.*W+36.*W*FQP+8.*W**2-36.*W**2*FQP+2.*W1**(-1)*LR+X1**(-
     1 2)-4.*X1**(-1)+X2**(-2)-4.*X2**(-1))
      CJ11=CJ11
     1 +CF*NC**(-1)*(-6.+6.*V-12.*V*W+12.*V*W**2-2.*V**2+4.*V**2*W-
     1 4.*V**2*W**2+12.*W-12.*W**2)*(LSQ-1.)
     1 +CF*LSQ*(12.-16.*V+16.*V*W-16.*V*W**2-V*X1**(-2)+2.*V*X1**(-1)-3
     1 .*V*X2**(-2)+10.*V*X2**(-1)+10.*V**2-12.*V**2*W+12.*V**2*W**2+3.
     1 *V**2*X2**(-2)-8.*V**2*X2**(-1)-4.*V**3+8.*V**3*W-8.*V**3*W**2-V
     1 **3*X2**(-2)+2.*V**3*X2**(-1)+4.*V**(-2)-8.*V**(-2)*W+8.*V**(-2)
     1 *W**2-8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W**2+4.*V1**(-1)-8.*V
     1 1**(-1)*W+8.*V1**(-1)*W**2-8.*W+8.*W**2+X1**(-2)-4.*X1**(-1)+X2*
     1 *(-2)-4.*X2**(-1))
C
       CJ5=
     1 +CF*NC**(-1)*LV*(-6.+6.*V-12.*V*W+12.*V*W**2-2.*V**2+4.*V**2*W-4
     1 .*V**2*W**2+12.*W-12.*W**2)
     1 +CF*LV*(12.-16.*V+16.*V*W-16.*V*W**2-V*X1**(-2)+2.*V*X1**(-1)-3.
     1 *V*X2**(-2)+10.*V*X2**(-1)+10.*V**2-12.*V**2*W+12.*V**2*W**2+3.*
     1 V**2*X2**(-2)-8.*V**2*X2**(-1)-4.*V**3+8.*V**3*W-8.*V**3*W**2-V*
     1 *3*X2**(-2)+2.*V**3*X2**(-1)+4.*V**(-2)-8.*V**(-2)*W+8.*V**(-2)*
     1 W**2-8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W**2+4.*V1**(-1)-8.*V1
     1 **(-1)*W+8.*V1**(-1)*W**2-8.*W+8.*W**2+X1**(-2)-4.*X1**(-1)+X2**
     1 (-2)-4.*X2**(-1))
C
       CJ8=
     1 +CF*NC**(-1)*LV1*(8.-12.*V+8.*V*W-8.*V*W**2+8.*V**2-4.*V**2*W+8.
     1 *V**2*W**2-4.*V**(-1))
     1 +CF*LV1*(-2.-2.*V+4.*V**2*W-8.*V**2*W**2+4.*V1**(-1)-8.*V1**(-1)
     1 *W+8.*V1**(-1)*W**2+8.*W-8.*W**2)
C
       CJ9=
     1 +CF*NC**(-1)*LW*(-2.+2.*V-8.*V*W-12.*V*W*FQP+12.*V*W**2*FQP+4.*V
     1 *W1**(-1)+6.*V*FQP+4.*V**2+8.*V**2*W+4.*V**2*W*FQP+4.*V**2*W**2-
     1 4.*V**2*W**2*FQP-6.*V**2*W1**(-1)-2.*V**2*FQP+8.*W+12.*W*FQP-12.
     1 *W**2*FQP-4.*W1**(-1)-6.*FQP)
     1 +CF*LW*(-8.*V+12.*V*W+12.*V*W*FQP-8.*V*W**2-12.*V*W**2*FQP-6.*V*
     1 FQP+8.*V**2-18.*V**2*W-8.*V**2*W*FQP+20.*V**2*W**2+8.*V**2*W**2*
     1 FQP-16.*V**2*W**3+4.*V**2*FQP-4.*V**3+8.*V**3*W+4.*V**3*W*FQP-8.
     1 *V**3*W**2-4.*V**3*W**2*FQP-2.*V**3*FQP-8.*V**(-2)*W*FQP+8.*V**(
     1 -2)*W**2*FQP+4.*V**(-2)*FQP+16.*V**(-1)*W*FQP-16.*V**(-1)*W**2*F
     1 QP-8.*V**(-1)*FQP+4.*V1**(-1)-8.*V1**(-1)*W-4.*V1**(-1)*W*FQP+8.
     1 *V1**(-1)*W**2+4.*V1**(-1)*W**2*FQP+2.*V1**(-1)*FQP+4.*W-12.*W*F
     1 QP-8.*W**2+12.*W**2*FQP+6.*FQP)
C
       CJ10=
     1 +CF*NC**(-1)*LW1*(-4.*V*W+12.*V*W*FQP+12.*V*W**2-12.*V*W**2*FQP-
     1 6.*V*FQP+2.*V**2-4.*V**2*W-4.*V**2*W*FQP+4.*V**2*W**2*FQP+2.*V**
     1 2*FQP+4.*W-12.*W*FQP-12.*W**2+12.*W**2*FQP+6.*FQP)
      CJ10=CJ10
     1 +CF*LW1*(12.-16.*V+12.*V*W-12.*V*W*FQP-16.*V*W**2+12.*V*W**2*FQP
     1 -V*X1**(-2)+2.*V*X1**(-1)-3.*V*X2**(-2)+10.*V*X2**(-1)+6.*V*FQP+
     1 10.*V**2-2.*V**2*W+8.*V**2*W*FQP-8.*V**2*W**2-8.*V**2*W**2*FQP+1
     1 6.*V**2*W**3+3.*V**2*X2**(-2)-8.*V**2*X2**(-1)-4.*V**2*FQP-4.*V*
     1 *3+8.*V**3*W-4.*V**3*W*FQP-8.*V**3*W**2+4.*V**3*W**2*FQP-V**3*X2
     1 **(-2)+2.*V**3*X2**(-1)+2.*V**3*FQP+4.*V**(-2)-8.*V**(-2)*W+8.*V
     1 **(-2)*W*FQP+8.*V**(-2)*W**2-8.*V**(-2)*W**2*FQP-4.*V**(-2)*FQP-
     1 8.*V**(-1)+16.*V**(-1)*W-16.*V**(-1)*W*FQP-16.*V**(-1)*W**2+16.*
     1 V**(-1)*W**2*FQP+8.*V**(-1)*FQP+4.*V1**(-1)-8.*V1**(-1)*W+4.*V1*
     1 *(-1)*W*FQP+8.*V1**(-1)*W**2-4.*V1**(-1)*W**2*FQP-2.*V1**(-1)*FQ
     1 P-4.*W+12.*W*FQP+8.*W**2-12.*W**2*FQP+X1**(-2)-4.*X1**(-1)+X2**(
     1 -2)-4.*X2**(-1)-6.*FQP)
C
       CJ10=CJ10
     1 +CF*(-12.*V*W+12.*V*W**2+6.*V
     1 +8.*V**2*W-8.*V**2*W**2
     1 -4.*V**2
     1 -4.*V**3*W+4.*V**3*W**2
     1 +2.*V**3+8.*V
     1 **(-2)*W-8.*V**(-2)*W**2-4.*V**(-2)
     1 -16.*V**(-1)*W+16.*
     1 V**(-1)*W**2+8.*V**(-1)+4.*V1*
     1 *(-1)*W-4.*V1**(-1)*W**2-2.*V1**(-1)
     1 +12.*W-12.*W**2
     1 -6.)
C
       CJ6=
     1 +CF*NC**(-1)*LX1*(-2.+2.*V+8.*V*W-4.*V**2-4.*V**2*W-4.*V**2*W**2
     1 )
     1 +CF*LX1*(2.+2.*V-4.*V*W+8.*V*W**2-2.*V**2*W+4.*V**2*W**2-4.*V1**
     1 (-1)+8.*V1**(-1)*W-8.*V1**(-1)*W**2-8.*W+8.*W**2)
C
       CJ7=
     1 +CF*NC**(-1)*LX2*(-2.+2.*V-8.*V*W-4.*V**2+4.*V**2*W-4.*V**2*W**2
     1 +2.*V**2*W1**(-1)+2.*W1**(-1))
     1 +CF*LX2*(6.-10.*V-4.*V*W-6.*V*X2**(-2)+20.*V*X2**(-1)+4.*V**2+2.
     1 *V**2*W+4.*V**2*W**2+6.*V**2*X2**(-2)-16.*V**2*X2**(-1)-2.*V**3*
     1 X2**(-2)+4.*V**3*X2**(-1)+2.*X2**(-2)-8.*X2**(-1))+0.
      CJ4=CJ5+CJ6+CJ7+CJ8+CJ9+CJ10+CJ11
      LV=LV+2.*LX2
      CFC=DQG*CF*X2**(-1)*(-2.*LV*V**3+8.*LV*V**2-10.*LV*V+
     . 4.*LV-2.*LW1*V**3+8.*LW1*V**2-10.*LW1*V+4.*LW1)+DQG*
     . CF*X2**(-2)*(LV*V**3-3.*LV*V**2+3.*LV*V-LV+LW1*V**3-
     . 3.*LW1*V**2+3.*LW1*V-LW1)+DQG*CF*(4.*LV*W**2*V**3-4.*
     . LV*W**2*V**2-4.*LV*W*V**3+4.*LV*W*V**2+2.*LV*V**3-6.
     . *LV*V**2+8.*LV*V-4.*LV+4.*LW1*W**2*V**3-4.*LW1*W**2*
     . V**2-4.*LW1*W*V**3+4.*LW1*W*V**2+2.*LW1*V**3-6.*LW1*
     . V**2+8.*LW1*V-4.*LW1)
      CJ4=CJ4+CFC
      LV=LV-2.*LX2
      RR3=2.*LV*X1**(-1)*V-4.*LV*X1**(-1)-LV*X1**(-2)*V+LV*
     . X1**(-2)+4.*LV*(-V+1.)**(-1)*W**2-4.*LV*(-V+1.)**(-1
     . )*W+2.*LV*(-V+1.)**(-1)-4.*LV*W**2*V-4.*LV*W**2+4.*
     . LV*W*V+4.*LV*W-2.*LV*V+2.*LV-2.*LX1*X1**(-1)*V+4.*
     . LX1*X1**(-1)+LX1*X1**(-2)*V-LX1*X1**(-2)-4.*LX1*(-V+
     . 1.)**(-1)*W**2+4.*LX1*(-V+1.)**(-1)*W-2.*LX1*(-V+1.)
     . **(-1)+4.*LX1*W**2*V+4.*LX1*W**2-4.*LX1*W*V-4.*LX1*W
     . +2.*LX1*V-2.*LX1+2.*LW1*X1**(-1)*V-4.*LW1*X1**(-1)-
     . LW1*X1**(-2)*V+LW1*X1**(-2)+4.*LW1*(-V+1.)**(-1)*W**
     . 2-4.*LW1*(-V+1.)**(-1)*W+2.*LW1*(-V+1.)**(-1)-4.*LW1
     . *W**2*V-4.*LW1*W**2+4.*LW1*W*V+4.*LW1*W-2.*LW1*V+2.*
     . LW1-68./3.*X1**(-1)*V**2+68.*X1**(-1)*V-136./3.*X1**
     . (-1)+68./3.*X1**(-2)*V**2-170./3.*X1**(-2)*V+34.*X1
     . **(-2)-34./3.*X1**(-3)*V**2+68./3.*X1**(-3)*V-34./3.
     . *X1**(-3)-68./3.*V+68./3.
      CJ4=CJ4+FGQ*RR3*(-CF)
C
      LV=1.
      RR4=2.*LV*X1**(-1)*V-4.*LV*X1**(-1)-LV*X1**(-2)*V+LV*
     . X1**(-2)+4.*LV*(-V+1.)**(-1)*W**2-4.*LV*(-V+1.)**(-1
     . )*W+2.*LV*(-V+1.)**(-1)-4.*LV*W**2*V-4.*LV*W**2+4.*
     . LV*W*V+4.*LV*W-2.*LV*V+2.*LV
      CJ4=CJ4+CF*RR4
      RETURN
      END
      FUNCTION CK4(VC,WC)
      implicit real*8(a-h,o-z)
      real*8 is,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/VAR/S,SS,PT
      EXTERNAL FQ1
      REAL*8 NC,LSQ,LV,LV1,LW,LW1,LX1,LX2
      COMMON/EC/FQA,FQB,FQC,FQE
      W=WC
      V=VC
      X22=(1.-GV)/(X3*(1.-V))
      Q2=FQA
      SH=X22*S
      LSQ=LOG(SH/Q2)
      LV=LOG(V)
      LV1=LOG(1.-V)
      V1=1.-V
      X1=1.-V*W
      X2=1.-V+V*W
      W1=1.-W
      LX1=LOG(X1)
      LX2=LOG(X2)
      LW=LOG(W)
      LW1=LOG(1.-W)
      CK11=
     1 +CF*(-2.+2.*V-12.*V*W+24.*V*W*FQP+20.*V*W**2-24.*V*W**2*FQP+2.*V
     1 *X1**(-2)-V*X1**(-1)+6.*V*X2**(-2)-11.*V*X2**(-1)-3.*V**2+4.*V**
     1 2*W-24.*V**2*W*FQP-6.*V**2*W**2+24.*V**2*W**2*FQP-6.*V**2*X2**(-
     1 2)+10.*V**2*X2**(-1)-4.*V**3*W+12.*V**3*W*FQP+4.*V**3*W**2-12.*V
     1 **3*W**2*FQP+2.*V**3*X2**(-2)-3.*V**3*X2**(-1)+4.*V1**(-1)*W-12.
     1 *V1**(-1)*W*FQP-4.*V1**(-1)*W**2+12.*V1**(-1)*W**2*FQP+4.*W-12.*
     1 W**2-2.*X1**(-2)+4.*X1**(-1)-2.*X2**(-2)+4.*X2**(-1))
     1 +CF*NC**(-1)*(2.*V**2-12.*V**2*W*FQP+12.*V**2*W**2*FQP)
     1 +CF*NC**(-1)*(LSQ-1.)*(2.*V**2-4.*V**2*W+4.*V**2*W**2)
     1 +CF*LSQ*(6.-14.*V+12.*V*W-12.*V*W**2-V*X1**(-2)+2.*V*X1**(-1)-3.
     1 *V*X2**(-2)+10.*V*X2**(-1)+10.*V**2-12.*V**2*W+12.*V**2*W**2+3.*
     1 V**2*X2**(-2)-8.*V**2*X2**(-1)-4.*V**3+8.*V**3*W-8.*V**3*W**2-V*
     1 *3*X2**(-2)+2.*V**3*X2**(-1)+4.*V1**(-1)-8.*V1**(-1)*W+8.*V1**(-
     1 1)*W**2+4.*W-4.*W**2+X1**(-2)-4.*X1**(-1)+X2**(-2)-4.*X2**(-1))
C
       CK5=
     1 +CF*NC**(-1)*LV*(-2.*V**2+4.*V**2*W-4.*V**2*W**2)
     1 +CF*LV*(6.-14.*V+12.*V*W-12.*V*W**2-V*X1**(-2)+2.*V*X1**(-1)-3.*
     1 V*X2**(-2)+10.*V*X2**(-1)+10.*V**2-12.*V**2*W+12.*V**2*W**2+3.*V
     1 **2*X2**(-2)-8.*V**2*X2**(-1)-4.*V**3+8.*V**3*W-8.*V**3*W**2-V**
     1 3*X2**(-2)+2.*V**3*X2**(-1)+4.*V1**(-1)-8.*V1**(-1)*W+8.*V1**(-1
     1 )*W**2+4.*W-4.*W**2+X1**(-2)-4.*X1**(-1)+X2**(-2)-4.*X2**(-1))
C
       CK8=
     1 +CF*NC**(-1)*LV1*(12.-8.*V-8.*V*W**2+2.*V*X1**(-1)+14.*V*X2**(-1
     1 )+4.*V**2-4.*V**2*W+8.*V**2*W**2-8.*V**2*X2**(-1)+2.*V**3*X2**(-
     1 1)-8.*V**(-1)+4.*V**(-1)*X1**(-1)+4.*V**(-1)*X2**(-1)-4.*X1**(-1
     1 )-12.*X2**(-1))
     1 +CF*LV1*(-2.-2.*V+8.*V*W-16.*V*W**2-4.*V**2*W+8.*V**2*W**2+4.*V1
     1 **(-1)-8.*V1**(-1)*W+8.*V1**(-1)*W**2+8.*W-8.*W**2)
C
       CK9=
     1 +CF*NC**(-1)*LW*(-4.+4.*V+8.*V*W+2.*V*X1**(-1)-14.*V*X2**(-1)-4.
     1 *V**2*W-4.*V**2*W*FQP+4.*V**2*W**2*FQP+8.*V**2*X2**(-1)+2.*V**2*
     1 FQP-2.*V**3*X2**(-1)+4.*V**(-1)*X1**(-1)-4.*V**(-1)*X2**(-1)-8.*
     1 W-4.*X1**(-1)+12.*X2**(-1))
     1 +CF*LW*(-8.*V+4.*V*W+8.*V*W*FQP-8.*V*W**2-8.*V*W**2*FQP-4.*V*FQP
     1 +8.*V**2-6.*V**2*W-8.*V**2*W*FQP-4.*V**2*W**2+8.*V**2*W**2*FQP+1
     1 6.*V**2*W**3+4.*V**2*FQP-4.*V**3+8.*V**3*W+4.*V**3*W*FQP-8.*V**3
     1 *W**2-4.*V**3*W**2*FQP-2.*V**3*FQP+4.*V1**(-1)-8.*V1**(-1)*W-4.*
     1 V1**(-1)*W*FQP+8.*V1**(-1)*W**2+4.*V1**(-1)*W**2*FQP+2.*V1**(-1)
     1 *FQP+12.*W-8.*W**2)
C
       CK10=
     1 +CF*NC**(-1)*LW1*(4.-4.*V-8.*V*W-2.*V*X1**(-1)+14.*V*X2**(-1)+2.
     1 *V**2+4.*V**2*W*FQP+4.*V**2*W**2-4.*V**2*W**2*FQP-8.*V**2*X2**(-
     1 1)-2.*V**2*FQP+2.*V**3*X2**(-1)-4.*V**(-1)*X1**(-1)+4.*V**(-1)*X
     1 2**(-1)+8.*W+4.*X1**(-1)-12.*X2**(-1))
     1 +CF*LW1*(6.-14.*V+16.*V*W-8.*V*W*FQP-12.*V*W**2+8.*V*W**2*FQP-V*
     1 X1**(-2)+2.*V*X1**(-1)-3.*V*X2**(-2)+10.*V*X2**(-1)+4.*V*FQP+10.
     1 *V**2-22.*V**2*W+8.*V**2*W*FQP+32.*V**2*W**2-8.*V**2*W**2*FQP-16
     1 .*V**2*W**3+3.*V**2*X2**(-2)-8.*V**2*X2**(-1)-4.*V**2*FQP-4.*V**
     1 3+8.*V**3*W-4.*V**3*W*FQP-8.*V**3*W**2+4.*V**3*W**2*FQP-V**3*X2*
     1 *(-2)+2.*V**3*X2**(-1)+2.*V**3*FQP+4.*V1**(-1)-8.*V1**(-1)*W+4.*
     1 V1**(-1)*W*FQP+8.*V1**(-1)*W**2-4.*V1**(-1)*W**2*FQP-2.*V1**(-1)
     1 *FQP-4.*W**2+X1**(-2)-4.*X1**(-1)+X2**(-2)-4.*X2**(-1))
       CK10=CK10
     1 +CF*V*(-8.*W+8.*W**2+4.+8.*V*W-8.*V*W**2-4.*V-4.*V**2*W)
     1 +CF*V*(4.*V**2*W**2+2.*V**2)+CF*V1**(-1)*(4.*W-4.*W**2-2.)
C
       CK6=
     1 +CF*LX1*(2.+2.*V-4.*V*W+8.*V*W**2+2.*V**2*W-4.*V**2*W**2-4.*V1**
     1 (-1)+8.*V1**(-1)*W-8.*V1**(-1)*W**2-8.*W+8.*W**2)
C
       CK7=
     1 +CF*(-4.*V*W*LX2-6.*V*X2**(-2)*LX2+20.*V*X2**(-1)*LX2-10.*V*LX2+
     1 6.*V**2*W*LX2-4.*V**2*W**2*LX2+6.*V**2*X2**(-2)*LX2-16.*V**2*X2*
     1 *(-1)*LX2+4.*V**2*LX2-2.*V**3*X2**(-2)*LX2+4.*V**3*X2**(-1)*LX2+
     1 2.*X2**(-2)*LX2-8.*X2**(-1)*LX2+6.*LX2)+0.
      CK4=CK5+CK6+CK7+CK8+CK9+CK10+CK11
      LV=LV+2.*LX2
      CFC=DQG*CF*X2**(-1)*(-2.*LV*V**3+8.*LV*V**2-10.*LV*V+
     . 4.*LV-2.*LW1*V**3+8.*LW1*V**2-10.*LW1*V+4.*LW1)+DQG*
     . CF*X2**(-2)*(LV*V**3-3.*LV*V**2+3.*LV*V-LV+LW1*V**3-
     . 3.*LW1*V**2+3.*LW1*V-LW1)+DQG*CF*(4.*LV*W**2*V**3-4.*
     . LV*W**2*V**2-4.*LV*W*V**3+4.*LV*W*V**2+2.*LV*V**3-6.
     . *LV*V**2+8.*LV*V-4.*LV+4.*LW1*W**2*V**3-4.*LW1*W**2*
     . V**2-4.*LW1*W*V**3+4.*LW1*W*V**2+2.*LW1*V**3-6.*LW1*
     . V**2+8.*LW1*V-4.*LW1)
      CK4=CK4+CFC
      LV=LV-2.*LX2
      RR3=2.*LV*X1**(-1)*V-4.*LV*X1**(-1)-LV*X1**(-2)*V+LV*
     . X1**(-2)+4.*LV*(-V+1.)**(-1)*W**2-4.*LV*(-V+1.)**(-1
     . )*W+2.*LV*(-V+1.)**(-1)-4.*LV*W**2*V-4.*LV*W**2+4.*
     . LV*W*V+4.*LV*W-2.*LV*V+2.*LV-2.*LX1*X1**(-1)*V+4.*
     . LX1*X1**(-1)+LX1*X1**(-2)*V-LX1*X1**(-2)-4.*LX1*(-V+
     . 1.)**(-1)*W**2+4.*LX1*(-V+1.)**(-1)*W-2.*LX1*(-V+1.)
     . **(-1)+4.*LX1*W**2*V+4.*LX1*W**2-4.*LX1*W*V-4.*LX1*W
     . +2.*LX1*V-2.*LX1+2.*LW1*X1**(-1)*V-4.*LW1*X1**(-1)-
     . LW1*X1**(-2)*V+LW1*X1**(-2)+4.*LW1*(-V+1.)**(-1)*W**
     . 2-4.*LW1*(-V+1.)**(-1)*W+2.*LW1*(-V+1.)**(-1)-4.*LW1
     . *W**2*V-4.*LW1*W**2+4.*LW1*W*V+4.*LW1*W-2.*LW1*V+2.*
     . LW1-68./3.*X1**(-1)*V**2+68.*X1**(-1)*V-136./3.*X1**
     . (-1)+68./3.*X1**(-2)*V**2-170./3.*X1**(-2)*V+34.*X1
     . **(-2)-34./3.*X1**(-3)*V**2+68./3.*X1**(-3)*V-34./3.
     . *X1**(-3)-68./3.*V+68./3.
      CK4=CK4+FGQ*RR3*(-CF)
C
      LV=1.
      RR4=2.*LV*X1**(-1)*V-4.*LV*X1**(-1)-LV*X1**(-2)*V+LV*
     . X1**(-2)+4.*LV*(-V+1.)**(-1)*W**2-4.*LV*(-V+1.)**(-1
     . )*W+2.*LV*(-V+1.)**(-1)-4.*LV*W**2*V-4.*LV*W**2+4.*
     . LV*W*V+4.*LV*W-2.*LV*V+2.*LV
      CK4=CK4+CF*RR4
      RETURN
      END

      SUBROUTINE DISTGAM(GPU,GPD,GPC,GPG,X,QQ2)
      implicit real*8(a-h,o-z)
C*****************************************************************************
C
C SUBROUTINE EMPRUNTEE A DIST.F CONTENU DANS HOME/CDC/MADRID
C
C*****************************************************************************
      DIMENSION Q(7),QQ(7)
      character*20 parm(20)
      real*8 h1pdf(20),h2pdf(20)
      common/mypdfsets/h1pdf,h2pdf
      real*8 UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL
      real * 4 delta,q2,q02
      COMMON/ANEW/DELTA
      COMMON/CONV/ IORD,ICONV,OWLAM,OWLAM2,RLAM,RLAM2
      COMMON/QS/Q2
      COMMON/Q000/Q02
      COMMON/GPHO/PAR1(30),PAR2(30),CAPOL(8,20,32),CAVDM(8,20,32)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      REAL*8 KA
c      EXTERNAL DIST,DIST2


      parm(1) = 'Nptype'
      parm(2) = 'Ngroup'
      parm(3) = 'Nset'

      sqrq2 = sqrt(qq2)

      call fonllpdfset(parm,h2pdf)

       call fonllstructm(dble(x),dble(sqrq2),upv,dnv,usea,dsea,str,
     #             chm,bot,top,gl)
       uplus = upv + usea
       dplus = dnv + dsea
       cplus = 2*chm

c......Questo per il set (3 6 3), cioe' AFG
c       uplus = upv 
c       dplus = dnv 
c       cplus = chm

       glu = gl

      goto 100
cccccccccccccccccccccccccccccccccccccc

c      KA=1.0
c      IORD=INT(PAR1(28))
c      FLAV=PAR1(25)
C
C      LES PARAMETRES SONT FIXES A PARTIR DU FICHIER POINTLIKE GRPOL
C     
c      DELTA=PAR1(29)
c      OWLAM=PAR1(1)
c      OWLAM2=OWLAM**2
c      Q02=PAR1(30)
c      FLAVOR=FLAV
c      Q2=QQ2

c      CALL DIST(X,Q)
c      ADD=Q(1)/FLAVOR
c      UPLUS=Q(5)+ADD
c      DPLUS=-Q(4)+ADD
c      SPLUS=-Q(6)+ADD
c      CPLUS=-Q(3)+ADD
c      SING=Q(1)
c      GLU=Q(7)         
c      CALL DIST2(X,QQ)
c      ADD2=QQ(1)/FLAVOR
c      UPLU2=QQ(5)+ADD2
c      DPLU2=-QQ(4)+ADD2
c      SPLU2=-QQ(6)+ADD2
c      CPLU2=-QQ(3)+ADD2
c      SING2=QQ(1)
c      GLU2=QQ(7)
c      UPLUS=UPLUS+UPLU2*KA
c      DPLUS=DPLUS+DPLU2*KA
c      SPLUS=SPLUS+SPLU2*KA
c      CPLUS=CPLUS+CPLU2*KA
c      SING=SING+SING2*KA
c      GLU=GLU+GLU2*KA

ccccccccccccccccccccccccccccccccccccccccccccc
 100   continue
c.....in termini degli output di structm in PDFLIB, vale
c.....uplus,dplus,cplus,glu = 2upv,2dnv,2chm,gl
      ALEMI=137.
      GPU=ALEMI*UPLUS/X/2.
c      gpu = 0.
      GPD=ALEMI*DPLUS/X/2.
c      gpd = 0.
      GPC=ALEMI*CPLUS/X/2.
c      gpc = 0.
      GPG=ALEMI*GLU/X
c      gpg = 0.
      RETURN
      END
      SUBROUTINE DISTPRO(GU,GD,GS,GC,GG,X,Q1S)
      implicit real*8(a-h,o-z)
      DIMENSION Q(7)
      integer ifbm
      character*20 parm(20)
      real*8 h1pdf(20),h2pdf(20)
      COMMON/HADR/IH
      common/mypdfsets/h1pdf,h2pdf
      real*8 UPV,DNV,SEA,STR,CHM,BOT,TOP,GL,usea,dsea,phfbm,assc,xmb,xmc
      COMMON/CONPR/ IORD,ICONV,OWLAM,OWLAM2,RLAM,RLAM2,DELTA
      COMMON/QSPR/Q2
      COMMON/GPRO/PAR(30),CAPR(8,20,32)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      common/w50510/iflprt
      common/fakepdf/ifbm
      common/asscale/assc
      common/hvqmass/xmb,xmc
      integer iret
      real*8 pdfeff,guo,gdo,gso,gco,ggo,dummy
c Check if pdf densities have been already evaluated
      pdfeff=1.d5*h1pdf(1)+1.d4*h1pdf(2)+h1pdf(3)
      call cachelookup(iret,x,q1s,pdfeff,guo,gdo,gso,gco,ggo,dummy)
      if(iret.eq.0)then
        gu=guo
        gd=gdo
        gs=gso
        gc=gco
        gg=ggo
        return
      endif
c
c If not, compute them and store their values
      Q2=Q1S
C  Q2 DEPENDENCE TURNED ON
      IORD=1
      DELTA=.08
      OWLAM=PAR(1)
      OWLAM2=OWLAM**2

      parm(1) = 'Nptype'
      parm(2) = 'Ngroup'
      parm(3) = 'Nset'


      sqrq2 = sqrt(q2)
      call fonllpdfset(parm,h1pdf)
      iflprt = 0
c The following is obsolete (original coding)            
c      CALL STRUCTF(dble(X),dble(sqrq2),UPV,DNV,SEA,STR,CHM,BOT,TOP,GL)
c new coding
      CALL FONLLSTRUCTM
     #  (dble(X),dble(sqrq2),UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
c end new coding
c This below is ridiculous... I take it away by force...P.N.
cc....if hadron is resolved photon there is no valence contribution
c      if(ih.eq.2) then
c         upv = 0.
c         dnv = 0.
c      endif

c old coding: see structf call earlier
c      gu = upv/x
c      gd = dnv/x
c      gs = sea/x
c new coding: this is a hack to make the program work
c with Q+Qbar production, in the general case with str#sea,
c and usea#dsea. It uses the fact that quark and antiquark
c cross sections are the same, and one can put all charge 2/3
c quark content in upv, and all charge 1/3 quark content in dv
      gu=(upv+2*usea)/x
      gd=(dnv+2*dsea+2*str)/x
      gs=0
c end new coding
      if(ifbm.eq.0) then
        gc = chm/x/assc
      elseif(ifbm.eq.1) then
        gc = phfbm(dble(x),dble(sqrq2),xmc)/x
      endif
      	
      gg = gl/x 
c       gg = 1./x

c      CALL DISTP(X,Q)
c      GU=(Q(1)-Q(3))/X
c      GD=(Q(2)-Q(3))/X
c      GS=Q(3)/X
c      GC=Q(5)/X
cc      GG=Q(7)/X 
c      gu = 0.
c      gd = 0.
c      gs = 0.
c      gc = 0.
c      gg = 0.
c Save the pdf densities 
c TEST ************ SEE IF WE CAN DO THIS:
c set to zero the see, set the valence equal to up+upbar, d+dbar,
c add s+sbar to d+dbar (they have the same electric charge)
c      gu=gu+2*gs
c      gd=gd+4*gs
c      gs=0
c THE TEST HAS WORKED OUT RIGHT!
      call cachestore(iret,x,q1s,pdfeff,gu,gd,gs,gc,gg,dummy)

77    RETURN
      END


      FUNCTION ECHE(YV)
      implicit real*8(a-h,o-z)
      COMMON/CS/CS1,CS2,CS3
      COMMON/VAR/S,SS,PT
      COMMON/CUSC/rusc,USC,PUSC
      COMMON/EC/FQA,FQB,FQC,FQE
      common/quarkmass/qm,xmq,ims
      
      PT2=PT**2
      
      if(ims.eq.1) then
          qm2 = xmq**2		! massive scale (transverse mass)
      else 
          qm2 = 0.
      endif
      
   1  FQA=USC*(PT2 + qm2)
      FQB=USC*(PT2 + qm2)
      FQC=USC*(PT2 + qm2)
      FQE=USC*(PT2 + qm2)
      IF(PUSC.GT.1.5) GO TO 2
      FQA=FQA*CS2/USC
      FQB=FQB*CS2/USC
      FQC=FQC*CS2/USC
    2  ECHE=0.
      RETURN
      END
      
      
c Phfonfra returns in xdcp4 the sum (divided by two) of charm and anticharm.
c Given the fact that in the rest of the code the short distance cross 
c section for the production of c and cbar are always summed, this amounts
c to computing the average cross section of heavy-flavoured hadron and
c heavy-flavoured antihadron. In this way, the only asymmetry which
c cancels in the sum is that due to the short-distance cross sections
c for the production of heavy flavours; the fragmentation functions
c are symmetric in the light-flavour sector at the NLO level
      SUBROUTINE PHFONFRA(XC4,IPI,QP24,XDUP4,XDUBP4,XDDP4,XDDBP4,XDSP4,
     # XDCP4,XDBP4,XDBBP4,XDTP4,XDTBP4,XDGP4)
c xc4 is called x3 in the rest of the code: is the fragmentation variable;
c however, it can be equal to one when computing the subtraction term 
c (see iflagmode). Thus, use x3frag instead
      implicit none
c iflagmode: 0 --> event
c            1 --> pole subtraction
c            2 --> pole compensation
      logical ini(3)
      integer iflagmode,k
      common/xiflagmode/iflagmode
      real*8 x3frag
      common/cx3frag/x3frag
      real*8 an, alpha,beta
      integer ipi
      real*8 XC4,QP24,XDUP4,XDUBP4,XDDP4,XDDBP4,XDSP4,XDCP4
     # ,xdcbp4,XDBP4,XDBBP4,XDTP4,XDTBP4,XDGP4
      real*8 XC,QP2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP,xdcbp
     # ,XDBP,XDBBP,XDTP,XDTBP,XDGP
c This arrays store the old values of x, q2 and densities
      real*8 QP24O(3),XDUP4O(3),XDUBP4O(3),XDDP4O(3),XDDBP4O(3),
     # XDSP4O(3),XDCP4O(3),XDCBP4O(3),XDBP4O(3),XDBBP4O(3),XDTP4O(3),
     # XDTBP4O(3),XDGP4O(3)
      real*8 x3frago(3)
c functions
      real*8 phfragmell
      data ini/.true.,.true.,.true./
c
      save x3frago,QP24O,XDUP4O,XDUBP4O,XDDP4O,XDDBP4O,
     # XDSP4O,XDCP4O,XDCBP4O,XDBP4O,XDBBP4O,XDTP4O,
     # XDTBP4O,XDGP4O,ini
c
      k=iflagmode+1
      if(.not.ini(k).and.qp24.eq.qp24o(k).and.
     #     (iflagmode.eq.2.or.x3frag.eq.x3frago(k)) ) then
         xdup4 = xdup4o(k)
         xdubp4 = xdubp4o(k)
         xddp4 = xddp4o(k)
         xddbp4 = xddbp4o(k)
         xdsp4 = xdsp4o(k)
         xdcp4 = xdcp4o(k)
         xdbp4 = xdbp4o(k)
         xdbbp4 = xdbbp4o(k)
         xdtp4 = xdtp4o(k)
         xdtbp4 = xdtbp4o(k)
         xdgp4 = xdgp4o(k)
         return
      endif
      ini(k)=.false.
      qp2 = qp24
      if(iflagmode.eq.2)then
c returns the first Mellin moment of the charm fragmentation function
         xdup4 = 0.
         xdubp4 = 0.
         xddp4 = 0.
         xddbp4 = 0.
         xdsp4 = 0.
         xdbp4 = 0.
         xdbbp4 = 0.  
         xdtp4 = 0.
         xdtbp4 = 0.  
         xdgp4 = 0.
c Here, divides by two because the main code assumes that anticharm
c fragments like charm, and thus sigma_c and sigma_cbar are summed.
c See the comment at the beginning of this routine
         xdcp4=phfragmell(qp2)/2.d0
         xdcbp4 = 0.
      else
         xc = x3frag
c first argument corresponds to ifragmode=0: It is changed inside
c fonfra. Using -1 causes fonfra to know it's being called from
c phfonfra, and therefore makes it to set the quark type hvqs = 'c'
c In this subroutine iflagmode is already handled properly
         call FONFRA(0,XC,IPI,QP2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP,
     #        xdcbp,XDBP,XDBBP,XDTP,XDTBP,XDGP)
         if(iflagmode.eq.0)then
            xdup4 = xdup
            xdubp4 = xdubp
            xddp4 = xddp
            xddbp4 = xddbp
            xdsp4 = xdsp
c.....because charm=anticharm in the kernel cross sections
c....changed on 29/4/98, MC
            xdcp4 = xdcp/2. + xdcbp/2.
            xdbp4 = xdbp
            xdbbp4 = xdbbp
            xdtp4 = xdtp
            xdtbp4 = xdtbp
            xdgp4 = xdgp
         elseif(iflagmode.eq.1)then
            xdup4 = 0.
            xdubp4 = 0.
            xddp4 = 0.
            xddbp4 = 0.
            xdsp4 = 0.
            xdcp4 = xdcp/2. + xdcbp/2.
            xdbp4 = 0.
            xdbbp4 = 0.  
            xdtp4 = 0.
            xdtbp4 = 0.  
            xdgp4 = 0.
         else
            write(*,*)'unknown iflagmode value = ',iflagmode
            stop
         endif
      endif
      qp24o(k) = qp24
      x3frago(k) = x3frag
      xdup4o(k) = xdup4
      xdubp4o(k) = xdubp4
      xddp4o(k) = xddp4
      xddbp4o(k) = xddbp4
      xdsp4o(k) = xdsp4
      xdcp4o(k) = xdcp4
      xdbp4o(k) = xdbp4
      xdbbp4o(k) = xdbbp4
      xdtp4o(k) = xdtp4
      xdtbp4o(k) = xdtbp4
      xdgp4o(k) = xdgp4
      end
      
      
      FUNCTION FQ1(E,F)
      implicit real*8(a-h,o-z)
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CS/CS1,CS2,CS3
      COMMON/VAR/S,SS,PT
      REAL*8 LAMP
      COMMON/CLA/LAMP
      VE=E
      VF=F
cc      FQ1=LAMP**2/.799814
correction by MC, 4/12/1996, see also corresponding correction in SIGMAPHOTON0
      FQ1=LAMP**2
      RETURN
      END

      SUBROUTINE QGAUS(FUNC,A,B,SS)                                           
      DIMENSION X(5),W(5)                                                     
      DATA X/.1488743389,.4333953941,.6794095682,.8650633666,.9739065285      
     */                                                                       
      DATA W/.2955242247,.2692667193,.2190863625,.1494513491,.0666713443      
     */                                                                       
      XM=0.5*(B+A)                                                            
      XR=0.5*(B-A)                                                            
      SS=0                                                                    
      DO 11 J=1,5                                                             
        DX=XR*X(J)                                                            
        SS=SS+W(J)*(FUNC(XM+DX)+FUNC(XM-DX))                                  
11    CONTINUE                                                                
      SS=XR*SS                                                                
      RETURN                                                                  
      END   
      
      FUNCTION GAUSS(F,A,B,EPS)
C.----------------------------------------------------------------------
C.
C.    GAUSS INTEGRAL OF THE FUNCTION F IN INTERVAL A,B
C.    LAST UPDATE: 12/03/87
C.
C.----------------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION W(12),X(12)
      EXTERNAL F
      DATA CONST/1.E-12/
      DATA W
     &/0.101228536290376, 0.222381034453374, 0.313706645877887,
     & 0.362683783378362, 0.027152459411754, 0.062253523938648,
     & 0.095158511682493, 0.124628971255534, 0.149595988816577,
     & 0.169156519395003, 0.182603415044924, 0.189450610455069/
      DATA X
     &/0.960289856497536, 0.796666477413627, 0.525532409916329,
     & 0.183434642495650, 0.989400934991650, 0.944575023073233,
     & 0.865631202387832, 0.755404408355003, 0.617876244402644,
     & 0.458016777657227, 0.281603550779259, 0.095012509837637/
C--
C--   INITIALISE
      DELTA=CONST*ABS(A-B)
      GAUSS=0.
      AA=A
C--
C--   ITERATION LOOP
   10 Y=B-AA
C--
C--   EPSILON REACHED ??
      IF (ABS(Y).LE.DELTA) RETURN
   20 BB=AA+Y
      C1=0.5*(AA+BB)
      C2=C1-AA
      S8=0.
      S16=0.
      DO 30 I=1,4
         U=X(I)*C2
   30 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 40 I=5,12
         U=X(I)*C2
   40 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF (ABS(S16-S8).GT.EPS*(1.0+ABS(S16))) GOTO 50
      GAUSS=GAUSS+S16
      AA=BB
      GOTO 10
   50 Y=0.5*Y
      IF (ABS(Y).GT.DELTA) GOTO 20
      WRITE (6,9000)
      GAUSS=0.
      RETURN
 9000 FORMAT(1H ,'****** GAUSS... TOO HIGH ACCURACY REQUIRED ******')
      END
 

       FUNCTION GAUSS1(F,A,B,EPS)
C.----------------------------------------------------------------------
C.
C.    GAUSS INTEGRAL OF THE FUNCTION F IN INTERVAL A,B
C.    LAST UPDATE: 12/03/87
C.
C.----------------------------------------------------------------------
      IMPLICIT REAL*4 (A-H,O-Z)
      DIMENSION W(12),X(12)
      EXTERNAL F
      DATA CONST/1.E-12/
      DATA W
     &/0.101228536290376, 0.222381034453374, 0.313706645877887,
     & 0.362683783378362, 0.027152459411754, 0.062253523938648,
     & 0.095158511682493, 0.124628971255534, 0.149595988816577,
     & 0.169156519395003, 0.182603415044924, 0.189450610455069/
      DATA X
     &/0.960289856497536, 0.796666477413627, 0.525532409916329,
     & 0.183434642495650, 0.989400934991650, 0.944575023073233,
     & 0.865631202387832, 0.755404408355003, 0.617876244402644,
     & 0.458016777657227, 0.281603550779259, 0.095012509837637/
C--
C--   INITIALISE
      DELTA=CONST*ABS(A-B)
      GAUSS1=0.
      AA=A
C--
C--   ITERATION LOOP
   10 Y=B-AA
C--
C--   EPSILON REACHED ??
      IF (ABS(Y).LE.DELTA) RETURN
   20 BB=AA+Y
      C1=0.5*(AA+BB)
      C2=C1-AA
      S8=0.
      S16=0.
      DO 30 I=1,4
         U=X(I)*C2
   30 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 40 I=5,12
         U=X(I)*C2
   40 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF (ABS(S16-S8).GT.EPS*(1.0+ABS(S16))) GOTO 50
      GAUSS1=GAUSS1+S16
      AA=BB
      GOTO 10
   50 Y=0.5*Y
      IF (ABS(Y).GT.DELTA) GOTO 20
      WRITE (6,9000)
      GAUSS1=0.
      RETURN
 9000 FORMAT(1H ,'****** GAUSS1... TOO HIGH ACCURACY REQUIRED ******')
      END
 
      FUNCTION PD1(XX)
      implicit real*8(a-h,o-z)
C     EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CE1,PD12,PD14,CE2,CE4
      COMMON/KAON/PKA
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/SOUS/PD12S,PD14S
      REAL*8 LAM
      COMMON/BK2/BO(2)
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      X3=XX
      X2B=(1.-GV)/(X3-AU)
      QB2=FQA
      ACB=PI/2.
      VB=AU/X3
      CALL DISTPRO(GU,GD,GS,GC,GG,X2B,QB2)
      GUA=GU
      GDA=GD
      GSA=GS
      CALL PHFONFRA(X3,IPI,QB2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      ALA=1.-AU/(X3-1.+GV)
      AL=LOG(ALA)
      CE21=CE1(VB,1.d0,2./3.d0,1./3.d0)+AL*CE2(VB,1.d0,2./3.d0,1./3.d0)
      CE12=CE1(VB,1.d0,1./3.d0,2./3.d0)+AL*CE2(VB,1.d0,1./3.d0,2./3.d0)
      CE11=CE1(VB,1.d0,1./3.d0,1./3.d0)+AL*CE2(VB,1.d0,1./3.d0,1./3.d0)
      CE22=CE1(VB,1.d0,2./3.d0,2./3.d0)+AL*CE2(VB,1.d0,2./3.d0,2./3.d0)
      IF (DIF.LT..5) GO TO 1
      DDIA=DDI
      PD1A=(2.*GUA*CE21-GDA*(CE12+CE11))*DDIA*X2B*ACB
      GO TO 2
   1  DFA=DF
      DNFA=DNF
      DSPA=DSP
      PD1A=(2.*CE21+CE22)*((GUA+GSA)*DFA+GSA*DNFA)+(2.*CE12+CE11)*
     !((GDA+GSA)*DNFA+GSA*DFA+2.*GSA*DSPA)+2.*GC*DCP*(2.*CE21+CE22)
      PD1A=PD1A+PKA*GSA*DKA*(2.*CE21-CE12-CE11)
      PD1A=PD1A*X2B*ACB
  2   CE21=CE2(VB,1.d0,2./3.d0,1./3.d0)
      CE12=CE2(VB,1.d0,1./3.d0,2./3.d0)
      CE11=CE2(VB,1.d0,1./3.d0,1./3.d0)
      CE22=CE2(VB,1.d0,2./3.d0,2./3.d0)
      IF (DIF.LT..5) GO TO 3
      DDIA=DDI
      PD12S=(2.*GUA*CE21-GDA*(CE12+CE11))*DDIA*X2B*ACB
      GO TO 4
 3    PD12S=(2.*CE21+CE22)*((GUA+GSA)*DFA+GSA*DNFA)+(2.*CE12+CE11)*
     !((GDA+GSA)*DNFA+GSA*DFA+2.*GSA*DSPA)+2.*GC*DCP*(2.*CE21+CE22)
      PD12S=PD12S+PKA*GSA*DKA*(2.*CE21-CE12-CE11)
      PD12S=PD12S*X2B*ACB
   4  BO(1)=AU/(X3-1.+GV)
      BO(2)=1.
      PD1B=QUADR1(8,1,PD12)
      PD1D=QUADR1(8,1,PD14)
      PD1=PD1A+PD1B+PD1D
      RETURN
      END
      FUNCTION PD12(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CE2
      COMMON/KAON/PKA
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/SOUS/PD12S,PD14S
      REAL*8 LAM
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
       AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      CE21=CE2(VV,WW,2./3.d0,1./3.d0)
      CE12=CE2(VV,WW,1./3.d0,2./3.d0)
      CE11=CE2(VV,WW,1./3.d0,1./3.d0)
      CE22=CE2(VV,WW,2./3.d0,2./3.d0)
      IF (DIF.LT..5) GO TO 1
      DDIB=DDI
      PD12N=(2.*GUB*CE21-GDB*(CE12+CE11))*DDIB*X2*AC
      GO TO 2
   1  DFB=DF
      DNFB=DNF
      DSPB=DSP
      PD12N=(2.*CE21+CE22)*((GUB+GSB)*DFB+GSB*DNFB)+(2.*CE12+CE11)*
     !((GDB+GSB)*DNFB+GSB*DFB+2.*GSB*DSPB)+2.*GC*DCP*(2.*CE21+CE22)
      PD12N=PD12N+PKA*GSB*DKA*(2.*CE21-CE12-CE11)
      PD12N=PD12N*X2*AC
   2  PD12=(PD12N/WW**2-PD12S)/(1.-WW)
      RETURN
      END
      FUNCTION PD14(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CE4
      COMMON/KAON/PKA
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/SOUS/PD12S,PD14S
      REAL*8 LAM
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      CE21=CE4(VV,WW,2./3.d0,1./3.d0)
      CE12=CE4(VV,WW,1./3.d0,2./3.d0)
      CE11=CE4(VV,WW,1./3.d0,1./3.d0)
      CE22=CE4(VV,WW,2./3.d0,2./3.d0)
      IF (DIF.LT..5) GO TO 1
      DDIB=DDI
      PD14N=(2.*GUB*CE21-GDB*(CE12+CE11))*DDIB*X2*AC
      GO TO 2
   1  DFB=DF
      DNFB=DNF
      DSPB=DSP
      PD14N=(2.*CE21+CE22)*((GUB+GSB)*DFB+GSB*DNFB)+(2.*CE12+CE11)*
     !((GDB+GSB)*DNFB+GSB*DFB+2.*GSB*DSPB)+2.*GC*DCP*(2.*CE21+CE22)
      PD14N=PD14N+PKA*GSB*DKA*(2.*CE21-CE12-CE11)
      PD14N=PD14N*X2*AC
  2   PD14=PD14N/WW**2
      RETURN
      END
      FUNCTION PD2(XX)
      implicit real*8(a-h,o-z)
      EXTERNAL FQ1,PD24
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK2/BO(2)
      X3=XX
      BO(1)=AU/(X3-1.+GV)
      BO(2)=1.
      PD2B=QUADR1(8,1,PD24)
      PD2=PD2B
      RETURN
      END
      FUNCTION PD24(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CF4
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      REAL*8 LAM
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GUB=GU
      GDB=GD
      GSB=GS
      CF21=CF4(VV,WW,2./3.d0,1./3.d0)
      CF12=CF4(VV,WW,1./3.d0,2./3.d0)
      CF11=CF4(VV,WW,1./3.d0,1./3.d0)
      CF2M1=CF4(VV,WW,2./3.d0,-1./3.d0)
      CF1M2=CF4(VV,WW,1./3.d0,-2./3.d0)
      CF22=CF4(VV,WW,2./3.d0,2./3.d0)
      CF2M2=CF4(VV,WW,2./3.d0,-2./3.d0)
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      IF (DIF.LT..5) GO TO 1
      DDIB=DDI
      PD24N=(GUB*(CF21-CF2M1)-GDB*(CF12-CF1M2))*DDIB*X2*AC
      GO TO 2
   1  DFB=DF
      DNFB=DNF
      DSPB=DSP
      DKAB=DKA
      CF1M1=CF4(VV,WW,1./3.d0,-1./3.d0)
      PD24N=CF21*((GUB+GSB+GC)*(DFB+DSPB)+(GSB+GC)*(DNFB+DSPB))
      PD24N=PD24N+CF2M1*((GUB+GSB+GC)*(DNFB+DSPB)+(GSB+GC)*(DFB+DSPB))
      PD24N=PD24N+CF12*((GDB+2.*GSB)*(DNFB+DCP)+2.*GSB*(DCP+DFB))
      PD24N=PD24N+CF1M2*((GDB+2.*GSB)*(DFB+DCP)+2.*GSB*(DNFB+DCP))
      PD24N=PD24N+CF11*((GDB+2.*GSB)*DSPB+GSB*(DFB+DNFB))
      PZ=PD24N+CF1M1*((GDB+2.*GSB)*DSPB+GSB*(DFB+DNFB))
      PD24N=PZ+(GUB+2.*GSB)*(CF2M2+CF22)*DCP+GC*(CF22+CF2M2)*(DFB+DNFB)
      PAUX=-CF21*(GUB+GSB)-CF2M1*GSB+CF12*(GDB+2.*GSB)
      PAUX=PAUX+2.*GSB*CF1M2-CF11*GSB-CF1M1*GSB
      PD24N=PD24N+DKAB*PAUX*PKA
      PD24N=PD24N*X2*AC
  2   PD24=PD24N/WW**2
      RETURN
      END
      FUNCTION PI1(XX)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CJ1,PI12,PI14,CJ2
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      REAL*8 LAM
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/BK2/BO(2)
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      X3=XX
      X2B=(1.-GV)/(X3-AU)
      QB2=FQA
      ACB=PI/2.
      VB=AU/X3
      CALL DISTPRO(GU,GD,GS,GC,GG,X2B,QB2)
      GUA=GU
      GDA=GD
      GSA=GS
      CALL PHFONFRA(X3,IPI,QB2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      ALA=1.-AU/(X3-1.+GV)
      AL=LOG(ALA)
      PI1A=CJ1(VB,1.d0)+AL*CJ2(VB,1.d0)
      IF (DIF.LT..5) GO TO 1
      DDIA=DDI
      GB=(4.*GUA-GDA)*X2B*DDIA*ACB/9.
      PI1A=PI1A*GB
      GO TO 2
   1  DFA=DF
      DNFA=DNF
      DSPA=DSP
      GB=((4.*GUA+5.*GSA)*DFA+(GDA+5.*GSA)*DNFA+2.*GSA*DSPA)*X2B*ACB/9.
      GB=GB+(8.*GC*DCP)*X2B*ACB/9.
      GB=GB+PKA*3.*GSA*DKA*X2B*ACB/9.
      PI1A=PI1A*GB
   2  BO(1)=AU/(X3-1.+GV)
      BO(2)=1.
      PI1B=QUADR1(8,1,PI12)
      PI1D=QUADR1(8,1,PI14)
      PI1=PI1A+PI1B+PI1D
      RETURN
      END
      FUNCTION PI12(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CJ2
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      REAL*8 LAM
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      IF(DIF.LT..5) GO TO 1
      DDIB=DDI
      G=(4.*GUB-GDB)*X2*DDIB*AC/9.
      PI12=(G*CJ2(VV,WW)/WW**2-GB*CJ2(VB,1.d0))/(1.-WW)
      GO TO 2
   1  G=((4.*GUB+5.*GSB)*DF+(GDB+5.*GSB)*DNF
     !+2.*GSB*DSP+8.*GC*DCP)*X2*AC/9.
      G=G+PKA*3.*GSB*DKA*X2*AC/9.
      PI12=(G*CJ2(VV,WW)/WW**2-GB*CJ2(VB,1.d0))/(1.-WW)
   2  RETURN
      END
      FUNCTION PI14(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CJ4
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      REAL*8 LAM
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      IF(DIF.LT..5) GO TO 1
      DDIB=DDI
      G=(4.*GUB-GDB)*X2*DDIB*AC/9.
      PI14=G*CJ4(VV,WW)/WW**2
      GO TO 2
   1  G=((4.*GUB+5.*GSB)*DF+(GDB+5.*GSB)*DNF
     !+2.*GSB*DSP+8.*GC*DCP)*X2*AC/9.
      G=G+PKA*3.*GSB*DKA*X2*AC/9.
      PI14=G*CJ4(VV,WW)/WW**2
   2  RETURN
      END
      FUNCTION PI2(XX)
      implicit real*8(a-h,o-z)
      EXTERNAL FQ1,PI24
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK2/BO(2)
      X3=XX
      BO(1)=AU/(X3-1.+GV)
      BO(2)=1.
      PI2B=QUADR1(8,1,PI24)
      PI2=PI2B
      RETURN
      END
      FUNCTION PI24(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CK4
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      IF(DIF.LT..5) GO TO 1
      DDIB=DDI
      G=-(4.*GUB-GDB)*X2*DDIB*AC/9.
      PI24=G*CK4(VV,WW)/WW**2
      GO TO 2
   1  G=((4.*GUB+5.*GSB)*DNF+(GDB+5.*GSB)*DF
     !+2.*GSB*DSP+8.*GC*DCP)*X2*AC/9.
      G=G+PKA*DKA*(4.*GUB-GDB+3.*GSB)*X2*AC/9.
      PI24=G*CK4(VV,WW)/WW**2
   2  RETURN
      END

      FUNCTION QUAX3nn(NBPT,NBINT,Y)
      implicit real*8(a-h,o-z)
      COMMON/BK1/BORNE(2)
      common/varint/xx3,xwv
      common/cx3frag/x3frag
      common/xiflagmode/iflagmode
      real*8 rcn
      common/crcn/rcn
      external y
c iflagmode: 0 --> event
c            1 --> pole subtraction
c            2 --> pole compensation
      integer iflagmode

      xjac=borne(2)-borne(1)
      x3=borne(1)+(borne(2)-borne(1))*xx3
c here, x3frag is the argument of the fragmentation function
      iflagmode=0
      x3frag=x3
      xev=y(x3)*xjac
c subtraction
c rcn is one plus the x power for the pole subtraction
      rcn=5
      iflagmode=1
      x3frag=xx3
      x3=1.d0
      xcnt=y(x3)*x3frag**(rcn-2)
c when iflagmode=2, the actual value of x3 and x3frag are irrelevant
      iflagmode=2
      xcomp=y(x3)
      quax3 = xev - xcnt + xcomp
      return
      END

      FUNCTION QUAX3(NBPT,NBINT,Y)
      implicit real*8(a-h,o-z)
      COMMON/BK1/BORNE(2)
      common/varint/xx3,xwv
      common/cx3frag/x3frag
      common/xiflagmode/iflagmode
      real*8 rcn
      common/crcn/rcn
      external y
c iflagmode: 0 --> event
c            1 --> pole subtraction
c            2 --> pole compensation
      integer iflagmode
      if(abs(borne(2)-1.d0).gt.1.d-8)then
         write(*,*)'ERROR: upper bound in QUAX3 =',borne(2)
         stop
      endif
c      xx=xx3**(0.8d0)
c      xjac0=0.8d0*xx/xx3
      xx=xx3
      xjac0=1.d0
      x3=borne(1)+(borne(2)-borne(1))*xx
      xjac=borne(2)-borne(1)
c here, x3frag is the argument of the fragmentation function
      iflagmode=0
      x3frag=x3
      xev=y(x3)*xjac
c subtraction
c rcn is one plus the x power for the pole subtraction
      rcn=5
      iflagmode=1
      x3frag=xx**(1-borne(1))
      x3=1.d0
      xcnt=y(x3)*x3frag**(rcn-2)*x3frag/xx*(1-borne(1))
c when iflagmode=2, the actual value of x3 and x3frag are irrelevant
      iflagmode=2
      xcomp=y(x3)
      quax3 = (xev - xcnt + xcomp)*xjac0
      return
      END

      FUNCTION QUAX3n(NBPT,NBINT,Y)
      implicit real*8(a-h,o-z)
      COMMON/BK1/BORNE(2)
      common/varint/xx3,xwv
      common/cx3frag/x3frag
      common/xiflagmode/iflagmode
      real*8 rcn
      common/crcn/rcn
      external y
c iflagmode: 0 --> event
c            1 --> pole subtraction
c            2 --> pole compensation
      integer iflagmode

      xjac=borne(2)-borne(1)
      x3=borne(1)+(borne(2)-borne(1))*xx3
c here, x3frag is the argument of the fragmentation function
      iflagmode=0
      x3frag=x3
      xev=y(x3)*xjac
c subtraction
c rcn is one plus the x power for the pole subtraction
      rcn=5
      iflagmode=1
      x3frag=x3
      x3=1.d0
      xcnt=y(x3)*x3frag**(rcn-2)
c extra range
      x3frag=borne(1)*xx3
      xjac=borne(1)
      xexr=y(x3)*x3frag**(rcn-2)*xjac
c when iflagmode=2, the actual value of x3 and x3frag are irrelevant
      iflagmode=2
      xcomp=y(x3)
      quax3 = xev - xcnt + (xcomp - xexr)
      return

      END


      FUNCTION QUADR1(NBPT,NBINT,Y)
      implicit real*8(a-h,o-z)
      parameter (small=1.d-6)
      COMMON/BK2/BORNE(2)
      common/varint/xx3,xwv
      external y
c
      xarg=borne(1)+(borne(2)-borne(1))*xwv
      if(xarg.eq.1.d0)then
         tmp=0.d0
      else
         xjac=borne(2)-borne(1)
         tmp=y(xarg)*xjac
      endif
      quadr1=tmp
      return
      END


      FUNCTION SB1(XX)
      implicit real*8(a-h,o-z)
c photon-quark ->quark contrib
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/KAON/PKA
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      X3=XX
      V=AU/X3
      X2B=(1.-GV)/(X3*(1.-V))
      Q2=FQA
      AS=PI
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      CALL DISTPRO(GU,GD,GS,GC,GG,X2B,Q2)
      GC=FLCA*GC
      IF(DIF.LT..5) GO TO 1
      SB1=(4.*GU-GD)*DDI/(9.*X3)
      GO TO 2
   1  SB1=(4.*GU+5.*GS)*DF/(X3*9.)
      SB1=SB1+(GD+5.*GS)*DNF/(X3*9.)
      SB1=SB1+(2.*GS*DSP+8.*GC*DCP)/(X3*9.)
      SB1=SB1+PKA*3.*GS*DKA/X3/9.
   2  SB1=SB1*2.*AS*(1.+(1.-V)**2)/(1.-V)
      SB1=SB1*CF
      RETURN
      END
      FUNCTION SB2(XX)
      implicit real*8(a-h,o-z)
c photon quark -> gluon contrib
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      EXTERNAL FQ1
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      X3=XX
      V=AU/X3
      X2B=(1.-GV)/(X3*(1.-V))
      Q2=FQB
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      CALL DISTPRO(GU,GD,GS,GC,GG,X2B,Q2)
      GC=GC*FLCA
      AS=PI
      SB2=(4.*GU+GD+10.*GS+2.*GS+8.*GC)/9.
      SB2=SB2*DGP/X3
      SB2=SB2*2.*AS*(1.+V**2)/V
      SB2=SB2*CF
      RETURN
      END
      FUNCTION SB3(XX)
      implicit real*8(a-h,o-z)
c gluon photon -> quark born contrib.
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/KAON/PKA
      COMMON/CST/LAM,PI,CAS
      COMMON/VA/T,U,GV,GW,AU,X3
      EXTERNAL FQ1
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      X3=XX
      V=AU/X3
      X2B=(1.-GV)/(X3*(1.-V))
      Q2=FQC
      CALL DISTPRO(GU,GD,GS,GC,GG,X2B,Q2)
      AS=PI
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      SB3=(5.*(DF+DNF)+2.*DSP+8.*DCP)*GG/(X3*9.)
      SB3=SB3-PKA*3.*GG*DKA/X3/9.
      SB3=SB3*(3./4.)*AS*(V**2+(1.-V)**2)/(V*(1.-V))
      SB3=SB3*CF
      RETURN
      END
      FUNCTION SC1(XX)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,SC12,SC13,SC14,CB1,CB2,CB3
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/BK2/BO(2)
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      X3=XX
      X2B=(1.-GV)/(X3-AU)
      QB2=FQA
      ACB=PI/2.
      VB=AU/X3
      CALL DISTPRO(GU,GD,GS,GC,GG,X2B,QB2)
      GUA=GU
      GDA=GD
      GSA=GS
      CALL PHFONFRA(X3,IPI,QB2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      ALA=1.-AU/(X3-1.+GV)
      AL=LOG(ALA)
      SC1A=CB1(VB,1.d0)+AL*CB2(VB,1.d0)+.5*AL**2*CB3(VB,1.d0)
      IF(DIF.LT..5) GO TO 1
      DDIA=DDI
      GB=(4.*GUA-GDA)*X2B*DDIA*ACB/9.
      SC1A=SC1A*GB
      GO TO 2
   1  DFA=DF
      DNFA=DNF
      DSPA=DSP
      GB=(4.*GUA+5.*GSA)*DFA/9.+(GDA+5.*GSA)*DNFA/9.+2.*GSA*DSPA/9.
      GB=GB+8.*GC*DCP/9.
      GB=GB+PKA*3.*GSA*DKA/9.
      GB=GB*X2B*ACB
      SC1A=SC1A*GB
   2  BO(1)=AU/(X3-1.+GV)
      BO(2)=1.
      SC1B=QUADR1(8,1,SC12)
      SC1C=QUADR1(8,1,SC13)
      SC1D=QUADR1(8,1,SC14)
      SC1=SC1A+SC1B+SC1C+SC1D
      SC1=SC1*CF
      RETURN
      END
      FUNCTION SC12(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CB2
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      IF(DIF.LT..5) GO TO 1
      DDIB=DDI
      G=(4.*GUB-GDB)*X2*DDIB*AC/9.
      SC12=(G*CB2(VV,WW)/WW**2-GB*CB2(VB,1.d0))/(1.-WW)
      GO TO 2
   1  G=((4.*GUB+5.*GSB)*DF+(GDB+5.*GSB)*DNF
     !+2.*GSB*DSP+8.*GC*DCP)*X2*AC/9.
      G=G+PKA*3.*GSB*DKA*X2*AC/9.
      SC12=(G*CB2(VV,WW)/WW**2-GB*CB2(VB,1.d0))/(1.-WW)
   2  RETURN
      END
      FUNCTION SC13(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CB3
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      IF(DIF.LT..5) GO TO 1
      DDIB=DDI
      G=(4.*GUB-GDB)*X2*DDIB*AC/9.
      SC13=(G*CB3(VV,WW)/WW**2-GB*CB3(VB,1.d0))*LOG(1.-WW)/(1.-WW)
      GO TO 2
   1  G=((4.*GUB+5.*GSB)*DF+(GDB+5.*GSB)*DNF
     !+2.*GSB*DSP+8.*GC*DCP)*X2*AC/9.
      G=G+PKA*3.*GSB*DKA*X2*AC/9.
      SC13=(G*CB3(VV,WW)/WW**2-GB*CB3(VB,1.d0))*LOG(1.-WW)/(1.-WW)
   2  RETURN
      END
      FUNCTION SC14(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CB4
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      IF(DIF.LT..5) GO TO 1
      DDIB=DDI
      G=(4.*GUB-GDB)*X2*DDIB*AC/9.
      SC14=G*CB4(VV,WW)/WW**2
      GO TO 2
   1  G=((4.*GUB+5.*GSB)*DF+(GDB+5.*GSB)*DNF
     !+2.*GSB*DSP+8.*GC*DCP)*X2*AC/9.
      G=G+PKA*3.*GSB*DKA*X2*AC/9.
      SC14=G*CB4(VV,WW)/WW**2
   2  RETURN
      END
      FUNCTION SC2(XX)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,SC22,SC23,SC24,CA1,CA2,CA3
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/BK2/BO(2)
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      X3=XX
      X2B=(1.-GV)/(X3-AU)
      QB2=FQB
      ACB=PI/2.
      VB=AU/X3
      CALL DISTPRO(GU,GD,GS,GC,GG,X2B,QB2)
      GC=GC*FLCA
      GUA=GU
      GDA=GD
      GSA=GS
      CALL PHFONFRA(X3,IPI,QB2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      ALA=1.-AU/(X3-1.+GV)
      AL=LOG(ALA)
      SC2A=CA1(VB,1.d0)+AL*CA2(VB,1.d0)+.5*AL**2*CA3(VB,1.d0)
      DGPA=DGP
      GB=(4.*GUA+GDA+10.*GSA+2.*GSA+8.*GC)*X2B*DGPA*ACB/9.
      SC2A=SC2A*GB
      BO(1)=AU/(X3-1.+GV)
      BO(2)=1.
      SC2B=QUADR1(8,1,SC22)
      SC2C=QUADR1(8,1,SC23)
      SC2D=QUADR1(8,1,SC24)
      SC2=SC2A+SC2B+SC2C+SC2D
      SC2=SC2*CF
      RETURN
      END
      FUNCTION SC22(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CA2
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQB
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GC=GC*FLCA
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      G=(4.*GUB+GDB+10.*GSB+2.*GSB+8.*GC)*X2*DGP*AC/9.
      SC22=(G*CA2(VV,WW)/WW**2-GB*CA2(VB,1.d0))/(1.-WW)
   2  RETURN
      END
      FUNCTION SC23(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CA3
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQB
       AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GC=GC*FLCA
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      G=(4.*GUB+GDB+10.*GSB+2.*GSB+8.*GC)*X2*DGP*AC/9.
      SC23=(G*CA3(VV,WW)/WW**2-GB*CA3(VB,1.d0))*LOG(1.-WW)/(1.-WW)
   2  RETURN
      END

      FUNCTION SC24(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      EXTERNAL FQ1,CA4
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQB
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GC=GC*FLCA
      GUB=GU
      GDB=GD
      GSB=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      G=(4.*GUB+GDB+10.*GSB+2.*GSB+8.*GC)*X2*DGP*AC/9.
      SC24=G*CA4(VV,WW)/WW**2
   2  RETURN
      END
      FUNCTION SCG(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL SC2
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SCG=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SCG=QUAX3(8,1,SC2)*AU**2
   2  CONTINUE
      RETURN
      END
      FUNCTION SCOD(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL SB1
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SCOD=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SCOD=QUAX3(12,1,SB1)*AU**2*(1.-GV)
   2  CONTINUE
      RETURN
      END
      FUNCTION SCOG(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL SB2
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SCOG=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SCOG=QUAX3(12,1,SB2)*AU**2*(1.-GV)
   2  CONTINUE
      RETURN
      END
      FUNCTION SCQ(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      EXTERNAL SC1
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY
      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SCQ=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SCQ=QUAX3(8,1,SC1)*AU**2
   2  CONTINUE
      RETURN
      END
      

      subroutine sigmaphotonwb(av,sd,ncall0,nitn0)
c Returns the cross section for photon-hadron collisions; the energy of
c the photons is distributed according to the WW function
      implicit none
      integer ncall0,nitn0
      real*8 xl,xu,acc
      integer ndim,ncall,itmx,nprn
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      real*8 av,sd,chi2
      integer j
      real*8 sigmaphoton2
      external sigmaphoton2
      ndim=3
      acc=0.001
      do j=1,ndim
         xl(j)=0
         xu(j)=1
      enddo
      ncall=ncall0
      itmx=nitn0
      call vegas(sigmaphoton2,av,sd,chi2)
      end


      function sigmaphoton2(xx,weight)
c Integrates wide-band photon cross section
      implicit none
      real*8 sigmaphoton2
      real*8 xx(*)
      real*8 weight
      real*8 xx3,xwv
      common/varint/xx3,xwv
      real*8 zmin1,zmax
      common/wwrange/zmin1,zmax
      real*8 z_ww,ansect,xjac,fww_ww
c
      xx3=xx(1)
      xwv=xx(2)
      z_ww=zmin1+(zmax-zmin1)*xx(3)**2
      xjac = (zmax-zmin1)*2*xx(3)
      call sigmaphoton0(z_ww,ansect)
      sigmaphoton2=fww_ww(z_ww)*ansect*xjac
      end


      subroutine sigmaphotonm(av,sd,ncall0,nitn0)
c Returns the cross section for monochromatic photon-hadron collisions
      implicit none
      integer ncall0,nitn0
      real*8 xl,xu,acc
      integer ndim,ncall,itmx,nprn
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      real*8 av,sd,chi2
      integer j
      real*8 sigmaphoton1
      external sigmaphoton1
      ndim=2
      acc=0.001
      do j=1,ndim
         xl(j)=0
         xu(j)=1
      enddo
      ncall=ncall0
      itmx=nitn0
      nprn=0
      call vegas(sigmaphoton1,av,sd,chi2)
      end


      function sigmaphoton1(xx,weight)
c Integrates monochromatic photon cross section
      implicit none
      real*8 sigmaphoton1
      real*8 xx(*)
      real*8 weight
      real*8 xx3,xwv
      common/varint/xx3,xwv
      real*8 z_ww,ansect
      xx3=xx(1)
      xwv=xx(2)
      z_ww=1.d0
      call sigmaphoton0(z_ww,ansect)
      sigmaphoton1=ansect
      end

      
**************************************** SIGMAPHOTON0 ********************
      SUBROUTINE SIGMAPHOTON0(ZZ,ANSECT)
      implicit real*8(a-h,o-z)
c Photon hadron cross section; zz is the photon energy fraction
c i.e. =1 for monochromatic photon, =z to be integrated with
c the Weizsaeker-Williams function for the electron.
c ANSECT is the return value
c 
      REAL*8 LAM,IS,NC,NF,LLD,LLQ,LLG,HO,LLR,LQ
      REAL*8 LAMP
      COMMON/CUSC/rusc,USC,PUSC
      COMMON/CST/LAM,PI,CAS
      COMMON/CLA/LAMP
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/VAR/S,SS,PT
      COMMON/KAON/PKA
      COMMON/FSM/FGG
      COMMON/CS/CS1,CS2,CS3
      COMMON/GLUN/DGG
      COMMON/CHARM/FLCA,IPI
      EXTERNAL SCOD,SLD,SCQ,SPDI,SPIA,SPIB,SCOG,SFU,SLLQ,SLLG,SCG
      EXTERNAL SFQ,SFG,SPDF,SLR,ALPS,FQ1,SLQ
c      EXTERNAL FOPTT,FOPTA
      EXTERNAL FOPTT
      COMMON/FOPT/RAE,RBE,RCE,RTE
      COMMON/GPRO/PAR(30),CAPR(8,20,32)
      COMMON/GPHO/PAR1(30),PAR2(30),CAPOL(8,20,32),CAVDM(8,20,32)
      REAL*8 AZ,BZ,EPS,ETA,FOPTT,SOR
      COMMON/NEWP/SCM,SSCM,PTM,YCM
      COMMON/FLAG/HO,PLLA
      common/quarkmass/qm,xmq,ims

      RX(X)=X*(2.+X*LOG(X/(1.+X)))
      CQCD=1.54
      PT=PTM
      IF(PUSC.LT.1.5) LMAX=12
      S=SCM*ZZ
      SS=SQRT(S)
      YA=YCM-LOG(ZZ)/2.
      
      raplim = dasinh(sqrt((s+qm**2)**2/4.d0
     #         /s/(ptm**2+qm**2) - 1.d0))
      if(abs(ya).gt.raplim) then
c          print*,' Out of kinematical limits, raplim = ',raplim,
c     #           '   ya = ',ya
          ansect = 0.
          return
      endif

      COEF=389.e6/137.d0
      COEP=COEF/PT**4

c....correction by MC, 4/12/1996 and then 13/12/96 for the addition of ims flag
cc      A=ALPS(PT**2*rusc)/PI
cc      print*,'alpha_s = ',aLPS(PT**2*rusc)
      if(ims.eq.1) then
           rensc2 = (xmq**2 + pt**2)*rusc
      else
           rensc2 = pt**2*rusc
      endif
      A=ALPS(rensc2)/PI
c      print*,'alpha_s = ',ALPS(rensc2)

      A2=A**2
c      print*,' before scod'
      COD=COEP*SCOD(YA)
c      print*,'cod in sigmaphoton0 ',cod

c This is idiotic
c      IF(abs(COD).LE.1E-15) then
c        bor_dir = 0.
c        bor_res = 0.
c        ansect = 0.
c        GO TO 3
c      endif

c.....LO resolved
      IF(PLLA.LT..5) GO TO 11
      LLD=COEP*SLD(YA)
      LLR=COEP*SLR(YA)
      GO TO 12



   11 LLD=0.
      LLR=0.
   12 IF(HO.LT..5) GO TO 13
      CQ=COEP*SCQ(YA)
      PDI=COEP*SPDI(YA)
      PDF=COEP*SPDF(YA)
      PIA=COEP*SPIA(YA)
      PIB=COEP*SPIB(YA )
      GO TO 14
   13 CQ=0.
      PDI=0.
      PDF=0.
      PIA=0.
      PIB=0.

   14 IF(DIF.GT..5) GO TO 19
      COG=COEP*SCOG(YA)
      FU=COEP*SFU(YA)
      

c.....LO resolved
      IF(PLLA.LT..5) GO TO 16
      LLQ=COEP*SLLQ(YA)
      LQ=COEP*SLQ(YA)
      LLG=COEP*SLLG(YA)
      GO TO 17




   16 LLQ=COEP*0.
      LLG=COEP*0.
      LQ=0.
   17 IF(HO.LT..5) GO TO 18
      CG=COEP*SCG(YA)
      FQ=COEP*SFQ(YA)
      FG=COEP*SFG(YA)
      GO TO 19
   18 CG=0.
      FQ=0.
      FG=0.
  19  GO TO 10
  10  BOR_dir = (COG+COD+FU)
      BOR_res = (LLR+LLD + LLQ+LQ+LLG)
      HO_dir  = (CQ+PDI+PDF+PIA+PIB+CG+FQ+FG)

c      xt = a*cqcd
c      ansect = (bor_dir/cqcd)*rx(xt) + ho_dir*a2

      IF(HO.LT..5) then

      ansect = bor_dir*a

      ansect = ansect + bor_res*a2

      else

      b0 = (33. - 2.*nf)/6.

cc....see eq. 4.6 of Nucl. Phys. B 286 (1987) 553 
cc      ansect = a*(bor_dir + .5*b0*log(PT**2*rusc/(lamp**2/.799814))
cc     #            *a*bor_dir) + ho_dir*a2

c....correction by MC, 4/12/1996, see also corresponding correcion in FQ1
      ansect = a*(bor_dir + .5*b0*log(rensc2/lamp**2)
     #            *a*bor_dir) + ho_dir*a2

      ansect = ansect + bor_res*a2

      endif

  3   CONTINUE

c      if(ho.eq.0) then
c	write(42,'(" born_dir = ",e12.6)') bor_dir*a
c	write(42,'(" born_res = ",e12.6)') bor_res*a2
c      else
c	write(42,'(" ho_dir  = ",e12.6)') ansect
c      endif

  2   CONTINUE
 1      CONTINUE
   7  FORMAT ( )
      RETURN
      END



      FUNCTION SF1(XX)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
c      EXTERNAL FQ1,SF12,SF13,SF14,CC1,CC2,CC3
      EXTERNAL FQ1,SF12,SF13,SF14,CC1,CC2,CC3
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/BK2/BO(2)
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      X3=XX
      X2B=(1.-GV)/(X3-AU)
      QB2=FQC
      ACB=PI/2.
      VB=AU/X3
      CALL DISTPRO(GU,GD,GS,GC,GG,X2B,QB2)
      ALA=1.-AU/(X3-1.+GV)
      AL=LOG(ALA)
      SF1A=CC1(VB,1.d0)+AL*CC2(VB,1.d0)+.5*AL**2*CC3(VB,1.d0)
      CALL PHFONFRA(X3,IPI,QB2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      DFA=DF
      DNFA=DNF
      DSPA=DSP
      GB=(5.*(DFA+DNFA)+2.*DSPA+8.*DCP)*GG*X2B*ACB/9.
      GB=GB-PKA*3.*DKA*GG*X2B*ACB/9.
      SF1A=SF1A*GB
   2  BO(1)=AU/(X3-1.+GV)
      BO(2)=1.
      SF1B=QUADR1(8,1,SF12)
      SF1C=QUADR1(8,1,SF13)
      SF1D=QUADR1(8,1,SF14)
      SF1=SF1A+SF1B+SF1C+SF1D
      RETURN
      END
      FUNCTION SF12(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
c      EXTERNAL FQ1,CC2
      EXTERNAL FQ1,CC2
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQC
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      DFB=DF
      DNFB=DNF
      DSPB=DSP
      G=(5.*(DFB+DNFB)+2.*DSPB+8.*DCP)*GG*X2*AC/9.
      G=G-PKA*3.*DKA*GG*X2*AC/9.
      SF12=(G*CC2(VV,WW)/WW**2-GB*CC2(VB,1.d0))/(1.-WW)
      RETURN
      END
      FUNCTION SF13(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
c      EXTERNAL FQ1,CC3
      EXTERNAL FQ1,CC3
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQC
       AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      DFB=DF
      DNFB=DNF
      DSPB=DSP
      G=(5.*(DFB+DNFB)+2.*DSPB+8.*DCP)*GG*X2*AC/9.
      G=G-PKA*3.*DKA*GG*X2*AC/9.
      SF13=(G*CC3(VV,WW)/WW**2-GB*CC3(VB,1.d0))*LOG(1.-WW)/(1.-WW)
      RETURN
      END
      FUNCTION SF14(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
c      EXTERNAL FQ1,CC4
      EXTERNAL FQ1,CC4
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQC
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      G=(5.*(DF+DNF)+2.*DSP+8.*DCP)*GG*X2*AC/9.
      G=G-PKA*3.*DKA*GG*X2*AC/9.
      SF14=G*CC4(VV,WW)/WW**2
      RETURN
      END
      FUNCTION SF2(XX)
      implicit real*8(a-h,o-z)
c      EXTERNAL FQ1,SF24
      EXTERNAL FQ1,SF24
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK2/BO(2)
      X3=XX
      BO(1)=AU/(X3-1.+GV)
      BO(2)=1.
      SF2B=QUADR1(8,1,SF24)
      SF2=SF2B
      RETURN
      END
      FUNCTION SF24(WA)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
c      EXTERNAL FQ1,CD4
      EXTERNAL FQ1,CD4
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/W/GB,GUA,GDA,GSA,DDIA,DFA,DNFA,DSPA,VB
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      WW=WA
      VV=AU/(WW*X3)
      X2=(1.-GV)/(X3*(1.-VV))
      Q2=FQA
      AC=PI/2.
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DGP=XDGP/X3
      DDI=.0
      DKA=.0       
      PDS=6.+FLCA*4.
      SF24=GG*DGP*CD4(VV,WW)*X2*AC*PDS/(9.*WW**2)
      RETURN
      END
      FUNCTION SFG(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL SF2
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SFG=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SFG=QUAX3(8,1,SF2)*AU**2
   2  CONTINUE
      RETURN
      END
      FUNCTION SFQ(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL SF1
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SFQ=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SFQ=QUAX3(8,1,SF1)*AU**2
   2  CONTINUE
      RETURN
      END
      FUNCTION SFU(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL SB3
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SFU=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SFU=QUAX3(12,1,SB3)*AU**2*(1.-GV)
   2  CONTINUE
      RETURN
      END
       
      FUNCTION SLD(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      EXTERNAL SLD1
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY
      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SLD=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SLD=QUAX3(24,1,SLD1)*AU*(1.-GV)
   2  CONTINUE
      RETURN
      END
      FUNCTION SLD1(XX)
      implicit real*8(a-h,o-z)
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK2/BO(2)
      EXTERNAL SLD2
      X3=XX
      BO(1)=AU/X3
      BO(2)=1.-(1.-GV)/X3
      SLD1=QUADR1(24,1,SLD2)
      RETURN
      END
      FUNCTION SLD2(VV)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
c      EXTERNAL FQ1
      EXTERNAL FQ1
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      V=VV
      V1=1.-V
      V2=V**2
      V12=V1**2
      X1=AU/(V*X3)
      X2=(1.-GV)/(V1*X3)
      Q2=FQA
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)      
      GC=GC*FLCA
      CALL DISTGAM(GPU,GPD,GPC,GPG,X1,Q2) 
      GPC=GPC*FLCA
      AS=PI
      PU=GPU
       PD=GPD
       UV=GU
       DV=GD
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0
       SI1=(4./9.)*(1.+V2)/V12
       SI1C=(4./9.)*(1.+V12)/V2
       SI2=(4./9.)*((1.+V2)/V12+(1.+V12)/V2)-(8./27.)/(V1*V)
       SI3=(4./9.)*((1.+V2)/V12+V2+V12)+(8./27.)*V2/V1
       SI3C=(4./9.)*((1.+V12)/V2+V2+V12)+(8./27)*V12/V
      SI6=(4./9.)*(V2+V12)
       IF(DIF.LT..5) GO TO 1
      SLD2=(2.*(PD*(2.*UV-DV)-PU*DV))*SI1C
     .+(PU*UV-PD*DV)*(SI2+SI3C-SI3)
      SLD2=SLD2*AS**2*DDI
      GO TO 2
   1   SE=GS
      DFA=DF
      DNFA=DNF
      DSPA=DSP
      DKAA=DKA*PKA
      UVS=UV+SE
      DVS=DV+SE
      SAA=(PU*(DV+4.*SE)+PD*(UV+4.*SE))*DFA+2.*SE*(UVS+DVS+2.*SE)*DSPA
      SAA=SAA+(PD*(UV+4.*SE)+PU*(DV+4.*SE))*DNFA
      SAA=SAA+(PU*(DV+4.*SE)-PD*(UV+4.*SE))*DKAA
      SAA=SAA+2.*GPC*DCP*(UV+DV+6.*SE)
c.....e' 6 o 4 il coeff. del mare??????????????????????? (vers. orig. 6)
c      SAA=SAA+2.*GPC*DCP*(UV+DV+4.*SE)
      SAB=(4.*PD*UVS+2.*(PU+PD)*SE)*DFA+4.*(PU+PD)*SE*DSPA
      SAB=SAB+(2.*(PU+PD)*DVS+2.*PD*SE)*DNFA
      SAB=SAB+2.*DKAA*SE*(PD-PU)
      SAB=SAB+4.*GC*(PU+PD+PD)*DCP
      SAC=(PU*UVS+PD*SE)*DFA+(PU*SE+PD*DVS)*DNFA+2.*PD*SE*DSPA
      SAC=SAC+(PU-PD)*SE*DKAA
      SAC=SAC+2.*GPC*GC*DCP
      SAD=(PU*SE+PD*DVS)*DFA+(PU*UVS+PD*SE)*DNFA+2.*PD*SE*DSPA
      SAD=SAD+(PU*UVS-PD*DVS)*DKAA
      SAD=SAD+2.*GPC*GC*DCP
      SLC=(PU*(UVS+SE)+PD*(DVS+SE))*(DFA+DNFA+DSPA)+4.*SE*(DFA+DNFA)*PD
      SLC=SLC+(PD*(DVS+SE)-PU*(UVS+SE))*DKAA
      SLC=SLC+2.*(PU*(UV+2.*SE)+PD*(DV+2.*SE)+PD*(2.*SE))*DCP
      SLD2=(SAA*SI1+SAB*SI1C+SAC*(SI2+SI3C)+SAD*SI3+SLC*SI6)*AS**2
   2  CONTINUE
      RETURN
      END
      FUNCTION SLLG(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      EXTERNAL SLLG1
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY
      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SLLG=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SLLG=QUAX3(12,1,SLLG1)*AU*(1.-GV)
   2  CONTINUE
      RETURN
      END
      FUNCTION SLLG1(XX)
      implicit real*8(a-h,o-z)
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK2/BO(2)
      EXTERNAL SLLG2
      X3=XX
      BO(1)=AU/X3
      BO(2)=1.-(1.-GV)/X3
      SLLG1=QUADR1(12,1,SLLG2)
      RETURN
      END
      FUNCTION SLLG2(VV)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
c      EXTERNAL FQ1
      EXTERNAL FQ1
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      V=VV
      V1=1.-V
      V2=V**2
      V12=V1**2
      X1=AU/(V*X3)
      X2=(1.-GV)/(V1*X3)
      Q2=FQB
      AS=PI
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      GC=GC*FLCA
      CALL DISTGAM(GPU,GPD,GPC,GPG,X1,Q2)
      PU=GPU                                                
      PD=GPD
       UV=GU
       DV=GD
      SE=GS
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      SI8=0.
      SI7=(32./27.)*(V/V1+V1/V)-(8./3.)*(V2+V12)
      DGA=DGP
      UVS=UV+SE
      DVS=DV+SE
      SLLG2=(2.*GPC*GC+PU*(UVS+SE)+PD*(DVS+SE)+2.*PD*SE)*SI7
      SLLG2=SLLG2*DGA*AS**2
      RETURN
      END
      FUNCTION SLLQ(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      EXTERNAL SLLQ1
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY
      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SLLQ=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SLLQ=QUAX3(12,1,SLLQ1)*AU*(1.-GV)
   2  CONTINUE
      RETURN
      END
      FUNCTION SLLQ1(XX)
      implicit real*8(a-h,o-z)
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK2/BO(2)
      EXTERNAL SLLQ2
      X3=XX
      BO(1)=AU/X3
      BO(2)=1.-(1.-GV)/X3
      SLLQ1=QUADR1(12,1,SLLQ2)
      RETURN
      END
      FUNCTION SLLQ2(VV)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
      COMMON/KAON/PKA
c      EXTERNAL FQ1
      EXTERNAL FQ1
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      V=VV
      V1=1.-V
      V2=V**2
      V12=V1**2
      X1=AU/(V*X3)
      X2=(1.-GV)/(V1*X3)
      Q2=FQC
      AS=PI
      CALL DISTGAM(GPU,GPD,GPC,GPG,X1,Q2)
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      PU=GPU
      PD=GPD
      GP=GG
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      SI4=(V2+1.)/V12+(4./9.)*(V2+1.)/V
      SI4C=0.
      DFA=DF
      DNFA=DNF
      DSPA=DSP
      SLA=((PU+PD)*(DFA+DNFA)+2.*PD*DSPA+2.*GPC*DCP)*GP
      SLA=SLA+PKA*(PU-PD)*DKA*GP
      SLC=0.
      SLLQ2=(SLA*SI4+SLC*SI4C)*AS**2
      RETURN
      END
      FUNCTION SLQ(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      EXTERNAL SLQ1,ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY
      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SLQ=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SLQ=QUAX3(12,1,SLQ1)*AU*(1.-GV)
   2  CONTINUE
      RETURN
      END
      FUNCTION SLQ1(XX)
      implicit real*8(a-h,o-z)
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK2/BO(2)
      EXTERNAL SLQ2
      X3=XX
      BO(1)=AU/X3
      BO(2)=1.-(1.-GV)/X3
      SLQ1=QUADR1(12,1,SLQ2)
      RETURN
      END
      FUNCTION SLQ2(VV)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
c      EXTERNAL FQ1
      EXTERNAL FQ1
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      V=VV
      V1=1.-V
      V2=V**2
      V12=V1**2
      X1=AU/(V*X3)
      X2=(1.-GV)/(V1*X3)
      Q2=FQA
      AS=PI
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      CALL DISTGAM(GPU,GPD,GPC,GPG,X1,Q2)
      PU=GPU
      PD=GPD
      GP=GG
      GPC=GPC*FLCA
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      SI4=0.
      SI4C=(V12+1.)/V2+(4./9.)*(1.+V12)/V1
      DGA=DGP
      SLA=0.
      SLC=2.*(GPC+PU+2.*PD)*GP*DGA
      SLQ2=(SLA*SI4+SLC*SI4C)*AS**2
      RETURN
      END
      FUNCTION SLR(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      EXTERNAL SLR1,ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY
      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SLR=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SLR=QUAX3(12,1,SLR1)*AU*(1.-GV)
   2  CONTINUE
      RETURN
      END
      FUNCTION SLR1(XX)
      implicit real*8(a-h,o-z)
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK2/BO(2)
      EXTERNAL SLR2
      X3=XX
      BO(1)=AU/X3
      BO(2)=1.-(1.-GV)/X3
      SLR1=QUADR1(12,1,SLR2)
      RETURN
      END
      FUNCTION SLR2(VV)
      implicit real*8(a-h,o-z)
C      EXTERNAL GU,GD,GS,GG,DDI,DF,DNF,DSP,DF0,DS0,DGP,DG0
      real*8 is,nc,lam,nf
      COMMON/PAR/FQQ,FQG,FGQ,FQP,DQQ,DQG,DGQ,FA,IS,DIF,PIP,PTL,CF,NC,NF
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/CST/LAM,PI,CAS
c      EXTERNAL FQ1
      EXTERNAL FQ1
      COMMON/EC/FQA,FQB,FQC,FQE
      COMMON/CHARM/FLCA,IPI
      V=VV
      V1=1.-V
      V2=V**2
      V12=V1**2
      X1=AU/(V*X3)
      X2=(1.-GV)/(V1*X3)
      Q2=FQE
      AS=PI
      CALL DISTPRO(GU,GD,GS,GC,GG,X2,Q2)
      CALL DISTGAM(GPU,GPD,GPC,GPG,X1,Q2)
      PG=GPG
      GP=GG
      UV=GU
      DV=GD
      CALL PHFONFRA(X3,IPI,Q2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     *,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      DF=XDUP/X3
      DNF=XDDP/X3
      DSP=XDSP/X3
      DGP=XDGP/X3
      DCP=FLCA*XDCP/X3
      DDI=.0
      DKA=.0       
      SI4C=(V12+1.)/V2+(4./9.)*(1.+V12)/V1
      SI5=(V2+V12)/(V*V1*6.)-(3./8.)*(V2+V12)
      SI4=(V2+1.)/V12+(4./9.)*(V2+1.)/V
      SI8=(9./2.)*(3.-V*V1+V/V12+V1/V2)
       IF(DIF.LT..5) GO TO 1
      SLR2=PG*(UV-DV)*SI4C
      SLR2=SLR2*AS**2*DDI
      GO TO 2
   1   SE=GS
      DFA=DF
      DNFA=DNF
      DSPA=DSP
      DGA=DGP
      UVS=UV+SE
      DVS=DV+SE
      SAE=PG*(2.*GC*DCP+(UVS+SE)*DFA+(DVS+SE)*DNFA+2.*SE*DSPA)
      SLB=2.*PG*GP*(DFA+DNFA+DSPA+DCP)
      SLF=PG*(UV+DV+6.*SE+2.*GC)*DGA
      SLG=PG*GP*DGA
      SLR2=(SAE*SI4C+SLB*SI5+SLF*SI4+SLG*SI8)*AS**2
c      print*,'slr2_photon',x1,pg
c      print*,'slr2_proton',x2,gp
   2  CONTINUE
      RETURN
      END
      FUNCTION SPDF(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL PD2
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SPDF=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SPDF=QUAX3(8,1,PD2)*AU**2
   2  CONTINUE
      RETURN
      END
      FUNCTION SPDI(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL PD1
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SPDI=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SPDI=QUAX3(8,1,PD1)*AU**2
   2  CONTINUE
      RETURN
      END
      FUNCTION SPIA(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL PI1
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SPIA=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SPIA=QUAX3(8,1,PI1)*AU**2
   2  CONTINUE
      RETURN
      END
      FUNCTION SPIB(YY)
      implicit real*8(a-h,o-z)
      COMMON/VAR/S,SS,PT
      COMMON/VA/T,U,GV,GW,AU,X3
      COMMON/BK1/B(2)
      common/quarkmass/qm,xmq,ims
      EXTERNAL PI2
      EXTERNAL ECHE
      EAX=ECHE(YY)
      EY=EXP(YY)
      T=-PT*SS/EY
      U=-PT*SS*EY

c   	teta = 2*atan(exp(-yy))
        amt = sqrt(qm**2 + pt**2)
        teta = atan(pt/amt/sinh(yy))
      ebeam = ss/2.
      pq = pt/sin(teta)
      eq = sqrt(pq**2 + qm**2)
      t = -2*ebeam*eq + 2*ebeam*pq*cos(teta) + qm**2
      u = -2*ebeam*eq - 2*ebeam*pq*cos(teta) + qm**2
      s = 4*ebeam**2

      TEST=S+T+U
      IF (TEST.GT..9) GO TO 1
      SPIB=0.
      GO TO 2
   1  GV=1.+T/S
      GW=-U/(T+S)
      AU=GV*GW
      B(1)=1.-GV+AU
      B(2)=1.
      SPIB=QUAX3(8,1,PI2)*AU**2
   2  CONTINUE
      RETURN
      END


      function phfragmell(qp2)
c returns the fragmentation function second moment
      implicit none
      real*8 qp2,phfragmell
      character*1 hvqs
      common/hvqtype/hvqs
      real*8 alam5qcd
      common/lamqcd/alam5qcd
      real*8 xmb,xmc
      common/hvqmass/xmb,xmc
      integer iloopas,iloopfr
      common/fonfrloop/iloopas,iloopfr
      real*8 ho,plla
      COMMON/FLAG/HO,PLLA
      real*8 cmu0
      common/cmu0/cmu0
      integer ifixed
      common/fixedorder/ifixed
      real*8 xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2
      integer nnf
      common/bndary/xxqsq0,xxqsqf,xxalf0,xxalff,xxsclf2,nnf
      real*8 rcn
      common/crcn/rcn
      complex * 16 ci
      integer k1,k2,init
      real*8 xpi
c local variables
      integer nf,lcoef,leadfl,j,icharm,icharmb
      real*8 q0,qf,sclf,xlam,xmq,xxf(14)
      complex * 16 cpsn(2,2),csn(2),cfl(20)
      complex * 16 cn,cpfm,cpfp,coefq,coefg
c functions
      real*8 alfas_p
c
      data ci/(0.d0,1.d0)/
      data k1,k2/2,1/,init/0/
      data xpi/3.141592653589793d0/
      xmq=xmc
      nnf = 4
      nf = nnf
      cn = rcn
      lcoef=0
      q0 = cmu0*xmq
      qf = sqrt(qp2)
      sclf = 1d0
      if(ifixed.eq.1)then
        q0=qf
      endif
      xlam = alam5qcd	! 'enter lambda_5'
c fill common block /bndary/
      xxsclf2 = sclf**2
      xxqsq0 = q0*q0
      xxqsqf = qf*qf*xxsclf2
      xxalf0 = alfas_p(xxqsq0,xlam,nnf)
      xxalff = alfas_p(xxqsqf,xlam,nnf)
c
      leadfl = iloopas		!'enter 1 for leading, 2 for next-to-leading'
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
      call parton(xxf,nf)
c the labeling convention for the flavour is given in the headers
c fonfra_old (b and bbar are 2 and 3, c and cbar are 4 and 5)
      icharm=4
      icharmb=5
      phfragmell = xxf(icharm)+xxf(icharmb)
      return
      end




      subroutine cachestore(iret,v1,v2,v3,v4,v5,v6,v7,v8,v9)
c saves v_i, i=1,..,9 up to 5 times
c if iret=0, values are already saved, and if iret=1 they are not
      implicit none
      integer iret
      real*8 v1,v2,v3,v4,v5,v6,v7,v8,v9
      integer n,m
      parameter (n=5,m=9)
      real*8 store(n,m)
      integer ilast,icurr,ii,jj
      save ilast,icurr
      if(ilast.lt.n) then
         ilast=ilast+1
         icurr=ilast
      else
         if(icurr.eq.n) then
            icurr=1
         else
            icurr=icurr+1
         endif
      endif
      store(icurr,1)=v1
      store(icurr,2)=v2
      store(icurr,3)=v3
      store(icurr,4)=v4
      store(icurr,5)=v5
      store(icurr,6)=v6
      store(icurr,7)=v7
      store(icurr,8)=v8
      store(icurr,9)=v9
      return
      entry cachelookup(iret,v1,v2,v3,v4,v5,v6,v7,v8,v9)
      iret=1
      do ii=1,ilast
         if( v1.eq.store(ii,1).and.
     #       v2.eq.store(ii,2).and.
     #       v3.eq.store(ii,3) )then
            iret=0
            jj=ii
         endif
      enddo
      if(iret.eq.0)then
         v4=store(jj,4)
         v5=store(jj,5)
         v6=store(jj,6)
         v7=store(jj,7)
         v8=store(jj,8)
         v9=store(jj,9)
      endif
      return
      entry cacheinit
      ilast=0
      end

