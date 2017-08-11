C****************************************************************************
C
C		 	EPPS16.f
C
C A fortran interface for the scale dependent nuclear modifications
C
C   		R_f^A(x,Q) = f_A(x,Q)/f_p(x,Q) 
C
C where f_A is the distribution of the parton flavour f for a PROTON in a
C nucleus A, and f_p is the corresponding parton distribution in the free proton.
C  
C When using this interface, please refer to:
C  
C K.J. Eskola, P. Paakkinen, H. Paukkunen and C.A. Salgado,
C "EPPS16: Nuclear parton distributions with LHC data."
C Published in EPJC
C Eprint: arXiv:1612.05741
C
C Questions & comments to:
C   hannu.paukkunen@jyu.fi
C   petja.paakkinen@jyu.fi
C   kari.eskola@jyu.fi
C   carlos.salgado@usc.es
C 
C ***************************************************************************
C Instructions:
C
C For given input values of
C
C     order (integer): dummy variable, retained for
C                      compatibility with the EPS09 routine
C
C
C     pset (integer):
C            1     = central fit
C            2,3   = error sets S{+1}, S{-1}
C            4,5   = error sets S{+2}, S{-2}
C            ...   ...
C            40,41 = error sets {S+20}, {S-20}
C
C     A (integer): atomic number
C     x (double precision) : Bjorken-x
C     Q (doubleprecision ) : scale in GeV
C
C the command
C
C   Call EPPS16(order, pset, A, x, Q, ruv, rdv, ru, rd, rs, rc, rb, rg)
C
C returns the bound proton nuclear corrections R_f^A(x,Q) (double precision)
C for
C	
C	ruv = up valence
C	rdv = down valence
C	ru  = up sea
C	rd  = down sea
C	rs  = strange
C	rc  = charm
C	rb  = bottom
C	rg  = gluons
C
C The nuclear corrections for bound neutrons can be obtained
C by the isospin symmetry, e.g. the total up quark distribution
C per nucleon in a nucleus A with Z protons is
C
C  u_A(x,Q) =    Z/A * [ruv*uV_p(x,Q) + ru*uSea_p(x,Q)] +
C            (A-Z)/A * [rdv*dV_p(x,Q) + rd*dSea_p(x,Q)]
C
C Note that the parametrization should only be applied at the
C kinematical domain
C
C             1e-7 <= x <= 1
C              1.3 <= Q <= 10000 GeV.
C
C Outside these boundaries the code stops and
C throws an error message.
C
C The data used by the program for required order
C and atomic number A, are stored in separate files
C
C   NLO: EPPS16NLOR_A
C
C which are assumed to be located in the working directory.
C
C The error bands for absolute cross-sections and for
C their nuclear ratios should be computed as explained
C in Secs. 2.5 and 4 of arXiv:0902.4154 [hep-ph]. For
C the absolute cross sections, both the errors in the
C free-proton PDFs f_p(x,Q) and the errors in
C the modifications R_f^A(x,Q) should be accounted for.
C For the nuclear ratios, it is usually sufficient to 
C account only for the errors in the modifications R_f^A(x,Q).
C
C *********************************************************
C *********************************************************/



      Subroutine EPPS16(order, pset, AAA, xxx, QQQ,
     &                   ruv, rdv, ru, rd, rs, rc, rb, rg)

      Implicit none
      Double precision :: ruv, rdv, ru, rd, rs, rc, rb, rg, QQQ, xxx
      Double precision :: LSTEP, x, Q, Q2, allvalues(1:41,1:8,0:50,0:80)
      Double precision :: x_i=0.0000001, arg(4), fu(4), res, fg(4)
      Double precision ::  result(9), dummy
      Double precision :: realQ, Q2min=1.69d0, Q2max=100000000.d0
      Double precision :: n_x, zero=0.d0, apu

      Character (Len=50) filenimi

      Integer :: Qsteps=30, r
      Integer :: xsteps=80, startline, lineno
      Integer :: k, p, t, Qpoint, xpoint, pset, iovar
      Integer :: Charmflag, Bottomflag
      Integer :: setnumber,j, A, openchannel, order, AAA
      Integer :: Alast = -10

      save Alast
      save allvalues

      Charmflag  = 0
      Bottomflag = 0

C *********************************************
C Stop if the set specifications are wrong ones
C *********************************************

      If (pset .LT. 1 .or. pset .GT. 41) then
      Write(*,*) 'In EPPS16: Wrong set!', pset
      Write(*,*) 'Central set: pset = 1'
      Write(*,*) 'Error sets : pset = 2...41'
      Stop
      End If

C ********************************
C Make sure not to change any
C specifications given by the user
C ********************************

      A  = AAA
      x  = xxx
      Q  = QQQ
      Q2 = Q*Q 

C *******************************
C Freeze x if it's < 10E-6 or > 1
C *******************************

      If (x .LT. x_i) Then
      Write(*,*) "In EPPS16: x<1E-7"
      stop
      End If
      If (x .GT. 1.d0) Then
      Write(*,*) "In EPPS16: x>1"
      stop
      End If

C ************************************
C Freeze Q^2 if it's < 1.69 or > 10E+6
C ************************************

      If (Q2 .LT. Q2min) Then
      Write(*,*) "In EPPS16: Q<1.3GeV"
      stop
      End If
      If (Q2 .GT. Q2max) Then
      Write(*,*) "In EPPS16: Q>10000GeV"
      stop
      End If

C If the set specifications have been changed, read the tables again

      If (A .NE. Alast) Then

C Read the table

        If (A < 10) Then
        Write(filenimi,'("EPPS16NLOR_", I1)'), A
        Else If (A < 100) Then
        Write(filenimi,'("EPPS16NLOR_", I2)'), A
        Else If (A < 1000) Then
        Write(filenimi,'("EPPS16NLOR_", I3)'), A
        End If

      Call NextUnitEPPS16(openchannel)

      OPEN (openchannel, file = filenimi, status='OLD', IOSTAT=iovar)

      If (iovar .NE. 0) Then
      Write(*,*) 'In EPPS16: Missing file: ',filenimi
      stop
      End If

      Do setnumber = 1, 41

      Do k = 0, Qsteps

      Read(openchannel,*) dummy

      Do t = 0,xsteps-1

      Read(openchannel,*) (allvalues(setnumber,p,k,t), p=1,8)

      End Do
      End Do

      End Do

      Close(openchannel)

      Alast     = A

      End If

C Find out the position in the loglog Q^2-grid

      realQ  = Qsteps * (log(log(Q2)/log(Q2min)))/
     &                  (log(log(Q2max)/log(Q2min)))
      Qpoint = Aint(realQ)

      If (Qpoint .LE. 0) Then
         Qpoint = 1
      End If
      If (Qpoint .GE. Qsteps-2) Then
         Qpoint = Qsteps-2
      End If

C Find the position in the x-grid

      LSTEP = (0.d0 - (Log(1.d0/x_i) + 2.d0*(1-x_i)))/(1.d0*xsteps)

      n_x  = ((Log(1.d0/x)+2.d0*(1-x))-(Log(1/x_i)+2.d0*(1-x_i)))/LSTEP
      xpoint = Aint(n_x)

      Do t=1,8

      If (t > 2 .and. t < 8) Then

       If (xpoint == 0) Then
       xpoint = 1
       Else If (xpoint > (xsteps-6)) Then
       xpoint = xsteps - 6
       End If

      Else

       If (xpoint == 0) Then
       xpoint = 1
       Else If (xpoint > (xsteps-4)) Then
       xpoint = xsteps - 4
       End If

      End If

C *********************
C Interpolate the grids 
C *********************

      Do k = 1, 4
      arg(k) = xpoint-2+k
      End Do

c For charm, don't use the point Q=1.3 for interpolation

      If (t == 6 .and. Qpoint == 1) Then
      Qpoint = 2
      Charmflag = 1
      End If

c For bottom, only use points Q>4.75 for interpolation

      If (t == 7 .and. Qpoint < 17 .and. Qpoint > 1) Then
      bottomflag = Qpoint;
      Qpoint = 17;
      End If

      Do j=1,4

      fu(1) = allvalues(pset,t,Qpoint-2+j,xpoint-1)
      fu(2) = allvalues(pset,t,Qpoint-2+j,xpoint)
      fu(3) = allvalues(pset,t,Qpoint-2+j,xpoint+1)
      fu(4) = allvalues(pset,t,Qpoint-2+j,xpoint+2)
      Call luoviEPPS16(fu,arg,4,n_x,res)
      fg(j) = res

      End Do

C *****************************************
C *****************************************

      arg(1) = Qpoint-1
      arg(2) = Qpoint
      arg(3) = Qpoint+1
      arg(4) = Qpoint+2

      Call luoviEPPS16(fg,arg,4,realQ,res)
  
      result(t) = res

c For charm, put back the original Qpoint

      If (Charmflag == 1) Then
      Qpoint = 1
      Charmflag = 0
      End If

c For bottom, put back the original Qpoint

      If (Bottomflag > 1) Then
      Qpoint = bottomflag
      bottomflag = 0
      End If

      End Do

      ruv = result(1)
      rdv = result(2)
      ru  = result(3)
      rd  = result(4)
      rs  = result(5)
      rc  = result(6)
      rb  = result(7)
      rg  = result(8)

c Put bottom to zero below the mass threshold

      If (Q < 4.75d0) Then
      rb  = 0.d0
      End If

200   Continue

      End Subroutine EPPS16

! ********************************
! Modified version of Cern Library
! interpolation routine E100
! ********************************

      SUBROUTINE luoviEPPS16(F,ARG,MMM,Z,SUM)

      Implicit none
      INTEGER :: MMM
      Double precision  :: F(MMM), ARG(MMM), COF(MMM), SUM, Z
      INTEGER :: M, MM, I, J, JNDEX, INDEX

      MM = MIN(MMM, 20)
      M = MM - 1
      DO 1780 I= 1, MM
      COF(I) = F(I)
 1780 Continue
      DO 1800 I= 1, M
      DO 1790 J= I, M
         JNDEX = MM - J
         INDEX = JNDEX + I
         COF(INDEX) = (COF(INDEX)-COF(INDEX-1))/(ARG(INDEX)-ARG(JNDEX))
 1790 CONTINUE
 1800 CONTINUE
      SUM = COF(MM)
      DO 1810 I= 1, M
         INDEX = MM - I
         SUM = (Z-ARG(INDEX))*SUM + COF(INDEX)
 1810 CONTINUE

      End SUBROUTINE luoviEPPS16

! **********************
! Find the open I/O unit
! **********************

      Subroutine NextUnitEPPS16(firstopen)

      Logical EX
      Integer firstopen, N

      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            firstopen = N
            Goto 20 
        Endif
10    Continue
      Stop ' There is no available I/O unit. '
20    Continue
      End Subroutine NextUnitEPPS16
