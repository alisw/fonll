C-CE PROGRAMME CONTIENT:
C----------- LES ELEMENTS DE MATRICE INTEGRES SUR L'ESPACE DE PHASE
C----------- LE CHANGEMENT DU SHEMA DE FACTORISATION
C----------- LES COMBINAISONS DE FONCTION DE STRUCTURE
C
C---CARREE DES ELEMENTS DE MATRICE ( INTEGRES SUR L'ESPACE DE PHASE)--
      DOUBLE PRECISION FUNCTION F0(V,X3)
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      COMMON/FCTD/GPPT,GPPC
      MU2=MU**2
      BX1=GV*GW/V/X3
      BX2=(1.-GV)/(1.-V)/X3
      SHD=BX1*BX2*GS
 
      IF (J0.EQ.1) THEN
      FBORN=GPPT*PI*CF/N/
     &SHD/V/(1.-V)*(V**2+1.)/(1.-V)**2
      FBORC=GPPC*PI*CF/N/
     &SHD/V/(1.-V)*(V**2-2.*V+2.)/V**2
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.2) THEN
      FBORN=1.E-30
      FBORC=1.E-30
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.3) THEN
      FBORN=GPPT*PI*CF/N/
     &SHD/V/(1.-V)*(V**2+1.)/(1.-V)**2
      FBORC=GPPC*PI*CF/N/
     &SHD/V/(1.-V)*(V**2-2.*V+2.)/V**2
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.4) THEN
      FBORN=1.E-30
      FBORC=1.E-30
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.5) THEN
      FBORN=GPPT*PI*CF/N/
     &SHD/V/(1.-V)*(2.*V**2-2.*V+1.)
      FBORC=GPPC*PI*CF/N/
     &SHD/V/(1.-V)*(2.*V**2-2.*V+1.)
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.6) THEN
      FBORN=GPPT*PI*2.*CF/(N**2)/
     &SHD/V/(1.-V)*(N*V**4-2.*N*V**3+4.*N*V**2+V**2-(3.*N+1.)*V+N)/
     &V**2/(1.-V)**2
      FBORC=GPPC*PI*2.*CF/(N**2)/
     &SHD/V/(1.-V)*(N*V**4-2.*N*V**3+4.*N*V**2+V**2-(3.*N+1.)*V+N)/
     &V**2/(1.-V)**2
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.7) THEN
      FBORN=1.E-30
      FBORC=1.E-30
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.8) THEN
      FBORN=1.E-30
      FBORC=1.E-30
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.9) THEN
      FBORN=1.E-30
      FBORC=1.E-30
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.10) THEN
      FBORN=1.E-30
      FBORC=1.E-30
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.11) THEN
      FBORN=GPPT*PI*2.*CF/(N**2)/
     &SHD/V/(1.-V)*(N*V**4-(3.*N+1.)*V**3+(4.*N+1.)*V**2-2.*N*V+N)/
     &(1.-V)**2
      FBORC=GPPC*PI*2.*CF/(N**2)/
     &SHD/V/(1.-V)*(N*V**4-(N-1.)*V**3+(N-2.)*V**2-(N-1.)*V+N)/
     &V**2
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.12) THEN
      FBORN=GPPT*PI*CF/(N**2)/
     &SHD/V/(1.-V)*(2.*V**2-2.*V+1.)*(2.*N**2*V**2-2.*N**2*V+N**2-1.)/
     &V/(1.-V)
      FBORC=GPPC*PI*CF/(N**2)/
     &SHD/V/(1.-V)*(2.*V**2-2.*V+1.)*(2.*N**2*V**2-2.*N**2*V+N**2-1.)/
     &V/(1.-V)
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.13) THEN
      FBORN=GPPT*PI/(2.*N**2)/
     &SHD/V/(1.-V)*(V**2+1.)*((N**2-1.)*V**2+2.*V+(N**2-1.))/
     &V/(1.-V)**2
      FBORC=GPPC*PI/(2.*N**2)/
     &SHD/V/(1.-V)*(V**2-2.*V+2.)*((N**2-1.)*V**2-2.*N**2*V+2.*N**2)/
     &V**2/(1.-V)
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.14) THEN
      FBORN=GPPT*PI/(2.*N**2)/
     &SHD/V/(1.-V)*(V**2-2.*V+2.)*((N**2-1.)*V**2-2.*N**2*V+2.*N**2)/
     &V**2/(1.-V)
      FBORC=GPPC*PI/(2.*N**2)/
     &SHD/V/(1.-V)*(V**2+1.)*((N**2-1.)*V**2+2.*V+(N**2-1.))/
     &V/(1.-V)**2
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.15) THEN
      FBORN=GPPT*PI*(4.*N**2)/VC/
     1SHD/V/(1.-V)*(3.-V*(1.-V)+V/(1.-V)**2+(1.-V)/V**2)
      FBORC=GPPC*PI*(4.*N**2)/VC/
     1SHD/V/(1.-V)*(3.-V*(1.-V)+V/(1.-V)**2+(1.-V)/V**2)
      F0=FBORN+FBORC
 
      ELSE IF (J0.EQ.16) THEN
      FBORN=GPPT*PI/(2.*N)/VC/
     1SHD/V/(1.-V)*(V**2+(1.-V)**2)*(2.*N**2*(V**2-V)+N**2-1.)/V/(1.-V)
      FBORC=GPPC*PI/(2.*N)/VC/
     1SHD/V/(1.-V)*(V**2+(1.-V)**2)*(2.*N**2*(V**2-V)+N**2-1.)/V/(1.-V)
      F0=FBORN+FBORC
 
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION AVWPL(W,V,S)
C     TERME EN (1-W)+
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      WPLUS=1.
C TERME EN WPLUS
 
      IF (J0.EQ.1) THEN
      AVWPL =(-2*LOG(1-V)*(V-1)*(V**2+1)*(4*V2+V1)-2*(V-1)*(V**2+1)*LOG(
     1   V)*(4*V2-5*V1)+2*LOG(S/MP**2)*(V-1)*(V**2+1)*V1+4*LOG(S/M**2)*(
     2   V-1)*(V**2+1)*V1+3*(V-1)*(V**2+1)*V1)*WPLUS/((V-1)**3*V)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=A0(V,S)/V/W*HQQW(W,V,1)+A0(V*W,S)/(1.-V)*HQQW(W,V,2)
     2+E0(W*V,S)*HGQW(W,V,2)/(1.-V)*N/VC
     2+A0(VZ,S)*HFQQW(W,V)/(1.-V+V*W)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.2) THEN
      AVWPL =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=(E0(V,S)/V/W*HGQW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQW(W,V,2))
     2*N/VC+HFGQW(W,V)/(1.-V+V*W)*(A0(VZ,S)+A0(1.-VZ,S))
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.3) THEN
      AVWPL =(2*(V-1)*(V**2+1)*LOG(V)*(8*V2+V1)-2*LOG(1-V)*(V-1)*(V**2+1
     1   )*(4*V2+V1)+2*LOG(S/MP**2)*(V-1)*(V**2+1)*V1+4*LOG(S/M**2)*(V-1
     2   )*(V**2+1)*V1+3*(V-1)*(V**2+1)*V1)*WPLUS/((V-1)**3*V)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=A0(V,S)/V/W*HQQW(W,V,1)+A0(V*W,S)/(1.-V)*HQQW(W,V,2)
     2+E0(W*V,S)*HGQW(W,V,2)/(1.-V)*N/VC
     2+A0(VZ,S)*HFQQW(W,V)/(1.-V+V*W)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.4) THEN
      AVWPL =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=(E0(V,S)/V/W*HGQW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQW(W,V,2))
     2*N/VC+HFGQW(W,V)/(1.-V+V*W)*(A0(VZ,S)+A0(1.-VZ,S))
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.5) THEN
      AVWPL =(2*(V-1)*(2*V**2-2*V+1)*LOG(V)*(8*V2+V1)-2*LOG(1-V)*(V-1)*(
     1   2*V**2-2*V+1)*(4*V2-3*V1)+2*LOG(S/MP**2)*(V-1)*(2*V**2-2*V+1)*V
     2   1+4*LOG(S/M**2)*(V-1)*(2*V**2-2*V+1)*V1+3*(V-1)*(2*V**2-2*V+1)*
     3   V1)*WPLUS/((V-1)*V)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=A2(1.-V,S)/V/W*HQQW(W,V,1)+A2(1.-V*W,S)/(1.-V)*HQQW(W,V,2)
     2+A2(1.-VZ,S)*HFQQW(W,V)/(1.-V+V*W)
     2+D1(VZ,S)*HFQGW(W,V)/(1.-V+V*W)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.6) THEN
      AVWPL =(4*(V-1)*LOG(V)*(6*(V-1)*V*V4-4*(V-1)*V*V3-4*(V**4-2*V**3+4
     1   *V**2-3*V+1)*V2+(3*V**4-2*V**3+6*V**2-3*V+1)*V1)+4*LOG(1-V)*(V-
     2   1)*(2*(V-1)*V*V4-4*(V-1)*V*V3-4*(V**4-2*V**3+4*V**2-3*V+1)*V2+(
     3   V**4-6*V**3+10*V**2-9*V+3)*V1)+4*LOG(S/MP**2)*(V-1)*(2*(V-1)*V*
     4   V4+(V**4-2*V**3+4*V**2-3*V+1)*V1)+8*LOG(S/M**2)*(V-1)*(2*(V-1)*
     5   V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1)+6*(V-1)*(2*(V-1)*V*V4+(V**4
     6   -2*V**3+4*V**2-3*V+1)*V1))*WPLUS/((V-1)**3*V**3)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=B0(V,S)/V/W*HQQW(W,V,1)+B0(V*W,S)/(1.-V)*HQQW(W,V,2)
     C+E0(1.-V,S)/W/V*N/VC*HGQW(W,V,1)
     1+E0(W*V,S)/(1.-V)*N/VC*HGQW(W,V,2)
     2+B0(VZ,S)*HFQQW(W,V)/(1.-V+V*W)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.7) THEN
      AVWPL =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=(E0(V,S)/V/W*HGQW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQW(W,V,2))
     2*N/VC+HFGQW(W,V)/(1.-V+V*W)*B0(VZ,S)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.8) THEN
      AVWPL =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=N/VC*D1(V,S)/V/W*HGQW(W,V,1)+A0(1.-V*W,S)/(1.-V)*HQGW(W,V,2)
     2/N*VC+HFQGW(W,V)/(1.-V+V*W)*E0(1.-VZ,S)+VC/N*HQGW(W,V,2)/(1.-V)*
     2A2(1.-V*W,S)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.9) THEN
      AVWPL =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
      ZWPL=N/VC*D1(VY,S)/V/W*HGQW(W,V,1)+A0(1.-V*W,S)/(1.-V)*HQGW(W,V,2)
     2/N*VC+HFQGW(W,V)/(1.-V+V*W)*E0(1.-VZ,S)+VC/N*HQGW(W,V,2)/(1.-V)*
     2A2(V*W,S)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.10) THEN
      AVWPL =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
      ZWPL=N/VC*D1(VY,S)/V/W*HGQW(W,V,1)+D0(1.-V*W,S)/(1.-V)*HQGW(W,V,2)
     2/N*VC+HFQGW(W,V)/(1.-V+V*W)*E0(1.-VZ,S)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.11) THEN
      AVWPL =-2*(2*LOG(V)*(2*(V-1)*V**2*V4+8*(V-1)*V**2*V3-8*(V**4-3*V**
     1   3+4*V**2-2*V+1)*V2+(-V**4+3*V**3-4*V**2+2*V-1)*V1)+2*LOG(1-V)*(
     2   2*(V-1)*V**2*V4-4*(V-1)*V**2*V3+4*(V**4-3*V**3+4*V**2-2*V+1)*V2
     3   +(-3*V**4+9*V**3-10*V**2+6*V-1)*V1)+2*LOG(S/MP**2)*(2*(V-1)*V**
     4   2*V4+(-V**4+3*V**3-4*V**2+2*V-1)*V1)+4*LOG(S/M**2)*(2*(V-1)*V**
     5   2*V4+(-V**4+3*V**3-4*V**2+2*V-1)*V1)+3*(2*(V-1)*V**2*V4+(-V**4+
     6   3*V**3-4*V**2+2*V-1)*V1))*WPLUS/((V-1)**2*V)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=D0(V,S)/V/W*HQQW(W,V,1)+D0(V*W,S)/(1.-V)*HQQW(W,V,2)
     C+E0(W*V,S)/(1.-V)*N/VC*HGQW(W,V,2)
     1+D0(VZ,S)*HFQQW(W,V)
     2/(1.-V+V*W)+D1(VZ,S)*HFQGW(W,V)/(1.-V+V*W)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.12) THEN
      AVWPL =(LOG(V)*(6*N**4*(V-1)*(2*V**2-2*V+1)*(2*CQ*(2*V**2-2*V+1)-1
     1   0*V**2+14*V-7)*VC+12*N**2*(V-1)*(V**2-V-CQ+4)*(2*V**2-2*V+1)*VC
     2   -6*(V-1)*(2*V**2-2*V+1)*VC)+LOG(1-V)*(-6*N**4*(V-1)*(2*V**2-2*V
     3   +1)*(2*V**2+2*V-1)*VC-12*N**2*(V-1)*(V**2-V-1)*(2*V**2-2*V+1)*V
     4   C+6*(V-1)*(2*V**2-2*V+1)*VC)+LOG(S/M**2)*(-12*N**4*(V-1)*(2*V**
     5   2-2*V+1)**2*VC+24*N**2*(V-1)*(V**2-V+1)*(2*V**2-2*V+1)*VC-12*(V
     6   -1)*(2*V**2-2*V+1)*VC)+LOG(S/MP**2)*(12*N**2*(V-1)*(2*V**2-2*V+
     7   1)*VC-12*N**4*(V-1)*(2*V**2-2*V+1)**2*VC)+N**2*(V-1)*(2*V**2-2*
     8   V+1)*(18*V**2-18*V+7)*VC+2*N**4*(V-1)*(2*V**2-2*V+1)**2*VC-4*GT
     9   R*N**3*(V-1)*(2*V**2-2*V+1)**2*VC+4*GTR*N*(V-1)*(2*V**2-2*V+1)*
     :   VC-9*(V-1)*(2*V**2-2*V+1)*VC)*WPLUS/(N**2*(V-1)**2*V**2)/3.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=D1(V,S)/V/W*HQQW(W,V,1)+D1(V*W,S)/(1.-V)*HQQW(W,V,2)
     C+2.*(2.*GTR-1.)*A2(VZ,S)/(1.-V+V*W)*HFGQW(W,V)
     C+E0(1.-V*W,S)/(1.-V)*N/VC*HGQW(W,V,2)
     1+N/VC*E0(V,S)/W/V*HGQW(W,V,1)+D1(VZ,S)*HFGGW(W,V)
     2/(1.-V+V*W)+(D0(1.-VZ,S)+D0(VZ,S))*HFGQW(W,V)/(1.-V+V*W)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.13) THEN
      AVWPL = (LOG(1-V)*(-12*N**2*(V-1)*(V**2+1)*(CQ*(V-1)**2-2*(
     1   V**2-3*V+1))*VC+12*(CQ-1)*N**4*(V-1)*(V**2+1)**2*VC+12*(
     2   V-1)**3*(V**2+1)*VC)+LOG(V)*(-6*N**4*(V-1)*(V**2+1)*(2*
     3   CQ*(V**2+1)-3*V**2-7)*VC+12*N**2*(V-1)*(-4*V**2+7*V+CQ*(
     4   V-1)**2-4)*(V**2+1)*VC+6*(V-1)**3*(V**2+1)*VC)+LOG(S/M*
     5   *2)*(-12*N**2*(V-1)*(V**2+1)*(2*V**2-3*V+2)*VC+18*N**4*(
     6   V-1)*(V**2+1)**2*VC+6*(V-1)**3*(V**2+1)*VC)+LOG(S/MP**2
     7   )*(-12*N**2*(V-1)*(V**2+1)*(V**2-V+1)*VC+6*N**4*(V-1)*(V
     8   **2+1)**2*VC+6*(V-1)**3*(V**2+1)*VC)-N**2*(V-1)*(V**2+1)
     9   *(7*V**2+4*V+7)*VC-2*N**4*(V-1)*(V**2+1)**2*VC+4*GTR*N**
     :   3*(V-1)*(V**2+1)**2*VC-4*GTR*N*(V-1)**3*(V**2+1)*VC+9*(V
     ;   -1)**3*(V**2+1)*VC)*WPLUS/(N**2*(V-1)**3*V**2)/3.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=E0(V,S)/V/W*HQQW(W,V,1)+E0(V*W,S)/(1.-V)*HGGW(W,V,2)
     C+2.*(2.*GTR-1.)*VC/N*A0(V*W,S)/(1.-V)*HQGW(W,V,2)
     C+(B0(V*W,S)+D0(V*W,S))/(1.-V)*VC/N*HQGW(W,V,2)
     1+N/VC*D1(V,S)/W/V*HGQW(W,V,1)+E0(VZ,S)*HFQQW(W,V)
     2/(1.-V+V*W)+E0(1.-VZ,S)*HFQGW(W,V)/(1.-V+V*W)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.14) THEN
      AVWPL = (LOG(1-V)*(4*CQ*N**2*(V-1)**2*V**2*(V**2-2*V+2)*VC-
     1   4*N**4*(V-1)**2*(V**2-2*V+2)*(CQ*(V**2-2*V+2)-2*(V-1)**2
     2   )*VC)+LOG(V)*(-4*N**2*(V-1)**2*(V**2-2*V+2)*(2*CQ*V**2-
     3   3*V**2+V-1)*VC+2*(4*CQ-9)*N**4*(V-1)**2*(V**2-2*V+2)**2*
     4   VC-2*(V-1)**2*V**2*(V**2-2*V+2)*VC)+LOG(S/M**2)*(4*N**2
     5   *(V-1)**2*(V**2-2*V+2)*(2*V**2-V+1)*VC-6*N**4*(V-1)**2*(
     6   V**2-2*V+2)**2*VC-2*(V-1)**2*V**2*(V**2-2*V+2)*VC)+LOG(
     7   S/MP**2)*(4*N**2*(V-1)**2*V**2*(V**2-2*V+2)*VC-4*N**4*(V
     8   -1)**2*(V**2-2*V+2)**2*VC))*WPLUS/(N**2*(V-1)**3*V**3)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=E0(1.-V,S)/V/W*HQQW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGGW(W,V,2)
     C+N/VC*F2(V,S)/W/V*
     1HGQW(W,V,1)+VC/N*D1(V*W,S)/(1.-V)*HQGW(W,V,2)+E0(VZ,S)*HFGQW(W,V)
     2/(1.-V+V*W)+E0(1.-VZ,S)*HFGGW(W,V)/(1.-V+V*W)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.15) THEN
      AVWPL = (-96*N**3*(V-1)*(V**2-V+1)**2*(2*CQ*(V**2-V+1)-4*V**
     1   2+5*V-5)*LOG(V)*VC+96*N**3*LOG(1-V)*(V-1)*(V**2-V+1)**
     2   2*(CQ*(V**2-V+1)-(V-1)**2)*VC+96*N**3*LOG(S/MP**2)*(V-1
     3   )*(V**2-V+1)**3*VC+192*N**3*LOG(S/M**2)*(V-1)*(V**2-V+1
     4   )**3*VC-88*N**3*(V-1)*(V**2-V+1)**3*VC+32*GTR*N**2*(V-1)
     5   *(V**2-V+1)**3*VC)*WPLUS/((V-1)**3*V**3)/3.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=F2(V,S)/V/W*HGGW(W,V,1)+F2(V*W,S)/(1.-V)*HGGW(W,V,2)
     C+4.*GTR*VC/N*(E0(VY,S)/W
     1/V*HQGW(W,V,1)+E0(V*W,S)/(1.-V)*HQGW(W,V,2))+F2(VZ,S)*HFGGW(W,V)
     2/(1.-V+V*W)+4.*GTR*D1(VZ,S)*HFGQW(W,V)/(1.-V+V*W)
      AVWPL=AVWPL+ZWPL
 
      ELSE IF (J0.EQ.16) THEN
      AVWPL = (LOG(V)*(4*N**4*(V-1)**2*(2*V**2-2*V+1)*(CQ*(2*V**2
     1   -2*V+1)-2*(3*V**2-4*V+2))*VC-4*(CQ-4)*N**2*(V-1)**2*(2*V
     2   **2-2*V+1)*VC)+LOG(1-V)*(4*CQ*N**2*(V-1)**2*(2*V**2-2*V
     3   +1)*VC-4*N**4*(V-1)**2*(2*V**2-2*V+1)*(CQ*(2*V**2-2*V+1)
     4   -2*(V-1)**2)*VC)+LOG(S/MP**2)*(-2*N**4*(V-1)**2*(2*V**2
     5   -2*V+1)**2*VC+4*N**2*(V-1)**2*(V**2-V+1)*(2*V**2-2*V+1)*
     6   VC-2*(V-1)**2*(2*V**2-2*V+1)*VC)+LOG(S/M**2)*(8*N**2*(V
     7   -1)**2*(2*V**2-2*V+1)*VC-8*N**4*(V-1)**2*(2*V**2-2*V+1)*
     8   *2*VC))*WPLUS/(N**2*(V-1)**3*V**2)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZWPL=D1(V,S)/V/W*HGGW(W,V,1)+D1(V*W,S)/(1.-V)*HGGW(W,V,2)
     C+VC/N*(E0(V,S)/W/V
     1*HQGW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HQGW(W,V,2))+D1(VZ,S)*HFQQW(W,V)
     2/(1.-V+V*W)+F2(VZ,S)*HFQGW(W,V)/(1.-V+V*W)
      AVWPL=ZWPL+AVWPL
 
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION AVDEL(V,S)
C     TERME EN DELTA(1-W)
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      DELW=1.
C    TERME EN DELTA
 
      IF (J0.EQ.1) THEN
      AVDEL =DELW*(3*LOG(1-V)*(V-1)**4*(V**2+1)*(16*GTR*(V**2+1)*V4+16*G
     1   TR*(V**2+1)*V3-4*(17*V**2-3*V+8)*V2+(-17*V**2+12*V-29)*V1)+12*L
     2   OG(S/MU**2)*(V-1)**4*(V**2+1)**2*(4*GTR*V4+4*GTR*V3-11*V2-11*V1
     3   )-80*GTR*(V-1)**4*(V**2+1)**2*(V4+V3)+2*(V-1)**4*(V**2+1)*((9*
     4   PI**2*V**2+170*V**2+27* PI**2+170)*V2+(24* PI**2*V**2+179*V**2+
     5   6* PI**2+179)*V1)+36*LOG(1-V)*(V-1)**4*(V**2+1)*LOG(V)*((3*V**2
     6   +1)*V2+(-5*V**2-3)*V1)+18*LOG(1-V)**2*(V-1)**4*(V**2+1)*((V**2+
     7   7)*V2+(3*V**2+1)*V1)+9*(V-1)**4*(V**2+1)*LOG(V)*(4*(V-1)*V2+(3*
     8   V**2-4*V+7)*V1)-18*(V-1)**4*(V**2+1)*(9*V**2+7)*LOG(V)**2*(V2-V
     9   1)+LOG(S/M**2)*(36*(V-1)**4*(V**2+1)**2*LOG(V)*V1-36*LOG(1-V)*(
     :   V-1)**4*(V**2+1)**2*V1+54*(V-1)**4*(V**2+1)**2*V1)+LOG(S/MP**2)
     ;   *(36*(V-1)**4*(V**2+1)**2*LOG(V)*V1+27*(V-1)**4*(V**2+1)**2*V1)
     <   )/((V-1)**6*V*(V**2+1))/18.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=A0(V,S)/V*HQQD(V,1)+A0(V,S)/(1.-V)*HQQD(V,2)
     2+E0(V,S)*HGQD(V,2)/(1.-V)*N/VC
     2+A0(V,S)*HFQQD(V)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.2) THEN
      AVDEL =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=(E0(V,S)/V*HGQD(V,1)+E0(1.-V,S)/(1.-V)*HGQD(V,2))
     2*N/VC+HFGQD(V)*(A0(V,S)+A0(1.-V,S))
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.3) THEN
      AVDEL =DELW*(3*LOG(1-V)*(V-1)**4*(V**2+1)*(16*GTR*(V**2+1)*V4+16*G
     1   TR*(V**2+1)*V3-4*(8*V**2-3*V+17)*V2+(-29*V**2+12*V-17)*V1)+12*L
     2   OG(S/MU**2)*(V-1)**4*(V**2+1)**2*(4*GTR*V4+4*GTR*V3-11*V2-11*V1
     3   )-80*GTR*(V-1)**4*(V**2+1)**2*(V4+V3)+2*(V-1)**4*(V**2+1)*(2*(1
     4   8* PI**2*V**2+85*V**2+85)*V2+(15* PI**2+179)*(V**2+1)*V1)+18*LO
     5   G(1-V)**2*(V-1)**4*(V**2+1)*((7*V**2+1)*V2+(V**2+3)*V1)-9*(V-1)
     6   **4*(V**2+1)*LOG(V)*(8*(V-1)*V2-3*(V**2+1)*V1)+36*(V-1)**4*(V**
     7   2+1)*(9*V**2+7)*LOG(V)**2*V2-144*LOG(1-V)*(V-1)**4*(V**2+1)*(3*
     8   V**2+2)*LOG(V)*V2+LOG(S/M**2)*(36*(V-1)**4*(V**2+1)**2*LOG(V)*V
     9   1-36*LOG(1-V)*(V-1)**4*(V**2+1)**2*V1+54*(V-1)**4*(V**2+1)**2*V
     :   1)+LOG(S/MP**2)*(36*(V-1)**4*(V**2+1)**2*LOG(V)*V1+27*(V-1)**4*
     ;   (V**2+1)**2*V1))/((V-1)**6*V*(V**2+1))/18.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=A0(V,S)/V*HQQD(V,1)+A0(V,S)/(1.-V)*HQQD(V,2)
     2+E0(V,S)*HGQD(V,2)/(1.-V)*N/VC
     2+A0(V,S)*HFQQD(V)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.4) THEN
      AVDEL =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=(E0(V,S)/V*HGQD(V,1)+E0(1.-V,S)/(1.-V)*HGQD(V,2))
     2*N/VC+HFGQD(V)*(A0(V,S)+A0(1.-V,S))
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.5) THEN
      AVDEL =DELW*(12*LOG(S/MU**2)*(V-1)**3*(2*V**2-2*V+1)**2*(4*GTR*V4+
     1   4*GTR*V3-11*V2-11*V1)-80*GTR*(V-1)**3*(2*V**2-2*V+1)**2*(V4+V3)
     2   -9*LOG(1-V)*(V-1)**3*(2*V**2-2*V+1)*(4*V*V2+(6*V**2-10*V+3)*V1)
     3   -9*(V-1)**3*(2*V**2-2*V+1)*LOG(V)*(8*(V-1)*V2-3*(2*V**2-2*V+1)*
     4   V1)-2*(V-1)**3*(2*V**2-2*V+1)**2*(2*(9* PI**2-85)*V2+(3* PI**2-
     5   179)*V1)-72*LOG(1-V)*(V-1)**3*(2*V**2-2*V+1)**2*LOG(V)*(3*V2-V1
     6   )+18*LOG(1-V)**2*(V-1)**3*(2*V-1)*(2*V**2-2*V+1)*(V2-V1)+36*(V-
     7   1)**3*(2*V**2-2*V+1)*(16*V**2-14*V+7)*LOG(V)**2*V2+LOG(S/M**2)*
     8   (36*(V-1)**3*(2*V**2-2*V+1)**2*LOG(V)*V1-36*LOG(1-V)*(V-1)**3*(
     9   2*V**2-2*V+1)**2*V1+54*(V-1)**3*(2*V**2-2*V+1)**2*V1)+LOG(S/MP*
     :   *2)*(36*(V-1)**3*(2*V**2-2*V+1)**2*LOG(V)*V1+27*(V-1)**3*(2*V**
     ;   2-2*V+1)**2*V1))/((V-1)**3*V*(2*V**2-2*V+1))/18.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=A2(1.-V,S)/V*HQQD(V,1)+A2(1.-V,S)/(1.-V)*HQQD(V,2)
     2+A2(1.-V,S)*HFQQD(V)
     2+D1(V,S)*HFQGD(V)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.6) THEN
      AVDEL =DELW*(LOG(S/M**2)*(36*(V-1)**4*(V**2+1)*(V**2-2*V+2)*LOG(V)
     1   *(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1)-36*LOG(1-V)*(V-1)
     2   **4*(V**2+1)*(V**2-2*V+2)*(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V
     3   +1)*V1)+54*(V-1)**4*(V**2+1)*(V**2-2*V+2)*(2*(V-1)*V*V4+(V**4-2
     4   *V**3+4*V**2-3*V+1)*V1))+LOG(S/MP**2)*(36*(V-1)**4*(V**2+1)*(V*
     5   *2-2*V+2)*LOG(V)*(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1)+2
     6   7*(V-1)**4*(V**2+1)*(V**2-2*V+2)*(2*(V-1)*V*V4+(V**4-2*V**3+4*V
     7   **2-3*V+1)*V1))+12*LOG(S/MU**2)*(V-1)**4*(V**2+1)*(V**2-2*V+2)*
     8   (4*GTR*(V**4-2*V**3+4*V**2-3*V+1)*V4-11*(V-1)*V*V4+4*GTR*(V**4-
     9   2*V**3+4*V**2-3*V+1)*V3-11*(V-1)*V*V3-11*(V**4-2*V**3+4*V**2-3*
     :   V+1)*V2+4*GTR*(V-1)*V*V2-11*(V**4-2*V**3+4*V**2-3*V+1)*V1)-80*G
     ;   TR*(V-1)**4*(V**2+1)*(V**2-2*V+2)*((V**4-2*V**3+4*V**2-3*V+1)*V
     <   4+(V**4-2*V**3+4*V**2-3*V+1)*V3+(V-1)*V*V2)+(V-1)**4*(V**2+1)*(
     =   V**2-2*V+2)*((V-1)*V*(18* PI**2*V**2-18* PI**2*V+15* PI**2+376)
     >   *V4-(V-1)*V*(18* PI**2*V**2-18* PI**2*V-45* PI**2-340)*V3+2*(9*
     ?    PI**2*V**4+170*V**4-18* PI**2*V**3-340*V**3+54* PI**2*V**2+680
     @   *V**2-45* PI**2*V-510*V+18* PI**2+170)*V2+2*(24* PI**2*V**4+179
     1   *V**4-48* PI**2*V**3-358*V**3+78* PI**2*V**2+716*V**2-54* PI**2
     2   *V-537*V+15* PI**2+179)*V1)+9*(V-1)**4*(V**2+1)*(V**2-2*V+2)*LO
     3   G(V)**2*((V-1)*V*(2*V**2-2*V+15)*V4-(V-1)*V*(2*V**2-2*V+11)*V3-
     4   2*(8*V**4-14*V**3+25*V**2-15*V+4)*V2+2*(6*V**4-6*V**3+13*V**2-7
     5   *V+2)*V1)-18*LOG(1-V)*(V-1)**4*(V**2+1)*(V**2-2*V+2)*LOG(V)*((V
     6   -1)*V*(2*V**2-2*V+3)*V4-(V-1)*V*(2*V**2-2*V+3)*V3-2*(V**2-V+2)*
     7   (3*V**2-3*V+1)*V2+2*V*(3*V**3-2*V**2+4*V-1)*V1)+9*LOG(1-V)**2*(
     8   V-1)**4*V*(V**2+1)*(V**2-2*V+2)*((V-1)*(2*V**2-2*V-1)*V4-(V-1)*
     9   (2*V**2-2*V-5)*V3+2*(2*V**2+V+1)*V2+2*(2*V**3-2*V**2+3*V-1)*V1)
     :   +3*(V-1)**4*(V**2+1)*(V**2-2*V+2)*LOG(V)*(8*GTR*(V-1)**2*(V**2-
     ;   2*V+2)*V4-2*(V-1)*V*(6*V-19)*V4+8*GTR*(V-1)**2*(V**2-2*V+2)*V3+
     <   4*(V-1)*V*(3*V-7)*V3-2*(V-1)*(17*V**3-51*V**2+53*V-22)*V2+8*GTR
     =   *(V-1)*V*V2+(5*V**4-14*V**3+26*V**2-9*V+1)*V1)+3*LOG(1-V)*(V-1)
     >   **4*(V**2+1)*(V**2-2*V+2)*(8*GTR*V**2*(V**2+1)*V4+2*(V-1)*V*(6*
     ?   V-5)*V4+8*GTR*V**2*(V**2+1)*V3-4*(V-1)*V*(3*V+4)*V3-2*V*(17*V**
     @   3+2*V+3)*V2+8*GTR*(V-1)*V*V2+(-13*V**4+30*V**3-58*V**2+33*V-9)*
     1   V1))/((V-1)**6*V**3*(V**2+1)*(V**2-2*V+2))/9.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=B0(V,S)/V*HQQD(V,1)+B0(V,S)/(1.-V)*HQQD(V,2)
     C+E0(1.-V,S)/V*N/VC*HGQD(V,1)
     1+E0(V,S)/(1.-V)*N/VC*HGQD(V,2)
     2+B0(V,S)*HFQQD(V)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.7) THEN
      AVDEL =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=(E0(V,S)/V*HGQD(V,1)+E0(1.-V,S)/(1.-V)*HGQD(V,2))
     2*N/VC+HFGQD(V)*B0(V,S)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.8) THEN
      AVDEL = 1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=N/VC*D1(V,S)/V*HGQD(V,1)+A0(1.-V,S)/(1.-V)*HQGD(V,2)
     2/N*VC+HFQGD(V)*E0(1.-V,S)+VC/N*HQGD(V,2)/(1.-V)*
     2A2(1.-V,S)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.9) THEN
      AVDEL = 1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ZDEL=N/VC*D1(1.-V,S)/V*HGQD(V,1)+A0(1.-V,S)/(1.-V)*HQGD(V,2)
     2/N*VC+HFQGD(V)*E0(1.-V,S)+VC/N*HQGD(V,2)/(1.-V)*
     2A2(V,S)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.10) THEN
      AVDEL = 1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ZDEL=N/VC*D1(1.-V,S)/V*HGQD(V,1)+D0(1.-V,S)/(1.-V)*HQGD(V,2)
     2/N*VC+HFQGD(V)*E0(1.-V,S)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.11) THEN
      AVDEL =-DELW*(2*(40*GTR*(V**4-3*V**3+4*V**2-2*V+1)*V4+4*(3* PI**2+
     1   47)*(V-1)*V**2*V4+40*GTR*(V**4-3*V**3+4*V**2-2*V+1)*V3+170*(V-1
     2   )*V**2*V3+(18* PI**2*V**4-170*V**4-54* PI**2*V**3+510*V**3+45*
     3   PI**2*V**2-680*V**2-36* PI**2*V+340*V+9* PI**2-170)*V2-40*GTR*(
     4   V-1)*V**2*V2+(3* PI**2*V**4-179*V**4-9* PI**2*V**3+537*V**3+3*
     5   PI**2*V**2-716*V**2-6* PI**2*V+358*V-6* PI**2-179)*V1)-12*LOG(S
     6   /MU**2)*(4*GTR*(V**4-3*V**3+4*V**2-2*V+1)*V4+11*(V-1)*V**2*V4+4
     7   *GTR*(V**4-3*V**3+4*V**2-2*V+1)*V3+11*(V-1)*V**2*V3-11*(V**4-3*
     8   V**3+4*V**2-2*V+1)*V2-4*GTR*(V-1)*V**2*V2-11*(V**4-3*V**3+4*V**
     9   2-2*V+1)*V1)-9*LOG(1-V)**2*((V-1)*(V**2+2*V-2)*V4-(V-1)*(5*V**2
     :   +2*V-2)*V3+2*V*(V**2+V+2)*V2-2*(V**3-3*V**2+2*V-2)*V1)-3*LOG(1-
     ;   V)*(8*GTR*(V**2+1)*V4+2*(V-1)*V*(5*V-6)*V4+8*GTR*(V**2+1)*V3+4*
     <   (V-1)*V*(4*V+3)*V3-2*(3*V**3+2*V**2+17)*V2-8*GTR*(V-1)*V**2*V2+
     =   (-9*V**4+33*V**3-58*V**2+30*V-13)*V1)+9*LOG(V)*(6*(V-1)*V**2*V4
     >   +4*(V-1)*(V**2-2*V+2)*V2-3*(V**4-3*V**3+4*V**2-2*V+1)*V1)+36*LO
     ?   G(1-V)*LOG(V)*(2*(V-1)*V**2*V4-8*(V-1)*V**2*V3+(6*V**4-18*V**3+
     @   27*V**2-12*V+7)*V2-(V-1)**2*(2*V**2-2*V+1)*V1)+9*LOG(S/MP**2)*(
     1   4*LOG(V)+3)*(2*(V-1)*V**2*V4+(-V**4+3*V**3-4*V**2+2*V-1)*V1)+18
     2   *LOG(S/M**2)*(2*LOG(V)-2*LOG(1-V)+3)*(2*(V-1)*V**2*V4+(-V**4+3*
     3   V**3-4*V**2+2*V-1)*V1)+36*LOG(V)**2*(8*(V-1)*V**2*V3+(-8*V**4+2
     4   3*V**3-30*V**2+14*V-7)*V2))/((V-1)**2*V)/9.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=D0(V,S)/V*HQQD(V,1)+D0(V,S)/(1.-V)*HQQD(V,2)
     C+E0(V,S)/(1.-V)*N/VC*HGQD(V,2)
     1+D0(V,S)*HFQQD(V)
     2+D1(V,S)*HFQGD(V)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.12) THEN
      AVDEL =DELW*(LOG(S/M**2)*(LOG(1-V)*(36*N**4*(V-1)**4*(2*V**2-2*V+1
     1   )**2*VC-72*N**2*(V-1)**4*(V**2-V+1)*(2*V**2-2*V+1)*VC+36*(V-1)*
     2   *4*(2*V**2-2*V+1)*VC)+LOG(V)*(-36*N**4*(V-1)**4*(2*V**2-2*V+1)*
     3   *2*VC+72*N**2*(V-1)**4*(V**2-V+1)*(2*V**2-2*V+1)*VC-36*(V-1)**4
     4   *(2*V**2-2*V+1)*VC)-54*N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC+108*N
     5   **2*(V-1)**4*(V**2-V+1)*(2*V**2-2*V+1)*VC-54*(V-1)**4*(2*V**2-2
     6   *V+1)*VC)+LOG(S/MP**2)*(LOG(V)*(72*N**2*(V-1)**4*(2*V**2-2*V+1)
     7   *VC-72*N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC)-66*N**4*(V-1)**4*(2*
     8   V**2-2*V+1)**2*VC+24*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)**2*VC+66*
     9   N**2*(V-1)**4*(2*V**2-2*V+1)*VC-24*GTR*N*(V-1)**4*(2*V**2-2*V+1
     :   )*VC)+LOG(V)**2*(18*N**4*(V-1)**4*(2*CQ*(2*V**2-2*V+1)**2-24*V*
     ;   *4+60*V**3-65*V**2+36*V-9)*VC-18*N**2*(V-1)**4*(V**3+2*CQ*(2*V*
     <   *2-2*V+1)-20*V**2+21*V-10)*VC-18*(V-1)**4*(V**2+1)*VC)+LOG(V)*(
     =   3*N**4*(V-1)**4*(52*V**4-74*V**3+26*V**2+14*V-5)*VC+6*N**2*(V-1
     >   )**4*(18*V**4-39*V**3+53*V**2-50*V+16)*VC-9*(V-1)**4*(8*V**2-14
     ?   *V+9)*VC-24*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)**2*VC+24*GTR*N*(V-
     @   1)**4*(2*V**2-2*V+1)*VC)+LOG(1-V)*(9*N**4*(V-1)**4*(12*V**4-34*
     1   V**3+28*V**2-12*V+3)*VC-18*N**2*(V-1)**4*(6*V**4-13*V**3+8*V**2
     2   -7*V+3)*VC+9*(V-1)**4*(4*V**2-10*V+3)*VC)+LOG(1-V)**2*(18*N**4*
     3   (V-1)**4*V*(4*V**2-5*V+2)*VC+18*N**2*(V-1)**5*(V**2-2*V+2)*VC-1
     4   8*(V-1)**4*(V**2-2*V+2)*VC)+LOG(S/MU**2)*(132*N**4*(V-1)**4*(2*
     5   V**2-2*V+1)**2*VC-48*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)**2*VC-132
     6   *N**2*(V-1)**4*(2*V**2-2*V+1)*VC+48*GTR*N*(V-1)**4*(2*V**2-2*V+
     7   1)*VC)+LOG(1-V)*LOG(V)*(-72*N**4*(V-1)**4*(2*V-1)*(2*V**2-2*V+1
     8   )*VC-72*N**2*(V-1)**4*(2*V**2-2*V+1)*VC)+2*N**4*(V-1)**4*(3*(4*
     9    PI**2-51)*CQ*(2*V**2-2*V+1)**2-36* PI**2*V**4-376*V**4+72* PI*
     :   *2*V**3+752*V**3-72* PI**2*V**2-725*V**2+36* PI**2*V+349*V-9* P
     ;   I**2-85)*VC+2*N**2*(V-1)**4*(60* PI**2*V**4+36*V**4-120* PI**2*
     <   V**3-72*V**3-3*(4* PI**2-51)*CQ*(2*V**2-2*V+1)+138* PI**2*V**2+
     =   269*V**2-78* PI**2*V-233*V+24* PI**2+103)*VC+(63*CQ+40)*GTR*N**
     >   3*(V-1)**4*(2*V**2-2*V+1)**2*VC-(63*CQ+40)*GTR*N*(V-1)**4*(2*V*
     ?   *2-2*V+1)*VC-6*(5* PI**2+6)*(V-1)**4*(2*V**2-2*V+1)*VC)/(N**2*(
     @   V-1)**5*V**2)/18.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=D1(V,S)/V*HQQD(V,1)+D1(V,S)/(1.-V)*HQQD(V,2)
     C+2.*(2.*GTR-1.)*A2(V,S)*HFGQD(V)
     C+E0(1.-V,S)/(1.-V)*N/VC*HGQD(V,2)
     1+N/VC*E0(V,S)/V*HGQD(V,1)+D1(V,S)*HFGGD(V)
     2+(D0(1.-V,S)+D0(V,S))*HFGQD(V)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.13) THEN
      AVDEL = DELW*(LOG(S/MP**2)*(LOG(V)*(-72*N**2*(V-1)**4*(V**
     1   2+1)*(V**2-V+1)*VC+36*N**4*(V-1)**4*(V**2+1)**2*VC+36*(V
     2   -1)**6*(V**2+1)*VC)-54*N**2*(V-1)**4*(V**2+1)*(V**2-V+1)
     3   *VC+27*N**4*(V-1)**4*(V**2+1)**2*VC+27*(V-1)**6*(V**2+1)
     4   *VC)+LOG(S/M**2)*(LOG(V)*(72*N**4*(V-1)**4*(V**2+1)**2
     5   *VC-72*N**2*(V-1)**6*(V**2+1)*VC)+LOG(1-V)*(72*N**2*(V-
     6   1)**6*(V**2+1)*VC-72*N**4*(V-1)**4*(V**2+1)**2*VC)-6*N**
     7   2*(V-1)**4*(V**2+1)*(20*V**2-31*V+20)*VC+93*N**4*(V-1)**
     8   4*(V**2+1)**2*VC-24*GTR*N**3*(V-1)**4*(V**2+1)**2*VC+24*
     9   GTR*N*(V-1)**6*(V**2+1)*VC+27*(V-1)**6*(V**2+1)*VC)+LOG
     :   (1-V)*LOG(V)*(36*N**4*(V-1)**4*(2*CQ*(V**2+1)**2-V**4-2
     ;   *V**3-5*V**2-4)*VC-36*N**2*(V-1)**4*(-2*V**4+6*V**3+2*CQ
     <   *(V-1)**2*(V**2+1)-5*V**2+9*V-4)*VC-36*(V-2)*(V-1)**6*V*
     =   VC)+LOG(V)**2*(-18*N**4*(V-1)**4*(2*CQ*(V**2+1)**2-2*V*
     >   *4-2*V**3-11*V**2-9)*VC+18*N**2*(V-1)**5*(-8*V**3+2*CQ*(
     ?   V-1)*(V**2+1)+8*V**2-9*V+10)*VC+18*(V-1)**6*(2*V**2-2*V+
     @   1)*VC)+LOG(V)*(-3*N**4*(V-1)**4*(13*V**4+38*V**2+6*V-5)
     1   *VC+6*N**2*(V-1)**4*(2*V**4-19*V**3+V**2+14*V-16)*VC+9*(
     2   V-1)**6*(3*V**2-4*V+9)*VC+24*GTR*N**3*(V-1)**4*(V**2+1)*
     3   *2*VC-24*GTR*N*(V-1)**6*(V**2+1)*VC)+LOG(1-V)**2*(-18*N
     4   **4*(V-1)**4*(V**2+1)*(2*CQ*(V**2+1)-(V+1)**2)*VC+18*N**
     5   2*(V-1)**6*(2*CQ*(V**2+1)-V)*VC+18*(V-1)**6*(3*V**2-4*V+
     6   3)*VC)+LOG(1-V)*(18*N**4*(V-1)**4*V*(V**2+10*V+1)*VC-18
     7   *N**2*(V-1)**4*V*(V**2+10*V+1)*VC+72*(V-1)**6*V*VC)+LOG
     8   (S/MU**2)*(-132*N**4*(V-1)**4*(V**2+1)**2*VC+48*GTR*N**3
     9   *(V-1)**4*(V**2+1)**2*VC+132*N**2*(V-1)**6*(V**2+1)*VC-4
     :   8*GTR*N*(V-1)**6*(V**2+1)*VC)+2*N**4*(V-1)**4*(6*CQ*(PI*
     ;   *2+3)*(V**2+1)**2+9*PI**2*V**4+85*V**4+18*PI**2*V**3+9*V
     <   **3+9*PI**2*V**2+188*V**2+9*V+85)*VC-2*N**2*(V-1)**4*(-1
     =   2*PI**2*V**4+103*V**4+18*PI**2*V**3-179*V**3+6*CQ*(PI**2
     >   +3)*(V-1)**2*(V**2+1)-15*PI**2*V**2+188*V**2-9*PI**2*V-1
     ?   79*V+6*PI**2+103)*VC+6*(V-1)**6*(5*PI**2*V**2+6*V**2-6*P
     @   I**2*V+2*PI**2+6)*VC-10*(3*CQ+4)*GTR*N**3*(V-1)**4*(V**2
     1   +1)**2*VC+10*(3*CQ+4)*GTR*N*(V-1)**6*(V**2+1)*VC)/(N**2*
     2   (V-1)**6*V**2)/18.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=E0(V,S)/V*HQQD(V,1)+E0(V,S)/(1.-V)*HGGD(V,2)
     C+2.*(2.*GTR-1.)*VC/N*A0(V,S)/(1.-V)*HQGD(V,2)
     C+(B0(V,S)+D0(V,S))/(1.-V)*VC/N*HQGD(V,2)
     1+N/VC*D1(V,S)/V*HGQD(V,1)+E0(V,S)*HFQQD(V)
     2+E0(1.-V,S)*HFQGD(V)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.14) THEN
      AVDEL = DELW*(LOG(S/M**2)*(LOG(1-V)*(24*N**4*(V-1)**4*(V**
     1   2-2*V+2)**2*VC-24*N**2*(V-1)**4*V**2*(V**2-2*V+2)*VC)+L
     2   OG(V)*(24*N**2*(V-1)**4*V**2*(V**2-2*V+2)*VC-24*N**4*(V-
     3   1)**4*(V**2-2*V+2)**2*VC)+2*N**2*(V-1)**4*(V**2-2*V+2)*(
     4   20*V**2-9*V+9)*VC-31*N**4*(V-1)**4*(V**2-2*V+2)**2*VC+8*
     5   GTR*N**3*(V-1)**4*(V**2-2*V+2)**2*VC-8*GTR*N*(V-1)**4*V*
     6   *2*(V**2-2*V+2)*VC-9*(V-1)**4*V**2*(V**2-2*V+2)*VC)+LOG
     7   (S/MP**2)*(LOG(V)*(24*N**2*(V-1)**4*V**2*(V**2-2*V+2)*V
     8   C-24*N**4*(V-1)**4*(V**2-2*V+2)**2*VC)-22*N**4*(V-1)**4*
     9   (V**2-2*V+2)**2*VC+8*GTR*N**3*(V-1)**4*(V**2-2*V+2)**2*V
     :   C+22*N**2*(V-1)**4*V**2*(V**2-2*V+2)*VC-8*GTR*N*(V-1)**4
     ;   *V**2*(V**2-2*V+2)*VC)+LOG(1-V)**2*(6*N**4*(V-1)**4*(2*
     <   CQ*(V**2-2*V+2)**2-2*V**4+10*V**3-21*V**2+20*V-8)*VC-6*N
     =   **2*(V-1)**4*V*(2*CQ*V*(V**2-2*V+2)-V-1)*VC-6*(V-1)**4*V
     >   **2*(2*V**2-2*V+1)*VC)+LOG(V)*(-6*N**2*(V-1)**4*(3*V**4
     ?   -8*V**3+2*V**2+12*V-6)*VC+3*N**4*(V-1)**4*(3*V**4-10*V**
     @   3-2*V**2+24*V-12)*VC+3*(V-1)**4*V**2*(3*V**2+2*V-2)*VC)+
     1   LOG(1-V)*LOG(V)*(12*N**2*(V-1)**4*V*((V-1)*(2*V**2-2*V
     2   +1)+2*CQ*V*(V**2-2*V+2))*VC-12*N**4*(V-1)**4*(2*CQ*(V**2
     3   -2*V+2)**2-4*V**4+18*V**3-37*V**2+36*V-16)*VC+12*(V-1)**
     4   4*V**2*(2*V**2-2*V+1)*VC)+LOG(V)**2*(-6*N**2*(V-1)**4*V
     5   **2*(4*CQ*(V**2-2*V+2)-2*V**2+5*V-5)*VC+6*N**4*(V-1)**4*
     6   (V**2-2*V+2)*(4*CQ*(V**2-2*V+2)-11*V**2+24*V-24)*VC-6*(V
     7   -1)**4*V**2*(3*V**2-2*V+2)*VC)+LOG(1-V)*(-6*N**2*(V-1)*
     8   *4*V*(2*V**2-7*V-1)*VC-6*(V-1)**4*V**2*(2*V+1)*VC+6*N**4
     9   *(V-1)**4*V*(2*V-5)*VC)+LOG(S/MU**2)*(44*N**4*(V-1)**4*
     :   (V**2-2*V+2)**2*VC-16*GTR*N**3*(V-1)**4*(V**2-2*V+2)**2*
     ;   VC-44*N**2*(V-1)**4*V**2*(V**2-2*V+2)*VC+16*GTR*N*(V-1)*
     <   *4*V**2*(V**2-2*V+2)*VC)+2*N**4*(V-1)**4*(CQ*(2*PI**2-57
     =   )*(V**2-2*V+2)**2-3*(2*PI**2*V**4+V**4-10*PI**2*V**3-5*V
     >   **3+21*PI**2*V**2+13*V**2-20*PI**2*V-16*V+8*PI**2+8))*VC
     ?   -2*N**2*(V-1)**4*V*(CQ*(2*PI**2-57)*V*(V**2-2*V+2)-3*(2*
     @   V**3-5*V**2+PI**2*V+5*V+PI**2))*VC-6*(V-1)**4*V**2*(2*PI
     1   **2*V**2+V**2-2*PI**2*V-2*V+PI**2+2)*VC+31*CQ*GTR*N**3*(
     2   V-1)**4*(V**2-2*V+2)**2*VC-31*CQ*GTR*N*(V-1)**4*V**2*(V*
     3   *2-2*V+2)*VC)/(N**2*(V-1)**5*V**3)/6.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=E0(1.-V,S)/V*HQQD(V,1)+E0(1.-V,S)/(1.-V)*HGGD(V,2)
     C+N/VC*F2(V,S)/V*
     1HGQD(V,1)+VC/N*D1(V,S)/(1.-V)*HQGD(V,2)+E0(V,S)*HFGQD(V)
     2+E0(1.-V,S)*HFGGD(V)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.15) THEN
      AVDEL = DELW*(LOG(S/M**2)*(288*N**3*(V-1)**4*(V**2-V+1)**3*
     1   LOG(V)*VC-288*N**3*LOG(1-V)*(V-1)**4*(V**2-V+1)**3*VC+
     2   528*N**3*(V-1)**4*(V**2-V+1)**3*VC-192*GTR*N**2*(V-1)**4
     3   *(V**2-V+1)**3*VC)+LOG(S/MP**2)*(288*N**3*(V-1)**4*(V**
     4   2-V+1)**3*LOG(V)*VC+264*N**3*(V-1)**4*(V**2-V+1)**3*VC-
     5   96*GTR*N**2*(V-1)**4*(V**2-V+1)**3*VC)+LOG(1-V)*LOG(V)
     6   *(72*N**3*(V-1)**4*(4*CQ*(V**2-V+1)**3-4*V**6+16*V**5-37
     7   *V**4+50*V**3-47*V**2+26*V-8)*VC+72*GTR*N**2*(V-1)**5*V*
     8   (2*V**2-2*V+1)*VC)+LOG(V)**2*(-72*N**3*(V-1)**4*(4*CQ*(
     9   V**2-V+1)**3-8*V**6+27*V**5-57*V**4+72*V**3-66*V**2+36*V
     :   -12)*VC-36*GTR*N**2*(V-1)**5*V**2*(V**2+V-1)*VC)+LOG(1-
     ;   V)**2*(36*GTR*N**2*(V-1)**6*V*(V**2-3*V+1)*VC-72*N**3*(V
     <   -1)**4*(2*CQ*(V**2-V+1)**3-2*V**6+7*V**5-14*V**4+16*V**3
     =   -14*V**2+7*V-2)*VC)+LOG(V)*(24*GTR*N**2*(V-1)**4*(V**2-
     >   V+1)*(4*V**4-3*V**3-V**2+8*V-4)*VC-24*N**3*(V-1)**4*(V**
     ?   2-V+1)*(11*V**4-15*V**3+4*V**2+22*V-11)*VC)+LOG(1-V)*(2
     @   4*N**3*(V-1)**4*V*(V**2-V+1)*(7*V**2+8*V+7)*VC-24*GTR*N*
     1   *2*(V-1)**4*V*(V**2-V+1)*(5*V**2-2*V+5)*VC)+LOG(S/MU**2
     2   )*(192*GTR*N**2*(V-1)**4*(V**2-V+1)**3*VC-528*N**3*(V-1)
     3   **4*(V**2-V+1)**3*VC)+4*N**3*(V-1)**4*(378*CQ*(V**2-V+1)
     4   **3+24*PI**2*V**6-134*V**6-72*PI**2*V**5+402*V**5+153*PI
     5   **2*V**4-831*V**4-186*PI**2*V**3+992*V**3+171*PI**2*V**2
     6   -831*V**2-90*PI**2*V+402*V+24*PI**2-134)*VC-4*GTR*N**2*(
     7   V-1)**4*(123*CQ*(V**2-V+1)**3-40*V**6+120*V**5+18*PI**2*
     8   V**4-294*V**4-36*PI**2*V**3+388*V**3+27*PI**2*V**2-294*V
     9   **2-9*PI**2*V+120*V-40)*VC)/((V-1)**6*V**3)/9.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=F2(V,S)/V*HGGD(V,1)+F2(V,S)/(1.-V)*HGGD(V,2)+VC*(E0(1.-V,S)
     1/V*HQGD(V,1)+E0(V,S)/(1.-V)*HQGD(V,2))*4.*GTR/N+F2(V,S)*HFGGD(V)
     2+4.*GTR*D1(V,S)*HFGQD(V)
      AVDEL=AVDEL+ZDEL
 
      ELSE IF (J0.EQ.16) THEN
      AVDEL = DELW*(LOG(S/M**2)*(LOG(1-V)*(24*N**4*(V-1)**4*(2*V
     1   **2-2*V+1)**2*VC-24*N**2*(V-1)**4*(2*V**2-2*V+1)*VC)+LO
     2   G(V)*(24*N**2*(V-1)**4*(2*V**2-2*V+1)*VC-24*N**4*(V-1)**
     3   4*(2*V**2-2*V+1)**2*VC)-44*N**4*(V-1)**4*(2*V**2-2*V+1)*
     4   *2*VC+16*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)**2*VC+44*N**2*
     5   (V-1)**4*(2*V**2-2*V+1)*VC-16*GTR*N*(V-1)**4*(2*V**2-2*V
     6   +1)*VC)+LOG(S/MP**2)*(LOG(V)*(-12*N**4*(V-1)**4*(2*V**
     7   2-2*V+1)**2*VC+24*N**2*(V-1)**4*(V**2-V+1)*(2*V**2-2*V+1
     8   )*VC-12*(V-1)**4*(2*V**2-2*V+1)*VC)-9*N**4*(V-1)**4*(2*V
     9   **2-2*V+1)**2*VC+18*N**2*(V-1)**4*(V**2-V+1)*(2*V**2-2*V
     :   +1)*VC-9*(V-1)**4*(2*V**2-2*V+1)*VC)+LOG(1-V)**2*(6*N**
     ;   4*(V-1)**4*(2*CQ*(2*V**2-2*V+1)**2-8*V**4+20*V**3-21*V**
     <   2+10*V-2)*VC+6*N**2*(V-1)**4*(V**2*(V+1)-2*CQ*(2*V**2-2*
     =   V+1))*VC-6*(V-1)**4*(V**2-2*V+2)*VC)+LOG(V)**2*(6*N**4*
     >   (V-1)**4*(2*CQ*(2*V**2-2*V+1)**2-24*V**4+60*V**3-65*V**2
     ?   +36*V-9)*VC-6*N**2*(V-1)**4*(V**3+2*CQ*(2*V**2-2*V+1)-20
     @   *V**2+21*V-10)*VC-6*(V-1)**4*(V**2+1)*VC)+LOG(1-V)*LOG
     1   (V)*(24*(CQ-2)*N**2*(V-1)**4*(2*V**2-2*V+1)*VC-24*N**4*(
     2   V-1)**4*(2*V**2-2*V+1)*(CQ*(2*V**2-2*V+1)-2*(V-1)**2)*VC
     3   )+LOG(S/MU**2)*(44*N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC-1
     4   6*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)**2*VC-44*N**2*(V-1)**
     5   4*(2*V**2-2*V+1)*VC+16*GTR*N*(V-1)**4*(2*V**2-2*V+1)*VC)
     6   +LOG(1-V)*(6*N**2*(V-1)**4*V*(V**2+7*V-2)*VC-6*N**4*(V-
     7   1)**4*V**2*(5*V-2)*VC-6*(V-1)**4*V*(V+2)*VC)+LOG(V)*(-6
     8   *N**2*(V-1)**5*(V**2-9*V+6)*VC+6*N**4*(V-1)**6*(5*V-3)*V
     9   C-6*(V-3)*(V-1)**5*VC)-2*N**4*(V-1)**4*(3*(8*V**4-16*V**
     :   3+13*V**2-5*V+1)+4*CQ*(PI**2+3)*(2*V**2-2*V+1)**2)*VC+2*
     ;   N**2*(V-1)**4*(3*(5*V**2-5*V+2)+4*CQ*(PI**2+3)*(2*V**2-2
     <   *V+1))*VC+20*CQ*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)**2*VC-2
     =   0*CQ*GTR*N*(V-1)**4*(2*V**2-2*V+1)*VC-6*(V-1)**4*(2*V**2
     >   -2*V+1)*VC)/(N**2*(V-1)**5*V**2)/6.D0
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ZDEL=D1(V,S)/V*HGGD(V,1)+D1(V,S)/(1.-V)*HGGD(V,2)+VC*(E0(V,S)
     1/V*HQGD(V,1)+E0(1.-V,S)/(1.-V)*HQGD(V,2))/N+D1(V,S)*HFQQD(V)
     2+F2(V,S)*HFQGD(V)
      AVDEL=AVDEL+ZDEL
 
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION AVLO(W,V,S)
C     TERME EN (LOG(1-W)/(1-W)+
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      LOPLUS=1.
C    TERME EN LOPLUS
 
      IF (J0.EQ.1) THEN
      AVLO = 4*LOPLUS*(V**2+1)*V1/((V-1)**2*V)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=A0(V,S)/V/W*HQQL(W,V,1)+A0(V*W,S)/(1.-V)*HQQL(W,V,2)
     2+E0(W*V,S)*HGQL(W,V,2)/(1.-V)*N/VC
     2+A0(VZ,S)*HFQQL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.2) THEN
      AVLO =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=(E0(V,S)/V/W*HGQL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQL(W,V,2))
     2*N/VC+HFGQL(W,V)/(1.-V+V*W)*(A0(VZ,S)+A0(1.-VZ,S))
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.3) THEN
      AVLO = 4*LOPLUS*(V**2+1)*V1/((V-1)**2*V)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=A0(V,S)/V/W*HQQL(W,V,1)+A0(V*W,S)/(1.-V)*HQQL(W,V,2)
     2+E0(W*V,S)*HGQL(W,V,2)/(1.-V)*N/VC
     2+A0(VZ,S)*HFQQL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.4) THEN
      AVLO =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=(E0(V,S)/V/W*HGQL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQL(W,V,2))
     2*N/VC+HFGQL(W,V)/(1.-V+V*W)*(A0(VZ,S)+A0(1.-VZ,S))
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.5) THEN
      AVLO = 4*LOPLUS*(2*V**2-2*V+1)*V1/V
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=A2(1.-V,S)/V/W*HQQL(W,V,1)+A2(1.-V*W,S)/(1.-V)*HQQL(W,V,2)
     2+A2(1.-VZ,S)*HFQQL(W,V)/(1.-V+V*W)
     2+D1(VZ,S)*HFQGL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.6) THEN
      AVLO =8*LOPLUS*(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1)/((V-1)
     1   **2*V**3)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=B0(V,S)/V/W*HQQL(W,V,1)+B0(V*W,S)/(1.-V)*HQQL(W,V,2)
     C+E0(1.-V,S)/W/V*N/VC*HGQL(W,V,1)
     1+E0(W*V,S)/(1.-V)*N/VC*HGQL(W,V,2)
     2+B0(VZ,S)*HFQQL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.7) THEN
      AVLO =1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=(E0(V,S)/V/W*HGQL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQL(W,V,2))
     2*N/VC+HFGQL(W,V)/(1.-V+V*W)*B0(VZ,S)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.8) THEN
      AVLO = 1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=N/VC*D1(V,S)/V/W*HGQL(W,V,1)+A0(1.-V*W,S)/(1.-V)*HQGL(W,V,2)
     2/N*VC+HFQGL(W,V)/(1.-V+V*W)*E0(1.-VZ,S)+VC/N*HQGL(W,V,2)/(1.-V)*
     2A2(1.-V*W,S)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.9) THEN
      AVLO = 1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
      ZVLO=N/VC*D1(VY,S)/V/W*HGQL(W,V,1)+A0(1.-V*W,S)/(1.-V)*HQGL(W,V,2)
     2/N*VC+HFQGL(W,V)/(1.-V+V*W)*E0(1.-VZ,S)+VC/N*HQGL(W,V,2)/(1.-V)*
     2A2(V*W,S)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.10) THEN
      AVLO = 1.D-30
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
      ZVLO=N/VC*D1(VY,S)/V/W*HGQL(W,V,1)+D0(1.-V*W,S)/(1.-V)*HQGL(W,V,2)
     2/N*VC+HFQGL(W,V)/(1.-V+V*W)*E0(1.-VZ,S)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.11) THEN
      AVLO =-8*LOPLUS*(2*(V-1)*V**2*V4+(-V**4+3*V**3-4*V**2+2*V-1)*V1)/(
     1   (V-1)**2*V)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=D0(V,S)/V/W*HQQL(W,V,1)+D0(V*W,S)/(1.-V)*HQQL(W,V,2)
     C+E0(W*V,S)/(1.-V)*N/VC*HGQL(W,V,2)
     1+D0(VZ,S)*HFQQL(W,V)
     2/(1.-V+V*W)+D1(VZ,S)*HFQGL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.12) THEN
      AVLO =LOPLUS*(4*N**2*(V-1)*(2*V**2-2*V+1)*(2*V**2-2*V-CQ+3)*VC+4*(
     1   CQ-2)*N**4*(V-1)*(2*V**2-2*V+1)**2*VC-4*(V-1)*(2*V**2-2*V+1)*VC
     2   )/(N**2*(V-1)**2*V**2)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=D1(V,S)/V/W*HQQL(W,V,1)+D1(V*W,S)/(1.-V)*HQQL(W,V,2)
     C+2.*(2.*GTR-1.)*A2(VZ,S)/(1.-V+V*W)*HFGQL(W,V)
     C+E0(1.-V*W,S)/(1.-V)*N/VC*HGQL(W,V,2)
     1+N/VC*E0(V,S)/W/V*HGQL(W,V,1)+D1(VZ,S)*HFGGL(W,V)
     2/(1.-V+V*W)+(D0(1.-VZ,S)+D0(VZ,S))*HFGQL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.13) THEN
      AVLO = LOPLUS*(-4*(CQ-2)*N**4*(V-1)*(V**2+1)**2*VC+4*N**2*(V
     1   -1)*(-3*V**2+4*V+CQ*(V-1)**2-3)*(V**2+1)*VC+4*(V-1)**3*(
     2   V**2+1)*VC)/(N**2*(V-1)**3*V**2)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=E0(V,S)/V/W*HQQL(W,V,1)+E0(V*W,S)/(1.-V)*HGGL(W,V,2)
     C+2.*(2.*GTR-1.)*VC/N*A0(V*W,S)/(1.-V)*HQGL(W,V,2)
     C+(B0(V*W,S)+D0(V*W,S))/(1.-V)*VC/N*HQGL(W,V,2)
     1+N/VC*D1(V,S)/W/V*HGQL(W,V,1)+E0(VZ,S)*HFQQL(W,V)
     2/(1.-V+V*W)+E0(1.-VZ,S)*HFQGL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.14) THEN
      AVLO = LOPLUS*(8*(CQ-2)*N**2*(V-1)**2*(V**2-2*V+2)**2*VC-8*(
     1   CQ-2)*(V-1)**2*V**2*(V**2-2*V+2)*VC)/((V-1)**3*V**3)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=E0(1.-V,S)/V/W*HQQL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGGL(W,V,2)
     C+N/VC*F2(V,S)/W/V*
     1HGQL(W,V,1)+VC/N*D1(V*W,S)/(1.-V)*HQGL(W,V,2)+E0(VZ,S)*HFGQL(W,V)
     2/(1.-V+V*W)+E0(1.-VZ,S)*HFGGL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.15) THEN
      AVLO = -32*(3*CQ-5)*LOPLUS*N**3*(V**2-V+1)**3*VC/((V-1)**2*V
     1   **3)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=F2(V,S)/V/W*HGGL(W,V,1)+F2(V*W,S)/(1.-V)*HGGL(W,V,2)
     C+4.*GTR*VC/N*(E0(VY,S)/W
     1/V*HQGL(W,V,1)+E0(V*W,S)/(1.-V)*HQGL(W,V,2))+F2(VZ,S)*HFGGL(W,V)
     2/(1.-V+V*W)+4.*GTR*D1(VZ,S)*HFGQL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ELSE IF (J0.EQ.16) THEN
      AVLO = LOPLUS*(8*(CQ-2)*N**2*(V-1)**2*(2*V**2-2*V+1)**2*VC-8
     1   *(CQ-2)*(V-1)**2*(2*V**2-2*V+1)*VC)/((V-1)**3*V**2)
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      AJOUT POUR LE SCHEMA DE FACTORIZATION  FINITE NEXT TO LEADING
C  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       VZ=V*W/(1.-V+V*W)
       VY=1.-V
       ZVLO=D1(V,S)/V/W*HGGL(W,V,1)+D1(V*W,S)/(1.-V)*HGGL(W,V,2)
     C+VC/N*(E0(V,S)/W/V
     1*HQGL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HQGL(W,V,2))+D1(VZ,S)*HFQQL(W,V)
     2/(1.-V+V*W)+F2(VZ,S)*HFQGL(W,V)/(1.-V+V*W)
      AVLO=AVLO+ZVLO
 
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION AVGO(W,V)
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      LOVW=LOG((1.-V*W)/(1.-V))/(1.-W)
      LOTVW=LOG(1.-V+V*W)/(1.-W)
      LOW=LOG(W)/(1.-W)
      LOPLUS=1.
C TERME REGULIER AVEC LES LOG (....) AVGO
 
      IF (J0.EQ.1) THEN
      AVGO =-2*((5*LOW*V**2+LOVW*V**2-2*LOTVW*V**2+3*LOW+7*LOVW+2*LOTVW)
     1   *V2+(-6*LOW*V**2+5*LOVW*V**2+2*LOTVW*V**2-4*LOW+3*LOVW+2*LOTVW)
     2   *V1)/((V-1)**2*V)
 
      ELSE IF (J0.EQ.2) THEN
      AVGO =1.D-30
 
      ELSE IF (J0.EQ.3) THEN
      AVGO =2*((10*LOW*V**2-7*LOVW*V**2-LOTVW*V**2+6*LOW-LOVW+LOTVW)*V2+
     1   (LOW*V**2-3*LOVW*V**2-LOTVW*V**2+LOW-5*LOVW-3*LOTVW)*V1)/((V-1)
     2   **2*V)
 
      ELSE IF (J0.EQ.4) THEN
      AVGO =1.D-30
 
      ELSE IF (J0.EQ.5) THEN
      AVGO =2*((16*LOW*V**2-8*LOTVW*V**2-12*LOW*V-2*LOVW*V+2*LOTVW*V+6*L
     1   OW+LOVW-LOTVW)*V2+(2*LOW*V**2-4*LOVW*V**2-8*LOTVW*V**2-2*LOW*V+
     2   6*LOVW*V+10*LOTVW*V+LOW-3*LOVW-5*LOTVW)*V1)/V
 
      ELSE IF (J0.EQ.6) THEN
      AVGO =2*((V-1)*V*(2*LOW*V**2-2*LOVW*V**2-2*LOW*V+2*LOVW*V+11*LOW-7
     1   *LOVW-8*LOTVW)*V4-(V-1)*V*(2*LOW*V**2-2*LOVW*V**2-2*LOW*V+2*LOV
     2   W*V+3*LOW+5*LOVW)*V3-2*V*(4*LOW*V**3-2*LOTVW*V**3-6*LOW*V**2+2*
     3   LOVW*V**2+4*LOTVW*V**2+9*LOW*V+LOVW*V-4*LOTVW*V-3*LOW+LOVW+2*LO
     4   TVW)*V2+2*(5*LOW*V**4-4*LOVW*V**4-2*LOTVW*V**4-8*LOW*V**3+6*LOV
     5   W*V**3+4*LOTVW*V**3+15*LOW*V**2-11*LOVW*V**2-8*LOTVW*V**2-10*LO
     6   W*V+7*LOVW*V+6*LOTVW*V+3*LOW-2*LOVW-2*LOTVW)*V1)/((V-1)**2*V**3
     7   )
 
      ELSE IF (J0.EQ.7) THEN
      AVGO =1.D-30
 
      ELSE IF (J0.EQ.8) THEN
      AVGO = 1.D-30
 
      ELSE IF (J0.EQ.9) THEN
      AVGO = 1.D-30
 
      ELSE IF (J0.EQ.10) THEN
      AVGO = 1.D-30
 
      ELSE IF (J0.EQ.11) THEN
      AVGO =-2*((V-1)*(4*LOW*V**2-7*LOVW*V**2-7*LOTVW*V**2+2*LOVW*V+2*LO
     1   TVW*V-2*LOVW-2*LOTVW)*V4+(V-1)*(16*LOW*V**2-5*LOVW*V**2-5*LOTVW
     2   *V**2-2*LOVW*V-2*LOTVW*V+2*LOVW+2*LOTVW)*V3-2*(8*LOW*V**4-4*LOT
     3   VW*V**4-22*LOW*V**3-LOVW*V**3+9*LOTVW*V**3+28*LOW*V**2-LOVW*V**
     4   2-7*LOTVW*V**2-12*LOW*V-2*LOVW*V+2*LOTVW*V+6*LOW)*V2-2*(LOW*V**
     5   4-2*LOVW*V**4-4*LOTVW*V**4-3*LOW*V**3+7*LOVW*V**3+13*LOTVW*V**3
     6   +4*LOW*V**2-11*LOVW*V**2-17*LOTVW*V**2-2*LOW*V+6*LOVW*V+10*LOTV
     7   W*V+LOW-4*LOVW-4*LOTVW)*V1)/((V-1)**2*V)
 
      ELSE IF (J0.EQ.12) THEN
      AVGO =(-2*N**4*(V-1)*(-4*CQ*LOTVW*(2*V**2-2*V+1)**2+12*LOW*V**4-8*
     1   LOVW*V**4-28*LOW*V**3+20*LOVW*V**3+29*LOW*V**2-21*LOVW*V**2+2*L
     2   OTVW*V**2-16*LOW*V+10*LOVW*V-2*LOTVW*V+4*LOW-2*LOVW+LOTVW)*VC+2
     3   *N**2*(V-1)*(4*LOW*V**4-8*LOVW*V**4-9*LOW*V**3+15*LOVW*V**3-4*C
     4   Q*LOTVW*(2*V**2-2*V+1)+18*LOW*V**2-17*LOVW*V**2+3*LOTVW*V**2-15
     5   *LOW*V+8*LOVW*V-3*LOTVW*V+6*LOW-2*LOVW+2*LOTVW)*VC-2*(V-1)*(3*L
     6   OW*V**2-5*LOVW*V**2-2*LOTVW*V**2-2*LOW*V+6*LOVW*V+2*LOTVW*V+2*L
     7   OW-4*LOVW-3*LOTVW)*VC)/(N**2*(V-1)**2*V**2)
 
      ELSE IF (J0.EQ.13) THEN
      AVGO = (-2*N**2*(V-1)*(4*LOW*V**4-2*LOVW*V**4-4*LOTVW*V**4-6
     1   *LOW*V**3+3*LOVW*V**3+3*LOTVW*V**3+9*LOW*V**2-2*LOVW*V**
     2   2-5*LOTVW*V**2-9*LOW*V+3*LOVW*V+6*LOW-2*LOVW-2*LOTVW)*VC
     3   +2*N**4*(V-1)*(LOW*V**4-3*LOVW*V**4-LOTVW*V**4+2*LOW*V**
     4   3-2*LOVW*V**3+5*LOW*V**2-6*LOVW*V**2-3*LOTVW*V**2-2*LOVW
     5   *V-2*LOTVW*V+4*LOW-3*LOVW-2*LOTVW)*VC+2*(V-1)**3*(3*LOW*
     6   V**2-3*LOVW*V**2-3*LOTVW*V**2-2*LOW*V+4*LOVW*V+2*LOTVW*V
     7   +2*LOW-3*LOVW-4*LOTVW)*VC)/(N**2*(V-1)**3*V**2)
 
      ELSE IF (J0.EQ.14) THEN
      AVGO = (-2*N**4*(V-1)**2*(-4*CQ*LOTVW*(V**2-2*V+2)**2+4*LOW*
     1   V**4-4*LOVW*V**4+LOTVW*V**4-18*LOW*V**3+18*LOVW*V**3-4*L
     2   OTVW*V**3+38*LOW*V**2-37*LOVW*V**2+7*LOTVW*V**2-40*LOW*V
     3   +36*LOVW*V-4*LOTVW*V+20*LOW-16*LOVW)*VC-2*N**2*(V-1)**2*
     4   (2*LOVW*V**4+3*LOW*V**3-4*LOVW*V**3-LOTVW*V**3+4*CQ*LOTV
     5   W*V**2*(V**2-2*V+2)-7*LOW*V**2+5*LOVW*V**2+8*LOW*V+LOVW*
     6   V-LOTVW*V-4*LOW)*VC-2*(V-1)**2*V**2*(4*LOW*V**2-2*LOVW*V
     7   **2-LOTVW*V**2-4*LOW*V+2*LOVW*V+4*LOW-LOVW-LOTVW)*VC)/(N
     8   **2*(V-1)**3*V**3)
 
      ELSE IF (J0.EQ.15) THEN
      AVGO = (16*N**3*(V-1)*(2*CQ*(LOW-2*LOTVW)*(V**2-V+1)**3+2*LO
     1   W*V**6-4*LOVW*V**6-7*LOW*V**5+13*LOVW*V**5+15*LOW*V**4-2
     2   6*LOVW*V**4+2*LOTVW*V**4-20*LOW*V**3+30*LOVW*V**3-4*LOTV
     3   W*V**3+20*LOW*V**2-26*LOVW*V**2+3*LOTVW*V**2-12*LOW*V+13
     4   *LOVW*V-LOTVW*V+4*LOW-4*LOVW)*VC-8*GTR*N**2*(V-1)**2*V*(
     5   LOW*V**3+LOVW*V**3+LOW*V**2-4*LOVW*V**2-LOTVW*V**2-LOW*V
     6   +4*LOVW*V+LOTVW*V-LOVW+LOTVW)*VC)/((V-1)**3*V**3)
 
      ELSE IF (J0.EQ.16) THEN
      AVGO = (-2*N**4*(V-1)**2*(2*CQ*LOW*(2*V**2-2*V+1)**2+8*LOW*V
     1   **4-16*LOVW*V**4-16*LOTVW*V**4-20*LOW*V**3+36*LOVW*V**3+
     2   32*LOTVW*V**3+21*LOW*V**2-37*LOVW*V**2-30*LOTVW*V**2-12*
     3   LOW*V+18*LOVW*V+14*LOTVW*V+3*LOW-4*LOVW-3*LOTVW)*VC-2*N*
     4   *2*(V-1)**2*(LOW*V**3+LOVW*V**3-2*CQ*LOW*(2*V**2-2*V+1)-
     5   8*LOW*V**2+5*LOVW*V**2+5*LOTVW*V**2+9*LOW*V-4*LOVW*V-5*L
     6   OTVW*V-4*LOW+2*LOVW+2*LOTVW)*VC-2*(V-1)**2*(LOW*V**2-LOV
     7   W*V**2-2*LOTVW*V**2+2*LOVW*V+2*LOTVW*V+LOW-2*LOVW-3*LOTV
     8   W)*VC)/(N**2*(V-1)**3*V**2)
 
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/EDF/J0
      IF (J0.EQ.1) THEN
      STRUV=STRUV1(W,V,X3,S)
      ELSE IF (J0.EQ.2) THEN
      STRUV=STRUV2(W,V,X3,S)
      ELSE IF (J0.EQ.3) THEN
      STRUV=STRUV3(W,V,X3,S)
      ELSE IF (J0.EQ.4) THEN
      STRUV=STRUV4(W,V,X3,S)
      ELSE IF (J0.EQ.5) THEN
      STRUV=STRUV5(W,V,X3,S)
      ELSE IF (J0.EQ.6) THEN
      STRUV=STRUV6(W,V,X3,S)
      ELSE IF (J0.EQ.7) THEN
      STRUV=STRUV7(W,V,X3,S)
      ELSE IF (J0.EQ.8) THEN
      STRUV=STRUV8(W,V,X3,S)
      ELSE IF (J0.EQ.9) THEN
      STRUV=STRUV9(W,V,X3,S)
      ELSE IF (J0.EQ.10) THEN
      STRUV=STRUV10(W,V,X3,S)
      ELSE IF (J0.EQ.11) THEN
      STRUV=STRUV11(W,V,X3,S)
      ELSE IF (J0.EQ.12) THEN
      STRUV=STRUV12(W,V,X3,S)
      ELSE IF (J0.EQ.13) THEN
      STRUV=STRUV13(W,V,X3,S)
      ELSE IF (J0.EQ.14) THEN
      STRUV=STRUV14(W,V,X3,S)
      ELSE IF (J0.EQ.15) THEN
      STRUV=STRUV15(W,V,X3,S)
      ELSE IF (J0.EQ.16) THEN
      STRUV=STRUV16(W,V,X3,S)
      ENDIF
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV1(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(V1*(-2*V**6*W**6+2*V**5*(V+1)*W**5-V**4*(V**2+6*V-
     1   1)*W**4+2*V**4*(V+7)*W**3-2*V**2*(3*V**2+V+7)*W**2+2*V*(V**2+5)
     2   *W-V**2-1)+CQ*V1*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+2*
     3   V**2*W**2+1)+2*(V-1)*V*V2*W*(3*V**3*W**3+V**2*W**2-V*(2*V+5)*W+
     4   3)+2*CQ*V*V2*W*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/
     5   ((V-1)**2*V*W*(V*W-1)**3)
 
      LV1 = -LOG(1-V)*(V1*(4*V**6*W**6-4*V**5*(V+1)*W**5+V**4*(2*V**2+9*
     1   V+1)*W**4-V**3*(7*V**2+11*V+10)*W**3+V**2*(8*V**2+13*V+15)*W**2
     2   -V*(7*V**2+3*V+10)*W+2*(V**2+1))+CQ*V1*(2*V**2*W**2-2*V*(V+1)*W
     3   +V**2+1)*(V**4*W**4+2*V**2*W**2+1)+2*CQ*V*V2*W*(V**2*W**2+1)*(2
     4   *V**2*W**2-2*V*(V+1)*W+V**2+1)+2*V*V2*W*(V**2*W**2-2*V*W+1)*(2*
     5   V**2*W**2+V**2+4*V-3))/((V-1)**2*V*W*(V*W-1)**3)
 
 
      LV = LOG(V)*(V1*(-2*V**6*W**6+2*V**5*(V+1)*W**5-V**3*(V**3+13*V**2
     1   -7*V+1)*W**4+V**2*(5*V**3+27*V**2-3*V+3)*W**3-V*(13*V**3+7*V**2
     2   +27*V+3)*W**2+(7*V**3-V**2+21*V+1)*W-2*(V**2+1))-2*V*V2*W*(2*V*
     3   *4*W**4-4*V**3*(2*V-1)*W**3+V**2*(3*V**2+10*V-1)*W**2-4*V*(V**2
     4   +V+3)*W+3*V**2-2*V+9)+CQ*V1*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V
     5   **4*W**4+2*V**2*W**2+1)+2*CQ*V*V2*W*(V**2*W**2+1)*(2*V**2*W**2-
     6   2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V*W-1)**3)
 
 
      LVW = 2*(V1*(2*V**3*W**3-2*(V-2)*V**2*W**2+V*(V**2+V+4)*W-V**2-1)+
     1   V*V2*W*(4*V*W-V+5))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)
     2   *LOG(1-V*W)/((V-1)**2*V*W*(V*W-1)**4)
 
 
      LTVW = 2*(V1*(2*V**2*W**2-2*(V-2)*V*W+V**2-V+4)+(1-V)*V2)*(V**4*W*
     1   *4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*LOG(V*W-V+1)/((V-1)**2*(V*W
     2   -1)**4)
 
 
      LW = -(V1*(V*(11*V**2-6*V+1)*W**2-(5*V**3+3*V**2+11*V+1)*W+3*(V**2
     1   +1))+2*V*V2*W*(2*V**2*W**2-5*(V-1)*V*W+3*V**2+V+6))*(V**3*W**3-
     2   3*V**2*W**2+3*V*W-1)*LOG(W)/((V-1)**2*V*W*(V*W-1)**4)
 
 
      CVC = (3*(V-1)*V1*(V*W-1)*(6*V**7*W**7-4*V**6*(3*V+1)*W**6+V**4*(7
     1   *V**3+6*V**2-5*V-2)*W**5-V**3*(V**4-V**3+13*V**2-37*V-4)*W**4+2
     2   *V**3*(7*V**2-22*V-2)*W**3-2*V*(7*V**4-25*V**3+17*V**2-3*V+2)*W
     3   **2-(8*V**4-27*V**3+26*V**2-11*V-2)*W-(V-1)*(V**2+4*V+1))+12*(V
     4   -1)*V2*W*(V*W-1)*(2*V**6*W**5-V**4*(6*V**2-3*V+1)*W**4+V**3*(6*
     5   V**3-3*V**2-V+2)*W**3-V**3*(4*V**3-16*V**2+27*V-11)*W**2-V*(4*V
     6   **4-10*V**3+3*V**2-V+2)*W+(V-1)*(2*V**3-6*V**2+3*V-1))-6*AL*(V-
     7   1)*V1*(V*W-1)*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*
     8   W**4+2*V**2*W**2+1)-8*CQ*(V-1)**2*V1*(V*W-V+1)*(V**2*W**2-2*V*W
     9   +1)*(V**4*W**4+2*V**2*W**2+1)-12*AL*(V-1)*V*V2*W*(V*W-1)*(V*W-V
     :   +1)*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)-16*CQ*(V-1)*
     ;   *2*V*V2*W*(V*W-V+1)*(V**2*W**2+1)*(V**2*W**2-2*V*W+1))/((V-1)**
     <   3*V*W*(V*W-1)**4*(V*W-V+1))/6.D0
 
 
      LM = LOG(S/M**2)*(-V1*(2*V**6*W**6-2*V**5*(V+1)*W**5+V**3*(V**3+2*
     1   V**2+4*V+1)*W**4-V**2*(2*V**3+10*V**2+5*V+3)*W**3+V*(5*V**3+8*V
     2   **2+10*V+3)*W**2-(4*V**3+2*V**2+9*V+1)*W+2*(V**2+1))-2*V*V2*W*(
     3   V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V
     4   *W-1)**3)
 
      LMP = -LOG(S/MP**2)*V1*(2*V**3*W**3-2*V**2*(2*V-3)*W**2+V*(3*V**2-
     1   8*V+9)*W-(V-1)*(V**2-V+4))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4
     2   *V*W+1)/((V-1)**2*(V*W-1)**4*(V*W-V+1))
 
      STRUV1=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV2(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+V**4*(16*V
     1   **4+24*V**3+11*V**2+1)*W**6-2*V**4*(8*V**4+12*V**3+21*V**2-V+2)
     2   *W**5+2*V**2*(8*V**6+4*V**5+29*V**4+9*V**3-6*V**2+V-1)*W**4-2*V
     3   **2*(6*V**6+2*V**5+11*V**4+30*V**3-16*V**2-2*V-3)*W**3+(4*V**8+
     4   12*V**7-17*V**6+54*V**5-5*V**4-32*V**3+V**2-2*V+1)*W**2-2*(V-1)
     5   *(4*V**6-2*V**5+7*V**4+6*V**3+6*V**2-4*V-1)*W+4*(V-1)**2*(V**2+
     6   1)*(V**2-V+2))-V1*(-4*V**8*W**8+4*V**7*(3*V+1)*W**7-V**4*(16*V*
     7   *4+2*V**3+17*V**2+1)*W**6+2*V**4*(8*V**4-11*V**3+25*V**2+4*V+2)
     8   *W**5-2*V**2*(8*V**6-19*V**5+20*V**4+33*V**3-10*V**2+V-1)*W**4+
     9   2*V**2*(6*V**6-9*V**5-11*V**4+63*V**3-16*V**2-6*V-3)*W**3-(4*V*
     :   *8+8*V**7-49*V**6+88*V**5+13*V**4-46*V**3-V**2-2*V+1)*W**2+2*(V
     ;   -1)*(4*V**6-6*V**5+4*V**4+11*V**3+9*V**2-5*V-1)*W-4*(V-1)**2*(V
     <   **2+1)*(V**2-2*V+3))+2*V2*(V*W-1)*(V*W-V+1)*(V**3*(6*V**2-3*V+1
     =   )*W**4-V**3*(8*V**2+9*V-5)*W**3+V*(6*V**4+9*V**3+15*V**2-13*V-1
     >   )*W**2-V*(9*V**3+V**2+11*V-13)*W+4*(V-1)*(V**2+1))+2*CQ*V2*(V*W
     ?   -1)*(V*W-V-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**3*W**3-V*(
     @   V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1)))/((V-1)**2*V**2*W**
     1   2*(V*W-1)**2*(V*W-V+1)**2)
 
 
      LV1 = LOG(1-V)*(2*V1*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**4*W**4-2*
     1   V**3*(V+2)*W**3+V**2*(V**2+2*V+9)*W**2-2*V*(V**2+2*V+3)*W+2*(V+
     2   1)*(V**2+1))+CQ*V1*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(
     3   V**4*W**4-4*V**3*W**3+8*V**2*W**2-8*V*W+4)-4*V2*(V*W-V-1)*(V**2
     4   *W**2-2*V*W+1)*(V**3*W**3-V**2*(V+1)*W**2+V*(3*V**2-2*V+1)*W-(V
     5   -1)*(V**2+1))-2*CQ*V2*(V*W-1)*(V*W-V+1)*(V**2*W**2-2*V*W+2)*(2*
     6   V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V**2*W**2*(V*W-1)**2*(
     7   V*W-V+1))
 
 
      LV = LOG(V)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+3*V**6*(V+1)*
     1   (5*V+3)*W**6-2*V**5*(6*V**3+12*V**2+17*V-1)*W**5+V**4*(9*V**4+6
     2   *V**3+46*V**2+14*V-15)*W**4-2*V**3*(3*V+5)*(V**4-2*V**3+6*V**2-
     3   2*V-1)*W**3+2*V**2*(V**6+2*V**5-7*V**4+20*V**3+3*V**2-18*V+3)*W
     4   **2-2*(V-1)*V*(2*V**5-3*V**4+4*V**3+4*V**2+6*V-5)*W+2*(V-1)**2*
     5   (V**2+1)*(V**2-2*V+3))-V1*(-4*V**8*W**8+4*V**7*(3*V+2)*W**7-V**
     6   4*(16*V**4+12*V**3+23*V**2+1)*W**6+2*V**4*(8*V**4-8*V**3+37*V**
     7   2+3*V+2)*W**5-2*V**2*(8*V**6-22*V**5+35*V**4+37*V**3-14*V**2+V-
     8   1)*W**4+2*V**2*(6*V**6-14*V**5-5*V**4+74*V**3-20*V**2-10*V-3)*W
     9   **3-(4*V**8+4*V**7-57*V**6+114*V**5+7*V**4-52*V**3-3*V**2-2*V+1
     :   )*W**2+2*(V-1)*(4*V**6-10*V**5+7*V**4+10*V**3+14*V**2-8*V-1)*W-
     ;   4*(V-1)**2*(V**2+1)*(V**2-3*V+4))-2*V2*(V*W-1)*(V*W-V+1)*(2*V**
     <   5*W**5-V**3*(11*V**2-2*V+1)*W**4+2*V**3*(4*V**2+11*V-5)*W**3-V*
     =   (V**4+16*V**3+24*V**2-20*V-1)*W**2-2*V*(V**4-4*V**3-2*V**2-8*V+
     >   9)*W+2*(V-3)*(V-1)*(V**2+1))+2*CQ*V2*(V*W-1)*(V**2*W**2-2*V*W+2
     ?   )*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**2*W**2-2*V*(V+1)*W+V**
     @   2+1))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
 
      LVW = -2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*
     1   (2*V**4*W**4-2*V**3*(V+3)*W**3+V**2*(V**2+3*V+12)*W**2-V*(V+3)*
     2   *2*W+4*(V**2+1))+V2*(-4*V**3*W**3+V**2*(5*V+7)*W**2-V*(3*V**2+6
     3   *V+7)*W+4*(V**2+1)))*LOG(1-V*W)/((V-1)**2*V**2*W**2*(V*W-1)**2*
     4   (V*W-V+1)**2)
 
 
      LTVW = -4*V1*(V**2*W**2-2*V*W+1)*(CQ*(V**2*W**2-2*V**2*W+V**2+1)*(
     1   V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)
     2   **4)-(V-1)*V*W*(V**3*W**3-V**2*(5*V-3)*W**2+(V-1)*V*(5*V+1)*W-(
     3   V-1)**2*(V+1)))*LOG(V*W-V+1)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*
     4   W-V+1)**2)
 
 
      LW = -(V**2*W**2-2*V*W+1)*(2*V1*(V*W-V+1)*(2*V**3*W**3-V**2*(V**2+
     1   7)*W**2+2*V*(V**3-V**2+3*V+3)*W-2*(V**2+1)*(V**2-V+2))+4*V2*(V*
     2   W-V-1)*(V**3*W**3-V**2*(3*V-1)*W**2+V*(V**2+2*V-1)*W-(V-1)*(V**
     3   2+1))-2*CQ*V*(V**2+1)*V2*(V*W-V+1)*(W**2-2*W+2)-CQ*(V**2+1)**2*
     4   V1*(V*W-V+1)*(W**2-2*W+2))*LOG(W)/((V-1)**2*V**2*W**2*(V*W-1)**
     5   2*(V*W-V+1))
 
 
      CVC = (3*V1*(8*V**7*W**7+V**4*(2*V**4-24*V**3-13*V**2+5*V-2)*W**6-
     1   V**4*(4*V**4-32*V**3-46*V**2+17*V-7)*W**5+2*V**2*(V**6-12*V**5-
     2   36*V**4+7*V**3+5*V**2-7*V+2)*W**4+2*V**2*(8*V**5+19*V**4+22*V**
     3   3-27*V**2+11*V-5)*W**3-(8*V**7+17*V**6+3*V**5+40*V**4-66*V**3+2
     4   1*V**2-9*V+2)*W**2+(V-1)*(16*V**5-11*V**4+34*V**3+2*V**2-6*V-3)
     5   *W-8*(V-1)**2*V*(V**2-V+2))+3*AL*V1*(V**2*W**2-2*(V-1)*V*W+(V-1
     6   )**2)*(2*V**6*W**6-2*V**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**
     7   2+1)*W**4-2*V*(V**5+3*V**4+10*V**3+20*V**2+V+1)*W**3+(2*V**6+4*
     8   V**5+13*V**4+24*V**3+36*V**2+4*V+1)*W**2-2*(2*V**5+V**4+8*V**3+
     9   6*V**2+10*V+1)*W+2*(V**2+1)*(V**2+3))+12*V*V2*(V*W-1)*(V**5*W**
     :   6-V**3*(4*V**2-V+1)*W**5+V**2*(5*V**3+2*V**2+3*V-1)*W**4-V*(5*V
     ;   **4-2*V**3+16*V**2-8*V-1)*W**3+(4*V**5-3*V**4+12*V**3-4*V**2-6*
     <   V+1)*W**2-(V-1)*(V**4+2*V**3+4*V**2+4*V-3)*W+(V-1)**2*(V**2+3))
     =   +4*CQ*V1*(V*W-1)*(V*W-V+1)*(V**5*(V+2)*W**5-V**4*(V**2+5*V-3)*W
     >   **4+V**2*(V**4+7*V**3-V**2-4*V+1)*W**3-V**2*(V**4+3*V**3+4*V**2
     ?   -5*V+1)*W**2+(V-1)*(4*V**4-6*V**3+13*V**2-4*V+1)*W-(V-1)**2*(3*
     @   V**2-6*V+7))-6*AL*V2*(V*W-1)*(V*W-V-1)*(V**2*W**2-2*(V-1)*V*W+(
     1   V-1)**2)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**
     2   2+1))-8*CQ*V2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)*
     3   *2)*((V-1)*V**2*W**2-V*(V**2+2*V-1)*W+2*(V-1)))/((V-1)**2*V**2*
     4   W**2*(V*W-1)**2*(V*W-V+1)**2)/3.D0
 
      LM = LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**6
     1   -2*V**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**
     2   5+3*V**4+10*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V*
     3   *3+36*V**2+4*V+1)*W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2
     4   *(V**2+1)*(V**2+3))-2*V2*(V*W-1)*(V*W-V-1)*(2*V**3*W**3-V*(V**2
     5   +4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1)))/((V-1)**2*V**2*W**2*(V
     6   *W-1)**2*(V*W-V+1)**2)
 
      LMP = 2*LOG(S/MP**2)*V1*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+V*
     1   *2+1)*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*
     2   W+(V-1)**4)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
      STRUV2=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV3(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(V1*(-2*V**6*W**6+2*V**5*(V+1)*W**5-V**4*(V**2+2*V+
     1   3)*W**4+2*V**3*(V**2+5*V+2)*W**3-2*V**2*(3*V**2+3*V+5)*W**2+2*V
     2   *(V**2+2*V+3)*W-V**2-1)+CQ*V1*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*
     3   (V**4*W**4+2*V**2*W**2+1)-2*(V-1)*V*V2*W*(3*V**3*W**3-7*V**2*W*
     4   *2+V*(2*V-1)*W+3)+2*CQ*V*V2*W*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V
     5   +1)*W+V**2+1))/((V-1)**2*V*W*(V*W-1)**3)
 
 
      LV1 = -LOG(1-V)*(V1*(4*V**6*W**6-4*V**5*(V+1)*W**5+V**4*(2*V**2+9*
     1   V+1)*W**4-V**3*(7*V**2+3*V+18)*W**3+V**2*(8*V**2-3*V+31)*W**2-V
     2   *(7*V**2-5*V+18)*W+2*(V**2+1))+CQ*V1*(2*V**2*W**2-2*V*(V+1)*W+V
     3   **2+1)*(V**4*W**4+2*V**2*W**2+1)+2*CQ*V*V2*W*(V**2*W**2+1)*(2*V
     4   **2*W**2-2*V*(V+1)*W+V**2+1)+2*V*V2*W*(V**2*W**2-2*V*W+1)*(2*V*
     5   *2*W**2+V**2-8*V+9))/((V-1)**2*V*W*(V*W-1)**3)
 
 
      LV = LOG(V)*(V1*(-2*V**6*W**6+2*V**5*(V+1)*W**5-V**3*(V**3+V**2+5*
     1   V+1)*W**4+V**2*(V**3+7*V**2+5*V+3)*W**3-V*(5*V**3+3*V**2+7*V+3)
     2   *W**2+(3*V**3+3*V**2+5*V+1)*W-2*(V**2+1))-2*V*V2*W*(2*V**4*W**4
     3   +2*V**3*(5*V-7)*W**3-V**2*(3*V**2+20*V-11)*W**2+2*V*(4*V**2+V+9
     4   )*W-3*V**2+4*V-15)+CQ*V1*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4
     5   *W**4+2*V**2*W**2+1)+2*CQ*V*V2*W*(V**2*W**2+1)*(2*V**2*W**2-2*V
     6   *(V+1)*W+V**2+1))/((V-1)**2*V*W*(V*W-1)**3)
 
 
      LVW = 2*(V1*(2*V**3*W**3-2*(V-2)*V**2*W**2+V*(V**2-V+6)*W-V**2-1)+
     1   V*V2*W*(4*V*W+5*V-1))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+
     2   1)*LOG(1-V*W)/((V-1)**2*V*W*(V*W-1)**4)
 
 
      LTVW = 2*(V1*(2*V**2*W**2-2*(V-2)*V*W+V**2-2*V+5)+2*(V-1)*V2)*(V**
     1   4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*LOG(V*W-V+1)/((V-1)**2*
     2   (V*W-1)**4)
 
 
      LW = -(V1*(V*(V**2+4*V+1)*W**2-(V+1)*(V**2+1)*W+3*(V**2+1))+2*V*V2
     1   *W*(2*V**2*W**2+10*(V-1)*V*W-3*V**2-2*V-9))*(V**3*W**3-3*V**2*W
     2   **2+3*V*W-1)*LOG(W)/((V-1)**2*V*W*(V*W-1)**4)
 
 
      CVC = (3*(V-1)*V1*(V*W-1)*(6*V**7*W**7-4*V**6*(3*V+1)*W**6+V**4*(7
     1   *V**3+2*V**2+3*V-6)*W**5-V**3*(V**4-9*V**3+21*V**2-29*V-12)*W**
     2   4-2*V**3*(4*V**3-19*V**2+34*V-2)*W**3-2*V*(7*V**4-21*V**3+13*V*
     3   *2-7*V+6)*W**2+(V+2)*(3*V**2-4*V+3)*W-(V-1)*(V**2+4*V+1))+12*(V
     4   -1)*V2*W*(V*W-1)*(2*V**6*W**5-V**4*(3*V**2+3*V-2)*W**4+V**3*(3*
     5   V**2+5*V-4)*W**3+V**3*(2*V**3-2*V**2-9*V+5)*W**2-V*(4*V**4-16*V
     6   **3+9*V**2+5*V-4)*W-2*(V-1)*(2*V-1)*(V**2-V+1))-6*AL*(V-1)*V1*(
     7   V*W-1)*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+2*
     8   V**2*W**2+1)-8*CQ*(V-1)**2*V1*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V*
     9   *4*W**4+2*V**2*W**2+1)-12*AL*(V-1)*V*V2*W*(V*W-1)*(V*W-V+1)*(V*
     :   *2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)-16*CQ*(V-1)**2*V*V2
     ;   *W*(V*W-V+1)*(V**2*W**2+1)*(V**2*W**2-2*V*W+1))/((V-1)**3*V*W*(
     <   V*W-1)**4*(V*W-V+1))/6.D0
 
 
      LM = LOG(S/M**2)*(-V1*(2*V**6*W**6-2*V**5*(V+1)*W**5+V**3*(V**3+2*
     1   V**2+4*V+1)*W**4-V**2*(2*V**3+10*V**2+5*V+3)*W**3+V*(5*V**3+8*V
     2   **2+10*V+3)*W**2-(4*V**3+2*V**2+9*V+1)*W+2*(V**2+1))-2*V*V2*W*(
     3   V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V
     4   *W-1)**3)
 
      LMP = -LOG(S/MP**2)*V1*(2*V**3*W**3-2*V**2*(2*V-3)*W**2+V*(3*V**2-
     1   8*V+9)*W-(V-1)*(V**2-V+4))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4
     2   *V*W+1)/((V-1)**2*(V*W-1)**4*(V*W-V+1))
 
      STRUV3=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV4(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+V**4*(16*V
     1   **4+24*V**3+11*V**2+1)*W**6-2*V**4*(8*V**4+12*V**3+21*V**2-V+2)
     2   *W**5+2*V**2*(8*V**6+4*V**5+29*V**4+9*V**3-6*V**2+V-1)*W**4-2*V
     3   **2*(6*V**6+2*V**5+11*V**4+30*V**3-16*V**2-2*V-3)*W**3+(4*V**8+
     4   12*V**7-17*V**6+54*V**5-5*V**4-32*V**3+V**2-2*V+1)*W**2-2*(V-1)
     5   *(4*V**6-2*V**5+7*V**4+6*V**3+6*V**2-4*V-1)*W+4*(V-1)**2*(V**2+
     6   1)*(V**2-V+2))-V1*(-4*V**8*W**8+4*V**7*(3*V+1)*W**7-V**4*(16*V*
     7   *4+12*V**3+7*V**2+1)*W**6+2*V**4*(8*V**4+2*V**3+19*V**2-3*V+2)*
     8   W**5-2*V**2*(8*V**6-6*V**5+27*V**4+13*V**3-10*V**2+V-1)*W**4+2*
     9   V**2*(6*V**6-4*V**5+5*V**4+46*V**3-24*V**2-2*V-3)*W**3-(4*V**8+
     :   8*V**7-29*V**6+74*V**5+7*V**4-56*V**3+9*V**2-2*V+1)*W**2+2*(V-1
     ;   )*(4*V**6-6*V**5+9*V**4+6*V**3+12*V**2-8*V-1)*W-4*(V-1)**2*(V**
     <   2+1)*(V**2-2*V+3))-2*V2*(V*W-1)*(V*W-V+1)*(V**3*(9*V**2-12*V-1)
     =   *W**4-4*V**3*(4*V**2-3*V-4)*W**3+V*(9*V**4-24*V**2-2*V+1)*W**2-
     >   2*V*(3*V**3-8*V**2-V+2)*W-4*(V-1)*(V**2+1))+2*CQ*V2*(V*W-1)*(V*
     ?   W-V-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**3*W**3-V*(V**2+4*
     @   V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1)))/((V-1)**2*V**2*W**2*(V*W-
     1   1)**2*(V*W-V+1)**2)
 
 
      LV1 = LOG(1-V)*(2*V1*(V**2*W**2-2*V*W+1)*(2*V**5*W**5-2*V**4*(2*V+
     1   1)*W**4+V**3*(3*V**2+4*V+5)*W**3-V**2*(V**3+5*V**2+7*V-1)*W**2+
     2   2*V*(3*V**3+V**2+V-1)*W-2*(V-1)*(V+1)*(V**2+1))+CQ*V1*(V*W-V+1)
     3   *(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4-4*V**3*W**3+8*V**2
     4   *W**2-8*V*W+4)-4*V2*(V*W-V-1)*(V**2*W**2-2*V*W+1)*(V**3*W**3-V*
     5   *2*(V+1)*W**2+2*V*(2*V-1)*W-(V-1)*(V**2+1))-2*CQ*V2*(V*W-1)*(V*
     6   W-V+1)*(V**2*W**2-2*V*W+2)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((
     7   V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1))
 
 
      LV = LOG(V)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+3*V**6*(V+1)*
     1   (5*V+3)*W**6-2*V**5*(6*V**3+12*V**2+17*V-1)*W**5+V**4*(9*V**4+6
     2   *V**3+46*V**2+14*V-15)*W**4-2*V**3*(3*V+5)*(V**4-2*V**3+6*V**2-
     3   2*V-1)*W**3+2*V**2*(V**6+2*V**5-7*V**4+20*V**3+3*V**2-18*V+3)*W
     4   **2-2*(V-1)*V*(2*V**5-3*V**4+4*V**3+4*V**2+6*V-5)*W+2*(V-1)**2*
     5   (V**2+1)*(V**2-2*V+3))-V1*(-4*V**8*W**8+4*V**7*(3*V+2)*W**7-V**
     6   4*(16*V**4+24*V**3+11*V**2+1)*W**6+2*V**4*(8*V**4+6*V**3+33*V**
     7   2-7*V+2)*W**5-2*V**2*(8*V**6-10*V**5+47*V**4+15*V**3-16*V**2+V-
     8   1)*W**4+2*V**2*(6*V**6-10*V**5+11*V**4+64*V**3-36*V**2-4*V-3)*W
     9   **3-(4*V**8+4*V**7-41*V**6+106*V**5+7*V**4-76*V**3+13*V**2-2*V+
     :   1)*W**2+2*(V-1)*(4*V**6-10*V**5+11*V**4+6*V**3+18*V**2-12*V-1)*
     ;   W-4*(V-1)**2*(V**2+1)*(V**2-3*V+4))-2*V2*(V*W-1)*(V*W-V+1)*(2*V
     <   **5*W**5+V**3*(7*V**2-16*V-1)*W**4-4*V**3*(4*V**2-4*V-5)*W**3+V
     =   *(11*V**4-4*V**3-24*V**2-4*V+1)*W**2-2*V*(V**4+2*V**3-8*V**2-2*
     >   V+3)*W+2*(V-3)*(V-1)*(V**2+1))+2*CQ*V2*(V*W-1)*(V**2*W**2-2*V*W
     ?   +2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**2*W**2-2*V*(V+1)*W+V
     @   **2+1))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
 
      LVW = -2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*
     1   (2*V**4*W**4-2*V**3*(V+3)*W**3+V**2*(V**2+4*V+11)*W**2-2*V*(V**
     2   2+3*V+4)*W+4*(V**2+1))-2*V2*(2*V**3*W**3-V**2*(V+5)*W**2+V*(3*V
     3   +5)*W-2*(V**2+1)))*LOG(1-V*W)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V
     4   *W-V+1)**2)
 
 
      LTVW = -4*(V**2*W**2-2*V*W+1)*(CQ*V1*(V**2*W**2-2*V**2*W+V**2+1)*(
     1   V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)
     2   **4)-3*(V-1)*V*V2*W*(V*W-V-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)+
     3   2*(V-1)**2*V**3*V1*(W-1)*W**2)*LOG(V*W-V+1)/((V-1)**2*V**2*W**2
     4   *(V*W-1)**2*(V*W-V+1)**2)
 
 
      LW = -(V**2*W**2-2*V*W+1)*(2*V1*(2*V**4*W**4-V**3*(V+1)*(V+3)*W**3
     1   +V**2*(3*V**3-V**2+13*V-3)*W**2-2*V*(2*V-1)*(V**3-V**2+3*V+1)*W
     2   +2*(V-1)*(V**2+1)*(V**2-V+2))+4*V2*(V*W-V-1)*(V**3*W**3-2*V**2*
     3   W**2+V*(V**2+2*V-1)*W-(V-1)*(V**2+1))-2*CQ*V*(V**2+1)*V2*(V*W-V
     4   +1)*(W**2-2*W+2)-CQ*(V**2+1)**2*V1*(V*W-V+1)*(W**2-2*W+2))*LOG(
     5   W)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1))
 
 
      CVC = (3*AL*V1*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**6*W**6-2*V**
     1   5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**5+3*V*
     2   *4+10*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V**3+36*
     3   V**2+4*V+1)*W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2*(V**2
     4   +1)*(V**2+3))+6*V*V2*(V*W-1)*(2*V**5*W**6-V**3*(V+1)*(5*V-1)*W*
     5   *5+V**2*(7*V**3+13*V**2-3*V+1)*W**4-V*(7*V**4+20*V**3-4*V**2-4*
     6   V+1)*W**3+(5*V**5+9*V**4+12*V**3-20*V**2+3*V-1)*W**2-2*(V-1)*(V
     7   **4+2*V**3+4*V**2+4*V-3)*W+2*(V-1)**2*(V**2+3))+3*V1*(V**2*W**2
     8   -2*(V-1)*V*W+(V-1)**2)*(8*V**5*W**5+V**2*(2*V**4-10*V**3-25*V**
     9   2+3*V-2)*W**4+V*(2*V**4+28*V**3+33*V**2-3*V+4)*W**3-(8*V**5+3*V
     :   **4+37*V**3+27*V**2+3*V+2)*W**2+(16*V**4-11*V**3+39*V**2+9*V+3)
     ;   *W-8*V*(V**2-V+2))+4*CQ*V1*(V*W-1)*(V*W-V+1)*(V**5*(V+2)*W**5-V
     <   **4*(V**2+5*V-3)*W**4+V**2*(V**4+7*V**3-V**2-4*V+1)*W**3-V**2*(
     =   V**4+3*V**3+4*V**2-5*V+1)*W**2+(V-1)*(4*V**4-6*V**3+13*V**2-4*V
     >   +1)*W-(V-1)**2*(3*V**2-6*V+7))-6*AL*V2*(V*W-1)*(V*W-V-1)*(V**2*
     ?   W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V
     @   *(V**2+3)*W-2*(V**2+1))-8*CQ*V2*(V**2*W**2-2*V*W+1)*(V**2*W**2-
     1   2*(V-1)*V*W+(V-1)**2)*((V-1)*V**2*W**2-V*(V**2+2*V-1)*W+2*(V-1)
     2   ))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)/3.D0
 
 
      LM = LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**6
     1   -2*V**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**
     2   5+3*V**4+10*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V*
     3   *3+36*V**2+4*V+1)*W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2
     4   *(V**2+1)*(V**2+3))-2*V2*(V*W-1)*(V*W-V-1)*(2*V**3*W**3-V*(V**2
     5   +4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1)))/((V-1)**2*V**2*W**2*(V
     6   *W-1)**2*(V*W-V+1)**2)
 
 
      LMP = 2*LOG(S/MP**2)*V1*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+V*
     1   *2+1)*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*
     2   W+(V-1)**4)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
      STRUV4=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV5(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V1*(-2*V**6*W**6+2
     1   *V**5*(2*V-1)*W**5-V**4*(6*V**2-8*V+3)*W**4+2*(V-1)*V**3*(8*V**
     2   2-9*V+2)*W**3-2*(V-1)**2*V**2*(11*V**2-13*V+5)*W**2+2*(V-1)**3*
     3   V*(6*V**2-8*V+3)*W-(V-1)**4*(2*V**2-2*V+1))+CQ*V1*(2*V**2*W**2-
     4   2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-
     5   1)**4)-2*(V-1)*V*V2*W*(3*V**3*W**3-7*(V-1)*V**2*W**2+(V-1)*V*(V
     6   +1)*W+3*(V-1)**3)+2*CQ*(V-1)*V*V2*W*(V**2*W**2+(V-1)**2)*(2*V**
     7   2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))/((V-1)*V*W*(V*W-1)**3*(V*W-
     8   V+1)**3)
 
 
      LV1 = -LOG(1-V)*(V**2*W**2-2*V*W+1)*(V**3*W**3-3*(V-1)*V**2*W**2+3
     1   *(V-1)**2*V*W-(V-1)**3)*(V1*(6*V**4*W**4-12*V**3*W**3+V**2*(6*V
     2   **2-4*V+7)*W**2-V*(20*V**3-48*V**2+50*V-19)*W+(V-1)*(20*V**2-33
     3   *V+18))+2*(V-1)*V2*(V*W-1)*(2*V**2*W**2-5*V*W+10*V**2-17*V+10))
     4   /((V-1)*(V*W-1)**3*(V*W-V+1)**4)
 
 
      LV = LOG(V)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V1*(-2*V**6*W**6+2*V*
     1   *5*(2*V-1)*W**5-V**3*(8*V**3-14*V**2+8*V-1)*W**4+(V-1)*V**2*(16
     2   *V**3-26*V**2+14*V-3)*W**3-(V-1)**2*V*(18*V**3-26*V**2+16*V-3)*
     3   W**2+(V-1)**3*(12*V**3-16*V**2+8*V-1)*W-2*(V-1)**4*(2*V**2-2*V+
     4   1))-2*(V-1)*V*V2*W*(2*V**4*W**4-2*V**3*(2*V-7)*W**3-V**2*(12*V*
     5   *2+2*V-11)*W**2+2*(V-1)*V*(14*V**2-19*V+9)*W-(V-1)**2*(14*V**2-
     6   26*V+15))+CQ*V1*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*
     7   W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)+2*CQ*(V-1)*V*V2*W*(V**2*W**
     8   2+(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))/((V-1)*V*
     9   W*(V*W-1)**3*(V*W-V+1)**3)
 
 
      LVW = 2*(V1*(2*V**2*W**2+2*(V-2)*V*W+4*V**2-8*V+5)+2*(V-1)*V2)*(V*
     1   *3*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-
     2   1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*LOG(1-V*W)/((V-1)*(V*W
     3   -1)**3*(V*W-V+1)**4)
 
 
      LTVW = 2*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(CQ*V1*(2*V**2*W**2-2*V*(
     1   2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4
     2   )+(V-1)*V*V2*W*(9*V**3*W**3-V**2*(8*V**2+3*V-9)*W**2+(V-1)**2*V
     3   *(16*V+3)*W-(V-1)**2*(8*V**2-7*V+1))+(V-1)*V*V1*W*(5*V**3*W**3-
     4   V**2*(8*V**2-9*V+5)*W**2+(V-1)*V*(16*V**2-25*V+13)*W-(V-1)**2*(
     5   8*V**2-11*V+7))+2*CQ*(V-1)*V*V2*W*(V**2*W**2+(V-1)**2)*(2*V**2*
     6   W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))*LOG(V*W-V+1)/((V-1)*V*W*(V*W-
     7   1)**3*(V*W-V+1)**3)
 
 
      LW = -(V1*(V*(6*V**2-6*V+1)*W**2-(2*V-1)*(2*V**2-2*V+1)*W+3*(V-1)*
     1   (2*V**2-2*V+1))+2*V*V2*W*(2*V**2*W**2+10*V*W-14*V**2+20*V-9))*(
     2   V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(
     3   V-1)**2*V*W-(V-1)**3)*LOG(W)/(V*W*(V*W-1)**3*(V*W-V+1)**4)
 
 
      CVC = (V1*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(6*V**7*W**7-4*V**6*(2*V+1
     1   )*W**6-(V-2)*V**4*(18*V**2-14*V+3)*W**5+V**3*(52*V**4-150*V**3+
     2   162*V**2-77*V+12)*W**4-2*(V-1)*V**3*(21*V**3-43*V**2+22*V+2)*W*
     3   *3+2*(V-1)**2*V*(2*V**4+12*V**3-32*V**2+17*V-6)*W**2+(V-1)**3*(
     4   6*V**4-16*V**3+19*V**2-11*V+6)*W-(V-1)**4*(6*V**2-6*V+1))+4*(V-
     5   1)*V2*W*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**6*W**5-V**4*(4*V**2
     6   +V-2)*W**4+(V-1)*V**3*(4*V**2+3*V-4)*W**3-V**3*(4*V**3-5*V**2-6
     7   *V+5)*W**2+(V-1)*V*(2*V**4+V**3-11*V+4)*W-2*(V-1)**2*(V+1)*(V**
     8   2-V+1)))/((V-1)*V*W*(V*W-1)**3*(V*W-V+1)**4)/2.D0
 
 
      LM = -LOG(S/M**2)*V1*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(
     1   V-1)**3)*(2*V**4*W**4+2*(V-3)*V**3*W**3+V*(6*V**3-14*V**2+12*V-
     2   1)*W**2-(2*V**4+2*V**3-8*V**2+6*V-1)*W+(V-1)*(2*V**2-2*V+1))/((
     3   V-1)*V*W*(V*W-1)*(V*W-V+1)**3)
 
 
      LMP = -LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V1*(2*V**6*W**6-2*V**5*(2
     1   *V-1)*W**5+V**4*(6*V**2-9*V+4)*W**4-(V-1)*V**3*(12*V**2-15*V+4)
     2   *W**3+(V-1)**2*V**2*(14*V**2-19*V+7)*W**2-(V-1)**3*V*(8*V**2-13
     3   *V+6)*W+(V-1)**4*(2*V**2-2*V+1))+2*(V-1)*V*V2*W*(V**2*W**2+(V-1
     4   )**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))/((V-1)*V*W*(V*W
     5   -1)**2*(V*W-V+1)**3)
 
      STRUV5=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV6(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(CQ*V1*(2*V**8*W**7-2*V**7*(V+1)*W**6+V**3*(2*V**5-
     1   5*V**4+17*V**3-16*V**2+12*V-4)*W**5-V**2*(2*V**6-7*V**5+13*V**4
     2   +8*V**3-24*V**2+28*V-12)*W**4+V*(2*V**7-4*V**6-V**5+25*V**4-32*
     3   V**3+16*V**2+12*V-12)*W**3-(6*V**7-24*V**6+43*V**5-27*V**4-10*V
     4   **3+32*V**2-12*V-4)*W**2+(6*V**6-28*V**5+63*V**4-72*V**3+41*V**
     5   2-8)*W-2*(V-1)*(V**2-2*V+2)**2)+V1*(-2*V**8*W**7+2*V**7*(V+1)*W
     6   **6-V**3*(2*V**5+3*V**4+7*V**3-14*V**2+12*V-4)*W**5+V**2*(2*V**
     7   6+V**5+5*V**4+30*V**3-46*V**2+36*V-12)*W**4-V*(2*V**7-4*V**6+21
     8   *V**5-39*V**4+80*V**3-62*V**2+36*V-12)*W**3+(6*V**7-24*V**6+63*
     9   V**5-99*V**4+108*V**3-50*V**2+12*V-4)*W**2-V*(6*V**5-28*V**4+69
     :   *V**3-96*V**2+79*V-28)*W+2*(V-1)*(V**2-2*V+2)**2)+2*CQ*V2*(2*V*
     ;   *7*W**6-V**3*(3*V**4-2*V**3+7*V**2-6*V+2)*W**5+V**2*(3*V**5-5*V
     <   **4+5*V**3+9*V**2-14*V+6)*W**4-V*(2*V**6-2*V**5-5*V**4+20*V**3-
     =   11*V**2-6*V+6)*W**3+(6*V**6-17*V**5+19*V**4+3*V**3-17*V**2+6*V+
     >   2)*W**2-2*(V-1)**2*(3*V+1)*(V**2-2*V+2)*W+2*(V-1)**2*(V**2-2*V+
     ?   2))+2*(V-1)*V2*W*(V**3*(3*V**3-2*V**2+4*V-2)*W**4-V**2*(3*V**4-
     @   14*V**3+8*V**2+8*V-6)*W**3+V*(7*V**4-44*V**3+36*V**2-6)*W**2-(9
     1   *V**4-42*V**3+40*V**2-8*V-2)*W+(V-1)*(3*V**2-10*V+4))+2*(V-1)*V
     2   *V4*W*(V**2*W**2-2*V*W+1)*(2*(V-1)*V**3*W**3-4*V**4*W**2+V*(2*V
     3   **3+2*V**2-2*V+7)*W-4*V**2+7*V-5)+2*(V-1)*V*V3*W*(V**2*W**2-2*V
     4   *W+1)*(2*V**3*W**3-2*(V-1)*V**2*W**2+V*(2*V**2+1)*W-2*V**3+6*V*
     5   *2-9*V+3))/((V-1)**2*V**3*W**2*(V*W-1)**3)
 
 
      LV1 = -LOG(1-V)*(2*V1*(2*V**8*W**7-2*V**7*(V+1)*W**6+V**5*(V**3+7*
     1   V**2-3*V+1)*W**5-V**3*(9*V**4-9*V**3+35*V**2-31*V+10)*W**4+V**2
     2   *(23*V**4-60*V**3+130*V**2-105*V+30)*W**3-V*(28*V**4-87*V**3+15
     3   8*V**2-119*V+30)*W**2+(15*V**4-48*V**3+77*V**2-52*V+10)*W-3*(V-
     4   1)*(V**2-2*V+2))+2*V2*(V**2*W**2-2*V*W+1)*(2*V**5*W**4-3*(V-1)*
     5   *2*V**3*W**3+V*(6*V**4-11*V**3+20*V**2-21*V+8)*W**2-(V-1)**2*(2
     6   *V**3+V**2-4*V+8)*W+2*(V-1)**2*(V**2-2*V+2))+2*(V-1)*V*V4*W*(V*
     7   W-1)*(2*V**5*W**4-2*(V-1)*V**3*(2*V+1)*W**3+V**2*(2*V**3+V-14)*
     8   W**2-V*(2*V**3-6*V**2+11*V-22)*W-2*(2*V**2-4*V+5))+CQ*V**2*V1*W
     9   *(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+2*V**2*W**2+1)+2*C
     :   Q*V**3*V2*W**2*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)+2
     ;   *(V-1)*V*V3*W*(V*(2*V**2-3*V+4)*W-2*(V-1)*(V**2-2*V+2))*(V**2*W
     <   **2-2*V*W+1))/((V-1)**2*V**3*W**2*(V*W-1)**3)
 
 
      LV = LOG(V)*(V1*(-2*V**8*W**7+2*V**7*(V+1)*W**6-V**3*(2*V**5+9*V**
     1   4+V**3-8*V**2+6*V-2)*W**5+V**2*(2*V**6-V**5+35*V**4-4*V**3-6*V*
     2   *2+12*V-6)*W**4-V*(2*V**7-4*V**6+13*V**5+23*V**4+22*V**2-6)*W**
     3   3+(6*V**7-24*V**6+53*V**5-41*V**4+26*V**3+22*V**2-12*V-2)*W**2-
     4   (6*V**6-28*V**5+65*V**4-76*V**3+49*V**2-6*V-6)*W+2*(V-1)*(V**2-
     5   2*V+2)**2)-2*V2*(2*V**7*W**6-V**3*(9*V**4-8*V**3+7*V**2-6*V+2)*
     6   W**5+V**2*(5*V**5-3*V**4+25*V**3-15*V**2-6*V+6)*W**4+V*(2*V**6-
     7   18*V**5+55*V**4-126*V**3+91*V**2-18*V-6)*W**3-(6*V**6-33*V**5+9
     8   1*V**4-155*V**3+113*V**2-30*V-2)*W**2+2*(V-1)**2*(3*V**3-7*V**2
     9   +12*V-6)*W-2*(V-1)**2*(V**2-2*V+2))+CQ*V**2*V1*W*(2*V**2*W**2-2
     :   *V*(V+1)*W+V**2+1)*(V**4*W**4+2*V**2*W**2+1)+4*(V-1)*V*V3*W*(V*
     ;   *2*W**2-2*V*W+1)*(V**3*W**2-4*V*W-(V-2)*(V**2-2*V+2))-4*(V-1)*V
     <   *V4*W*(V**2*W**2-2*V*W+1)*(V*(V**2+1)*W**2-(5*V+1)*W-V**3+4*V**
     =   2-6*V+5)+2*CQ*V**3*V2*W**2*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)
     >   *W+V**2+1))/((V-1)**2*V**3*W**2*(V*W-1)**3)
 
 
      LVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V1*(2*V**5*W*
     1   *4-2*(V-2)*V**4*W**3+V**3*(V**2+2*V+3)*W**2-(3*V**4-8*V**3+17*V
     2   **2-14*V+4)*W+2*(V-1)*(V**2-2*V+2))+V*V2*W*(4*V**3*W**2-(V-5)*V
     3   **2*W-2*(V-1)**2)+(V-1)*V*V4*W*(4*V**2*W-4*V**2+7*V-7)+(V-1)*V*
     4   V3*W*(2*V**2*W-2*V**2+3*V-5))*LOG(1-V*W)/((V-1)**2*V**3*W**2*(V
     5   *W-1)**4)
 
 
      LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V1*(2*V**5*W
     1   **4-2*(V-2)*V**4*W**3+V**3*(2*V**2-3*V+5)*W**2-(V-1)**2*(2*V**3
     2   -V**2-4*V+4)*W+2*(V-1)**3*(V**2-2*V+2))+4*(V-1)*V*V4*W*(V**3*W*
     3   *2-V**2*(2*V-3)*W+(V-1)*(V**2-2*V+2))-(V-1)*V*V2*W*(V**2*W-(V-4
     4   )*(V-1)))*LOG(V*W-V+1)/((V-1)**2*V**3*W**2*(V*W-1)**4)
 
 
      LW = -(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V2*(2*V**5*W**4-(V-1)*V**
     1   3*(6*V-1)*W**3+V**2*(3*V**3+2*V**2+4*V+1)*W**2-2*(V-1)**2*V*(V*
     2   *2-2*V+2)*W+2*(V-1)**2*(V**2-2*V+2))+2*V1*(V*(V**5+4*V**4-3*V**
     3   3+3*V**2-3*V+1)*W**3-(2*V**6-2*V**5+2*V**4+17*V**3-12*V**2+2*V+
     4   1)*W**2+(2*V**6-8*V**5+20*V**4-30*V**3+37*V**2-23*V+5)*W-2*(V-1
     5   )*(V**2-2*V+2)**2)+2*(V-1)*V*V4*W*(2*(V-1)*V**3*W**3-2*V*(2*V**
     6   3-3*V**2-1)*W**2+(2*V**4-2*V**3+V**2-10*V-2)*W-2*(V**3-4*V**2+6
     7   *V-6))+2*(V-1)*V**2*V3*W**2*(2*V**2*W**2-2*(V-1)*V*W+V+4)-2*CQ*
     8   (V-1)**2*(V**2-2*V+2)*V2*(V*W-1)*(W**2-2*W+2)+CQ*(V-1)*(V**2-2*
     9   V+2)**2*V1*(V*W-1)*(W**2-2*W+2))*LOG(W)/((V-1)**2*V**3*W**2*(V*
     :   W-1)**4)
 
 
      CVC = (3*(V-1)*V1*(V*W-1)*(6*V**9*W**8-4*V**8*(3*V+1)*W**7+V**4*(8
     1   *V**5+V**4+19*V**3-62*V**2+64*V-24)*W**6-(V-2)*V**3*(8*V**5-22*
     2   V**4+75*V**3-35*V**2-22*V+24)*W**5+V**3*(12*V**6-43*V**5+55*V**
     3   4+62*V**3-344*V**2+368*V-144)*W**4-V*(6*V**8+4*V**7-93*V**6+245
     4   *V**5-304*V**4+90*V**3+200*V**2-200*V+48)*W**3+(18*V**8-60*V**7
     5   -3*V**6+289*V**5-551*V**4+503*V**3-198*V**2-16*V+24)*W**2-(V-1)
     6   *(18*V**6-66*V**5+59*V**4+56*V**3-119*V**2+100*V-36)*W+2*(V-1)*
     7   *2*(3*V**4-10*V**3+8*V**2+2))+12*(V-1)*V2*(V*W-1)*(2*V**8*W**7-
     8   V**4*(10*V**4-23*V**3+35*V**2-24*V+6)*W**6+V**3*(12*V**5-25*V**
     9   4+13*V**3+28*V**2-36*V+12)*W**5-V**3*(10*V**5-25*V**4-11*V**3+9
     :   8*V**2-98*V+30)*W**4+V*(2*V**7+3*V**6-27*V**5+24*V**4+62*V**3-1
     ;   16*V**2+66*V-12)*W**3-(V-1)*(6*V**6-10*V**5-3*V**4+22*V**3-7*V*
     <   *2-12*V+6)*W**2+(V-1)**2*(V**2-2*V+2)*(6*V**2-V-3)*W-2*(V-1)**3
     =   *(V**2-2*V+2))-6*AL*(V-1)*V1*(V*W-1)*(V*W-V+1)*(2*V**8*W**7-2*V
     >   **7*(V+1)*W**6+V**3*(2*V**5-5*V**4+17*V**3-16*V**2+12*V-4)*W**5
     ?   -V**2*(2*V**6-7*V**5+13*V**4+8*V**3-24*V**2+28*V-12)*W**4+V*(2*
     @   V**7-4*V**6-V**5+25*V**4-32*V**3+16*V**2+12*V-12)*W**3-(6*V**7-
     1   24*V**6+43*V**5-27*V**4-10*V**3+32*V**2-12*V-4)*W**2+(6*V**6-28
     2   *V**5+63*V**4-72*V**3+41*V**2-8)*W-2*(V-1)*(V**2-2*V+2)**2)-12*
     3   AL*(V-1)*V2*(V*W-1)*(V*W-V+1)*(2*V**7*W**6-V**3*(3*V**4-2*V**3+
     4   7*V**2-6*V+2)*W**5+V**2*(3*V**5-5*V**4+5*V**3+9*V**2-14*V+6)*W*
     5   *4-V*(2*V**6-2*V**5-5*V**4+20*V**3-11*V**2-6*V+6)*W**3+(6*V**6-
     6   17*V**5+19*V**4+3*V**3-17*V**2+6*V+2)*W**2-2*(V-1)**2*(3*V+1)*(
     7   V**2-2*V+2)*W+2*(V-1)**2*(V**2-2*V+2))-8*CQ*(V-1)**2*V1*W*(V*W-
     8   V+1)*(V**2*W**2-2*V*W+1)*(V**6*W**4+V**2*(V**4-4*V**3+10*V**2-8
     9   *V+4)*W**2-2*V*(V**2-2*V+2)**2*W+V**4-4*V**3+9*V**2-8*V+4)+12*(
     :   V-1)**2*V*V4*W*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**5*W**4-V**2*(2
     ;   *V**3+2*V**2+V+4)*W**3+V*(V**4+4*V**3-6*V**2+7*V+8)*W**2-(2*V**
     <   4+3*V**3-16*V**2+11*V+4)*W+V**3+5*V**2-13*V+5)-16*CQ*(V-1)**2*V
     =   2*W*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**5*W**3-(V-1)*V**2*(V**2-2
     >   *V+2)*W**2+V*(2*V**3-5*V**2+8*V-4)*W-(V-1)*(V**2-2*V+2))-12*(V-
     ?   1)**3*V*V3*W**2*(V*W-V+1)*(V**2*(4*V-3)*W**2-2*V*(2*V**2+4*V-3)
     @   *W+V**2+4*V-3)*(V**2*W**2-2*V*W+1))/((V-1)**3*V**3*W**2*(V*W-1)
     1   **4*(V*W-V+1))/6.D0
 
 
      LM = LOG(S/M**2)*(-V1*(2*V**8*W**7-2*V**7*(V+1)*W**6+V**3*(2*V**5-
     1   V**4+11*V**3-8*V**2+6*V-2)*W**5-V**2*(2*V-1)*(V**5-V**4+8*V**3+
     2   10*V**2-4*V+6)*W**4+V*(2*V**7-4*V**6+9*V**5+11*V**4+20*V**3-18*
     3   V**2+12*V-6)*W**3-(6*V**7-24*V**6+57*V**5-63*V**4+68*V**3-26*V*
     4   *2-2)*W**2+(6*V**6-28*V**5+71*V**4-96*V**3+83*V**2-30*V-2)*W-2*
     5   (V-1)*(V**2-2*V+2)*(V**2-2*V+3))-2*V2*(2*V**7*W**6-V**3*(3*V**4
     6   -2*V**3+7*V**2-6*V+2)*W**5+V**2*(3*V**5-5*V**4+5*V**3+9*V**2-14
     7   *V+6)*W**4-V*(2*V**6-2*V**5-5*V**4+20*V**3-11*V**2-6*V+6)*W**3+
     8   (6*V**6-17*V**5+19*V**4+3*V**3-17*V**2+6*V+2)*W**2-2*(V-1)**2*(
     9   3*V+1)*(V**2-2*V+2)*W+2*(V-1)**2*(V**2-2*V+2))-4*(V-1)*V*V4*W*(
     :   V*W-1)*(V**2*W**3-V*(3*V+2)*W**2-(V**2-6*V-1)*W-V**2+2*V-3))/((
     ;   V-1)**2*V**3*W**2*(V*W-1)**3)
 
 
      LMP = -2*LOG(S/MP**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*
     1   (V1*(V**6*W**5-V**5*(2*V-3)*W**4+V**4*(2*V**2-5*V+5)*W**3-(V-1)
     2   *V*(V+1)*(2*V-1)*(V**2-2*V+2)*W**2+(V-1)**3*(2*V**3-3*V**2+2)*W
     3   -(V-1)**4*(V**2-2*V+2))+2*(V-1)*V*V4*W*(V*W-V+1)*(V**3*W**2-V**
     4   2*(2*V-3)*W+(V-1)*(V**2-2*V+2)))/((V-1)**2*V**3*W**2*(V*W-1)**4
     5   *(V*W-V+1))
 
      STRUV6=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV7(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+V**4*(16*V
     1   **4+24*V**3+11*V**2+1)*W**6-2*V**4*(8*V**4+12*V**3+21*V**2-V+2)
     2   *W**5+2*V**2*(8*V**6+4*V**5+29*V**4+9*V**3-6*V**2+V-1)*W**4-2*V
     3   **2*(6*V**6+2*V**5+11*V**4+30*V**3-16*V**2-2*V-3)*W**3+(4*V**8+
     4   12*V**7-17*V**6+54*V**5-5*V**4-32*V**3+V**2-2*V+1)*W**2-2*(V-1)
     5   *(4*V**6-2*V**5+7*V**4+6*V**3+6*V**2-4*V-1)*W+4*(V-1)**2*(V**2+
     6   1)*(V**2-V+2))-V1*(-4*V**8*W**8+4*V**7*(3*V+1)*W**7-V**4*(16*V*
     7   *4+2*V**3+17*V**2+1)*W**6+2*V**4*(8*V**4-11*V**3+25*V**2+4*V+2)
     8   *W**5-2*V**2*(8*V**6-19*V**5+20*V**4+33*V**3-10*V**2+V-1)*W**4+
     9   2*V**2*(6*V**6-9*V**5-11*V**4+63*V**3-16*V**2-6*V-3)*W**3-(4*V*
     :   *8+8*V**7-49*V**6+88*V**5+13*V**4-46*V**3-V**2-2*V+1)*W**2+2*(V
     ;   -1)*(4*V**6-6*V**5+4*V**4+11*V**3+9*V**2-5*V-1)*W-4*(V-1)**2*(V
     <   **2+1)*(V**2-2*V+3))+2*V2*(V*W-1)*(V*W-V+1)*(V**3*(6*V**2-3*V+1
     =   )*W**4-V**3*(8*V**2+9*V-5)*W**3+V*(6*V**4+9*V**3+15*V**2-13*V-1
     >   )*W**2-V*(9*V**3+V**2+11*V-13)*W+4*(V-1)*(V**2+1))+2*CQ*V2*(V*W
     ?   -1)*(V*W-V-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**3*W**3-V*(
     @   V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1))-4*CQ*(V-1)*V*V4*W*(
     1   V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-
     2   2*V**2*W+V**2+1)+4*(V-1)*V4*W*(V*W-1)*(V*W-V-1)*(V**2*W**2-2*(V
     3   -1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)+4*(V-1)*V*V3*W*(V
     4   *W-1)*(V*W-V+1)*(V*W**2-W+V-1)*(V**2*W**2-2*V**2*W+V**2+1))/((V
     5   -1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
 
      LV1 = LOG(1-V)*(2*V1*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**4*W**4-2*
     1   V**3*(V+2)*W**3+V**2*(V**2+2*V+9)*W**2-2*V*(V**2+2*V+3)*W+2*(V+
     2   1)*(V**2+1))+CQ*V1*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(
     3   V**4*W**4-4*V**3*W**3+8*V**2*W**2-8*V*W+4)-4*V2*(V*W-V-1)*(V**2
     4   *W**2-2*V*W+1)*(V**3*W**3-V**2*(V+1)*W**2+V*(3*V**2-2*V+1)*W-(V
     5   -1)*(V**2+1))-2*CQ*V2*(V*W-1)*(V*W-V+1)*(V**2*W**2-2*V*W+2)*(2*
     6   V**2*W**2-2*V*(V+1)*W+V**2+1)+4*(V-1)*V**2*V4*W**2*(V*W-1)*(V*W
     7   -V+1)*(V**2*W**2-2*V**2*W+V**2+1)+4*(V-1)**2*V3*W*(V*W-1)*(V*W-
     8   V-1)*(V**2*W**2-2*V**2*W+V**2+1))/((V-1)**2*V**2*W**2*(V*W-1)**
     9   2*(V*W-V+1))
 
 
      LV = LOG(V)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+3*V**6*(V+1)*
     1   (5*V+3)*W**6-2*V**5*(6*V**3+12*V**2+17*V-1)*W**5+V**4*(9*V**4+6
     2   *V**3+46*V**2+14*V-15)*W**4-2*V**3*(3*V+5)*(V**4-2*V**3+6*V**2-
     3   2*V-1)*W**3+2*V**2*(V**6+2*V**5-7*V**4+20*V**3+3*V**2-18*V+3)*W
     4   **2-2*(V-1)*V*(2*V**5-3*V**4+4*V**3+4*V**2+6*V-5)*W+2*(V-1)**2*
     5   (V**2+1)*(V**2-2*V+3))-V1*(-4*V**8*W**8+4*V**7*(3*V+2)*W**7-V**
     6   4*(16*V**4+12*V**3+23*V**2+1)*W**6+2*V**4*(8*V**4-8*V**3+37*V**
     7   2+3*V+2)*W**5-2*V**2*(8*V**6-22*V**5+35*V**4+37*V**3-14*V**2+V-
     8   1)*W**4+2*V**2*(6*V**6-14*V**5-5*V**4+74*V**3-20*V**2-10*V-3)*W
     9   **3-(4*V**8+4*V**7-57*V**6+114*V**5+7*V**4-52*V**3-3*V**2-2*V+1
     :   )*W**2+2*(V-1)*(4*V**6-10*V**5+7*V**4+10*V**3+14*V**2-8*V-1)*W-
     ;   4*(V-1)**2*(V**2+1)*(V**2-3*V+4))-2*V2*(V*W-1)*(V*W-V+1)*(2*V**
     <   5*W**5-V**3*(11*V**2-2*V+1)*W**4+2*V**3*(4*V**2+11*V-5)*W**3-V*
     =   (V**4+16*V**3+24*V**2-20*V-1)*W**2-2*V*(V**4-4*V**3-2*V**2-8*V+
     >   9)*W+2*(V-3)*(V-1)*(V**2+1))+2*CQ*V2*(V*W-1)*(V**2*W**2-2*V*W+2
     ?   )*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**2*W**2-2*V*(V+1)*W+V**
     @   2+1)-4*CQ*(V-1)*V*V4*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V
     1   *W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)+4*(V-1)*V*V4*W*(V*W-2)
     2   *(V*W-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V
     3   **2+1)+4*(V-1)*V*V3*W*(V*W-1)*(V*W-V+1)*(V*W+V-1)*(V**2*W**2-2*
     4   V**2*W+V**2+1))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
 
      LVW = -2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*
     1   (2*V**4*W**4-2*V**3*(V+3)*W**3+V**2*(V**2+3*V+12)*W**2-V*(V+3)*
     2   *2*W+4*(V**2+1))+V2*(-4*V**3*W**3+V**2*(5*V+7)*W**2-V*(3*V**2+6
     3   *V+7)*W+4*(V**2+1)))*LOG(1-V*W)/((V-1)**2*V**2*W**2*(V*W-1)**2*
     4   (V*W-V+1)**2)
 
 
      LTVW = -4*(V**2*W**2-2*V*W+1)*(CQ*V1*(V**2*W**2-2*V**2*W+V**2+1)*(
     1   V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)
     2   **4)-(V-1)*V*V1*W*(V**3*W**3-V**2*(5*V-3)*W**2+(V-1)*V*(5*V+1)*
     3   W-(V-1)**2*(V+1))+2*CQ*(V-1)*V*V4*W*(V**2*W**2-2*(V-1)*V*W+(V-1
     4   )**2)*(V**2*W**2-2*V**2*W+V**2+1))*LOG(V*W-V+1)/((V-1)**2*V**2*
     5   W**2*(V*W-1)**2*(V*W-V+1)**2)
 
 
      LW = (-2*V1*(V*W-1)*(V*W-V+1)*(2*V**3*W**3-V**2*(V**2+7)*W**2+2*V*
     1   (V**3-V**2+3*V+3)*W-2*(V**2+1)*(V**2-V+2))-4*V2*(V*W-1)*(V*W-V-
     2   1)*(V**3*W**3-V**2*(3*V-1)*W**2+V*(V**2+2*V-1)*W-(V-1)*(V**2+1)
     3   )+4*(V-1)**2*V4*W*(V*W-1)*(V*W-V+1)*(V**2*W**2-2*V**2*W+V**2+1)
     4   -4*(V-1)*V*V3*W**2*(V*W-V-1)*(V**2*W**2-2*V**2*W+V**2+1)+2*CQ*V
     5   *(V**2+1)*V2*(V*W-1)*(V*W-V+1)*(W**2-2*W+2)+CQ*(V**2+1)**2*V1*(
     6   V*W-1)*(V*W-V+1)*(W**2-2*W+2))*LOG(W)/((V-1)**2*V**2*W**2*(V*W-
     7   1)*(V*W-V+1))
 
 
      CVC = (3*V1*(8*V**7*W**7+V**4*(2*V**4-24*V**3-13*V**2+5*V-2)*W**6-
     1   V**4*(4*V**4-32*V**3-46*V**2+17*V-7)*W**5+2*V**2*(V**6-12*V**5-
     2   36*V**4+7*V**3+5*V**2-7*V+2)*W**4+2*V**2*(8*V**5+19*V**4+22*V**
     3   3-27*V**2+11*V-5)*W**3-(8*V**7+17*V**6+3*V**5+40*V**4-66*V**3+2
     4   1*V**2-9*V+2)*W**2+(V-1)*(16*V**5-11*V**4+34*V**3+2*V**2-6*V-3)
     5   *W-8*(V-1)**2*V*(V**2-V+2))+3*AL*V1*(V**2*W**2-2*(V-1)*V*W+(V-1
     6   )**2)*(2*V**6*W**6-2*V**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**
     7   2+1)*W**4-2*V*(V**5+3*V**4+10*V**3+20*V**2+V+1)*W**3+(2*V**6+4*
     8   V**5+13*V**4+24*V**3+36*V**2+4*V+1)*W**2-2*(2*V**5+V**4+8*V**3+
     9   6*V**2+10*V+1)*W+2*(V**2+1)*(V**2+3))+12*V*V2*(V*W-1)*(V**5*W**
     :   6-V**3*(4*V**2-V+1)*W**5+V**2*(5*V**3+2*V**2+3*V-1)*W**4-V*(5*V
     ;   **4-2*V**3+16*V**2-8*V-1)*W**3+(4*V**5-3*V**4+12*V**3-4*V**2-6*
     <   V+1)*W**2-(V-1)*(V**4+2*V**3+4*V**2+4*V-3)*W+(V-1)**2*(V**2+3))
     =   +4*CQ*V1*(V*W-1)*(V*W-V+1)*(V**5*(V+2)*W**5-V**4*(V**2+5*V-3)*W
     >   **4+V**2*(V**4+7*V**3-V**2-4*V+1)*W**3-V**2*(V**4+3*V**3+4*V**2
     ?   -5*V+1)*W**2+(V-1)*(4*V**4-6*V**3+13*V**2-4*V+1)*W-(V-1)**2*(3*
     @   V**2-6*V+7))+12*(V-1)*V*V4*W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(
     1   V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)-6*
     2   AL*V2*(V*W-1)*(V*W-V-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**
     3   3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1))+24*CQ*(V-
     4   1)*V*V4*W*(V**2*W**2-2*V*W+1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V
     5   -1)**2*V*W-(V-1)**3)-8*CQ*V2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(
     6   V-1)*V*W+(V-1)**2)*((V-1)*V**2*W**2-V*(V**2+2*V-1)*W+2*(V-1)))/
     7   ((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)/3.D0
 
 
      LM = LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**6
     1   -2*V**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**
     2   5+3*V**4+10*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V*
     3   *3+36*V**2+4*V+1)*W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2
     4   *(V**2+1)*(V**2+3))-2*V2*(V*W-1)*(V*W-V-1)*(2*V**3*W**3-V*(V**2
     5   +4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1)))/((V-1)**2*V**2*W**2*(V
     6   *W-1)**2*(V*W-V+1)**2)
 
 
      LMP = 2*LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+V**2+
     1   1)*(V1*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V
     2   *W+(V-1)**4)+2*(V-1)*V*V4*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2))/(
     3   (V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
      STRUV7=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV8(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V
     1   -2)*W**7+V**4*(52*V**4-50*V**3+17*V**2-4*V+1)*W**6-2*V**4*(42*V
     2   **4-59*V**3+30*V**2-7*V+2)*W**5+2*V**2*(44*V**6-64*V**5+15*V**4
     3   +25*V**3-16*V**2+5*V-1)*W**4-2*V**2*(28*V**6-22*V**5-60*V**4+11
     4   4*V**3-71*V**2+20*V-3)*W**3+(16*V**8+40*V**7-204*V**6+280*V**5-
     5   150*V**4+12*V**3+15*V**2-6*V+1)*W**2-2*(V-1)*(16*V**6-28*V**5+6
     6   *V**4+30*V**3-29*V**2+10*V-1)*W+4*(V-1)**2*(2*V**2-3*V+2)*(2*V*
     7   *2-2*V+1))+V1*(-4*V**8*W**8+4*V**7*(4*V-1)*W**7-V**4*(36*V**4-3
     8   0*V**3+13*V**2-4*V+1)*W**6+2*V**4*(28*V**4-39*V**3+22*V**2-5*V+
     9   2)*W**5-2*V**2*(32*V**6-46*V**5+V**4+37*V**3-20*V**2+5*V-1)*W**
     :   4+2*V**2*(24*V**6-20*V**5-66*V**4+130*V**3-79*V**2+20*V-3)*W**3
     ;   -(16*V**8+32*V**7-204*V**6+292*V**5-138*V**4-12*V**3+23*V**2-6*
     <   V+1)*W**2+2*(V-1)*(16*V**6-32*V**5+4*V**4+46*V**3-43*V**2+14*V-
     =   1)*W-4*(V-1)**2*(2*V**2-4*V+3)*(2*V**2-2*V+1))-2*(V-1)*V2*(V*W-
     >   1)*(V*W-V+1)*(V**3*(4*V**2-14*V+1)*W**4-4*V**3*(3*V**2-11*V+4)*
     ?   W**3+V*(16*V**4-50*V**3+24*V**2+2*V-1)*W**2-2*(V-1)*V*(4*V**3-4
     @   *V**2-5*V+2)*W+4*(V-1)**2*(2*V**2-2*V+1))-2*CQ*(V-1)*V2*(V*W-2*
     1   V+1)*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**3*W**3-V*(6*V**2-6*V+1
     2   )*W**2+2*V*(4*V**2-6*V+3)*W-2*(V-1)*(2*V**2-2*V+1)))/((V-1)*V**
     3   2*W**2*(V*W-1)**4*(V*W-V+1)**2)
 
 
      LV1 = -2*LOG(1-V)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)
     1   **2)*(CQ*V1*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-V**3*W
     2   **3+V**2*W**2-V*W+1)+2*V*V1*(W-1)*(V*W-1)*(V**4*W**4-V**3*(2*V-
     3   1)*W**3+V**2*(2*V**2-3*V+2)*W**2-V**2*W+2*V**2-2*V+1)-2*(V-1)*V
     4   2*(V*W-1)*(V**4*W**4-2*V**3*(V+1)*W**3+V**3*(2*V+3)*W**2-V*(8*V
     5   **2-9*V+5)*W+2*V**2-2*V+1))/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+
     6   1)**2)
 
 
      LV = LOG(V)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V-2)
     1   *W**7+3*V**6*(2*V-1)*(8*V-3)*W**6-2*V**5*(34*V**3-43*V**2+14*V+
     2   1)*W**5+V**4*(60*V**4-80*V**3-2*V**2+46*V-15)*W**4-2*V**3*(8*V-
     3   5)*(2*V**4-6*V**2+6*V-1)*W**3+2*V**2*(4*V**6+12*V**5-64*V**4+88
     4   *V**3-42*V**2+3)*W**2-2*(V-1)*V*(8*V**5-16*V**4+2*V**3+22*V**2-
     5   19*V+5)*W+2*(V-1)**2*(2*V**2-4*V+3)*(2*V**2-2*V+1))+V1*(-4*V**8
     6   *W**8+4*V**7*(5*V-2)*W**7-V**4*(52*V**4-50*V**3+17*V**2-4*V+1)*
     7   W**6+2*V**4*(42*V**4-59*V**3+24*V**2-V+2)*W**5-2*V**2*(44*V**6-
     8   64*V**5-9*V**4+59*V**3-26*V**2+5*V-1)*W**4+2*V**2*(28*V**6-22*V
     9   **5-98*V**4+180*V**3-101*V**2+22*V-3)*W**3-(16*V**8+40*V**7-260
     :   *V**6+380*V**5-178*V**4-16*V**3+27*V**2-6*V+1)*W**2+2*(V-1)*(16
     ;   *V**6-36*V**5+2*V**4+62*V**3-57*V**2+18*V-1)*W-4*(V-1)**2*(2*V*
     <   *2-5*V+4)*(2*V**2-2*V+1))+2*(V-1)*V2*(V*W-1)*(V*W-V+1)*(2*V**5*
     =   W**5-V**3*(10*V**2-18*V+1)*W**4+4*V**3*(5*V**2-14*V+5)*W**3-V*(
     >   20*V**4-60*V**3+30*V**2-1)*W**2+2*V*(4*V**4-8*V**3-4*V**2+10*V-
     ?   3)*W-2*(V-1)*(2*V-3)*(2*V**2-2*V+1))-2*CQ*(V-1)*V2*(V*W-V+1)*(V
     @   **2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W*
     1   *2-2*V*(2*V-1)*W+2*V**2-2*V+1))/((V-1)*V**2*W**2*(V*W-1)**4*(V*
     2   W-V+1)**2)
 
 
      LVW = 4*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*
     1   V**2*W**2-4*V*W+1)*(V1*(V**4*W**4-V**3*(2*V-1)*W**3+V**2*(2*V**
     2   2-4*V+3)*W**2+V*(2*V**2-4*V+1)*W+2*V**2-2*V+1)-3*(V-1)*V*V2*W*(
     3   V*W-2*V+1))*LOG(1-V*W)/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2
     4   )
 
 
      LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*V1*(2*V**
     1   2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3
     2   +8*(V-1)**2*V**2*W**2-8*(V-1)**3*V*W+4*(V-1)**4)-2*CQ*(V-1)*V2*
     3   (V*W-V+1)*(V**2*W**2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W**2-2*V*(
     4   2*V-1)*W+2*V**2-2*V+1)-2*(V-1)*V*V2*W*(V*W-V+1)*(V**2*W**2-V*(2
     5   *V-1)*W+(V-1)**2)+(V-1)**2*V**2*V1*W*((2*V-1)*W-2*(V-1)))*LOG(V
     6   *W-V+1)/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
 
 
      LW = (V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W
     1   -1)*(2*V1*(2*(V-1)*V**4*W**4-V**3*(2*V-1)*(4*V-3)*W**3+V**2*(12
     2   *V**3-16*V**2+4*V+3)*W**2-2*V*(V+1)*(4*V**3-8*V**2+6*V-1)*W+2*(
     3   2*V**2-3*V+2)*(2*V**2-2*V+1))+4*(V-1)*V2*(V*W-2*V+1)*(V**3*W**3
     4   -2*(V-1)*V**2*W**2+V*(2*V**2-1)*W-2*V**2+2*V-1)-2*CQ*(V-1)*V*(2
     5   *V**2-2*V+1)*V2*(V*W-1)*(W**2-2*W+2)-CQ*(2*V**2-2*V+1)**2*V1*(V
     6   *W-1)*(W**2-2*W+2))*LOG(W)/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1
     7   )**2)
 
 
      CVC = (3*AL*V1*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2
     1   )*(2*V**6*W**6-2*V**5*(2*V+1)*W**5+V**2*(8*V**3-4*V**2+4*V-1)*W
     2   **4+2*V*(2*V**2-V+1)*(2*V**3-2*V**2-2*V+1)*W**3-(8*V**6-16*V**4
     3   +16*V**3-10*V**2+1)*W**2+2*(V-1)*(8*V**4-4*V**3+2*V**2+2*V-1)*W
     4   -4*(V-1)*V*(2*V**2-2*V+1))-6*(V-1)*V*V2*(V*W-V+1)*(V**2*W**2-2*
     5   V*W+1)*(2*V**5*W**6-V**3*(2*V-1)*(4*V+1)*W**5+V**2*(18*V**3-10*
     6   V**2-1)*W**4-V*(20*V**4-4*V**3-10*V**2+1)*W**3+(8*V**5+20*V**4-
     7   40*V**3+18*V**2-2*V+1)*W**2-2*(8*V**4-10*V**3-2*V**2+8*V-3)*W+2
     8   *(V-1)*(4*V**2-6*V+3))-3*V1*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-
     9   4*V*W+1)*(8*(V-1)*V**5*W**5-V**2*(32*V**4-59*V**3+28*V**2-5*V+2
     :   )*W**4+(V-1)*V*(64*V**4-101*V**3+48*V**2-13*V+4)*W**3-(V-1)*(80
     ;   *V**5-180*V**4+156*V**3-59*V**2+13*V-2)*W**2+(V-1)**2*(56*V**4-
     <   106*V**3+84*V**2-21*V+3)*W-8*(V-1)**3*V*(2*V**2-3*V+2))-6*AL*(V
     =   -1)*V*(2*V**2-2*V+1)*V2*(W**2-2*W+2)*(V**2*W**2-2*(V-1)*V*W+(V-
     >   1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)-8*CQ*(V-1)*V
     ?   *(2*V**2-2*V+1)*V2*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**
     @   4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)-4*CQ*(2*V**2-2*V+1)**2*V1*W*
     1   (V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*
     2   W**2-4*V*W+1))/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)/3.D0
 
 
      LM = -LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**
     1   6-2*V**5*(2*V+1)*W**5+V**2*(8*V**4-8*V**3+12*V**2-4*V+1)*W**4-2
     2   *V*(4*V**5-2*V**4+6*V**2-3*V+1)*W**3+(8*V**6-8*V**4+16*V**3-2*V
     3   **2+1)*W**2-2*(4*V**2-2*V+1)*(2*V**3-2*V**2+V+1)*W+4*(V**2-V+1)
     4   *(2*V**2-2*V+1))+2*(V-1)*V*(2*V**2-2*V+1)*V2*(W**2-2*W+2)*(V**2
     5   *W**2-2*V*W+1))/((V-1)*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
 
      LMP = -LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2*(
     1   V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V1*(V**2*W**
     2   2-2*(V-1)*V*W+2*(V-1)**2)-2*(V-1)*V2*(V*W-V+1))/((V-1)*V**2*W**
     3   2*(V*W-1)**2*(V*W-V+1)**2)
 
      STRUV8=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV9(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V
     1   -2)*W**7+V**4*(52*V**4-50*V**3+17*V**2-4*V+1)*W**6-2*V**4*(42*V
     2   **4-59*V**3+30*V**2-7*V+2)*W**5+2*V**2*(44*V**6-64*V**5+15*V**4
     3   +25*V**3-16*V**2+5*V-1)*W**4-2*V**2*(28*V**6-22*V**5-60*V**4+11
     4   4*V**3-71*V**2+20*V-3)*W**3+(16*V**8+40*V**7-204*V**6+280*V**5-
     5   150*V**4+12*V**3+15*V**2-6*V+1)*W**2-2*(V-1)*(16*V**6-28*V**5+6
     6   *V**4+30*V**3-29*V**2+10*V-1)*W+4*(V-1)**2*(2*V**2-3*V+2)*(2*V*
     7   *2-2*V+1))+V1*(-4*V**8*W**8+4*V**7*(4*V-1)*W**7-V**4*(36*V**4-4
     8   0*V**3+23*V**2-4*V+1)*W**6+2*V**4*(28*V**4-59*V**3+49*V**2-12*V
     9   +2)*W**5-2*V**2*(32*V**6-79*V**5+54*V**4+17*V**3-20*V**2+5*V-1)
     :   *W**4+2*V**2*(24*V**6-46*V**5-23*V**4+121*V**3-91*V**2+24*V-3)*
     ;   W**3-(16*V**8+16*V**7-196*V**6+354*V**5-232*V**4+38*V**3+13*V**
     <   2-6*V+1)*W**2+2*(V-1)*(16*V**6-40*V**5+26*V**4+23*V**3-31*V**2+
     =   11*V-1)*W-4*(V-1)**2*(2*V**2-4*V+3)*(2*V**2-2*V+1))-2*(V-1)*V2*
     >   (V*W-1)*(V*W-V+1)*(V**3*(4*V**2+V+1)*W**4-V**3*(12*V**2+V-5)*W*
     ?   *3+V*(16*V**4+4*V**3-30*V**2+17*V-1)*W**2-(V-1)*V*(8*V**3+16*V*
     @   *2-28*V+13)*W+4*(V-1)**2*(2*V**2-2*V+1))-2*CQ*(V-1)*V2*(V*W-2*V
     1   +1)*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**3*W**3-V*(6*V**2-6*V+1)
     2   *W**2+2*V*(4*V**2-6*V+3)*W-2*(V-1)*(2*V**2-2*V+1)))/((V-1)*V**2
     3   *W**2*(V*W-1)**4*(V*W-V+1)**2)
 
 
      LV1 = -2*LOG(1-V)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)
     1   **2)*(CQ*V1*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-V**3*W
     2   **3+V**2*W**2-V*W+1)+V*V1*(V**2*W**2-2*V*W+1)*(2*V**3*W**4-2*V*
     3   *2*(3*V-2)*W**3+V*(8*V**2-11*V+5)*W**2-(2*V-1)*(2*V**2-3*V+3)*W
     4   +2*(2*V**2-2*V+1))-(V-1)*V2*(V*W-1)*(V**2*W**2+2*V*W+1)*(2*V**2
     5   *W**2-V*(4*V-1)*W+2*(2*V**2-2*V+1)))/((V-1)*V**2*W**2*(V*W-1)**
     6   4*(V*W-V+1)**2)
 
 
      LV = LOG(V)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V-2)
     1   *W**7+3*V**6*(2*V-1)*(8*V-3)*W**6-2*V**5*(34*V**3-43*V**2+14*V+
     2   1)*W**5+V**4*(60*V**4-80*V**3-2*V**2+46*V-15)*W**4-2*V**3*(8*V-
     3   5)*(2*V**4-6*V**2+6*V-1)*W**3+2*V**2*(4*V**6+12*V**5-64*V**4+88
     4   *V**3-42*V**2+3)*W**2-2*(V-1)*V*(8*V**5-16*V**4+2*V**3+22*V**2-
     5   19*V+5)*W+2*(V-1)**2*(2*V**2-4*V+3)*(2*V**2-2*V+1))+V1*(-4*V**8
     6   *W**8+4*V**7*(5*V-2)*W**7-V**4*(52*V**4-62*V**3+29*V**2-4*V+1)*
     7   W**6+2*V**4*(42*V**4-83*V**3+58*V**2-11*V+2)*W**5-2*V**2*(44*V*
     8   *6-102*V**5+57*V**4+29*V**3-24*V**2+5*V-1)*W**4+2*V**2*(28*V**6
     9   -50*V**5-48*V**4+166*V**3-115*V**2+28*V-3)*W**3-(16*V**8+24*V**
     :   7-252*V**6+452*V**5-298*V**4+56*V**3+11*V**2-6*V+1)*W**2+2*(V-1
     ;   )*(16*V**6-44*V**5+26*V**4+34*V**3-41*V**2+14*V-1)*W-4*(V-1)**2
     <   *(2*V**2-5*V+4)*(2*V**2-2*V+1))+2*(V-1)*V2*(V*W-1)*(V*W-V+1)*(2
     =   *V**5*W**5-V**3*(10*V**2+1)*W**4+2*V**3*(10*V**2-V-5)*W**3-V*(2
     >   0*V**4-42*V**2+24*V-1)*W**2+2*V*(4*V**4+4*V**3-28*V**2+28*V-9)*
     ?   W-2*(V-1)*(2*V-3)*(2*V**2-2*V+1))-2*CQ*(V-1)*V2*(V*W-V+1)*(V**2
     @   *W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W**2-
     1   2*V*(2*V-1)*W+2*V**2-2*V+1))/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V
     2   +1)**2)
 
 
      LVW = 4*V1*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-V**3*(2*V-1
     1   )*W**3+V**2*(2*V**2-5*V+4)*W**2+V*(4*V**2-7*V+2)*W+2*V**2-2*V+1
     2   )*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*LOG(1-V*W)/((V-1)
     3   *V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
 
 
      LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*V1*(2*V**
     1   2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3
     2   +8*(V-1)**2*V**2*W**2-8*(V-1)**3*V*W+4*(V-1)**4)-(V-1)*V*V1*W*(
     3   V**3*W**3-V**2*(4*V-3)*W**2+(V-1)*V*(3*V-2)*W+(V-1)**2)-2*CQ*(V
     4   -1)*V2*(V*W-V+1)*(V**2*W**2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W**
     5   2-2*V*(2*V-1)*W+2*V**2-2*V+1)+(V-1)*V*V2*W*(V*W-V+1)*(V**2*W**2
     6   -V*(5*V-4)*W+(V-1)*(4*V-1)))*LOG(V*W-V+1)/((V-1)*V**2*W**2*(V*W
     7   -1)**4*(V*W-V+1)**2)
 
 
      LW = (V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W
     1   -1)*(2*V1*(V*W-1)*(2*(V-1)*V**3*W**3-V**2*(8*V**2-14*V+7)*W**2+
     2   2*V*(6*V**3-14*V**2+12*V-3)*W-2*(2*V**2-3*V+2)*(2*V**2-2*V+1))+
     3   4*(V-1)*V2*(V*W-2*V+1)*(V**3*W**3-V**2*(2*V+1)*W**2+V*(2*V**2-1
     4   )*W-2*V**2+2*V-1)-2*CQ*(V-1)*V*(2*V**2-2*V+1)*V2*(V*W-1)*(W**2-
     5   2*W+2)-CQ*(2*V**2-2*V+1)**2*V1*(V*W-1)*(W**2-2*W+2))*LOG(W)/((V
     6   -1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
 
 
      CVC = (-3*V1*(V**2*W**2-2*V*W+1)*(8*(V-1)*V**7*W**7-V**4*(32*V**4-
     1   43*V**3+10*V**2-3*V+2)*W**6+V**4*(64*V**4-101*V**3+37*V**2-11*V
     2   +7)*W**5-2*V**2*(40*V**6-66*V**5+25*V**4-3*V**3+5*V-2)*W**4+2*(
     3   V-1)*V**2*(28*V**5-V**4-43*V**3+33*V**2-14*V+5)*W**3-(V-1)*(16*
     4   V**7+56*V**6-156*V**5+124*V**4-26*V**3-9*V**2+5*V-2)*W**2+(V-1)
     5   **2*(32*V**5-24*V**4-26*V**3+52*V**2-21*V+3)*W-8*(V-1)**3*V*(2*
     6   V**2-3*V+2))+3*AL*V1*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W
     7   +(V-1)**2)*(2*V**6*W**6-2*V**5*(2*V+1)*W**5+V**2*(8*V**3-4*V**2
     8   +4*V-1)*W**4+2*V*(2*V**2-V+1)*(2*V**3-2*V**2-2*V+1)*W**3-(8*V**
     9   6-16*V**4+16*V**3-10*V**2+1)*W**2+2*(V-1)*(8*V**4-4*V**3+2*V**2
     :   +2*V-1)*W-4*(V-1)*V*(2*V**2-2*V+1))-12*(V-1)*V*V2*(V*W-V+1)*(V*
     ;   *2*W**2-2*V*W+1)*(V**5*W**6-V**3*(4*V**2-V+1)*W**5+V**2*(9*V**3
     <   -5*V**2+1)*W**4-V*(10*V**4-2*V**3-14*V**2+12*V-1)*W**3+(4*V**5+
     =   10*V**4-26*V**3+18*V**2-V-1)*W**2-(8*V**4-10*V**3-2*V**2+8*V-3)
     >   *W+(V-1)*(4*V**2-6*V+3))-6*AL*(V-1)*V*(2*V**2-2*V+1)*V2*(W**2-2
     ?   *W+2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6
     @   *V**2*W**2-4*V*W+1)-8*CQ*(V-1)*V*(2*V**2-2*V+1)*V2*W*(V**2*W**2
     1   -2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W
     2   +1)-4*CQ*(2*V**2-2*V+1)**2*V1*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2
     3   )*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1))/((V-1)*V**2*W**2
     4   *(V*W-1)**4*(V*W-V+1)**2)/3.D0
 
 
      LM = -LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**
     1   6-2*V**5*(2*V+1)*W**5+V**2*(8*V**4-8*V**3+12*V**2-4*V+1)*W**4-2
     2   *V*(4*V**5-2*V**4+6*V**2-3*V+1)*W**3+(8*V**6-8*V**4+16*V**3-2*V
     3   **2+1)*W**2-2*(4*V**2-2*V+1)*(2*V**3-2*V**2+V+1)*W+4*(V**2-V+1)
     4   *(2*V**2-2*V+1))+2*(V-1)*V*(2*V**2-2*V+1)*V2*(W**2-2*W+2)*(V**2
     5   *W**2-2*V*W+1))/((V-1)*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
 
      LMP = -LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2*(
     1   V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V1*(V**2*W**
     2   2-2*(V-1)*V*W+2*(V-1)**2)-2*(V-1)*V2*(V*W-V+1))/((V-1)*V**2*W**
     3   2*(V*W-1)**2*(V*W-V+1)**2)
 
      STRUV9=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV10(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V
     1   **7*(5*V-2)*W**7+V**4*(52*V**4-50*V**3+17*V**2-4*V+1)*W*
     2   *6-2*V**4*(42*V**4-59*V**3+30*V**2-7*V+2)*W**5+2*V**2*(4
     3   4*V**6-64*V**5+15*V**4+25*V**3-16*V**2+5*V-1)*W**4-2*V**
     4   2*(28*V**6-22*V**5-60*V**4+114*V**3-71*V**2+20*V-3)*W**3
     5   +(16*V**8+40*V**7-204*V**6+280*V**5-150*V**4+12*V**3+15*
     6   V**2-6*V+1)*W**2-2*(V-1)*(16*V**6-28*V**5+6*V**4+30*V**3
     7   -29*V**2+10*V-1)*W+4*(V-1)**2*(2*V**2-3*V+2)*(2*V**2-2*V
     8   +1))+V1*(-4*V**8*W**8+4*V**7*(4*V-1)*W**7-V**4*(36*V**4-
     9   40*V**3+23*V**2-4*V+1)*W**6+2*V**4*(28*V**4-59*V**3+49*V
     :   **2-12*V+2)*W**5-2*V**2*(32*V**6-79*V**5+54*V**4+17*V**3
     ;   -20*V**2+5*V-1)*W**4+2*V**2*(24*V**6-46*V**5-23*V**4+121
     <   *V**3-91*V**2+24*V-3)*W**3-(16*V**8+16*V**7-196*V**6+354
     =   *V**5-232*V**4+38*V**3+13*V**2-6*V+1)*W**2+2*(V-1)*(16*V
     >   **6-40*V**5+26*V**4+23*V**3-31*V**2+11*V-1)*W-4*(V-1)**2
     ?   *(2*V**2-4*V+3)*(2*V**2-2*V+1))-2*(V-1)*V2*(V*W-1)*(V*W-
     @   V+1)*(V**3*(4*V**2+V+1)*W**4-V**3*(12*V**2+V-5)*W**3+V*(
     1   16*V**4+4*V**3-30*V**2+17*V-1)*W**2-(V-1)*V*(8*V**3+16*V
     2   **2-28*V+13)*W+4*(V-1)**2*(2*V**2-2*V+1))-2*CQ*(V-1)*V2*
     3   (V*W-2*V+1)*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**3*W**3-V
     4   *(6*V**2-6*V+1)*W**2+2*V*(4*V**2-6*V+3)*W-2*(V-1)*(2*V**
     5   2-2*V+1))+4*CQ*V*V4*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(
     6   V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)-4*(
     7   V-1)*V4*W*(V*W-2*V+1)*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**
     8   2*W**2-2*V**2*W+2*V**2-2*V+1)-4*(V-1)*V*V3*W*(V*W-1)*(V*
     9   W-V+1)*(V*W**2-(V-1)*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V
     :   +1))/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
 
      LV1 = -2*LOG(1-V)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V
     1   *W+(V-1)**2)*(CQ*V1*(V*W-V+1)*(V**2*W**2-2*V**2*W+2*V**2
     2   -2*V+1)*(V**4*W**4-V**3*W**3+V**2*W**2-V*W+1)+V*V1*(V*W-
     3   V+1)*(V**2*W**2-2*V*W+1)*(2*V**3*W**4-2*V**2*(3*V-2)*W**
     4   3+V*(8*V**2-11*V+5)*W**2-(2*V-1)*(2*V**2-3*V+3)*W+2*(2*V
     5   **2-2*V+1))-(V-1)*V2*(V*W-1)*(V*W-V+1)*(V**2*W**2+2*V*W+
     6   1)*(2*V**2*W**2-V*(4*V-1)*W+2*(2*V**2-2*V+1))+2*CQ*V*V4*
     7   W*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V*
     8   *2-2*V+1)+2*V*V4*W*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**2*W
     9   **2-2*V**2*W+2*V**2-2*V+1)-2*(V-1)*V3*W*(V*W-1)*(V*(V+1)
     :   *W-V+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1))/((V-1)*V**2*W
     ;   **2*(V*W-1)**4*(V*W-V+1)**3)
 
      LV = LOG(V)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7
     1   *(5*V-2)*W**7+3*V**6*(2*V-1)*(8*V-3)*W**6-2*V**5*(34*V**
     2   3-43*V**2+14*V+1)*W**5+V**4*(60*V**4-80*V**3-2*V**2+46*V
     3   -15)*W**4-2*V**3*(8*V-5)*(2*V**4-6*V**2+6*V-1)*W**3+2*V*
     4   *2*(4*V**6+12*V**5-64*V**4+88*V**3-42*V**2+3)*W**2-2*(V-
     5   1)*V*(8*V**5-16*V**4+2*V**3+22*V**2-19*V+5)*W+2*(V-1)**2
     6   *(2*V**2-4*V+3)*(2*V**2-2*V+1))+V1*(-4*V**8*W**8+4*V**7*
     7   (5*V-2)*W**7-V**4*(52*V**4-62*V**3+29*V**2-4*V+1)*W**6+2
     8   *V**4*(42*V**4-83*V**3+58*V**2-11*V+2)*W**5-2*V**2*(44*V
     9   **6-102*V**5+57*V**4+29*V**3-24*V**2+5*V-1)*W**4+2*V**2*
     :   (28*V**6-50*V**5-48*V**4+166*V**3-115*V**2+28*V-3)*W**3-
     ;   (16*V**8+24*V**7-252*V**6+452*V**5-298*V**4+56*V**3+11*V
     <   **2-6*V+1)*W**2+2*(V-1)*(16*V**6-44*V**5+26*V**4+34*V**3
     =   -41*V**2+14*V-1)*W-4*(V-1)**2*(2*V**2-5*V+4)*(2*V**2-2*V
     >   +1))+2*(V-1)*V2*(V*W-1)*(V*W-V+1)*(2*V**5*W**5-V**3*(10*
     ?   V**2+1)*W**4+2*V**3*(10*V**2-V-5)*W**3-V*(20*V**4-42*V**
     @   2+24*V-1)*W**2+2*V*(4*V**4+4*V**3-28*V**2+28*V-9)*W-2*(V
     1   -1)*(2*V-3)*(2*V**2-2*V+1))-2*CQ*(V-1)*V2*(V*W-V+1)*(V**
     2   2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2*(V-1)**2)*(2*V*
     3   *2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)+4*CQ*V*V4*W*(V**2*W*
     4   *2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-
     5   2*V**2*W+2*V**2-2*V+1)-4*V*V4*W*(V*W-2*(V-1))*(V*W-V+1)*
     6   (V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)-4*
     7   (V-1)*V*V3*W*(V*W-1)*(V*W+1)*(V*W-V+1)*(V**2*W**2-2*V**2
     8   *W+2*V**2-2*V+1))/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)*
     9   *2)
 
      LVW = 4*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**
     1   3)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V1*(V**4
     2   *W**4-V**3*(2*V-1)*W**3+V**2*(2*V**2-5*V+4)*W**2+V*(4*V*
     3   *2-7*V+2)*W+2*V**2-2*V+1)+2*V*V4*W*(V**2*W**2-2*V**2*W+2
     4   *V**2-2*V+1))*LOG(1-V*W)/((V-1)*V**2*W**2*(V*W-1)**4*(V
     5   *W-V+1)**3)
 
      LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*V1
     1   *(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4-4*(
     2   V-1)*V**3*W**3+8*(V-1)**2*V**2*W**2-8*(V-1)**3*V*W+4*(V-
     3   1)**4)-(V-1)*V*V1*W*(V**3*W**3-V**2*(4*V-3)*W**2+(V-1)*V
     4   *(3*V-2)*W+(V-1)**2)-2*CQ*(V-1)*V2*(V*W-V+1)*(V**2*W**2-
     5   2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**
     6   2-2*V+1)+(V-1)*V*V2*W*(V*W-V+1)*(V**2*W**2-V*(5*V-4)*W+(
     7   V-1)*(4*V-1)))*LOG(V*W-V+1)/((V-1)*V**2*W**2*(V*W-1)**4
     8   *(V*W-V+1)**2)
 
      LW = -(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W*
     1   *2+3*V*W-1)*(-2*V1*(V*W-1)*(V*W-V+1)*(2*(V-1)*V**3*W**3-
     2   V**2*(8*V**2-14*V+7)*W**2+2*V*(6*V**3-14*V**2+12*V-3)*W-
     3   2*(2*V**2-3*V+2)*(2*V**2-2*V+1))-4*(V-1)*V2*(V*W-2*V+1)*
     4   (V*W-V+1)*(V**3*W**3-V**2*(2*V+1)*W**2+V*(2*V**2-1)*W-2*
     5   V**2+2*V-1)+4*V4*W*(V*W-1)*(V*W-V+1)*(V**2*W**2-2*V**2*W
     6   +2*V**2-2*V+1)-4*(V-1)*V*V3*W**2*(V*W-2*V+1)*(V**2*W**2-
     7   2*V**2*W+2*V**2-2*V+1)+2*CQ*(V-1)*V*(2*V**2-2*V+1)*V2*(V
     8   *W-1)*(V*W-V+1)*(W**2-2*W+2)+CQ*(2*V**2-2*V+1)**2*V1*(V*
     9   W-1)*(V*W-V+1)*(W**2-2*W+2))*LOG(W)/((V-1)*V**2*W**2*(V
     :   *W-1)**4*(V*W-V+1)**3)
 
      CVC = (AL*V1*(6*V**11*W**11-30*V**11*W**10+(54*V**11+30*V**
     1   10-30*V**9+12*V**8-3*V**7)*W**9+(-18*V**11-144*V**10+120
     2   *V**9-36*V**8+3*V**7+3*V**6)*W**8+(-84*V**11+258*V**10-9
     3   6*V**9-108*V**8+129*V**7-60*V**6+9*V**5)*W**7+(144*V**11
     4   -132*V**10-306*V**9+540*V**8-375*V**7+111*V**6+9*V**5-9*
     5   V**4)*W**6+(-96*V**11-180*V**10+744*V**9-678*V**8+126*V*
     6   *7+270*V**6-264*V**5+90*V**4-9*V**3)*W**5+(24*V**11+264*
     7   V**10-480*V**9-120*V**8+834*V**7-900*V**6+474*V**5-78*V*
     8   *4-27*V**3+9*V**2)*W**4+(-96*V**10-96*V**9+816*V**8-1140
     9   *V**7+642*V**6+24*V**5-300*V**4+195*V**3-48*V**2+3*V)*W*
     :   *3+(144*V**9-336*V**8+60*V**7+540*V**6-774*V**5+528*V**4
     ;   -177*V**3+3*V**2+15*V-3)*W**2+(-96*V**8+384*V**7-612*V**
     <   6+480*V**5-150*V**4-60*V**3+84*V**2-36*V+6)*W+24*V**7-12
     =   0*V**6+252*V**5-288*V**4+192*V**3-72*V**2+12*V)+V2*((12*
     >   V**10-12*V**11)*W**10+(72*V**11-84*V**10+24*V**9-12*V**8
     ?   )*W**9+(-216*V**11+276*V**10-60*V**9-12*V**8+12*V**7)*W*
     @   *8+(384*V**11-420*V**10-264*V**9+516*V**8-252*V**7+36*V*
     1   *6)*W**7+(-396*V**11+72*V**10+1416*V**9-1752*V**8+684*V*
     2   *7+12*V**6-36*V**5)*W**6+(216*V**11+552*V**10-2256*V**9+
     3   1764*V**8+516*V**7-1224*V**6+468*V**5-36*V**4)*W**5+(-48
     4   *V**11-600*V**10+1188*V**9+840*V**8-3576*V**7+3132*V**6-
     5   936*V**5-36*V**4+36*V**3)*W**4+(192*V**10+240*V**9-2256*
     6   V**8+3384*V**7-1428*V**6-864*V**5+996*V**4-276*V**3+12*V
     7   **2)*W**3+(-288*V**9+720*V**8+108*V**7-2160*V**6+2928*V*
     8   *5-1680*V**4+360*V**3+24*V**2-12*V)*W**2+(192*V**8-840*V
     9   **7+1416*V**6-1032*V**5+84*V**4+348*V**3-204*V**2+36*V)*
     :   W-48*V**7+264*V**6-612*V**5+768*V**4-552*V**3+216*V**2-3
     ;   6*V)+V1*((24*V**10-24*V**11)*W**10+(120*V**11-129*V**10+
     <   6*V**9-9*V**8+6*V**7)*W**9+(-288*V**11+288*V**10+60*V**9
     =   -12*V**8-18*V**7-6*V**6)*W**8+(432*V**11-315*V**10-372*V
     >   **9+201*V**8-36*V**7+72*V**6-18*V**5)*W**7+(-408*V**11-5
     ?   4*V**10+1200*V**9-888*V**8+318*V**7-180*V**6+18*V**5+18*
     @   V**4)*W**6+(216*V**11+594*V**10-1902*V**9+1245*V**8-48*V
     1   **7-243*V**6+240*V**5-126*V**4+18*V**3)*W**5+(-48*V**11-
     2   600*V**10+1056*V**9+600*V**8-2340*V**7+2040*V**6-912*V**
     3   5+204*V**4+18*V**3-18*V**2)*W**4+(192*V**10+240*V**9-205
     4   2*V**8+3048*V**7-1653*V**6-132*V**5+591*V**4-300*V**3+72
     5   *V**2-6*V)*W**3+(-288*V**9+720*V**8-24*V**7-1794*V**6+26
     6   28*V**5-1704*V**4+510*V**3-36*V**2-18*V+6)*W**2+(192*V**
     7   8-840*V**7+1458*V**6-1146*V**5+147*V**4+432*V**3-324*V**
     8   2+90*V-9)*W-48*V**7+264*V**6-624*V**5+816*V**4-624*V**3+
     9   264*V**2-48*V)+AL*V4*(12*V**10*W**10+(-60*V**10-12*V**9)
     :   *W**9+(132*V**10+72*V**9-24*V**8)*W**8+(-156*V**10-204*V
     ;   **9+96*V**8+24*V**7)*W**7+(96*V**10+336*V**9-180*V**8-72
     <   *V**7)*W**6+(-24*V**10-288*V**9+60*V**8+252*V**7-72*V**6
     =   )*W**5+(96*V**9+192*V**8-480*V**7+252*V**6-72*V**5+24*V*
     >   *4)*W**4+(-144*V**8+192*V**7+60*V**6-180*V**5+96*V**4-24
     ?   *V**3)*W**3+(96*V**7-288*V**6+336*V**5-204*V**4+72*V**3-
     @   12*V**2)*W**2+(-24*V**6+96*V**5-156*V**4+132*V**3-60*V**
     1   2+12*V)*W)+V4*(-12*V**10*W**10+(36*V**10+36*V**9)*W**9+(
     2   -36*V**10-144*V**9)*W**8+(12*V**10+180*V**9+144*V**8-96*
     3   V**7)*W**7+(-72*V**9-324*V**8+144*V**7+72*V**6)*W**6+(18
     4   0*V**8+180*V**7-360*V**6+72*V**5)*W**5+(-240*V**7+180*V*
     5   *6+144*V**5-96*V**4)*W**4+(180*V**6-324*V**5+144*V**4)*W
     6   **3+(-72*V**5+180*V**4-144*V**3+36*V**2)*W**2+(12*V**4-3
     7   6*V**3+36*V**2-12*V)*W)+AL*V2*((-12*V**11+24*V**10-18*V*
     8   *9+6*V**8)*W**9+(60*V**11-108*V**10+66*V**9-12*V**8-6*V*
     9   *7)*W**8+(-132*V**11+168*V**10+30*V**9-150*V**8+102*V**7
     :   -18*V**6)*W**7+(156*V**11-36*V**10-426*V**9+516*V**8-228
     ;   *V**7+18*V**5)*W**6+(-96*V**11-216*V**10+744*V**9-492*V*
     <   *8-156*V**7+360*V**6-162*V**5+18*V**4)*W**5+(24*V**11+26
     =   4*V**10-444*V**9-264*V**8+960*V**7-756*V**6+198*V**5+36*
     >   V**4-18*V**3)*W**4+(-96*V**10-96*V**9+816*V**8-1044*V**7
     ?   +360*V**6+270*V**5-294*V**4+90*V**3-6*V**2)*W**3+(144*V*
     @   *9-336*V**8+36*V**7+588*V**6-738*V**5+372*V**4-48*V**3-2
     1   4*V**2+6*V)*W**2+(-96*V**8+384*V**7-600*V**6+432*V**5-84
     2   *V**4-84*V**3+60*V**2-12*V)*W+24*V**7-120*V**6+252*V**5-
     3   288*V**4+192*V**3-72*V**2+12*V)+CQ*V2*((-16*V**11+32*V**
     4   10-24*V**9+8*V**8)*W**8+(48*V**11-80*V**10+40*V**9-8*V**
     5   7)*W**7+(-48*V**11+168*V**9-216*V**8+120*V**7-24*V**6)*W
     6   **6+(16*V**11+112*V**10-312*V**9+256*V**8-48*V**7-48*V**
     7   6+24*V**5)*W**5+(-64*V**10+32*V**9+288*V**8-544*V**7+432
     8   *V**6-168*V**5+24*V**4)*W**4+(96*V**9-288*V**8+288*V**7-
     9   48*V**6-120*V**5+96*V**4-24*V**3)*W**3+(-64*V**8+272*V**
     :   7-480*V**6+456*V**5-248*V**4+72*V**3-8*V**2)*W**2+(16*V*
     ;   *7-80*V**6+168*V**5-192*V**4+128*V**3-48*V**2+8*V)*W)+CQ
     <   *V1*((-16*V**11+32*V**10-32*V**9+16*V**8-4*V**7)*W**8+(4
     =   8*V**11-80*V**10+64*V**9-16*V**8-4*V**7+4*V**6)*W**7+(-4
     >   8*V**11+144*V**9-240*V**8+180*V**7-72*V**6+12*V**5)*W**6
     ?   +(16*V**11+112*V**10-304*V**9+320*V**8-140*V**7-12*V**6+
     @   36*V**5-12*V**4)*W**5+(-64*V**10+32*V**9+256*V**8-560*V*
     1   *7+560*V**6-312*V**5+96*V**4-12*V**3)*W**4+(96*V**9-288*
     2   V**8+336*V**7-144*V**6-72*V**5+120*V**4-60*V**3+12*V**2)
     3   *W**3+(-64*V**8+272*V**7-512*V**6+560*V**5-384*V**4+164*
     4   V**3-40*V**2+4*V)*W**2+(16*V**7-80*V**6+176*V**5-224*V**
     5   4+180*V**3-92*V**2+28*V-4)*W))/((V-1)*V**2*W**2*(V*W-1)*
     6   *4*(V*W-V+1)**3)/3.D0
 
      LM = -LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*
     1   V**6*W**6-2*V**5*(2*V+1)*W**5+V**2*(8*V**4-8*V**3+12*V**
     2   2-4*V+1)*W**4-2*V*(4*V**5-2*V**4+6*V**2-3*V+1)*W**3+(8*V
     3   **6-8*V**4+16*V**3-2*V**2+1)*W**2-2*(4*V**2-2*V+1)*(2*V*
     4   *3-2*V**2+V+1)*W+4*(V**2-V+1)*(2*V**2-2*V+1))+4*V*V4*W*(
     5   V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)+2*(
     6   V-1)*V*(2*V**2-2*V+1)*V2*(W**2-2*W+2)*(V**2*W**2-2*V*W+1
     7   ))/((V-1)*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
      LMP = -LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)
     1   *V*W+2*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1
     2   )*(V1*(V**2*W**2-2*(V-1)*V*W+2*(V-1)**2)-2*(V-1)*V2*(V*W
     3   -V+1))/((V-1)*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
 
      STRUV10=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV11(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(CQ*V1*(2*V**10*W**9-4*V**8*(V**2+2*V-2)*W**8+3*V**
     1   8*(2*V**2+3*V-3)*W**7-V**6*(8*V**4+V**3-15*V**2+28*V-14)*W**6+3
     2   *V**6*(2*V**4+V**3-10*V**2+18*V-9)*W**5-V**4*(4*V**6-3*V**5-7*V
     3   **4+26*V**3-28*V**2+18*V-6)*W**4+V**4*(2*V**6-17*V**4+45*V**3-5
     4   0*V**2+33*V-11)*W**3-(V-1)*V**2*(6*V**6-18*V**5+17*V**4+5*V**2-
     5   6*V+2)*W**2+(V-1)**2*V**2*(6*V**4-20*V**3+25*V**2-10*V+5)*W-2*(
     6   V-1)**3*(V**4-3*V**3+4*V**2-2*V+1))+V1*(-2*V**10*W**9+4*V**8*(V
     7   **2+2*V-2)*W**8-V**8*(6*V**2+11*V-11)*W**7+V**6*(16*V**4-21*V**
     8   3+19*V**2+4*V-2)*W**6-V**6*(22*V**4-39*V**3+24*V**2+30*V-15)*W*
     9   *5+3*V**4*(4*V**6+V**5-27*V**4+62*V**3-56*V**2+30*V-10)*W**4-V*
     :   *4*(2*V**6+24*V**5-93*V**4+143*V**3-84*V**2+15*V-5)*W**3+(V-1)*
     ;   V**2*(6*V**6+6*V**5-65*V**4+128*V**3-89*V**2+30*V-10)*W**2-(V-1
     <   )**2*V**2*(6*V**4-12*V**3+V**2+22*V-11)*W+2*(V-1)**3*(V**4-3*V*
     =   *3+4*V**2-2*V+1))+2*CQ*V**2*V2*W*(2*V**6*(V**2-2*V+2)*W**7-4*V*
     >   *6*(V**2-V+1)*W**6+2*V**4*(2*V**4-V**3+4*V**2-6*V+3)*W**5-V**4*
     ?   (4*V**4-6*V**3+13*V**2-14*V+7)*W**4+2*V**2*(V**6+V**5-8*V**4+12
     @   *V**3-V**2-6*V+2)*W**3-(V-1)*V**2*(6*V**4-11*V**3-V**2+24*V-12)
     1   *W**2+2*(V-1)**2*(3*V**4-7*V**3+8*V**2-2*V+1)*W-(V-1)**3*(2*V**
     2   2-3*V+3))-2*(V-1)*V**2*V2*W*(3*V**6*W**6-V**4*(7*V**2+4*V-4)*W*
     3   *5+V**4*(V**2+29*V-29)*W**4+3*V**2*(V**4-6*V**3-6*V**2+24*V-12)
     4   *W**3-(V-1)*V**2*(9*V**2-43*V+43)*W**2+(V-1)**2*(7*V**2-16*V+16
     5   )*W-3*(V-1)**3)+2*(V-1)*V**2*V4*(V**2*W**2-2*V*W+1)*(V**2*W**2-
     6   2*(V-1)*V*W+(V-1)**2)*(2*V**2*W**4-V**2*W**3-4*(V**2-2*V+2)*W**
     7   2+2*(V-1))+2*(V-1)*V**2*V3*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(
     8   V-1)*V*W+(V-1)**2)*(2*V**2*W**3-3*V**2*W**2+4*(V-1)*W-2*(V**2-2
     9   *V+2)))/((V-1)**2*V*W*(V*W-1)**3*(V*W-V+1)**3)
 
 
      LV1 = -LOG(1-V)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)*
     1   *3)*(2*V1*(V**7*(3*V-1)*W**7-4*V**6*(4*V-3)*W**6+V**5*(3*V**3-2
     2   *V**2+25*V-20)*W**5-V**4*(10*V**4-27*V**3+46*V**2-20*V+5)*W**4+
     3   V**3*(30*V**4-98*V**3+147*V**2-92*V+31)*W**3-V**2*(30*V**4-103*
     4   V**3+151*V**2-89*V+21)*W**2+V*(10*V**4-33*V**3+47*V**2-23*V+1)*
     5   W-(V-1)*(V**2+1))-2*(V-1)*V**2*V4*(V*W-1)*(6*V**3*W**5-4*V**2*W
     6   **4-V*(7*V**2+2*V+2)*W**3+(22*V**2-15*V+4)*W**2-(2*V**2+11*V-11
     7   )*W+2*(V-1))+CQ*V1*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(
     8   V**4*W**4+2*V**2*W**2+1)+2*V*V2*W*(V**2*W**2-2*V*W+1)*(2*V**3*(
     9   V**2-2*V+2)*W**3-(V-1)*V**2*(7*V-5)*W**2+V*(10*V**4-37*V**3+60*
     :   V**2-55*V+24)*W-(V-1)*(10*V**3-26*V**2+19*V-1))+2*CQ*V*V2*W*(V*
     ;   W-V+1)*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)-2*(V-1)*V
     <   **2*V3*W*(2*(V-1)*V*W**2+(3*V**2-4*V+4)*W+(V-1)*(2*V-5))*(V**2*
     =   W**2-2*V*W+1))/((V-1)**2*V*W*(V*W-1)**3*(V*W-V+1)**4)
 
 
      LV = LOG(V)*(CQ*V1*(2*V**10*W**9-4*V**8*(V**2+2*V-2)*W**8+3*V**8*(
     1   2*V**2+3*V-3)*W**7-V**6*(8*V**4+V**3-15*V**2+28*V-14)*W**6+3*V*
     2   *6*(2*V**4+V**3-10*V**2+18*V-9)*W**5-V**4*(4*V**6-3*V**5-7*V**4
     3   +26*V**3-28*V**2+18*V-6)*W**4+V**4*(2*V**6-17*V**4+45*V**3-50*V
     4   **2+33*V-11)*W**3-(V-1)*V**2*(6*V**6-18*V**5+17*V**4+5*V**2-6*V
     5   +2)*W**2+(V-1)**2*V**2*(6*V**4-20*V**3+25*V**2-10*V+5)*W-2*(V-1
     6   )**3*(V**4-3*V**3+4*V**2-2*V+1))+V1*(-2*V**10*W**9+4*V**8*(V**2
     7   +2*V-2)*W**8-V**6*(8*V**4+3*V**3-V**2-4*V+2)*W**7+V**6*(16*V**4
     8   -23*V**3+17*V**2+12*V-6)*W**6-3*V**4*(6*V**6-9*V**5-2*V**4+24*V
     9   **3-17*V**2+6*V-2)*W**5+V**4*(12*V**6-9*V**5-47*V**4+136*V**3-1
     :   28*V**2+72*V-24)*W**4-V**2*(4*V**8+12*V**7-75*V**6+137*V**5-90*
     ;   V**4+9*V**3+25*V**2-24*V+6)*W**3+(V-1)*V**2*(12*V**6-24*V**5-7*
     <   V**4+72*V**3-61*V**2+30*V-10)*W**2-(V-1)**2*(12*V**6-36*V**5+41
     =   *V**4-8*V**3-V**2+6*V-2)*W+4*(V-1)**3*(V**4-3*V**3+4*V**2-2*V+1
     >   ))-2*V**2*V2*W*(2*V**6*(V**2-2*V+2)*W**7-4*(V-2)**2*V**6*W**6-2
     ?   *V**4*(6*V**4-17*V**3+42*V**2-50*V+25)*W**5+V**4*(28*V**4-86*V*
     @   *3+175*V**2-178*V+89)*W**4-2*V**2*(7*V**6+V**5-54*V**4+164*V**3
     1   -227*V**2+174*V-58)*W**3+(V-1)*V**2*(42*V**4-117*V**3+173*V**2-
     2   112*V+56)*W**2-6*(V-1)**2*(7*V**4-21*V**3+30*V**2-18*V+9)*W+(V-
     3   1)**3*(14*V**2-37*V+37))+2*CQ*V**2*V2*W*(2*V**6*(V**2-2*V+2)*W*
     4   *7-4*V**6*(V**2-V+1)*W**6+2*V**4*(2*V**4-V**3+4*V**2-6*V+3)*W**
     5   5-V**4*(4*V**4-6*V**3+13*V**2-14*V+7)*W**4+2*V**2*(V**6+V**5-8*
     6   V**4+12*V**3-V**2-6*V+2)*W**3-(V-1)*V**2*(6*V**4-11*V**3-V**2+2
     7   4*V-12)*W**2+2*(V-1)**2*(3*V**4-7*V**3+8*V**2-2*V+1)*W-(V-1)**3
     8   *(2*V**2-3*V+3))+4*(V-1)*V**2*V3*W*(V**2*W**2-2*V*W+1)*(V**2*W*
     9   *2-2*(V-1)*V*W+(V-1)**2)*(4*V**2*W**3-2*V**2*W**2-2*(3*V**2-4*V
     :   +4)*W-V**2+8*V-8)-4*(V-1)*V**2*V4*(V**2*W**2-2*V*W+1)*(V**2*W**
     ;   2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**3-4*(V-1)*W**2+3*(V-1)*W-V+1))
     <   /((V-1)**2*V*W*(V*W-1)**3*(V*W-V+1)**3)
 
 
      LVW = 2*(V1*(2*V**4*W**3+2*(V-2)**2*V**2*W**2+V*(4*V**3-11*V**2+12
     1   *V+1)*W-V**2-1)+V*V2*W*(4*V*W+2*V**2+V+1)-(V-1)*V**2*V4*W*(4*W+
     2   3)-(V-1)*V**2*V3*W*(4*W+3))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-
     3   4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V
     4   -1)**3*V*W+(V-1)**4)*LOG(1-V*W)/((V-1)**2*V*W*(V*W-1)**4*(V*W-V
     5   +1)**4)
 
 
      LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(2*V*V1*W*(V*
     1   *5*W**5-V**4*(4*V-5)*W**4+V**3*(9*V**2-21*V+14)*W**3-(V-1)*V**2
     2   *(4*V**3-3*V**2-7*V+12)*W**2+(V-1)**2*V*(8*V**3-18*V**2+13*V+3)
     3   *W-(V-1)**3*(4*V**3-9*V**2+8*V-1))+CQ*(V-1)*V1*(2*V**2*W**2-2*V
     4   *(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)*
     5   *4)-(V-1)*V**2*V4*W*(4*W+3)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1
     6   )**2*V*W-(V-1)**3)-(V-1)*V**2*V3*W*(4*W+3)*(V**3*W**3-3*(V-1)*V
     7   **2*W**2+3*(V-1)**2*V*W-(V-1)**3)+(V-1)*V*V2*W*(V**3*(9*V-7)*W*
     8   *3-(V-1)*V**2*(8*V**2+3*V-3)*W**2+(V-1)**2*V*(16*V**2-13*V+3)*W
     9   -(V-1)**3*(8*V**2-7*V+3))+2*CQ*(V-1)**2*V*V2*W*(V**2*W**2+(V-1)
     :   **2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))*LOG(V*W-V+1)/((V
     ;   -1)**2*V*W*(V*W-1)**4*(V*W-V+1)**3)
 
 
      LW = -2*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W*
     1   *2+3*(V-1)**2*V*W-(V-1)**3)*((V-1)*V**2*V4*(2*V**2*W**4-3*V**2*
     2   W**3+4*V**2*W**2-2*(2*V**2-V+1)*W+4*(V-1))+V1*(V**2*(3*V**4-9*V
     3   **3+10*V**2-2*V+1)*W**3-2*V**2*(V**4-2*V**3+V**2+2*V-1)*W**2+(3
     4   *V**6-10*V**5+13*V**4-5*V**3+3*V-1)*W-3*(V-1)*(V**4-3*V**3+4*V*
     5   *2-2*V+1))+V**2*V2*W*(2*V**2*(V**2-2*V+2)*W**3+8*(V-1)*V**2*W**
     6   2-2*(7*V**4-24*V**3+43*V**2-38*V+19)*W+(V-1)*(14*V**2-31*V+31))
     7   -(V-1)*V**2*V3*W*(6*V**2*W**3-V**2*W**2-12*(V**2-V+1)*W+12*(V-1
     8   )))*LOG(W)/((V-1)**2*V*W*(V*W-1)**4*(V*W-V+1)**4)
 
 
      CVC = (3*(V-1)*V1*(V*W-1)*(V*W-V+1)*(6*V**10*W**9-8*V**8*(V**2+4*V
     1   -3)*W**8-V**6*(18*V**4-121*V**3+109*V**2-24*V+12)*W**7+V**6*(52
     2   *V**4-193*V**3+207*V**2-100*V+62)*W**6-V**4*(42*V**6-75*V**5-88
     3   *V**4+282*V**3-183*V**2+92*V-36)*W**5+V**4*(4*V**6+95*V**5-385*
     4   V**4+462*V**3-24*V**2-202*V+54)*W**4+V**2*(6*V**8-48*V**7+45*V*
     5   *6+293*V**5-876*V**4+965*V**3-495*V**2+152*V-36)*W**3-(V-1)*V**
     6   2*(18*V**6-102*V**5+229*V**4-256*V**3+149*V**2-38*V+18)*W**2+(V
     7   -1)**2*(18*V**6-76*V**5+131*V**4-114*V**3+75*V**2-28*V+12)*W-2*
     8   (V-1)**3*(3*V**4-9*V**3+10*V**2-2*V+1))+12*(V-1)*V2*W*(V*W-1)*(
     9   V*W-V+1)*(2*V**8*(V**2-2*V+2)*W**7-V**6*(4*V**4-3*V**3-V**2+8*V
     :   -4)*W**6+V**6*(4*V**4-V**3-5*V**2+12*V-6)*W**5-V**4*(4*V**6-5*V
     ;   **5-5*V**4+8*V**3+26*V**2-36*V+12)*W**4+V**4*(V**2+2*V-2)*(2*V*
     <   *4-V**3-17*V**2+36*V-18)*W**3-(V-1)*V**2*(6*V**6-10*V**5-7*V**4
     =   +22*V**3+19*V**2-36*V+12)*W**2+(V-1)**2*V**2*(6*V**4-13*V**3+27
     >   *V**2-28*V+14)*W-2*(V-1)**3*(V**2-V+1)*(V**2+2*V-2))-12*(V-1)**
     ?   2*V**2*V4*(V**2*W**2-2*V*W+1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V
     @   -1)**2*V*W-(V-1)**3)*(2*V**3*W**5-2*V**2*(2*V+3)*W**4+2*V*(V**2
     1   +9*V-3)*W**3+(V**3-14*V**2+6*V+2)*W**2-(2*V**2-9*V+8)*W+V-1)-6*
     2   AL*(V-1)*V1*(V*W-1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4
     3   +2*V**2*W**2+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W*
     4   *2-4*(V-1)**3*V*W+(V-1)**4)-8*CQ*(V-1)**2*V1*(V**2*W**2-2*V*W+1
     5   )*(V**4*W**4+2*V**2*W**2+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1
     6   )**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)-12*AL*(V-1)*V*V2*W*(V*W
     7   -1)*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4-4
     8   *(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)-
     9   16*CQ*(V-1)**2*V*V2*W*(V**2*W**2+1)*(V**2*W**2-2*V*W+1)*(V**4*W
     :   **4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)
     ;   **4)+12*(V-1)**3*V**2*V3*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-
     <   1)*V*W+(V-1)**2)*(2*V**2*W**3-6*V**2*W**2+(V**2+14*V-14)*W-4*(V
     =   -1)))/((V-1)**3*V*W*(V*W-1)**4*(V*W-V+1)**4)/6.D0
 
 
      LM = -LOG(S/M**2)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1
     1   )**3)*(V1*(2*V**7*W**6+2*V**5*(V**2-7*V+4)*W**5+V**3*(6*V**4-23
     2   *V**3+46*V**2-23*V+2)*W**4-2*V**2*(V**5+6*V**4-25*V**3+42*V**2-
     3   17*V+3)*W**3+2*V*(3*V**5-14*V**3+29*V**2-8*V+3)*W**2-2*(3*V**5-
     4   6*V**4+4*V**3+4*V**2+2*V+1)*W+2*V**4-6*V**3+9*V**2-4*V+3)-4*(V-
     5   1)*V**2*V4*(V*W-1)*(V**2*W**4+2*(V-1)*V*W**3-(V**2+6*V-2)*W**2+
     6   (2*V+3)*W-1)+2*V*V2*W*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V*
     7   *2+1))/((V-1)**2*V*W*(V*W-1)**3*(V*W-V+1)**3)
 
 
      LMP = -LOG(S/MP**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V
     1   1*(2*V**7*W**6-2*V**5*(2*V**2+V-4)*W**5+V**4*(6*V**3-2*V**2-19*
     2   V+19)*W**4-4*(V-1)*V**3*(3*V**3-4*V**2-2*V+6)*W**3+2*(V-1)**2*V
     3   **2*(7*V**3-14*V**2+8*V+5)*W**2-2*(V-1)**3*V*(4*V**3-10*V**2+9*
     4   V-1)*W+(V-1)**5*(2*V**2-2*V+1))-4*(V-1)*V**2*V4*W*(V*W-V+1)*(V*
     5   *2*W**3-(V-2)*V*W**2-2*(V-1)*(V+1)*W+2*(V-1)**2)+2*(V-1)**2*V*V
     6   2*W*(V**2*W**2+(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+
     7   1))/((V-1)**2*V*W*(V*W-1)**4*(V*W-V+1)**3)
 
      STRUV11=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV12(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(-N**3*VC*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V
     1   -1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(4*V**8*W**8-4*V*
     2   *7*(3*V+2)*W**7+V**4*(16*V**4+24*V**3+11*V**2+1)*W**6-2*V**4*(8
     3   *V**4+12*V**3+21*V**2-V+2)*W**5+2*V**2*(8*V**6+4*V**5+29*V**4+9
     4   *V**3-6*V**2+V-1)*W**4-2*V**2*(6*V**6+2*V**5+11*V**4+30*V**3-16
     5   *V**2-2*V-3)*W**3+(4*V**8+12*V**7-17*V**6+54*V**5-5*V**4-32*V**
     6   3+V**2-2*V+1)*W**2-2*(V-1)*(4*V**6-2*V**5+7*V**4+6*V**3+6*V**2-
     7   4*V-1)*W+4*(V-1)**2*(V**2+1)*(V**2-V+2))-4*V**8*W**8+4*V**7*(3*
     8   V+1)*W**7-V**4*(16*V**4+12*V**3+7*V**2+1)*W**6+2*V**4*(8*V**4+2
     9   *V**3+19*V**2-3*V+2)*W**5-2*V**2*(8*V**6-6*V**5+27*V**4+13*V**3
     :   -10*V**2+V-1)*W**4+2*V**2*(6*V**6-4*V**5+5*V**4+46*V**3-24*V**2
     ;   -2*V-3)*W**3-(4*V**8+8*V**7-29*V**6+74*V**5+7*V**4-56*V**3+9*V*
     <   *2-2*V+1)*W**2+2*(V-1)*(4*V**6-6*V**5+9*V**4+6*V**3+12*V**2-8*V
     =   -1)*W-4*(V-1)**2*(V**2+1)*(V**2-2*V+3))+N*VC*(V*W-1)*(V**4*W**4
     >   -4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4
     ?   )*(CQ*(4*V**8*W**8-4*V**7*(3*V+1)*W**7+V**4*(16*V**4+10*V**3+3*
     @   V**2-2*V+1)*W**6-2*V**4*(8*V**4+V**3+7*V**2-6*V+2)*W**5+2*V**2*
     1   (8*V**6-7*V**5+10*V**4-8*V**3-V**2+3*V-1)*W**4-2*(V-1)*V**2*(6*
     2   V**5+V**4-4*V**3+8*V**2-8*V+3)*W**3+(V-1)**2*(4*V**6+16*V**5-11
     3   *V**4+6*V**3-2*V**2-2*V+1)*W**2-2*(V-1)**3*(4*V**4+2*V**3+V**2-
     4   1)*W+4*(V-1)**4*(V**2+1))-4*V**8*W**8+4*V**7*(3*V+1)*W**7-V**4*
     5   (16*V**4+30*V**3-17*V**2-2*V+1)*W**6+2*V**4*(8*V**4+27*V**3-5*V
     6   **2-20*V+2)*W**5-2*V**2*(8*V**6+19*V**5+24*V**4-48*V**3-V**2+3*
     7   V-1)*W**4+2*(V-1)*V**2*(6*V**5+11*V**4+38*V**3+16*V**2-16*V+3)*
     8   W**3-(V-1)*(4*V**7+12*V**6+13*V**5+29*V**4-8*V**3-20*V**2+3*V-1
     9   )*W**2+2*(V-1)**2*(V**2+1)*(4*V**3-2*V**2+5*V+1)*W-4*(V-1)**4*(
     :   V**2+1))+2*N**4*(V-1)*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V*
     ;   *3*W**3-3*V**2*W**2+3*V*W-1)*(2*V**8*W**7-V**7*(6*V-5)*W**6+V**
     <   6*(10*V**2-16*V+11)*W**5-V**4*(30*V**4-67*V**3+67*V**2-31*V+4)*
     =   W**4+2*(V-1)*V**3*(31*V**4-73*V**3+79*V**2-43*V+8)*W**3-2*CQ*(V
     >   **2*W-(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*(
     ?   V-1)*V*(2*V**2-2*V+1)*W+(V-1)**2)-(V-1)**2*V**2*(58*V**4-151*V*
     @   *3+180*V**2-101*V+24)*W**2+(V-1)**3*V*(22*V**4-62*V**3+81*V**2-
     1   53*V+16)*W-(V-1)**4*(V**2-V+1)*(2*V**2-3*V+4))-2*N**2*(V-1)*VC*
     2   W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V
     3   **2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**7*W**6-2*V**5*(V**2-V-1
     4   )*W**5-CQ*(V*W-1)*(V**5*W**4-(V-1)*V**4*W**3-2*V**2*(3*V**2-5*V
     5   +1)*W**2-(V-1)*V*(V**3-8*V**2+9*V-4)*W+(V-1)**2*(V**3-2*V**2+5*
     6   V-2))+V**4*(4*V**3-13*V**2-V+3)*W**4-V**3*(12*V**4-42*V**3+29*V
     7   **2-9*V+4)*W**3+V**2*(10*V**5-35*V**4+23*V**3-V**2+3*V-4)*W**2-
     8   (V-1)*V*(2*V**5-27*V**3+32*V**2-19*V+6)*W+(V-1)**2*(2*V**4-10*V
     9   **3+11*V**2-12*V+5))-2*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**
     :   4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**
     ;   4)*(V**6*W**5-V**4*(4*V**2+V-1)*W**4+2*V**3*(V+2)*(3*V**2-3*V+1
     <   )*W**3+CQ*V*(V*W-1)*(V**2*W**2+(V-1)*V*W+(V-1)**2)*(V**2*W**2-2
     =   *V**2*W+V**2+1)-2*V**2*(2*V**4+4*V**3-5*V**2+2*V-2)*W**2+(V-1)*
     >   V**2*(V**3+5*V**2-V+5)*W-(V-1)**2*(V**3+V**2-V+1))-4*(CQ-1)*GTR
     ?   *N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*(
     @   V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*V**
     1   2*W**2+3*V*W-1)+4*(CQ-1)*GTR*N*(V-1)**2*V**2*VC*W**2*(V**2*W**2
     2   +(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W
     3   +V**2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1))/(N**2*(V-1)**2*V**2*W
     4   **2*(V*W-1)**3*(V*W-V+1)**6)
 
 
      LV1 = LOG(1-V)*(N**4*(V-1)*VC*W*(V*W-1)*(V**5*W**5-5*(V-1)*V**4*W*
     1   *4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(
     2   V-1)**5)*(12*V**7*W**6-36*V**6*W**5+4*V**5*(3*V**2-3*V+13)*W**4
     3   -2*V**3*(12*V**4-15*V**3+16*V**2+11*V+1)*W**3+V**2*(50*V**4-98*
     4   V**3+94*V**2-21*V+2)*W**2-V*(28*V**4-62*V**3+56*V**2-15*V-2)*W+
     5   (V-1)*(2*V**3-3*V**2+2))-2*N**2*(V-1)*VC*W*(V*W-1)*(6*V**6*W**5
     6   +V**5*(2*V-19)*W**4+V**3*(8*V**3-21*V**2+39*V-1)*W**3-V**2*(15*
     7   V**3-28*V**2+36*V+3)*W**2+V*(8*V**3-9*V**2+9*V+7)*W-V**3-V**2+2
     8   *V-3)*(V**6*W**6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-
     9   1)**3*V**3*W**3+15*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)+
     :   N**3*VC*(V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W
     ;   **3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(2*(V**2*W**
     <   2-2*V*W+1)*(2*V**5*W**5-2*V**4*(2*V+1)*W**4+V**3*(3*V**2+4*V+5)
     =   *W**3-V**2*(V**3+5*V**2+7*V-1)*W**2+2*V*(3*V**3+V**2+V-1)*W-2*(
     >   V-1)*(V+1)*(V**2+1))+CQ*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2
     ?   +1)*(V**4*W**4-4*V**3*W**3+8*V**2*W**2-8*V*W+4))+(V-1)*VC*W*(V*
     @   W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1
     1   )**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(2*V**6*W**5+2*V**4*(V*
     2   *2-3*V+1)*W**4-2*V**3*(2*V**3+3*V**2-7*V+3)*W**3+V**2*(2*V-1)*(
     3   2*V**3+4*V**2+3*V-6)*W**2-V*(8*V**4-8*V**2+7*V-6)*W+(V-1)*(V+2)
     4   *(4*V**2-5*V+4))-N*V*VC*W*(V*W-1)*(2*(V**2*W**2-2*V*W+1)*(2*V**
     5   4*W**4-4*V**4*W**3+V**2*(3*V**2+1)*W**2-(V-1)*V*(V**2+4*V-1)*W+
     6   4*(V-1)**2*(V+1))+CQ*V*W*(V*W-V+1)*(V**2*W**2-2*V*W+2)*(2*V**2*
     7   W**2-2*V*(V+1)*W+V**2+1))*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)
     8   **2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))/(
     9   N**2*(V-1)**2*V**2*W**2*(V*W-1)**3*(V*W-V+1)**6)
 
 
      LV = LOG(V)*(-N**3*VC*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)
     1   **2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(4*V**8*W**8-4*V**7*
     2   (3*V+2)*W**7+3*V**6*(V+1)*(5*V+3)*W**6-2*V**5*(6*V**3+12*V**2+1
     3   7*V-1)*W**5+V**4*(9*V**4+6*V**3+46*V**2+14*V-15)*W**4-2*V**3*(3
     4   *V+5)*(V**4-2*V**3+6*V**2-2*V-1)*W**3+2*V**2*(V**6+2*V**5-7*V**
     5   4+20*V**3+3*V**2-18*V+3)*W**2-2*(V-1)*V*(2*V**5-3*V**4+4*V**3+4
     6   *V**2+6*V-5)*W+2*(V-1)**2*(V**2+1)*(V**2-2*V+3))-4*V**8*W**8+4*
     7   V**7*(3*V+2)*W**7-V**4*(16*V**4+24*V**3+11*V**2+1)*W**6+2*V**4*
     8   (8*V**4+6*V**3+33*V**2-7*V+2)*W**5-2*V**2*(8*V**6-10*V**5+47*V*
     9   *4+15*V**3-16*V**2+V-1)*W**4+2*V**2*(6*V**6-10*V**5+11*V**4+64*
     :   V**3-36*V**2-4*V-3)*W**3-(4*V**8+4*V**7-41*V**6+106*V**5+7*V**4
     ;   -76*V**3+13*V**2-2*V+1)*W**2+2*(V-1)*(4*V**6-10*V**5+11*V**4+6*
     <   V**3+18*V**2-12*V-1)*W-4*(V-1)**2*(V**2+1)*(V**2-3*V+4))+N*VC*(
     =   V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1
     >   )**3*V*W+(V-1)**4)*(CQ*(4*V**8*W**8-4*V**7*(3*V+1)*W**7+V**6*(1
     ?   5*V**2+12*V+1)*W**6-2*V**5*(6*V**3+5*V**2+3*V-2)*W**5+V**4*(9*V
     @   **4-2*V**3+12*V**2-10*V-1)*W**4-2*(V-1)*V**3*(3*V**4+V**3-V**2+
     1   5*V-2)*W**3+2*(V-1)**2*V**2*(V**4+4*V**3-3*V**2+2)*W**2-2*(V-1)
     2   **3*V*(2*V**3+V**2-1)*W+2*(V-1)**4*(V**2+1))-4*V**8*W**8+4*V**7
     3   *(3*V+1)*W**7-V**4*(2*V-1)*(8*V**3+21*V**2-1)*W**6+2*V**4*(8*V*
     4   *4+29*V**3-V**2-26*V+2)*W**5-2*V**2*(8*V**6+17*V**5+34*V**4-52*
     5   V**3-5*V**2+3*V-1)*W**4+2*(V-1)*V**2*(6*V**5+9*V**4+36*V**3+28*
     6   V**2-20*V+3)*W**3-(V-1)*(4*V**7+12*V**6+5*V**5+33*V**4+8*V**3-3
     7   2*V**2+3*V-1)*W**2+2*(V-1)**2*(4*V**5-2*V**4+7*V**3-V**2+7*V+1)
     8   *W-4*(V-1)**4*(V**2+1))+N**4*(V-1)*VC*W*(V**2*W**2-2*(V-1)*V*W+
     9   (V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(4*V**8*W**7-4*V**7*(
     :   3*V-2)*W**6+V**4*(24*V**4-32*V**3+20*V**2-4*V+1)*W**5-V**3*(80*
     ;   V**5-196*V**4+212*V**3-120*V**2+33*V-4)*W**4+2*(V-1)*V**2*(82*V
     <   **5-222*V**4+268*V**3-162*V**2+41*V-3)*W**3-4*CQ*(V**2*W-(V-1)*
     =   *2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*(V-1)*V*(2*V*
     >   *2-2*V+1)*W+(V-1)**2)-2*(V-1)**2*V*(78*V**5-228*V**4+296*V**3-1
     ?   80*V**2+49*V-2)*W**2+(V-1)**3*(64*V**5-192*V**4+264*V**3-176*V*
     @   *2+57*V-1)*W-(V-1)**4*(8*V**4-20*V**3+32*V**2-24*V+13))-2*N**2*
     1   (V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(
     2   V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**7*W**6-V**6*(6
     3   *V-7)*W**5-CQ*(V*W-1)*(V**5*W**4-(V-1)*V**4*W**3-2*V**2*(3*V**2
     4   -5*V+1)*W**2-(V-1)*V*(V**3-8*V**2+9*V-4)*W+(V-1)**2*(V**3-2*V**
     5   2+5*V-2))+V**3*(2*V-1)*(5*V**3-9*V**2+V-1)*W**4-(V-1)*V**2*(14*
     6   V**4-26*V**3+10*V**2-11*V+1)*W**3+V*(12*V**6-35*V**5+32*V**4-13
     7   *V**3+8*V**2-4*V-1)*W**2-(V-1)**2*(4*V**5+V**4-16*V**3+14*V**2-
     8   10*V+1)*W+(V-1)**2*(4*V**4-13*V**3+17*V**2-16*V+7))-(V-1)*VC*W*
     9   (V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*
     :   (V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(6*V**5*W**4-V**2*(
     ;   14*V**3-2*V+1)*W**3+2*CQ*V*(V**2*W**2+(V-1)*V*W+(V-1)**2)*(V**2
     <   *W**2-2*V**2*W+V**2+1)+V*(12*V**4+6*V**3-18*V**2+7*V-2)*W**2-(V
     =   -1)*(2*V**4+10*V**3-2*V**2+5*V-1)*W-(V-1)**2*(2*V**3-2*V**2+4*V
     >   -1))-4*(CQ-1)*GTR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+(V-1)**
     ?   2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)
     @   *(V**3*W**3-3*V**2*W**2+3*V*W-1)+4*(CQ-1)*GTR*N*(V-1)**2*V**2*V
     1   C*W**2*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V
     2   **2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1))/(N**
     3   2*(V-1)**2*V**2*W**2*(V*W-1)**3*(V*W-V+1)**6)
 
 
      LVW = (-2*N**3*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V**4*W**4-2*V
     1   **3*(V+3)*W**3+V**2*(V**2+4*V+11)*W**2-2*V*(V**2+3*V+4)*W+4*(V*
     2   *2+1))*(V**6*W**6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V
     3   -1)**3*V**3*W**3+15*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)
     4   -2*N**4*(V-1)*VC*W*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(4*V**4*W**3
     5   +V**3*(4*V-9)*W**2+V**2*(8*V**2-16*V+13)*W-V**2+2*V-2)*(V**6*W*
     6   *6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**
     7   3+15*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)+2*N**2*(V-1)*V
     8   C*W*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(4*V**4*W**3+2*V**3*(2*V-3)
     9   *W**2+V**2*(8*V**2-17*V+12)*W-(V-2)*(V-1))*(V**6*W**6-6*(V-1)*V
     :   **5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(V-1)**
     ;   4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)+2*N*V*VC*W*(V**3*W**3-3*V*
     <   *2*W**2+3*V*W-1)*(2*V**3*W**3-2*V**2*(V+1)*W**2+V*(V+1)**2*W-2*
     =   (V-1)*(V+1))*(V**6*W**6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4
     >   -20*(V-1)**3*V**3*W**3+15*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-
     ?   1)**6)+2*(V-1)*VC*W*(V**3*W**2-V**2*(V+4)*W+2*(V**2-V+2))*(V**3
     @   *W**3-3*V**2*W**2+3*V*W-1)*(V**6*W**6-6*(V-1)*V**5*W**5+15*(V-1
     1   )**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(V-1)**4*V**2*W**2-6*(V
     2   -1)**5*V*W+(V-1)**6))*LOG(1-V*W)/(N**2*(V-1)**2*V**2*W**2*(V*W-
     3   1)**3*(V*W-V+1)**6)
 
 
      LTVW = (-2*N**4*(V-1)*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*
     1   W**3-3*V**2*W**2+3*V*W-1)*(V**7*W**6+2*V**6*(4*V-5)*W**5-V**4*(
     2   9*V**3-12*V**2+2*V-1)*W**4+4*(V-1)*V**3*(3*V**2-1)*W**3+4*CQ*(V
     3   **2*W-(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*(
     4   V-1)*V*(2*V**2-2*V+1)*W+(V-1)**2)-(V-1)**2*V**2*(9*V**3-6*V**2+
     5   7*V-6)*W**2+2*(V-1)**3*V*(4*V**3-5*V**2+7*V-2)*W+(V-1)**4*(V+1)
     6   *(V**2-V+1))+2*N**2*(V-1)*VC*W*(V**3*W**3-3*V**2*W**2+3*V*W-1)*
     7   (V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*
     8   W+(V-1)**4)*(2*CQ*(V**5*W**4-(V-1)*V**4*W**3-2*V**2*(3*V**2-5*V
     9   +1)*W**2-(V-1)*V*(V**3-8*V**2+9*V-4)*W+(V-1)**2*(V**3-2*V**2+5*
     :   V-2))+2*V**5*W**4-V**4*(9*V-5)*W**3+2*V**2*(7*V**3-6*V**2+1)*W*
     ;   *2-(V-1)*V*(9*V**3-2*V**2+V+4)*W+2*(V-1)**2*(V+1)*(V**2-V+1))+4
     <   *N*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W
     =   **3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(V**2*W**
     >   2-2*V**2*W+V**2+1)*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**
     ?   2-(V-1)**3*V*W+(V-1)**4)+(V-1)*V*W*(3*V**3*W**3-V**2*(7*V-1)*W*
     @   *2+(V-1)*V*(7*V+3)*W-3*(V-1)**2*(V+1)))-4*N**3*VC*(V**3*W**3-3*
     1   V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2
     2   *W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(V**2*W**2-2*V**2*W+V**2+1)*
     3   (V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1
     4   )**4)+2*(V-1)**2*V**3*(W-1)*W**2)-2*(V-1)*VC*W*(V**3*W**3-3*V**
     5   2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W*
     6   *2-4*(V-1)**3*V*W+(V-1)**4)*(-3*V**5*W**4+V**4*(6*V+1)*W**3+2*C
     7   Q*V*(V**2*W**2+(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)-
     8   V**2*(6*V**3-3*V**2+2*V+3)*W**2+(V-1)*V*(6*V**3-3*V**2-5*V+6)*W
     9   -(V-1)**2*(V+1)*(3*V**2-2*V+3))-8*(CQ-1)*GTR*N**3*(V-1)**2*V**2
     :   *VC*W**2*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*
     ;   (V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1)+8*(
     <   CQ-1)*GTR*N*(V-1)**2*V**2*VC*W**2*(V**2*W**2+(V-1)**2)*(V**2*W*
     =   *2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3
     >   -3*V**2*W**2+3*V*W-1))*LOG(V*W-V+1)/(N**2*(V-1)**2*V**2*W**2*(V
     ?   *W-1)**3*(V*W-V+1)**6)
 
 
      LW = (2*N**2*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**3*(4*V**2-5*V+2)*W
     1   **3-V*(2*V**4-2*V**3+5*V**2-2*V+1)*W**2+(6*V**5-11*V**4+12*V**3
     2   -5*V**2+4*V+1)*W-(2*V**2-V+2)*(3*V**2-5*V+3))*(V**6*W**6-6*(V-1
     3   )*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(V-1
     4   )**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)-N**3*VC*(V**3*W**3-3*V*
     5   *2*W**2+3*V*W-1)*(2*(2*V**4*W**4-V**3*(V+1)*(V+3)*W**3+V**2*(3*
     6   V**3-V**2+13*V-3)*W**2-2*V*(2*V-1)*(V**3-V**2+3*V+1)*W+2*(V-1)*
     7   (V**2+1)*(V**2-V+2))-CQ*(V**2+1)**2*(V*W-V+1)*(W**2-2*W+2))*(V*
     8   *5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**
     9   2*W**2+5*(V-1)**4*V*W-(V-1)**5)-(V-1)*VC*W*(V**2*W**2-2*V*W+1)*
     :   (4*V**5*W**4-V**2*(4*V**3+10*V**2+1)*W**3+2*V**2*(V**3+5*V**2+2
     ;   *V+1)*W**2+(2*V**5-8*V**4+V**2-5*V+1)*W-(V-1)*(2*V**3-6*V**2+8*
     <   V-5))*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-
     =   1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)-N*(V-1)*VC*(V**3*W**3-
     >   3*V**2*W**2+3*V*W-1)*(2*(V**3*(V+3)*W**3-V**2*(3*V**2+5)*W**2+2
     ?   *(V-1)*V*(2*V**2-V+1)*W-2*(V-1)**2*(V**2+1))+CQ*(V-1)*(V**2+1)*
     @   (V*W-V+1)*(W**2-2*W+2))*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**
     1   2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)-N**4
     2   *(V-1)*VC*W*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V**4*W**3-V*(12*
     3   V**4-14*V**3+6*V**2-4*V+1)*W**2+(24*V**5-68*V**4+84*V**3-50*V**
     4   2+12*V-1)*W-(V-1)*(12*V**4-24*V**3+28*V**2-18*V+7))*(V**5*W**5-
     5   5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5
     6   *(V-1)**4*V*W-(V-1)**5))*LOG(W)/(N**2*(V-1)**2*V**2*W**2*(V*W-1
     7   )**3*(V*W-V+1)**6)
 
 
      CVC = (-N**4*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W
     1   +(V-1)**2)*(6*V**9*W**8-6*V**8*(V+1)*W**7-2*V**6*(45*V**3-88*V*
     2   *2+42*V-3)*W**6+V**5*(250*V**4-678*V**3+688*V**2-291*V+16)*W**5
     3   -V**4*(234*V**5-754*V**4+884*V**3-376*V**2+5*V-6)*W**4+2*(V-1)*
     4   V**3*(33*V**5-24*V**4-138*V**3+181*V**2-46*V+8)*W**3+2*(V-1)**2
     5   *V**2*(7*V**5-66*V**4+129*V**3-107*V**2+83*V-7)*W**2-(V-1)**3*V
     6   **2*(6*V**4-4*V**3-38*V**2+27*V-19)*W+(V-1)**4*(6*V**4-12*V**3+
     7   8*V**2+9*V+2))+N**2*(V-1)*VC*W*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*
     8   W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(6*V**8*W**7
     9   +6*(V-5)*V**7*W**6-V**4*(84*V**4-174*V**3+79*V**2-15*V+6)*W**5+
     :   V**4*(120*V**4-204*V**3+69*V**2+19)*W**4-V**2*(42*V**6+54*V**5-
     ;   303*V**4+296*V**3-117*V**2+42*V-12)*W**3-V**2*(6*V**6-102*V**5+
     <   245*V**4-232*V**3+156*V**2-95*V+26)*W**2+(V-1)*(12*V**6-66*V**5
     =   +121*V**4-103*V**3+71*V**2-21*V+6)*W-(V-1)**2*(6*V**4-6*V**3-16
     >   *V**2+27*V-7))-3*N*VC*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V
     ?   -1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(4*V**7*W**7+V**4*(2*
     @   V**4-16*V**3+3*V**2+V-2)*W**6-V**4*(4*V**4-20*V**3-6*V**2+17*V-
     1   7)*W**5+2*V**2*(V**6-6*V**5-2*V**4+7*V**3-V**2-3*V+2)*W**4+2*(V
     2   -1)*V**2*(4*V**4+3*V**3-11*V**2+2*V+5)*W**3-(V-1)**2*(4*V**5+13
     3   *V**4-11*V**3-11*V**2-V+2)*W**2+(V-1)**3*(8*V**3+V**2-3)*W-4*(V
     4   -1)**4*V)+3*(V-1)*VC*W*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(
     5   V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**7*W**6-(V-2)*(
     6   V-1)*V**4*(2*V-1)*W**5-V**4*(19*V**2-17*V+7)*W**4-V**2*(2*V**5-
     7   29*V**4+30*V**3+V**2-18*V+4)*W**3+V**2*(2*V**5-13*V**4-5*V**3+5
     8   0*V**2-47*V+10)*W**2-(V-1)*(4*V**5-21*V**4+17*V**3+9*V**2-9*V+2
     9   )*W+(V-1)**2*(2*V**3-8*V**2+6*V-3))+3*AL*N**3*VC*(V*W-1)*(V**6*
     :   W**6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W
     ;   **3+15*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)*(2*V**6*W**6
     <   -2*V**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**
     =   5+3*V**4+10*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V*
     >   *3+36*V**2+4*V+1)*W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2
     ?   *(V**2+1)*(V**2+3))-3*AL*N*VC*(V*W-1)*(V**6*W**6-6*(V-1)*V**5*W
     @   **5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(V-1)**4*V**
     1   2*W**2-6*(V-1)**5*V*W+(V-1)**6)*(2*V**6*W**6-2*V**5*(V+3)*W**5+
     2   V**2*(2*V**4+2*V**3+11*V**2-2*V+1)*W**4-2*V*(V**5+2*V**3+3*V**2
     3   -V+1)*W**3+(V**2+1)*(2*V**4-3*V**2+2*V+1)*W**2-2*(V-1)**2*(2*V+
     4   1)*(V**2+1)*W+2*(V-1)**2*(V**2+1))+3*N**3*VC*(V*W-1)*(8*V**5*W*
     5   *5+V**2*(2*V**4-10*V**3-25*V**2+3*V-2)*W**4+V*(2*V**4+28*V**3+3
     6   3*V**2-3*V+4)*W**3-(8*V**5+3*V**4+37*V**3+27*V**2+3*V+2)*W**2+(
     7   16*V**4-11*V**3+39*V**2+9*V+3)*W-8*V*(V**2-V+2))*(V**6*W**6-6*(
     8   V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(
     9   V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)-4*GTR*N*(V-1)*VC*W*(
     :   V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)
     ;   *(V**6*W**5+V**4*(V**2-2*V-1)*W**4+2*(V-1)*V**3*(14*V**2-15*V+2
     <   )*W**3-2*(V-1)**2*V**2*(14*V**2+V+3)*W**2-(V-1)**3*V*(V**2+V-4)
     =   *W-(V-1)**4*(V**2+1))+4*CQ*N**3*VC*(V**2*W**2-2*V*W+1)*(V**5*(V
     >   +2)*W**5-V**4*(V**2+5*V-3)*W**4+V**2*(V**4+7*V**3-V**2-4*V+1)*W
     ?   **3-V**2*(V**4+3*V**3+4*V**2-5*V+1)*W**2+(V-1)*(4*V**4-6*V**3+1
     @   3*V**2-4*V+1)*W-(V-1)**2*(3*V**2-6*V+7))*(V**5*W**5-5*(V-1)*V**
     1   4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V
     2   *W-(V-1)**5)-4*CQ*N*VC*(V**2*W**2-2*V*W+1)*(V**5*(V+2)*W**5-V**
     3   4*(V**2+3*V-1)*W**4+(V-1)*V**2*(V**3+4*V**2+V-1)*W**3-(V-1)**2*
     4   V**2*(V**2+3*V+3)*W**2+(V-1)**3*(4*V**2+1)*W-3*(V-1)**4)*(V**5*
     5   W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W
     6   **2+5*(V-1)**4*V*W-(V-1)**5)-12*CQ*N**2*(V-1)*V*VC*W*(V**2*W**2
     7   +(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**5*W**5
     8   -5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+
     9   5*(V-1)**4*V*W-(V-1)**5)+12*CQ*(V-1)*V*VC*W*(V**2*W**2+(V-1)*V*
     :   W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**5*W**5-5*(V-1)*
     ;   V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**
     <   4*V*W-(V-1)**5)-4*GTR*N**3*(V-1)*VC*W*(V**2*W**2-2*(V-1)*V*W+(V
     =   -1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**6*(2*V-3)*W**5+V**4
     >   *(4*V**4-12*V**3+6*V**2+3*V+1)*W**4-(V-1)*V**3*(12*V**4-32*V**3
     ?   +45*V**2-37*V+4)*W**3+(V-1)**2*V**2*(12*V**4-36*V**3+53*V**2-11
     @   *V+6)*W**2-2*(V-1)**3*V*(2*V**4-7*V**3+5*V**2-6*V+2)*W+(V-1)**4
     1   *(V**2+1))+24*CQ*GTR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+(V-1
     2   )**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W
     3   **2+3*(V-1)**2*V*W-(V-1)**3)-24*CQ*GTR*N*(V-1)**2*V**2*VC*W**2*
     4   (V**2*W**2+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**3*W**3
     5   -3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3))/(N**2*(V-1)**2*V**
     6   2*W**2*(V*W-1)**3*(V*W-V+1)**6)/3.D0
 
 
      LM = LOG(S/M**2)*(N**3*VC*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(
     1   V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**6*W**6-2*V**5*
     2   (V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**5+3*V**4
     3   +10*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V**3+36*V*
     4   *2+4*V+1)*W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2*(V**2+1
     5   )*(V**2+3))-N*VC*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**
     6   2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**6*W**6-2*V**5*(V+3)*
     7   W**5+V**2*(2*V**4+2*V**3+11*V**2-2*V+1)*W**4-2*V*(V**5+2*V**3+3
     8   *V**2-V+1)*W**3+(V**2+1)*(2*V**4-3*V**2+2*V+1)*W**2-2*(V-1)**2*
     9   (2*V+1)*(V**2+1)*W+2*(V-1)**2*(V**2+1))+N**4*(V-1)*VC*W*(V*W-1)
     :   *(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V
     ;   *W+(V-1)**4)*(4*V**6*W**5+4*(V-4)*V**5*W**4+V**2*(12*V**4-32*V*
     <   *3+40*V**2-4*V+1)*W**3-V*(4*V**5+16*V**4-40*V**3+40*V**2-5*V+2)
     =   *W**2+(8*V**5-4*V**4-4*V**3+7*V**2+2*V+1)*W-4*V**4+8*V**3-9*V**
     >   2+6*V-3)-2*N**2*(V-1)*VC*W*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3
     ?   +6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**6*W**5+2*(
     @   V-4)*V**5*W**4+V**2*(6*V**4-16*V**3+22*V**2-3*V+1)*W**3-V*(2*V*
     1   *5+8*V**4-18*V**3+22*V**2-3*V+2)*W**2+(4*V**5-2*V**4+V**3+2*V**
     2   2+3*V+1)*W-2*V**4+4*V**3-6*V**2+5*V-3)+(V-1)*VC*W*(V*W-1)*(V**2
     3   *(4*V**2-2*V+1)*W**3-V*(4*V**3+4*V**2-V+2)*W**2+(6*V**3-3*V**2+
     4   4*V+1)*W-3*V**2+4*V-3)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*
     5   V**2*W**2-4*(V-1)**3*V*W+(V-1)**4))/(N**2*(V-1)**2*V**2*W**2*(V
     6   *W-1)**3*(V*W-V+1)**4)
 
 
      LMP = LOG(S/MP**2)*(4*N**4*(V-1)*VC*W*(V**3*W**3-3*V**2*W**2+3*V*W
     1   -1)*(V**8*W**7-V**7*(3*V-2)*W**6+V**6*(5*V**2-8*V+5)*W**5-V**4*
     2   (11*V**4-26*V**3+25*V**2-10*V+1)*W**4+(V-1)**2*V**3*(19*V**3-26
     3   *V**2+19*V-4)*W**3-(V-1)**2*V**2*(17*V**4-44*V**3+49*V**2-26*V+
     4   6)*W**2+(V-1)**4*V*(7*V**3-12*V**2+11*V-4)*W-(V-1)**4*(V**2-V+1
     5   )**2)-2*N**2*(V-1)*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*
     6   W**3-3*V**2*W**2+3*V*W-1)*(2*V**6*W**5-V**5*(6*V-5)*W**4+V**4*(
     7   8*V**2-13*V+7)*W**3-2*(V-1)*V**2*(4*V**3-6*V**2+6*V-1)*W**2+(V-
     8   1)*V*(6*V**4-15*V**3+20*V**2-13*V+4)*W-(V-1)**4*(2*V**2-V+2))+2
     9   *N**3*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V
     :   **2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-(V-1)*V**3*W*
     ;   *3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4)-2*N*VC*(V**2*W**2-
     <   2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*
     =   V**2*W**2+3*V*W-1)*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**
     >   2-(V-1)**3*V*W+(V-1)**4)+2*(V-1)*V*VC*W*(V**2*W**2-2*(V-1)*V*W+
     ?   (V-1)**2)*(V**2*W**2+(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V*
     @   *2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1)+4*GTR*N**3*(V-1)**2*V**2*
     1   VC*W**2*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*
     2   W**3-3*V**2*W**2+3*V*W-1)-4*GTR*N*(V-1)**2*V**2*VC*W**2*(V**2*W
     3   **2+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*V**2*W**
     4   2+3*V*W-1))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)**3*(V*W-V+1)**4)
 
      STRUV12=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV13(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(N**3*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V
     1   -1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*
     2   (4*V**8*W**8-4*V**7*(5*V-2)*W**7+V**4*(52*V**4-50*V**3+1
     3   7*V**2-4*V+1)*W**6-2*V**4*(42*V**4-59*V**3+30*V**2-7*V+2
     4   )*W**5+2*V**2*(44*V**6-64*V**5+15*V**4+25*V**3-16*V**2+5
     5   *V-1)*W**4-2*V**2*(28*V**6-22*V**5-60*V**4+114*V**3-71*V
     6   **2+20*V-3)*W**3+(16*V**8+40*V**7-204*V**6+280*V**5-150*
     7   V**4+12*V**3+15*V**2-6*V+1)*W**2-2*(V-1)*(16*V**6-28*V**
     8   5+6*V**4+30*V**3-29*V**2+10*V-1)*W+4*(V-1)**2*(2*V**2-3*
     9   V+2)*(2*V**2-2*V+1))-4*V**8*W**8+4*V**7*(4*V-1)*W**7-V**
     :   4*(36*V**4-30*V**3+13*V**2-4*V+1)*W**6+2*V**4*(28*V**4-3
     ;   9*V**3+22*V**2-5*V+2)*W**5-2*V**2*(32*V**6-46*V**5+V**4+
     <   37*V**3-20*V**2+5*V-1)*W**4+2*V**2*(24*V**6-20*V**5-66*V
     =   **4+130*V**3-79*V**2+20*V-3)*W**3-(16*V**8+32*V**7-204*V
     >   **6+292*V**5-138*V**4-12*V**3+23*V**2-6*V+1)*W**2+2*(V-1
     ?   )*(16*V**6-32*V**5+4*V**4+46*V**3-43*V**2+14*V-1)*W-4*(V
     @   -1)**2*(2*V**2-4*V+3)*(2*V**2-2*V+1))-N*(V-1)**2*VC*(V**
     1   2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*V*
     2   *2*W**2-4*V*W+1)*(CQ*(4*V**8*W**8-4*V**7*(4*V-1)*W**7+V*
     3   *4*(28*V**4-14*V**3+3*V**2-2*V+1)*W**6-2*V**4*(12*V**4-5
     4   *V**3+V**2-2*V+2)*W**5+2*V**2*(4*V**6+6*V**5-5*V**4+2*V*
     5   *3-V**2+3*V-1)*W**4-2*V**2*(6*V**5+2*V**3-6*V**2+7*V-3)*
     6   W**3+(12*V**6-10*V**4+2*V**3+3*V**2-4*V+1)*W**2-2*(V-1)*
     7   (6*V**4-5*V**2+4*V-1)*W+4*(V-1)**2*(2*V**2-2*V+1))-4*V**
     8   8*W**8+4*V**7*(4*V-1)*W**7-V**4*(28*V**4+6*V**3-17*V**2-
     9   2*V+1)*W**6+2*V**4*(12*V**4+35*V**3-53*V**2+12*V+2)*W**5
     :   -2*V**2*(4*V**6+72*V**5-111*V**4+42*V**3-V**2+3*V-1)*W**
     ;   4+2*V**2*(58*V**5-86*V**4+20*V**3+18*V**2-V-3)*W**3-(32*
     <   V**7-4*V**6-124*V**5+178*V**4-98*V**3+23*V**2-4*V+1)*W**
     =   2+2*(V-1)*(2*V**2-2*V+1)*(8*V**3-11*V**2+8*V-1)*W-4*(V-1
     >   )**2*(2*V**2-2*V+1))-2*N**4*(V-1)*VC*W*(V**2*W**2-2*V*W+
     ?   1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(
     @   V-1)**3*V*W+(V-1)**4)*(2*V**8*W**7-V**7*(V+5)*W**6+V**6*
     1   (5*V**2-6*V+11)*W**5-V**4*(3*V**4+10*V**3-2*V**2+15*V+4)
     2   *W**4+2*V**3*(2*V**4+12*V**3-2*V**2+11*V+8)*W**3-2*CQ*(V
     3   **2*W-1)*(V**2*W**2-2*V*W+1)*((V-1)**2*V**2*W**2-2*V*(V*
     4   *2+1)*W+(V-1)**2)-V**2*(10*V**4-2*V**3+21*V**2+5*V+24)*W
     5   **2+V*(4*V**4-5*V**3+18*V**2-11*V+16)*W-(V**2-V+1)*(3*V*
     6   *2-5*V+4))+2*N**2*(V-1)*VC*W*(V**3*W**3-3*(V-1)*V**2*W**
     7   2+3*(V-1)**2*V*W-(V-1)**3)*(V**4*W**4-4*V**3*W**3+6*V**2
     8   *W**2-4*V*W+1)*(2*V**7*W**6+2*V**5*(V**2-3*V+1)*W**5-CQ*
     9   (V-1)*(V*W-V+1)*(V**5*W**4-V**4*W**3+2*(V-1)*V**2*(V**2-
     :   3*V-1)*W**2+V*(2*V**3-2*V**2+3*V-4)*W+2*V**3-2*V**2-V+2)
     ;   -V**4*(7*V**3-6*V**2-8*V+3)*W**4+V**3*(6*V**4+5*V**3-26*
     <   V**2+7*V-4)*W**3-V**2*(4*V**5+2*V**3-29*V**2+17*V-4)*W**
     =   2+V*(6*V**5-4*V**4-15*V**3+16*V**2-11*V+6)*W-(V-1)*(4*V*
     >   *4-4*V**3-5*V**2+8*V-5))+2*(V-1)**2*VC*W*(V**3*W**3-3*(V
     ?   -1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V**4*W**4-4*V**3
     @   *W**3+6*V**2*W**2-4*V*W+1)*(V**6*W**5-V**4*(4*V**2+V-1)*
     1   W**4+2*V**3*(3*V-2)*(V**2+V+1)*W**3+CQ*V*(V*W-V+1)*(V**2
     2   *W**2+V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)-2*V**2*(V
     3   **4+8*V**3-11*V**2+6*V-2)*W**2+V**2*(10*V**3-18*V**2+14*
     4   V-5)*W-(V-1)*(2*V**3-2*V**2+2*V-1))+4*(CQ-1)*GTR*N**3*(V
     5   -1)**2*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2-2*V*W+1)*(V
     6   **2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*
     7   W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)-4*(CQ
     8   -1)*GTR*N*(V-1)**2*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2
     9   -2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-4
     :   *(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-
     ;   1)**4))/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**6*(V*W-V+1)**4
     <   )
 
      LV1 = LOG(1-V)*(2*N*(V-1)**2*VC*(V**4*W**4-4*V**3*W**3+6*V
     1   **2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**
     2   2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*(V*W-1)*(V**5*W*
     3   *5-2*V**5*W**4+2*V**3*(V**2-2*V+2)*W**3+2*V**3*(2*V-3)*W
     4   **2-2*V*(2*V-3)*(2*V**2-2*V+1)*W-2*V**2+2*V-1)+CQ*(V**2*
     5   W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-V**3*W**3+V**2*W*
     6   *2-V*W+1))-2*N**4*(V-1)*VC*W*(V**4*W**4-4*V**3*W**3+6*V*
     7   *2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2
     8   *V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*((V*W-1)*(4*V**5*W**
     9   4-4*(V-1)*V**4*W**3+V**3*(5*V**2-6*V+9)*W**2-V*(V+1)*(V*
     :   *3-V**2+3*V+1)*W+(V+1)**2*(V**2-V+1))+2*CQ*(V**2*W-1)*((
     ;   V-1)**2*V**2*W**2-2*V*(V**2+1)*W+(V-1)**2))+2*N**2*(V-1)
     <   *VC*W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)*
     =   *3)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*((V**2*W
     >   **2-2*V*W+1)*(4*V**5*W**4-7*(V-1)*V**4*W**3+V**2*(4*V**3
     ?   -8*V**2+13*V-1)*W**2+V*(4*V**4-15*V**3+9*V**2-4*V+2)*W-(
     @   V-1)*(4*V**4-14*V**3+14*V**2-9*V+1))+CQ*(V-1)*(V*W-V+1)*
     1   (V**5*W**4-V**4*W**3+2*(V-1)*V**2*(V**2-3*V-1)*W**2+V*(2
     2   *V**3-2*V**2+3*V-4)*W+2*V**3-2*V**2-V+2))-2*(V-1)**2*VC*
     3   W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*
     4   (V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*((V*W-1)*(V*
     5   *5*W**4-(V-1)**2*V**3*W**3-(V-1)*V**2*(5*V**2-5*V+6)*W**
     6   2+V*(4*V**4-7*V**3+2*V**2+V-1)*W-(V-1)*(4*V**3-12*V**2+1
     7   1*V-4))+CQ*V*(V*W-V+1)*(V**2*W**2+V*W+1)*(V**2*W**2-2*V*
     8   *2*W+2*V**2-2*V+1))-2*N**3*(V-1)**2*VC*(V**4*W**4-4*V**3
     9   *W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+
     :   6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(V**2*
     ;   W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-V**3*W**3+V**2*W*
     <   *2-V*W+1)+2*V*(W-1)*(V*W-1)*(V**4*W**4-V**3*(2*V-1)*W**3
     =   +V**2*(2*V**2-3*V+2)*W**2-V**2*W+2*V**2-2*V+1))-4*CQ*GTR
     >   *N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2-2*V
     ?   *W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-
     @   1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**
     1   4)+4*CQ*GTR*N*(V-1)**2*V**2*VC*W**2*(V**2*W**2+1)*(V**2*
     2   W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W*
     3   *4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W
     4   +(V-1)**4))/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**6*(V*W-V+1
     5   )**4)
 
      LV = LOG(V)*(N**3*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)
     1   **2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*(4*
     2   V**8*W**8-4*V**7*(5*V-2)*W**7+3*V**6*(2*V-1)*(8*V-3)*W**
     3   6-2*V**5*(34*V**3-43*V**2+14*V+1)*W**5+V**4*(60*V**4-80*
     4   V**3-2*V**2+46*V-15)*W**4-2*V**3*(8*V-5)*(2*V**4-6*V**2+
     5   6*V-1)*W**3+2*V**2*(4*V**6+12*V**5-64*V**4+88*V**3-42*V*
     6   *2+3)*W**2-2*(V-1)*V*(8*V**5-16*V**4+2*V**3+22*V**2-19*V
     7   +5)*W+2*(V-1)**2*(2*V**2-4*V+3)*(2*V**2-2*V+1))-4*V**8*W
     8   **8+4*V**7*(5*V-2)*W**7-V**4*(52*V**4-50*V**3+17*V**2-4*
     9   V+1)*W**6+2*V**4*(42*V**4-59*V**3+24*V**2-V+2)*W**5-2*V*
     :   *2*(44*V**6-64*V**5-9*V**4+59*V**3-26*V**2+5*V-1)*W**4+2
     ;   *V**2*(28*V**6-22*V**5-98*V**4+180*V**3-101*V**2+22*V-3)
     <   *W**3-(16*V**8+40*V**7-260*V**6+380*V**5-178*V**4-16*V**
     =   3+27*V**2-6*V+1)*W**2+2*(V-1)*(16*V**6-36*V**5+2*V**4+62
     >   *V**3-57*V**2+18*V-1)*W-4*(V-1)**2*(2*V**2-5*V+4)*(2*V**
     ?   2-2*V+1))-N*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)
     @   *(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*(4*V**8
     1   *W**8-4*V**7*(4*V-1)*W**7+V**6*(28*V**2-14*V+1)*W**6-2*V
     2   **5*(12*V**3-5*V**2-3*V+2)*W**5+V**4*(8*V**4+12*V**3-24*
     3   V**2+14*V-1)*W**4-2*V**3*(6*V**4-6*V**3+2*V**2+3*V-2)*W*
     4   *3+2*V**2*(4*V**4-6*V**3+9*V**2-8*V+2)*W**2-2*(V-1)*V*(2
     5   *V**3+2*V**2-3*V+1)*W+2*(V-1)**2*(2*V**2-2*V+1))-4*V**8*
     6   W**8+4*V**7*(4*V-1)*W**7-V**4*(V+1)*(28*V**3-18*V**2-3*V
     7   +1)*W**6+2*V**4*(12*V**4+43*V**3-67*V**2+18*V+2)*W**5-2*
     8   V**2*(4*V**6+82*V**5-137*V**4+62*V**3-5*V**2+3*V-1)*W**4
     9   +2*V**2*(62*V**5-100*V**4+30*V**3+22*V**2-5*V-3)*W**3-(3
     :   2*V**7-4*V**6-144*V**5+230*V**4-142*V**3+35*V**2-4*V+1)*
     ;   W**2+2*(V-1)*(16*V**5-42*V**4+56*V**3-37*V**2+12*V-1)*W-
     <   4*(V-1)**2*(2*V**2-2*V+1))-N**4*(V-1)*VC*W*(V**2*W**2-2*
     =   V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2
     >   -4*(V-1)**3*V*W+(V-1)**4)*(4*V**8*W**7-4*V**7*(V+2)*W**6
     ?   +V**4*(9*V**4+14*V**2+1)*W**5-V**3*(5*V**5+20*V**4+10*V*
     @   *3+28*V**2+13*V+4)*W**4+2*V**2*(4*V**5+23*V**4-2*V**3+28
     1   *V**2+26*V+3)*W**3-4*CQ*(V**2*W-1)*(V**2*W**2-2*V*W+1)*(
     2   (V-1)**2*V**2*W**2-2*V*(V**2+1)*W+(V-1)**2)-2*V*(13*V**5
     3   -10*V**4+30*V**3+4*V**2+39*V+2)*W**2+(16*V**5-31*V**4+68
     4   *V**3-42*V**2+52*V+1)*W-9*V**4+24*V**3-38*V**2+28*V-13)+
     5   2*N**2*(V-1)*VC*W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**
     6   2*V*W-(V-1)**3)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W
     7   +1)*(2*V**7*W**6+(V-7)*V**6*W**5-CQ*(V-1)*(V*W-V+1)*(V**
     8   5*W**4-V**4*W**3+2*(V-1)*V**2*(V**2-3*V-1)*W**2+V*(2*V**
     9   3-2*V**2+3*V-4)*W+2*V**3-2*V**2-V+2)-V**3*(V+1)*(4*V**3-
     :   10*V**2+2*V-1)*W**4+V**2*(12*V**4-35*V**3+17*V**2-7*V-1)
     ;   *W**3-V*(V**6-4*V**5+14*V**4-41*V**3+27*V**2-10*V+1)*W**
     <   2+(6*V**5-24*V**4+24*V**3-16*V**2+5*V+1)*W-(V-1)*(V**4+V
     =   **3-11*V**2+12*V-7))+(V-1)**2*VC*W*(V**4*W**4-4*V**3*W**
     >   3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V
     ?   -1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(6*V**5*W**4-V
     @   **2*(13*V**3+V**2+V-1)*W**3+2*CQ*V*(V**2*W**2+V*W+1)*(V*
     1   *2*W**2-2*V**2*W+2*V**2-2*V+1)+V*(5*V**4+17*V**3-9*V**2+
     2   V-2)*W**2-(14*V**4-17*V**3+7*V**2-V-1)*W-3*V**3+3*V**2-V
     3   -1)+4*(CQ-1)*GTR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+1
     4   )*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*
     5   (V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1
     6   )**3*V*W+(V-1)**4)-4*(CQ-1)*GTR*N*(V-1)**2*V**2*VC*W**2*
     7   (V**2*W**2+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*
     8   V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2
     9   *W**2-4*(V-1)**3*V*W+(V-1)**4))/(N**2*(V-1)**3*V**2*W**2
     :   *(V*W-1)**6*(V*W-V+1)**4)
 
      LVW = (4*N**3*(V-1)**2*VC*(V**4*W**4-V**3*(2*V-1)*W**3+V**2
     1   *(2*V**2-4*V+3)*W**2+V*(2*V**2-4*V+1)*W+2*V**2-2*V+1)*(V
     2   **4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)*
     3   *3*V*W+(V-1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*
     4   V**3*W**3+15*V**2*W**2-6*V*W+1)-4*N*(V-1)**2*VC*(V**4*W*
     5   *4-V**3*(2*V-1)*W**3+V**3*(2*V-1)*W**2-V*(4*V**2-5*V+2)*
     6   W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*
     7   V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**6*W**6-6*V**5*W**
     8   5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1)-2*N**2
     9   *(V-1)*VC*W*(4*V**4*W**3-4*V**3*W**2-V**2*(3*V-7)*W+2*(2
     :   *V**4-5*V**3+3*V**2-V-1))*(V**4*W**4-4*(V-1)*V**3*W**3+6
     ;   *(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**6*W**6-
     <   6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W
     =   +1)+2*N**4*(V-1)*VC*W*(4*V**4*W**3-V**3*(3*V-7)*W**2+2*V
     >   **2*(3*V**2-4*V+7)*W-(V**2-V+1)*(2*V**2-V+3))*(V**4*W**4
     ?   -4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(
     @   V-1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**
     1   3+15*V**2*W**2-6*V*W+1)-2*(V-1)**2*VC*W*(V**3*W**2-V**2*
     2   (3*V-1)*W+4*V**3-8*V**2+8*V-3)*(V**4*W**4-4*(V-1)*V**3*W
     3   **3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**6*
     4   W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-
     5   6*V*W+1))*LOG(1-V*W)/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**
     6   6*(V*W-V+1)**4)
 
      LTVW = (2*N**3*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)
     1   *(CQ*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4
     2   -4*(V-1)*V**3*W**3+8*(V-1)**2*V**2*W**2-8*(V-1)**3*V*W+4
     3   *(V-1)**4)+(V-1)**2*V**2*W*((2*V-1)*W-2*(V-1)))*(V**6*W*
     4   *6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*
     5   V*W+1)-2*N**2*(V-1)*VC*W*(4*V**4*W**3-2*(V-3)*V**3*W**2+
     6   V**2*(3*V**2-7*V+12)*W+(V-2)*(V-1)**2)*(V**4*W**4-4*(V-1
     7   )*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4
     8   )*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V*
     9   *2*W**2-6*V*W+1)+2*N**4*(V-1)*VC*W*(4*V**4*W**3-V**3*(5*
     :   V-9)*W**2+V**2*(5*V**2-10*V+13)*W-(V-1)**2*(V**2-2*V+2))
     ;   *(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-
     <   1)**3*V*W+(V-1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-
     =   20*V**3*W**3+15*V**2*W**2-6*V*W+1)-2*(V-1)**2*VC*W*(V**3
     >   *W**2-V**2*(5*V-4)*W+2*(V-1)*(2*V**2-3*V+2))*(V**4*W**4-
     ?   4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V
     @   -1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3
     1   +15*V**2*W**2-6*V*W+1)-2*N*(V-1)**2*V*VC*W*(V**2*W**2-2*
     2   (V-1)*V*W+(V-1)**2)*((V-1)*(2*V**3*W**3-2*V**2*(3*V-2)*W
     3   **2+(V-1)*V*(8*V-5)*W-2*(V-1)**2*(2*V-1))+CQ*V*W*(V**2*W
     4   **2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2
     5   *V**2-2*V+1))*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**
     6   3*W**3+15*V**2*W**2-6*V*W+1))*LOG(V*W-V+1)/(N**2*(V-1)*
     7   *3*V**2*W**2*(V*W-1)**6*(V*W-V+1)**4)
 
      LW = (-2*N**2*(V-1)*VC*W*(V**3*(V**2+V+2)*W**3-V*(4*V**4-6*
     1   V**3+5*V**2-2*V+1)*W**2+(7*V**5-19*V**4+31*V**3-21*V**2+
     2   9*V-1)*W-(V-1)*(V**2-V+3)*(3*V**2-3*V+2))*(V**3*W**3-3*(
     3   V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V**6*W**6-6*V**
     4   5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1)+N
     5   **3*(V-1)**2*VC*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*
     6   V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*(2*(V-1)*V**4*W**4
     7   -V**3*(2*V-1)*(4*V-3)*W**3+V**2*(12*V**3-16*V**2+4*V+3)*
     8   W**2-2*V*(V+1)*(4*V**3-8*V**2+6*V-1)*W+2*(2*V**2-3*V+2)*
     9   (2*V**2-2*V+1))-CQ*(2*V**2-2*V+1)**2*(V*W-1)*(W**2-2*W+2
     :   ))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*
     ;   W-1)+(V-1)**2*VC*W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)*
     <   *2*V*W-(V-1)**3)*(4*V**5*W**4-V**2*(15*V**3-13*V**2+3*V-
     =   1)*W**3+2*V**2*(9*V**3-12*V**2+5*V-1)*W**2-(9*V**5-20*V*
     >   *4+17*V**3-9*V**2+1)*W+(V-1)*(V**3-5*V**2+7*V-5))*(V**5*
     ?   W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)+N*(V
     @   -1)**2*VC*(2*(V**3*(4*V-3)*W**3-V**2*(8*V**2-10*V+5)*W**
     1   2+2*V*(2*V**2-V+1)*W-2*(2*V**2-2*V+1))+CQ*(2*V**2-2*V+1)
     2   *(V*W-1)*(W**2-2*W+2))*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V
     3   -1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**5*W**5-5*V
     4   **4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)+N**4*(V-1)*V
     5   C*W*(2*(V-1)*V**4*W**3-V*(V**4+10*V**3+1)*W**2+(V**5+7*V
     6   **4-4*V**3+12*V**2+7*V+1)*W-5*V**4+6*V**3-16*V**2+10*V-7
     7   )*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V
     8   -1)**3*V*W+(V-1)**4)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3
     9   -10*V**2*W**2+5*V*W-1))*LOG(W)/(N**2*(V-1)**3*V**2*W**2
     :   *(V*W-1)**6*(V*W-V+1)**4)
 
      CVC = (N**4*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-
     1   1)*V*W+(V-1)**2)*(6*V**10*W**9-6*V**9*(3*V+2)*W**8+2*V**
     2   7*(10*V**3+2*V**2-6*V-3)*W**7-V**6*(23*V**4-97*V**3+145*
     3   V**2-167*V-10)*W**6+2*V**5*(14*V**5-113*V**4+297*V**3-37
     4   8*V**2+106*V+5)*W**5-V**4*(13*V**6-159*V**5+668*V**4-107
     5   2*V**3+475*V**2+11*V+22)*W**4-2*V**3*(14*V**6-164*V**5+4
     6   52*V**4-474*V**3+226*V**2-62*V+1)*W**3-V**2*(78*V**6-406
     7   *V**5+835*V**4-877*V**3+503*V**2-113*V-14)*W**2-2*(V-1)*
     8   V*(14*V**5-56*V**4+79*V**3-60*V**2+18*V-1)*W-(V-1)**2*(1
     9   3*V**4-33*V**3+35*V**2-11*V+2))-N**2*(V-1)*VC*W*(V**2*W*
     :   *2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W
     ;   **2-4*V*W+1)*(6*V**8*W**7-6*V**7*(4*V-1)*W**6+V**4*(20*V
     <   **4+11*V**3-70*V**2+9*V-6)*W**5+V**4*(4*V**4-70*V**3+207
     =   *V**2-64*V+19)*W**4-V**2*(2*V**6-42*V**5+111*V**4+76*V**
     >   3-123*V**2+30*V-12)*W**3-V**2*(4*V**6+15*V**5-51*V**4-46
     ?   *V**3+119*V**2-61*V+26)*W**2+(V-1)*(20*V**6-58*V**5+82*V
     @   **4-91*V**3+68*V**2-15*V+6)*W-(V-1)**2*(4*V**4-3*V**3-V*
     1   *2+13*V-7))+3*N*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)
     2   **2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(4*(V-1
     3   )*V**7*W**7-V**4*(12*V**4-15*V**3+6*V**2-7*V+2)*W**6+V**
     4   4*(12*V**4-9*V**3-3*V**2-11*V+7)*W**5-2*V**2*(2*V**6+4*V
     5   **5-13*V**4+13*V**3-14*V**2+9*V-2)*W**4+2*(V-1)*V**2*(3*
     6   V**4-7*V**3+25*V**2-22*V+5)*W**3+(V-1)*(4*V**5-36*V**4+3
     7   0*V**3+5*V**2-9*V+2)*W**2+(V-1)**2*(6*V**3+8*V**2-9*V+3)
     8   *W-4*(V-1)**3*V)+3*AL*N**3*(V-1)**2*VC*(V**4*W**4-4*V**3
     9   *W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+
     :   6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**6*W*
     ;   *6-2*V**5*(2*V+1)*W**5+V**2*(8*V**3-4*V**2+4*V-1)*W**4+2
     <   *V*(2*V**2-V+1)*(2*V**3-2*V**2-2*V+1)*W**3-(8*V**6-16*V*
     =   *4+16*V**3-10*V**2+1)*W**2+2*(V-1)*(8*V**4-4*V**3+2*V**2
     >   +2*V-1)*W-4*(V-1)*V*(2*V**2-2*V+1))-3*N**3*(V-1)**2*VC*(
     ?   V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(8*(V-1)*V**5*W**5-V**2*
     @   (32*V**4-59*V**3+28*V**2-5*V+2)*W**4+(V-1)*V*(64*V**4-10
     1   1*V**3+48*V**2-13*V+4)*W**3-(V-1)*(80*V**5-180*V**4+156*
     2   V**3-59*V**2+13*V-2)*W**2+(V-1)**2*(56*V**4-106*V**3+84*
     3   V**2-21*V+3)*W-8*(V-1)**3*V*(2*V**2-3*V+2))*(V**6*W**6-6
     4   *V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+
     5   1)-4*CQ*N**3*(V-1)**2*(2*V**2-2*V+1)**2*VC*W*(V**4*W**4-
     6   4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V
     7   -1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3
     8   +15*V**2*W**2-6*V*W+1)+4*CQ*N*(V-1)**2*(2*V**2-2*V+1)*VC
     9   *W*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(
     :   V-1)**3*V*W+(V-1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**
     ;   4-20*V**3*W**3+15*V**2*W**2-6*V*W+1)-4*GTR*N**3*(V-1)*VC
     <   *W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V
     =   -1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*((V-3)*(V-1)*V
     >   **6*W**5-V**4*(2*V**4-13*V**3+21*V**2-7*V+1)*W**4-V**3*(
     ?   8*V**4-37*V**3+42*V**2-21*V-4)*W**3-V**2*(24*V**4-61*V**
     @   3+56*V**2-13*V+6)*W**2-2*V*(4*V**4-7*V**3+V**2+2*V-2)*W-
     1   (V-1)**2*(2*V**2-2*V+1))-3*(V-1)**2*VC*W*(V**3*W**3-3*(V
     2   -1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V**4*W**4-4*V**3
     3   *W**3+6*V**2*W**2-4*V*W+1)*(2*V**6*W**5+V**3*(2*V**3-V**
     4   2-V-2)*W**4-(V-1)*V**2*(7*V**3+8*V+2)*W**3+V*(3*V**5-8*V
     5   **4+20*V**3-13*V**2-6*V+2)*W**2+(2*V**5-7*V**4-V**3+12*V
     6   **2-2*V-2)*W+(V-1)*(3*V**3-7*V**2+5*V-3))-3*AL*N*(V-1)**
     7   2*VC*W*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4
     8   *W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*
     9   V*W+(V-1)**4)*(2*V**6*W**5-2*V**5*(2*V+1)*W**4+V**2*(4*V
     :   **4+2*V**2+2*V-1)*W**3-2*V*(2*V**4-2*V**3+2*V**2+V-1)*W*
     ;   *2-(4*V**3-8*V**2+2*V+1)*W+2*(V-1)*(2*V**2-1))+4*GTR*N*(
     <   V-1)**3*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3
     =   *W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**
     >   6*W**5-V**4*(2*V**2-4*V+1)*W**4+2*V**3*(V**2+11*V+2)*W**
     ?   3-2*V**2*(18*V**2-7*V+3)*W**2+V*(2*V**2-7*V+4)*W-2*V**2+
     @   2*V-1)-6*AL*N**2*(V-1)**2*V*VC*W*(V**2*W**2+V*W+1)*(V**2
     1   *W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-4*V**3*W**3+6*V*
     2   *2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2
     3   *V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)+6*AL*(V-1)**2*V*VC*W
     4   *(V**2*W**2+V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V*
     5   *4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*(V
     6   -1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)*
     7   *4)+12*AL*GTR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+1)*(
     8   V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V*
     9   *4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**
     :   3*V*W+(V-1)**4)-12*AL*GTR*N*(V-1)**2*V**2*VC*W**2*(V**2*
     ;   W**2+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2
     <   *V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-
     =   4*(V-1)**3*V*W+(V-1)**4))/(N**2*(V-1)**3*V**2*W**2*(V*W-
     >   1)**6*(V*W-V+1)**4)/3.D0
 
      LM = LOG(S/M**2)*(-N**4*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)*
     1   *2)*(4*V**8*W**7-4*V**7*(V+2)*W**6+V**4*(9*V**4-8*V**3+2
     2   2*V**2+1)*W**5-V**3*(5*V**5+12*V**4+6*V**3+32*V**2+5*V+4
     3   )*W**4+2*V**2*(2*V**5+19*V**4+4*V**3+20*V**2+10*V+3)*W**
     4   3-2*V*(7*V**5+2*V**4+20*V**3+8*V**2+15*V+2)*W**2+(4*V**5
     5   +9*V**4+16*V**3-2*V**2+20*V+1)*W-5*V**4+8*V**3-14*V**2+8
     6   *V-5)-N**3*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-
     7   1)*V*W+(V-1)**2)*(2*V**6*W**6-2*V**5*(2*V+1)*W**5+V**2*(
     8   8*V**4-8*V**3+12*V**2-4*V+1)*W**4-2*V*(4*V**5-2*V**4+6*V
     9   **2-3*V+1)*W**3+(8*V**6-8*V**4+16*V**3-2*V**2+1)*W**2-2*
     :   (4*V**2-2*V+1)*(2*V**3-2*V**2+V+1)*W+4*(V**2-V+1)*(2*V**
     ;   2-2*V+1))+N*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V
     <   -1)*V*W+(V-1)**2)*(2*V**6*W**6-2*V**5*(2*V+1)*W**5+V**2*
     =   (4*V**4+6*V**2-2*V+1)*W**4-2*V*(2*V**4+2*V**3+2*V**2-V+1
     >   )*W**3+(2*V+1)*(4*V**3+1)*W**2-2*(6*V**3-2*V**2+V+1)*W+4
     ?   *(2*V**2-2*V+1))+2*N**2*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W
     @   **2-2*(V-1)*V*W+(V-1)**2)*(2*V**6*W**5-V**5*(V+5)*W**4+V
     1   **2*(3*V**4-2*V**3+9*V**2-V+1)*W**3-V*(V**5+7*V**4-6*V**
     2   3+9*V**2+V+2)*W**2+(4*V**5-3*V**4+8*V**3-3*V**2+5*V+1)*W
     3   -V**4+V**3-5*V**2+4*V-3)-(V-1)*VC*W*(V**2*W**2-2*V*W+1)*
     4   (V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**5*W**4-V**2*(3*V*
     5   *3-V**2-V+1)*W**3+(V-1)*V*(3*V**3-6*V**2-V-2)*W**2+(6*V*
     6   *4-9*V**3+3*V**2-V-1)*W+3*V**3-3*V**2+V+1)-4*GTR*N**3*(V
     7   -1)*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2-2*(V-1)*V*W+(V
     8   -1)**2)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)+4*GTR*N*(V-1)*
     9   V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)*
     :   *2)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1))/(N**2*(V-1)**2*V*
     ;   *2*W**2*(V*W-1)**4*(V*W-V+1)**2)
 
      LMP = LOG(S/MP**2)*(-N**4*VC*W*(V**4*W**4-4*V**3*W**3+6*V*
     1   *2*W**2-4*V*W+1)*(4*V**6*W**5-4*V**5*(3*V-4)*W**4+8*V**4
     2   *(2*V**2-5*V+4)*W**3-2*(V-1)*V**2*(7*V**3-15*V**2+17*V-1
     3   )*W**2+(V-1)**2*V*(7*V**3-14*V**2+19*V-4)*W-(V-1)**4*(V*
     4   *2-2*V+2))+2*N**2*VC*W*(V**4*W**4-4*V**3*W**3+6*V**2*W**
     5   2-4*V*W+1)*(2*V**6*W**5-2*V**5*(3*V-4)*W**4+V**4*(9*V**2
     6   -22*V+17)*W**3-(V-1)*V**2*(10*V**3-21*V**2+21*V-2)*W**2+
     7   (V-1)**2*V*(6*V**3-13*V**2+15*V-4)*W-(V-1)**4*(V**2-2*V+
     8   2))-N**3*(V-1)*VC*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+
     9   1)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**
     :   4-4*(V-1)*V**3*W**3+8*(V-1)**2*V**2*W**2-8*(V-1)**3*V*W+
     ;   4*(V-1)**4)-(V-1)**2*VC*W*(2*V**4*W**3-2*V**2*(3*V**2-3*
     <   V+1)*W**2+(V-1)*V*(5*V**2-7*V+4)*W-(V-1)**2*(V**2-2*V+2)
     =   )*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)+N*(V-1)*V*
     >   *2*VC*W**2*(V**2*W**2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W*
     ?   *2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4-4*V**3*W**3+6*
     @   V**2*W**2-4*V*W+1))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)**4*
     1   (V*W-V+1)**2)
 
      STRUV13=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV14(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W
     1   **2-2*(V-1)*V*W+(V-1)**2)*(CQ*(4*V**12*W**10-4*V**10*(5*
     2   V**2-2*V+2)*W**9+4*V**6*(13*V**6-13*V**5+19*V**4-14*V**3
     3   +12*V**2-6*V+2)*W**8-4*V**6*(23*V**6-37*V**5+67*V**4-62*
     4   V**3+36*V**2-6*V+2)*W**7+V**4*(120*V**8-227*V**7+440*V**
     5   6-378*V**5+45*V**4+240*V**3-192*V**2+96*V-24)*W**6-V**6*
     6   (112*V**6-178*V**5+275*V**4-24*V**3-413*V**2+510*V-170)*
     7   W**5+V**2*(64*V**10-7*V**9-194*V**8+750*V**7-1245*V**6+1
     8   020*V**5-228*V**4-240*V**3+240*V**2-120*V+24)*W**4-V**2*
     9   (16*V**10+96*V**9-373*V**8+838*V**7-1161*V**6+1004*V**5-
     :   596*V**4+368*V**3-272*V**2+120*V-24)*W**3+(V-1)*(48*V**1
     ;   0-48*V**9-29*V**8+325*V**7-662*V**6+793*V**5-563*V**4+20
     <   8*V**3+8*V**2-40*V+8)*W**2-(V-1)**2*(48*V**8-128*V**7+27
     =   5*V**6-392*V**5+457*V**4-358*V**3+194*V**2-64*V+16)*W+16
     >   *(V-1)**3*(V**6-3*V**5+7*V**4-11*V**3+13*V**2-9*V+3))-8*
     ?   V**12*W**10+20*V**10*(2*V**2-V+1)*W**9-4*V**6*(26*V**6-3
     @   1*V**5+34*V**4-8*V**3+9*V**2-6*V+2)*W**8+2*V**6*(92*V**6
     1   -171*V**5+184*V**4-14*V**3-23*V**2+36*V-12)*W**7-V**4*(2
     2   40*V**8-529*V**7+492*V**6+166*V**5-337*V**4+372*V**3-236
     3   *V**2+96*V-24)*W**6+V**4*(224*V**8-444*V**7+113*V**6+786
     4   *V**5-607*V**4-12*V**3+452*V**2-384*V+96)*W**5-V**2*(128
     5   *V**10-73*V**9-726*V**8+1934*V**7-1699*V**6+552*V**5+432
     6   *V**4-672*V**3+348*V**2-120*V+24)*W**4+V**2*(32*V**10+17
     7   6*V**9-951*V**8+1822*V**7-1601*V**6+976*V**5-932*V**4+12
     8   40*V**3-1210*V**2+600*V-120)*W**3-(V-1)*(96*V**10-160*V*
     9   *9-191*V**8+1091*V**7-1554*V**6+1303*V**5-565*V**4+64*V*
     :   *3+44*V**2-40*V+8)*W**2+(V-1)**2*(96*V**8-336*V**7+621*V
     ;   **6-650*V**5+573*V**4-432*V**3+368*V**2-192*V+48)*W-32*(
     <   V-1)**3*(V**2-2*V+2)*(V**2-V+1)**2)+2*N**2*(V-1)*V**2*VC
     =   *W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*
     >   (CQ*(2*V**10*W**9-2*V**8*(5*V**2-2*V+2)*W**8+V**8*(22*V*
     ?   *2-13*V+13)*W**7-V**6*(26*V**4-11*V**3+15*V**2-8*V+4)*W*
     @   *6+2*V**6*(8*V**4+5*V**3+3*V**2-16*V+8)*W**5-V**4*(4*V**
     1   6+21*V**5+10*V**4-72*V**3+61*V**2-30*V+10)*W**4+(V-1)*V*
     2   *4*(9*V**4+32*V**3-35*V**2+6*V-3)*W**3-(V-1)**2*V**2*(13
     3   *V**4+15*V**3-19*V**2+8*V-4)*W**2+(V-1)**3*V**2*(9*V**2+
     4   2*V-2)*W-(V-1)**4*(3*V**2-2*V+2))-(V*W-1)*(V*W-V+1)*(4*V
     5   **8*W**7-4*V**6*(4*V**2-3*V+3)*W**6+V**6*(28*V**2-31*V+3
     6   1)*W**5-V**4*(24*V**4-15*V**3-V**2+32*V-16)*W**4+V**4*(8
     7   *V**4+20*V**3-41*V**2+42*V-21)*W**3-2*(V-1)*V**2*(8*V**4
     8   +5*V**3-16*V**2+22*V-11)*W**2+4*(V-1)**2*V**2*(4*V**2-3*
     9   V+3)*W-2*(V-1)**3*(2*V**2-V+1)))-(V-1)**2*V**2*VC*W*(V**
     :   2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(4*V**8
     ;   *W**8+CQ*(2*V**8*W**7-2*V**6*(V+2)*(3*V-2)*W**6+V**6*(7*
     <   V**2+27*V-27)*W**5-V**4*(4*V**4+37*V**3-43*V**2+12*V-6)*
     =   W**4+V**4*(V**4+25*V**3-27*V**2+4*V-2)*W**3-(V-1)*V**4*(
     >   7*V**2+10*V-10)*W**2+(V-1)**2*V**2*(7*V**2+V-1)*W-(V-1)*
     ?   *3*(3*V**2-2*V+2))-18*V**8*W**7+4*V**6*(8*V**2+5*V-5)*W*
     @   *6-V**6*(29*V**2+63*V-63)*W**5+V**4*(14*V**4+79*V**3-73*
     1   V**2-12*V+6)*W**4-V**4*(3*V**4+49*V**3-27*V**2-44*V+22)*
     2   W**3+(V-1)*V**2*(V**2+2*V-2)*(13*V**2+2*V-2)*W**2-(V-1)*
     3   *2*V**2*(13*V**2+5*V-5)*W+(V-1)**3*(5*V**2-2*V+2)))/(N**
     4   2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
 
      LV1 = LOG(1-V)*(N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**4*W*
     1   *4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W
     2   +(V-1)**4)*(CQ*(V-1)*(V*W-V+1)*(4*V**8*W**7-8*V**7*(V+1)
     3   *W**6+4*V**6*(2*V**2+3*V+2)*W**5-4*V**3*(4*V**4+3*V**3+5
     4   *V**2-8*V+4)*W**4+V**2*(16*V**4+4*V**3+69*V**2-112*V+48)
     5   *W**3-2*V*(4*V**4-V**3+50*V**2-72*V+24)*W**2+(2*V**4-2*V
     6   **3+65*V**2-80*V+16)*W-16*(V-1))+2*(V**2*W**2-2*V*W+1)*(
     7   4*V**8*W**6-16*V**6*(V**2-V+1)*W**5+V**5*(28*V**3-47*V**
     8   2+48*V+7)*W**4-2*V**2*(12*V**6-26*V**5+36*V**4-3*V**3-4*
     9   V**2+13*V-4)*W**3+V**3*(8*V**5-11*V**4+29*V**3+7*V+7)*W*
     :   *2-(8*V**7-3*V**6+7*V**5+17*V**4-37*V**3+50*V**2-34*V+8)
     ;   *W+8*(V-1)*(V+1)*(V**2-V+1)**2))-2*N**2*(V-1)*V*VC*W*(V*
     <   *2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2
     =   *V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(V-1)*V*(V*W-V+1
     >   )*(2*V**6*W**6-4*V**5*(V+1)*W**5+V**4*(4*V**2+6*V+5)*W**
     ?   4-V**3*(8*V**2+8*V+3)*W**3+2*V**2*(5*V**2+V+2)*W**2-V*(6
     @   *V**2-2*V+3)*W+2*V**2-2*V+1)+(V**3*W**3-3*V**2*W**2+3*V*
     1   W-1)*(4*V**6*W**4-4*V**4*(3*V**2-3*V+2)*W**3+V**3*(16*V*
     2   *3-33*V**2+37*V-12)*W**2-(V-1)*V*(8*V**4-14*V**3+16*V**2
     3   -V+1)*W+(V-1)**4*(V+1)))+(V-1)**2*V**2*VC*W*(V**2*W**2-2
     4   *V*W+1)*(2*(V**2*W-(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*
     5   W-1)+CQ*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(2*V**2*W**2-2
     6   *V*W+1))*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W
     7   **3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))/(N**
     8   2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
 
      LV = LOG(V)*(-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2
     1   -2*(V-1)*V*W+(V-1)**2)*(CQ*(4*V**12*W**10-4*V**10*(5*V**
     2   2-2*V+2)*W**9+4*V**10*(11*V**2-7*V+7)*W**8-4*V**6*(13*V*
     3   *6-7*V**5+7*V**4+8*V**3-24*V**2+24*V-8)*W**7+V**6*(32*V*
     4   *6+13*V**5+8*V**4+22*V**3-171*V**2+192*V-64)*W**6-V**4*(
     5   8*V**8+38*V**7+35*V**6-160*V**5+211*V**4-426*V**3+590*V*
     6   *2-384*V+96)*W**5+(V-1)*V**4*(17*V**6+87*V**5-107*V**4+1
     7   36*V**3-308*V**2+288*V-96)*W**4-(V-1)**2*V**2*(35*V**6+8
     8   4*V**5-188*V**4+304*V**3-392*V**2+288*V-96)*W**3+(V-1)**
     9   3*V**4*(59*V**2-45*V+45)*W**2-(V-1)**4*(51*V**4-98*V**3+
     :   130*V**2-64*V+32)*W+16*(V-1)**5*(V**2-2*V+2))-8*V**12*W*
     ;   *10+16*V**10*(3*V**2-2*V+2)*W**9-V**6*(144*V**6-199*V**5
     <   +219*V**4-44*V**3+32*V**2-12*V+4)*W**8+4*V**6*(68*V**6-1
     =   31*V**5+135*V**4+2*V**3-26*V**2+30*V-10)*W**7-V**4*(344*
     >   V**8-741*V**7+573*V**6+460*V**5-552*V**4+420*V**3-196*V*
     ?   *2+48*V-12)*W**6+2*V**4*(144*V**8-265*V**7-62*V**6+729*V
     @   **5-486*V**4-39*V**3+321*V**2-264*V+66)*W**5-V**2*(144*V
     1   **10-34*V**9-1130*V**8+2671*V**7-2073*V**6+537*V**5+437*
     2   V**4-600*V**3+240*V**2-60*V+12)*W**4+2*V**2*(16*V**10+10
     3   4*V**9-569*V**8+1019*V**7-730*V**6+331*V**5-437*V**4+712
     4   *V**3-718*V**2+360*V-72)*W**3-2*(V-1)*(48*V**10-80*V**9-
     5   165*V**8+730*V**7-997*V**6+846*V**5-422*V**4+108*V**3-12
     6   *V**2-10*V+2)*W**2+2*(V-1)**2*(48*V**8-184*V**7+338*V**6
     7   -351*V**5+309*V**4-233*V**3+199*V**2-104*V+26)*W-16*(V-1
     8   )**3*(V**2-V+1)**2*(2*V**2-5*V+5))+2*N**2*(V-1)*VC*W*(V*
     9   *2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(CQ*V*
     :   *2*(2*V**10*W**9-2*V**8*(5*V**2-2*V+2)*W**8+V**8*(22*V**
     ;   2-13*V+13)*W**7-V**6*(26*V**4-11*V**3+15*V**2-8*V+4)*W**
     <   6+2*V**6*(8*V**4+5*V**3+3*V**2-16*V+8)*W**5-V**4*(4*V**6
     =   +21*V**5+10*V**4-72*V**3+61*V**2-30*V+10)*W**4+(V-1)*V**
     >   4*(9*V**4+32*V**3-35*V**2+6*V-3)*W**3-(V-1)**2*V**2*(13*
     ?   V**4+15*V**3-19*V**2+8*V-4)*W**2+(V-1)**3*V**2*(9*V**2+2
     @   *V-2)*W-(V-1)**4*(3*V**2-2*V+2))-(V*W-1)*(V*W-V+1)*(4*V*
     1   *10*W**7-4*V**8*(4*V**2-3*V+3)*W**6+V**4*(28*V**6-30*V**
     2   5+27*V**4+8*V**3-9*V**2+6*V-2)*W**5-V**4*(24*V**6-15*V**
     3   5-2*V**4+40*V**3-35*V**2+18*V-6)*W**4+V**2*(8*V**8+19*V*
     4   *7-35*V**6+32*V**5-12*V**4-16*V**3+24*V**2-16*V+4)*W**3-
     5   (V-1)*V**2*(16*V**6+13*V**5-37*V**4+56*V**3-48*V**2+24*V
     6   -8)*W**2+(V-1)**2*(18*V**6-13*V**5+14*V**4-5*V**2+6*V-2)
     7   *W-(V-1)**3*(6*V**4-5*V**3+7*V**2-4*V+2)))-(V-1)**2*V**2
     8   *VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**
     9   2)*(8*V**8*W**8+CQ*(2*V**8*W**7-2*V**6*(V+2)*(3*V-2)*W**
     :   6+V**6*(7*V**2+27*V-27)*W**5-V**4*(4*V**4+37*V**3-43*V**
     ;   2+12*V-6)*W**4+V**4*(V**4+25*V**3-27*V**2+4*V-2)*W**3-(V
     <   -1)*V**4*(7*V**2+10*V-10)*W**2+(V-1)**2*V**2*(7*V**2+V-1
     =   )*W-(V-1)**3*(3*V**2-2*V+2))-V**6*(35*V**2-2*V+2)*W**7+2
     >   *V**6*(31*V**2+12*V-12)*W**6-3*V**4*(19*V**4+30*V**3-32*
     ?   V**2+4*V-2)*W**5+2*V**6*(14*V**2+61*V-61)*W**4-V**2*(6*V
     @   **6+80*V**5-49*V**4-68*V**3+49*V**2-18*V+6)*W**3+2*(V-1)
     1   *V**2*(11*V**4+23*V**3-25*V**2+4*V-2)*W**2-2*(V-1)**2*(1
     2   1*V**4+2*V**3-3*V**2+2*V-1)*W+4*(V-1)**3*(2*V**2-V+1)))/
     3   (N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
 
      LVW = (-2*N**4*(V-1)*VC*(4*V**6*W**4-8*V**4*(V**2-V+2)*W**3
     1   +V**3*(8*V**3-17*V**2+26*V+15)*W**2+(6*V**5-24*V**4-11*V
     2   **3+33*V**2-52*V+16)*W+16*(V**2-V+1)**2)*(V**5*W**5-5*V*
     3   *4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*
     4   (V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W
     5   **2+5*(V-1)**4*V*W-(V-1)**5)+2*N**2*(V-1)*V*VC*W*(4*V**5
     6   *W**3-2*V**3*(4*V**2-5*V+5)*W**2+V**3*(8*V**2-15*V+15)*W
     7   -(V-1)*(2*V**3-V**2-V+1))*(V**5*W**5-5*V**4*W**4+10*V**3
     8   *W**3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4
     9   +10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*
     :   V*W-(V-1)**5)+2*(V-1)**2*V**2*VC*W*(2*V**2*W**2-3*V**2*W
     ;   +2*V**2-2*V+1)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V*
     <   *2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**
     =   2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**
     >   5))*LOG(1-V*W)/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W
     ?   -V+1)**5)
 
      LTVW = (-2*N**4*(V-1)*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(
     1   V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)
     2   *(CQ*(4*V**8*W**7-8*V**7*(2*V-1)*W**6+4*V**6*(7*V**2-7*V
     3   +2)*W**5-4*(V-1)*V**3*(8*V**4-5*V**3+5*V**2-8*V+4)*W**4+
     4   (V-1)**2*V**2*(25*V**4+2*V**3+21*V**2-80*V+48)*W**3-2*(V
     5   -1)**3*V*(5*V**4+21*V**3-22*V**2-24*V+24)*W**2+(V-1)**4*
     6   (V**4+48*V**3-79*V**2+16*V+16)*W-16*(V-1)**7)-(V-1)*V*W*
     7   (V**5*W**4-V**3*(3*V**3+3*V**2-4)*W**3+(V-1)*V**2*(2*V**
     8   3+5*V**2+15*V-12)*W**2+(V-1)**2*V*(5*V**3-3*V**2-16*V+12
     9   )*W-4*(V-1)**3*(V+1)*(V**2-V+1)))+2*N**2*(V-1)*V*VC*W*(V
     :   **2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**5*W**5-5*V**4*W**4+10
     ;   *V**3*W**3-10*V**2*W**2+5*V*W-1)*(2*CQ*V*(2*V**6*W**6-4*
     <   V**5*(2*V-1)*W**5+V**4*(15*V**2-16*V+5)*W**4-(V-1)*V**3*
     =   (19*V**2-14*V+3)*W**3+2*(V-1)**2*V**2*(8*V**2-5*V+2)*W**
     >   2-(V-1)**3*V*(7*V**2-4*V+3)*W+(V-1)**4*(V**2+1))+(V-1)*(
     ?   V*W-V+1)*(2*V**5*W**4-V**4*(5*V-6)*W**3+V**2*(11*V**3-14
     @   *V**2-4*V+1)*W**2-(V-1)*V*(9*V**3-7*V**2+2)*W+(V-1)**2*(
     1   V+1)*(V**2-V+1)))+2*(V-1)**2*V**2*VC*W*(V**2*W**2-2*(V-1
     2   )*V*W+(V-1)**2)*(V*W*(2*V**4*W**4-V**3*(7*V-4)*W**3+V**2
     3   *(10*V**2-13*V+5)*W**2-(V-1)*V*(7*V**2-7*V+2)*W+(V-1)**2
     4   *(2*V**2-V+1))-CQ*(V-1)*(V**2*W**2-2*V**2*W+V**2+1)*(2*V
     5   **2*W**2-2*(V-1)*V*W+(V-1)**2))*(V**5*W**5-5*V**4*W**4+1
     6   0*V**3*W**3-10*V**2*W**2+5*V*W-1))*LOG(V*W-V+1)/(N**2*(
     7   V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
 
      LW = (-(V-1)**2*V**2*VC*W*(4*V**2*W**2-(9*V**2-2*V+2)*W+5*V
     1   **2-6*V+6)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W
     2   **2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V*
     3   *3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)-N
     4   **4*(V-1)*VC*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)
     5   *(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-
     6   1)**3*V*W+(V-1)**4)*(4*V**6*(2*V**2-3*V+3)*W**5-V**2*(40
     7   *V**6-89*V**5+97*V**4-12*V**3-4*V**2+12*V-4)*W**4+2*V**2
     8   *(40*V**6-104*V**5+121*V**4-14*V**3-43*V**2+60*V-20)*W**
     9   3-8*CQ*(V**2-V+1)**3*(V*W-1)*(V*W-V+1)*(W**2-2*W+2)-(80*
     :   V**8-223*V**7+308*V**6-174*V**5+101*V**4-28*V**3+28*V**2
     ;   -16*V+4)*W**2+(32*V**8-64*V**7+33*V**6+92*V**5-85*V**4-5
     <   4*V**3+186*V**2-144*V+36)*W-16*(V-1)*(V**2-V+1)**2*(2*V*
     =   *2-3*V+3))-2*N**2*(V-1)**2*VC*W*(V**2*(V**2-2*V+2)*(V**2
     >   -V+1)*W**3-(V-1)*V**2*(V**2+4*V-4)*W**2+(V**2+2*V-2)*(V*
     ?   *4-V**3+2*V**2-2*V+1)*W-(V-1)*(3*V**4-4*V**3+6*V**2-4*V+
     @   2))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W*
     1   *4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W
     2   +(V-1)**4))*LOG(W)/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*
     3   (V*W-V+1)**5)
 
      CVC = (3*N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-
     1   1)*V*W+(V-1)**2)*(8*V**12*W**10-8*V**10*(3*V**2+2*V-2)*W
     2   **9+V**6*(24*V**6+91*V**5-141*V**4+148*V**3-194*V**2+144
     3   *V-48)*W**8+2*V**6*(12*V**6-125*V**5+230*V**4-306*V**3+3
     4   93*V**2-288*V+96)*W**7-V**4*(128*V**8-494*V**7+927*V**6-
     5   1036*V**5+799*V**4+66*V**3-694*V**2+576*V-144)*W**6+V**4
     6   *(192*V**8-522*V**7+773*V**6-186*V**5-1129*V**4+2676*V**
     7   3-2908*V**2+1728*V-432)*W**5-V**2*(128*V**10-107*V**9-42
     8   2*V**8+1971*V**7-3646*V**6+4107*V**5-2461*V**4+72*V**3+1
     9   062*V**2-720*V+144)*W**4+V**2*(32*V**10+192*V**9-891*V**
     :   8+1964*V**7-2315*V**6+1082*V**5+1366*V**4-3208*V**3+2962
     ;   *V**2-1440*V+288)*W**3-(V-1)*(96*V**10-96*V**9-217*V**8+
     <   1119*V**7-1982*V**6+2191*V**5-1393*V**4+280*V**3+290*V**
     =   2-240*V+48)*W**2+(V-1)**2*(96*V**8-256*V**7+427*V**6-302
     >   *V**5+3*V**4+312*V**3-328*V**2+192*V-48)*W-32*(V-1)**3*(
     ?   V**2-V+1)**3)-6*AL*N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V*
     @   *5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)
     1   **3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(4*(V-1)*V**8*W**
     2   7-8*(V-1)*V**7*(V+1)*W**6+4*V**3*(6*V**5-11*V**4+11*V**3
     3   -12*V**2+6*V-2)*W**5+4*V**2*(4*V**7-10*V**6+10*V**5+5*V*
     4   *4-15*V**3+24*V**2-14*V+6)*W**4-V*(16*V**8-40*V**6+128*V
     5   **5-113*V**4+81*V**3+16*V**2-24*V+24)*W**3+2*(24*V**8-48
     6   *V**7+72*V**6-29*V**5-5*V**4+46*V**3-24*V**2+12*V+4)*W**
     7   2-(48*V**7-128*V**6+238*V**5-236*V**4+173*V**3-47*V**2+1
     8   6)*W+16*(V**2-V+1)**3)-6*N**2*(V-1)*V**2*VC*W*(V**3*W**3
     9   -3*V**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-
     :   1)**2*V*W-(V-1)**3)*(4*V**8*W**7-8*V**6*(V**2+V-1)*W**6+
     ;   V**4*(4*V**4+39*V**3-60*V**2+42*V-21)*W**5-3*(V-1)*V**4*
     <   (15*V**2-11*V+11)*W**4+3*(V-1)*V**2*(7*V**4-6*V**3-8*V**
     =   2+28*V-14)*W**3-(V-1)*V**2*(7*V**4-18*V**3+12*V**2+12*V-
     >   6)*W**2+(V-1)**3*(10*V**2-21*V+21)*W-(V-1)**3*(3*V**2-5*
     ?   V+5))+3*(V-1)**2*V**2*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**
     @   2-2*(V-1)*V*W+(V-1)**2)*(V**6*(11*V**2-4*V+4)*W**7-2*V**
     1   6*(15*V**2-11*V+11)*W**6+V**4*(30*V**4-23*V**3+11*V**2+2
     2   4*V-12)*W**5-V**4*(14*V**4-23*V**3+89*V**2-132*V+66)*W**
     3   4+V**2*(3*V**6-23*V**5+144*V**4-254*V**3+157*V**2-36*V+1
     4   2)*W**3+(V-1)*V**2*(5*V**4-44*V**3+54*V**2-20*V+10)*W**2
     5   +(V-1)**2*(5*V**4-5*V**3+V**2+8*V-4)*W-(V-1)**3*(V**2-2*
     6   V+2))+12*AL*N**2*(V-1)**2*V**2*VC*W*(V**2*W**2-2*V*W+1)*
     7   (V**2*W**2-V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(2*V
     8   **2*W**2-2*V*W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)*
     9   *2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)*
     :   *5)-6*AL*(V-1)**2*V**2*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W*
     ;   *2-2*V**2*W+2*V**2-2*V+1)*(2*V**2*W**2-2*V*W+1)*(V**5*W*
     <   *5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V
     =   **2*W**2+5*(V-1)**4*V*W-(V-1)**5)+4*CQ*N**4*(V-1)*VC*W*(
     >   V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(12
     ?   *V**6*W**4-24*(V-1)*V**5*W**3+8*V**2*(2*V**6-6*V**5+15*V
     @   **4-20*V**3+15*V**2-6*V+2)*W**2-4*(V-1)*V*(8*V**6-24*V**
     1   5+51*V**4-62*V**3+51*V**2-24*V+8)*W+(V-1)**2*(16*V**6-48
     2   *V**5+99*V**4-118*V**3+99*V**2-48*V+16))*(V**5*W**5-5*V*
     3   *4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)-24*CQ*N**2*(V
     4   -1)*V**2*VC*W*(V**2*W**2-(V-1)*V*W+(V-1)**2)*(2*V**2*W**
     5   2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(
     6   V-1)**2*V*W-(V-1)**3)*(V**5*W**5-5*V**4*W**4+10*V**3*W**
     7   3-10*V**2*W**2+5*V*W-1)+12*CQ*(V-1)**3*V**2*VC*W*(2*V**2
     8   *W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*(V-1)*V**2*W**2
     9   +3*(V-1)**2*V*W-(V-1)**3)*(V**5*W**5-5*V**4*W**4+10*V**3
     :   *W**3-10*V**2*W**2+5*V*W-1))/(N**2*(V-1)**3*V**3*W**2*(V
     ;   *W-1)**5*(V*W-V+1)**5)/6.D0
 
      LM = LOG(S/M**2)*(N**4*VC*(V**3*W**3-3*(V-1)*V**2*W**2+3*(
     1   V-1)**2*V*W-(V-1)**3)*(4*V**9*W**7-4*V**7*(2*V**2+V+4)*W
     2   **6+V**3*(16*V**6-15*V**5+55*V**4+20*V**3+32*V**2-12*V+4
     3   )*W**5-V**2*(16*V**7-3*V**6+26*V**5+77*V**4+24*V**3+108*
     4   V**2-40*V+12)*W**4+V*(16*V**8+11*V**6+80*V**5+54*V**4+23
     5   *V**3+164*V**2-48*V+12)*W**3-(48*V**8-96*V**7+195*V**6-1
     6   28*V**5+213*V**4-96*V**3+164*V**2-24*V+4)*W**2+(48*V**7-
     7   128*V**6+275*V**5-297*V**4+319*V**3-177*V**2+108*V-4)*W-
     8   16*(V**2-V+1)**2*(V**2-V+2))-2*N**2*VC*W*(V*W-1)*(V**3*W
     9   **3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(2*V**8*W
     :   **5-4*V**6*(V**2+1)*W**4+V**2*(4*V**6+V**5-V**4+15*V**3-
     ;   9*V**2+6*V-2)*W**3-V*(V**2-V+1)*(7*V**4+5*V**3-2*V**2+6*
     <   V-4)*W**2+(10*V**6-19*V**5+22*V**4-12*V**3+3*V**2+2*V-2)
     =   *W-(V-1)*(3*V**4-5*V**3+6*V**2-4*V+2))+(V-1)*V**2*VC*W*(
     >   V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V*
     ?   *3*(V**2+2)*W**4-V**2*(V**3+5*V**2-2*V+6)*W**3+V*(7*V**3
     @   -3*V**2+3*V+6)*W**2-(7*V**3-7*V**2+6*V+2)*W+3*V**2-4*V+3
     1   ))/(N**2*(V-1)**2*V**3*W**2*(V*W-1)**3*(V*W-V+1)**3)
 
      LMP = LOG(S/MP**2)*(N**4*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1
     1   )*(4*V**9*W**7-4*V**7*(7*V**2-9*V+4)*W**6+4*V**6*(23*V**
     2   3-56*V**2+51*V-16)*W**5-4*(V-1)*V**3*(45*V**5-110*V**4+1
     3   01*V**3-27*V**2-8*V+4)*W**4+(V-1)**2*V**2*(224*V**5-535*
     4   V**4+498*V**3-123*V**2-80*V+48)*W**3-2*(V-1)**3*V*(88*V*
     5   *5-203*V**4+205*V**3-70*V**2-24*V+24)*W**2+(V-1)**4*(80*
     6   V**5-175*V**4+208*V**3-111*V**2+16*V+16)*W-16*(V-1)**5*(
     7   V**2-V+1)**2)-2*N**2*V**2*VC*W*(V*W-V+1)*(V**3*W**3-3*V*
     8   *2*W**2+3*V*W-1)*(2*V**6*W**5-4*V**4*(2*V**2-2*V+1)*W**4
     9   +V**3*(2*V-1)*(7*V**2-10*V+7)*W**3-6*(V-1)*V**2*(2*V-1)*
     :   (V**2-V+1)*W**2+2*(V-1)**2*V*(2*V**3-V**2+2*V+1)*W-(V-1)
     ;   **3*(V**2+1))+(V-1)**2*V**2*VC*W*(V**2*W**2-2*V**2*W+V**
     <   2+1)*(2*V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**
     =   2*W**2+3*V*W-1))/(N**2*(V-1)**2*V**3*W**2*(V*W-1)**3*(V*
     >   W-V+1)**3)
 
      STRUV14=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV15(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(-16*N**5*(V-1)*VC*(V*W-1)*(V*W-V+1)*(4*V**
     1   14*W**13-2*V**12*(10*V**2+V-1)*W**12+V**8*(52*V**6-V**5+
     2   7*V**4-14*V**3+12*V**2-6*V+2)*W**11-V**8*(92*V**6-36*V**
     3   5+71*V**4-80*V**3+65*V**2-30*V+10)*W**10+V**6*(110*V**8-
     4   59*V**7+75*V**6-22*V**5-22*V**4+62*V**3-58*V**2+32*V-8)*
     5   W**9-V**6*(90*V**8-14*V**7-115*V**6+350*V**5-437*V**4+40
     6   4*V**3-284*V**2+128*V-32)*W**8+V**4*(62*V**10+V**9-166*V
     7   **8+340*V**7-183*V**6-30*V**5+122*V**4-168*V**3+132*V**2
     8   -60*V+12)*W**7-2*CQ*(V**2*W-V+1)*(V**2*W**2-2*V*W+1)*(V*
     9   *2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**6*(V**2-2*V+2)*W**6-V*
     :   *4*(V**4-2*V**3+4*V**2-4*V+2)*W**5-2*V**4*(V**4-4*V**3+5
     ;   *V**2-2*V+1)*W**4+2*V**2*(V**6-4*V**5+7*V**4-8*V**3+9*V*
     <   *2-6*V+2)*W**3-2*(V-1)*V**2*(V**4-4*V**3+5*V**2-2*V+1)*W
     =   **2-(V-1)**2*(V**4-2*V**3+4*V**2-4*V+2)*W+(V-1)**3*(V**2
     >   -2*V+2))-V**4*(42*V**10-14*V**9-88*V**8+45*V**7+481*V**6
     ?   -937*V**5+975*V**4-784*V**3+466*V**2-180*V+36)*W**6+V**2
     @   *(20*V**12+31*V**11-214*V**10+344*V**9-127*V**8-2*V**7-1
     1   66*V**6+328*V**5-370*V**4+280*V**3-144*V**2+48*V-8)*W**5
     2   -V**2*(4*V**12+50*V**11-190*V**10+241*V**9-22*V**8-71*V*
     3   *7-221*V**6+592*V**5-739*V**4+570*V**3-290*V**2+96*V-16)
     4   *W**4+(V-1)*(16*V**12+14*V**11-206*V**10+535*V**9-696*V*
     5   *8+679*V**7-569*V**6+436*V**5-301*V**4+150*V**3-52*V**2+
     6   12*V-2)*W**3-(V-1)**2*(24*V**10-58*V**9+32*V**8+109*V**7
     7   -202*V**6+189*V**5-77*V**4+15*V**2-10*V+2)*W**2+(V-1)**3
     8   *(16*V**8-58*V**7+113*V**6-128*V**5+115*V**4-78*V**3+54*
     9   V**2-24*V+6)*W-4*(V-1)**4*(V**2-2*V+2)*(V**2-V+1)**2)-8*
     :   GTR*N**2*(V-1)**2*VC*W*(V*W-1)*(V*W-V+1)*(CQ*(2*V**8*(2*
     ;   V**4-3*V**3+5*V**2-4*V+2)*W**10-V**8*(19*V**4-30*V**3+50
     <   *V**2-40*V+20)*W**9+V**6*(42*V**6-65*V**5+83*V**4-20*V**
     =   3-30*V**2+48*V-16)*W**8-2*V**6*(27*V**6-34*V**5+10*V**4+
     >   80*V**3-120*V**2+96*V-32)*W**7+2*V**4*(21*V**8-13*V**7-3
     ?   3*V**6+126*V**5-136*V**4+54*V**3+38*V**2-48*V+12)*W**6-V
     @   **4*(19*V**8+16*V**7-66*V**6+80*V**5+82*V**4-348*V**3+45
     1   2*V**2-288*V+72)*W**5+V**2*(V**2+4*V-4)*(4*V**8+9*V**7-3
     2   7*V**6+70*V**5-66*V**4+26*V**3+10*V**2-16*V+4)*W**4-2*(V
     3   -1)*V**2*(5*V**8+15*V**7-49*V**6+68*V**5-18*V**4-64*V**3
     4   +96*V**2-64*V+16)*W**3+2*(V-1)**2*(8*V**8-4*V**7-2*V**6+
     5   21*V**5-31*V**4+19*V**3+3*V**2-8*V+2)*W**2-(V-1)**3*(10*
     6   V**6-5*V**5+7*V**4-10*V**2+12*V-4)*W+(V-1)**4*(4*V**4-5*
     7   V**3+7*V**2-4*V+2))-2*(V**8*(3*V**4-3*V**3+5*V**2-4*V+2)
     8   *W**10-2*V**8*(8*V**4-7*V**3+12*V**2-10*V+5)*W**9+V**6*(
     9   38*V**6-24*V**5+33*V**4-10*V**3-15*V**2+24*V-8)*W**8-V**
     :   6*(50*V**6-3*V**5-19*V**4+76*V**3-118*V**2+96*V-32)*W**7
     ;   +V**4*(38*V**8+43*V**7-75*V**6+98*V**5-122*V**4+54*V**3+
     <   38*V**2-48*V+12)*W**6-V**4*(16*V**8+60*V**7-43*V**6-44*V
     =   **5+83*V**4-174*V**3+226*V**2-144*V+36)*W**5+V**2*(3*V**
     >   10+36*V**9+24*V**8-162*V**7+216*V**6-238*V**5+182*V**4-4
     ?   0*V**3-50*V**2+40*V-8)*W**4-(V-1)*V**2*(9*V**8+47*V**7-5
     @   8*V**6+20*V**5+11*V**4-70*V**3+98*V**2-64*V+16)*W**3+(V-
     1   1)**2*(14*V**8+16*V**7-21*V**6+19*V**5-30*V**4+19*V**3+3
     2   *V**2-8*V+2)*W**2-(V-1)**3*(9*V**6+3*V**5-3*V**4+2*V**3-
     3   6*V**2+6*V-2)*W+(V-1)**4*(3*V**4-2*V**3+3*V**2-2*V+1)))+
     4   8*(CQ-1)*GTR*N**4*(V-1)**2*VC*W*(V*W-1)*(V*W-V+1)*(2*V**
     5   8*(V**4-2*V**3+4*V**2-4*V+2)*W**10-V**8*(9*V**4-20*V**3+
     6   40*V**2-40*V+20)*W**9+V**6*(20*V**6-47*V**5+73*V**4-36*V
     7   **3-22*V**2+48*V-16)*W**8-2*V**6*(13*V**6-28*V**5+20*V**
     8   4+48*V**3-104*V**2+96*V-32)*W**7+2*V**4*(10*V**8-17*V**7
     9   -5*V**6+84*V**5-130*V**4+72*V**3+32*V**2-48*V+12)*W**6-V
     :   **4*(9*V**8-8*V**7-10*V**6+52*V**5+6*V**4-240*V**3+416*V
     ;   **2-288*V+72)*W**5+V**2*(2*V**10+5*V**9+3*V**8-78*V**7+2
     <   62*V**6-442*V**5+390*V**4-112*V**3-92*V**2+80*V-16)*W**4
     =   -2*(V-1)*V**2*(2*V**8+7*V**7-35*V**6+64*V**5-36*V**4-40*
     >   V**3+88*V**2-64*V+16)*W**3+2*(V-1)**2*(4*V**8-7*V**7+6*V
     ?   **6+12*V**5-29*V**4+22*V**3+2*V**2-8*V+2)*W**2-(V-1)**3*
     @   (4*V**6-3*V**5+7*V**4-4*V**3-8*V**2+12*V-4)*W+(V-1)**4*(
     1   2*V**4-3*V**3+5*V**2-4*V+2))+8*GTR*(V-1)**2*V**2*VC*W*(V
     2   **3*W**3-3*V**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W*
     3   *2+3*(V-1)**2*V*W-(V-1)**3)*(2*CQ*(V**4*(V**2-V+1)*W**6-
     4   3*V**4*(V**2-V+1)*W**5+2*V**2*(2*V**4-2*V**3+V**2+2*V-1)
     5   *W**4-V**2*(3*V**4-3*V**3-V**2+8*V-4)*W**3+(V**6-V**4+V*
     6   *3+2*V**2-3*V+1)*W**2-(V-1)*(V**4+V**3-2*V**2+2*V-1)*W+(
     7   V-1)**2*(V**2-V+1))-V**4*(3*V**2-2*V+2)*W**6+2*V**4*(5*V
     8   **2-4*V+4)*W**5-2*V**2*(V**2-V+1)*(7*V**2+2*V-2)*W**4+2*
     9   V**2*(5*V**4-2*V**3-2*V**2+8*V-4)*W**3-(V**2+V-1)**2*(3*
     :   V**2-2*V+2)*W**2+2*(V-1)*V**2*(V+1)*(2*V-1)*W-(V-1)**2*(
     ;   3*V**2-2*V+2)))/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W
     <   -V+1)**5)
 
 
      LV1 = LOG(1-V)*(-16*N**5*(V-1)*VC*(V**3*W**3-3*V**2*W**2+3
     1   *V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**
     2   2-4*(V-1)**3*V*W+(V-1)**4)*((V*W-1)*(4*V**8*W**7-10*V**8
     3   *W**6+2*V**5*(8*V**3-5*V**2+6*V+1)*W**5-V**5*(14*V**3-15
     4   *V**2+18*V+5)*W**4+V**2*(6*V**6-6*V**5+6*V**4+15*V**3-8*
     5   V**2+5*V-2)*W**3-V**3*(2*V**5-3*V**4+9*V**3-2*V**2+5*V-1
     6   )*W**2+(2*V-1)*(V**2-V+1)*(V**4+V**3+V**2+3*V-2)*W-2*(V-
     7   1)*(V+1)*(V**2-V+1)**2)+2*CQ*((V-1)*W-1)*(V*W-V+1)*(V*W+
     8   V-1)*((V-1)*V*W+1)*(V**2*W-1)*(V**2*W-V+1))-4*GTR*N**4*(
     9   V-1)**2*V*VC*W*(V*W-1)*(CQ*V*(V**2*W**2-2*V**2*W+2*V**2-
     :   2*V+1)*(V**4*W**4+2*V**2*W**2+1)-2*(V**2*W-(V-1)**2)*(V*
     ;   *4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1))*(V**5*W**5-5*(
     <   V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W*
     =   *2+5*(V-1)**4*V*W-(V-1)**5)+8*GTR*N**2*(V-1)**2*V*VC*W*(
     >   V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2
     ?   -4*(V-1)**3*V*W+(V-1)**4)*((V**3*W**3-3*V**2*W**2+3*V*W-
     @   1)*(2*V**5*W**4-V**4*(3*V+2)*W**3+V**2*(5*V**3-3*V**2+5*
     1   V-2)*W**2-2*V**5*W+(V-1)*(V+1)*(2*V**2-3*V+2))+CQ*V*(V*W
     2   -V+1)*(V**2*W**2+1)*(V**2*W**2-V*W+1)*(V**2*W**2-2*V**2*
     3   W+2*V**2-2*V+1))-4*GTR*(V-1)**2*V**2*VC*W*(V**3*W**3-3*V
     4   **2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**
     5   2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*(V*W-1)*(2*V**3*
     6   W**4-4*V**3*W**3+V*(V**3+2*V**2-V+2)*W**2-V*(V**3-V+2)*W
     7   +(V-1)**3)+CQ*(V*W-V+1)*(V**2*W**2+1)*(V**2*W**2-2*V**2*
     8   W+2*V**2-2*V+1)))/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V
     9   *W-V+1)**5)
 
 
      LV = LOG(V)*(-16*N**5*(V-1)*VC*(V*W-1)*(V*W-V+1)*(4*V**14*
     1   W**13-22*V**14*W**12+2*V**8*(31*V**6-5*V**5+8*V**4-7*V**
     2   3+6*V**2-3*V+1)*W**11-V**8*(116*V**6-53*V**5+83*V**4-70*
     3   V**3+60*V**2-30*V+10)*W**10+2*V**6*(74*V**8-41*V**7+40*V
     4   **6+7*V**5-20*V**4+31*V**3-29*V**2+16*V-4)*W**9-V**6*(13
     5   2*V**8-34*V**7-137*V**6+432*V**5-473*V**4+398*V**3-282*V
     6   **2+128*V-32)*W**8+V**4*(92*V**10+8*V**9-261*V**8+514*V*
     7   *7-265*V**6-36*V**5+124*V**4-168*V**3+132*V**2-60*V+12)*
     8   W**7-2*CQ*(V**2*W-V+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*
     9   (V-1)*V*W+(V-1)**2)*(V**6*(V**2-2*V+2)*W**6-V**4*(V**4-2
     :   *V**3+4*V**2-4*V+2)*W**5-2*V**4*(V**4-4*V**3+5*V**2-2*V+
     ;   1)*W**4+2*V**2*(V**6-4*V**5+7*V**4-8*V**3+9*V**2-6*V+2)*
     <   W**3-2*(V-1)*V**2*(V**4-4*V**3+5*V**2-2*V+1)*W**2-(V-1)*
     =   *2*(V**4-2*V**3+4*V**2-4*V+2)*W+(V-1)**3*(V**2-2*V+2))-V
     >   **4*(54*V**10+19*V**9-228*V**8+241*V**7+430*V**6-999*V**
     ?   5+1005*V**4-792*V**3+468*V**2-180*V+36)*W**6+V**2*(22*V*
     @   *12+56*V**11-293*V**10+388*V**9+23*V**8-242*V**7-30*V**6
     1   +280*V**5-358*V**4+280*V**3-144*V**2+48*V-8)*W**5-V**2*(
     2   4*V**12+56*V**11-198*V**10+171*V**9+214*V**8-359*V**7-41
     3   *V**6+508*V**5-703*V**4+560*V**3-288*V**2+96*V-16)*W**4+
     4   (V-1)*(16*V**12+18*V**11-242*V**10+597*V**9-719*V**8+661
     5   *V**7-549*V**6+424*V**5-298*V**4+150*V**3-52*V**2+12*V-2
     6   )*W**3-2*(V-1)**2*(12*V**10-31*V**9+11*V**8+76*V**7-133*
     7   V**6+127*V**5-61*V**4+10*V**3+5*V**2-5*V+1)*W**2+2*(V-1)
     8   **3*(V**2-V+1)*(8*V**6-24*V**5+31*V**4-17*V**3+16*V**2-9
     9   *V+3)*W-2*(V-1)**4*(V**2-V+1)**2*(2*V**2-5*V+5))+4*GTR*N
     :   **4*(V-1)**2*VC*W*(V*W-1)*(V*W-V+1)*(CQ*V**2*(2*V**10*W*
     ;   *10-8*V**10*W**9+V**8*(19*V**2-18*V+18)*W**8-4*V**8*(7*V
     <   **2-12*V+12)*W**7+4*V**6*(6*V**4-14*V**3+25*V**2-22*V+11
     =   )*W**6-4*V**6*(3*V**4-9*V**3+29*V**2-40*V+20)*W**5+V**4*
     >   (3*V**6-6*V**5+56*V**4-144*V**3+182*V**2-132*V+44)*W**4-
     ?   4*(V-1)*V**4*(V**4+5*V**3-17*V**2+24*V-12)*W**3+2*(V-1)*
     @   *2*V**2*(5*V**4-8*V**3+17*V**2-18*V+9)*W**2-4*(V-1)**3*V
     1   **2*(V**2+2*V-2)*W+(V-1)**4*(3*V**2-2*V+2))-2*(2*V**8*(V
     2   **4-2*V**3+4*V**2-4*V+2)*W**10-V**8*(9*V**4-20*V**3+40*V
     3   **2-40*V+20)*W**9+V**6*(20*V**6-47*V**5+73*V**4-36*V**3-
     4   22*V**2+48*V-16)*W**8-2*V**6*(13*V**6-28*V**5+20*V**4+48
     5   *V**3-104*V**2+96*V-32)*W**7+2*V**4*(10*V**8-17*V**7-5*V
     6   **6+84*V**5-130*V**4+72*V**3+32*V**2-48*V+12)*W**6-V**4*
     7   (9*V**8-8*V**7-10*V**6+52*V**5+6*V**4-240*V**3+416*V**2-
     8   288*V+72)*W**5+V**2*(2*V**10+5*V**9+3*V**8-78*V**7+262*V
     9   **6-442*V**5+390*V**4-112*V**3-92*V**2+80*V-16)*W**4-2*(
     :   V-1)*V**2*(2*V**8+7*V**7-35*V**6+64*V**5-36*V**4-40*V**3
     ;   +88*V**2-64*V+16)*W**3+2*(V-1)**2*(4*V**8-7*V**7+6*V**6+
     <   12*V**5-29*V**4+22*V**3+2*V**2-8*V+2)*W**2-(V-1)**3*(4*V
     =   **6-3*V**5+7*V**4-4*V**3-8*V**2+12*V-4)*W+(V-1)**4*(2*V*
     >   *4-3*V**3+5*V**2-4*V+2)))-8*GTR*N**2*(V-1)**2*VC*W*(V*W-
     ?   1)*(V*W-V+1)*(CQ*V**2*(2*V**10*W**10-9*V**10*W**9+V**8*(
     @   21*V**2-10*V+10)*W**8-2*V**8*(15*V**2-14*V+14)*W**7+2*V*
     1   *6*(13*V**4-15*V**3+25*V**2-20*V+10)*W**6-V**6*(13*V**4-
     2   10*V**3+44*V**2-68*V+34)*W**5+V**4*(3*V**6+8*V**5+16*V**
     3   4-68*V**3+84*V**2-60*V+20)*W**4-2*(V-1)*V**4*(3*V**4+9*V
     4   **3-19*V**2+20*V-10)*W**3+2*(V-1)**2*V**2*(5*V**4-V**3+6
     5   *V**2-10*V+5)*W**2-(V-1)**3*V**2*(6*V**2+5*V-5)*W+(V-1)*
     6   *4*(3*V**2-2*V+2))-2*V**8*(2*V**4-3*V**3+5*V**2-4*V+2)*W
     7   **10+V**8*(23*V**4-30*V**3+50*V**2-40*V+20)*W**9-V**6*(5
     8   8*V**6-61*V**5+79*V**4-20*V**3-30*V**2+48*V-16)*W**8+8*V
     9   **6*(10*V**6-5*V**5-V**4+20*V**3-30*V**2+24*V-8)*W**7-4*
     :   V**4*(16*V**8+9*V**7-29*V**6+57*V**5-65*V**4+27*V**3+19*
     ;   V**2-24*V+6)*W**6+V**4*(29*V**8+80*V**7-82*V**6-16*V**5+
     <   130*V**4-348*V**3+452*V**2-288*V+72)*W**5-V**2*(6*V**10+
     =   59*V**9+15*V**8-234*V**7+392*V**6-482*V**5+366*V**4-80*V
     >   **3-100*V**2+80*V-16)*W**4+2*(V-1)*V**2*(9*V**8+36*V**7-
     ?   56*V**6+40*V**5-4*V**4-64*V**3+96*V**2-64*V+16)*W**3-2*(
     @   V-1)**2*(14*V**8+7*V**7-11*V**6+17*V**5-29*V**4+19*V**3+
     1   3*V**2-8*V+2)*W**2+(V-1)**3*(18*V**6-V**5+3*V**4-10*V**2
     2   +12*V-4)*W-(V-1)**4*(6*V**4-5*V**3+7*V**2-4*V+2))+4*GTR*
     3   (V-1)**2*V**2*VC*W*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**3
     4   *W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(CQ*(2*
     5   V**6*W**6-6*V**6*W**5+3*V**4*(3*V**2-2*V+2)*W**4-4*V**4*
     6   (2*V**2-3*V+3)*W**3+V**2*(3*V**4-4*V**3+10*V**2-12*V+6)*
     7   W**2-2*(V-1)*V**2*(V**2+3*V-3)*W+(V-1)**2*(3*V**2-2*V+2)
     8   )-4*(V**4*(V**2-V+1)*W**6-3*V**4*(V**2-V+1)*W**5+2*V**2*
     9   (V+1)*(2*V-1)*(V**2-V+1)*W**4-V**2*(V+2)*(3*V-2)*(V**2-V
     :   +1)*W**3+(V**2+V-1)*(V**4+2*V**3-3*V**2+2*V-1)*W**2-(V-1
     ;   )*(2*V**4+3*V**3-4*V**2+2*V-1)*W+(V-1)**2*(2*V**2-V+1)))
     <   )/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
 
 
      LVW = (16*N**5*(V-1)*VC*(V**5*W**5-5*V**4*W**4+10*V**3*W**3
     1   -10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(
     2   V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(
     3   V-1)**5)*(4*V**6*W**5-4*V**5*(V+1)*W**4+V**4*(8*V**2-11*
     4   V+19)*W**3-2*V**3*(V**2+V+4)*W**2-(V**2-V+1)*(2*V**3-5*V
     5   **2-9*V+4)*W-4*(V**2-V+1)**2)-16*GTR*N**2*(V-1)**2*V**2*
     6   VC*W*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**5*W**5-5*V**4
     7   *W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V
     8   -1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**
     9   2+5*(V-1)**4*V*W-(V-1)**5)+8*GTR*(V-1)**2*V**2*VC*W*(V**
     :   2*W**2-2*(V-1)*V*W+2*V**2-2*V+1)*(V**5*W**5-5*V**4*W**4+
     ;   10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V*
     <   *4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V
     =   -1)**4*V*W-(V-1)**5)-8*GTR*N**4*(V-1)**2*V*VC*W*(V**2*W-
     >   (V-1)**2)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W*
     ?   *2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**
     @   3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))*
     1   LOG(1-V*W)/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)
     2   **5)
 
 
      LTVW = (16*N**5*(V-1)*VC*(V*W-V+1)*(V**5*W**5-5*V**4*W**4+1
     1   0*V**3*W**3-10*V**2*W**2+5*V*W-1)*((V-1)*V*W*(V**7*W**6+
     2   2*(V-2)*V**6*W**5-V**4*(9*V**3-16*V**2+6*V-1)*W**4+4*(V-
     3   1)*V**3*(3*V**3-3*V**2+2*V-1)*W**3-(V-1)**2*V**2*(9*V**3
     4   -6*V**2+11*V-6)*W**2+2*(V-1)**3*V*(V**3-V**2+4*V-2)*W+(V
     5   -1)**4*(V+1)*(V**2-V+1))+4*CQ*(W-V+1)*(V*W+1)*(V*W+(V-1)
     6   **2)*(V**2*W-(V-1)**2)*(V**2*W-V+1)*(V**2*W**2-2*(V-1)*V
     7   *W+(V-1)**2))+8*GTR*N**4*(V-1)**2*V*VC*W*(V*W-V+1)*(V**5
     8   *W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(-V
     9   **7*W**6+V**6*(V+1)*W**5+CQ*V*(V**2*W**2-2*V**2*W+V**2+1
     :   )*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)+V**4*(V**3-4
     ;   *V**2+2*V-1)*W**4-2*(V-2)*(V-1)**3*V**3*W**3+(V-1)**2*V*
     <   *2*(V**3-6*V**2+7*V-6)*W**2+(V-1)**4*V*(V**2+V-4)*W-(V-1
     =   )**4*(V+1)*(V**2-V+1))+8*GTR*(V-1)**2*V**2*VC*W*(CQ*(V**
     >   2*W**2+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)+4*(V-1)*V**
     ?   2*(W-1)*W)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(
     @   V-1)**3)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**
     1   2+5*V*W-1)-16*GTR*N**2*(V-1)**2*V**2*VC*W*(V*W-V+1)*(V**
     2   2*W**2-2*V**2*W+V**2+1)*(CQ*(V**2*W**2+(V-1)**2)*(V**2*W
     3   **2-(V-1)*V*W+(V-1)**2)-(V-1)*V*W*(3*V**2*W**2-4*(V-1)*V
     4   *W+3*(V-1)**2))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V
     5   **2*W**2+5*V*W-1))*LOG(V*W-V+1)/(N**2*(V-1)**3*V**3*W**
     6   2*(V*W-1)**5*(V*W-V+1)**5)
 
 
      LW = (16*N**5*(V-1)*VC*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4
     1   *V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**
     2   2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**6*(V**2-V+1)*W**6-V**6*
     3   (8*V**2-13*V+13)*W**5+V**4*(14*V**4-29*V**3+30*V**2-2*V+
     4   1)*W**4-2*V**2*(7*V**6-15*V**5+12*V**4+8*V**3-9*V**2+6*V
     5   -2)*W**3+V**2*(10*V**6-24*V**5+26*V**4-3*V**3-V**2+3*V-1
     6   )*W**2-4*(V**2-V+1)*(V**6-V**5-V**4+3*V**3+V**2-3*V+1)*W
     7   +2*(V-1)*(V**2-V+1)**2*(2*V**2-3*V+3))-4*GTR*N**4*(V-1)*
     8   *2*VC*W*(CQ*(V**2-2*V+2)**2*(2*W**2-2*W+1)-2*V**2*(V**2*
     9   W-V+1))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2
     :   +5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*
     ;   W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)+8*GT
     <   R*N**2*(V-1)**2*VC*W*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-
     =   4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W*
     >   *2-4*(V-1)**3*V*W+(V-1)**4)*(V**2*(2*V**4*W**4-V**2*(5*V
     ?   **2-2*V+2)*W**3+3*V**2*(V**2+V-1)*W**2-2*(V**4+V**2-2*V+
     @   1)*W+(V-1)*(2*V**2-V+1))+CQ*(V**2-2*V+2)*(V**2-V+1)*(V*W
     1   -1)*(V*W-V+1)*(2*W**2-2*W+1))-4*GTR*(V-1)**2*V**2*VC*W*(
     2   V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*
     3   (V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1
     4   )**4)*(2*(V**4*W**4-V**2*(V**2+2*V-2)*W**3+5*(V-1)*V**2*
     5   W**2-2*(V-1)*(2*V**2-V+1)*W+2*(V-1)*(V**2-V+1))+CQ*(V**2
     6   -2*V+2)*(V*W-1)*(V*W-V+1)*(2*W**2-2*W+1)))*LOG(W)/(N**2
     7   *(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
 
 
      CVC = (-8*N**5*(V-1)*VC*(V*W-1)*(V*W-V+1)*((V-1)*V**8*(17*V
     1   **4-50*V**3+90*V**2-80*V+40)*W**10-2*(V-1)*V**8*(34*V**4
     2   -114*V**3+213*V**2-198*V+99)*W**9+V**6*(11*V**8+103*V**7
     3   -571*V**6+1145*V**5-935*V**4-13*V**3+751*V**2-640*V+160)
     4   *W**8-V**6*(33*V**8+71*V**7-602*V**6+806*V**5+869*V**4-3
     5   296*V**3+4048*V**2-2528*V+632)*W**7+V**4*(33*V**10+81*V*
     6   *9-642*V**8+796*V**7+1121*V**6-3554*V**5+3350*V**4-416*V
     7   **3-1696*V**2+1200*V-240)*W**6-V**4*(11*V**10+123*V**9-6
     8   51*V**8+1137*V**7-719*V**6+743*V**5-3309*V**4+6872*V**3-
     9   7028*V**2+3540*V-708)*W**5+(V-1)*V**2*(61*V**10-124*V**9
     :   +13*V**8-168*V**7+1657*V**6-3402*V**5+3178*V**4-792*V**3
     ;   -1002*V**2+800*V-160)*W**4-(V-1)**2*V**2*(84*V**8-166*V*
     <   *7-183*V**6+818*V**5-397*V**4-888*V**3+1752*V**2-1248*V+
     =   312)*W**3+(V-1)**3*(86*V**8-285*V**7+382*V**6-24*V**5-37
     >   3*V**4+350*V**3+70*V**2-160*V+40)*W**2-(V-1)**4*(29*V**6
     ?   -47*V**5+87*V**4-42*V**3-74*V**2+114*V-38)*W+17*(V-1)**5
     @   *(V**2-V+1)**2)+8*GTR*N**4*(V-1)*VC*(V*W-1)*(V*W-V+1)*((
     1   V-1)*V**8*(7*V**4-34*V**3+78*V**2-88*V+44)*W**10-(V-1)*V
     2   **8*(31*V**4-156*V**3+372*V**2-432*V+216)*W**9+V**6*(4*V
     3   **8+53*V**7-368*V**6+922*V**5-1015*V**4+172*V**3+764*V**
     4   2-704*V+176)*W**8-2*V**6*(6*V**8+23*V**7-194*V**6+314*V*
     5   *5+257*V**4-1460*V**3+2092*V**2-1376*V+344)*W**7+V**4*(1
     6   2*V**10+45*V**9-381*V**8+602*V**7+694*V**6-3226*V**5+367
     7   0*V**4-640*V**3-1820*V**2+1320*V-264)*W**6-V**4*(4*V**10
     8   +51*V**9-351*V**8+876*V**7-952*V**6+892*V**5-3060*V**4+6
     9   976*V**3-7504*V**2+3840*V-768)*W**5+(V-1)*V**2*(23*V**10
     :   -95*V**9+287*V**8-834*V**7+2174*V**6-3702*V**5+3362*V**4
     ;   -768*V**3-1128*V**2+880*V-176)*W**4-2*(V-1)**2*V**2*(9*V
     <   **8+20*V**7-213*V**6+422*V**5-133*V**4-564*V**3+972*V**2
     =   -672*V+168)*W**3+(V-1)**3*(34*V**8-138*V**7+173*V**6+72*
     >   V**5-347*V**4+250*V**3+122*V**2-176*V+44)*W**2+(V-1)**4*
     ?   (2*V**6-29*V**5+33*V**4-48*V**3+124*V**2-120*V+40)*W+(V-
     @   1)**5*(V**2-V+1)*(7*V**2-4*V+4))+24*GTR*N**2*(V-1)**2*VC
     1   *(V*W-1)*(V*W-V+1)*(2*(V-1)*V**8*(5*V**2-6*V+6)*W**10+V*
     2   *8*(V**4-46*V**3+106*V**2-120*V+60)*W**9-V**6*(4*V**6-95
     3   *V**5+175*V**4-112*V**3-64*V**2+144*V-48)*W**8+6*V**6*(V
     4   **6-18*V**5+14*V**4+40*V**3-100*V**2+96*V-32)*W**7-2*V**
     5   4*(2*V**8-49*V**7+35*V**6+142*V**5-320*V**4+198*V**3+102
     6   *V**2-144*V+36)*W**6+V**4*(V**8-80*V**7+166*V**6-136*V**
     7   5+194*V**4-756*V**3+1260*V**2-864*V+216)*W**5+(V-1)*V**2
     8   *(41*V**8-88*V**7+218*V**6-436*V**5+610*V**4-336*V**3-11
     9   2*V**2+192*V-48)*W**4-2*(V-1)*V**2*(5*V**8-11*V**7+49*V*
     :   *6-80*V**5+2*V**4+180*V**3-284*V**2+192*V-48)*W**3+2*(V-
     ;   1)**2*(4*V**8+2*V**7-6*V**6-11*V**5+47*V**4-33*V**3-17*V
     <   **2+24*V-6)*W**2-(V-1)**3*(10*V**6-29*V**5+35*V**4-24*V*
     =   *3+42*V**2-36*V+12)*W-(V-1)**5*V**2)-24*AL*GTR*N**2*(V-1
     >   )**2*VC*(V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2
     ?   *V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5
     @   )*(V**4*(3*V**4-6*V**3+10*V**2-8*V+4)*W**6-V**3*(4*V**5+
     1   3*V**4-14*V**3+32*V**2-28*V+16)*W**5+V**2*(3*V**6+5*V**5
     2   -4*V**4+30*V**2-32*V+24)*W**4-2*V*(3*V**6+V**5-3*V**4+10
     3   *V**3-4*V+8)*W**3+(10*V**6-12*V**5+11*V**4+10*V**3-10*V*
     4   *2+8*V+4)*W**2-(6*V**5-10*V**4+15*V**3-6*V**2+4)*W+3*V**
     5   4-5*V**3+6*V**2-4*V+2)+12*AL*GTR*N**4*(V-1)**2*VC*(V*W-1
     6   )*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*
     7   (V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(V**4*(3*V**
     8   4-8*V**3+16*V**2-16*V+8)*W**6-4*V**3*(V**5-4*V**3+12*V**
     9   2-14*V+8)*W**5+V**2*(3*V**6+2*V**5-9*V**4+8*V**3+36*V**2
     :   -64*V+48)*W**4-4*V*(V**6-2*V**4+8*V**3-4*V**2-4*V+8)*W**
     ;   3+(10*V**6-20*V**5+21*V**4+8*V**3-24*V**2+16*V+8)*W**2-4
     <   *(V**5-3*V**4+6*V**3-4*V**2+2)*W+3*V**4-6*V**3+9*V**2-8*
     =   V+4)-24*GTR*(V-1)**2*V**2*VC*(V**2*W**2+V**2-2*V+2)*(V**
     >   5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V
     ?   **5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1
     @   )**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)+12*AL*GTR*(V-1)*
     1   *2*V**2*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**2*(3*V**2
     2   -4*V+4)*W**4-4*V*(V**3-V+2)*W**3+(3*V**4-2*V**2+4*V+4)*W
     3   **2-2*(V**3+2)*W+3*V**2-4*V+3)*(V**5*W**5-5*(V-1)*V**4*W
     4   **4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)*
     5   *4*V*W-(V-1)**5)-24*CQ*GTR*(V-1)**2*V**2*VC*(V**2*W**2+(
     6   V-1)**2)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W*
     7   *2-4*(V-1)**3*V*W+(V-1)**4)*(V**5*W**5-5*V**4*W**4+10*V*
     8   *3*W**3-10*V**2*W**2+5*V*W-1)-24*CQ*GTR*N**4*(V-1)**2*V*
     9   *2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4+2*(V-1
     :   )**2*V**2*W**2+(V-1)**4)*(V**5*W**5-5*V**4*W**4+10*V**3*
     ;   W**3-10*V**2*W**2+5*V*W-1)+48*CQ*GTR*N**2*(V-1)**2*V**2*
     <   VC*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)
     =   *(V**2*W**2-(V-1)*V*W+(V-1)**2)*(V**5*W**5-5*V**4*W**4+1
     >   0*V**3*W**3-10*V**2*W**2+5*V*W-1))/(N**2*(V-1)**3*V**3*W
     ?   *(V*W-1)**5*(V*W-V+1)**5)/3.D0
 
 
      LM = LOG(S/M**2)*(-32*N**5*VC*(V**4*W**4-4*(V-1)*V**3*W**3
     1   +6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**10*W*
     2   *9-V**9*(V+5)*W**8+V**4*(3*V**6-2*V**5+21*V**4-7*V**3+6*
     3   V**2-3*V+1)*W**7-V**3*(V**7+9*V**6-9*V**5+47*V**4-22*V**
     4   3+21*V**2-11*V+4)*W**6+V**2*(V**8+18*V**6-12*V**5+53*V**
     5   4-18*V**3+24*V**2-14*V+6)*W**5-V*(V**9+V**8+20*V**6-3*V*
     6   *5+27*V**4+10*V**3+6*V**2-6*V+4)*W**4+(4*V**9-6*V**8+13*
     7   V**7+2*V**6+15*V**5-2*V**4+25*V**3-6*V**2+V+1)*W**3-(6*V
     8   **8-14*V**7+30*V**6-27*V**5+34*V**4-17*V**3+18*V**2-3*V+
     9   1)*W**2+V*(V**2-V+1)*(4*V**4-7*V**3+13*V**2-7*V+8)*W-(V*
     :   *2-V+1)**2*(V**2-V+2))+8*GTR*N**2*(V-1)*VC*W*(V**4*W**4-
     ;   4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V
     <   -1)**4)*(V**4*(3*V**4-6*V**3+10*V**2-8*V+4)*W**6-V**3*(4
     =   *V**5+3*V**4-14*V**3+32*V**2-28*V+16)*W**5+V**2*(3*V**6+
     >   5*V**5-4*V**4+30*V**2-32*V+24)*W**4-2*V*(3*V**6+V**5-3*V
     ?   **4+10*V**3-4*V+8)*W**3+(10*V**6-12*V**5+11*V**4+10*V**3
     @   -10*V**2+8*V+4)*W**2-(6*V**5-10*V**4+15*V**3-6*V**2+4)*W
     1   +3*V**4-5*V**3+6*V**2-4*V+2)-4*GTR*N**4*(V-1)*VC*W*(V**4
     2   *W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*
     3   V*W+(V-1)**4)*(V**4*(3*V**4-8*V**3+16*V**2-16*V+8)*W**6-
     4   4*V**3*(V**5-4*V**3+12*V**2-14*V+8)*W**5+V**2*(3*V**6+2*
     5   V**5-9*V**4+8*V**3+36*V**2-64*V+48)*W**4-4*V*(V**6-2*V**
     6   4+8*V**3-4*V**2-4*V+8)*W**3+(10*V**6-20*V**5+21*V**4+8*V
     7   **3-24*V**2+16*V+8)*W**2-4*(V**5-3*V**4+6*V**3-4*V**2+2)
     8   *W+3*V**4-6*V**3+9*V**2-8*V+4)-4*GTR*(V-1)*V**2*VC*W*(V*
     9   *2*W**2-2*V*W+1)*(V**2*(3*V**2-4*V+4)*W**4-4*V*(V**3-V+2
     :   )*W**3+(3*V**4-2*V**2+4*V+4)*W**2-2*(V**3+2)*W+3*V**2-4*
     ;   V+3)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4
     <   *(V-1)**3*V*W+(V-1)**4))/(N**2*(V-1)**2*V**3*W**2*(V*W-1
     =   )**4*(V*W-V+1)**4)
 
 
      LMP = LOG(S/MP**2)*(-32*N**5*VC*(V**4*W**4-4*V**3*W**3+6*V
     1   **2*W**2-4*V*W+1)*(V**10*W**9-V**9*(6*V-5)*W**8+V**8*(18
     2   *V**2-31*V+15)*W**7-(V-1)*V**7*(35*V**2-57*V+30)*W**6+(V
     3   -1)*V**4*(47*V**5-120*V**4+114*V**3-41*V**2-2*V+1)*W**5-
     4   (V-1)**3*V**3*(45*V**4-66*V**3+39*V**2+3*V-4)*W**4+(V-1)
     5   **3*V**2*(32*V**5-75*V**4+72*V**3-24*V**2-8*V+6)*W**3-(V
     6   -1)**4*V*(17*V**5-38*V**4+40*V**3-16*V**2-2*V+4)*W**2+(V
     7   -1)**5*(V**2-V+1)*(6*V**3-7*V**2+3*V+1)*W-(V-1)**6*(V**2
     8   -V+1)**2)-4*GTR*N**4*(V-1)*V**2*VC*W*(V**2*W**2-2*V**2*W
     9   +V**2+1)*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)*(V**4
     :   *W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)+8*GTR*N**2*(V-1)*
     ;   V**2*VC*W*(V**2*W**2+(V-1)**2)*(V**2*W**2-(V-1)*V*W+(V-1
     <   )**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**4*W**4-4*V**3*W**3
     =   +6*V**2*W**2-4*V*W+1)-4*GTR*(V-1)*V**2*VC*W*(V**2*W**2+(
     >   V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V
     ?   **2*W+V**2+1)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1
     @   ))/(N**2*(V-1)**2*V**3*W**2*(V*W-1)**4*(V*W-V+1)**4)
 
      STRUV15=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION STRUV16(W,V,X3,S)
C     CALCUL DES TERMES CONSTANTS
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/FONLLSCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
 
      LW1 = LOG(1-W)*(-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W
     1   **2-2*(V-1)*V*W+(V-1)**2)*(CQ*(16*V**11*W**11-16*V**10*(
     2   5*V-1)*W**10+V**6*(224*V**5-141*V**4+23*V**3+4*V**2+2)*W
     3   **9-V**6*(384*V**5-323*V**4-28*V**3+79*V**2+8)*W**8+V**4
     4   *(416*V**7-280*V**6-431*V**5+470*V**4-88*V**3+5*V**2+6*V
     5   -6)*W**7-V**4*(320*V**7-120*V**6-770*V**5+829*V**4-170*V
     6   **3-37*V**2+18*V-18)*W**6+V**2*(192*V**9-35*V**8-679*V**
     7   7+670*V**6+129*V**5-358*V**4+122*V**3-15*V**2-12*V+6)*W*
     8   *5-V**2*(80*V**9+83*V**8-710*V**7+903*V**6-244*V**5-203*
     9   V**4+98*V**3+9*V**2-24*V+12)*W**4+(V-1)*(16*V**10+144*V*
     :   *9-413*V**8+265*V**7+223*V**6-298*V**5+67*V**4+29*V**3-1
     ;   5*V**2-4*V+2)*W**3-(V-1)**2*(48*V**8-293*V**6+522*V**5-3
     <   20*V**4+60*V**3+7*V**2+2*V-2)*W**2+(V-1)**3*(12*V**2-13*
     =   V+5)*(4*V**4-5*V**3-V**2+5*V+1)*W-4*(V-1)**4*(2*V**2-4*V
     >   +3)*(2*V**2-2*V+1))-32*V**11*W**11+16*V**9*(8*V**2+V-1)*
     ?   W**10-V**6*(288*V**5-51*V**4+13*V**3-44*V**2+20*V-2)*W**
     @   9+V**6*(448*V**5-335*V**4+262*V**3-243*V**2+116*V-8)*W**
     1   8-V**4*(448*V**7-398*V**6-13*V**5+294*V**4-276*V**3+189*
     2   V**2-66*V+6)*W**7+V**4*(320*V**7-206*V**6-732*V**5+1539*
     3   V**4-1246*V**3+621*V**2-210*V+18)*W**6-V**2*(224*V**9-23
     4   9*V**8-489*V**7+952*V**6-29*V**5-768*V**4+590*V**3-267*V
     5   **2+72*V-6)*W**5+V**2*(128*V**9-83*V**8-666*V**7+1167*V*
     6   *6-82*V**5-1091*V**4+948*V**3-425*V**2+124*V-12)*W**4-(V
     7   -1)*(32*V**10+208*V**9-815*V**8+903*V**7-147*V**6-96*V**
     8   5-197*V**4+229*V**3-107*V**2+24*V-2)*W**3+(V-1)**2*(96*V
     9   **8-64*V**7-439*V**6+938*V**5-648*V**4+172*V**3+V**2-10*
     :   V+2)*W**2-(V-1)**3*(96*V**6-240*V**5+199*V**4+26*V**3-76
     ;   *V**2+18*V+9)*W+8*(V-1)**4*(2*V**2-4*V+3)*(2*V**2-2*V+1)
     <   )+2*N**2*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**5*W**5-5*(V-1)
     =   *V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5
     >   *(V-1)**4*V*W-(V-1)**5)*(CQ*(V**3*(3*V**4-3*V**3+4*V**2-
     ?   2*V+2)*W**6-V**2*(4*V**5+5*V**4-5*V**3+10*V**2-4*V+6)*W*
     @   *5+V*(3*V**6+7*V**5+4*V**4-3*V**3+9*V**2+6)*W**4-(9*V**6
     1   +V**5+7*V**4-8*V**3+7*V**2+4*V+2)*W**3+(13*V**5-13*V**4+
     2   17*V**3-16*V**2+7*V+2)*W**2-(9*V**4-13*V**3+12*V**2-9*V+
     3   3)*W+2*(V-1)*(2*V**2-2*V+1))-(V*W-1)*(4*V**3*(V**3+V**2-
     4   V+1)*W**5-2*V**2*(3*V**4+2*V**3+5*V**2-2*V+4)*W**4+V*(4*
     5   V**5+8*V**4-2*V**3+5*V**2+9*V+4)*W**3-V*(16*V**4-22*V**3
     6   +23*V**2-15*V+14)*W**2+(16*V**4-28*V**3+25*V**2-14*V+5)*
     7   W-4*(V-1)*(2*V**2-2*V+1)))-(V-1)**2*VC*W*(V**2*W**2-2*V*
     8   W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-
     9   10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(CQ*(V**3
     :   *(3*V**3-2*V**2+2*V-2)*W**5-2*V**2*(2*V**4+3*V**3-2*V**2
     ;   +2*V-3)*W**4+V*(3*V**5+7*V**4+4*V**3-V**2-6)*W**3-(7*V**
     <   5+3*V**4+V**3+V**2-4*V-2)*W**2+(7*V**4-5*V**3+3*V**2-V-2
     =   )*W-(V-1)*(V**2+1))-V**3*(5*V**3-2*V**2+2*V-2)*W**5+2*V*
     >   *2*(4*V**4+6*V**3-2*V**2+2*V-3)*W**4-V*(5*V**5+19*V**4+1
     ?   2*V**3-V**2-6)*W**3+(13*V**5+15*V**4+9*V**3+V**2-4*V-2)*
     @   W**2-(13*V**4-V**3+9*V**2-V-2)*W+(3*V-1)*(V**2+1)))/(N**
     1   2*(V-1)**3*V**2*W**2*(V*W-1)**5*(V*W-V+1)**5)
 
      LV1 = LOG(1-V)*(N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**4*W*
     1   *4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W
     2   +(V-1)**4)*(2*(V**2*W**2-2*V*W+1)*(16*V**7*W**7-40*V**7*
     3   W**6+2*V**4*(32*V**3-19*V**2+7*V+4)*W**5-V**4*(56*V**3-5
     4   7*V**2+23*V+14)*W**4+V**2*(24*V**5-20*V**4-9*V**3+31*V**
     5   2-4*V-6)*W**3-V**2*(8*V**5-11*V**4+17*V**3-21*V**2+23*V-
     6   12)*W**2+(V-1)*(8*V**5-3*V**4+3*V**2-2*V+2)*W-4*(V-1)**2
     7   *V*(2*V**2-2*V+1))+CQ*(V-1)*(V*W-V+1)*(V**6*(16*V-15)*W*
     8   *6+2*V**5*(8*V**2-41*V+30)*W**5-V**3*(46*V**3-150*V**2+8
     9   7*V+4)*W**4+4*V**2*(10*V**3-30*V**2+12*V+3)*W**3+4*V*(8*
     :   V**2+V-3)*W**2-4*(4*V**3-2*V**2+3*V-1)*W+4*(2*V**2-2*V+1
     ;   )))-2*N**2*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**5*W**5-5*(V-
     <   1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2
     =   +5*(V-1)**4*V*W-(V-1)**5)*(CQ*(V-1)*(V**6*W**6-V**5*(2*V
     >   +3)*W**5+2*V**3*(V**3+2*V**2+3*V-1)*W**4-V**2*(6*V**3+4*
     ?   V**2+9*V-6)*W**3+V*(10*V**3-2*V**2+11*V-6)*W**2-2*(4*V**
     @   3-2*V**2+3*V-1)*W+2*(2*V**2-2*V+1))+(2*V**2*(5*V-1)*W**3
     1   -V**2*(7*V+1)*W**2+V*(V**3+4*V**2-3*V+2)*W-4*(V-1)*(2*V*
     2   *2-2*V+1))*(V**3*W**3-3*V**2*W**2+3*V*W-1))+(V-1)**2*VC*
     3   W*(V**2*W**2-2*V*W+1)*(2*(V-1)*(V**2-2*V+2)*(V**3*W**3-3
     4   *V**2*W**2+3*V*W-1)+CQ*V**2*W*(V**2*W**2-2*V*W+2)*(V**2*
     5   W**2-2*V**2*W+2*V**2-2*V+1))*(V**5*W**5-5*(V-1)*V**4*W**
     6   4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4
     7   *V*W-(V-1)**5))/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**5*(V*W
     8   -V+1)**5)
 
      LV = LOG(V)*(-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2
     1   -2*(V-1)*V*W+(V-1)**2)*(CQ*(16*V**11*W**11-16*V**10*(5*V
     2   -1)*W**10+V**9*(224*V**2-143*V+23)*W**9-V**8*(384*V**3-3
     3   31*V**2-28*V+63)*W**8+V**6*(416*V**5-293*V**4-437*V**3+4
     4   50*V**2-100*V+4)*W**7-V**6*(320*V**5-131*V**4-788*V**3+8
     5   25*V**2-206*V-12)*W**6+(V-1)*V**4*(192*V**6+152*V**5-548
     6   *V**4+127*V**3+226*V**2-113*V+12)*W**5-(V-1)*V**4*(80*V*
     7   *6+162*V**5-560*V**4+341*V**3+97*V**2-119*V+15)*W**4+4*(
     8   V-1)*V**2*(4*V**8+36*V**7-104*V**6+64*V**5+56*V**4-78*V*
     9   *3+19*V**2+7*V-3)*W**3-4*(V-1)**2*V**2*(12*V**6-74*V**4+
     :   130*V**3-81*V**2+14*V+2)*W**2+4*(V-1)**3*(12*V**6-28*V**
     ;   5+18*V**4+12*V**3-15*V**2+3*V+1)*W-4*(V-1)**4*(2*V**2-4*
     <   V+3)*(2*V**2-2*V+1))-32*V**11*W**11+144*V**11*W**10-V**6
     =   *(352*V**5-109*V**4+71*V**3-44*V**2+20*V-2)*W**9+V**6*(5
     >   76*V**5-435*V**4+316*V**3-185*V**2+80*V-8)*W**8-V**4*(62
     ?   4*V**7-552*V**6+5*V**5+446*V**4-380*V**3+189*V**2-66*V+6
     @   )*W**7+V**4*(480*V**7-340*V**6-858*V**5+1885*V**4-1470*V
     1   **3+627*V**2-198*V+18)*W**6-V**2*(304*V**9-219*V**8-927*
     2   V**7+1698*V**6-539*V**5-674*V**4+606*V**3-267*V**2+72*V-
     3   6)*W**5+V**2*(144*V**9+7*V**8-1102*V**7+1815*V**6-368*V*
     4   *5-1315*V**4+1234*V**3-539*V**2+144*V-12)*W**4-(V-1)*(32
     5   *V**10+240*V**9-913*V**8+873*V**7+243*V**6-586*V**5+47*V
     6   **4+181*V**3-107*V**2+24*V-2)*W**3+(V-1)**2*(96*V**8-64*
     7   V**7-553*V**6+1186*V**5-808*V**4+176*V**3+35*V**2-22*V+2
     8   )*W**2-(V-1)**3*(96*V**6-272*V**5+225*V**4+64*V**3-126*V
     9   **2+36*V+9)*W+8*(V-1)**4*(2*V**2-5*V+4)*(2*V**2-2*V+1))+
     :   2*N**2*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**5*W**5-5*(V-1)*V
     ;   **4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(
     <   V-1)**4*V*W-(V-1)**5)*(CQ*(V-1)*(V**6*W**6-V**5*(2*V+3)*
     =   W**5+2*V**3*(V**3+2*V**2+3*V-1)*W**4-V**2*(6*V**3+4*V**2
     >   +9*V-6)*W**3+V*(10*V**3-2*V**2+11*V-6)*W**2-2*(4*V**3-2*
     ?   V**2+3*V-1)*W+2*(2*V**2-2*V+1))-(V*W-1)*(V**3*(3*V**3+5*
     @   V**2-4*V+4)*W**5-2*V**2*(2*V**4+3*V**3+5*V**2-2*V+4)*W**
     1   4+V*(3*V**2+3*V+1)*(V**3+2*V**2-3*V+4)*W**3-V*(16*V**4-2
     2   0*V**3+21*V**2-15*V+14)*W**2+(17*V**4-29*V**3+26*V**2-15
     3   *V+5)*W-4*(V-1)*(2*V**2-2*V+1)))-(V-1)**2*VC*W*(V**2*W**
     4   2-2*V*W+1)*(-V**3*(3*V**3-2*V**2+2*V-2)*W**5+2*V**2*(2*V
     5   **4+3*V**3-2*V**2+2*V-3)*W**4-V*(3*V**5+7*V**4+4*V**3-V*
     6   *2-6)*W**3+CQ*V**2*W*(V**2*W**2-2*V*W+2)*(V**2*W**2-2*V*
     7   *2*W+2*V**2-2*V+1)+(7*V**5+3*V**4+V**3+V**2-4*V-2)*W**2-
     8   (7*V**4-5*V**3+3*V**2-V-2)*W+(V-1)*(V**2+1))*(V**5*W**5-
     9   5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2
     :   *W**2+5*(V-1)**4*V*W-(V-1)**5))/(N**2*(V-1)**3*V**2*W**2
     ;   *(V*W-1)**5*(V*W-V+1)**5)
 
      LVW = (-2*N**4*(V-1)*VC*(V**5*W**5-5*V**4*W**4+10*V**3*W**3
     1   -10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(
     2   V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(
     3   V-1)**5)*(16*V**5*W**5-16*V**4*(V+1)*W**4+V**3*(32*V**2-
     4   43*V+43)*W**3-V**2*(10*V**2-7*V+13)*W**2-(6*V**4-21*V**3
     5   +12*V**2+3*V-4)*W-4*(V-1)*(2*V**2-2*V+1))+2*N**2*(V-1)*V
     6   C*(V**3*(V+7)*W**3-V**2*(2*V**2+V+5)*W**2+(2*V**4+V**3-V
     7   +2)*W-4*(V-1)*(2*V**2-2*V+1))*(V**5*W**5-5*V**4*W**4+10*
     8   V**3*W**3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*
     9   W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)
     :   **4*V*W-(V-1)**5)-2*(V-1)**2*VC*W*(V**3*W**2-2*V**3*W+(2
     ;   *V-1)*(V**2-V+2))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10
     <   *V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1
     =   )**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1
     >   )**5))*LOG(1-V*W)/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**5*(
     ?   V*W-V+1)**5)
 
      LTVW = (-2*N**4*(V-1)*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(
     1   V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)
     2   *(8*CQ*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**6*W*
     3   *6-3*(V-1)*V**5*W**5+6*(V-1)**2*V**4*W**4-7*(V-1)**3*V**
     4   3*W**3+6*(V-1)**4*V**2*W**2-3*(V-1)**5*V*W+(V-1)**6)+(V-
     5   1)*W*(5*V**6*W**5+V**5*(5*V-7)*W**4-(V-1)*V**3*(16*V**3-
     6   8*V**2+8*V-3)*W**3+(V-1)**2*V**2*(32*V**3-40*V**2+36*V-9
     7   )*W**2-(V-1)**3*V*(16*V**3-27*V**2+27*V-9)*W-(V-1)**4*(5
     8   *V**2-3*V+3)))+2*N**2*(V-1)**2*VC*W*(V**3*W**2-(V-3)*V**
     9   2*W+(V-2)*(V-1)**2)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-
     :   10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V
     ;   -1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V
     <   -1)**5)-2*(V-1)**2*VC*W*(V**3*W**2+V**2*W+(V-1)*(V**2-2*
     =   V+3))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5
     >   *V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W*
     ?   *3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))*LOG(
     @   V*W-V+1)/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**5*(V*W-V+1)**
     1   5)
 
      LW = (-N**4*(V-1)*VC*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V
     1   *W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-
     2   4*(V-1)**3*V*W+(V-1)**4)*(2*(8*V**5*(V**2-V+1)*W**6-V**5
     3   *(32*V**2-51*V+43)*W**5+2*V**3*(28*V**4-56*V**3+47*V**2-
     4   4*V-1)*W**4-V**2*(56*V**5-117*V**4+67*V**3+45*V**2-38*V+
     5   3)*W**3+V*(40*V**6-94*V**5+54*V**4+55*V**3-74*V**2+29*V-
     6   6)*W**2-(V-1)*(16*V**6-16*V**5-30*V**4+73*V**3-45*V**2+7
     7   *V+3)*W+4*(V-1)**2*(2*V**2-3*V+2)*(2*V**2-2*V+1))-CQ*(V*
     8   *2+1)**2*W*(V*W-1)*(V*W-V+1)*(2*W**2-2*W+1))+(V-1)**2*VC
     9   *W*(2*(V**3*W**2+V**2*W+V-1)+CQ*(V-1)*(V**2+1)*(2*W**2-2
     :   *W+1))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+
     ;   5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W
     <   **3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)-2*N**
     =   2*(V-1)*VC*(V**4*(V+7)*W**4-V**3*(V+3)**2*W**3+CQ*(V**2+
     >   1)*(V**2-V+1)*W*(V*W-1)*(2*W**2-2*W+1)+V*(8*V**3+V**2-V+
     ?   4)*W**2-(10*V**4-15*V**3+16*V**2-11*V+4)*W+4*(V-1)*(2*V*
     @   *2-2*V+1))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(
     1   V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-
     2   1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))*LOG(W)/(N**2*
     3   (V-1)**3*V**2*W**2*(V*W-1)**5*(V*W-V+1)**5)
 
      CVC = (-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1
     1   )*V*W+(V-1)**2)*(32*V**11*W**11-32*V**10*(5*V-1)*W**10+V
     2   **6*(416*V**5-245*V**4+59*V**3-78*V**2+52*V-12)*W**9-V**
     3   6*(704*V**5-641*V**4+96*V**3-177*V**2+192*V-46)*W**8+V**
     4   4*(832*V**7-906*V**6-165*V**5+228*V**4-123*V**3+314*V**2
     5   -192*V+36)*W**7-V**4*(704*V**7-714*V**6-964*V**5+1507*V*
     6   *4-1102*V**3+965*V**2-522*V+102)*W**6+V**2*(416*V**9-161
     7   *V**8-1817*V**7+2436*V**6-1023*V**5+48*V**4+341*V**3-456
     8   *V**2+228*V-36)*W**5-V**2*(160*V**9+235*V**8-1864*V**7+2
     9   185*V**6-180*V**5-1327*V**4+1396*V**3-943*V**2+396*V-66)
     :   *W**4+(V-1)*(32*V**10+288*V**9-793*V**8+13*V**7+1340*V**
     ;   6-1476*V**5+683*V**4-49*V**3-134*V**2+76*V-12)*W**3-(V-1
     <   )**2*(96*V**8-667*V**6+992*V**5-496*V**4-62*V**3+125*V**
     =   2-46*V+10)*W**2+(V-1)**3*(96*V**6-224*V**5+53*V**4+220*V
     >   **3-233*V**2+66*V-10)*W-8*(V-1)**4*(4*V**4-12*V**3+12*V*
     ?   *2-4*V-1))-2*N**2*(V-1)*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)*
     @   *2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*((V-1)*V**6*(3*V**2-
     1   6*V+10)*W**8-(V-1)*V**5*(10*V**3-24*V**2+41*V-10)*W**7+V
     2   **4*(15*V**5-60*V**4+125*V**3-77*V**2-19*V+20)*W**6-V**3
     3   *(15*V**6-80*V**5+198*V**4-119*V**3-80*V**2+94*V-20)*W**
     4   5+V**2*(10*V**7-55*V**6+143*V**5-52*V**4-145*V**3+123*V*
     5   *2-10*V-10)*W**4-(V-1)*V*(3*V**7-15*V**6+62*V**5-17*V**4
     6   -30*V**3-2*V**2+25*V-10)*W**3+(V-1)**2*V*(23*V**4-27*V**
     7   3+54*V**2-39*V+13)*W**2-(V-1)**3*(7*V**4-18*V**3+39*V**2
     8   -20*V+8)*W+4*(V-1)**4)-2*AL*N**4*(V-1)*VC*(V**2*W**2-2*V
     9   *W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3
     :   -10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(V**3*(3
     ;   *V**4-V**3+4*V**2+2)*W**6-2*V**2*(2*V**5+4*V**4+6*V**2+V
     <   +3)*W**5+V*(3*V**6+10*V**5+11*V**4+3*V**3+13*V**2+6*V+6)
     =   *W**4-(11*V**6+6*V**5+12*V**4+7*V**2+6*V+2)*W**3+(19*V**
     >   5-14*V**4+18*V**3-8*V**2+3*V+2)*W**2-(17*V**4-24*V**3+18
     ?   *V**2-8*V+1)*W+4*(V-1)*(2*V**2-2*V+1))+4*AL*N**2*(V-1)*V
     @   C*(V**2*W**2-2*V*W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V
     1   -1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V
     2   -1)**5)*(V**3*(3*V**4-3*V**3+4*V**2-2*V+2)*W**6-V**2*(4*
     3   V**5+5*V**4-5*V**3+10*V**2-4*V+6)*W**5+V*(3*V**6+7*V**5+
     4   4*V**4-V**3+7*V**2+6)*W**4-(9*V**6+V**5+5*V**4+V**2+4*V+
     5   2)*W**3+(13*V**5-13*V**4+13*V**3-6*V**2+V+2)*W**2-(9*V**
     6   4-13*V**3+10*V**2-5*V+1)*W+2*(V-1)*(2*V**2-2*V+1))+(V-1)
     7   **2*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**
     8   3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**4*(V
     9   **3+2*V**2-12*V+12)*W**6-V**3*(V+1)*(V**3+6*V**2-26*V+24
     :   )*W**5+V**3*(V**4+11*V**3+6*V**2-56*V+56)*W**4-V*(V**6+5
     ;   *V**5+20*V**4-22*V**3-22*V**2+60*V-24)*W**3+(5*V**6-23*V
     <   **5+78*V**4-85*V**3+40*V**2+4*V-12)*W**2+(V-1)*(5*V**4-8
     =   *V**3-29*V**2+24*V-10)*W+(V-1)**2*(3*V**2-2*V+6))-2*AL*(
     >   V-1)**2*VC*W*(V**2*W**2-2*V*W+1)*(V**3*(3*V**3-2*V**2+2*
     ?   V-2)*W**5-2*V**2*(2*V**4+3*V**3-2*V**2+2*V-3)*W**4+V*(3*
     @   V**5+7*V**4+4*V**3-V**2-6)*W**3-(7*V**5+3*V**4+V**3+V**2
     1   -4*V-2)*W**2+(7*V**4-5*V**3+3*V**2-V-2)*W-(V-1)*(V**2+1)
     2   )*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*
     3   (V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))/(N**2*(V-1)
     4   **3*V**2*W**2*(V*W-1)**5*(V*W-V+1)**5)/2.D0
 
      LM = LOG(S/M**2)*(N**4*VC*(V**3*W**3-3*(V-1)*V**2*W**2+3*(
     1   V-1)**2*V*W-(V-1)**3)*(16*V**8*W**8-16*V**7*(V+4)*W**7+V
     2   **3*(48*V**5-45*V**4+207*V**3-44*V**2+20*V-2)*W**6-2*V**
     3   2*(8*V**6+50*V**5-76*V**4+176*V**3-56*V**2+29*V-3)*W**5+
     4   V*(16*V**7-13*V**6+138*V**5-149*V**4+263*V**3-71*V**2+54
     5   *V-6)*W**4-(16*V**8-37*V**6+134*V**5-72*V**4+64*V**3+19*
     6   V**2+14*V-2)*W**3+(48*V**7-96*V**6+83*V**5+6*V**4-10*V**
     7   3-4*V**2+23*V-2)*W**2-(48*V**6-128*V**5+165*V**4-112*V**
     8   3+42*V**2-8*V+1)*W+8*(V-1)*(V**2-V+1)*(2*V**2-2*V+1))-2*
     9   N**2*VC*(V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*
     :   V*W-(V-1)**3)*(V**3*(3*V**3+5*V**2-4*V+4)*W**5-2*V**2*(2
     ;   *V**4+5*V**3+3*V**2-2*V+4)*W**4+V*(3*V**5+11*V**4+2*V**3
     <   +3*V**2+5*V+4)*W**3-V*(14*V**4-10*V**3+9*V**2-3*V+6)*W**
     =   2+(3*V**2-3*V+1)*(5*V**2-2*V+1)*W-4*(V-1)*(2*V**2-2*V+1)
     >   )+(V-1)*VC*W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W
     ?   -(V-1)**3)*(V**3*(3*V**3-2*V**2+2*V-2)*W**5-2*V**2*(2*V*
     @   *4+3*V**3-2*V**2+2*V-3)*W**4+V*(3*V**5+7*V**4+4*V**3-V**
     1   2-6)*W**3-(7*V**5+3*V**4+V**3+V**2-4*V-2)*W**2+(7*V**4-5
     2   *V**3+3*V**2-V-2)*W-(V-1)*(V**2+1)))/(N**2*(V-1)**2*V**2
     3   *W**2*(V*W-1)**3*(V*W-V+1)**3)
 
      LMP = LOG(S/MP**2)*(N**4*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1
     1   )*(16*V**8*W**8-16*V**7*(5*V-4)*W**7+V**6*(208*V**2-351*
     2   V+151)*W**6-(V-1)*V**5*(352*V**2-559*V+230)*W**5+2*(V-1)
     3   **2*V**3*(204*V**3-307*V**2+123*V+1)*W**4-2*(V-1)**3*V**
     4   2*(168*V**3-235*V**2+90*V+3)*W**3+(V-1)**4*V*(200*V**3-2
     5   55*V**2+95*V+6)*W**2-(V-1)**5*(80*V**3-95*V**2+38*V+2)*W
     6   +8*(V-1)**6*(2*V**2-2*V+1))-2*N**2*(V-1)*VC*W*(V**3*W**3
     7   -3*V**2*W**2+3*V*W-1)*(V**6*W**5-V**5*(2*V-3)*W**4-(V-1)
     8   *V**3*(4*V**3-7*V**2+10*V-2)*W**3+(V-1)**2*V**2*(8*V**3-
     9   15*V**2+17*V-6)*W**2-(V-1)**3*V*(4*V**3-10*V**2+13*V-6)*
     :   W-(V-1)**4*(V**2-2*V+2))+(V-1)*VC*W*(V**2*W**2-2*(V-1)*V
     ;   *W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**3-
     <   (V-2)*V**3*W**2+(V-1)*V*(V**2-5*V+2)*W-(V-1)**2*(V**2-2*
     =   V+2)))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)**3*(V*W-V+1)**3)
 
      STRUV16=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
      RETURN
      END
 
C---CHANGEMENT DU SHEMA DE FACTORISATION-------------------------------
C    THIS PROGRAM ALLOWS TO MOVE FROM A FACTORIZATION SCHEME TO ANOTHER
      DOUBLE PRECISION FUNCTION HQQD(V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      UN=1.D0
      IF(J.EQ.1) THEN
      HQQD=CQQD(UN)
      ELSE
      HQQD=(1.-V)/V*(CQQD(UN)+DLOG(V/(1.-V))*CQQW(UN)+DLOG((1.-V)/V)**2/
     12.*CQQL(UN))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HGQD(V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      UN=1.D0
      IF(J.EQ.1) THEN
      HGQD=CGQD(UN)
      ELSE
      HGQD=(1.-V)/V*(CGQD(UN)+DLOG(V/(1.-V))*CGQW(UN)+DLOG((1.-V)/V)**2/
     12.*CGQL(UN))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HQGD(V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      UN=1.D0
      IF(J.EQ.1) THEN
      HQGD=CQGD(UN)
      ELSE
      HQGD=(1.-V)/V*(CQGD(UN)+DLOG(V/(1.-V))*CQGW(UN)+DLOG((1.-V)/V)**2/
     12.*CQGL(UN))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HGGD(V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      UN=1.D0
      IF(J.EQ.1) THEN
      HGGD=CGGD(UN)
      ELSE
      HGGD=(1.-V)/V*(CGGD(UN)+DLOG(V/(1.-V))*CGGW(UN)+DLOG((1.-V)/V)**2/
     12.*CGGL(UN))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HQQW(W,V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      IF(J.EQ.1) THEN
      HQQW=CQQW(W)
      ELSE
      HQQW=(1.-V*W)/V*(CQQW((1.-V)/(1.-V*W))+DLOG(V/(1.-V*W))*CQQL((1.
     1-V)/(1.-V*W)))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HQGW(W,V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      IF(J.EQ.1) THEN
      HQGW=CQGW(W)
      ELSE
      HQGW=(1.-V*W)/V*(CQGW((1.-V)/(1.-V*W))+DLOG(V/(1.-V*W))*CQGL((1.
     1-V)/(1.-V*W)))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HGGW(W,V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      IF(J.EQ.1) THEN
      HGGW=CGGW(W)
      ELSE
      HGGW=(1.-V*W)/V*(CGGW((1.-V)/(1.-V*W))+DLOG(V/(1.-V*W))*CGGL((1.
     1-V)/(1.-V*W)))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HGQW(W,V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      IF(J.EQ.1) THEN
      HGQW=CGQW(W)
      ELSE
      HGQW=(1.-V*W)/V*(CGQW((1.-V)/(1.-V*W))+DLOG(V/(1.-V*W))*CGQL((1.
     1-V)/(1.-V*W)))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HGQL(W,V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      IF(J.EQ.1) THEN
      HGQL=CGQL(W)
      ELSE
      HGQL=(1.-V*W)/V*CGQL((1.-V)/(1.-V*W))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HGGL(W,V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      IF(J.EQ.1) THEN
      HGGL=CGGL(W)
      ELSE
      HGGL=(1.-V*W)/V*CGGL((1.-V)/(1.-V*W))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HQQL(W,V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      IF(J.EQ.1) THEN
      HQQL=CQQL(W)
      ELSE
      HQQL=(1.-V*W)/V*CQQL((1.-V)/(1.-V*W))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION HQGL(W,V,J)
      IMPLICIT REAL*8 (A-H,L-Z)
      IF(J.EQ.1) THEN
      HQGL=CQGL(W)
      ELSE
      HQGL=(1.-V*W)/V*CQGL((1.-V)/(1.-V*W))
      ENDIF
      RETURN
      END
C    POUR LA FRAGMENTATION
      DOUBLE PRECISION FUNCTION HFQQL(W,V)
      IMPLICIT REAL*8 (A-H,L-Z)
      HFQQL=FQQL(1.-V+V*W)/V
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFQQW(W,V)
      IMPLICIT REAL*8 (A-H,L-Z)
      HFQQW=1./V*(FQQW(1.-V+V*W)+DLOG(V)*FQQL(1.-V+V*W))
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFQQD(V)
      IMPLICIT REAL*8 (A-H,L-Z)
      UN=1.D0
      HFQQD=1./V*(FQQD(UN)+DLOG(V)*FQQW(UN)+DLOG(V)**2/2.*FQQL(UN))
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFQGL(W,V)
      IMPLICIT REAL*8 (A-H,L-Z)
      HFQGL=FQGL(1.-V+V*W)/V
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFQGW(W,V)
      IMPLICIT REAL*8 (A-H,L-Z)
      HFQGW=1./V*(FQGW(1.-V+V*W)+DLOG(V)*FQGL(1.-V+V*W))
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFQGD(V)
      IMPLICIT REAL*8 (A-H,L-Z)
      UN=1.D0
      HFQGD=1./V*(FQGD(UN)+DLOG(V)*FQGW(UN)+DLOG(V)**2/2.*FQGL(UN))
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFGGL(W,V)
      IMPLICIT REAL*8 (A-H,L-Z)
      HFGGL=FGGL(1.-V+V*W)/V
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFGGW(W,V)
      IMPLICIT REAL*8 (A-H,L-Z)
      HFGGW=1./V*(FGGW(1.-V+V*W)+DLOG(V)*FGGL(1.-V+V*W))
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFGGD(V)
      IMPLICIT REAL*8 (A-H,L-Z)
      UN=1.D0
      HFGGD=1./V*(FGGD(UN)+DLOG(V)*FGGW(UN)+DLOG(V)**2/2.*FGGL(UN))
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFGQL(W,V)
      IMPLICIT REAL*8 (A-H,L-Z)
      HFGQL=FGQL(1.-V+V*W)/V
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFGQW(W,V)
      IMPLICIT REAL*8 (A-H,L-Z)
      HFGQW=1./V*(FGQW(1.-V+V*W)+DLOG(V)*FGQL(1.-V+V*W))
      RETURN
      END
      DOUBLE PRECISION FUNCTION HFGQD(V)
      IMPLICIT REAL*8 (A-H,L-Z)
      UN=1.D0
      HFGQD=1./V*(FGQD(UN)+DLOG(V)*FGQW(UN)+DLOG(V)**2/2.*FGQL(UN))
      RETURN
      END
      DOUBLE PRECISION FUNCTION A(S,T,U)
      IMPLICIT REAL*8 (A-H,L-Z)
      N=3.D0
      VC=(N**2-1.)
      A=2.*VC*(S**2+U**2)/T**2
      RETURN
      END
      DOUBLE PRECISION FUNCTION B(S,T,U)
      IMPLICIT REAL*8 (A-H,L-Z)
      N=3.D0
      VC=(N**2-1.)
      B=-4.*VC/N*S**2/U/T
      RETURN
      END
      DOUBLE PRECISION FUNCTION C(S,T,U)
      IMPLICIT REAL*8 (A-H,L-Z)
      N=3.D0
      VC=(N**2-1.)
      C=2.*VC/N*(VC/U/T-2.*N**2/S**2)*(T**2+U**2)
      RETURN
      END
      DOUBLE PRECISION FUNCTION D(S,T,U)
      IMPLICIT REAL*8 (A-H,L-Z)
      N=3.D0
      VC=(N**2-1.)
      D=16.*VC*N**2*(3.-U*T/S**2-U*S/T**2-S*T/U**2)
      RETURN
      END
      DOUBLE PRECISION FUNCTION A0(X,S)
      IMPLICIT REAL*8 (A-H,L-Z)
C     QJ+QK-->QJ+QK
      T=-S*(1.-X)
      U=-S*X
      A0=A(S,T,U)
      RETURN
      END
      DOUBLE PRECISION FUNCTION A2(X,S)
      IMPLICIT REAL*8 (A-H,L-Z)
C     Q+QB --> Q+QB DIFFERENT FLAVORS
      T=-S*(1.-X)
      U=-S*X
      A2=A(T,S,U)
      RETURN
      END
      DOUBLE PRECISION FUNCTION B0(X,S)
      IMPLICIT REAL*8 (A-H,L-Z)
C     Q+Q --> Q+Q SAME FLAVOR
      T=-S*(1.-X)
      U=-S*X
      B0=A(S,T,U)+A(S,U,T)+B(S,T,U)
      RETURN
      END
      DOUBLE PRECISION FUNCTION D0(X,S)
      IMPLICIT REAL*8 (A-H,L-Z)
C     Q+QB ->Q+QB SAME FLAVOR
      T=-S*(1.-X)
      U=-S*X
      D0=A(U,T,S)+A(U,S,T)+B(U,T,S)
      RETURN
      END
      DOUBLE PRECISION FUNCTION D1(X,S)
      IMPLICIT REAL*8 (A-H,L-Z)
C     Q+QB -->G+G    AND ALSO G+G --> Q+QB
      T=-S*(1.-X)
      U=-S*X
      D1=C(S,T,U)
      RETURN
      END
      DOUBLE PRECISION FUNCTION E0(X,S)
      IMPLICIT REAL*8 (A-H,L-Z)
C     Q+G -->Q+G
      T=-S*(1.-X)
      U=-S*X
      E0=-C(T,S,U)
      RETURN
      END
      DOUBLE PRECISION FUNCTION F2(X,S)
      IMPLICIT REAL*8 (A-H,L-Z)
C     G+G -->G+G
      T=-S*(1.-X)
      U=-S*X
      F2=D(S,T,U)
      RETURN
      END
 
