      FUNCTION HQHDpnM0(T1,RO,NL)
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER NL
      DATA PI/3.141 592 653 589 793/
      DD = 0
      HQHDpnM0 = DD
      RETURN
      END
      FUNCTION HQHLpnM0(TX,T1,RO)
      IMPLICIT DOUBLE PRECISION (A-Z)
      T2 = -TX-T1+1
      T22 = 1/(TX+T1)
      PP = 4.0D0*(-T2-T1+1)**2/(3.0D0*T1*(1-T2)**2)
      PP = 4.0D0*(-T2-T1+1)*(2*T2+T1-1)/(3.0D0*T1*(1-T2))+PP
      PP = (-8.0D0)*(-T2-T1+1)*(2*T2-1)/(3.0D0*(T1-1)*T1)+PP
      PP = 8.0D0*(-T2-T1+1)*T2/(3.0D0*(T1-1)**2)+PP
      TMP0 = -T1**2*T22+2*T2+2*T1-2
      TMP0 = (-8.0D0)*(T2-1)*(T2+T1-1)*(2*T2**2-2*T2+1)*T22**2*TMP0/(3.0
     1   D0*T1**2*T2)
      PP = TMP0+PP
      HQHLpnM0 = PP
      RETURN 
      END 
      FUNCTION HQHPpnM0(TX,T1,RO)
      IMPLICIT DOUBLE PRECISION (A-Z)
      XLRO = LOG(RO/4.0D0)
      T2 = -TX-T1+1
      T22 = 1/(TX+T1)
      PP = (-4.0D0)*(-T2-T1+1)**2*XLRO/(3.0D0*T1*(1-T2)**2)
      PP = (-4.0D0)*(-T2-T1+1)*(2*T2+T1-1)*XLRO/(3.0D0*T1*(1-T2))+PP
      PP = 8.0D0*(-T2-T1+1)*(2*T2-1)*XLRO/(3.0D0*(T1-1)*T1)+PP
      PP = (-8.0D0)*(-T2-T1+1)*T2*XLRO/(3.0D0*(T1-1)**2)+PP
      TMP0 = -T1**2*T22+2*T2+2*T1-2
      TMP0 = (-8.0D0)*(-T2-T1+1)*(T2-1)*(2*T2**2-2*T2+1)*T22**2*TMP0*XLR
     1   O/(3.0D0*T1**2*T2)
      PP = TMP0+PP
      PP = 8.0D0*(-T2-T1+1)*(4*T2**2+2*T1*T2-4*T2+T1**2-2*T1+2)*LOG(-T2/
     1   (T1-1))/(3.0D0*T1**2*T2)+PP
      PP = (-1.6D1)*(-T2-T1+1)**2/(3.0D0*(T1-1)*T1)+PP
      PP = 8.0D0*LOG(T1/(1-T2))*(-T2-T1+1)*(4*T2**2+2*T1*T2-4*T2+T1**2-2
     1   *T1+2)/(3.0D0*T1**2*T2)+PP
      PP = (-4.0D0)*(T2+T1-1)*(2*T2+T1-1)/(3.0D0*T1*(1-T2))+PP
      PP = (-4.0D0)*(T2+T1-1)*(4*T2**3+4*T1*T2**2-12*T2**2+T1**2*T2-8*T1
     1   *T2+12*T2-2*T1**2+4*T1-4)/(3.0D0*T1**2*(1-T2)**2)+PP
      PP = 8.0D0*(T2+T1-1)*(2*T1*T2+T2+T1**2-1)/(3.0D0*(T1-1)**2*T1)+PP
      PP = 4.0D0*(T2+T1-1)*(8*T2**2+8*T1*T2-8*T2-T1**2)/(3.0D0*T1**2*T2)
     1   +PP
      HQHPpnM0 = PP
      RETURN 
      END 
