      FUNCTION HQHDpcM0(T1,RO,NL)
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER NL
      DATA PI/3.141 592 653 589 793/
      DD = 0
      HQHDpcM0 = DD
      RETURN
      END
      FUNCTION HQHLpcM0(TX,T1,RO)
      IMPLICIT DOUBLE PRECISION (A-Z)
      T2 = -TX-T1+1
      T11 = 1/(1-T1)
      PP = 8.0D0*(-T2-T1+1)**2/(3.0D0*(T1-1)**2)
      PP = 8.0D0*(-T2-T1+1)*T2/(3.0D0*(-T2-T1))+PP
      PP = (-4.0D0)*(-T2-T1+1)*T2*(T2+T1)/(3.0D0*(-T2-T1)**3)+PP
      PP = (-8.0D0)*(-T2-T1+1)*(2*T2+T1-1)/(3.0D0*(T1-1))+PP
      PP = (-8.0D0)*(2*T1**2-2*T1+1)*(T2+T1-1)*(2*T11**2*T2**2-2*T11*T2+
     1   1)/(3.0D0*T2)+PP
      HQHLpcM0 = PP
      RETURN 
      END 
      FUNCTION HQHPpcM0(TX,T1,RO)
      IMPLICIT DOUBLE PRECISION (A-Z)
      XLRO = LOG(RO/4.0D0)
      FOUR = 4
      XLOG4 = LOG(FOUR)
      T2 = -TX-T1+1
      T11 = 1/(1-T1)
      PP = (-8.0D0)*(-T2-T1+1)*T2*(XLRO+XLOG4-LOG(4*((-T2-T1+1)**2-2*(-T
     1   2-T1+1)+1)))/(3.0D0*(-T2-T1))
      PP = 4.0D0*(-T2-T1+1)*T2*(T2+T1)*(XLRO+XLOG4-LOG(4*((-T2-T1+1)**2-
     1   2*(-T2-T1+1)+1)))/(3.0D0*(-T2-T1)**3)+PP
      PP = (-8.0D0)*(-T2-T1+1)**2*XLRO/(3.0D0*(T1-1)**2)+PP
      PP = 8.0D0*(-T2-T1+1)*(2*T2+T1-1)*XLRO/(3.0D0*(T1-1))+PP
      PP = (-8.0D0)*(2*T1**2-2*T1+1)*(-T2-T1+1)*(2*T11**2*T2**2-2*T11*T2
     1   +1)*XLRO/(3.0D0*T2)+PP
      PP = 8.0D0*(-T2-T1+1)*(4*T2**2+4*T1*T2-2*T2+2*T1**2-2*T1+1)*LOG(-T
     1   2/(T1-1))/(3.0D0*T2)+PP
      PP = (-8.0D0)*(-T2-T1+1)*(3*T2**2+4*T1*T2-2*T2+T1**2-T1)/(3.0D0*T2
     1   )+PP
      PP = (-8.0D0)*(-T2-T1+1)*(2*T2+T1)/(3.0D0*(T1-1))+PP
      PP = 4.0D0*(T2+T1-1)*(2*T2**3+4*T1*T2**2-2*T2**2+2*T1**2*T2-T1*T2+
     1   T2+T1**2-T1)/(3.0D0*(-T2-T1)**2)+PP
      PP = 8.0D0*(T2+T1-1)*(T1*T2+2*T2+T1-1)/(3.0D0*(T1-1)**2)+PP
      PP = PP
      HQHPpcM0 = PP
      RETURN 
      END 
