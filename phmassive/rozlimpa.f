      FUNCTION HQHLpaM0(TX,T1,RO)
      IMPLICIT DOUBLE PRECISION (A-Z)
      T2 = -TX-T1+1
      PP = (-4.0D0)*(-T2-T1+1)*(2*T2**2+3*T1*T2+T1**2-T1)/(3.0D0*T1*(-T2
     1   -T1))
      PP = 4.0D0*(-T2-T1+1)*T2*(T2+T1)**2/(3.0D0*T1*(-T2-T1)**3)+PP
      PP = 4.0D0*(-T2-T1+1)*(T2-1)*(2*T2+T1-1)/(3.0D0*T1*(1-T2))+PP
      PP = 8.0D0*(-T2-T1+1)*(3*T1*T2+3*T2+T1**2-2)/(3.0D0*(T1-1)*T1)+PP
      PP = (-8.0D0)*(-T2-T1+1)*(T1*T2+T2-1)/(3.0D0*(T1-1)*T1)+PP
      PP = 1.6D1*(T2+T1-1)*(4*T2**2+6*T1*T2-4*T2+3*T1**2-4*T1+2)/(3.0D0*
     1   (T1-1)*T1)+PP
      PP = PP
      HQHLpaM0 = PP
      RETURN 
      END 
      FUNCTION HQHPpaM0(TX,T1,RO)
      IMPLICIT DOUBLE PRECISION (A-Z)
      XLRO = LOG(RO/4.0D0)
      FOUR = 4
      XLOG4 = LOG(FOUR)
      T2 = -TX-T1+1
      PP = XLRO+XLOG4-LOG(-4*(T1-1)*(1-T2))
      PP = (-8.0D0)*PP*(-T2-T1+1)*(4*T2**2+4*T1*T2-6*T2+2*T1**2-4*T1+3)/
     1   (3.0D0*(T1-1)*T1)
      PP = 4.0D0*(-T2-T1+1)*(2*T2**2+3*T1*T2+T1**2-T1)*(XLRO+XLOG4-LOG(4
     1   *((-T2-T1+1)**2-2*(-T2-T1+1)+1)))/(3.0D0*T1*(-T2-T1))+PP
      PP = (-4.0D0)*(-T2-T1+1)*T2*(T2+T1)**2*(XLRO+XLOG4-LOG(4*((-T2-T1+
     1   1)**2-2*(-T2-T1+1)+1)))/(3.0D0*T1*(-T2-T1)**3)+PP
      PP = (-4.0D0)*(-T2-T1+1)*(T2-1)*(2*T2+T1-1)*XLRO/(3.0D0*T1*(1-T2))
     1   +PP
      PP = (-8.0D0)*(-T2-T1+1)*(3*T1*T2+3*T2+T1**2-2)*XLRO/(3.0D0*(T1-1)
     1   *T1)+PP
      PP = 8.0D0*(-T2-T1+1)*(T1*T2+T2-1)*XLRO/(3.0D0*(T1-1)*T1)+PP
      TMP0 = -2*XLRO-2*XLOG4+LOG(16*T1**2/(T1**2-2*T1+1))
      TMP0 = 4.0D0*(T2+T1-1)*(4*T2**2+6*T1*T2-4*T2+3*T1**2-4*T1+2)*TMP0/
     1   (3.0D0*(T1-1)*T1)
      PP = TMP0+PP
      PP = 8.0D0*(-T2-T1+1)*(8*T2**2+6*T1*T2-6*T2+3*T1**2-4*T1+3)*LOG(-T
     1   2/(T1-1))/(3.0D0*(T1-1)*T1)+PP
      PP = 8.0D0*LOG(T1/(1-T2))*(-T2-T1+1)*(4*T2**2+2*T1*T2-4*T2+T1**2-2
     1   *T1+2)/(3.0D0*(T1-1)*T1)+PP
      PP = 8.0D0*(-T2-T1+1)*(2*T2+T1-1)/(3.0D0*(T1-1)*T1)+PP
      PP = (-4.0D0)*(T2+T1-1)*(T2+T1)*(T1*T2+T2+T1**2-T1)/(3.0D0*T1*(-T2
     1   -T1)**2)+PP
      PP = (-2.0D0)*(T2+T1-1)*(4*T2**2+3*T1*T2-5*T2+T1+1)/(3.0D0*T1*(1-T
     1   2))+PP
      PP = (-8.0D0)*(T1+1)*(T2+T1-1)*(T1*T2+T2+T1-1)/(3.0D0*(T1-1)**2*T1
     1   )+PP
      PP = 2.0D0*(T2+T1-1)*(8*T2-T1**2+4*T1-3)/(3.0D0*(T1-1)*T1)+PP
      HQHPpaM0 = PP
      RETURN 
      END 
