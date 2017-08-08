
C----------- LES COMBINAISONS DE FONCTION DE STRUCTURE
C*===================================================================
C
C                 STRUCTURE FUNCTION AT SCALE Q2
C
C==================================================================  */
c the lines commented with * were used by P.N. as an alternative
c (easier to read) implementation, that was tested for identity
c with the original one. A bug in the original implementation
c was then fixed. The lines have been left here, for ease in cases
c of future extensions.
      SUBROUTINE STRU(IFRAGMODE
     # ,XA,XB,XC,IH1,IH2,Q2,QP2,GPPV1,GPPV2,GPPV3,GPPV4
     # ,GPPV5,GPPV6,GPPV7,GPPV8,GPPV9,GPPV10,GPPV11,GPPV12,GPPV13
     # ,GPPV14,GPPV15,GPPV16)
c      IMPLICIT REAL*8(A-H,M-Z)
      real * 8
     #  XA,XB,XC,Q2,QP2,GPPV1,GPPV2,GPPV3,GPPV4
     # ,GPPV5,GPPV6,GPPV7,GPPV8,GPPV9,GPPV10,GPPV11,GPPV12,GPPV13
     # ,GPPV14,GPPV15,GPPV16
      integer IFRAGMODE,ih1,ih2
      integer ipion
      COMMON/PIONCM/IPION
      integer iasyflag
      common/ciasyflag/iasyflag
      integer ip,ipb,ipi
      real * 8 XUHA,XUBHA,XDHA,XDBHA,XSHA,
     # XCHA,XBHA,XTHA,XGPROA,XUHB,XUBHB,XDHB,XDBHB,XSHB,
     # XCHB,XBHB,XTHB,XGPROB,XDUP,XDUBP,XDDP,XDDBP,XDSP,
     # XDCP,xdcbp,XDBP,XDBBP,XDTP,XDTBP,XDGP,GPPU1,GPPU2,GPPU3,GPPU4
     # ,GPPU5,GPPU6,GPPU7,GPPU8,GPPU9,GPPU10,GPPU11,GPPU12,GPPU13
     # ,GPPU14,GPPU15,GPPU16
*      real * 8 fa(-6:6),fb(-6:6),fd(-6:6),tmp
*      integer j,k
      IF (IPION.EQ.0) THEN
ccc        CALL FRAEVOL(XC,QP2)
      ENDIF
C  /* HADRON 0 =P, H 1 =PBAR */
      IP=IH1
      IPB=IH2
      IPI=IPION
c      CALL FONSTRU(XA,IP,Q2/xc**2,XUHA,XUBHA,XDHA,XDBHA,XSHA,
      CALL FONSTRU(XA,IP,Q2,XUHA,XUBHA,XDHA,XDBHA,XSHA,
     # XCHA,XBHA,XTHA,XGPROA)
c      CALL FONSTRU(XB,IPB,Q2/xc**2,XUHB,XUBHB,XDHB,XDBHB,XSHB,
      CALL FONSTRU(XB,IPB,Q2,XUHB,XUBHB,XDHB,XDBHB,XSHB,
     # XCHB,XBHB,XTHB,XGPROB)

c      CALL FONFRA(IFRAGMODE,XC,IPI,QP2/xc**2,XDUP,XDUBP,XDDP,XDDBP,XDSP,
      CALL FONFRA(IFRAGMODE,XC,IPI,QP2,XDUP,XDUBP,XDDP,XDDBP,XDSP,
     #  XDCP,xdcbp,XDBP,XDBBP,XDTP,XDTBP,XDGP)
      if(iasyflag.eq.0) then
         xdcp=(xdcp+xdcbp)/2
         xdcbp=xdcp
         xdbp=(xdbp+xdbbp)/2
         xdbbp=xdbp
      endif
c The following is incorrect; better leave it out
c      IF (XA.LE.1.D-5.OR.XB.LE.1.D-5) THEN
c        WRITE (15,*) 'XA OR XB OUT OF THEIR RANGE'
c      ELSE
*         fa(0)=xgproa
*         fa(1)=xuha
*         fa(-1)=xubha
*         fa(2)=xdha
*         fa(-2)=xdbha
*         fa(3)=xsha
*         fa(-3)=xsha
*         fa(4)=xcha
*         fa(-4)=xcha
*         fa(5)=xbha
*         fa(-5)=xbha
*         fa(6)=xtha
*         fa(-6)=xtha
* c
*         fb(0)=xgprob
*         fb(1)=xuhb
*         fb(-1)=xubhb
*         fb(2)=xdhb
*         fb(-2)=xdbhb
*         fb(3)=xshb
*         fb(-3)=xshb
*         fb(4)=xchb
*         fb(-4)=xchb
*         fb(5)=xbhb
*         fb(-5)=xbhb
*         fb(6)=xthb
*         fb(-6)=xthb
* c
*         fd(0)=xdgp
*         fd(1)=xdup
*         fd(-1)=xdubp
*         fd(2)=xddp
*         fd(-2)=xddbp
*         fd(3)=xdsp
*         fd(-3)=xdsp
*         fd(4)=xdcp
*         fd(-4)=xdcbp
*         fd(5)=xdbp
*         fd(-5)=xdbbp
*         fd(6)=xdtp
*         fd(-6)=xdtbp
c
c this seems q qprime -> q + qbar qbarprime -> qbar
*         tmp=0
*         do j=1,6
*            do k=1,6
*               if(k.ne.j) then
*                  tmp=tmp+fa(j)*fb(k)*fd(j)+fa(-j)*fb(-k)*fd(-j)
*               endif
*            enddo
*         enddo
      GPPU1=XUHA*(XDHB+XSHB+XCHB+XBHB)*XDUP
     &+XDHA*(XUHB+XSHB+XCHB+XBHB)*XDDP+
     &XSHA*(XUHB+XDHB+XCHB+XBHB)*XDSP+
     &XCHA*(XUHB+XDHB+XSHB+XBHB)*XDCP+
     &XBHA*(XUHB+XDHB+XSHB+XCHB)*XDBP+
     &(XUHA*XDUP+XDHA*XDDP+XSHA*XDSP+
     &XCHA*XDCP+XBHA*XDBP)*XTHB+XTHA*
     &(XUHB+XDHB+XSHB+XCHB+XBHB)*XDTP
     &+XUBHA*(XDBHB+XSHB+XCHB+XBHB)*XDUBP
     &+XDBHA*(XUBHB+XSHB+XCHB+XBHB)*XDDBP+
     &XSHA*(XUBHB+XDBHB+XCHB+XBHB)*XDSP+
     &XCHA*(XUBHB+XDBHB+XSHB+XBHB)*xdcbp+
     &XBHA*(XUBHB+XDBHB+XSHB+XCHB)*XDBBP+
     &(XUBHA*XDUBP+XDBHA*XDDBP+XSHA*XDSP+
     &XCHA*xdcbp+XBHA*XDBBP)*XTHB+XTHA*
     &(XUBHB+XDBHB+XSHB+XCHB+XBHB)*XDTBP
      GPPV1=GPPU1/XA/XB/XC**3
*       call chk(tmp,gppu1)
C*********************************************************************
c seems q qprime -> g + cc (qprime>q)
*         tmp=0
*         do j=1,6
*            do k=j+1,6
*               tmp=tmp+fa(j)*fb(k)*fd(0)+fa(-j)*fb(-k)*fd(0)
*            enddo
*         enddo
      GPPU2=(XUHA*(XDHB+XSHB+XCHB+XBHB+XTHB)
     &+XDHA*(XSHB+XCHB+XBHB+XTHB)+
     &XSHA*(XCHB+XBHB+XTHB)+
     &XCHA*(XBHB+XTHB)+
     &XBHA*(XTHB)
     &+XUBHA*(XDBHB+XSHB+XCHB+XBHB+XTHB)
     &+XDBHA*(XSHB+XCHB+XBHB+XTHB)+
     &XSHA*(XCHB+XBHB+XTHB)+
     &XCHA*(XBHB+XTHB)+
     &XBHA*(XTHB))*XDGP
      GPPV2=GPPU2/XA/XB/XC**3
*       call chk(tmp,gppu2)
C*********************************************************************
c This seems q + qbarprime -> q
*         tmp=0
*         do j=1,6
*            do k=1,6
*               if(k.ne.j) then
*                  tmp=tmp+fa(j)*fb(-k)*fd(j)+fa(-j)*fb(k)*fd(-j)
*               endif
*            enddo
*         enddo
      GPPU3=XUHA*(XDBHB+XSHB+XCHB+XBHB)*XDUP
     &+XDHA*(XUBHB+XSHB+XCHB+XBHB)*XDDP+
     &XSHA*(XUBHB+XDBHB+XCHB+XBHB)*XDSP+
     &XCHA*(XUBHB+XDBHB+XSHB+XBHB)*XDCP+
     &XBHA*(XUBHB+XDBHB+XSHB+XCHB)*XDBP+
     &(XUHA*XDUP+XDHA*XDDP+XSHA*XDSP+
     &XCHA*XDCP+XBHA*XDBP)*XTHB+XTHA*
     &(XUBHB+XDBHB+XSHB+XCHB+XBHB)*XDTP
     &+XUBHA*(XDHB+XSHB+XCHB+XBHB)*XDUBP
     &+XDBHA*(XUHB+XSHB+XCHB+XBHB)*XDDBP+
     &XSHA*(XUHB+XDHB+XCHB+XBHB)*XDSP+
     &XCHA*(XUHB+XDHB+XSHB+XBHB)*xdcbp+
     &XBHA*(XUHB+XDHB+XSHB+XCHB)*XDBBP+
     &(XUBHA*XDUBP+XDBHA*XDDBP+XSHA*XDSP+
     &XCHA*xdcbp+XBHA*XDBBP)*XTHB+XTHA*
     &(XUHB+XDHB+XSHB+XCHB+XBHB)*XDTBP
      GPPV3=GPPU3/XA/XB/XC**3
*       call chk(tmp,gppu3)
C*********************************************************************
c q + qbarprime -> g
*         tmp=0
*         do j=1,6
*            do k=1,6
*               if(k.ne.j) then
*                  tmp=tmp+fa(j)*fb(-k)*fd(0)
*               endif
*            enddo
*         enddo
      GPPU4=(XUHA*(XDBHB+XSHB+XCHB+XBHB+XTHB)
     &+XDHA*(XUBHB+XSHB+XCHB+XBHB+XTHB)+
     &XSHA*(XUBHB+XDBHB+XCHB+XBHB+XTHB)+
     &XCHA*(XUBHB+XDBHB+XSHB+XBHB+XTHB)+
     &XBHA*(XUBHB+XDBHB+XSHB+XCHB+XTHB)+
     &XTHA*(XUBHB+XDBHB+XSHB+XCHB+XBHB))
     &*XDGP
      GPPV4=GPPU4/XA/XB/XC**3
*       call chk(tmp,gppu4)
C*********************************************************************
c q qbar -> qprime + cc
*         tmp=0
*         do j=1,6
*            do k=1,6
*               if(k.ne.j) then
*                  tmp=tmp+fa(j)*fb(-j)*fd(k)+fa(-j)*fb(j)*fd(-k)
*               endif
*            enddo
*         enddo
      GPPU5=(XDHA*XDBHB+XSHA*XSHB+XCHA*
     &XCHB+XBHA*XBHB+XTHA*XTHB)*XDUP+(XUHA
     &*XUBHB+XSHA*XSHB+XBHA*XBHB+XCHA*XCHB
     &+XTHA*XTHB)*XDDP+(XUHA*XUBHB+XDHA
     &*XDBHB+XCHA*XCHB+XBHA*XBHB+XTHA*XTHB)
     &*XDSP+(XUHA*XUBHB+XDHA*XDBHB+
     &XSHA*XSHB+XBHA*XBHB+XTHA*XTHB)*XDCP+
     &(XUHA*XUBHB+XDHA*XDBHB+XSHA*XSHB
     &+XCHA*XCHB+XTHA*XTHB)*XDBP+(XUHA*
     &XUBHB+XDHA*XDBHB+XSHA*XSHB+XCHA*
     &XCHB+XBHA*XBHB)*XDTP+
     &(XDBHA*XDHB+XSHA*XSHB+XCHA*
     &XCHB+XBHA*XBHB+XTHA*XTHB)*XDUBP+(XUBHA
     &*XUHB+XSHA*XSHB+XBHA*XBHB+XCHA*XCHB
     &+XTHA*XTHB)*XDDBP+(XUBHA*XUHB+XDBHA
     &*XDHB+XCHA*XCHB+XBHA*XBHB+XTHA*XTHB)
     &*XDSP+(XUBHA*XUHB+XDBHA*XDHB+
     &XSHA*XSHB+XBHA*XBHB+XTHA*XTHB)*xdcbp+
     &(XUBHA*XUHB+XDBHA*XDHB+XSHA*XSHB
     &+XCHA*XCHB+XTHA*XTHB)*XDBBP+(XUBHA*
     &XUHB+XDBHA*XDHB+XSHA*XSHB+XCHA*
     &XCHB+XBHA*XBHB)*XDTBP
      GPPV5=GPPU5/XA/XB/XC**3
*       call chk(tmp,gppu5)
C*********************************************************************
c q q -> q  + cc ( 1/2)
*         tmp=0
*         do j=1,6
*            tmp=tmp+fa(j)*fb(j)*fd(j)+fa(-j)*fb(-j)*fd(-j)
*         enddo
      GPPU6=XUHA*XUHB*XDUP+XDHA*XDHB
     &*XDDP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP
     &+XUBHA*XUBHB*XDUBP+XDBHA*XDBHB
     &*XDDBP+XSHA*XSHB*XDSP+XCHA*XCHB*xdcbp
     &+XBHA*XBHB*XDBP+XBHA*XBHB*XDBBP+
     &XTHA*XTHB*XDTP+XTHA*XTHB*XDTBP
      GPPV6=GPPU6/XA/XB/2.D0/XC**3
*       call chk(tmp,gppu6)
C*********************************************************************
c q q -> g + cc (1/2)
*         tmp=0
*         do j=1,6
*            tmp=tmp+fa(j)*fb(j)*fd(0)+fa(-j)*fb(-j)*fd(0)
*         enddo
      GPPU7=(XUHA*XUHB+XDHA*XDHB
     &+XSHA*XSHB+XCHA*XCHB+XBHA*XBHB*2.
     &+XUBHA*XUBHB+XDBHA*XDBHB
     &+XSHA*XSHB+XCHA*XCHB+2.*XTHA*XTHB)*XDGP
      GPPV7=GPPU7/XA/XB/2.D0/XC**3
*       call chk(tmp,gppu7)
C*********************************************************************
c q  g  -> qprime
*         tmp=0
*         do j=1,6
*            do k=1,6
*               if(k.ne.j) then
*                  tmp=tmp+fa(j)*fb(0)*fd(k)+fa(-j)*fb(0)*fd(-k)
*               endif
*            enddo
*         enddo
      GPPU8=((XDHA+XSHA+XCHA+XBHA+XTHA)*XDUP+
     &(XUHA+XSHA+XCHA+XBHA+XTHA)*XDDP+
     &(XUHA+XDHA+XCHA+XBHA+XTHA)*XDSP+
     &(XUHA+XDHA+XSHA+XBHA+XTHA)*XDCP+
     &(XUHA+XDHA+XSHA+XCHA+XTHA)*XDBP+
     &(XUHA+XDHA+XSHA+XCHA+XBHA)*XDTP+
c should be the antiquark in XD: P.Nason, 28-11-2002
     &(XDBHA+XSHA+XCHA+XBHA+XTHA)*XDUbP+
     &(XUBHA+XSHA+XCHA+XBHA+XTHA)*XDDbP+
     &(XUBHA+XDBHA+XCHA+XBHA+XTHA)*XDSP+
     &(XUBHA+XDBHA+XSHA+XBHA+XTHA)*XDCbP+
     &(XUBHA+XDBHA+XSHA+XCHA+XTHA)*XDBbP+
     &(XUBHA+XDBHA+XSHA+XCHA+XBHA)*XDTbP
     &)*XGPROB
      GPPV8=GPPU8/XA/XB/XC**3
*       call chk(tmp,gppu8)
C*********************************************************************
c q g -> qbarprime
*         tmp=0
*         do j=1,6
*            do k=1,6
*               if(k.ne.j) then
*                  tmp=tmp+fa(j)*fb(0)*fd(-k)+fa(-j)*fb(0)*fd(k)
*               endif
*            enddo
*         enddo
      GPPU9=((XDHA+XSHA+XCHA+XBHA+XTHA)*XDUBP+
     &(XUHA+XSHA+XCHA+XBHA+XTHA)*XDDBP+
     &(XUHA+XDHA+XCHA+XBHA+XTHA)*XDSP+
     &(XUHA+XDHA+XSHA+XBHA+XTHA)*xdcbp+
     &(XUHA+XDHA+XSHA+XCHA+XTHA)*XDBBP+
     &(XUHA+XDHA+XSHA+XCHA+XBHA)*XDTBP+
c should be the quark in XD: P.Nason, 28-11-2002
     &(XDBHA+XSHA+XCHA+XBHA+XTHA)*XDUP+
     &(XUBHA+XSHA+XCHA+XBHA+XTHA)*XDDP+
     &(XUBHA+XDBHA+XCHA+XBHA+XTHA)*XDSP+
     &(XUBHA+XDBHA+XSHA+XBHA+XTHA)*xdcp+
     &(XUBHA+XDBHA+XSHA+XCHA+XTHA)*XDBP+
     &(XUBHA+XDBHA+XSHA+XCHA+XBHA)*XDTP
     &)*XGPROB
      GPPV9=GPPU9/XA/XB/XC**3
*       call chk(tmp,gppu9)
C*********************************************************************
c q g -> qbar
*         tmp=0
*         do j=1,6
*            tmp=tmp+fa(j)*fb(0)*fd(-j)+fa(-j)*fb(0)*fd(j)
*         enddo
      GPPU10=(XUHA*XDUBP+XDHA*XDDBP+XSHA*
     &XDSP+XCHA*xdcbp+XBHA*XDBBP+XTHA*XDTBP+
     &XUBHA*XDUP+XDBHA*XDDP+XSHA*XDSP+
     &XCHA*XDCP+XBHA*XDBP+XTHA*XDTP
     &)*XGPROB
      GPPV10=GPPU10/XA/XB/XC**3
*       call chk(tmp,gppu10)
C*********************************************************************
c q qbar -> q
*         tmp=0
*         do j=1,6
*            tmp=tmp+fa(j)*fb(-j)*fd(j)+fa(-j)*fb(j)*fd(-j)
*         enddo
      GPPU11=XUHA*XUBHB*XDUP+XDHA*XDBHB
     &*XDDP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP
     &+XBHA*XBHB*XDBP+XTHA*XTHB*XDTP+
     &XUBHA*XUHB*XDUBP+XDBHA*XDHB
     &*XDDBP+XSHA*XSHB*XDSP+XCHA*XCHB*xdcbp
     &+XBHA*XBHB*XDBBP+XTHA*XTHB*XDTBP
      GPPV11=GPPU11/XA/XB/XC**3
*       call chk(tmp,gppu11)
C*********************************************************************
c q qbar -> g
*         tmp=0
*         do j=1,6
*            tmp=tmp+fa(j)*fb(-j)*fd(0)
*         enddo
      GPPU12=(XUHA*XUBHB+XDHA*XDBHB+XSHA*
     &XSHB+XCHA*XCHB+XBHA*XBHB+XTHA*XTHB)*XDGP
      GPPV12=GPPU12/XA/XB/XC**3
*       call chk(tmp,gppu12)
C*********************************************************************
c q g -> q
*         tmp=0
*         do j=1,6
*            tmp=tmp+fa(j)*fb(0)*fd(j)+fa(-j)*fb(0)*fd(-j)
*         enddo
      GPPU13=(XUHA*XDUP+XUBHA*XDUBP+XDHA*
     &XDDP+XDBHA*XDDBP+2.*XSHA*XDSP+
     &XCHA*(XDCP+xdcbp)+XBHA*(XDBP+XDBBP)+XTHA*(
     &XDTP+XDTBP))*XGPROB
      GPPV13=GPPU13/XA/XB/XC**3
*       call chk(tmp,gppu13)
C*********************************************************************
c q g -> g
*         tmp=0
*         do j=1,6
*            tmp=tmp+fa(j)*fb(0)*fd(0)+fa(-j)*fb(0)*fd(0)
*         enddo
      GPPU14=(XUHA+XUBHA+XDHA+XDBHA+2.*XSHA
     &+2.*XCHA+2.*XBHA+2.*XTHA)*XGPROB*XDGP
      GPPV14=GPPU14/XA/XB/XC**3
*       call chk(tmp,gppu14)
C*********************************************************************
c g g -> g
      GPPU15=XGPROA*XGPROB*XDGP/2.
      GPPV15=GPPU15/XA/XB/XC**3
C*********************************************************************
c g g -> q
*         tmp=0
*         do j=1,6
*            tmp=tmp+fa(0)*fb(0)*fd(j)+fa(0)*fb(0)*fd(-j)
*         enddo      
*         tmp=tmp/2
      GPPU16=XGPROA*XGPROB*(XDUP+XDUBP+XDDP+XDDBP
     &+XDSP*2.+XDCP+xdcbp+XDBP+XDBBP+XDTP+XDTBP)/2.
      GPPV16=GPPU16/XA/XB/XC**3
*       call chk(tmp,gppu16)
c      ENDIF
      RETURN
      END



*      subroutine chk(tmp,gppu16)
*      real * 8 tmp,gppu16
*      if(.not.(tmp.eq.0.and.gppu16.eq.0).and.
*     #    abs(tmp-gppu16)/abs(tmp+gppu16).gt.1d-10) then
*         write(*,*) tmp,gppu16
*         pause
*      endif
*      end

