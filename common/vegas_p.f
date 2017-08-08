      block data vegas0
      implicit  real*8 (a-h,o-z)
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/bveg3/alph,ndmx,mds
      data ndmx/50/,alph/1.5d0/,mds/1/
      data ncall/10000/,itmx/10/,nprn/ 0/,acc/-1.d0/,
     1    xl/1.d-3, 1.d-3,1.d-3,1.d-3,1.d-3,1.d-3,1.d-3,1.d0,1.d0,1.d0/,
     2    xu/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/
      data xi/500*1.d0/
      end
C
C
C      NCALL IS THE NUMBER OF CALLS TO VEGAS.
C      NPRN >  0 VEGAS PRINTS THE RESULTS OF EACH ITERATION.
C      NPRN 0 VEGAS PRINTS NOTHING.
C      NPRN < 0 VEGAS PRINTS ALL.
C      XL(I) IS LOWER INTEGRATION LIMIT ON I TH AXIS.
C       XU(I) IS UPPER INTEGRATION LIMIT ON I THE AXIS.
c
      subroutine vegas(fxn,avgi,sd,chi2a)
c
c   routine performs n dim monte carlo inte
c written by p lepage
c
      implicit none
      real * 8 fxn,avgi,sd,chi2a
      real * 8 xl,xu,acc
      integer ndim,ncall,itmx,nprn
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      real * 8 xi,si,si2,swgt,schi
      integer ndo,it
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      real * 8 alph
      integer ndmx,mds
      common/bveg3/alph,ndmx,mds
      real * 8 calls,ti,tsi
      common/bveg4/calls,ti,tsi
      integer num,num2
      common/seed/num,num2
      real * 8 d(50,10),di(50,10),xin(50),r(50),dx(10),dt(10),
     1     x(10)
      integer kg(10),ia(10)
      real * 8 rand(10)
      real * 8 one,dxg,dv2g,xnd,xjac,rc,
     #          f,f2,fb,f2b,wgt,ti2,xo,xn,dr
      integer j,nd,ng,npg,k,ndm,i
      data one/1.d0/
c ndmx = maximum subdivisions
c 
      num=1
c       num2 is irrelevant
      ndo=1
      do j=1,ndim
         xi(1,j)=one
      enddo
c
      entry vegas1(fxn,avgi,sd,chi2a)
c     initialises  cumulative  variables but not grid
      it=0
      si=0.
      si2=si
      swgt=si
      schi=si
      entry vegas2(fxn,avgi,sd,chi2a)
c     no initialisation
      nd=ndmx
      ng=1
      if(mds.eq.0)go to 2
      ng=(ncall/2.)**(1./ndim)
      mds=1
      if((2*ng-ndmx).lt.0)go to 2
      mds=-1
      npg=ng/ndmx+1
      nd=ng/npg
      ng=npg*nd
 2    k=ng**ndim
      npg=ncall/k
      if(npg.lt.2)npg=2
      calls=npg*k
      dxg=one/ng
      dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-one)
      xnd=nd
      ndm=nd-1
      dxg=dxg*xnd
      xjac=one/calls
      do j=1,ndim
         dx(j)=xu(j)-xl(j)
         xjac=xjac*dx(j)
      enddo
c
c     rebin preserving bin density
c
      if(nd.eq.ndo)go to 8
      rc=ndo/xnd
      do J=1,ndim
         k=0
         xn=0.
         dr=xn
         i=k
 4       k=k+1
         dr=dr+one
         xo=xn
         xn=xi(k,j)
 5       if(rc.gt.dr)go to 4
         i=i+1
         dr=dr-rc
         xin(i)=xn-(xn-xo)*dr
         if(i.lt.ndm)go to 5
         do i=1,ndm
            xi(i,j)=xin(i)
         enddo
         xi(nd,j)=one
      enddo
      ndo=nd
c
 8    if(nprn.ne.0)write(6,200)ndim,calls,it,itmx,acc
     1     ,mds,nd,(xl(j),xu(j),j=1,ndim)
c
      entry vegas3(fxn,avgi,sd,chi2a)
c     main integration loop
 9    it=it+1
      ti=0.
      tsi=ti
      do j=1,ndim
         kg(j)=1
         do i=1,nd
            d(i,j)=ti
            di(i,j)=ti
         enddo
      enddo
c
 11   fb=0.
      f2b=fb
      k=0
 12   k=k+1
      call randa(ndim,rand)
      wgt=xjac
      do j=1,ndim
         xn=(kg(j)-rand(j))*dxg+one
         ia(j)=xn
         if(ia(j).gt.1)go to 13
         xo=xi(ia(j),j)
         rc=(xn-ia(j))*xo
         go to 14
 13      xO=xi(ia(j),j)-xi(ia(j)-1,j)
         rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
 14      x(j)=xl(j)+rc*dx(j)
         wgt=wgt*xo*xnd
      enddo
c
      f=wgt
      f=f*fxn(x,wgt)
      f2=f*f
      fb=fb+f
      f2b=f2b+f2
      do j=1,ndim
         di(ia(j),j)=di(ia(j),j)+f
         if(mds.ge.0)d(ia(j),J)=d(ia(j),J)+f2
      enddo
      if(k.lt.npg) go to 12
c
 888  FORMAT(1X,'F',G14.6,'F2',G14.6,'FB',G14.6,'F2B',G14.6)
      f2b= sqrt(f2b*      NPG)
      f2b=(f2b-fb)*(f2b+fb)
 1661 FORMAT(1X,'F2B',G14.6,'NPG',  I10)
      ti=ti+fb
      tsi=tsi+f2b
 33   FORMAT(1X,'TSI',G14.6,'F2B',G14.6)
      if(mds.ge.0)go to 18
      do j=1,ndim
         d(ia(j),j)=d(ia(j),j)+f2b
      enddo
 18   k=ndim
 19   kg(k)=mod(kg(k),ng)+1
      if(kg(k).ne.1)go to 11
      k=k-1
      if(k.gt.0)go to 19
c
c final results for this iteration
c
      tsi=tsi*dv2g
      ti2=ti*ti
 88   format(1x,'tsi',g14.6)
      if(abs(tsi).gt.0) then
         wgt=ti2/tsi
         si=si+ti*wgt
         si2=si2+ti2
         swgt=swgt+wgt
         schi=schi+ti2*wgt
 995     FORMAT(1X,'SWGT',G14.6,'SI2',G14.6)
         avgi=si/swgt
         sd=swgt*it/si2
         chi2a=sd*(schi/swgt-avgi*avgi)/(it-.999d0)
         sd=dsqrt(one/sd)
      else
         write(*,*) '***VEGAS WARNING***: zero error!'
         write(*,*) ' we guess that integral is exact '
         avgi=ti
         sd=0
         chi2a=-1
         return
      endif
c
      if(nprn.eq.0)go to 21
      tsi=dsqrt(tsi)
      write(6,201)it,ti,tsi,avgi,sd,chi2a
      if(nprn.ge.0)go to 21
      do j=1,ndim
         write(6,202) j,(xi(i,j),di(i,j),d(i,j),i=1,nd)
      enddo
c
c      refine grid
c
 21   continue
      if(sd.eq.0) return
      do j=1,ndim
         xo=d(1,j)
         xn=d(2,j)
         d(1,j)=(xo+xn)/2.
         dt(j)=d(1,j)
         do i=2,ndm
            d(i,j)=xo+xn
            xo=xn
            xn=d(i+1,j)
            d(i,j)=(d(i,j)+xn)/3.
            dt(j)=dt(j)+d(i,j)
         enddo
         d(nd,j)=(xn+xo)/2.
         dt(j)=dt(j)+d(nd,j)
      enddo
c
      do j=1,ndim
         rc=0.
         do i=1,nd
            r(i)=0.
            if(d(i,j).le.0.)go to 24
            xo=dt(j)/d(i,j)
            r(i)=((xo-one)/xo/dlog(xo))**alph
            rc=rc+r(i)
 24         continue
         enddo
         rc=rc/xnd
         k=0
         xn=0.
         dr=xn
         i=k
 25      k=k+1
         dr=dr+r(k)
         xo=xn
         xn=xi(k,j)
 26      if(rc.gt.dr)go to 25
         i=i+1
         dr=dr-rc
         xin(i)=xn-(xn-xo)*dr/r(k)
         if(i.lt.ndm)go to 26
         do i=1,ndm
            xi(i,j)=xin(i)
         enddo
         xi(nd,j)=one
      enddo
c
      if(it.lt.itmx.and.acc*dabs(avgi).lt.sd)go to 9
 200  format(1X,'0input parameters for vegas:  ndim=',i3,
     1     '   ncall=',f8.0/28x,'  it=',i5,'    itmx=',i5/28x,
     2     '  acc=',g9.3/28x,'  mds=',i3,'     nd=',i4/28x,
     3     '  (xl,xu)=',(t40,'( ',g12.6,' , ',g12.6,' )'))
 201  format(///' integration by vegas' / '0iteration no.',i5,
     1     ':  integral=',g14.8/21x,'std dev =',g10.4 /
     2     ' accumulated results:   integral=',g14.8/
     3     24x,'std dev =',g10.4 / 24x,'chi**2 per it''n =',g10.4)
 202  format(1X,'0data for axis',i2,/,' ',6x,'x',7x,'  delt i ',
     1     2x,'conv','ce   ',11x,'x',7x,'  delt i ',2x,'conv','ce  '
     2     ,11x,'x',7x,'   delt i ',2x,'conv','CE  ',/,
     3     (1X,' ',3g12.4,5x,3g12.4,5x,3g12.4))
      return
      entry vegas4(fxn,avgi,sd,chi2a)
      avgi=si/swgt
      sd=swgt*it/si2
      chi2a=sd*(schi/swgt-avgi*avgi)/(it-.999d0)
      sd=dsqrt(one/sd)
      if(nprn.ne.0) write(6,201)it,0.d0,0.d0,avgi,sd,chi2a
      return
      end

c        subroutine save(ndim)
c        implicit real*8 (a-h,o-z)
c       implicit integer*4 (i-n)
c        common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
c
c      stores vegas data   (unit 7) for later initialisation
c
c        write(7,200) ndo,it,si,si2,swgt,schi,
c     1       ((xi(i,j),i=1,ndo),j=1,ndim)
c        return
c        entry restr(ndim)
c
c         enters initialisation data for vegas
c
c        read(7,200) ndo,it,si,si2,swgt,schi,
c     1    ((xi(i,j),i= 1,ndo),j=1,ndim)
c 200    format(2i8,4z16/(5z16))
c        return
c        end

      
      FUNCTION RANDOM(SEED)
*     -----------------
* Ref.: K. Park and K.W. Miller, Comm. of the ACM 31 (1988) p.1192
* Use seed = 1 as first value.
*
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION MINV,RANDOM
      SAVE
      PARAMETER(M=2147483647,A=16807,Q=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = SEED/Q
      LO = MOD(SEED,Q)
      SEED = A*LO - R*HI
      IF(SEED.LE.0) SEED = SEED + M
      RANDOM = SEED*MINV
      END

      subroutine randa(n,rand)
      implicit double precision (a-h,o-z)
      COMMON/SEED/NUM,NUM2
      common/caso/caso(5)
      dimension rand(10)
      do 1 i=1,n
      rand(i)=random(NUM)
1     continue
      do 2 i=1,5
      caso(i)=random(NUM)
2     continue
      return
      end
