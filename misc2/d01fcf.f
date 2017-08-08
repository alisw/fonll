c...this fake d01fcf subroutine actually performs the integration by
c...calling VEGAS. NB. ncalls and itmx have been tuned to this specific
c...problem.
      SUBROUTINE D01FCF(NDIMi,A,B,MINPTS,MAXPTS,F,EPS,ACCi,LENWRK,
     *                  WRKSTR,FINVAL,IFAIL)
      implicit none
      integer minpts,maxpts,lenwrk,ifail,ndim,ndimi
      real*8 A(NDIMi), B(NDIMi), WRKSTR(LENWRK)
      real*8 EPS, FINVAL,acci
      integer k
      real * 8 tmp
      real * 8 f,accum,accerr,epsint,coeff,res,err
      real * 8 xl,xu,acc,xi,si,si2,swgt,schi,avgi,sd,chi2a,alph
      integer ncall,itmx,nprn,ndo,it,ndmx,mds
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/bveg3/alph,ndmx,mds
      real * 8 calls,ti,tsi
      common/bveg4/calls,ti,tsi
      real * 8 avgip
      integer j
      external f
      
      ndim=ndimi
      acc=eps
      coeff=1d0
      
      do j=1,2
         xl(j) = a(j)
         xu(j) = b(j)
      enddo

      ncall=1000
      itmx=1
      call vegas(f,avgi,sd,chi2a)
      finval = avgi
      acci   = sd/finval
c      print*,finval,acci*finval
      itmx=10
      ncall=4000
      call vegas1(f,avgi,sd,chi2a)
      finval = avgi
      acci   = sd/finval     
c      print*,finval,acci*finval

      end
