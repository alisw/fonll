      subroutine e01baf(m,x,y,lambda,c,lck,wrk,lwrk,ifail)
      implicit none
      integer m,lck,lwrk,ifail,k
      real * 8 x(m),y(m),
     #       wrk(lwrk)
      real * 8 lambda(lck),c(lck)
      do k=1,m
         lambda(k)=x(k)
         c(k)=y(k)
      enddo
      ifail=0
      end

      subroutine e02bbf(n,lambda,c,x,s,ifail)
      implicit none
      integer n,ifail
      real * 8 lambda(n),c(n)
      real * 8 x,s
      real * 8 xtmp,dfint
      integer na(1)
      na(1)=n-4
      xtmp=x
      s=dfint(1,xtmp,na,lambda,c)
      ifail = 0
      end

      subroutine dovegas(ndimi,f,maxpts,accum,accerr,epsint,
     #     coeff,res,err)
      implicit none
      integer maxpts,ndimi,k
      real * 8 tmp
      real * 8 f,accum,accerr,epsint,coeff,res,err
      real * 8 xl,xu,acc,xi,si,si2,swgt,schi
     # ,avgi,sd,chi2a,alph
      integer ndim,ncall,itmx,nprn,ndo,it,ndmx,mds
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/bveg3/alph,ndmx,mds
      real * 8 calls,ti,tsi
      common/bveg4/calls,ti,tsi
      real * 8 avgip
      integer j
      external f
c alph controls the rebinning in vegas
      alph=1.5
      ndim=ndimi
      do k=1,ndim
         xl(k)=0
         xu(k)=1
      enddo
      acc=-1
c      ncall = maxpts/5
      if(ndim.ge.3) then
         ncall = 4000
      else
         ncall = 1000
      endif


c... MC, 27/4/2005. Let VEGAS _really_ do the job 
c changed on 6/7/2007: now use epsint as given in input
      acc = epsint
c
      nprn=1
      alph=1.5
      if(ndim.ge.3) then
c	 acc = 2d-2
         itmx = 1
         ncall = 4000
         call vegas(f,avgi,sd,chi2a)
         itmx = 10
         ncall = 12000
         call vegas1(f,avgi,sd,chi2a)
      else
c	 acc = 1d-3
         itmx = 20
         ncall = 4000
         call vegas(f,avgi,sd,chi2a)
      endif
      res = avgi*coeff
      err = sd*coeff
      return
c........      



c... MC, 26/4/2005. Let VEGAS do the job (well, sort of...)
      nprn=1
      alph=1.5
      if(ndim.ge.3) then
	 acc = 2d-2
         ncall = 10000
         itmx = 1
         call vegas(f,avgi,sd,chi2a)
         itmx = 10
         ncall = 5000
         call vegas1(f,avgi,sd,chi2a) ! inizialize cumulants, not grid
	 if(dabs(sd/avgi).gt.acc) then
  	   alph=0	! freeze the grid
           itmx = 20    ! NB. previous iterations still count, so it's 10 here
           ncall = 10000
           call vegas2(f,avgi,sd,chi2a) ! don't initialize anything
         endif
      else
	 acc = 1d-3
         itmx = 1
         ncall = 4000
         call vegas(f,avgi,sd,chi2a) 
         itmx = 20
	 ncall = 4000
         call vegas1(f,avgi,sd,chi2a) ! inizialize cumulants, not grid
      endif
      res = avgi*coeff
      err = sd*coeff
      return
c........      
      
      itmx=1
      avgip=0
      do j=1,15
         if(j.eq.1) then
            call vegas(f,avgi,sd,chi2a)
         else
            call vegas2(f,avgi,sd,chi2a)
         endif
         avgip = avgip+ti
         write(*,*) ' vegas estimate:', avgi,' +- ',sd
         write(*,*) ' alternative estimate:',avgip/it
         if((accum+coeff*avgi).ne.0.d0)then
           tmp = abs(sqrt((coeff*sd)**2)/(accum+coeff*avgi))
         else
c this is arbitrary; freeze the grid (shouldn't make any difference)
           tmp=1.d8
         endif
         if(tmp.gt.epsint) then
c after the 5th iteration freeze the grid
            if(it.gt.5) alph=0
         else
            goto 998
         endif
      enddo
 998  continue
      avgip=avgip/it
      if(abs(avgi-avgip).gt.2*sd) then
         write(*,*) ' WARNING! alternative estimate differs'//
     # 'from vegas estimate! increase ncalls!'
         sd=abs(avgi-avgip)
         avgi=avgip
      endif
      res = avgi*coeff
      err = sd*coeff
      end


