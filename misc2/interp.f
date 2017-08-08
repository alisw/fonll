c...program to read a table in pt and y, include a smearing function, 
c...and convolute the pt distribution with a non-perturbative
c...fragmentation function
c...It can also perform integrations over one or both variables of
c...the table
c...MC 1998-2002

      implicit none
      real*8 pt(1000),y(1000),
     #       sum(1000),ptvec(100),sumvec(100),
     #       sumint(1000)
      real*8 x,sqrts,convol,dgauss,norm,sumintvec(100)
      real*8 splinesumint(10000),yvec(100)
      real*8 splinesum(10000),spline(10000)
      real*8 par1,par2,qm,enh1,enh2,br
      integer iwhatnp
      real*8 smear,tiny
      integer i,ntot,npt,ny,k
      common/int/sumvec,ptvec,x,sqrts,npt
      common/int2/spline,yvec,ny
      common/nppar/par1,par2,iwhatnp
      common/kinematics/qm,enh1,enh2
      external convol,dnp
      real*8 ptmin,ptmax,ymin,ymax
      integer ptpoints,ypoints,icross
      common/tables/ptmin,ptmax,ymin,ymax,ptpoints,ypoints,icross
      data tiny /1d-8/
      logical logfit
      common/whatfit/logfit
c      data logfit /.true./ 
      data logfit /.false./ 
      integer iwhatdpt
      common/whatdpt/iwhatdpt
      integer iunit,itm
      data iunit /37/,itm /53/
      real*8 mult
      common/symmrap/mult
      character*30 tablefile
      integer iformat
      real*8 dummy,mv,ml,rs
      integer jprefix
      character*30 prefix
      real*8 csmear
      data csmear /5d0/
            
      write(*,*) ' enter prefix for this run'
      read(*,*) prefix
      do jprefix=len(prefix),1,-1
         if(prefix(jprefix:jprefix).ne.' ') goto 9
      enddo
 9    continue
      open(unit=itm,file=prefix(1:jprefix)//'-interp.log')
      write(itm,*) prefix(1:jprefix)

c.....0->only interpolate, 1->\int dydpt, 2->\int dy, 3->\int dpt
      write(*,*) 'Run type:'
      write(*,*) '0 -> only interpolate'
      write(*,*) '1 -> int dydpt'
      write(*,*) '2 -> int dy'
      write(*,*) '3 -> int dpt'
      read(*,*) icross
      write(itm,*) icross,'  ',
     # '! icross:0->only interpol.,1->int dydpt,2->int dy,3->int dpt'

c.....output as: 1 -> d/dpt, 2 -> d/dpt^2
      write(*,*) 'output type: 1 -> d/dpt, 2 -> d/dpt^2'
      read(*,*) iwhatdpt
      write(itm,*) iwhatdpt,'   ! output type: 1 -> d/dpt, 2 -> d/dpt^2'
      
      write(*,*) 'ptmin,ptmax,ptpoints'
      read(*,*) ptmin,ptmax,ptpoints
      write(itm,*) ptmin,ptmax,ptpoints,'    ! ptmin,ptmax,ptpoints' 

      write(*,*) 'ymin,ymax,ypoints'
      read(*,*) ymin,ymax,ypoints
      write(itm,*) ymin,ymax,ypoints,'    ! ymin,ymax,ypoints' 
      
c.....iwhatpt: 0 -> no non pert. corr, 1 -> kartvelishvili, 2 -> peterson
c.....br -> branching ratio
c.....par1, par2: params. for the non-perturbative FF
c.....br -> branching ratio
c.....par1, par2: params. for the non-perturbative FF
      write(*,*) 'Non-perturbative FF: iwhatnp, br, par1, par2'
      write(*,*) 'iwhatpt = 0 -> no non pert. corr'
      write(*,*) 'iwhatpt = 1 -> Kartvelishvili' 
      write(*,*) 'iwhatpt = 2 -> Peterson'
      write(*,*) 'br : branching ration from quark to hadron'
      write(*,*) 'par1, par2: parameters for non-perturbative FF'
      read(*,*) iwhatnp,br,par1,par2
      write(itm,*) iwhatnp,br,par1,par2,
     #'  ! iwhatnp,br,par1,par2'
     
c.....quark mass, hadron 1 energy, hadron 2 energy (lab frame)      
      write(*,*) 'quark mass, hadr. 1 energy, hadr. 2 energy (labframe)'     
      read(*,*) qm,enh1,enh2
      
c....NB: when integrating over rapidity above y=0 only, 
c....negative rapidities can be accounted for by
c....multiplying by two the result!
      write(*,*) 'When integrating over rapidity above y=0 only,'
      write(*,*) 'negative rapidities can be accounted for by'
      write(*,*) 'multiplying by mult=2 the result!'
      write(*,*) 'Enter mult (usually 1 or 2):'
      read(*,*) mult
      write(itm,*) mult,'    ! mult'
      
c....table filename 
      write(*,*) 'Enter table format and table filename'
      write(*,*) 'Format: 0 -> old format (pt,et,y,mv,ml,rs,sum)'
      write(*,*) '        1 -> new format (pt,y, smeared sum)'
      read(*,*) iformat, tablefile
      do k=len(tablefile),1,-1
         if(tablefile(k:k).ne.' ') goto 11
      enddo
 11   continue
      write(itm,*) 
     #       iformat,' ',tablefile(1:k),' ! table format and filename'
      
      close(itm)

c....end of inputs

c...read in the data. The table should be pt, y, dsigma/dpt^2/dy
c...y should be the fastest rotating variable     
      open(unit=iunit,file=tablefile(1:k))      
      do i = 1,1000
c...read old format tables
        if(iformat.eq.0) then
           read(iunit,*,end=50) pt(i),dummy,y(i),mv,ml,rs,sum(i)
	   sum(i) = mv + (rs-ml)*pt(i)**2/(pt(i)**2 + csmear**2*qm**2)
c...read new format tables
	elseif(iformat.eq.1) then
           read(iunit,*,end=50) pt(i),y(i),sum(i)
	else
	   write(*,*) 'Table format not implemented'
	   stop
	endif
        if(sum(i).le.0.d0) then
	   sum(i) = tiny
	endif
      enddo
 50   close(iunit)

c...determine total number of points in the table    
      ntot = i-1
         
c...determine number of points in y and pt	 
      i = 2
 60   if(pt(i).ne.pt(1)) then
         ny = i-1
      else
         i = i+1
	 goto 60
      endif
 
      npt = ntot/ny

c....filling up pt vectors for FO and FONLL
      do i=1,ny
         yvec(i) = y(i)
         do k=1,npt
	   ptvec(k) = pt((k-1)*ny+1+i-1)
c...logs or weights used for easier interpolation
           if(logfit) then
	      sumvec(k) = log(sum((k-1)*ny+1+i-1))
           else
              sumvec(k) = sum((k-1)*ny+1+i-1)*smear(ptvec(k))
           endif
           splinesum(k + (i-1)*npt) = sumvec(k)
	 enddo
c....calculating the convolution
         if(iwhatnp.gt.0) then	 
	 do k=1,npt
           sqrts = 2*ptvec(npt)
c           norm = dgauss(d,0.d0,1.d0,1d-3)
           norm = 1.
	   x = 2*ptvec(k)/sqrts
	   sumintvec(k) = dgauss(convol,x,1.d0,1d-3)/norm/2./ptvec(k)*br
	   if(sumintvec(k).le.0d0) then
	       sumintvec(k) = tiny
	   endif
	   sumint((k-1)*ny+1+i-1) = sumintvec(k)
c...logs or weights used for easier interpolation
           if(logfit) then
	     splinesumint(k + (i-1)*npt) = log(sumintvec(k))
           else
             splinesumint(k + (i-1)*npt) = sumintvec(k)*smear(ptvec(k))
           endif
         enddo
	 endif
      enddo
      
c....rewrite out the numbers we have read, plus the convolution
      open(15,file=prefix(1:jprefix)//'-inputs.dat')
      write(15,*) '(* dsigma/dpt^2/dy (pb/GeV^2)'
      write(15,*),'(*   pt        y         FONLL    FONLL + conv'          
      do i=1,ntot
         write(15,'(1x,2(f9.4,1x),1x,2(e11.5,2x))') 
     #                       pt(i),y(i),sum(i),sumint(i)
      enddo
      close(15)
      
c.....write out an interpolated point or table
c.....FONLL, pure quark
      open(98,file=prefix(1:jprefix)//'-fonll.dat')
      call writetable(iwhatdpt,98,splinesum)
      close(98)
      if(iwhatnp.gt.0) then
c.....FONLL convoluted with non-pert. fragm. function (see convol)
         open(99,file=prefix(1:jprefix)//'-fonllconv.dat')
         call writetable(iwhatdpt,99,splinesumint)
         close(99)
      endif
         
      if(icross.eq.0) stop
            
c...make integrals and write topdrawer file

      open(45,file=prefix(1:jprefix)//'.top')
      write(45,*) 'set scale y log'
      write(45,*) 'set order x y dummy'

c.....fill the spline vector
c.....the FIXED ORDER (mv)
c      do i=1,10000
c         spline(i) = splinemv(i)
c      enddo      
c.....do the integrals
c      write(45,*) '(* FIXED ORDER'
c      call integrali(ptvec(npt))
c      write(45,*) 'join 1 dots'


c.....fill the spline vector
c.....the NON-convoluted one (i.e. the FONLL matched calculation)
      do i=1,10000
         spline(i) = splinesum(i)
      enddo      
c.....do the integrals
      write(45,*) '(* FONLL, NOT convoluted'
      call integrali(ptvec(npt))
      write(45,*) 'join 1 dashes'


c.....fill the spline vector
      if(iwhatnp.gt.0) then
c.....the convoluted one
      do i=1,10000
         spline(i) = splinesumint(i)
      enddo      
c.....do the integrals
      write(45,*) '(* FONLL + convolution'
      write(45,*) '(* iwhatnp = ',iwhatnp,', par1 = ',par1,', br = ',br
      call integrali(ptvec(npt))
      write(45,*) 'join 1'
      endif
      
      close(45)
      
      end

      subroutine writetable(iwhat,outunit,splineinput)
      implicit none
      real*8 splineinput(10000)
      real*8 spline(10000),yvec(100)
      integer ny,outunit,npt,nrap,i,k,iwhat
      real*8 dpt,dy,dpty,xx(2),res
      common/int2/spline,yvec,ny
      real*8 ptmin,ptmax,ymin,ymax
      integer ptpoints,ypoints,icross
      common/tables/ptmin,ptmax,ymin,ymax,ptpoints,ypoints,icross
      
c.....fill the spline vector
      do i=1,10000
         spline(i) = splineinput(i)
      enddo      
      
      npt = ptpoints
      nrap = ypoints
            
      if(iwhat.eq.2) then      
       write(outunit,*)'    p_T       y      dsigma/dpt^2/dy (pb/GeV^2)'
      elseif(iwhat.eq.1) then
       write(outunit,*)'    p_T       y      dsigma/dpt/dy (pb/GeV)'
      else
          write(*,*) 'iwhat=?? in writetable'
          stop
      endif
      
      write(outunit,*) ' '

      if(npt.eq.1) then
         dpt = 0.d0
      else
         dpt = (ptmax-ptmin)/float(npt-1)
      endif
      do i=1,npt
          xx(1) = ptmin + dpt*float(i-1)  
          if(nrap.eq.1) then
              dy = 0.d0
          else
              dy = (ymax-ymin)/float(nrap-1)
          endif   
          do k=1,nrap
              xx(2) = ymin + dy*float(k-1)
	      res = dpty(2,xx)
	      if(iwhat.eq.2) then         ! convert ds/dpt -> ds/dpt2
	          res = dpty(2,xx)/2./xx(1)
	      endif
              write(outunit,'(1x,2(f11.6,1x),1x,e13.6)')
     #                       xx(1), xx(2), res
          enddo
c          write (outunit,*) 'join'
      enddo
      
      end

      subroutine integrali(ptlim)
      implicit real*8 (a-h,o-z)
      common/kinematics/qm,enh1,enh2
      real*8 ptmin,ptmax,ymin,ymax
      integer ptpoints,ypoints,icross
      common/tables/ptmin,ptmax,ymin,ymax,ptpoints,ypoints,icross
      integer iwhatdpt
      common/whatdpt/iwhatdpt
      real*8 mult
      common/symmrap/mult


      if(ptmax.gt.ptlim) then
         ptmax = ptlim
      endif

      npt = ptpoints
      nrap = ypoints
            
      write(45,*)'(* icross = ',icross

      if(icross.eq.1) then       
c....integration over rapidity and pt, for total cross sections
      write(45,*)'(* sigma [pb], pt>ptmin, 2 x ymin<y<ymax'
      write(45,*)'(* ptmin,ptmax = ',ptmin,ptmax
      write(45,*)'(* ymin,ymax = ',ymin,ymax, '[but x mult]'
            
      if(npt.gt.1) then
          dpt = (ptmax-ptmin)/float(npt-1)
      else
          dpt=0.
      endif

      do i=1,npt
          pt = ptmin + dpt*float(i-1)

        call integrate(pt,ptlim,ymin,ymax,result,error)
c...to account for negative rapidities (mult must be set to two)
	  result = result*mult
	  error = error*mult

c          write(*,'(1x,f12.4,2(2x,e15.4))') pt, result
          write(45,'(1x,f12.4,2(2x,e15.4))') pt, result
      enddo

      elseif(icross.eq.2) then
c....integration over rapidity, result is ds/dpt^2, as a function of pt
      if(iwhatdpt.eq.2) then
         write(45,*)'(* dsigma/dpt^2 [pb/GeV^2]'
      elseif(iwhatdpt.eq.1) then
         write(45,*)'(* dsigma/dpt [pb/GeV]'
      endif
      write(45,*)'(* ymin,ymax = ',ymin,ymax, '[but x mult]'

      if(npt.gt.1) then
          dpt = (ptmax-ptmin)/float(npt-1)
      else
          dpt=0.
      endif
              
      do i=1,npt
          pt = ptmin + dpt*float(i-1)

          call integratey(pt,ymin,ymax,result,error)
c...to account for negative rapidities (mult must be set to two)
	  result = result*mult
	  error = error*mult
          
	  if(iwhatdpt.eq.2) then
c           write(*,'(1x,f12.4,2(2x,e15.4))') pt, result
           write(45,'(1x,f12.4,2(2x,e15.4))') pt, result  
	  elseif(iwhatdpt.eq.1) then
c...dsigma/dpt
c           write(*,'(1x,f12.4,2(2x,e15.6))') pt, result*2*pt,error*2*pt
           write(45,'(1x,f12.4,2(2x,e15.6))') pt, result*2*pt,error*2*pt
          endif
      enddo


      elseif(icross.eq.3) then
c.....integration over pt, result as a function of rapidity
      write(45,*)'(* dsigma/dy [pb]'
      write(45,*)'(* ptmin,ptmax = ',ptmin,ptmax

      if(nrap.gt.1) then
          dy = (ymax-ymin)/float(nrap-1)
      else
          dy=0.
      endif
      
      do i=1,nrap
          rap = ymin + dy*float(i-1)

          call integratept(rap,ptmin,ptmax,res,err)
          
c          write(*,'(1x,f12.4,2(2x,e15.4))') rap, res
          write(45,'(1x,f12.4,2(2x,e15.4))') rap, res
          
      enddo
      
      endif

      end


      
      real*8 function convol(z)
      implicit none
      real*8 z,sumvec(100),ptvec(100),pt,ddivdif,dsigdpt,x,sqrts
      real*8 dnp,smear
      integer npt
      common/int/sumvec,ptvec,x,sqrts,npt
      logical logfit
      common/whatfit/logfit
            
      pt = sqrts*x/2.

      if (logfit) then
c.....interpolated value converted back from log      
         dsigdpt = 2*pt/z*exp(ddivdif(sumvec,ptvec,npt,pt/z,4))
      else
c...weights used
         dsigdpt = 2*pt/z*ddivdif(sumvec,ptvec,npt,pt/z,4)/smear(pt/z)
      endif
      
      convol = 1./z*dnp(z)*dsigdpt
                  
      end
      
      real*8 function dnp(z)
      implicit none
      real*8 z,par1,par2
      real*8 kart,peterson
      integer iwhatnp
      common/nppar/par1,par2,iwhatnp

      if(iwhatnp.eq.1) then
         dnp = kart(z)
      elseif(iwhatnp.eq.2) then
         dnp = peterson(z)
      else
         write(*,*) 'iwhatnp = ',iwhatnp,': undefined np function'
	 stop
      endif
      
      end

      real*8 function kart(z)
      implicit none
      real*8 z,alpha,anorm,par1,par2
      integer iwhatnp
      common/nppar/par1,par2,iwhatnp

      alpha = par1      
      anorm = (1d0+alpha)*(2d0+alpha)
      kart = (1d0 - z)* z**alpha   
      kart = kart*anorm
      
      end
      
      real*8 function peterson(z)
      implicit none
      real*8 z,pet,eps,anorm,root,par1,par2
      integer ini,iwhatnp
      common/nppar/par1,par2,iwhatnp
      data ini/0/
      save ini,anorm
      
      eps = par1

      if(ini.eq.0) then
        root = sqrt(4*eps-eps**2)
        anorm = (eps**2 - 6*eps + 4.)/(4.-eps)/root*
     #        (atan(eps/root) + atan((2.-eps)/root)) + 0.5*log(eps) +
     #        1./(4.-eps)      
        ini = 1
      endif
      pet = z*(1.-z)**2/((1.-z)**2 + eps*z)**2
      peterson = pet/anorm

      end

 
 
 
      function raplim(rap,pt)
      implicit real*8 (a-h,o-z)
      common/kinematics/qm,enh1,enh2

      shift = 0.
      if(enh1.eq.enh2) then
           rpcm = rap
      else
           shift = log(enh1/enh2)/2.
           rpcm = rap - shift
      endif

      gs = 4.*enh1*enh2
	
c      raplim = dasinh(sqrt((gs+qm**2)**2/4.
c     #         /gs/(pt**2+qm**2) - 1d0))
c..forse era sbagliato?? MC, 9/12/2002
      raplim = dasinh(sqrt(gs/4./(pt**2+qm**2) - 1d0))
      if(abs(rpcm).gt.raplim) then
c          write(*,*) ' Out of kinematical limits'
c          write(*,*) ' pt = ',pt,'    rap = ',rap
c          write(*,*) ' raplim = ',-raplim+shift, raplim+shift
          raplim = 0.d0
      else
          raplim = 1.d0
      endif

      end

c.....integrates over both pt and y
      subroutine integrate(ptmin,ptmax,ymin,ymax,result,error)
      PARAMETER        (indim=2,imaxpts=2000*indim,
     +               iLENWRK=(indim+2)*(1+imaxpts/(2**indim
     +               +2*indim*indim+2*indim+1)))
      parameter (iwk = (2*indim+3)*(1+imaxpts/(2**indim +
     #                2*indim*(indim+1)+1))/2)
      implicit real*8 (a-h,o-z)
*     .. Local Arrays ..
      dimension A2(indim), B2(indim), WRKSTR(iLENWRK)
      real*8 wk(iwk)
*     .. External Subroutines ..
      EXTERNAL         dpty,dptyvegas
      
      A2(1) = ptmin
      B2(1) = ptmax
      a2(2) = ymin
      b2(2) = ymax
   
      EPS = 5d-3
      iMINPTS = 0
      IFAIL = -1

      result = 0d0
c..call to the NAG library routing
c      CALL D01FCF(iNDIM,A2,B2,iMINPTS,iMAXPTS,dpty,EPS,relerr,
c     #            iLENWRK,WRKSTR,result,IFAIL)
c      error = result*relerr
c..call to VEGAS
c....variable change: y = 1/pT^6
      a2(1) = (ptmin**(-6+1))/(-6+1)
      B2(1) = (ptmax**(-6+1))/(-6+1)
      CALL D01FCF(iNDIM,A2,B2,iMINPTS,iMAXPTS,dptyvegas,EPS,relerr,
     #            iLENWRK,WRKSTR,result,IFAIL)
      error = result*relerr
c..call to DADMUL, CERNLIB routine
c      call dadmul(dpty,indim,a2,b2,iminpts,imaxpts,eps,
c     #            wk,iwk,result,relerr,nfnevl,ifail)
c      error = result*relerr
c      print*,result,error,nfnevl,ifail
       
      end

c.....integrates over y and returns dsigma/dpt^2
      subroutine integratey(ptin,ymin,ymax,result,error)
      implicit real*8 (a-h,o-z)
*     .. External Subroutines ..
      EXTERNAL         dpty1
      common/ptvalue/pt
      
      pt=ptin
      
      A = ymin
      B = ymax
   
      EPSrel= 1d-2
      IFAIL = 0

      result = 0d0
c      result = d01ahf(a,b,epsrel,npts,relerr,dpty1,-1,ifail)
      result = dgauss(dpty1,a,b,epsrel)
      error = result*relerr

      end
	
c.....integrates over pt and returns dsigma/dy	
      subroutine integratept(yin,ptmin,ptmax,result,error)
      implicit real*8 (a-h,o-z)
*     .. External Subroutines ..
      EXTERNAL         dpty2
      common/yvalue/y
      
      y=yin

      A = ptmin 
      B = ptmax
   
      EPSrel= 1d-2
      IFAIL = 0

      result = 0d0
c      result = d01ahf(a,b,epsrel,npts,relerr,dpty2,-1,ifail)
      result = dgauss(dpty2,a,b,epsrel)
      error = result*relerr

      end
	
c.....returns dsigma/dpt/dy at a given pt
      function dpty2(pt)
      implicit real*8 (a-h,o-z)
      real*8 xx(2)
      common/yvalue/yval

      xx(1) = pt
      xx(2) = yval

      dpty2 = dpty(2,xx)

      end

c.....returns dsigma/dpt^2/dy at a given yval
      function dpty1(yval)
      implicit real*8 (a-h,o-z)
      real*8 xx(2)
      common/ptvalue/pt

      xx(1) = pt
      xx(2) = yval
            
      dpty1 = dpty(2,xx)/2./xx(1)

      end


c.....returns the interpolated dsigma/dpt/dy differential distribution
      function dpty(indim,xx)
      implicit real*8 (a-h,o-z)
c      dimension xx(indim)
      dimension xx(*)
      real*8 spline(10000),ptvec(100),yvec(100),sumvec(100)
      real*8 wk(3,100,100),xi,yi,zi
      common/int/sumvec,ptvec,x,sqrts,npt
      common/int2/spline,yvec,ny
      logical logfit
      common/whatfit/logfit

      pt = xx(1)
      yval = xx(2)
      
      if (raplim(yval,pt).eq.0) then
         dpty = 0d0
	 return
      endif	 

      xi = pt
      yi = yval

*     Evaluate the spline.
      call RGBI3P(1,npt,ny,ptvec,yvec,spline,1,XI,YI, ZI,IER, WK)

      if(logfit) then
c.....from log scale back to real values
         zi = exp(zi)
      else
c.....weights used
         zi = zi/smear(xi)
      endif
      
c      print*,xi,yi,zi

c.....from dpt2 to dpt
      dpty = zi*2*xx(1)
            
      end
 
c....dpty in VEGAS version
      real*8 function dptyvegas(xx,wgt)
      implicit none
      real*8 xx(2),wgt,dpty

c...variable change and jacobian      
      xx(1) = ((-6+1)*xx(1))**(1./(-6+1))
      dptyvegas = dpty(2,xx) * xx(1)**6
      
      end
 

      real*8 function smear(x)
      implicit none
      real*8 x
      
      smear = x**6 + 10.
c      smear = (x+5.)**6
      
      end


*
* $Id: ddivdif.F,v 1.1.1.1 1996/02/15 17:48:36 mclareni Exp $
*
* $Log: ddivdif.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:36  mclareni
* Kernlib
*
*
c#include "kernnum/pilot.h"
      FUNCTION ddivdif(F,A,NN,X,MM)
      implicit real*8 (a-h,o-z)
      real*8 A(NN),F(NN),T(20),D(20)
      LOGICAL EXTRA
      LOGICAL MFLAG,RFLAG
      DATA MMAX/10/
C
C  TABULAR INTERPOLATION USING SYMMETRICALLY PLACED ARGUMENT POINTS.
C
C  START.  FIND SUBSCRIPT IX OF X IN ARRAY A.
      IF( (NN.LT.2) .OR. (MM.LT.1) ) GOTO 20
      N=NN
      M=MIN0(MM,MMAX,N-1)
      MPLUS=M+1
      IX=0
      IY=N+1
      IF(A(1).GT.A(N)) GOTO 4
C     (SEARCH INCREASING ARGUMENTS.)
    1    MID=(IX+IY)/2
         IF(X.GE.A(MID)) GOTO 2
            IY=MID
            GOTO 3
C        (IF TRUE.)
    2       IX=MID
    3    IF(IY-IX.GT.1) GOTO 1
         GOTO 7
C     (SEARCH DECREASING ARGUMENTS.)
    4    MID=(IX+IY)/2
         IF(X.LE.A(MID)) GOTO 5
            IY=MID
            GOTO 6
C        (IF TRUE.)
    5       IX=MID
    6    IF(IY-IX.GT.1) GOTO 4
C
C  COPY REORDERED INTERPOLATION POINTS INTO (T(I),D(I)), SETTING
C  *EXTRA* TO TRUE IF M+2 POINTS TO BE USED.
    7 NPTS=M+2-MOD(M,2)
      IP=0
      L=0
      GOTO 9
    8    L=-L
         IF(L.GE.0) L=L+1
    9    ISUB=IX+L
         IF((1.LE.ISUB).AND.(ISUB.LE.N)) GOTO 10
C        (SKIP POINT.)
            NPTS=MPLUS
            GOTO 11
C        (INSERT POINT.)
   10       IP=IP+1
            T(IP)=A(ISUB)
            D(IP)=F(ISUB)
   11    IF(IP.LT.NPTS) GOTO 8
      EXTRA=NPTS.NE.MPLUS
C
C  REPLACE D BY THE LEADING DIAGONAL OF A DIVIDED-DIFFERENCE TABLE, SUP-
C  PLEMENTED BY AN EXTRA LINE IF *EXTRA* IS TRUE.
      DO 14 L=1,M
         IF(.NOT.EXTRA) GOTO 12
            ISUB=MPLUS-L
            D(M+2)=(D(M+2)-D(M))/(T(M+2)-T(ISUB))
   12    I=MPLUS
         DO 13 J=L,M
            ISUB=I-L
            D(I)=(D(I)-D(I-1))/(T(I)-T(ISUB))
            I=I-1
   13    CONTINUE
   14 CONTINUE
C
C  EVALUATE THE NEWTON INTERPOLATION FORMULA AT X, AVERAGING TWO VALUES
C  OF LAST DIFFERENCE IF *EXTRA* IS TRUE.
      SUM=D(MPLUS)
      IF(EXTRA) SUM=0.5*(SUM+D(M+2))
      J=M
      DO 15 L=1,M
         SUM=D(J)+(X-T(J))*SUM
         J=J-1
   15 CONTINUE
      ddivdif=SUM
      RETURN
C
   20 stop
      CALL KERMTR('E105.1',LGFILE,MFLAG,RFLAG)
      IF(MFLAG) THEN
         IF(LGFILE.EQ.0) THEN
            IF(MM.LT.1) WRITE(*,101) MM
            IF(NN.LT.2) WRITE(*,102) NN
         ELSE
            IF(MM.LT.1) WRITE(LGFILE,101) MM
            IF(NN.LT.2) WRITE(LGFILE,102) NN
         ENDIF
      ENDIF
      IF(.NOT.RFLAG) CALL ABEND
      RETURN
  101 FORMAT( 7X, 'FUNCTION ddivdif ... M =',I6,' IS LESS THAN 1')
  102 FORMAT( 7X, 'FUNCTION ddivdif ... N =',I6,' IS LESS THAN 2')
      END
