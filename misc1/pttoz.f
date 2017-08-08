

      subroutine ptpoints(ebeam1,zmax1,ebeam2,zmax2,
     # ffact,y,xm,ptmax,iflagmode,z,pt)
      implicit none
      real * 8 ebeam1,zmax1,ebeam2,zmax2,ffact,y,xm,ptmax
      integer iflagmode
      real * 8 pt,z
c common
      integer jptseqfl
      common/cjptseqfl/jptseqfl
c local
      real * 8 q2min
      parameter (q2min=2d0)
      real * 8 sh,xt
      real * 8 epart1,epart2,ptmin,xmt,xmtmax,xmtmin,
     # ey1max,y1max,ypart,ycm,ey
      data ptmin/0.05d0/
      if(jptseqfl.ne.0) then
         if(iflagmode.eq.0) then
            pt=z
         else
            z=pt
         endif
         return
      endif
c the factorization scale is xmt*ffact. With this choice
c of xmtmin it is never below sqrt(q2min) (to stay within pdf
c boundaries) and always above xm (to avoid the coulomb singularity)
c zmax1, zmax2 are the maximum energy fraction allowed in beam1, beam2.
c In practice, if beam1(2) is an electron, zmax1(2) is the upper bound
c of the WW variable z. In all other cases they are 1.
      epart1=ebeam1*zmax1
      epart2=ebeam2*zmax2
      xmtmin=max(sqrt(ptmin**2+xm**2),sqrt(q2min)/ffact)
      sh=(epart2+epart1)**2-(epart2-epart1)**2
      ypart=0.5*log(epart1/epart2)
      xmt=xm
c Upper bound for y; obtained by imposing E_quark<sqrt(shat)/2
c shat is the CM energy of the part1-part2 system
      xt=2*xmt/sqrt(sh)
      ey1max=1/xt*(1+sqrt(1-xt**2))
      y1max=log(ey1max)
c      write(*,*) '( y limit,',-y1max+ypart,'<y<',y1max+ypart
      ycm=y-ypart
      if(abs(ycm).lt.y1max) then
c generate a pt sequence
         ey=exp(ycm)
c Upper limit for m_T, obtained by imposing E_quark<sqrt(shat)/2
         xmtmax=sqrt(sh)/(2*cosh(ycm))
         xmtmax=min( xmtmax,sqrt(ptmax**2+xm**2) )
         if(xmtmax.le.xmtmin)then
            write(*,*)'y is out of range'
            stop
         else
            if(iflagmode.eq.0) then
               xmt=xmtmin*exp(z*log(xmtmax/xmtmin))
               pt=sqrt(xmt**2-xm**2)
            else
               xmt=sqrt(pt**2+xm**2)
               z=log(xmt/xmtmin)/log(xmtmax/xmtmin)
            endif
         endif
      else
         write(*,*)'y is out of range'
         stop
      endif
      end
