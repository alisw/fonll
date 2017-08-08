c This file contains all the relevant cross section formulae for heavy
c quark photoproduction.
c Common blocks :
c    /process/prc           (character * 2)
c    /nl/nl
c    /schhad/schhad         (character * 2)
c    /schpho/schpho         (character * 2)
c    /betfac/betfac,delta
c Functions called:
c    ddilog  (cernlib)
c
      subroutine ppsv(sveh2,svel2,sveleh,s,t,xm2,xmur2,xmufph2,xmufh2)
      implicit real * 8 (a-h,o-z)
      character * 2 prc
      common/process/prc
      common/nl/nl
      if(prc.eq.'pg')then
        csir =  (33.d0-2.d0*nl)/3 * log(xmur2/xmufh2)
        sveh2 = pgqq2(s,t,xm2,xmufh2,nl)+csir*pgborn(s,t,xm2)
        svel2 = 0
        sveleh = 0
      elseif(prc.eq.'pq')then
        call pqqq2(sveh2,svel2,sveleh,s,t,xm2,nl)
      else
        write(*,*)'PPSV: non existent process ',prc
        stop
      endif
      end

      function ppcolp(y,s,q1q,x,xm2,xlmude)
      implicit real * 8 (a-h,o-z)
      character * 2 prc
      common/nl/nl
      common/process/prc
      if( y .eq. 1 .and. prc.eq.'pq' ) then
         ppcolp = pqcolp1(s,q1q,x,xm2,xlmude,nl)
      elseif( y .eq. -1) then
         if( prc.eq.'pg' ) then
            ppcolp = pgcolp2(s,q1q,x,xm2,xlmude,nl)
         elseif( prc.eq.'pq' ) then
            ppcolp = pqcolp2(s,q1q,x,xm2,xlmude,nl)
         endif
      else
         write(*,*) 'error in ppcolp: prc=',prc,' y=',y
         stop
      endif
      end

      function ppcoll(y,s,q1q,x,xm2)
      implicit real * 8 (a-h,o-z)
      character * 2 prc
      common/nl/nl
      common/process/prc
      if( y .eq. 1 .and. prc.eq.'pq' ) then
         ppcoll = pqcoll1(s,q1q,x,xm2,nl)
      elseif( y .eq. -1) then
         if( prc.eq.'pg' ) then
            ppcoll = pgcoll2(s,q1q,x,xm2,nl)
         elseif( prc.eq.'pq' ) then
            ppcoll = pqcoll2(s,q1q,x,xm2,nl)
         endif
      else
         write(*,*) 'error in ppcoll: prc=',prc,' y=',y
         stop
      endif
      end

      subroutine fpp
     # (feh2,fel2,feleh,s0,x,y,xm20,q1q0,q2q0,w1h,w2h,cth2)
      implicit real * 8 (a-h,o-z)
      character * 2 prc
      common/process/prc
      s=1
      xm2=xm20/s0
      q1q=q1q0/s0
      q2q=q2q0/s0
      if(prc.eq.'pg') then
         feh2 = fpg(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
         fel2 = 0
         feleh = 0
      elseif(prc.eq.'pq') then
         call fpq(feh2,fel2,feleh,
     #            s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
      else
         write(*,*) 'FPP: non existent subprocess',prc
         stop
      endif
      end

c 
c Contributo di born al processo p + g ----> Q + Qbar + X
c
c
c         (born)    
c  d sigma       =  g^2 (e_h)^2 pgborn dphi_2
c         pg        
c
      function pgborn(s,t,m2)
      implicit real * 8 (a-z)

      vtf = 1/2.d0
      u = -t-s
      born = 2*(u**2+4*m2*s*(1-m2*s/(t*u))+t**2)*vtf/(t*u)
      pgborn = born/(2*s)
      return 
      end 
c
c 
c Contributo di born al processo q + qbar ----> Q + Qbar + X
c
c
c         (born)    
c  d sigma       =  g^4 bornqq dphi_2
c         qq        
c
c
      function bornqq(s,t,m2)
      implicit real * 8 (a-z)

      vtf = 1/2.d0
      vda = 8
      vdf = 3
      u = -t-s
      born = (8*s**2-16*(t*u-m2*s))*vda*vtf**2/(s**2*vdf**2*4.d0)
      bornqq = born/(2*s)
      return 
      end 
c
c
c
c Contributo collineare al processo p + g ----> Q + Qbar + X
c
c All'ordine minimo X = g e si ha la fattorizzazione della funzione
c di splitting P_gg(x) di Altarelli Parisi
c
c d sigma_pg (c-) = g^4 (e_h)^2 [ 1/(1-x)_rho pgcolp
c                      +        [log(1-x)/(1-x)]_rho * pgcoll ] d phi_2x
c
c xlmude = log(s/xmu2) + log(delta/2)
c
c
       function pgcolp2(s,q1q,x,m2,xlmude,nl)

       implicit real * 8 (a-z)
       integer nl
       character * 2 schhad
       common/schhad/schhad
       data pi/3.141 592 653 589 793/

       vca = 3
       pgcolp2 = vca/(4*pi**2)*xlmude*(x+(1-x)**2/x+x*(1-x)**2)
       if(schhad.eq.'DI')then
           pgcolp2 = pgcolp2 - xkpgg(x,nl)/(8*pi**2)
       elseif(schhad.ne.'MS')then
           write(6,*)'scheme ',schhad,'not known'
           stop
       endif
       pgcolp2 = pgcolp2*pgborn(x*s,q1q,m2)
       return
       end


       function pgcoll2(s,q1q,x,m2,nl)

       implicit real * 8 (a-z)
       integer nl
       character * 2 schhad
       common/schhad/schhad
       data pi/3.141 592 653 589 793/

       vca = 3
       pgcoll2 = 2*vca/(4*pi**2)*(x+(1-x)**2/x+x*(1-x)**2)
       if(schhad.eq.'DI')then
           pgcoll2 = pgcoll2 - xklgg(x,nl)/(8*pi**2)
       elseif(schhad.ne.'MS')then
           write(6,*)'scheme ',schhad,'not known'
           stop
       endif
       pgcoll2 = pgcoll2*pgborn(x*s,q1q,m2)
       return
       end
c
c
c Contributo collineare al processo p + q ----> Q + Qbar + X
c
c All'ordine minimo X = q e si ha la fattorizzazione della funzione
c di splitting P_qp(x) di Altarelli Parisi
c
c d sigma_pq (c+) = g^4 (e_l)^2 [ 1/(1-x)_rho pqcolp1
c                     +         [log(1-x)/(1-x)]_rho * pqcoll1 ] d phi_2x
c
c xlmude = log(s/xmu2) + log(delta/2)
c
c
       function pqcolp1(s,q2q,x,m2,xlmude,nl)

       implicit real * 8 (a-z)
       integer nl
       character * 2 schpho
       common/schpho/schpho
       data pi/3.141 592 653 589 793/

       vdf = 3
       vnc = 3
       vtf = 1/2.d0
       pqcolp1 = vdf/(8*pi**2) *
     #   ( xlmude*(x**2 + (1-x)**2) + 2*x*(1-x) )*(1-x)
       if(schpho.eq.'DI')then
c should be xkpqp = vnc/vtf xkpqg
           pqcolp1 = pqcolp1 - vnc/vtf*xkpqg(x,nl)/(8*pi**2)
       elseif(schpho.ne.'MS')then
           write(6,*)'scheme ',schpho,'not known'
           stop
       endif
       pqcolp1 = pqcolp1 * bornqq(s*x,q2q,m2)
       return
       end


       function pqcoll1(s,q2q,x,m2,nl)

       implicit real * 8 (a-z)
       integer nl
       character * 2 schpho
       common/schpho/schpho
       data pi/3.141 592 653 589 793/

       vdf = 3
       vnc = 3
       vtf = 1/2.d0
       pqcoll1 = vdf/(8*pi**2) *
     #   2*(x**2 + (1-x)**2) * (1-x)
       if(schpho.eq.'DI')then
c should be xklqp = vnc/vtf xklqg
           pqcoll1 = pqcoll1 - vnc/vtf*xklqg(x,nl)/(8*pi**2)
       elseif(schpho.ne.'MS')then
           write(6,*)'scheme ',schpho,'not known'
           stop
       endif
       pqcoll1 = pqcoll1 * bornqq(s*x,q2q,m2)
       return
       end

c
c Contributo collineare al processo p + q ----> Q + Qbar + X
c
c All'ordine minimo X = q e si ha la fattorizzazione della funzione
c di splitting P_gq(x) di Altarelli Parisi
c
c d sigma_pq (c-) = g^4 (e_h)^2 [ 1/(1-x)_rho pqcolp2
c                     +         [log(1-x)/(1-x)]_rho * pqcoll2 ] d phi_2x
c
c xlmude = log(s/xmu2) + log(delta/2)
c
c
       function pqcolp2(s,q1q,x,m2,xlmude,nl)

       implicit real * 8 (a-z)
       integer nl
       character * 2 schhad
       common/schhad/schhad
       data pi/3.141 592 653 589 793/

       cf = 4.d0/3.d0
       pqcolp2 = cf/(8*pi**2) *
     #   ( xlmude*(1 + (1-x)**2)/x + x )*(1-x)
       if(schhad.eq.'DI')then
           pqcolp2 = pqcolp2 - xkpgq(x,nl)/(8*pi**2)
       elseif(schhad.ne.'MS')then
           write(6,*)'scheme ',schhad,'not known'
           stop
       endif
       pqcolp2 = pqcolp2 * pgborn(s*x,q1q,m2)
       return
       end


       function pqcoll2(s,q1q,x,m2,nl)

       implicit real * 8 (a-z)
       integer nl
       character * 2 schhad
       common/schhad/schhad
       data pi/3.141 592 653 589 793/

       cf = 4.d0/3.d0
       pqcoll2 = cf/(8*pi**2) *
     #   2*(1 + (1-x)**2)/x * (1-x)
       if(schhad.eq.'DI')then
           pqcoll2 = pqcoll2 - xklgq(x,nl)/(8*pi**2)
       elseif(schhad.ne.'MS')then
           write(6,*)'scheme ',schhad,'not known'
           stop
       endif
       pqcoll2 = pqcoll2 * pgborn(s*x,q1q,m2)
       return
       end
c
c 
c Contributo due corpi al processo  p + q ----> Q + Qbar + X
c
c Effetto spurio dovuto unicamente al cambio di schema
c
c
      subroutine pqqq2(sveh2,svel2,sveleh,s,t,m2,nl)
      implicit double precision (a-z)
      integer nl
      character * 2 schpho
      character * 2 schhad
      common/betfac/betfac,delta
      common/schpho/schpho
      common/schhad/schhad
      data pi/3.141 592 653 589 793/
                                            
      if(      (schhad.ne.'MS'.and.schhad.ne.'DI')
     #     .or.(schpho.ne.'MS'.and.schpho.ne.'DI')  )  then
           write(6,*)'pqqq2: scheme ',schhad,schpho,'not known'
           stop
      endif
      sveh2 = 0
      svel2 = 0
      sveleh = 0
      vnc = 3
      vtf = 1/2.d0
      ro = 4*m2/s
      b = dsqrt(1-ro)
      lb = log(b*betfac)
      if(schhad.eq.'DI')then
           one = 1
           xk2 = xkdgq(nl) + 2*xkpgq(one,nl)*lb + 2*xklgq(one,nl)*lb**2
           sveh2 = -xk2*16*pi**2*pgborn(s,t,m2)/(8*pi**2)
      endif
      if(schpho.eq.'DI')then
c should be xkdqp = vnc/vtf xkdqg
           one = 1
           xk1 = xkdqg(nl) + 2*xkpqg(one,nl)*lb + 2*xklqg(one,nl)*lb**2
           svel2 = -xk1*vnc/vtf*16*pi**2*bornqq(s,t,m2)/(8*pi**2)
      endif
      return
      end  
c
c
c Contributo soft virtuale al processo  p + g ----> Q + Qbar + X
c
c         (sv)   g^4 (e_h)^2
c  d sigma    = -------------   pgqq2 dphi_2
c         pg      (4 pi)^2
c
c
      function pgqq2(s,t,m2,mu2,nl)
      implicit double precision (a-z)
      integer nl
      character * 2 schhad
      common/betfac/betfac,delta
      common/schhad/schhad
      data pi/3.141 592 653 589 793/

      ro = 4*m2/s
      t1 = -t/s
      zg = 1
      zeh = 1
      vca = 3
      vcf = 4
      vcf = vcf/ 3.d0
      vtf = 1
      vtf = vtf/ 2.d0
      nlf = nl
      t2 = 1-t1
      b = dsqrt(1-ro)
      lp = (b+1)/ 2.d0
      lm = (1-b)/ 2.d0
      at = s*t1
      aw = s*t2
      vlm2 = dlog(m2/mu2)
      vltm = dlog(at/m2)
      vlpm = dlog(lp/lm)
      vlsm = dlog(s/m2)
      vlsmu = dlog(s/mu2)
      vlwm = dlog(aw/m2)
      vlbl = dlog(b/lm)
      vdw = ddilog((aw-m2)/aw)-vlwm**2/ 2.d0
      vdt = ddilog((at-m2)/at)-vltm**2/ 2.d0
      vdmp = ddilog(-lm/lp)
      vdmb = vlbl**2/ 2.d0+ddilog(-lm/b)
      auinv = 1/(m2-aw)
      atinv = 1/(m2-at)
      softt1 = ddilog(1-2*t1/(b+1))+ddilog(1-2*t1/(1-b))+log(2*t1/(1-b))
     1   *log(2*t1/(b+1))
      softt2 = ddilog(1-2*t2/(b+1))+ddilog(1-2*t2/(1-b))+log(2*t2/(1-b))
     1   *log(2*t2/(b+1))
      softb = ddilog(2*b/(b+1))-ddilog(-2*b/(1-b))
      lt1 = log(t1)
      lt2 = log(t2)
      lb = log(b*betfac)
      ss = -2*lb*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**
     1   2)*vca*vlsmu*vtf*zeh**2*zg**4/(s*(t1-1)**2*t1**2)
      ss = (8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2)*vls
     1   mu*vtf*(4*nlf*vtf-11*vca)*zeh**2*zg**4/( 12.d0*s*(t1-1)**2*t1**
     2   2)+ss
      ss = ss-3*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2
     1   )*vca*vlsm**2*vtf*zeh**2*zg**4/( 4.d0*s*(t1-1)**2*t1**2)
      ss = (ro-2)*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro*
     1   *2)*(2*vcf-vca)*vlpm*vlsm*vtf*zeh**2*zg**4/( 4.d0*b*s*(t1-1)**2
     2   *t1**2)+ss
      ss = ss-lt2*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro*
     1   *2)*vca*vlsm*vtf*zeh**2*zg**4/( 2.d0*s*(t1-1)**2*t1**2)
      ss = ss-lt1*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro*
     1   *2)*vca*vlsm*vtf*zeh**2*zg**4/( 2.d0*s*(t1-1)**2*t1**2)
      ss = ss-2*lb*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro
     1   **2)*vca*vlsm*vtf*zeh**2*zg**4/(s*(t1-1)**2*t1**2)
      ss = ss-vlsm*vtf*(32*nlf*t1**4*vtf-64*nlf*t1**3*vtf+16*nlf*ro*t1**
     1   2*vtf+48*nlf*t1**2*vtf-16*nlf*ro*t1*vtf-16*nlf*t1*vtf+4*nlf*ro*
     2   *2*vtf-96*t1**4*vcf+192*t1**3*vcf-48*ro*t1**2*vcf-144*t1**2*vcf
     3   +48*ro*t1*vcf+48*t1*vcf-12*ro**2*vcf-184*t1**4*vca+368*t1**3*vc
     4   a-140*ro*t1**2*vca-228*t1**2*vca+140*ro*t1*vca+44*t1*vca-35*ro*
     5   *2*vca)*zeh**2*zg**4/( 12.d0*s*(t1-1)**2*t1**2)
      ss = ss-(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2)*
     1   vca*vlpm**2*vtf*zeh**2*zg**4/( 4.d0*s*(t1-1)**2*t1**2)
      ss = lb*(ro-2)*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+
     1   ro**2)*(2*vcf-vca)*vlpm*vtf*zeh**2*zg**4/(b*s*(t1-1)**2*t1**2)+
     2   ss
      ss = ss-(8*ro*t1**4*vcf-16*ro*t1**3*vcf+8*ro**2*t1**2*vcf+8*t1**2*
     1   vcf-8*ro**2*t1*vcf+8*ro*t1*vcf-8*t1*vcf+2*ro**3*vcf-2*ro**2*vcf
     2   -4*ro*t1**4*vca+8*t1**4*vca+8*ro*t1**3*vca-16*t1**3*vca-4*ro**2
     3   *t1**2*vca+4*ro*t1**2*vca+8*t1**2*vca+4*ro**2*t1*vca-8*ro*t1*vc
     4   a-ro**3*vca+2*ro**2*vca)*vlpm*vtf*zeh**2*zg**4/( 2.d0*b*s*(t1-1
     5   )**2*t1**2)
      ss = ss-softt2*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+
     1   ro**2)*vca*vtf*zeh**2*zg**4/( 2.d0*s*(t1-1)**2*t1**2)
      ss = ss-softt1*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+
     1   ro**2)*vca*vtf*zeh**2*zg**4/( 2.d0*s*(t1-1)**2*t1**2)
      ss = ss-(ro-2)*softb*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1
     1   -4*t1+ro**2)*(2*vcf-vca)*vtf*zeh**2*zg**4/( 4.d0*b*s*(t1-1)**2*
     2   t1**2)
      ss = pi**2*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**
     1   2)*vca*vtf*zeh**2*zg**4/( 6.d0*s*(t1-1)**2*t1**2)+ss
      ss = ss-2*lb*lt2*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t
     1   1+ro**2)*vca*vtf*zeh**2*zg**4/(s*(t1-1)**2*t1**2)
      ss = lt2*(2*t1**2-2*t1+ro)**2*vca*vtf*zeh**2*zg**4/(s*(t1-1)**2*t1
     1   **2)+ss
      ss = ss-2*lb*lt1*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t
     1   1+ro**2)*vca*vtf*zeh**2*zg**4/(s*(t1-1)**2*t1**2)
      ss = lt1*(2*t1**2-2*t1+ro)**2*vca*vtf*zeh**2*zg**4/(s*(t1-1)**2*t1
     1   **2)+ss
      ss = ss-4*lb**2*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1
     1   +ro**2)*vca*vtf*zeh**2*zg**4/(s*(t1-1)**2*t1**2)
      ss = 4*lb*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**2
     1   )*vcf*vtf*zeh**2*zg**4/(s*(t1-1)**2*t1**2)+ss
      ss = 2*vtf*(4*nlf*t1**4*vtf-8*nlf*t1**3*vtf+4*nlf*ro*t1**2*vtf+4*n
     1   lf*t1**2*vtf-4*nlf*ro*t1*vtf+nlf*ro**2*vtf-12*t1**4*vcf+24*t1**
     2   3*vcf-12*ro*t1**2*vcf-12*t1**2*vcf+12*ro*t1*vcf-3*ro**2*vcf-17*
     3   t1**4*vca+34*t1**3*vca-20*ro*t1**2*vca-17*t1**2*vca+20*ro*t1*vc
     4   a-5*ro**2*vca)*zeh**2*zg**4/( 3.d0*s*(t1-1)**2*t1**2)+ss
      dd = -auinv*pi*(4*t1+ro-4)*(8*t1**4-16*t1**3+12*ro*t1**2+8*t1**2-1
     1   2*ro*t1+3*ro**2)*(4*vcf-vca)*vlwm**2*vtf*zeh**2*zg**4/( 8.d0*(t
     2   1-1)**3*t1)
      dd = pi*(16*ro*t1**3*vcf+32*t1**3*vcf+8*ro**2*t1**2*vcf-16*ro*t1**
     1   2*vcf-96*t1**2*vcf-16*ro**2*t1*vcf+8*ro*t1*vcf+96*t1*vcf+12*ro*
     2   *2*vcf-8*ro*vcf-32*vcf-4*t1**4*vca-8*ro*t1**3*vca-4*ro**2*t1**2
     3   *vca+16*ro*t1**2*vca+24*t1**2*vca+8*ro**2*t1*vca-14*ro*t1*vca-3
     4   2*t1*vca-5*ro**2*vca+6*ro*vca+12*vca)*vlwm**2*vtf*zeh**2*zg**4/
     5   (s*(t1-1)**3*t1)+dd
      dd = 4*pi*(4*t1**4-8*t1**3+ro**2*t1**2+2*ro*t1**2+6*t1**2-ro**2*t1
     1   -2*ro*t1-2*t1+ro**2)*vca*vltm*vlwm*vtf*zeh**2*zg**4/(s*(t1-1)**
     2   2*t1**2)+dd
      dd = 2*pi*(4*t1**3+2*ro**2*t1**2+2*ro*t1**2-12*t1**2-6*ro**2*t1-2*
     1   ro*t1+16*t1-ro**3+6*ro**2-8)*(2*vcf-vca)*vlpm*vlwm*vtf*zeh**2*z
     2   g**4/(b*s*(t1-1)**2*t1)+dd
      dd = auinv**2*pi*s*(192*t1**5*vcf+32*ro*t1**4*vcf-512*t1**4*vcf+19
     1   2*ro*t1**3*vcf+512*t1**3*vcf+48*ro**2*t1**2*vcf-384*ro*t1**2*vc
     2   f-384*t1**2*vcf+24*ro**2*t1*vcf+64*ro*t1*vcf+320*t1*vcf+12*ro**
     3   3*vcf-72*ro**2*vcf+96*ro*vcf-128*vcf-32*t1**5*vca-8*ro*t1**4*vc
     4   a-32*t1**4*vca-64*ro*t1**3*vca+256*t1**3*vca-12*ro**2*t1**2*vca
     5   +96*ro*t1**2*vca-256*t1**2*vca-12*ro**2*t1*vca+32*ro*t1*vca+32*
     6   t1*vca-3*ro**3*vca+24*ro**2*vca-56*ro*vca+32*vca)*vlwm*vtf*zeh*
     7   *2*zg**4/( 8.d0*(t1-1)**2*t1)+dd
      dd = dd-auinv*pi*(32*t1**4*vcf+16*ro*t1**3*vcf-40*t1**3*vcf+4*ro**
     1   2*t1**2*vcf+4*ro*t1**2*vcf+4*ro**2*t1*vcf-20*ro*t1*vcf-8*t1*vcf
     2   +2*ro**3*vcf-4*ro**2*vcf+16*vcf+4*t1**4*vca-32*t1**3*vca-2*ro*t
     3   1**2*vca+48*t1**2*vca-6*ro*t1*vca-16*t1*vca-ro**2*vca+8*ro*vca-
     4   4*vca)*vlwm*vtf*zeh**2*zg**4/((t1-1)**2*t1)
      dd = dd-4*pi*(2*t1**2-2*t1+ro)**2*vca*vlwm*vtf*zeh**2*zg**4/(s*(t1
     1   -1)**2*t1**2)
      dd = atinv*pi*(4*t1-ro)*(8*t1**4-16*t1**3+12*ro*t1**2+8*t1**2-12*r
     1   o*t1+3*ro**2)*(4*vcf-vca)*vltm**2*vtf*zeh**2*zg**4/( 8.d0*(t1-1
     2   )*t1**3)+dd
      dd = dd-pi*(16*ro*t1**3*vcf+32*t1**3*vcf-8*ro**2*t1**2*vcf-32*ro*t
     1   1**2*vcf+24*ro*t1*vcf-4*ro**2*vcf+4*t1**4*vca-8*ro*t1**3*vca-16
     2   *t1**3*vca+4*ro**2*t1**2*vca+8*ro*t1**2*vca-6*ro*t1*vca+ro**2*v
     3   ca)*vltm**2*vtf*zeh**2*zg**4/(s*(t1-1)*t1**3)
      dd = 2*pi*(4*t1**3-2*ro**2*t1**2-2*ro*t1**2-2*ro**2*t1+2*ro*t1+4*t
     1   1+ro**3-2*ro**2)*(2*vcf-vca)*vlpm*vltm*vtf*zeh**2*zg**4/(b*s*(t
     2   1-1)*t1**2)+dd
      dd = atinv**2*pi*s*(192*t1**5*vcf-32*ro*t1**4*vcf-448*t1**4*vcf+32
     1   0*ro*t1**3*vcf+384*t1**3*vcf-48*ro**2*t1**2*vcf-384*ro*t1**2*vc
     2   f+120*ro**2*t1*vcf-12*ro**3*vcf-32*t1**5*vca+8*ro*t1**4*vca+128
     3   *t1**4*vca-80*ro*t1**3*vca-64*t1**3*vca+12*ro**2*t1**2*vca+80*r
     4   o*t1**2*vca-28*ro**2*t1*vca+3*ro**3*vca)*vltm*vtf*zeh**2*zg**4/
     5   ( 8.d0*(t1-1)*t1**2)+dd
      dd = atinv*pi*(32*t1**4*vcf-16*ro*t1**3*vcf-88*t1**3*vcf+4*ro**2*t
     1   1**2*vcf+52*ro*t1**2*vcf+72*t1**2*vcf-12*ro**2*t1*vcf-36*ro*t1*
     2   vcf+2*ro**3*vcf+4*ro**2*vcf+4*t1**4*vca+8*t1**3*vca-2*ro*t1**2*
     3   vca-8*t1**2*vca+6*ro*t1*vca-ro**2*vca)*vltm*vtf*zeh**2*zg**4/((
     4   t1-1)*t1**2)+dd
      dd = dd-4*pi*(2*t1**2-2*t1+ro)**2*vca*vltm*vtf*zeh**2*zg**4/(s*(t1
     1   -1)**2*t1**2)
      dd = dd-pi*(8*t1**4-16*t1**3-6*ro**2*t1**2+2*ro*t1**2+20*t1**2+6*r
     1   o**2*t1-2*ro*t1-12*t1-ro**3+2*ro**2)*(2*vcf-vca)*vlpm*vlsm*vtf*
     2   zeh**2*zg**4/(b*s*(t1-1)**2*t1**2)
      dd = pi*(ro-1)*(4*t1**2-4*t1-ro-2)*(2*vcf-vca)*vlpm**2*vtf*zeh**2*
     1   zg**4/(b*s*(t1-1)*t1)+dd
      dd = dd-pi*(4*t1**2-4*t1-ro**2+6)*(2*vcf-vca)*vlpm**2*vtf*zeh**2*z
     1   g**4/(s*(t1-1)*t1)
      dd = 4*b*pi*(2*vcf-vca)*vlpm*vtf*zeh**2*zg**4/(s*(t1-1)*t1)+dd
      dd = 2*pi*(4*ro*t1**4-16*t1**4-8*ro*t1**3+32*t1**3+4*ro**2*t1**2-4
     1   *ro*t1**2-18*t1**2-4*ro**2*t1+8*ro*t1+2*t1+ro**3-2*ro**2)*(2*vc
     2   f-vca)*vlpm*vtf*zeh**2*zg**4/(b*s*(t1-1)**2*t1**2)+dd
      dd = dd-pi*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-4*t1+ro**
     1   2)*vlm2*vtf*(4*nlf*vtf-11*vca)*zeh**2*zg**4/( 3.d0*s*(t1-1)**2*
     2   t1**2)
      tmp0 = -auinv*pi*(4*t1+ro-4)*(32*t1**4*vcf-64*t1**3*vcf+48*ro*t1**
     1   2*vcf+32*t1**2*vcf-48*ro*t1*vcf+12*ro**2*vcf-8*t1**4*vca+16*t1*
     2   *3*vca-12*ro*t1**2*vca+12*ro*t1*vca-16*t1*vca-3*ro**2*vca+8*vca
     3   )*vdw*vtf*zeh**2*zg**4/( 8.d0*(t1-1)**3*t1)
      dd = tmp0+dd
      dd = pi*(16*ro*t1**3*vcf+32*t1**3*vcf+8*ro**2*t1**2*vcf-16*ro*t1**
     1   2*vcf-96*t1**2*vcf-16*ro**2*t1*vcf+8*ro*t1*vcf+96*t1*vcf+12*ro*
     2   *2*vcf-8*ro*vcf-32*vcf-12*t1**4*vca-8*ro*t1**3*vca+24*t1**3*vca
     3   -2*ro**2*t1**2*vca+12*ro*t1**2*vca+4*ro**2*t1*vca-6*ro*t1*vca-2
     4   4*t1*vca-3*ro**2*vca+2*ro*vca+12*vca)*vdw*vtf*zeh**2*zg**4/(s*(
     5   t1-1)**3*t1)+dd
      dd = atinv*pi*(4*t1-ro)*(8*t1**4-16*t1**3+12*ro*t1**2+8*t1**2-12*r
     1   o*t1+3*ro**2)*(4*vcf-vca)*vdt*vtf*zeh**2*zg**4/( 8.d0*(t1-1)*t1
     2   **3)+dd
      dd = dd-pi*(16*ro*t1**3*vcf+32*t1**3*vcf-8*ro**2*t1**2*vcf-32*ro*t
     1   1**2*vcf+24*ro*t1*vcf-4*ro**2*vcf+12*t1**4*vca-8*ro*t1**3*vca-2
     2   4*t1**3*vca+2*ro**2*t1**2*vca+12*ro*t1**2*vca+4*t1**2*vca-6*ro*
     3   t1*vca+ro**2*vca)*vdt*vtf*zeh**2*zg**4/(s*(t1-1)*t1**3)
      dd = 2*pi*(8*t1**4-16*t1**3-6*ro**2*t1**2+2*ro*t1**2+20*t1**2+6*ro
     1   **2*t1-2*ro*t1-12*t1-ro**3+2*ro**2)*(2*vcf-vca)*vdmp*vtf*zeh**2
     2   *zg**4/(b*s*(t1-1)**2*t1**2)+dd
      dd = dd-2*pi*(ro-2)*(8*t1**4-16*t1**3+4*ro*t1**2+12*t1**2-4*ro*t1-
     1   4*t1+ro**2)*(2*vcf-vca)*vdmb*vtf*zeh**2*zg**4/(b*s*(t1-1)**2*t1
     2   **2)
      dd = pi**3*(32*ro*t1**4-56*t1**4-64*ro*t1**3+112*t1**3+10*ro**2*t1
     1   **2+18*ro*t1**2-76*t1**2-10*ro**2*t1+14*ro*t1+20*t1+3*ro**3-6*r
     2   o**2)*(2*vcf-vca)*vtf*zeh**2*zg**4/( 6.d0*b*s*(t1-1)**2*t1**2)+
     3   dd
      tmp0 = 128*pi**2*t1**5*vcf+384*t1**5*vcf+32*pi**2*ro*t1**4*vcf-384
     1   *pi**2*t1**4*vcf+384*t1**4*vcf+128*pi**2*ro*t1**3*vcf
      tmp0 = 480*ro*t1**3*vcf+384*pi**2*t1**3*vcf-1920*t1**3*vcf+48*pi**
     1   2*ro**2*t1**2*vcf-48*ro**2*t1**2*vcf-352*pi**2*ro*t1**2*vcf+576
     2   *ro*t1**2*vcf-128*pi**2*t1**2*vcf-384*t1**2*vcf+384*ro**2*t1*vc
     3   f+192*pi**2*ro*t1*vcf-2592*ro*t1*vcf+3072*t1*vcf+12*pi**2*ro**3
     4   *vcf-48*pi**2*ro**2*vcf-336*ro**2*vcf+1536*ro*vcf-1536*vcf-32*p
     5   i**2*t1**5*vca-8*pi**2*ro*t1**4*vca+96*pi**2*t1**4*vca-1152*t1*
     6   *4*vca-32*pi**2*ro*t1**3*vca-288*ro*t1**3*vca-96*pi**2*t1**3*vc
     7   a+3072*t1**3*vca-12*pi**2*ro**2*t1**2*vca+88*pi**2*ro*t1**2*vca
     8   +192*ro*t1**2*vca+32*pi**2*t1**2*vca-2304*t1**2*vca-72*ro**2*t1
     9   *vca-48*pi**2*ro*t1*vca+480*ro*t1*vca-3*pi**2*ro**3*vca+12*pi**
     :   2*ro**2*vca+72*ro**2*vca-384*ro*vca+384*vca+tmp0
      tmp0 = -auinv*pi*tmp0*vtf*zeh**2*zg**4/( 48.d0*(t1-1)**3*t1)
      dd = tmp0+dd
      dd = atinv*pi*(128*pi**2*t1**5*vcf+384*t1**5*vcf-32*pi**2*ro*t1**4
     1   *vcf-256*pi**2*t1**4*vcf-2304*t1**4*vcf+256*pi**2*ro*t1**3*vcf+
     2   480*ro*t1**3*vcf+128*pi**2*t1**3*vcf+3456*t1**3*vcf-48*pi**2*ro
     3   **2*t1**2*vcf+48*ro**2*t1**2*vcf-224*pi**2*ro*t1**2*vcf-2016*ro
     4   *t1**2*vcf+96*pi**2*ro**2*t1*vcf+288*ro**2*t1*vcf-12*pi**2*ro**
     5   3*vcf-32*pi**2*t1**5*vca+8*pi**2*ro*t1**4*vca+64*pi**2*t1**4*vc
     6   a+1152*t1**4*vca-64*pi**2*ro*t1**3*vca-288*ro*t1**3*vca-32*pi**
     7   2*t1**3*vca-768*t1**3*vca+12*pi**2*ro**2*t1**2*vca+56*pi**2*ro*
     8   t1**2*vca+480*ro*t1**2*vca-24*pi**2*ro**2*t1*vca-72*ro**2*t1*vc
     9   a+3*pi**2*ro**3*vca)*vtf*zeh**2*zg**4/( 48.d0*(t1-1)*t1**3)+dd
      tmp0 = 64*nlf*t1**6*vtf-192*nlf*t1**5*vtf+64*nlf*ro*t1**4*vtf+192*
     1   nlf*t1**4*vtf-128*nlf*ro*t1**3*vtf-64*nlf*t1**3*vtf+16*nlf*ro**
     2   2*t1**2*vtf+64*nlf*ro*t1**2*vtf-16*nlf*ro**2*t1*vtf-48*pi**2*t1
     3   **6*vcf-576*t1**6*vcf+144*pi**2*t1**5*vcf+1728*t1**5*vcf-4*pi**
     4   2*ro**2*t1**4*vcf-48*pi**2*ro*t1**4*vcf-288*ro*t1**4*vcf-184*pi
     5   **2*t1**4*vcf-2496*t1**4*vcf+8*pi**2*ro**2*t1**3*vcf+96*pi**2*r
     6   o*t1**3*vcf
      tmp0 = 576*ro*t1**3*vcf+128*pi**2*t1**3*vcf+2112*t1**3*vcf-12*pi**
     1   2*ro**2*t1**2*vcf-96*ro**2*t1**2*vcf-72*pi**2*ro*t1**2*vcf-384*
     2   ro*t1**2*vcf-40*pi**2*t1**2*vcf-768*t1**2*vcf+8*pi**2*ro**2*t1*
     3   vcf+96*ro**2*t1*vcf+24*pi**2*ro*t1*vcf+96*ro*t1*vcf-4*pi**2*ro*
     4   *2*vcf+64*pi**2*t1**6*vca+16*t1**6*vca-192*pi**2*t1**5*vca-48*t
     5   1**5*vca+10*pi**2*ro**2*t1**4*vca+24*pi**2*ro*t1**4*vca-200*ro*
     6   t1**4*vca+232*pi**2*t1**4*vca+336*t1**4*vca-20*pi**2*ro**2*t1**
     7   3*vca-48*pi**2*ro*t1**3*vca+400*ro*t1**3*vca-144*pi**2*t1**3*vc
     8   a-592*t1**3*vca+20*pi**2*ro**2*t1**2*vca-44*ro**2*t1**2*vca+30*
     9   pi**2*ro*t1**2*vca-176*ro*t1**2*vca+40*pi**2*t1**2*vca+288*t1**
     :   2*vca-10*pi**2*ro**2*t1*vca+44*ro**2*t1*vca-6*pi**2*ro*t1*vca-2
     ;   4*ro*t1*vca+pi**2*ro**2*vca+tmp0
      tmp0 = -pi*tmp0*vtf*zeh**2*zg**4/( 6.d0*s*(t1-1)**3*t1**3)
      dd = tmp0+dd
c Correction from difference with RK Ellis
      pnmrk = 8*vca*vtf*zeh**2*zg**4
     # *pi*3*(t1**2+t2**2+ro*(1-ro/(4*t1*t2)))/(s*t1*t2)
      dd = dd - pnmrk
      pgqq2 = ss+dd/( 4.d0*pi)
      if(schhad.eq.'DI')then
           one = 1
           xk = xkdgg(nl) + 2*xkpgg(one,nl)*lb + 2*xklgg(one,nl)*lb**2
           xk = xk*16*pi**2         
           pgqq2 = pgqq2 - xk*pgborn(s,t,m2)/(8*pi**2)
      elseif(schhad.ne.'MS')then
           write(6,*)'scheme ',schhad,'not known'
           stop
      endif
      return
      end
c
       function fpg(s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)
c
c      d sigma_pg (f) = N g^4 eh^2 / (64 pi^2 s) beta_x d costh1 d th2 dy dx
c                    1/(1-x)_rho ( 1/(1-y)_+ + 1/(1+y)_+ ) fpg
c
c           N = 1 / (4 pi)^2
c
       implicit real * 8(a-h,o-z)
       tiny = 0.1d-10
       tf = 1/2.d0
       cf = 4.d0/3.d0
       ca = 3
       tk = - (1-x)*(1-y)*s/2
       uk = - (1-x)*(1+y)*s/2
       if(1-x.le.tiny)then
          q2c=-s-q2q
          t = q1q
          u = q2c
          born = 2*(u**2+4*xm2*s*(1-xm2*s/(t*u))+t**2)*tf/(t*u)
          p13 = -q1q/2
          p23 = -q2c/2
          p12 = s/2
C         --------------------------------------------------------------
C         Fattori iconali moltiplicati per 4*tk*uk
c         p1.k = -tk/2
c         p2.k = -uk/2
c         p3.k = w1/2
c         p4.k = w2/2
c
          p14 = p23
          p24 = p13
c
          e23 = 16*(1-y)*p23/w1h
          e24 = 16*(1-y)*p24/w2h
          e33 = 16*(1-y)*(1+y)* xm2/w1h**2
          e44 = 16*(1-y)*(1+y)* xm2/w2h**2
          e34 = 16*(1-y)*(1+y)* (s/2-xm2)/(w1h*w2h)
          sum = born*(  cf * (-e33-e44+2*e23+2*e24)
     #       + 2*(cf-ca/2) * ( e34-e23-e24)         )
          fpg = 1/(2*s)*sum
       elseif(1+y.le.tiny)then
          q2c = -s-uk-q2q
          sx = s*x
          t =   q1q
          u = - sx - t
          azidep = - (tf*4)*(t*u/(sx*xm2)-1)*
     #               (2*cth2**2-1)*xm2**2*sx**2/(t**2*u**2)
          born = 2*(u**2+4*xm2*sx*(1-xm2*sx/(t*u))+t**2)*tf/(t*u)
          sum = - born*(8*tk)*2*ca*(x/(1-x)+(1-x)/x+x*(1-x))
     #          + azidep*(8*tk)*4*ca*(1-x)/x
          fpg = 1/(2*sx)*sum
       else
          s2 = s+tk+uk
          q1c=-s-tk-q1q
          q2c=-s-uk-q2q
          w1 =q2q-q1q-tk
          w2 =q1q-q2q-uk
c
          p12 = s/2
          p13 = q1q/2
          p14 = q1c/2
          p15 = tk/2
          p23 = q2c/2
          p24 = q2q/2
          p25 = uk/2
          p34 = (s2-2*xm2)/2
          p35 = w1/2
          p45 = w2/2
          ans = fpg1(xm2,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45)
          ans = 4*tk*uk*ans/(2*s*4*8)
          fpg = ans
      endif
      return
      end
C
C       Author:
C
C               R. K. Ellis 
C               Fermilab, November 1986
C
C       Modified by S. Frixione, P. Nason and G. Ridolfi, 21-7-92
C
C ---
c +++
c       This function subroutine, FPG1, calculates the
c       invariant matrix element squared for the process:
c
c       Q(-P1) + QBAR(-P2) --> G(P3) + G(P4) + PHOTON(P5)
c
c       summed over initial and final spins and colours including
c       masses for the quark antiquark pair.
c       No averaging is performed for initial spins or colours.
c
c       where
c
c       1.      P1...P5 are the four momenta of the
c               partons defined as above and satisfying
c                       P1 + P2 + P3 + P4 + P5 = 0
c               Thus in any physical process the incoming partons
c               will have negative energy.
c
c       2.      QM2 is the mass squared of the Quark and AntiQuark lines
c               with
c                       P1**2 = P2**2 = QM**2
c
c       Example:
c
c       invariant matrix element squared for the physical process,
c
c               Q(K1) + QBAR(K2) --> G(P3) + G(P4) + PH(P5)
c
c       (where K1 and K2 have positive energies)
c       averaged over initial and summed over final spins is given
c       by,
c               G**4 * 1.0/4.0/XN**2 *
c                      FPG1(QM2,P45,P25,P15,P35,P24,P14,P34,P12,P23,P13)
c
c       where G is the QCD coupling constant, EQH is the quark charge,
c       XN is the number of colours and where
c
c               P1 = -K1, and P2 = -K2
c
c ---
      function fpg1(qm2,p45,p25,p15,p35,p24,p14,p34,p12,p23,p13)
      implicit real * 8 (a-z)
      parameter         (xn = 3.d0)
      parameter         (xv = 8.d0)

      qm4 = qm2**2
      qm6 = qm2**3
      s   = 2 * (p12+qm2)

      res=
     & +(p13*p23*(p13**2+p23**2)+p14*p24*(p14**2+p24**2)
     & +p15*p25*(p15**2+p25**2))/p13/p23/p14/p24/p15/p25
     & *(-s+2*xn**2*(p13*p24+p23*p14)/p34)

      res=res+qm6*2*xv*
     & (+1/p13**2/p24**2+1/p13**2/p25**2
     & +1/p14**2/p23**2+1/p14**2/p25**2
     & +1/p15**2/p23**2+1/p15**2/p24**2)
     & +4*xv*qm6
     & *( +1/p13**2/p24/p25+1/p14**2/p23/p25
     & +1/p23**2/p14/p15+1/p24**2/p13/p15
     & +1/p15/p25*(1/p13/p24+1/p14/p23))
     & -4*qm6*(1/p15**2/p23/p24+1/p25**2/p14/p13)
     & -4*qm6*(p13*p23+p14*p24+p15*p25+p15*p23
     & +p15*p24+p25*p13+p25*p14)/p13/p23/p14/p24/p15/p25

      res=res+4*xv*qm4*(+1/p15**2*(1/p24+1/p23)
     & +1/p25**2*(1/p14+1/p13)+1/p14**2*(1/p23+1/p25)
     & +1/p13**2*(1/p24+1/p25)+1/p23**2*(1/p14+1/p15)
     & +1/p24**2*(1/p13+1/p15))
     & + qm4
     & * ( +s*(p13**2+p23**2+p14**2+p24**2+p15**2+p25**2)
     & -4*(p13*p14*p25+p23*p24*p15)  
     & -4*(p13*p24+p14*p23)*(p15+p25)
     & +4*s*(p15*p25+p14*p24+p13*p23))/p13/p23/p14/p24/p15/p25
     & +8*qm4*xn**2/p34*(1/p15**2+1/p25**2)

      res=res+qm4*xn**2*(32/p34*p13*p14*p23*p24
     & -2/p34*p13*p24*(p23**2+p14**2+p15**2+p25**2)
     & -2/p34*p14*p23*(p13**2+p24**2+p15**2+p25**2)+2/p34
     & *(p13*p24**2+p14*p23**2-p24*p13**2-p23*p14**2)*(p25-p15)
     & -2/p34*p13*p23*p24*p25-2/p34*p13*p14*p15*p23
     & -2/p34*p13*p14*p15*p24-2/p34*p14*p23*p24*p25
     & +2/p34*p13*p23*(p24**2+p14**2)
     & +2/p34*p14*p24*(p13**2+p23**2)
     & -8/p34*(p13*p24+p14*p23)*p15*p25
     & +2*(p13*p24**2+p13**2*p24+p14*p23**2+p14**2*p23)
     & +6*(p13*p24*p25+p13*p15*p24+p14*p23*p25+p14*p15*p23)
     & +8*(p13*p23*p24+p13*p14*p23+p13*p14*p24+p14*p23*p24)
     &  )/p13/p23/p14/p24/p15/p25
      res=res+qm2*xn**2*(
     & -4/p34*(p14*p15*p23*p25*(p15+p25)+(p13+p14)*p15*p23*p24*p25
     &             +p13*p15*p24*p25*(p15+p25)+p13*p14*p15*(p23+p24)*p25)
     & +4/p34*(p13*p23*p24**2*p25+p13*p14**2*p15*p23
     & +p13**2*p14*p15*p24+p14*p23**2*p24*p25)
     & -2/p34*(p13*p23*p24*p25**2+p13*p14*p15*p23**2
     & +p13*p14*p15*p24**2+p13*p14*p15**2*p23
     & +p13*p14*p15**2*p24+p13**2*p23*p24*p25
     & +p14*p23*p24*p25**2+p14**2*p23*p24*p25)
     & +24/p34*(p15+p25)*p13*p14*p23*p24
     & )/p13/p23/p14/p24/p15/p25
      res=res+qm2*xn**2*(
     & -2*(p23/p24+p24/p23)/p15**2-2*(p13/p14+p14/p13)/p25**2
     & -2*(p15/p13+p13/p15)/p24**2-2*(p14/p15+p15/p14)/p23**2
     & -2*(p23/p25+p25/p23)/p14**2-2*(p24/p25+p25/p24)/p13**2)
      res=res+qm2*xn**2*(+6*p13*p14*p23*p24
     & +27.d0/2*p13*p14*p25*(p23+p24)+27.d0/2*p15*p23*p24*(p13+p14)
     & +16*(p13*p24+p14*p23)*p15*p25 
     & -(p13*p24+ p14*p23)*(p25**2+p15**2)
     & -p14**2*p15*p23-p13**2*p15*p24-p13*p24**2*p25-p14*p23**2*p25 
     & -p13*p15*p24**2-p13**2*p24*p25-p14*p15*p23**2-p14**2*p23*p25
     & +1.5d0*(p13*p23*p24**2+p14*p24*p23**2
     & +p13*p23*p14**2+p14*p24*p13**2)
     & -0.5d0*p13*p24*(p23**2+p14**2)-0.5d0*p14*p23*(p13**2+p24**2)
     & +0.5d0*p13*p14*(p23**2+p24**2)+0.5d0*(p13**2+p14**2)*p23*p24
     & +0.5d0*(p13+p14)*p23*p24*p25+0.5d0*p13*p14*p15*(p23+p24)
     & )/p13/p23/p14/p24/p15/p25
      res=res+qm2*(
     & +p15*p25*(p13*p14+p23*p24)+p13*p14*p25**2+p23*p24*p15**2
     & +p13*p14**2*p25+p15*p23*p24**2+p13**2*p14*p25+p15*p23**2*p24
     & +(p13*p24+p14*p23)*(p25**2+p15**2)
     & +p14*p15*p23**2+p13**2*p24*p25
     & +p13*p15*p24**2+p14**2*p23*p25+p13**2*p15*p24+p14*p23**2*p25
     & +p14**2*p15*p23+p13*p24**2*p25
     & -p13*p23*(p24**2+p14**2)-p14*p24*(p23**2+p13**2)
     & -p13*p24*(p23**2+p14**2)-p14*p23*(p24**2+p13**2)
     & -p13*p14*p25*(p23+p24)-p23*p24*p15*(p13+p14)
     & -p13*p14*(p23**2+p24**2)-p23*p24*(p13**2+p14**2)
     & -2*p15*p25*(p13*p24+p14*p23)-4*p15*p25*(p13*p23+p14*p24))
     & /p13/p23/p14/p24/p15/p25
      res=res-2*qm2*(p14*p15*(p24**2+p25**2)+p24*p25*(p14**2+p15**2)
     & +p13*p15*(p23**2+p25**2)+p23*p25*(p13**2+p15**2))
     & /p13/p23/p14/p24/p15/p25
      res=res
     & +6*qm2*(1/p24/p25+1/p14/p15)
     & +6*qm2*(1/p23/p24+1/p13/p14)
     & +6*qm2*(1/p23/p25+1/p13/p15)
     & -12*qm2*(1/p14/p24+1/p13/p23+1/p15/p25 )
     & -16*qm2*(1/p14/p23+1/p13/p24)
     & -16*qm2*(1/p14/p25+1/p15/p24)
     & -16*qm2*(1/p13/p25+1/p15/p23)
     & +2*qm2*(1/p13**2*(p25/p24+p24/p25)+1/p14**2*(p25/p23+p23/p25)
     & +1/p15**2*(p24/p23+p23/p24)+1/p23**2*(p15/p14+p14/p15)
     & +1/p24**2*(p15/p13+p13/p15)+1/p25**2*(p14/p13+p13/p14))
     & +4*xn**2*qm2*((p23**2+p24**2)/p15**2/p34**2
     & +(p13**2+p14**2)/p25**2/p34**2
     & -(p13-p14)*(p23-p24)/p34**2/p15/p25-1/p34**2)

      fpg1 = xv/xn*res

      return
      end
c 
c Process pq
c
C
C       Author:
C
C               R. K. Ellis 
C               Fermilab, November 1986
C
C       Modified by S. Frixione, P. Nason and G. Ridolfi, 21-7-92
C
C ---
      subroutine fpq1(fpqeh2,fpqel2,fpqeleh,
     #    qm2,p45,p15,p25,p35,p14,p24,p34,p12,p13,p23)
      implicit real * 8 (a-z)
      qm4 = qm2**2
      qm6 = qm2**3
      s   = 2 * (p12+qm2)
      fpqeh2=
     & (p13**2+p23**2+p14**2+p24**2+qm2*s)*p12/p15/p25/s/p34
     & +0.5d0*qm2*(p35**2+p45**2)/p34**2/p15/p25
     & +qm2*(p13**2+p14**2+p23**2+p24**2-s*p34)/s/p34**2*(1/p15+1/p25)
     & -0.5d0*qm2/p34**2*((p13**2+p14**2+qm2*p34)/p25**2
     &                 +(p23**2+p24**2+qm2*p34)/p15**2)
      fpqeh2 = -fpqeh2
      fpqel2=
     & (p13**2+p23**2+p14**2+p24**2+2*qm2*p34)/s/p35/p45
     & +2*qm2*(p35**2+p45**2)/s**2/p35/p45
      fpqel2 = -fpqel2
      fpqeleh=
     & (p13**2+p23**2+p14**2+p24**2+qm2*(p34+0.5d0*s))/s/p34
     & *(p13/p15/p35+p24/p25/p45-p14/p15/p45-p23/p25/p35)
     & +2*qm2*((p23-p24)/p15-(p13-p14)/p25)/s/p34
      end

c
       subroutine fpq(feh2,fel2,feleh,
     #     s,x,y,xm2,q1q,q2q,w1h,w2h,cth2)

       implicit real * 8(a-z)
       xm = sqrt(xm2)
       tiny = .1e-4
       vtf = 1/2.d0
       vdf = 3
       vda = 8
       cf = 4.d0/3.d0
       ca = 3
       tk = - (1-x)*(1-y)*s/2
       uk = - (1-x)*(1+y)*s/2
       if (x.eq.1)then
          feh2 = 0
          fel2 = 0
          feleh = 0
       elseif(1+y.le.tiny)then
          sx = s*x
          t =   q1q
          u = - sx - t
          azidep = - (vtf*4)*(t*u/(sx*xm2)-1)*
     #               (2*cth2**2-1)*xm2**2*sx**2/(t**2*u**2)
          born = 2*(u**2+4*xm2*sx*(1-xm2*sx/(t*u))+t**2)*vtf/(t*u)
          sum = - born*(8*tk)*cf*(1+(1-x)**2)/x
     #          + azidep*(8*tk)*4*cf*(1-x)/x
          feh2 = 1/(2*sx)*sum
          fel2 = 0
          feleh = 0
       elseif(1-y.le.tiny)then
          sx = s*x
          t = q2q
          feh2 = 0
          fel2 = -8*uk*vdf*(x**2+(1-x)**2)*bornqq(sx,t,xm2)
          feleh = 0
       else
          s2 = s+tk+uk
          q1c=-s-tk-q1q
          q2c=-s-uk-q2q
          w1 =q2q-q1q-tk
          w2 =q1q-q2q-uk
c
          p12 = s/2
          p13 = q1q/2
          p14 = q1c/2
          p15 = tk/2
          p23 = q2c/2
          p24 = q2q/2
          p25 = uk/2
          p34 = (s2-2*xm2)/2
          p35 = w1/2
          p45 = w2/2
          call fpq1(fpqeh2,fpqel2,fpqeleh,
     #    xm2,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45)
          norm = vda/(2.d0*s*vdf)
          feh2 = 4*tk*uk*norm*fpqeh2
          fel2 = 4*tk*uk*norm*fpqel2
          feleh = 4*tk*uk*norm*fpqeleh
      endif
      return
      end
c
c
c- Single inclusive cross sections for heavy quark photoproduction
c
c
      function hqh0pg(t1,ro)
      implicit double precision (a-z)
      t2=1-t1
      hqh0pg=( t2**2 + t1**2 + ro*(1-ro/(4*t1*t2)) )/(t1*t2)
      end

      function phhqh0qa(t1,ro)
      implicit double precision (a-z)
      t2=1-t1
      phhqh0qa=2*(2*t1**2+2*t2**2+ro)/9
      return
      end

      function hqhlpg(tx,t1,ro)
      implicit double precision (a-z)
      hqhlpg = - 2*hqbppg(tx,t1,ro)
      return
      end

      function hqhlpn(tx,t1,ro)
      implicit double precision (a-z)
      hqhlpn = - 2*hqbppn(tx,t1,ro)
      return
      end

      function hqhlpc(tx,t1,ro)
      implicit double precision (a-z)
      hqhlpc = - 2*hqbppc(tx,t1,ro)
      return
      end

      function hqbdpg(t1,ro)
      implicit double precision (a-z)
      hqbdpg = 6*log(t1)*hqh0pg(t1,ro)
      end

      function hqbppg(tx,t1,ro)
      implicit double precision (a-z)
      pgg(x) = 6*(x+(1-x)**2*(1/x+x))
      t2 = 1-t1-tx
      hqbppg = -(1-t2)/t1*hqh0pg(1-t2,ro*(1-t2)/t1)*pgg(t1/(1-t2))
      return
      end

      function hqbppn(tx,t1,ro)
      implicit double precision (a-z)
      parameter (cf=4.d0/3)
      t2 = 1-t1-tx
      x = t1/(1-t2)
      pgq = cf*(1+(1-x)**2)/x * (1-x)
      hqbppn = -1/x*hqh0pg(1-t2,ro/x )*pgq
      return
      end

      function hqbppc(tx,t1,ro)
      implicit double precision (a-z)
      t2 = 1-t1-tx
      x = t2/(1-t1)
      pqp = 3*(x**2+(1-x)**2)*(1-x)
      hqbppc = -1/x*phhqh0qa(1-t1,ro*(1-t1)/t2)*pqp
      return
      end

      function hqhppg(tx,t1,ro)
      implicit double precision (a-z)
      if(tx.eq.0)then
      t2 = 1-t1
      b = dsqrt(1-ro)
      vlsm = dlog(4/ro)
      srl12 = dlog(t1/t2)
      srlg1 = dlog((b+1)/(1-b))/b
      pp = -3*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+ro**2)*vlsm/(
     1   (t1-1)**2*t1**2)
      pp = pp-(ro-2)*srlg1*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+
     1   ro**2)/( 12.d0*(t1-1)**2*t1**2)
      pp = 3*srl12*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+ro**2)/(
     1    2.d0*(t1-1)**2*t1**2)+pp
      pp = 4*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+ro**2)/( 3.d0*
     1   (t1-1)**2*t1**2)+pp
      hqhppg = pp
      return 
      else 
      t2 = -tx-t1+1
      t11 = 1/(1-t1)
      t22 = 1/(tx+t1)
      dro = 1/(4*tx+ro)
      b = dsqrt((1-tx)**2-ro)
      dlam2 = 1/(tx**2-2*tx-ro+1)
      vlsm = dlog(4/ro)
      rlgxro = dlog((4*tx+ro)/ro)/tx
      rlg12 = dlog(t1/(1-t2))/tx
      rlg21 = dlog(t2/(1-t1))/tx
      rlg22 = dlog(t2/(1-t2))
      rr = t1*(t1-ro*(1-t2))
      d2435 = 1/rr
      rr = dsqrt(rr)
      rl3524 = 2*dlog((t1+rr)/(t1-rr))/rr
      rr = t2*(t2-ro*(1-t1))
      d1435 = 1/rr
      rr = dsqrt(rr)
      rl3514 = 2*dlog((t2+rr)/(t2-rr))/rr
      rr = tx**2+ro*(1-t1)*(1-t2)
      rr = dsqrt(rr)
      rl1424 = 2*dlog((tx+rr)/(rr-tx))/(rr*tx)
      rl35 = dlog((2*tx+ro-2*b-2)/(2*tx+ro+2*b-2))/b
      pp = t1*(-t2-t1+1)*t22**2+(-t2-t1+1)/t1+t1/(-t2-t1+1)
      pp = -3*pp*(-t2-t1+1)*(-4/t22-4*ro/(t1*t22**2)+ro**2/(t1**2*t22**2
     1   )+12/t22**2+4*ro/(t1*t22**3)-16/t22**3+8/t22**4)*t22**2*vlsm/(t
     2   1*t2**2)
      tmp0 = -8*((8*t1+36)*t2**4+(25*t1**2+93*t1-117)*t2**3+(18*t1**3+99
     1   *t1**2-246*t1+135)*t2**2+(52*t1**3-169*t1**2+181*t1-63)*t2-18*t
     2   1**3+45*t1**2-36*t1+9)+4*ro*((8*t1+34)*t2**3+(10*t1**2+84*t1-85
     3   )*t2**2+(56*t1**2-148*t1+68)*t2-38*t1**2+56*t1-17)+ro**2*(-(4*t
     4   1+52)*t2**2-(52*t1-93)*t2+46*t1-41)+2*ro**3*(t2-1)
      tmp0 = -rlgxro*(t2+t1-1)*t22**3*tmp0/( 48.d0*t1*t2)
      pp = tmp0+pp
      tmp0 = -32*((3*t1**2-20*t1-13)*t2**2+(-20*t1**2-6*t1+26)*t2-13*t1*
     1   *2+26*t1-13)*tx+8*ro*((t1**2+18*t1+47)*t2**2+(18*t1**2+76*t1-94
     2   )*t2+47*t1**2-94*t1+47)-12*ro**2*((t1+6)*t2+6*t1-6)+3*ro**3
      tmp0 = dro*(t2+t1-1)*(t2+2*t1-1)*t22**3*tmp0/( 12.d0*t1*t2*(4*tx+r
     1   o))
      pp = tmp0+pp
      pp = pp-4*(t2+t1-1)**2*(ro-2*t1*t2)*t22**3/( 3.d0*t1)
      tmp0 = -ro*(144*t2**5+(80*t1-432)*t2**4+(-180*t1**2-315*t1+432)*t2
     1   **3+(-148*t1**3-270*t1**2+658*t1-144)*t2**2+(-316*t1**3+1106*t1
     2   **2-711*t1)*t2+360*t1**3-648*t1**2+288*t1)*tx-4*t1*((44*t1+19)*
     3   t2**4+(78*t1**2+22*t1-48)*t2**3+(50*t1**3+68*t1**2-236*t1+39)*t
     4   2**2+(70*t1**3-302*t1**2+242*t1-10)*t2-72*t1**3+144*t1**2-72*t1
     5   )*tx
      tmp0 = ro**2*(56*t2**4+(-48*t1-112)*t2**3+(-48*t1**2-40*t1+128)*t2
     1   **2+(-122*t1**2+303*t1-144)*t2+144*t1**2-216*t1+72)*tx+tmp0+ro*
     2   *3*(5*t2**4+(-8*t1-10)*t2**3+(-4*t1**2-9*t1+23)*t2**2+(-17*t1**
     3   2+53*t1-36)*t2+18*t1**2-36*t1+18)
      tmp0 = -t22**2*tmp0/( 6.d0*t1**2*t2**2*(4*tx+ro))
      pp = tmp0+pp
      tmp0 = -4*t1*(53*t2**3+(106*t1-113)*t2**2+(71*t1**2-149*t1+76)*t2+
     1   34*t1**3-69*t1**2+58*t1-16)+2*ro*(44*t2**3+(86*t1-88)*t2**2+(48
     2   *t1**2-87*t1+44)*t2+6*t1**3-24*t1**2)-2*ro**2*(10*t2**2+(14*t1-
     3   10)*t2-5*t1**2-2*t1)+ro**3*(t2+t1)
      tmp0 = -rlgxro*(t2+t1-1)*t22*tmp0/( 24.d0*t1**2*t2)
      pp = tmp0+pp
      pp = 2*d2435*(ro-2*t1)*(-2*t1**2-2*ro*t1+ro**2)*(t2+t1-1)**2*t22/(
     1    3.d0*t1)+pp
      tmp0 = -2*ro*t1*((5*t1-8)*t2**5+(9*t1**2-31*t1+32)*t2**4+(6*t1**3-
     1   30*t1**2+63*t1-48)*t2**3+(2*t1**4-15*t1**3+37*t1**2-53*t1+32)*t
     2   2**2+(-4*t1**4+17*t1**3-20*t1**2+16*t1-8)*t2+4*t1**4-8*t1**3+4*
     3   t1**2)
      tmp0 = tmp0+t1**3*(3*t2**5+(16*t1-8)*t2**4+(16*t1**2-51*t1+6)*t2**
     1   3+(8*t1**3-56*t1**2+70*t1)*t2**2+(-20*t1**3+72*t1**2-51*t1-1)*t
     2   2+16*t1**3-32*t1**2+16*t1)+ro**2*((4*t1**2-4)*t2**4+(3*t1**3-8*
     3   t1**2-8*t1+16)*t2**3+(-3*t1**3+24*t1-24)*t2**2+(t1**4+8*t1**2-2
     4   4*t1+16)*t2-4*t1**2+8*t1-4)
      tmp0 = -3*t22*tmp0/( 2.d0*t1**4*(t2-1)*t2**2)
      pp = tmp0+pp
      pp = pp-12*(ro-2*t1)*(t2**2+2*t1*t2-2*t2+2*t1**2-2*t1+1)/(t1**3*t2
     1   **2*t22)
      pp = pp-12*(t2**2+2*t1*t2-2*t2+2*t1**2-2*t1+1)/(t1**2*t2**2*t22**2
     1   )
      tmp0 = ro*((8*t1**2-20*t1-52)*t2**2+(-8*t1**3-12*t1**2-32*t1+52)*t
     1   2-7*t1**3+13*t1**2-5*t1-1)-2*t1*t2*((36*t1-52)*t2**2+(20*t1**2-
     2   72*t1+52)*t2+t1**3-3*t1**2+3*t1-1)-2*ro**2*((4*t1**2-2*t1-10)*t
     3   2-2*t1**2-3*t1+5)+4*ro**3*(t1-1)
      tmp0 = rlgxro*t11**3*(t2+t1-1)*tmp0/( 12.d0*t1*t2)
      pp = tmp0+pp
      tmp0 = 16*((3*t1**2+4*t1-1)*t2**2+(4*t1**2-6*t1+2)*t2-t1**2+2*t1-1
     1   )*tx+8*ro*((t1**2+9*t1+2)*t2**2+(9*t1**2-5*t1-4)*t2+2*t1**2-4*t
     2   1+2)-6*ro**2*((2*t1+3)*t2+3*t1-3)+3*ro**3
      tmp0 = -2*dro*t11**3*(t2+t1-1)*(2*t2+t1-1)*tmp0/( 3.d0*t1*t2*(4*tx
     1   +ro))
      pp = tmp0+pp
      pp = pp-4*t11**3*(t2+t1-1)**2*(ro-2*t1*t2)/( 3.d0*t2)
      tmp0 = -2*t2*((104*t1+52)*t2**3+(132*t1**2-76*t1-104)*t2**2+(53*t1
     1   **3-77*t1**2-29*t1+53)*t2+t1**3-3*t1**2+3*t1-1)+ro*t2*((74*t1+6
     2   8)*t2**2+(90*t1**2-54*t1-76)*t2+32*t1**3-63*t1**2+22*t1+9)-4*ro
     3   **2*((6*t1+4)*t2**2+(6*t1**2-4*t1-3)*t2+2*t1**3-4*t1**2+2*t1)+2
     4   *ro**3*((t1+2)*t2+t1**2-t1)
      tmp0 = t11**2*(t2+t1-1)*tmp0/( 3.d0*t1*t2**2*(4*tx+ro))
      pp = tmp0+pp
      tmp0 = 2*t2*(88*t2**3+(80*t1-132)*t2**2+(34*t1**2-68*t1+58)*t2+8*t
     1   1**3-23*t1**2+22*t1-7)+ro*(12*t2**3+(60*t1-12)*t2**2+(22*t1**2-
     2   21*t1-9)*t2-8*t1**3+16*t1**2-8*t1)+ro**2*(-14*t2**2-(31*t1-7)*t
     3   2-8*t1**2+8*t1)+4*ro**3*(t2+t1)
      tmp0 = rlgxro*t11*(t2+t1-1)*tmp0/( 12.d0*t1*t2**2)
      pp = tmp0+pp
      pp = pp-d1435*t11*(ro-2*t2)*(t2+t1-1)**2*(16*t2**2+16*ro*t2+ro**2)
     1   /( 12.d0*t2)
      pp = 3*t11*(ro-2*t2)*(t2+t1-1)**2*(ro-2*t1*t2)/(t1**2*t2**2)+pp
      pp = 3*rlg22*(2*t2*(t2-t1-1)*(2*t2**2+(t1-2)*t2+t1**2-t1+1)+ro**2*
     1   (t2**2+(2*t1-1)*t2+t1**2-t1+1)-2*ro*t2*((t1-1)*t2+t1**2+1))/( 2
     2   .d0*t1*(t2-1)*t2**2)+pp
      tmp0 = ro**2*((2*t1**2-2)*t2**3+(2*t1**3-2*t1**2-4*t1+6)*t2**2+(-2
     1   *t1**4+2*t1**3-6*t1**2+8*t1-6)*t2-t1**5+3*t1**4-5*t1**3+6*t1**2
     2   -4*t1+2)-2*ro*t1*t2*((2*t1-4)*t2**3+(6*t1**2-12*t1+12)*t2**2+(5
     3   *t1**3-16*t1**2+18*t1-12)*t2+t1**4-6*t1**3+10*t1**2-8*t1+4)
      tmp0 = tmp0-2*t1**2*t2*(2*t2+t1-2)**2*(2*t2**2+(t1-2)*t2+t1**2-t1+
     1   1)
      tmp0 = -3*rlg21*(t2+t1-1)*tmp0/( 2.d0*t1**4*(t2-1)*t2**2)
      pp = tmp0+pp
      tmp0 = ro**2*((2*t1**2-2)*t2**3+(2*t1**3-2*t1**2-4*t1+6)*t2**2+(-2
     1   *t1**4+2*t1**3-6*t1**2+8*t1-6)*t2-t1**5+3*t1**4-5*t1**3+6*t1**2
     2   -4*t1+2)-2*ro*t1*t2*((2*t1-4)*t2**3+(6*t1**2-12*t1+12)*t2**2+(5
     3   *t1**3-16*t1**2+18*t1-12)*t2+t1**4-6*t1**3+10*t1**2-8*t1+4)
      tmp0 = tmp0-2*t1**2*t2*(2*t2+t1-2)**2*(2*t2**2+(t1-2)*t2+t1**2-t1+
     1   1)
      tmp0 = -3*rlg12*(t2+t1-1)*tmp0/( 2.d0*t1**4*(t2-1)*t2**2)
      pp = tmp0+pp
      tmp0 = -2*ro*t1*(2*t2**3+(2*t1-4)*t2**2+(-6*t1-1)*t2+3*t1+3)-8*t1*
     1   *2*(2*t2**2+(2*t1-2)*t2-t1)+ro**2*(t2-1)*(3*t2-t1-3)
      tmp0 = -d2435**2*rl3524*ro**3*(t2+t1-1)*tmp0/ 24.d0
      pp = tmp0+pp
      tmp0 = -d2435*rl3524*ro*(t2+t1-1)*(-4*t1**2*(t1+1)*tx+2*ro*t1*((5*
     1   t1-1)*t2+t1**2-3*t1+1)-2*ro**2*t1*((t1+1)*t2+t1)+ro**3*((t1-1)*
     2   t2+1))/( 6.d0*t1)
      pp = tmp0+pp
      tmp0 = 2*ro*((8*t1+8)*t2**4+(16*t1**2-63*t1-16)*t2**3+(8*t1**3-96*
     1   t1**2+93*t1+8)*t2**2+(-25*t1**3+61*t1**2-38*t1)*t2-t1**3+t1**2)
      tmp0 = tmp0+2*ro**2*(8*t2**4+(24*t1-8)*t2**3+(24*t1**2-24*t1)*t2**
     1   2+(8*t1**3-16*t1**2+8*t1)*t2+t1**3+t1**2+t1)+ro**3*(t2-1)*((8*t
     2   1-8)*t2**2+(8*t1**2-15*t1+9)*t2+t1**2)-4*t1*(t1**2+1)*(8*t2**2+
     3   (8*t1-8)*t2+t1)
      tmp0 = -rl3524*tmp0/( 24.d0*t1*(t2-1)*t2)
      pp = tmp0+pp
      tmp0 = 2*ro*t2*((2*t1**2-6*t1+3)*t2+2*t1**3-4*t1**2-t1+3)+8*t2**2*
     1   ((2*t1-1)*t2+2*t1**2-2*t1)+ro**2*(t1-1)*(t2-3*t1+3)
      tmp0 = d1435**2*rl3514*ro**3*(t2+t1-1)*tmp0/ 24.d0
      pp = tmp0+pp
      tmp0 = -32*t2**2*(t2+1)*tx+16*ro*t2*(t2**2+(5*t1-3)*t2-t1+1)-2*ro*
     1   *2*t2*((8*t1-1)*t2-t1)+ro**3*((8*t1-9)*t2+t1-1)
      tmp0 = -d1435*rl3514*ro*(t2+t1-1)*tmp0/( 48.d0*t2)
      pp = tmp0+pp
      pp = pp-rl3514*(2*ro*((8*t1-8)*t2**3+(16*t1**2-71*t1+17)*t2**2+(8*
     1   t1**3-64*t1**2+65*t1-9)*t2-t1**3+t1**2)+2*ro**2*(8*t2**3+(24*t1
     2   -1)*t2**2+(24*t1**2-18*t1-1)*t2+8*t1**3-18*t1**2+18*t1-9)+ro**3
     3   *((8*t1-8)*t2**2+(8*t1**2-15*t1+9)*t2+t1**2)-4*t2*(8*t2+8*t1-9)
     4   *(t2**2+1))/( 24.d0*t1*t2)
      pp = dlam2*rl35*(t2-t1)*(t2+t1-1)*(t2+t1)*(ro-2*t1*(t2+t1))/( 12.d
     1   0*t1*t2)+pp
      pp = pp-rl35*(2*(9*t2**3+(-t1-10)*t2**2+(6*t1-9*t1**2)*t2+t1**3-2*
     1   t1**2)+ro*(-2*t2**2+3*t2-2*t1**2+3*t1)+ro**2*(t2+t1))/( 12.d0*t
     2   1*t2)
      tmp0 = -4*(t2**2+(2*t1-2)*t2+t1**2-2*t1+2)*tx+2*ro*((t1+1)*t2**2+(
     1   t1**2-2*t1-1)*t2+t1**2-t1)+2*ro**2*(t2**2+(2*t1-3)*t2+t1**2-3*t
     2   1+3)+ro**3*(t1-1)*(t2-1)
      tmp0 = rl1424*(t2+t1-1)**2*(8*t2+1)*tmp0/( 24.d0*t1*(t2-1)*t2)
      pp = tmp0+pp
      tmp0 = 6*t1*t2*(t2+t1-2)*tx+ro*(t2+t1)*(t2**2-t2+t1**2-t1)
      tmp0 = 16*dro*(t2+t1-1)**2*tmp0/( 3.d0*t1**2*t2**2*(4*tx+ro))
      pp = tmp0+pp
      tmp0 = ro-2*((t1+1)*t2+t1**2-t1)
      tmp0 = dlam2*(t2-t1)*(t2+t1-1)*(t2+t1)*tmp0/( 3.d0*t1*t2*(4*tx+ro)
     1   )
      pp = tmp0+pp
      pp = d2435**2*ro**2*(-4*ro*t1*(t2**3+(2*t1-2)*t2**2+t1**2*t2+2*t1*
     1   *2-t1+1)+ro**2*(3*t2**2+(2*t1-6)*t2+3*t1**2-2*t1+3)+4*t1**2*(2*
     2   t2**2+(4*t1-2)*t2+2*t1**2-2*t1+1))/ 6.d0+pp
      pp = 2*d2435*(ro**2-2*ro-2)*(2*t1*(t2+t1)+ro*(t2-t1-1))/ 3.d0+pp
      pp = d1435**2*ro**2*(-4*ro*t2*((t1+2)*t2**2+(2*t1**2-1)*t2+t1**3-2
     1   *t1**2+1)+ro**2*(3*t2**2+(2*t1-2)*t2+3*t1**2-6*t1+3)+4*t2**2*(2
     2   *t2**2+(4*t1-2)*t2+2*t1**2-2*t1+1))/ 6.d0+pp
      pp = pp-2*d1435*(ro**2-2*ro-2)*(ro*(t2-t1+1)-2*t2*(t2+t1))/ 3.d0
      tmp0 = -4*t1**2*(t2-1)*(144*t2**4+(313*t1-288)*t2**3+(182*t1**2-34
     1   1*t1+72)*t2**2+(-118*t1**3+91*t1**2-117*t1+144)*t2-72*t1**4+144
     2   *t1**3-216*t1**2+144*t1-72)*tx
      tmp0 = tmp0+ro*t1*((20*t1-1296)*t2**5+(259*t1**2-3948*t1+5472)*t2*
     1   *4+(540*t1**3-4876*t1**2+12660*t1-8928)*t2**3+(592*t1**4-2754*t
     2   1**3+10287*t1**2-14492*t1+6912)*t2**2+(200*t1**5-704*t1**4+2990
     3   *t1**3-6966*t1**2+6696*t1-2448)*t2-56*t1**5+112*t1**4-776*t1**3
     4   +1296*t1**2-936*t1+288)
      tmp0 = tmp0-2*ro**2*((73*t1**2-162*t1-126)*t2**4+(139*t1**3-470*t1
     1   **2+144*t1+504)*t2**3+(76*t1**4-330*t1**3+379*t1**2+540*t1-756)
     2   *t2**2+(10*t1**5+40*t1**4+29*t1**3+360*t1**2-864*t1+504)*t2+44*
     3   t1**5-116*t1**4+162*t1**3-342*t1**2+342*t1-126)
      tmp0 = tmp0+ro**3*((23*t1**2-63)*t2**3+(19*t1**3-23*t1**2-126*t1+1
     1   89)*t2**2+(-4*t1**4+17*t1**3-81*t1**2+252*t1-189)*t2+22*t1**4-3
     2   6*t1**3+81*t1**2-126*t1+63)
      tmp0 = -tmp0/( 6.d0*t1**4*(t2-1)*t2**2*(4*tx+ro))
      pp = tmp0+pp
      hqhppg = pp
      return 
      endif 
      end 

      function hqhdpg(t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      data pi/3.141 592 653 589 793/
      t2 = 1-t1
      b = dsqrt(1-ro)
      lp = (b+1)/ 2.d0
      lm = (1-b)/ 2.d0
      at = t1
      aw = t2
      vltm = dlog(4*at/ro)
      vlpm = dlog(lp/lm)
      vlsm = dlog(4/ro)
      vlwm = dlog(4*aw/ro)
      vlbl = dlog(b/lm)
      vdw = ddilog((aw-ro/ 4.d0)/aw)-vlwm**2/ 2.d0
      vdt = ddilog((at-ro/ 4.d0)/at)-vltm**2/ 2.d0
      vdmp = ddilog(-lm/lp)
      vdmb = vlbl**2/ 2.d0+ddilog(-lm/b)
      auinv = 1/(ro/ 4.d0-aw)
      atinv = 1/(ro/ 4.d0-at)
      srlgpr = dlog((1-b)*t2/((b+1)*t1))*dlog((b+1)*t2/((1-b)*t1))
      srl22p = ddilog(1-(b+1)*t2/((1-b)*t1))
      srl22m = ddilog(1-(1-b)*t2/((b+1)*t1))
      srlg1 = dlog((b+1)/(1-b))/b
      srl2l = (ddilog(-4*b/(1-b)**2)-ddilog(4*b/(b+1)**2))/b
      sl3525 = srlgpr+b**2*srlg1**2+srl22p+srl22m
      dd = -(8*(t1-1)**2*(16*t1**2-14*t1+5)+ro**2*(8*t1**2-16*t1+15)+4*r
     1   o*(t1-1)**2*(4*t1-11))*vlwm**2/( 48.d0*(t1-1)**3*t1)
      dd = 3*(2*(t1-1)*t1*(2*t1**2-2*t1+1)+ro**2*(t1**2-t1+1)+2*ro*(t1-1
     1   )*t1)*vltm*vlwm/( 2.d0*(t1-1)**2*t1**2)+dd
      dd = dd-3*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+ro**2)*vlsm
     1   *vlwm/( 2.d0*(t1-1)**2*t1**2)
      dd = (-4*(t1-1)*(t1**2-2*t1+2)-2*ro**2*(t1**2-3*t1+3)-2*ro*(t1-1)*
     1   t1+ro**3)*vlpm*vlwm/( 12.d0*b*(t1-1)**2*t1)+dd
      dd = auinv**2*(32*(t1-1)**3*(15*t1**2-28*t1+7)+8*ro*(t1-1)**2*(7*t
     1   1**2+38*t1-15)+12*ro**2*(t1-1)*(7*t1+6)+21*ro**3)*vlwm/( 192.d0
     2   *(t1-1)**2*t1)+dd
      dd = dd-auinv*(4*(t1-1)**2*(41*t1**2-30*t1+7)+2*ro*(t1-1)*(32*t1**
     1   2+31*t1-36)+ro**2*(16*t1**2+16*t1-25)+8*ro**3)*vlwm/( 24.d0*(t1
     2   -1)**2*t1)
      dd = dd-(8*t1**2*(16*t1**2-18*t1+7)+ro**2*(8*t1**2+7)-4*ro*t1**2*(
     1   4*t1+7))*vltm**2/( 48.d0*(t1-1)*t1**3)
      dd = 3*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+ro**2)*vlsm*vl
     1   tm/( 2.d0*(t1-1)**2*t1**2)+dd
      dd = dd-(-2*ro**2*(t1**2+t1+1)+4*t1*(t1**2+1)-2*ro*(t1-1)*t1+ro**3
     1   )*vlpm*vltm/( 12.d0*b*(t1-1)*t1**2)
      dd = dd-atinv**2*(8*ro*t1**2*(7*t1**2-70*t1+102)-160*t1**3*(3*t1**
     1   2-4*t1+6)+12*ro**2*t1*(7*t1-19)+21*ro**3)*vltm/( 192.d0*(t1-1)*
     2   t1**2)
      dd = atinv*(4*t1**2*(41*t1**2-70*t1+54)-2*ro*t1*(32*t1**2-95*t1+45
     1   )+ro**2*(16*t1**2-48*t1+7)+8*ro**3)*vltm/( 24.d0*(t1-1)*t1**2)+
     2   dd
      dd = dd-3*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+ro**2)*vlsm
     1   **2/( 2.d0*(t1-1)**2*t1**2)
      dd = dd-(-4*(t1-1)*t1*(10*t1**2-10*t1+7)+2*ro*(t1-1)*t1*(8*t1**2-8
     1   *t1-5)+2*ro**2*(7*t1**2-7*t1-3)+3*ro**3)*vlpm*vlsm/( 24.d0*b*(t
     2   1-1)**2*t1**2)
      dd = 4*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+ro**2)*vlsm/( 
     1   3.d0*(t1-1)**2*t1**2)+dd
      dd = dd-b*(ro-2*(2*t1**2-2*t1-1))*vlpm**2/( 24.d0*(t1-1)*t1)
      dd = dd-(ro**2-2*(2*t1**2-2*t1+3))*vlpm**2/( 24.d0*(t1-1)*t1)
      dd = dd-b*vlpm/( 6.d0*(t1-1)*t1)
      dd = (-(t1-1)*t1*(28*t1**2-28*t1+15)+8*ro*(t1-1)*t1*(2*t1**2-2*t1-
     1   1)+4*ro**2*(2*t1**2-2*t1-1)+2*ro**3)*vlpm/( 6.d0*b*(t1-1)**2*t1
     2   **2)+dd
      dd = (-16*(t1-1)**2*(17*t1**2-16*t1+7)-4*ro*(t1-1)**2*(4*t1+7)+7*r
     1   o**2*(2*t1-3)*(2*t1-1))*vdw/( 48.d0*(t1-1)**3*t1)+dd
      dd = (-16*t1**2*(17*t1**2-18*t1+8)+4*ro*t1**2*(4*t1-11)+7*ro**2*(2
     1   *t1-1)*(2*t1+1))*vdt/( 48.d0*(t1-1)*t1**3)+dd
      dd = (2*ro**2*(3*t1**2-3*t1-1)-4*(t1-1)*t1*(2*t1**2-2*t1+3)-2*ro*(
     1   t1-1)*t1+ro**3)*vdmp/( 12.d0*b*(t1-1)**2*t1**2)+dd
      dd = (ro-2)*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+ro**2)*vd
     1   mb/( 12.d0*b*(t1-1)**2*t1**2)+dd
      dd = dd-(ro-2)*srl2l*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+
     1   ro**2)/( 48.d0*(t1-1)**2*t1**2)
      dd = dd-3*sl3525*(4*(t1-1)*t1*(2*t1**2-2*t1+1)+4*ro*(t1-1)*t1+ro**
     1   2)/( 4.d0*(t1-1)**2*t1**2)
      dd = dd-pi**2*(2*ro*(t1-1)*t1*(16*t1**2-16*t1-7)-4*(t1-1)*t1*(14*t
     1   1**2-14*t1+5)+2*ro**2*(5*t1**2-5*t1-3)+3*ro**3)/( 144.d0*b*(t1-
     2   1)**2*t1**2)
      dd = dd-auinv*(56*ro*(t1-1)**2*(pi**2*(t1**2+6*t1)-12*t1+48)+32*(t
     1   1-1)**3*(7*pi**2*t1**2+48*t1**2-132*t1+84)+12*ro**2*(t1-1)*(pi*
     2   *2*(7*t1+7)-16*t1+58)+21*pi**2*ro**3)/( 1152.d0*(t1-1)**3*t1)
      dd = dd-atinv*(-32*t1**3*(pi**2*(7*t1**2-14*t1+7)+48*t1**2+36*t1+2
     1   16)+8*ro*t1**2*(pi**2*(7*t1**2-56*t1+49)+84*t1+468)+12*ro**2*t1
     2   *(pi**2*(7*t1-14)-16*t1-42)+21*pi**2*ro**3)/( 1152.d0*(t1-1)*t1
     3   **3)
      dd = dd-(ro**2*(pi**2*(74*t1**4-148*t1**3+96*t1**2-22*t1-7)-96*t1*
     1   *2+96*t1)-6*ro*(t1-1)*t1*(pi**2*(20*t1**2-20*t1+7)+36*t1**2-36*
     2   t1+28)+8*(t1-1)**2*t1**2*(pi**2*(12*t1**2-12*t1+7)-96*t1**2+96*
     3   t1-174))/( 144.d0*(t1-1)**3*t1**3)
      hqhdpg = dd
      return
      end

      function hqhppn(tx,t1,ro)
      implicit double precision (a-z)
      if(tx.eq.0)then
      pp = 0
      hqhppn = pp
      return 
      else 
      t2 = -tx-t1+1
      t11 = 1/(1-t1)
      t22 = 1/(tx+t1)
      vlsm = dlog(4/ro)
      rlgxro = dlog((4*tx+ro)/ro)/tx
      rlg12 = dlog(t1/(1-t2))/tx
      rlg21 = dlog(t2/(1-t1))/tx
      pp = (-t2-t1+1)*(t1*t22+2*(1-t2)/t1-2)
      pp = -2*pp*(8*(1-t2)**4+4*ro*(1-t2)**3/t1-16*(1-t2)**3-4*ro*(1-t2)
     1   **2/t1+ro**2*(1-t2)**2/t1**2+12*(1-t2)**2-4*(1-t2))*t22**2*vlsm
     2   /( 3.d0*t1*t2**2)
      pp = pp-4*rlgxro*(t2+t1-1)**3*t22**2/( 3.d0*t1)
      tmp0 = -2*ro*t1*t2*((12*t1+16)*t2**4+(12*t1**2-4*t1-64)*t2**3+(2*t
     1   1**3-t1**2-60*t1+96)*t2**2+(5*t1**3-35*t1**2+84*t1-64)*t2-8*t1*
     2   *3+24*t1**2-32*t1+16)*tx
      tmp0 = 4*ro**2*(t2-1)*((2*t1**2+2*t1+2)*t2**3-6*t2**2+(t1**3+t1**2
     1   -6*t1+6)*t2+t1**3-3*t1**2+4*t1-2)*tx+tmp0
      tmp0 = -8*t1**2*t2**2*(4*t2**3+(4*t1-12)*t2**2+(t1**2-8*t1+12)*t2-
     1   2*t1**2+4*t1-4)*tx**2+tmp0+ro**3*(t2-1)*((t1**2+2)*t2**3+(-t1**
     2   2+4*t1-6)*t2**2+(3*t1**2-8*t1+6)*t2+t1**3-3*t1**2+4*t1-2)
      tmp0 = -2*t22**2*tmp0/( 3.d0*t1**4*t2**2*(4*tx+ro))
      pp = tmp0+pp
      pp = rlgxro*(t2+t1-1)**2*t22*(4*ro*tx+4*t1*(2*t2+t1-1)+ro**2)/( 3.
     1   d0*t1**2)+pp
      pp = pp-t22*(4*ro*t1*t2*((5*t1-8)*t2**2+(t1**2-13*t1+16)*t2+8*t1-8
     1   )*tx-4*t1**3*t2**2*(2*t2+t1-1)*tx+ro**2*((5*t1**2-8)*t2**3+(t1*
     2   *3-5*t1**2-16*t1+24)*t2**2+(-8*t1**2+32*t1-24)*t2+8*t1**2-16*t1
     3   +8))/( 3.d0*t1**4*t2**2)
      pp = pp-4*rlgxro*t11**2*(t2+t1-1)**2*(ro-2*t1*t2)/( 3.d0*t1)
      pp = 8*t11**2*(t2+t1-1)**2*(ro*(t1+2)-4*((2*t1+1)*t2+t1**2-1))/( 3
     1   .d0*t1*(4*tx+ro))+pp
      pp = 8*rlgxro*t11*(t2+t1-1)**2*(2*t2-1)/( 3.d0*t1)+pp
      pp = 2*t11*(ro-4*t2)*(t2+t1-1)**2*(ro-2*t1*t2)/( 3.d0*t1**2*t2**2)
     1   +pp
      tmp0 = -4*t1**2*t2*(4*t2**2+(2*t1-4)*t2+t1**2-2*t1+2)+ro**2*((t1**
     1   2-2)*t2+t1**2-2*t1+2)-4*ro*t1*t2*((t1-2)*t2+t1**2-2*t1+2)
      tmp0 = -2*rlg21*(t2+t1-1)**2*tmp0/( 3.d0*t1**4*t2**2)
      pp = tmp0+pp
      tmp0 = -4*t1**2*t2*(4*t2**2+(2*t1-4)*t2+t1**2-2*t1+2)+ro**2*((t1**
     1   2-2)*t2+t1**2-2*t1+2)-4*ro*t1*t2*((t1-2)*t2+t1**2-2*t1+2)
      tmp0 = -2*rlg12*(t2+t1-1)**2*tmp0/( 3.d0*t1**4*t2**2)
      pp = tmp0+pp
      pp = pp-(16*t1**2*t2*(8*t2**2+(8*t1-8)*t2-t1**2)*tx**2+4*ro**2*((6
     1   *t1**2-18*t1-18)*t2**2+(2*t1**3-18*t1**2-18*t1+36)*t2+t1**3-19*
     2   t1**2+36*t1-18)*tx-4*ro*t1*t2*((2*t1-72)*t2**2+(8*t1**2-146*t1+
     3   144)*t2+5*t1**3-76*t1**2+144*t1-72)*tx+ro**3*((3*t1**2-18)*t2**
     4   2+(t1**3-36*t1+36)*t2+t1**3-19*t1**2+36*t1-18))/( 3.d0*t1**4*t2
     5   **2*(4*tx+ro))
      hqhppn = pp
      return 
      endif 
      end
 
      function hqhppc(tx,t1,ro)
      implicit double precision (a-z)
      if(tx.eq.0)then
      pp = 0
      hqhppc = pp
      return 
      else 
      t2 = -tx-t1+1
      t11 = 1/(1-t1)
      b = dsqrt((1-tx)**2-ro)
      dlam2 = 1/(tx**2-2*tx-ro+1)
      vlsm = dlog(4/ro)
      rlgxro = dlog((4*tx+ro)/ro)/tx
      rlg21 = dlog(t2/(1-t1))/tx
      rl34 = dlog((2*tx**2+2*b*tx-2*tx-ro)/(2*tx**2-2*b*tx-2*tx-ro))/(b*
     1   tx)
      pp = 2*t11**2*t2**2-2*t11*t2+1
      pp = 4*pp*(ro*(1-t1)/t2+4*t1**2-4*t1+2)*(-t2-t1+1)*vlsm/( 3.d0*t2)
      pp = pp-8*rlgxro*t11**2*(t2+t1-1)**3/ 3.d0
      pp = 8*t11**2*(t2+t1-1)*(4*((t1+2)*t2+t1-1)*tx+ro*(3*t2+2*t1-2))/(
     1    3.d0*(4*tx+ro))+pp
      pp = pp-4*rlgxro*t11*(t2+t1-1)**2*(ro-2*(2*t2+t1-1))/ 3.d0
      pp = 4*t11*(t2+t1-1)*(ro-2*(2*t2+t1))/ 3.d0+pp
      tmp0 = ro*(2*t2+t1-1)-2*t2*(4*t2**2+(4*t1-2)*t2+2*t1**2-2*t1+1)
      tmp0 = -4*rlg21*(t2+t1-1)**2*tmp0/( 3.d0*t2**2)
      pp = tmp0+pp
      tmp0 = -2*ro*((t1+1)*t2+t1**2-t1+1)+4*t2*(t2+t1)+ro**2
      tmp0 = -dlam2*rl34*(t2+t1-1)**2*tmp0/ 3.d0
      pp = tmp0+pp
      pp = pp-4*rl34*(ro-2*t2)*(t2+t1-1)**2/ 3.d0
      tmp0 = 4*(2*t2**3+(4*t1-2)*t2**2+(2*t1**2-t1+1)*t2+t1**2-t1)-2*ro*
     1   (t2**2+(t1+2)*t2+2*t1-1)+ro**2
      tmp0 = -4*dlam2*(t2+t1-1)**2*tmp0/( 3.d0*(4*tx+ro))
      pp = tmp0+pp
      pp = pp-4*(t2+t1-1)*(ro*(2*t2+t1-1)-2*t2*(3*t2**2+(4*t1-2)*t2+t1**
     1   2-t1))/( 3.d0*t2**2)
      hqhppc = pp
      return 
      endif 
      end 

      function hqhppa(tx,t1,ro)
      implicit double precision (a-z)
      if(tx.eq.0)then
      pp = 0
      hqhppa = pp
      return 
      else 
      t2 = -tx-t1+1
      t11 = 1/(1-t1)
      t22 = 1/(tx+t1)
      b = dsqrt((1-tx)**2-ro)
      dlam2 = 1/(tx**2-2*tx-ro+1)
      rlgxro = dlog((4*tx+ro)/ro)/tx
      rlg12 = dlog(t1/(1-t2))/tx
      rlg21 = dlog(t2/(1-t1))/tx
      rlgro = dlog(4*(1-t1)*(1-t2)/ro)
      rr = t1*(t1*tx**2+ro*(1-t1))
      rr = dsqrt(rr)
      rl3414 = 2*dlog((t1*tx+rr)/(rr-t1*tx))/(rr*tx)
      rl34 = dlog((2*tx**2+2*b*tx-2*tx-ro)/(2*tx**2-2*b*tx-2*tx-ro))/(b*
     1   tx)
      pp = ro-2*(2*t2+t1-1)
      pp = -2*pp*rlgxro*(t2-1)*(t2+t1-1)**2*t22/( 3.d0*t1)
      pp = 2*(t2+t1-1)*t22*(-4*(4*t2**2+(3*t1-5)*t2+t1+1)*tx+ro*(-8*t2**
     1   2-(5*t1-11)*t2+t1-3)+ro**2*(t2-1))/( 3.d0*t1*(4*tx+ro))+pp
      tmp0 = ro*t1-2*(t1-1)*((t1+1)*t2-1)
      tmp0 = 4*rlgxro*t11**2*(t2+t1-1)**2*tmp0/( 3.d0*t1)
      pp = tmp0+pp
      tmp0 = ro-2*((t1+1)*t2+t1-1)
      tmp0 = -16*(t1+1)*t11**2*(t2+t1-1)**2*tmp0/( 3.d0*t1*(4*tx+ro))
      pp = tmp0+pp
      pp = 4*rlgxro*t11*(t2+t1-1)**2*(ro*(t1+1)-2*((3*t1+3)*t2+t1**2-2))
     1   /( 3.d0*t1)+pp
      pp = pp-4*t11*(t2+t1-1)*(ro-2*(2*t2+t1-1))/( 3.d0*t1)
      pp = pp-4*rlgro*(t2+t1-1)*(2*(4*t2**2+(4*t1-6)*t2+2*t1**2-4*t1+3)+
     1   ro)/( 3.d0*(t1-1)*t1)
      tmp0 = 2*t1*(8*t2**2+(6*t1-6)*t2+3*t1**2-4*t1+3)+ro*(t1-4)*(t1-1)
      tmp0 = 4*rlg21*(t2+t1-1)**2*tmp0/( 3.d0*(t1-1)*t1**2)
      pp = tmp0+pp
      tmp0 = 2*(4*t2**2+(2*t1-4)*t2+t1**2-2*t1+2)+ro*t1
      tmp0 = 4*rlg12*(t2+t1-1)**2*tmp0/( 3.d0*(t1-1)*t1)
      pp = tmp0+pp
      tmp0 = 2*t1*(4*t2**2+(6*t1-4)*t2+3*t1**2-4*t1+2)*tx+ro*((t1**2+4*t
     1   1-4)*t2+t1**3+t1**2-2*t1)
      tmp0 = -2*rl3414*(t2+t1-1)**2*tmp0/( 3.d0*(t1-1)*t1)
      pp = tmp0+pp
      pp = pp-2*dlam2*rl34*(t2+t1-1)**2*(t2+t1)*(ro-2*t2*(t2+t1))/( 3.d0
     1   *t1)
      tmp0 = ro*(t2+t1+2)-2*(2*t2**2+3*t1*t2+t1**2-t1)
      tmp0 = 2*rl34*(t2+t1-1)**2*tmp0/( 3.d0*t1)
      pp = tmp0+pp
      tmp0 = ro-2*((t1+1)*t2+t1**2-t1)
      tmp0 = -8*dlam2*(t2+t1-1)**2*(t2+t1)*tmp0/( 3.d0*t1*(4*tx+ro))
      pp = tmp0+pp
      pp = 2*(t2+t1-1)*(4*(8*t2-t1**2+4*t1-3)*tx+ro*(-(8*t1-24)*t2-5*t1*
     1   *2+24*t1-19)+ro**2*(t1-3))/( 3.d0*(t1-1)*t1*(4*tx+ro))+pp
      hqhppa = pp
      return 
      endif 
      end 
c
c Correction terms for change of subtraction scheme in single inclusive.
c Assume that the change of structure functions is given by:
c xkdij(), xkpij(x), xklij(x).
c
      function cthdpg(t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      data one/1.d0/
      xkd1 = xkdgg(nl)
      xkp1 = xkpgg(one,nl)
      xkl1 = xklgg(one,nl)
      xlg1 = log(t1)
      tmp  = -(xkd1-xlg1*xkp1+xlg1**2*xkl1/2) * hqh0pg(t1,ro)
      cthdpg = tmp
      end

      function cthppg(tx,t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      t2 = 1 - t1 - tx
      r2 = t2/(1-t1)
      r1 = t1/(1-t2)
      xkp1 = xkpgg(r1,nl)
      xkl1 = xklgg(r1,nl)
      tmp = -(xkp1-log(1-t2)*xkl1)*hqh0pg(t2,ro/r1)/r1
      cthppg = tmp
      end

      function cthlpg(tx,t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      t2 = 1 - t1 - tx
      r2 = t2/(1-t1)
      r1 = t1/(1-t2)
      xkl1 = xklgg(r1,nl)
      tmp = -xkl1*hqh0pg(t2,ro/r1)/r1
      cthlpg = tmp
      end

      function cthdpn(t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      data one/1.d0/
      xkd1 = xkdgq(nl)
      xkp1 = xkpgq(one,nl)
      xkl1 = xklgq(one,nl)
      xlg1 = log(t1)
      tmp  = -(xkd1-xlg1*xkp1+xlg1**2*xkl1/2) * hqh0pg(t1,ro)
      cthdpn = tmp
      return
      end

      function cthppn(tx,t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      t2 = 1 - t1 - tx
      r2 = t2/(1-t1)
      r1 = t1/(1-t2)
      xkp1 = xkpgq(r1,nl)
      xkl1 = xklgq(r1,nl)
      tmp = -(xkp1-log(1-t2)*xkl1)*hqh0pg(t2,ro/r1)/r1
      cthppn = tmp
      return
      end

      function cthlpn(tx,t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      t2 = 1 - t1 - tx
      r1 = t1/(1-t2)
      xkl1 = xklgq(r1,nl)
      tmp = -xkl1*hqh0pg(t2,ro/r1)/r1
      cthlpn = tmp
      return
      end

      function cthdpc(t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      data one/1.d0/
      vnc = 3
      vtf = 1/2.d0
      xkd2 = xkdqg(nl)
      xkp2 = xkpqg(one,nl)
      xkl2 = xklqg(one,nl)
      t2   = 1 - t1
      xlg2 = log(t2)
c      tmp  = -(xkd2-xlg2*xkp2+xlg2**2*xkl2/2) * phhqh0qa(t1,ro)
      tmp  = -vnc/vtf*(xkd2-xlg2*xkp2+xlg2**2*xkl2/2) * phhqh0qa(t1,ro)
      cthdpc = tmp
      end

      function cthppc(tx,t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      vnc = 3
      vtf = 1/2.d0
      t2 = 1 - t1 - tx
      r2 = t2/(1-t1)
      xkp2 = xkpqg(r2,nl)
      xkl2 = xklqg(r2,nl)
c      tmp = -(xkp2-log(1-t1)*xkl2)*phhqh0qa(t1,ro/r2)/r2
      tmp = -vnc/vtf*(xkp2-log(1-t1)*xkl2)*phhqh0qa(t1,ro/r2)/r2
      cthppc = tmp
      return
      end

      function cthlpc(tx,t1,ro,nl)
      implicit double precision (a-z)
      integer nl
      vnc = 3
      vtf = 1/2.d0
      t2 = 1 - t1 - tx
      r2 = t2/(1-t1)
      xkl2 = xklqg(r2,nl)
c      tmp = -xkl2*phhqh0qa(t1,ro/r2)/r2
      tmp = -vnc/vtf*xkl2*phhqh0qa(t1,ro/r2)/r2
      cthlpc = tmp
      end
c
c
c Funzioni di cambio di schema MS ----> DI
c
c
      function xkpqq(x,nl)
      implicit double precision (a-z)
      parameter (fot=4/3.d0)
      integer nl
      xkpqq = fot*(-1.5d0-(1+x**2)*log(x)+(1-x)*(3+2*x))
      return
      end

      function xkdqq(nl)
      implicit double precision (a-z)
      parameter (fot=4/3.d0)
      integer nl
      data pi/3.141 592 653 589 793/
      xkdqq = -fot*(4.5d0 + pi**2/3)
      return
      end

      function xklqq(x,nl)
      implicit double precision (a-z)
      parameter (fot=4/3.d0)
      integer nl
      xklqq = fot*(1+x**2)
      return
      end

      function xkpqg(x,nl)
      implicit double precision (a-z)
      integer nl
      xkpqg = (1-x)*(-(x**2+(1-x)**2)*log(x)+8*x*(1-x)-1)/2
      return
      end

      function xkdqg(nl)
      implicit double precision (a-z)
      integer nl
      xkdqg = 0
      return
      end

      function xklqg(x,nl)
      implicit double precision (a-z)
      integer nl
      xklqg = (1-x)*(x**2+(1-x)**2)/2
      return
      end

      function xkdgg(nl)
      implicit double precision (a-z)
      integer nl
      xkdgg = - 2 * nl * xkdqg(nl)
      return
      end

      function xkpgg(x,nl)
      implicit double precision (a-z)
      integer nl
      xkpgg = - 2 * nl * xkpqg(x,nl)
      return
      end

      function xklgg(x,nl)
      implicit double precision (a-z)
      integer nl
      xklgg = - 2 * nl * xklqg(x,nl)
      return
      end

      function xkdgq(nl)
      implicit double precision (a-z)
      integer nl
      xkdgq = - xkdqq(nl)
      return
      end

      function xkpgq(x,nl)
      implicit double precision (a-z)
      integer nl
      xkpgq = - xkpqq(x,nl)
      return
      end

      function xklgq(x,nl)
      implicit double precision (a-z)
      integer nl
      xklgq = - xklqq(x,nl)
      return
      end

c
c Fit to heavy quark partonic total production cross section.
c The cross section is:
c
c  A)  p+g --> Q+Qb
c   sipg(ro) = 
c   as^2/M^2*(eQ^2*c0pg(ro)+g^2*eQ^2*(f1pg(ro)+fbpg(ro)*log(mu^2/M^2))
c
c  B)  p+q --> Q+Qb
c   sipq(ro) = as^2/M^2*(eQ^2*g^2*(c1pq(ro)+cbpq(ro)*log(mu^2/M^2)) +
c                        eq^2*g^2*(d1pq(ro)+dbpq(ro)*log(mu^2/M^2)) );
c
c All above is in ms_bar. for q+qb-->Q+Qb we also give the
c cross section to be convoluted with the parton model densities
c from deep inelastic scattering:
c
c  D)  p+q --> Q+Qb
c   sipq(ro) = as^2/M^2*(eQ^2*g^2*(cppq(ro)+cbpq(ro)*log(mu^2/M^2)) +
c                        eq^2*g^2*(dppq(ro)+dbpq(ro)*log(mu^2/M^2)) );
c
c with g^2 = 4*pi*as, ro = 4*M^2/s, s is the partonic CM energy,
c M is the heavy quark mass.
c The symbols are: g for gluon, p for photon, q for light quark,
c qb for light antiquark, Q for heavy quark, Qb for heavy antiquark.
c The ratio of the fits to the numerically integrated results is
c well within 1% of accuracy.
c
      function f0qq(ro)
      implicit double precision(a-h,o-z)
      data xk1/1.163552835d-1/
      b = dsqrt(1-ro)
      f0qq = b*xk1*ro*(ro+2)
      return
      end

      function c0pg(ro)
      implicit double precision(a-h,o-z)
      data pi/3.141592654d0/
      b = dsqrt(1-ro)
      c0pg = pi/4*ro*((3-b**4)*log((b+1)/(1-b))+b*(-4+2*b**2))
      return
      end

      function cbpg(ro)
      implicit double precision(a-h,o-z)
      data pi/3.141592654d0/,dlog2/.6931471806d0/
      b = dsqrt(1-ro)
      xlpm = log((b+1)/(1-b))
      h1   = log((1+b)/2)**2-log((1-b)/2)**2+
     #       2*(ddilog((1+b)/2)-ddilog((1-b)/2))
      h2   = ddilog(2*b/(1+b))-ddilog(-2*b/(1-b))
      c0pg = pi/4*ro*((3-b**4)*xlpm+b*(-4+2*b**2))
      cbpg = pi/4*( -b/6*(112-250*ro+285*ro**2) 
     #  - ro/4*(72*(1-ro)+47*ro**2)*xlpm + 3*ro*(ro*(ro+2)-4)*h1
     #  - 6*ro*(ro*(ro-2)-2)*h2   ) + 6*c0pg*log(ro/4/b**2)
      cbpg = cbpg/8/pi**2
      return
      end

      function cbpq(ro)
      implicit double precision(a-h,o-z)
      data pi/3.141592654d0/
      b = dsqrt(1-ro)
      xlpm = log((b+1)/(1-b))
      h1   = log((1+b)/2)**2-log((1-b)/2)**2+
     #       2*(ddilog((1+b)/2)-ddilog((1-b)/2))
      cbpq = pi/3*(ro*(ro-2)*h1-4*b/9*(14-29*ro+3*ro**2)
     #       -2*ro/3*(3+ro**2)*xlpm)
      cbpq = cbpq/8/pi**2
      return
      end

      function dbpq(ro)
      implicit double precision(a-h,o-z)
      data pi/3.141 592 653 589 793d0/
      b = dsqrt(1-ro)
      xlpm = log((b+1)/(1-b))
      dbpq = pi/27*ro*((9*ro**2/2-6)*xlpm+b*(16-13*ro))
      dbpq = dbpq/8/pi**2
      return
      end

      function c1pg(ro)
      implicit double precision (a-h,o-z)
      data pi/3.141 592 653 589 793d0/
      data
     # a0 / 0.5493001968106041  d0/
     #,a1 / -0.04893965650111131d0/
     #,a2 / -1.387054355623665  d0/
     #,a3 / -3.203219176270447  d0/
     #,a4 / -1.839152562092848  d0/
     #,a5 / 1.132469363107884   d0/
     #,a6 / -0.6690506398471837 d0/
     #,a7 / -2.304832667341383  d0/
      b=sqrt(1-ro)
      xlgb2=log(8*b**2)
      xlgro=log(ro)
      tmp=ro*(6*b*xlgb2**2-30*b*xlgb2-pi**2/6.0)/pi/16.0+41.0*b/(36.0*pi
     1   )
      c1pg = b*ro*(a5*ro*xlgro**2+a6*xlgro**2+a4*ro*xlgro+a7*xlgro+a1*b*
     1   *2*xlgb2**2+a2*b**2*xlgb2+a3*b**2+a0)+tmp
      return 
      end 

      function c1pq(ro)
      implicit double precision (a-h,o-z)
      data pi/3.141 592 653 589 793d0/
      data
     # a0 / 0.06936784403775915 d0/
     #,a1 / -2.978407164169949  d0/
     #,a2 / 0.09065877200667612 d0/
     #,a3 / -0.2954977301247192 d0/
     #,a4 / -1.500457878130202  d0/
     #,a5 / 0.5010419401425481  d0/
     #,a6 / -1.309177419861434  d0/
     #,a7 / -0.3495119637625088 d0/
      b=sqrt(1-ro)
      xlgb=log(b)
      xlgro=log(ro)
      c1pq = ro*(a5*b*ro*xlgro**2+a7*b*xlgro**2+a4*b*ro*xlgro+a6*b*xlgro
     1   +a2*b**5*xlgb+a0*b**3*xlgb+a3*b**5+a1*b**3)+41.0*b**3/(81.0*pi)
      return 
      end 

      function d1pq(ro)
      implicit double precision (a-h,o-z)
      data pi/3.141 592 653 589 793d0/
      data
     #  a0 / 0.03530077461522012  d0/
     # ,a1 / -0.05381769977456701 d0/
     # ,a2 / -0.01615128875446254 d0/
     # ,a3 / 0.03291557700487208  d0/
     # ,a4 / -0.03949469158293252 d0/
     # ,a5 / 0.0004413852878932756d0/
     # ,a6 / -0.003753695895391341d0/
     # ,a7 / 0.008976644490046603 d0/
      b=sqrt(1-ro)
      xlgb=log(b)
      xlgro=log(ro)
      d1pq = ro*(a5*b*ro*xlgro**2+a7*b*xlgro**2+a4*b*ro*xlgro+a6*b*xlgro
     1   +a2*b**5*xlgb+a0*b**3*xlgb+a3*b**5+a1*b**3)
      return 
      end 
c
c-------------------------------------------------------------------------
c 
       function cpgg(ro)
       implicit real*8 (a-h,o-z)
       common/nl/nl
       data pi/3.141 592 653 589 793/
       data
     #  a1 / -.385522052D-02/,
     #  a2 / 0.113890896D+00/,
     #  a3 / -.675542043D-01/,
     #  a4 / -.217390046D+00/,
     #  a5 / -.820328537D-01/,
     #  a6 / -.111661820D-01/,
     #  a7 / -.351533910D-02/,
     #  a8 / 0.924159776D-01/,
     #  a9  / -.816968116D-02/,
     #  a10 / 0.751969048D+00/,
     #  a11 / -.617323899D-03/,
     #  a12 / -.585222400D+00/,
     #  a13 / 0.123517429D+00/,
     #  a14 / 0.712403132D-01/

       vtf =1/2.d0
       b=sqrt(1-ro)
       xlgb=log(b)
       xlgro=log(ro)
       cpgg = a1*ro**3+a2*ro**3*b**2+a3*ro**3*b**4+a4*b**2+
     #         a5*b**2*xlgro+a6*b**2*xlgro**2+a7*b**2*xlgro**3+
     #         a8*ro*b**2*xlgro**2+a9*ro*b**2*xlgro**3+
     #         a10*ro*b**2*xlgb+a11*xlgb+a12*b**2*xlgb+
     #         a13*b**2*xlgb**2+a14*b**2*xlgb**3
       cpgg = 2*nl*vtf*cpgg*ro*b**2
       return
       end


       function cpgq(ro)
       implicit real*8 (a-h,o-z)
       external cpgqdel,f0pg
       data pi/3.141 592 653 589 793/
       data
     #  a1 / 0.301758208D+00/,
     #  a2 / 0.518109439D+00/,
     #  a3 / 0.320434185D+00/,
     #  a4 / 0.530286582D+00/,
     #  a5 / 0.115467870D-01/,
     #  a6 / 0.409453571D-01/,
     #  a7 / -.384938863D-02/,
     #  a8 / -.135069167D+00/,
     #  a15 / -.312940377D-02/,
     #  a16 / -.950855211D-01/,
     #  a17 / -.211402755D-01/

       b=sqrt(1-ro)
       xlgb=log(b)
       xlgro=log(ro)
       cpgq = a1*ro**3+a2*ro**3*b**2+a3*ro**3*b**4+a4*b**2+
     #         a5*b**2*xlgro+a6*b**2*xlgro**2+a7*b**2*xlgro**3+
     #         a8*ro*b**2*xlgro**2+a15*xlgb**2/b+a16*xlgb/b+a17/b
       cpgq = cpgq*ro*b**2
       cpgq = cpgq+f0pg(ro)*cpgqdel(ro)
       return
       end


       function cpgqdel(ro)
       implicit real*8 (a-h,o-z)
       data pi/3.141 592 653 589 793/

       pi2 = pi*pi
       vcf = 4/3.d0
       be = sqrt(1-ro)
       cpgqdel = vcf*(9/2.d0+pi2/3.d0+3*log(be)-4*log(be)**2)
       cpgqdel = -1/(8*pi2)*cpgqdel
       return
       end       


       function f0pg(ro)
       implicit real*8 (a-h,o-z)
       data pi/3.141 592 653 589 793/

       pi2 = pi*pi
       be = sqrt(1-ro)
       f0pg = pi/4.d0*
     #   (ro*(3-be**4)*log((1+be)/(1-be))-4*ro*be+2*ro*be**3)
       return
       end       


       function cqqp(ro)
       implicit real*8 (a-h,o-z)
       data pi/3.141 592 653 589 793/

       pi2 = pi*pi
       vnc = 3.d0
       vtf = 1/2.d0
       b = sqrt(1-ro)
       cqqp = (1-3*ro**2/4.d0)*( log(ro/(4*b**2))*log((1+b)/(1-b))
     #                           +ddilog(2*b/(1+b))-ddilog(-2*b/(1-b)) )
     #        +b*(13*ro/6.d0-8/3.d0)*log(ro/(4*b**2))
     #        -3*ro**2*log((1+b)/(1-b))+b*(115*ro/9.d0-70/9.d0)
       cqqp = cqqp*pi*ro/(27*8*pi2)
       cqqp = vnc/vtf*cqqp
       return
       end

