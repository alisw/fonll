C isigm       for choice of cross-section
c     isigm=1 ==> dsigma/dy/dpt2
c     isigm=2 ==> e*dsigma/d3p
c     isigm=3 ==> dsigma/dy/dpt
c iloopas=1,2 for 1 or 2 loop in alphas (not implemented!)
c iloopfr=1,2 for leading or nl fragmentation function
c ihior=1,2   for born and NLO partonic cross sections
      integer ias2term,ias3term,iwhichsfh
      common/hdrsiord/ias2term,ias3term,iwhichsfh
      integer iloopas,iloopfr
      common/fonfrloop/iloopas,iloopfr
c this not to do change of scale in structure functions, for testing
      integer istrsc
      common/test0/istrsc
c this is to use alternative expansion point in evmat
      integer ialtev
      common/altevc/ialtev
