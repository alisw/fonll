      implicit none
      integer nptype,ngroup,nset,ih,ndns
      real * 4 fx(-6:6)
      real * 8 dxpdf(-6:6)
c used by pdfset
      real * 8 val(20)
      character * 20 parm(20)
c
      real * 8 x,q,tmp
      integer j
 1    continue
      write(*,*) ' enter nptype,ngroup,nset (pdl) and ndns (mlm),x,q'
      read(*,*) nptype,ngroup,nset,ndns,x,q
      if(nptype.eq.1) then
         ih=1
      elseif(nptype.eq.2) then
         ih=3
      elseif(nptype.eq.3) then
         ih=4
      endif
      parm(1) = 'NPTYPE'
      val (1) =  nptype
      parm(2) = 'NGROUP'
      val (2) =  ngroup
      parm(3) = 'NSET'
      val (3) =  nset
      call pdfset(parm,val)
      call pftopdg(x,q,dxpdf)
      tmp     = dxpdf(-1)
      dxpdf(-1) = dxpdf(-2)
      dxpdf(-2) = tmp
      dxpdf(1)=dxpdf(-1)
      dxpdf(2)=dxpdf(-2)
      if(nptype.eq.3.and.ngroup.eq.6.and.nset.eq.3) then
         do j=-6,6
            if(j.ne.0)dxpdf(j)=dxpdf(j)/2.d0
         enddo
      endif
c      call mlmpdf(ndns,ih,sngl(q**2),sngl(x),fx,6)
      do j=-6,6
         write(*,*) dxpdf(j)/x
      enddo
      goto 1
      end
