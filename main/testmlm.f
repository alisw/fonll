      implicit none
      integer nptype,ngroup,nset,ih,ndns
      real * 4 fx(-6:6)
      real * 8 dxpdf(-6:6)
c used by pdfset
      real * 8 val(20)
      character * 20 parm(20)
c
      real * 8 x,q
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
c      call pdfset(parm,val)
c      call pftopdg(x,q,dxpdf)
      call mlmpdf(ndns,ih,sngl(q**2),sngl(x),fx,6)
      do j=-6,6
         write(*,*) fx(j)
      enddo
      goto 1
      end
