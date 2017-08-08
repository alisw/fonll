      subroutine pdftype
      call prntsf
      end

      subroutine setlam5(xlam,sche)
      implicit none
      real * 8 xlam
      integer mode
      common/modec/mode
      integer ih,iret
      character * 2 sche
      ih=1
      call pdfpar(mode,ih,xlam,sche,iret)
      if(iret.ne.0) then
         write(*,*) ' error: set ',mode,'does not exist in mlm'
         stop
      endif
      write(*,*) ' lambda(pdf)=',xlam
      write(*,*) ' scheme(pdf)=',sche
      end

      subroutine selectpdf(nptype, ngroup, nset)
      implicit none
      integer nptype, ngroup, nset
c select the structure function set
      integer mode
      common/modec/mode
      mode=nset
      end

      subroutine setlam52(mode,ih,xlam,sche)
      implicit none
      real * 8 xlam
      integer mode,ih,iret
      character * 2 sche
c
      call pdfpar(mode,ih,xlam,sche,iret)
      if(iret.ne.0) then
         write(*,*) ' error: set ',mode,'does not exist in mlm'
         stop
      endif
      write(*,*) ' lambda(pdf)=',xlam
      write(*,*) ' scheme(pdf)=',sche
      end

      subroutine
     # fxab(beam,proc,nf,xmu2,x,xfg,xfu,xub,xfd,xdb,xfs,xfc,xfb)
c returns the densities times x
      implicit none
c      implicit double precision (a-h,o-z)
      real * 8 xmu2,x,xfg,xfu,xub,xfd,xdb,xfs,xfc,xfb
      integer nf
      character * 3 beam,proc*2
c Global
      real * 8 xlramu,ffact
      integer irunsc,istrsc
      common/rensca/xlramu,ffact,irunsc,istrsc
      integer mode
      common/modec/mode
c Local
      integer ih
      real * 8 tmp
      real * 8 xq2,xx
      real * 8 fx(-5:5)
c
      xfg = 0
      xfu = 0
      xfd = 0
      xub = 0
      xdb = 0
      xfs = 0
      xfc = 0
      xfb = 0
      call sybnumb(beam,ih)
      xq2=xmu2
      xx=x
      if(istrsc.eq.1) then
         xq2=xq2/ffact**2
      endif
      call fonllmlmpdf(mode,ih,xq2,xx,fx,5)
      if(proc.ne.'qq') xfg = fx(0)*x
      if(proc.ne.'gg') then
         xfu = fx(1)*x
         xfd = fx(2)*x
         xub = fx(-1)*x
         xdb = fx(-2)*x
         xfs = fx(3)*x
         if(nf.ge.4) xfc = fx(4)*x
         if(nf.ge.5) xfb = fx(5)*x
      endif
      end

      subroutine
     # fxab2(beam,mode,proc,nf,xmu2,x,xfg,xfu,xub,xfd,xdb,xfs,xfc,xfb)
c returns the densities times x
      implicit none
c      implicit double precision (a-h,o-z)
      real * 8 xmu2,x,xfg,xfu,xub,xfd,xdb,xfs,xfc,xfb
      integer nf
      character * 3 beam,proc*2
c Global
      real * 8 xlramu,ffact
      integer irunsc,istrsc
      common/rensca/xlramu,ffact,irunsc,istrsc
      integer mode
c Local
      integer ih
      real * 8 tmp
      character sign*1,part*2
      real * 8 xq2,xx
      real * 8 fx(-5:5)
c
      sign = beam(3:3)
      part = beam(1:2)
      xfg = 0
      xfu = 0
      xfd = 0
      xub = 0
      xdb = 0
      xfs = 0
      xfc = 0
      xfb = 0
c      if(part.eq.'pr'.or.part.eq.'nu') then
c         ih=1
      if(part.eq.'pr') then
         ih=1
      elseif(part.eq.'nu') then
         ih=0
      elseif(part.eq.'pi') then
         ih=3
      elseif(part.eq.'ph') then
         ih=4
      elseif(part.eq.'el') then
         ih=5
      endif
      xq2=xmu2
      xx=x
      if(istrsc.eq.1) then
         xq2=xq2/ffact**2
      endif
      call fonllmlmpdf(mode,ih,xq2,xx,fx,5)
      if(proc.ne.'qq') xfg = fx(0)*x
      if(proc.ne.'gg') then
         xfu = fx(1)*x
         xfd = fx(2)*x
         xub = fx(-1)*x
         xdb = fx(-2)*x
         xfs = fx(3)*x
         if(nf.ge.4) xfc = fx(4)*x
         if(nf.ge.5) xfb = fx(5)*x
      endif
      if(part.eq.'nu') then
         xfu = (xfu+xfd)/2
         xfd = xfu
         xub = (xub+xdb)/2
         xdb = xub
      endif
      if( sign.eq.'+'.or.
     #   (sign.eq.' '.and.(part.eq.'ph'.or.part.eq.'el')) ) then
      elseif(sign.eq.'-') then
         tmp  = xfu
         xfu  = xub
         xub  = tmp
         tmp  = xfd
         xfd  = xdb
         xdb  = tmp
      else
         write(6,*)'error: unknown charge label', sign
         stop
      endif
      return
      end
