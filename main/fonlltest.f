c Driving program for FONLL calculations
      program fonll
      implicit none
c include version string
      include 'version.h'
c prefix for generated files
      character * 30 prefix
c jprefix is last non-blank in prefix
      integer jprefix
      integer ibeam1,ibeam2,nptype1,nptype2,ngroup1,ngroup2,nset1,nset2
     # ,icalctype
      real * 8 ener1,ener2,xmq,xlam,ffact,fren,pt,ylab
      real * 8 hdmassive,hdmassless,hdresummed,hdreserr,sigmafonll,
     #         phmassive,phmassless,phresummed,phreserr
c ww common
      integer itypeww
      real * 8 xmuww,thcww,zminww,zmaxww,elphen
      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
c Unit used for temporary files. Keep closed
c before calling fonll0
      integer itm
      data itm/57/
      prefix='test1fonll'
      jprefix=10
c
      ibeam1=1
      ener1=920
      nptype1=0
      ngroup1=0
      nset1=131
      ibeam2=5
      ener2=27.5
      nptype2=0
      ngroup2=0
      nset2=42
      xmq=1.5
      xlam=-1
      ffact=2
      fren=0.5
      pt=20
      ylab=1
      icalctype=1
      itypeww=1
      xmuww=1
      zminww=0.2
      zmaxww=0.8
c End of input
      call fonll0(icalctype,prefix,jprefix,
     #  ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
      open(11,file=prefix(1:jprefix)//'.outlog',access='append')
c wind file to end
c      call toend(11)
      write(11,'(a)') versionstring
      write(11,111)ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,sigmafonll,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr
      write(11,*)
      close(11)
c
      ibeam1=5
      ener1=27.5
      nptype1=0
      ngroup1=0
      nset1=42
      ibeam2=1
      ener2=920
      nptype2=0
      ngroup2=0
      nset2=131
      xmq=1.5
      xlam=-1
      ffact=2
      fren=0.5
      pt=20
      ylab=-1
      icalctype=1
      itypeww=1
      xmuww=1
      zminww=0.2
      zmaxww=0.8
c End of input
      call fonll0(icalctype,prefix,jprefix,
     #  ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
      open(11,file=prefix(1:jprefix)//'.outlog',access='append')
c wind file to end
c      call toend(11)
      write(11,111)ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,sigmafonll,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr
      write(11,*)
      close(11)
c
      ibeam1=4
      ener1=200
      nptype1=0
      ngroup1=0
      nset1=42
      ibeam2=0
      ener2=0
      nptype2=0
      ngroup2=0
      nset2=131
      xmq=1.5
      xlam=-1
      ffact=0.5
      fren=2
      pt=5
      ylab=2
      icalctype=1
c End of input
      call fonll0(icalctype,prefix,jprefix,
     #  ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
      open(11,file=prefix(1:jprefix)//'.outlog',access='append')
c wind file to end
c      call toend(11)
      write(11,111)ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,sigmafonll,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr
      write(11,*)
      close(11)
c Bottom at Tevatron
      ibeam1=1
      ener1=900
      nptype1=0
      ngroup1=0
      nset1=131
      ibeam2=-1
      ener2=900
      nptype2=0
      ngroup2=0
      nset2=131
      xmq=5
      xlam=-1
      ffact=2
      fren=2
      pt=20
      ylab=1
      icalctype=1
c End of input
      call fonll0(icalctype,prefix,jprefix,
     #  ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
      open(11,file=prefix(1:jprefix)//'.outlog',access='append')
c wind file to end
c      call toend(11)
      write(11,111)ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,sigmafonll,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr
      write(11,*)
      close(11)
c Bottom at Tevatron (980+980)
      ibeam1=1
      ener1=980
      nptype1=0
      ngroup1=0
      nset1=131
      ibeam2=-1
      ener2=980
      nptype2=0
      ngroup2=0
      nset2=131
      xmq=5
      xlam=-1
      ffact=2
      fren=2
      pt=20
      ylab=1
      icalctype=1
c End of input
      call fonll0(icalctype,prefix,jprefix,
     #  ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
      open(11,file=prefix(1:jprefix)//'.outlog',access='append')
c wind file to end
c      call toend(11)
      write(11,111)ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,sigmafonll,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr
      write(11,*)
      close(11)
c charm at Tevatron
      ibeam1=1
      ener1=900
      nptype1=0
      ngroup1=0
      nset1=131
      ibeam2=-1
      ener2=900
      nptype2=0
      ngroup2=0
      nset2=131
      xmq=1.5
      xlam=-1
      ffact=2
      fren=1
      pt=10
      ylab=1
      icalctype=1
c End of input
      call fonll0(icalctype,prefix,jprefix,
     #  ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
      open(11,file=prefix(1:jprefix)//'.outlog',access='append')
c wind file to end
c      call toend(11)
      write(11,111)ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,sigmafonll,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr
      write(11,*)
      close(11)
c charm at fixed target with pion beams
      ibeam1=3
      ener1=800
      nptype1=0
      ngroup1=0
      nset1=32
      ibeam2=0
      ener2=0
      nptype2=0
      ngroup2=0
      nset2=131
      xmq=1.5
      xlam=-1
      ffact=1.5
      fren=2
      pt=3
      ylab=2
      icalctype=1
c End of input
      call fonll0(icalctype,prefix,jprefix,
     #  ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr,sigmafonll)
      open(11,file=prefix(1:jprefix)//'.outlog',access='append')
c wind file to end
c      call toend(11)
      write(11,111)ibeam1,ener1, nptype1, ngroup1, nset1,
     #  ibeam2,ener2, nptype2, ngroup2, nset2,
     #  xmq,xlam,ffact,fren,pt,ylab,sigmafonll,
     #  hdmassive,hdmassless,phmassive,phmassless,
     #  hdresummed,hdreserr,phresummed,phreserr
      write(11,*)
      close(11)
c
 111  format('beam1=',i2,',e1=',d12.6,',pdf1=',i1,',',i1,',',i3,
     #      ',beam2=',i2,',e2=',d12.6,',pdf2=',i1,',',i1,',',i3,/,
     #       'mq=',d12.6,',l5=',d12.6,',ff=',d12.6,',fr=',d12.6,/,
     #       'pt,y,fonll',3(1x,d12.6),/,
     #       'hdmv=',d12.6,',hdml=',d12.6,',phmv=',d12.6,',phml=',d12.6,
     # / ,'hdrs=',d12.4,'+-',d8.2,',phrs=',d12.6,'+-',d8.2)
      end
