      implicit none
      integer n
      parameter        (n=6)
      real * 8 fx(-6:6)
      real * 8 ee
c ww common
      integer itypeww
      real * 8 xmuww,thcww,zminww,zmaxww,elphen
      common/fonllww/xmuww,thcww,zminww,zmaxww,elphen,itypeww
      integer j,k,l
      real * 8 xs,xs1
      itypeww=1
      xmuww=1
      zminww=0.2
      zmaxww=0.8
      xs=0.1
      call fonllmlmpdf(1,5,1.d0,xs,fx,6)
      do k=0,n
         call getarrs(k,0,xs,ee)
         call fonllmlmpdf(1,5,1.d0,xs,fx,6)
         do j=-6,6
            call getarrs(k,j,xs,ee)
            if(abs(xs*fx(j)-ee).ne.0) then
               write(*,*) ' error',xs*fx(j)/ee,k,j
               read(*,*)
            endif
         enddo
      enddo
      write(*,*) ' done 1'
      write(11,*) 'set scale x log y log'
      do k=0,n-1
         do l=1,10
            call getarrs(k,0,xs,ee)
            call getarrs(k+1,0,xs1,ee)
            xs=xs+(l-1)*(xs1-xs)/10
            call fonllmlmpdf(1,5,1.d0,xs,fx,6)
            write(11,*) xs,fx(0)
         enddo
      enddo
      write(11,*) 'join'
      write(11,*) 'join'
      write(11,*) ' set symbol 5O'
      do k=0,n
         call getarrs(k,0,xs,ee)
         write(11,*) xs,ee/xs
      enddo
      write(11,*) 'plot'
      end

