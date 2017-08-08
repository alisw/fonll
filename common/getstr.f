
      subroutine getstr(str)
c------------------------------------------------------------
c Turns to lowercase, and replaces ascii tabs with blanks.
c
      character * (*) str
      logical init
      data init/.true./, it/9/
c acii tab is 9
      if(init) then
         ia = ichar('A')
         iz = ichar('Z')
         ilo = ichar('z') - iz
         init=.false.
      endif
      read(55,'(a)') str
      do j=1,len(str)
      i = ichar(str(j:j))
      if(i.eq.it) then
        str(j:j)=' '
      elseif(i.ge.ia.and.i.le.iz) then
        str(j:j)=char(i+ilo)
      endif
      enddo
      return
      end

