      subroutine MASKVALUE(
     & mask
     & ,imasklo0
     & ,imaskhi0
     & ,test
     & ,itestlo0
     & ,itesthi0
     & ,ibxlo0
     & ,ibxhi0
     & ,val
     & ,onoff
     & )
      implicit none
      integer imasklo0
      integer imaskhi0
      REAL*8 mask(
     & imasklo0:imaskhi0)
      integer itestlo0
      integer itesthi0
      integer test(
     & itestlo0:itesthi0)
      integer ibxlo0
      integer ibxhi0
      integer val
      integer onoff
      integer i0
      if (onoff .eq. 1) then
      do i0 = ibxlo0,ibxhi0
         if (test(i0) .eq. val) then
            mask(i0) = (1.0d0)
         else
            mask(i0) = 0
         endif
      enddo
      else
      do i0 = ibxlo0,ibxhi0
         if (test(i0) .eq. val) then
            mask(i0) = 0
         else
            mask(i0) = (1.0d0)
         endif
      enddo
      endif
      return
      end
