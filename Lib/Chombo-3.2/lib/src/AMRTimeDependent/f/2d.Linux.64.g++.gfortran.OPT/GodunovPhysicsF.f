      subroutine FLUXDIFFF(
     & diff
     & ,idifflo0,idifflo1
     & ,idiffhi0,idiffhi1
     & ,ndiffcomp
     & ,F
     & ,iFlo0,iFlo1
     & ,iFhi0,iFhi1
     & ,nFcomp
     & ,idir
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndiffcomp
      integer idifflo0,idifflo1
      integer idiffhi0,idiffhi1
      REAL*8 diff(
     & idifflo0:idiffhi0,
     & idifflo1:idiffhi1,
     & 0:ndiffcomp-1)
      integer nFcomp
      integer iFlo0,iFlo1
      integer iFhi0,iFhi1
      REAL*8 F(
     & iFlo0:iFhi0,
     & iFlo1:iFhi1,
     & 0:nFcomp-1)
      integer idir
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer i0,i1
      integer c2fLo0,c2fLo1
      integer c2fHi0,c2fHi1
      integer iv
      c2fLo0= 0*CHF_ID(0, idir)
      c2fLo1= 0*CHF_ID(1, idir)
      c2fHi0= 1*CHF_ID(0, idir)
      c2fHi1= 1*CHF_ID(1, idir)
      do iv = 0,ndiffcomp - 1
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            diff(i0,i1, iv) =
     & F(i0 +c2fHi0,i1 +c2fHi1, iv) -
     & F(i0 +c2fLo0,i1 +c2fLo1, iv)
      enddo
      enddo
      enddo
      return
      end
