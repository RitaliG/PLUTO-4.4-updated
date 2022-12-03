#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine FLUXDIFFF(
     &           diff
     &           ,idifflo0
     &           ,idiffhi0
     &           ,ndiffcomp
     &           ,F
     &           ,iFlo0
     &           ,iFhi0
     &           ,nFcomp
     &           ,idir
     &           ,iboxlo0
     &           ,iboxhi0
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ndiffcomp
      integer idifflo0
      integer idiffhi0
      REAL_T diff(
     &           idifflo0:idiffhi0,
     &           0:ndiffcomp-1)
      integer nFcomp
      integer iFlo0
      integer iFhi0
      REAL_T F(
     &           iFlo0:iFhi0,
     &           0:nFcomp-1)
      integer idir
      integer iboxlo0
      integer iboxhi0
      integer i0
      integer c2fLo0
      integer c2fHi0
      integer iv
      
      c2fLo0= 0*CHF_ID(0, idir)

      
      c2fHi0= 1*CHF_ID(0, idir)

      do iv = 0,ndiffcomp - 1
         
      do i0 = iboxlo0,iboxhi0

            diff(i0, iv) =
     &        F(i0 +c2fHi0, iv) -
     &        F(i0 +c2fLo0, iv)
         
      enddo
      enddo
      return
      end
