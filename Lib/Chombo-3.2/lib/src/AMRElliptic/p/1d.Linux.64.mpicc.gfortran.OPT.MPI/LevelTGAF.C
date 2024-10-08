#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine RESGHOSTBC(
     &           state
     &           ,istatelo0
     &           ,istatehi0
     &           ,nstatecomp
     &           ,ibcBoxlo0
     &           ,ibcBoxhi0
     &           ,idir
     &           ,side
     &           ,ncomp
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nstatecomp
      integer istatelo0
      integer istatehi0
      REAL_T state(
     &           istatelo0:istatehi0,
     &           0:nstatecomp-1)
      integer ibcBoxlo0
      integer ibcBoxhi0
      integer idir
      integer side
      integer ncomp
      integer nc
      integer ii
      integer i
      REAL_T nearval, farval
      
      ii = side*CHF_ID(0,idir)
      do nc = 0, ncomp-1
         
      do i = ibcBoxlo0,ibcBoxhi0

           nearval = state(i-ii,nc)
           farval  = state(i-2*ii,nc)
           state(i,nc) = two*nearval - farval
         
      enddo
      enddo
      return
      end
      subroutine HORESGHOSTBC(
     &           state
     &           ,istatelo0
     &           ,istatehi0
     &           ,nstatecomp
     &           ,ibcBoxlo0
     &           ,ibcBoxhi0
     &           ,idir
     &           ,side
     &           ,ncomp
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nstatecomp
      integer istatelo0
      integer istatehi0
      REAL_T state(
     &           istatelo0:istatehi0,
     &           0:nstatecomp-1)
      integer ibcBoxlo0
      integer ibcBoxhi0
      integer idir
      integer side
      integer ncomp
      integer nc
      integer ii
      integer i
      REAL_T nearval, midval, farval
      
      ii = side*CHF_ID(0,idir)
      do nc = 0, ncomp-1
         
      do i = ibcBoxlo0,ibcBoxhi0

           nearval = state(i-ii,nc)
           midval  = state(i-2*ii,nc)
           farval  = state(i-3*ii,nc)
           state(i,nc) = three*(nearval - midval) + farval
         
      enddo
      enddo
      return
      end
