      subroutine RESGHOSTBC(
     & state
     & ,istatelo0,istatelo1
     & ,istatehi0,istatehi1
     & ,nstatecomp
     & ,ibcBoxlo0,ibcBoxlo1
     & ,ibcBoxhi0,ibcBoxhi1
     & ,idir
     & ,side
     & ,ncomp
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL*8 state(
     & istatelo0:istatehi0,
     & istatelo1:istatehi1,
     & 0:nstatecomp-1)
      integer ibcBoxlo0,ibcBoxlo1
      integer ibcBoxhi0,ibcBoxhi1
      integer idir
      integer side
      integer ncomp
      integer nc
      integer ii,jj
      integer i,j
      REAL*8 nearval, farval
      ii = side*CHF_ID(0,idir)
      jj = side*CHF_ID(1,idir)
      do nc = 0, ncomp-1
      do j = ibcBoxlo1,ibcBoxhi1
      do i = ibcBoxlo0,ibcBoxhi0
           nearval = state(i-ii,j-jj,nc)
           farval = state(i-2*ii,j-2*jj,nc)
           state(i,j,nc) = (2.0d0)*nearval - farval
      enddo
      enddo
      enddo
      return
      end
      subroutine HORESGHOSTBC(
     & state
     & ,istatelo0,istatelo1
     & ,istatehi0,istatehi1
     & ,nstatecomp
     & ,ibcBoxlo0,ibcBoxlo1
     & ,ibcBoxhi0,ibcBoxhi1
     & ,idir
     & ,side
     & ,ncomp
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nstatecomp
      integer istatelo0,istatelo1
      integer istatehi0,istatehi1
      REAL*8 state(
     & istatelo0:istatehi0,
     & istatelo1:istatehi1,
     & 0:nstatecomp-1)
      integer ibcBoxlo0,ibcBoxlo1
      integer ibcBoxhi0,ibcBoxhi1
      integer idir
      integer side
      integer ncomp
      integer nc
      integer ii,jj
      integer i,j
      REAL*8 nearval, midval, farval
      ii = side*CHF_ID(0,idir)
      jj = side*CHF_ID(1,idir)
      do nc = 0, ncomp-1
      do j = ibcBoxlo1,ibcBoxhi1
      do i = ibcBoxlo0,ibcBoxhi0
           nearval = state(i-ii,j-jj,nc)
           midval = state(i-2*ii,j-2*jj,nc)
           farval = state(i-3*ii,j-3*jj,nc)
           state(i,j,nc) = (3.0d0)*(nearval - midval) + farval
      enddo
      enddo
      enddo
      return
      end
