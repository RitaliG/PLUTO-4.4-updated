      subroutine AVERAGECODIM(
     & coarse
     & ,icoarselo0,icoarselo1
     & ,icoarsehi0,icoarsehi1
     & ,ncoarsecomp
     & ,fine
     & ,ifinelo0,ifinelo1
     & ,ifinehi0,ifinehi1
     & ,nfinecomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,refRatio
     & ,refFactor
     & ,refDim
     & ,ibreflo0,ibreflo1
     & ,ibrefhi0,ibrefhi1
     & )
      implicit none
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & 0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & 0:nfinecomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer refRatio
      integer refFactor
      integer refDim
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer var
      integer ic0,ic1
      integer ip0,ip1
      integer ii0,ii1
      REAL*8 refScale,coarseSum
      refScale = (1.0d0)/(refFactor**refDim)
      do var = 0, ncoarsecomp - 1
      do ic1 = iboxlo1,iboxhi1
      do ic0 = iboxlo0,iboxhi0
               ip0 = ic0*refRatio
               ip1 = ic1*refRatio
            coarseSum = 0
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0
               coarseSum = coarseSum + fine(ip0+ii0,ip1+ii1,var)
      enddo
      enddo
            coarse(ic0,ic1,var) = coarseSum*refScale
      enddo
      enddo
      enddo
      return
      end
      subroutine AVERAGECODIMHARMONIC(
     & coarse
     & ,icoarselo0,icoarselo1
     & ,icoarsehi0,icoarsehi1
     & ,ncoarsecomp
     & ,fine
     & ,ifinelo0,ifinelo1
     & ,ifinehi0,ifinehi1
     & ,nfinecomp
     & ,iboxlo0,iboxlo1
     & ,iboxhi0,iboxhi1
     & ,refRatio
     & ,refFactor
     & ,refDim
     & ,ibreflo0,ibreflo1
     & ,ibrefhi0,ibrefhi1
     & )
      implicit none
      integer ncoarsecomp
      integer icoarselo0,icoarselo1
      integer icoarsehi0,icoarsehi1
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & icoarselo1:icoarsehi1,
     & 0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0,ifinelo1
      integer ifinehi0,ifinehi1
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & 0:nfinecomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer refRatio
      integer refFactor
      integer refDim
      integer ibreflo0,ibreflo1
      integer ibrefhi0,ibrefhi1
      integer var
      integer ic0,ic1
      integer ip0,ip1
      integer ii0,ii1
      REAL*8 refScale,coarseSum
      refScale = (1.0d0)/(refFactor**refDim)
      do var = 0, ncoarsecomp - 1
      do ic1 = iboxlo1,iboxhi1
      do ic0 = iboxlo0,iboxhi0
               ip0 = ic0*refRatio
               ip1 = ic1*refRatio
            coarseSum = 0
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0
               coarseSum = coarseSum + (1.0d0)/
     & fine(ip0+ii0,ip1+ii1,var)
      enddo
      enddo
            coarse(ic0,ic1,var) = (1.0d0)/(coarseSum*refScale)
      enddo
      enddo
      enddo
      return
      end
