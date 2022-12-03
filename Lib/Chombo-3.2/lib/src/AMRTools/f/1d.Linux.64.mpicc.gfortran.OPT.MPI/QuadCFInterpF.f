      subroutine QUADINTERP(
     & phi
     & ,iphilo0
     & ,iphihi0
     & ,nphicomp
     & ,phistar
     & ,iphistarlo0
     & ,iphistarhi0
     & ,nphistarcomp
     & ,iboxlo0
     & ,iboxhi0
     & ,ihilo
     & ,h
     & ,idir
     & ,scomp
     & ,ecomp
     & ,nref
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0,
     & 0:nphicomp-1)
      integer nphistarcomp
      integer iphistarlo0
      integer iphistarhi0
      REAL*8 phistar(
     & iphistarlo0:iphistarhi0,
     & 0:nphistarcomp-1)
      integer iboxlo0
      integer iboxhi0
      integer ihilo
      REAL*8 h
      integer idir
      integer scomp
      integer ecomp
      integer nref
      integer i0
      integer ii0
      integer n
      REAL*8 x, pa, pb, ps, a, b, frac, denom, xsquared
      REAL*8 mult, invh
      frac = (2.0d0)/(h*h)
      denom = nref*nref + 4*nref + 3
      mult = frac/denom
      invh = (1.0d0) / h
      x = (2.0d0) * h
      xsquared = (4.0d0) * h*h
      ii0= ihilo*CHF_ID(0, idir)
      do n=scomp,ecomp
      do i0 = iboxlo0,iboxhi0
            pa = phi(i0 -2*ii0,n)
            pb = phi(i0 -ii0,n)
            ps = phistar(i0 +ii0,n)
            a = mult*((2.0d0)*ps + (nref+1)*pa - (nref+3)*pb)
            b = (pb-pa)*invh - a*h
            phi(i0,n) = xsquared*a + b*x + pa
      enddo
      enddo
      return
      end
      subroutine PHISTAR(
     & fPhiStar
     & ,ifPhiStarlo0
     & ,ifPhiStarhi0
     & ,nfPhiStarcomp
     & ,iregionlo0
     & ,iregionhi0
     & ,phic
     & ,iphiclo0
     & ,iphichi0
     & ,nphiccomp
     & ,coarslope
     & ,icoarslopelo0
     & ,icoarslopehi0
     & ,ncoarslopecomp
     & ,coarcurva
     & ,icoarcurvalo0
     & ,icoarcurvahi0
     & ,ncoarcurvacomp
     & ,coarmixed
     & ,icoarmixedlo0
     & ,icoarmixedhi0
     & ,ncoarmixedcomp
     & ,dxf
     & ,ivar
     & ,dir
     & ,sign
     & ,nRef
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfPhiStarcomp
      integer ifPhiStarlo0
      integer ifPhiStarhi0
      REAL*8 fPhiStar(
     & ifPhiStarlo0:ifPhiStarhi0,
     & 0:nfPhiStarcomp-1)
      integer iregionlo0
      integer iregionhi0
      integer nphiccomp
      integer iphiclo0
      integer iphichi0
      REAL*8 phic(
     & iphiclo0:iphichi0,
     & 0:nphiccomp-1)
      integer ncoarslopecomp
      integer icoarslopelo0
      integer icoarslopehi0
      REAL*8 coarslope(
     & icoarslopelo0:icoarslopehi0,
     & 0:ncoarslopecomp-1)
      integer ncoarcurvacomp
      integer icoarcurvalo0
      integer icoarcurvahi0
      REAL*8 coarcurva(
     & icoarcurvalo0:icoarcurvahi0,
     & 0:ncoarcurvacomp-1)
      integer ncoarmixedcomp
      integer icoarmixedlo0
      integer icoarmixedhi0
      REAL*8 coarmixed(
     & icoarmixedlo0:icoarmixedhi0,
     & 0:ncoarmixedcomp-1)
      REAL*8 dxf
      integer ivar
      integer dir
      integer sign
      integer nRef
      REAL*8 xf1, xc1, xf2, xc2, x1, x2, dxc
      REAL*8 aa, update1, update2, update3
      integer i0
      integer ii0
      integer ir0
      integer ic(0:1 -1)
      integer ivf(0:1 -1)
      integer YOU(1:2, 0:2), you1, you2
      data YOU / 1, 2, 0, 2, 0, 1 /
      dxc = nRef * dxf
      you1 = YOU(1,dir)
      you2 = YOU(2,dir)
      ii0= sign*CHF_ID(0, dir)
      do ir0 = iregionlo0,iregionhi0
         ic(0)=ir0/nRef
         ivf(0)=ir0
         i0=ir0+ii0
         xf1 = (ivf(you1)+(0.500d0))*dxf
         xc1 = ( ic(you1)+(0.500d0))*dxc
         xf2 = (ivf(you2)+(0.500d0))*dxf
         xc2 = ( ic(you2)+(0.500d0))*dxc
         x1 = xf1-xc1
         x2 = xf2-xc2
         aa= phic(ic(0),ivar)
         update1=x1*coarslope(ic(0),you1) +
     & (0.500d0)*x1*x1*coarcurva(ic(0),you1)
         update2=x2*coarslope(ic(0),you2) +
     & (0.500d0)*x2*x2*coarcurva(ic(0),you2)
         update3=x1*x2*coarmixed(ic(0),0)
         fPhiStar(i0,ivar) = aa+update1+update2+update3
      enddo
      return
      end
