#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine QUADINTERP(
     &           phi
     &           ,iphilo0,iphilo1
     &           ,iphihi0,iphihi1
     &           ,nphicomp
     &           ,phistar
     &           ,iphistarlo0,iphistarlo1
     &           ,iphistarhi0,iphistarhi1
     &           ,nphistarcomp
     &           ,iboxlo0,iboxlo1
     &           ,iboxhi0,iboxhi1
     &           ,ihilo
     &           ,h
     &           ,idir
     &           ,scomp
     &           ,ecomp
     &           ,nref
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphicomp
      integer iphilo0,iphilo1
      integer iphihi0,iphihi1
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           0:nphicomp-1)
      integer nphistarcomp
      integer iphistarlo0,iphistarlo1
      integer iphistarhi0,iphistarhi1
      REAL_T phistar(
     &           iphistarlo0:iphistarhi0,
     &           iphistarlo1:iphistarhi1,
     &           0:nphistarcomp-1)
      integer iboxlo0,iboxlo1
      integer iboxhi0,iboxhi1
      integer ihilo
      REAL_T h
      integer idir
      integer scomp
      integer ecomp
      integer nref
      integer i0,i1
      integer ii0,ii1
      integer n
      real_t x, pa, pb, ps, a, b, frac, denom, xsquared
      real_t  mult, invh
      frac = two/(h*h)
      denom =  nref*nref + 4*nref + 3
      mult = frac/denom
      invh = one / h
      x = two * h
      xsquared = four * h*h
      
      ii0= ihilo*CHF_ID(0, idir)

      ii1= ihilo*CHF_ID(1, idir)

      do n=scomp,ecomp
         
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

            pa = phi(i0 -2*ii0,i1 -2*ii1,n)
            pb = phi(i0 -ii0,i1 -ii1,n)
            ps = phistar(i0 +ii0,i1 +ii1,n)
            a  = mult*(two*ps + (nref+1)*pa - (nref+3)*pb)
            b  = (pb-pa)*invh - a*h
            phi(i0,i1,n) = xsquared*a + b*x + pa
         
      enddo
      enddo
      enddo
      return
      end
      subroutine PHISTAR(
     &           fPhiStar
     &           ,ifPhiStarlo0,ifPhiStarlo1
     &           ,ifPhiStarhi0,ifPhiStarhi1
     &           ,nfPhiStarcomp
     &           ,iregionlo0,iregionlo1
     &           ,iregionhi0,iregionhi1
     &           ,phic
     &           ,iphiclo0,iphiclo1
     &           ,iphichi0,iphichi1
     &           ,nphiccomp
     &           ,coarslope
     &           ,icoarslopelo0,icoarslopelo1
     &           ,icoarslopehi0,icoarslopehi1
     &           ,ncoarslopecomp
     &           ,coarcurva
     &           ,icoarcurvalo0,icoarcurvalo1
     &           ,icoarcurvahi0,icoarcurvahi1
     &           ,ncoarcurvacomp
     &           ,coarmixed
     &           ,icoarmixedlo0,icoarmixedlo1
     &           ,icoarmixedhi0,icoarmixedhi1
     &           ,ncoarmixedcomp
     &           ,dxf
     &           ,ivar
     &           ,dir
     &           ,sign
     &           ,nRef
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfPhiStarcomp
      integer ifPhiStarlo0,ifPhiStarlo1
      integer ifPhiStarhi0,ifPhiStarhi1
      REAL_T fPhiStar(
     &           ifPhiStarlo0:ifPhiStarhi0,
     &           ifPhiStarlo1:ifPhiStarhi1,
     &           0:nfPhiStarcomp-1)
      integer iregionlo0,iregionlo1
      integer iregionhi0,iregionhi1
      integer nphiccomp
      integer iphiclo0,iphiclo1
      integer iphichi0,iphichi1
      REAL_T phic(
     &           iphiclo0:iphichi0,
     &           iphiclo1:iphichi1,
     &           0:nphiccomp-1)
      integer ncoarslopecomp
      integer icoarslopelo0,icoarslopelo1
      integer icoarslopehi0,icoarslopehi1
      REAL_T coarslope(
     &           icoarslopelo0:icoarslopehi0,
     &           icoarslopelo1:icoarslopehi1,
     &           0:ncoarslopecomp-1)
      integer ncoarcurvacomp
      integer icoarcurvalo0,icoarcurvalo1
      integer icoarcurvahi0,icoarcurvahi1
      REAL_T coarcurva(
     &           icoarcurvalo0:icoarcurvahi0,
     &           icoarcurvalo1:icoarcurvahi1,
     &           0:ncoarcurvacomp-1)
      integer ncoarmixedcomp
      integer icoarmixedlo0,icoarmixedlo1
      integer icoarmixedhi0,icoarmixedhi1
      REAL_T coarmixed(
     &           icoarmixedlo0:icoarmixedhi0,
     &           icoarmixedlo1:icoarmixedhi1,
     &           0:ncoarmixedcomp-1)
      REAL_T dxf
      integer ivar
      integer dir
      integer sign
      integer nRef
      REAL_T xf1, xc1, xf2, xc2, x1, x2, dxc
      REAL_T aa, update1, update2, update3
      integer i0,i1
      integer ii0,ii1
      integer ir0,ir1
      integer ic(0:CH_SPACEDIM-1)
      integer ivf(0:CH_SPACEDIM-1)
      integer YOU(1:2, 0:2), you1, you2
      data YOU / 1, 2, 0, 2, 0, 1 /
#if CH_SPACEDIM > 3
      call MAYDAY_ERROR()
#else
      dxc = nRef * dxf
      you1 = YOU(1,dir)
      you2 = YOU(2,dir)
      
      ii0= sign*CHF_ID(0, dir)

      ii1= sign*CHF_ID(1, dir)

      
      do ir1 = iregionlo1,iregionhi1
      do ir0 = iregionlo0,iregionhi0

         
         ic(0)=ir0/nRef
         ic(1)=ir1/nRef
         
         ivf(0)=ir0
         ivf(1)=ir1
         
         i0=ir0+ii0
         i1=ir1+ii1
         xf1 = (ivf(you1)+half)*dxf
         xc1 = (  ic(you1)+half)*dxc
         xf2 = (ivf(you2)+half)*dxf
         xc2 = (  ic(you2)+half)*dxc
         x1 = xf1-xc1
         x2 = xf2-xc2
         aa= phic(ic(0),ic(1),ivar)
         update1=x1*coarslope(ic(0),ic(1),you1) +
     &        half*x1*x1*coarcurva(ic(0),ic(1),you1)
         update2=x2*coarslope(ic(0),ic(1),you2) +
     &        half*x2*x2*coarcurva(ic(0),ic(1),you2)
         update3=x1*x2*coarmixed(ic(0),ic(1),0)
         fPhiStar(i0,i1,ivar) = aa+update1+update2+update3
      
      enddo
      enddo
#endif
      return
      end
