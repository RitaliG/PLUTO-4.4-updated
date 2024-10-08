#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine NODEINTERPMG_GETWEIGHTS(
     &           nref
     &           ,ibreflo0
     &           ,ibrefhi0
     &           ,wtcrnr
     &           ,iwtcrnrlo0
     &           ,iwtcrnrhi0
     &           ,nwtcrnrcomp
     &           )

      implicit none
      integer nref
      integer ibreflo0
      integer ibrefhi0
      integer nwtcrnrcomp
      integer iwtcrnrlo0
      integer iwtcrnrhi0
      REAL_T wtcrnr(
     &           iwtcrnrlo0:iwtcrnrhi0,
     &           0:nwtcrnrcomp-1)
      integer iref0
      integer ib0
      integer nvwt
      integer ibmax0
      REAL_T refinv, wt
      REAL_T fraci0
      REAL_T wti0
      refinv = one / nref
      nvwt = 0
      
      do iref0 = ibreflo0,ibrefhi0

         
         call maxb(iref0, ibmax0) 
         
         fraci0 = iref0 * refinv 
         
         do ib0 = 0, ibmax0 
            wt = one
            
            call wtside(ib0, fraci0, wti0)
            wt = wt * wti0 
            wtcrnr( ib0, nvwt ) = wt
         
         end do 
         nvwt = nvwt + 1
      
      enddo
      return
      end
      subroutine NODEINTERPMG(
     &           fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,crse
     &           ,icrselo0
     &           ,icrsehi0
     &           ,ncrsecomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,nref
     &           ,ibreflo0
     &           ,ibrefhi0
     &           ,wtcrnr
     &           ,iwtcrnrlo0
     &           ,iwtcrnrhi0
     &           ,nwtcrnrcomp
     &           )

      implicit none
      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer ncrsecomp
      integer icrselo0
      integer icrsehi0
      REAL_T crse(
     &           icrselo0:icrsehi0,
     &           0:ncrsecomp-1)
      integer iregionlo0
      integer iregionhi0
      integer nref
      integer ibreflo0
      integer ibrefhi0
      integer nwtcrnrcomp
      integer iwtcrnrlo0
      integer iwtcrnrhi0
      REAL_T wtcrnr(
     &           iwtcrnrlo0:iwtcrnrhi0,
     &           0:nwtcrnrcomp-1)
      integer iref0, icrse0
      integer ifine0, ib0;
      integer var, ncomp, nvwt
      integer ibmax0
      integer icmin0
      integer icmax0;
      REAL_T csum, refinv
      ncomp = nfinecomp
      if (ncomp .ne. ncrsecomp) then
         print *, 'nodeinterpmg: fine and crse incompatible'
         call MAYDAY_ERROR()
      endif
      refinv = one / nref
      
      icmin0 = iregionlo0 
      nvwt = 0
      
      do iref0 = ibreflo0,ibrefhi0

         
         call maxb(iref0, ibmax0) 
         
         icmax0 = iregionhi0 + (1-ibmax0) 
         
         do icrse0 = icmin0, icmax0 
            
            ifine0 = nref*icrse0 + iref0 
            do var = 0, ncomp-1
               csum = 0
               
               do ib0 = 0, ibmax0 
                  csum = csum + wtcrnr( ib0, nvwt ) *
     &              crse( icrse0+ib0, var)
               
               end do 
               fine ( ifine0, var ) = csum +
     &              fine ( ifine0, var )
            end do
         
         end do 
         nvwt = nvwt + 1
      
      enddo
      return
      end
      subroutine WTSIDE(
     &           i
     &           ,frac
     &           ,wt
     &           )

      implicit none
      integer i
      REAL_T frac
      REAL_T wt
      if (i .eq. 0) then
         wt = one - frac
      else
         wt = frac
      endif
      return
      end
      subroutine MAXB(
     &           iref
     &           ,ibmax
     &           )

      implicit none
      integer iref
      integer ibmax
      if (iref .eq. 0) then
         ibmax = 0
      else
         ibmax = 1
      endif
      return
      end
