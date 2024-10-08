#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine INTERPFACECONSTANT(
     &           fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,coarse
     &           ,icoarselo0
     &           ,icoarsehi0
     &           ,ncoarsecomp
     &           ,iblo0
     &           ,ibhi0
     &           ,ref_ratio
     &           ,ibreflo0
     &           ,ibrefhi0
     &           ,dir
     &           )

      implicit none
      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           0:ncoarsecomp-1)
      integer iblo0
      integer ibhi0
      integer ref_ratio
      integer ibreflo0
      integer ibrefhi0
      integer dir
      integer var
      integer ic0
      integer ifine0
      integer ii0
      do var = 0, ncoarsecomp - 1
         
      do ic0 = iblo0,ibhi0

            
      do ii0 = ibreflo0,ibrefhi0

            
               ifine0 = ic0*ref_ratio + ii0
               fine(ifine0,var) = coarse(ic0,var)
            
      enddo
         
      enddo
      end do
      return
      end
      subroutine INTERPLINEARFACE(
     &           fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,slope
     &           ,islopelo0
     &           ,islopehi0
     &           ,nslopecomp
     &           ,iblo0
     &           ,ibhi0
     &           ,dir
     &           ,ref_ratio
     &           ,ibreffacelo0
     &           ,ibreffacehi0
     &           )

      implicit none
      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer nslopecomp
      integer islopelo0
      integer islopehi0
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           0:nslopecomp-1)
      integer iblo0
      integer ibhi0
      integer dir
      integer ref_ratio
      integer ibreffacelo0
      integer ibreffacehi0
      integer  ic 0
      integer  ifine 0
      integer  ii 0
      integer var, id
      REAL_T dxf
      do var = 0, nfinecomp - 1
         
      do ic0 = iblo0,ibhi0

              
      do ii0 = ibreffacelo0,ibreffacehi0

              
                  ifine0 = ic0*ref_ratio + ii0
              
                  if (dir .eq. 0) then
                      id = ii0
                  endif
                  dxf = -half + ( (id+half) / ref_ratio )
                  fine( ifine0,var) =
     &                 fine( ifine0,var) +
     &                 dxf * slope (   ic 0, var )
              
      enddo
          
      enddo
      end do
      return
      end
      subroutine INTERPLINEARINTERIORFACE(
     &           fine
     &           ,ifinelo0
     &           ,ifinehi0
     &           ,nfinecomp
     &           ,ibcoarselo0
     &           ,ibcoarsehi0
     &           ,ref_ratio
     &           ,facedir
     &           ,iinteriorrefboxlo0
     &           ,iinteriorrefboxhi0
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfinecomp
      integer ifinelo0
      integer ifinehi0
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           0:nfinecomp-1)
      integer ibcoarselo0
      integer ibcoarsehi0
      integer ref_ratio
      integer facedir
      integer iinteriorrefboxlo0
      integer iinteriorrefboxhi0
      integer ic0
      integer ifine0
      integer ii0
      integer iloface0
      integer ihiface0
      integer var, id
      REAL_T dxf, diff
      REAL_T loval, hival
      do var=0, nfinecomp -1
         
      do ic0 = ibcoarselo0,ibcoarsehi0

            
      do ii0 = iinteriorrefboxlo0,iinteriorrefboxhi0

            
              ifine0 = ic0*ref_ratio + ii0
              
              iloface0 = ic0*ref_ratio + (1-CHF_ID(0,facedir))*ii0
              
              ihiface0 = iloface0 + ref_ratio*CHF_ID(0,facedir)
              
              if (facedir .eq. 0) then
                 id = ii0
              endif
              dxf = float(id)/ref_ratio
              diff = fine(ihiface0,var)
     &                -fine(iloface0,var)
              fine( ifine0,var) =
     &            fine(iloface0,var)
     &           +dxf*diff
            
      enddo
          
      enddo
       enddo
       return
       end
      subroutine INTERPLIMITFACE(
     &           islope
     &           ,iislopelo0
     &           ,iislopehi0
     &           ,nislopecomp
     &           ,jslope
     &           ,ijslopelo0
     &           ,ijslopehi0
     &           ,njslopecomp
     &           ,kslope
     &           ,ikslopelo0
     &           ,ikslopehi0
     &           ,nkslopecomp
     &           ,lslope
     &           ,ilslopelo0
     &           ,ilslopehi0
     &           ,nlslopecomp
     &           ,mslope
     &           ,imslopelo0
     &           ,imslopehi0
     &           ,nmslopecomp
     &           ,nslope
     &           ,inslopelo0
     &           ,inslopehi0
     &           ,nnslopecomp
     &           ,state
     &           ,istatelo0
     &           ,istatehi0
     &           ,nstatecomp
     &           ,iblo0
     &           ,ibhi0
     &           ,ibnlo0
     &           ,ibnhi0
     &           ,ivalidBoxlo0
     &           ,ivalidBoxhi0
     &           ,normaldir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nislopecomp
      integer iislopelo0
      integer iislopehi0
      REAL_T islope(
     &           iislopelo0:iislopehi0,
     &           0:nislopecomp-1)
      integer njslopecomp
      integer ijslopelo0
      integer ijslopehi0
      REAL_T jslope(
     &           ijslopelo0:ijslopehi0,
     &           0:njslopecomp-1)
      integer nkslopecomp
      integer ikslopelo0
      integer ikslopehi0
      REAL_T kslope(
     &           ikslopelo0:ikslopehi0,
     &           0:nkslopecomp-1)
      integer nlslopecomp
      integer ilslopelo0
      integer ilslopehi0
      REAL_T lslope(
     &           ilslopelo0:ilslopehi0,
     &           0:nlslopecomp-1)
      integer nmslopecomp
      integer imslopelo0
      integer imslopehi0
      REAL_T mslope(
     &           imslopelo0:imslopehi0,
     &           0:nmslopecomp-1)
      integer nnslopecomp
      integer inslopelo0
      integer inslopehi0
      REAL_T nslope(
     &           inslopelo0:inslopehi0,
     &           0:nnslopecomp-1)
      integer nstatecomp
      integer istatelo0
      integer istatehi0
      REAL_T state(
     &           istatelo0:istatehi0,
     &           0:nstatecomp-1)
      integer iblo0
      integer ibhi0
      integer ibnlo0
      integer ibnhi0
      integer ivalidBoxlo0
      integer ivalidBoxhi0
      integer normaldir
      integer   i 0, var
      integer   ii 0
      integer   in 0
      REAL_T statemax, statemin, deltasum, etamax, etamin, eta
      REAL_T tempone, tempzero
      tempone = one
      tempzero = 0
      do var = 0, nislopecomp - 1
         
      do i0 = iblo0,ibhi0

             statemax = state ( i0, var )
             statemin = state ( i0, var )
             
      do ii0 = ibnlo0,ibnhi0

             
                 in0 = i0 + ii0
                 if (
                 
     &                in0 .ge. ivalidBoxlo0 .and.
     &                in0 .le. ivalidBoxhi0 
     &                )
     &        then
                    statemax = max ( statemax, state(in0,var))
                    statemin = min ( statemin, state(in0,var))
                 endif
             
      enddo
             deltasum = half * (
                
     &            (1-CHF_ID(normaldir,0))*abs(islope(i0,var))
     &            )
             if ( deltasum .gt. 0 ) then
                etamax = ( statemax - state ( i0, var ) )
     &               / deltasum
                etamin = ( state ( i0, var ) - statemin )
     &               / deltasum
                eta = max ( min ( etamin, etamax, tempone ), tempzero )
                
                islope ( i0, var ) =
     &               eta * islope ( i0, var ) 
             end if
         
      enddo
      end do
      return
      end
