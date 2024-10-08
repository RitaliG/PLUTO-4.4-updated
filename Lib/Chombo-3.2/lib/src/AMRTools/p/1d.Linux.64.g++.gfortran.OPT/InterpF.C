#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine INTERPCONSTANT(
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
     &           ,dx
     &           ,stretch
     &           ,xbeg
     &           ,ixbeghi0
     &           ,ibreflo0
     &           ,ibrefhi0
     &           ,geometry
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
      REAL_T dx
      REAL_T stretch
      integer ixbeghi0
      REAL_T xbeg(
     &           0:ixbeghi0)
      integer ibreflo0
      integer ibrefhi0
      integer geometry
      integer var
      integer ic0
      integer if0
      integer ii0
      REAL_T  xlf0
      REAL_T  xrf0
      REAL_T  volume
      if (geometry .eq. 1) then
      do var = 0, ncoarsecomp - 1
         
      do ic0 = iblo0,ibhi0

            
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               fine(if0,var) = coarse(ic0,var)
            
      enddo
         
      enddo
      end do
      endif
      if ((geometry .eq. 2) .or. (geometry .eq. 5)) then
      do var = 0, ncoarsecomp - 1
         
      do ic0 = iblo0,ibhi0

            
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               volume = abs(xbeg(0)/dx+if0+half)
               fine(if0,var) = coarse(ic0,var)*volume
            
      enddo
         
      enddo
      end do
      endif
      if (geometry .eq. 3) then
      do var = 0, ncoarsecomp - 1
         
      do ic0 = iblo0,ibhi0

            
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               xlf0 = xbeg(0)+if0*dx
               xrf0 = xlf0 + dx
               volume = (xrf0*xrf0*xrf0-xlf0*xlf0*xlf0)*third
               fine(if0,var) = coarse(ic0,var)*volume
            
      enddo
         
      enddo
      end do
      endif
      if (geometry .eq. 4) then
      do var = 0, ncoarsecomp - 1
         
      do ic0 = iblo0,ibhi0

            
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               xlf0 = if0*dx
               xrf0 = xlf0 + dx
               volume = (exp(three*xrf0)-exp(three*xlf0))*third
               fine(if0,var) = coarse(ic0,var)*volume
            
      enddo
         
      enddo
      end do
      endif
      if (geometry .eq. 6) then
      do var = 0, ncoarsecomp - 1
         
      do ic0 = iblo0,ibhi0

            
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               xlf0 = if0*dx
               xrf0 = xlf0 + dx
               volume = (exp(two*xrf0)-exp(two*xlf0))*half
               fine(if0,var) = coarse(ic0,var)*volume
            
      enddo
         
      enddo
      end do
      endif
      return
      end
      subroutine INTERPCENTRALSLOPE(
     &           slope
     &           ,islopelo0
     &           ,islopehi0
     &           ,nslopecomp
     &           ,state
     &           ,istatelo0
     &           ,istatehi0
     &           ,nstatecomp
     &           ,iblo0
     &           ,ibhi0
     &           ,dir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopecomp
      integer islopelo0
      integer islopehi0
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           0:nslopecomp-1)
      integer nstatecomp
      integer istatelo0
      integer istatehi0
      REAL_T state(
     &           istatelo0:istatehi0,
     &           0:nstatecomp-1)
      integer iblo0
      integer ibhi0
      integer dir
      integer i0
      integer ii0
      integer var
      
      ii0=CHF_ID(0, dir)

      do var = 0, nstatecomp - 1
         
      do i0 = iblo0,ibhi0

          slope (i0,var) = half * (
     &        state (i0+ii0,var) -
     &        state (i0-ii0,var) )
          
      enddo
       end do
      return
      end
      subroutine INTERPHISIDESLOPE(
     &           slope
     &           ,islopelo0
     &           ,islopehi0
     &           ,nslopecomp
     &           ,state
     &           ,istatelo0
     &           ,istatehi0
     &           ,nstatecomp
     &           ,iblo0
     &           ,ibhi0
     &           ,dir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopecomp
      integer islopelo0
      integer islopehi0
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           0:nslopecomp-1)
      integer nstatecomp
      integer istatelo0
      integer istatehi0
      REAL_T state(
     &           istatelo0:istatehi0,
     &           0:nstatecomp-1)
      integer iblo0
      integer ibhi0
      integer dir
      integer i0
      integer ii0
      integer var
      
      ii0=CHF_ID(0, dir)

      do var = 0, nstatecomp - 1
         
      do i0 = iblo0,ibhi0

          slope (i0,var) =
     &          state ( i0+ii0, var)
     &        - state ( i0, var)
          
      enddo
       enddo
      return
      end
      subroutine INTERPLOSIDESLOPE(
     &           slope
     &           ,islopelo0
     &           ,islopehi0
     &           ,nslopecomp
     &           ,state
     &           ,istatelo0
     &           ,istatehi0
     &           ,nstatecomp
     &           ,iblo0
     &           ,ibhi0
     &           ,dir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopecomp
      integer islopelo0
      integer islopehi0
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           0:nslopecomp-1)
      integer nstatecomp
      integer istatelo0
      integer istatehi0
      REAL_T state(
     &           istatelo0:istatehi0,
     &           0:nstatecomp-1)
      integer iblo0
      integer ibhi0
      integer dir
      integer i0, var
      integer ii0
      
      ii0=CHF_ID(0, dir)

      do var = 0, nstatecomp - 1
         
      do i0 = iblo0,ibhi0

         slope (i0,var) =
     &        state (  i 0, var) -
     &        state (  i0-ii0, var)
          
      enddo
       end do
      return
      end
      subroutine INTERPLIMIT(
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
     &           ,state
     &           ,istatelo0
     &           ,istatehi0
     &           ,nstatecomp
     &           ,ibcoarselo0
     &           ,ibcoarsehi0
     &           ,ibnlo0
     &           ,ibnhi0
     &           ,iphysdomainlo0
     &           ,iphysdomainhi0
     &           )

      implicit none
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
      integer nstatecomp
      integer istatelo0
      integer istatehi0
      REAL_T state(
     &           istatelo0:istatehi0,
     &           0:nstatecomp-1)
      integer ibcoarselo0
      integer ibcoarsehi0
      integer ibnlo0
      integer ibnhi0
      integer iphysdomainlo0
      integer iphysdomainhi0
      integer  i0, var
      integer  ii0
      integer  in0
      REAL_T statemax, statemin, deltasum,  eta
      do var = 0, nislopecomp - 1
         
      do i0 = ibcoarselo0,ibcoarsehi0

             statemax = state ( i0, var )
             statemin = state ( i0, var )
             
      do ii0 = ibnlo0,ibnhi0

             
                 in0 = i0 + ii0
                 if (
                 
     &                in0 .ge. iphysdomainlo0 .and.
     &                in0 .le. iphysdomainhi0 
     &                ) then
                    statemax = max ( statemax, state(in0,var))
                    statemin = min ( statemin, state(in0,var))
                 endif
             
      enddo
             deltasum = half * (
                
     &            abs ( islope ( i0, var ) )
     &            )
#if CH_SPACEDIM > 3
                call MAYDAY_ERROR()
#endif
              eta = min(statemax - state(i0,var),
     &                  state(i0,var) - statemin)
              if( eta .le. 1.e-9*abs(statemax) ) then
                 eta = zero
              else
              if (deltasum .gt. eta) then
                eta = eta/deltasum
              else
                eta = one
              endif
              endif
              
              islope ( i0, var ) =
     &             eta * islope ( i0, var ) 
#if CH_SPACEDIM > 3
              call MAYDAY_ERROR()
#endif
         
      enddo
      end do
      return
      end
      subroutine INTERPLINEAR(
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
     &           ,dx
     &           ,stretch
     &           ,xbeg
     &           ,ixbeghi0
     &           ,ibreflo0
     &           ,ibrefhi0
     &           ,geometry
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
      REAL_T dx
      REAL_T stretch
      integer ixbeghi0
      REAL_T xbeg(
     &           0:ixbeghi0)
      integer ibreflo0
      integer ibrefhi0
      integer geometry
      integer ic0
      integer if0
      integer ii0
      integer var, id
      REAL_T  xlf0
      REAL_T  xrf0
      REAL_T  xlc0
      REAL_T  xrc0
      REAL_T dxf, volume
      if (geometry .eq. 1) then 
      do var = 0, nfinecomp - 1
          
      do ic0 = iblo0,ibhi0

              
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
              
                  if (dir .eq. 0) then
                      id = ii0
                  endif
                  dxf = -half + ( (id+half) / ref_ratio )
                  fine( if0,var) =
     &                 fine( if0,var) +
     &                 dxf * slope (   ic 0, var )
              
      enddo
          
      enddo
      end do
      end if
      if ((geometry .eq. 2) .or. (geometry .eq. 5)) then
      do var = 0, nfinecomp - 1
          
      do ic0 = iblo0,ibhi0

               xlc0 = xbeg(0)+ic0*ref_ratio*dx
               xrc0 = xlc0 + ref_ratio*dx;
              
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  volume = abs(xbeg(0)/dx+if0+half)
              
                  if (dir .eq. 0) then
                      xlf0 = xbeg(0)+if0*dx
                      xrf0 = xlf0 + dx
                      dxf = half*(xlf0*xlf0+xrf0*xrf0-xlc0*xlc0-xrc0*xrc0)/
     &                           (xrc0*xrc0-xlc0*xlc0)
                  endif
                  fine( if0,var) =
     &                 fine( if0,var) +
     &                 dxf * slope (   ic 0, var )*volume
              
      enddo
          
      enddo
      end do
      endif
      if (geometry .eq. 3) then
      do var = 0, nfinecomp - 1
          
      do ic0 = iblo0,ibhi0

              
               xlc0 = xbeg(0)+ic0*ref_ratio*dx
               xrc0 = xlc0 + ref_ratio*dx
              
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  xlf0 = xbeg(0)+if0*dx
                  xrf0 = xlf0 + dx
                  volume = (xrf0*xrf0*xrf0-xlf0*xlf0*xlf0)*third
              
                  if (dir .eq. 0) then
                      dxf = half*(xlf0*xlf0*xlf0+xrf0*xrf0*xrf0-xlc0*xlc0*xlc0-xrc0*xrc0*xrc0)/
     &                           (xrc0*xrc0*xrc0-xlc0*xlc0*xlc0)
                  endif
                  fine( if0,var) =
     &                 fine( if0,var) +
     &                 dxf * slope (   ic 0, var )*volume
              
      enddo
          
      enddo
      end do
      endif
      if (geometry .eq. 4) then
      do var = 0, nfinecomp - 1
          
      do ic0 = iblo0,ibhi0

              
               xlc0 = ic0*ref_ratio*dx
               xrc0 = xlc0 + ref_ratio*dx
              
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  xlf0 = if0*dx
                  xrf0 = xlf0 + dx
                  volume = (exp(three*xrf0)-exp(three*xlf0))*third
              
                  if (dir .eq. 0) then
                      dxf = half*(exp(three*xlf0)+exp(three*xrf0)-exp(three*xlc0)-exp(three*xrc0))/
     &                           (exp(three*xrc0)-exp(three*xlc0))
                  endif
                  fine( if0,var) =
     &                 fine( if0,var) +
     &                 dxf * slope (   ic 0, var )*volume
              
      enddo
          
      enddo
      end do
      endif
      if (geometry .eq. 6) then
      do var = 0, nfinecomp - 1
          
      do ic0 = iblo0,ibhi0

               xlc0 = ic0*ref_ratio*dx
               xrc0 = xlc0 + ref_ratio*dx;
              
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  xlf0 = if0*dx
                  xrf0 = xlf0 + dx
                  volume = (exp(two*xrf0)-exp(two*xlf0))*half
              
                  if (dir .eq. 0) then
                      dxf = half*(exp(two*xlf0)+exp(two*xrf0)-exp(two*xlc0)-exp(two*xrc0))/
     &                           (exp(two*xrc0)-exp(two*xlc0))
                  endif
                  fine( if0,var) =
     &                 fine( if0,var) +
     &                 dxf * slope (   ic 0, var )*volume
              
      enddo
          
      enddo
      end do
      endif
      return
      end
      subroutine INTERPHOMO_OLD(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,x1
     &           ,dxCrse
     &           ,idir
     &           ,ihilo
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T x1
      REAL_T dxCrse
      integer idir
      integer ihilo
      REAL_T x2, denom, idenom, x, xsquared, m1, m2
      REAL_T q1, q2
      REAL_T pa, pb, a, b
      INTEGER ncomp,  n
      INTEGER ii0
      INTEGER i0
      x2 = half*(three*x1+dxCrse)
      denom = one-((x1+x2)/x1)
      idenom = one/(denom)
      x = two*x1
      xsquared = x*x
      m1 = one/(x1*x1)
      m2 = one/(x1*(x1-x2))
      q1 = one/(x1-x2)
      q2 = x1+x2
      ihilo = ihilo*(-1)
      ncomp = nphicomp
      
      ii0= ihilo*CHF_ID(0, idir)

      do n = 0, ncomp-1
          
      do i0 = iregionlo0,iregionhi0

          pa=phi(i0+2*ii0,n)
          pb=phi(i0+ii0,n)
          a=((pb-pa)*m1 - (pb)*m2)*idenom
          b=(pb)*q1 - a*q2
          phi(i0,n) = a*xsquared + b*x + pa
          
      enddo
      enddo
      return
      end
      subroutine INTERPHOMO(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,x1
     &           ,dxCrse
     &           ,idir
     &           ,ihilo
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T x1
      REAL_T dxCrse
      integer idir
      integer ihilo
      REAL_T c1, c2
      REAL_T pa, pb
      INTEGER ncomp,  n
      INTEGER ii0
      INTEGER i0
      c1=two*(dxCrse-x1)/(dxCrse+x1)
      c2=   -(dxCrse-x1)/(dxCrse+three*x1)
      ihilo = ihilo*(-1)
      ncomp = nphicomp
      
      ii0= ihilo*CHF_ID(0, idir)

      do n = 0, ncomp-1
          
      do i0 = iregionlo0,iregionhi0

          pa=phi(i0+ii0,n)
          pb=phi(i0+2*ii0,n)
          phi(i0,n) = c1*pa + c2*pb
          
      enddo
      enddo
      return
      end
      subroutine INTERPHOMOLINEAR(
     &           phi
     &           ,iphilo0
     &           ,iphihi0
     &           ,nphicomp
     &           ,iregionlo0
     &           ,iregionhi0
     &           ,x1
     &           ,dxCrse
     &           ,idir
     &           ,ihilo
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           0:nphicomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL_T x1
      REAL_T dxCrse
      integer idir
      integer ihilo
      INTEGER ncomp,  n
      INTEGER ii0
      INTEGER i0
      REAL_T pa, factor
      ihilo = ihilo*(-1)
      ncomp = nphicomp
      
      ii0= ihilo*CHF_ID(0, idir)

      factor = one - two*x1/(x1+dxCrse)
      do n = 0, ncomp-1
          
      do i0 = iregionlo0,iregionhi0

          pa=phi(i0+ii0,n)
          phi(i0,n) = factor*pa
          
      enddo
      enddo
      return
      end
