#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
#define SMALLNUMBER 1.0d-9
      subroutine AVEFACESCALTOFACEVECT(
     &           facevect
     &           ,ifacevectlo0,ifacevectlo1,ifacevectlo2
     &           ,ifacevecthi0,ifacevecthi1,ifacevecthi2
     &           ,nfacevectcomp
     &           ,facescal
     &           ,ifacescallo0,ifacescallo1,ifacescallo2
     &           ,ifacescalhi0,ifacescalhi1,ifacescalhi2
     &           ,facedir
     &           ,vectdir
     &           ,idcalcfacelo0,idcalcfacelo1,idcalcfacelo2
     &           ,idcalcfacehi0,idcalcfacehi1,idcalcfacehi2
     &           ,ioffboxlo0,ioffboxlo1,ioffboxlo2
     &           ,ioffboxhi0,ioffboxhi1,ioffboxhi2
     &           )

      implicit none
      integer nfacevectcomp
      integer ifacevectlo0,ifacevectlo1,ifacevectlo2
      integer ifacevecthi0,ifacevecthi1,ifacevecthi2
      REAL_T facevect(
     &           ifacevectlo0:ifacevecthi0,
     &           ifacevectlo1:ifacevecthi1,
     &           ifacevectlo2:ifacevecthi2,
     &           0:nfacevectcomp-1)
      integer ifacescallo0,ifacescallo1,ifacescallo2
      integer ifacescalhi0,ifacescalhi1,ifacescalhi2
      REAL_T facescal(
     &           ifacescallo0:ifacescalhi0,
     &           ifacescallo1:ifacescalhi1,
     &           ifacescallo2:ifacescalhi2)
      integer facedir
      integer vectdir
      integer idcalcfacelo0,idcalcfacelo1,idcalcfacelo2
      integer idcalcfacehi0,idcalcfacehi1,idcalcfacehi2
      integer ioffboxlo0,ioffboxlo1,ioffboxlo2
      integer ioffboxhi0,ioffboxhi1,ioffboxhi2
      integer i0,i1,i2
      integer ioff0,ioff1,ioff2
      integer numpts
      real_t wt, tot
      if (facedir .eq. vectdir) then
         
      do i2 = idcalcfacelo2,idcalcfacehi2
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0

            facevect(i0,i1,i2, vectdir) = facescal(i0,i1,i2)
         
      enddo
      enddo
      enddo
      else
         numpts = 2**CH_SPACEDIM
         wt = one / (numpts * one)
         
      do i2 = idcalcfacelo2,idcalcfacehi2
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0

            tot = zero
            
      do ioff2 = ioffboxlo2,ioffboxhi2
      do ioff1 = ioffboxlo1,ioffboxhi1
      do ioff0 = ioffboxlo0,ioffboxhi0

               tot = tot + facescal(i0 +ioff0,i1 +ioff1,i2 +ioff2)
            
      enddo
      enddo
      enddo
            facevect(i0,i1,i2, vectdir) = wt * tot
         
      enddo
      enddo
      enddo
      endif
      return
      end
      subroutine AVESCALTOFACE(
     &           facescal
     &           ,ifacescallo0,ifacescallo1,ifacescallo2
     &           ,ifacescalhi0,ifacescalhi1,ifacescalhi2
     &           ,cellscal
     &           ,icellscallo0,icellscallo1,icellscallo2
     &           ,icellscalhi0,icellscalhi1,icellscalhi2
     &           ,idir
     &           ,idcalcfacelo0,idcalcfacelo1,idcalcfacelo2
     &           ,idcalcfacehi0,idcalcfacehi1,idcalcfacehi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ifacescallo0,ifacescallo1,ifacescallo2
      integer ifacescalhi0,ifacescalhi1,ifacescalhi2
      REAL_T facescal(
     &           ifacescallo0:ifacescalhi0,
     &           ifacescallo1:ifacescalhi1,
     &           ifacescallo2:ifacescalhi2)
      integer icellscallo0,icellscallo1,icellscallo2
      integer icellscalhi0,icellscalhi1,icellscalhi2
      REAL_T cellscal(
     &           icellscallo0:icellscalhi0,
     &           icellscallo1:icellscalhi1,
     &           icellscallo2:icellscalhi2)
      integer idir
      integer idcalcfacelo0,idcalcfacelo1,idcalcfacelo2
      integer idcalcfacehi0,idcalcfacehi1,idcalcfacehi2
      integer i0,i1,i2
      integer f2cLo0,f2cLo1,f2cLo2
      integer f2cHi0,f2cHi1,f2cHi2
      
      f2cLo0= -1*CHF_ID(0, idir)

      f2cLo1= -1*CHF_ID(1, idir)

      f2cLo2= -1*CHF_ID(2, idir)

      
      f2cHi0= 0*CHF_ID(0, idir)

      f2cHi1= 0*CHF_ID(1, idir)

      f2cHi2= 0*CHF_ID(2, idir)

      
      do i2 = idcalcfacelo2,idcalcfacehi2
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0

         facescal(i0,i1,i2) =  half * (
     &     cellscal(i0 +f2cLo0,i1 +f2cLo1,i2 +f2cLo2) +
     &     cellscal(i0 +f2cHi0,i1 +f2cHi1,i2 +f2cHi2) )
      
      enddo
      enddo
      enddo
      return
      end
      subroutine AVECELLTOFACE(
     &           facevel
     &           ,ifacevello0,ifacevello1,ifacevello2
     &           ,ifacevelhi0,ifacevelhi1,ifacevelhi2
     &           ,cellvel
     &           ,icellvello0,icellvello1,icellvello2
     &           ,icellvelhi0,icellvelhi1,icellvelhi2
     &           ,ncellvelcomp
     &           ,idir
     &           ,idcalcfacelo0,idcalcfacelo1,idcalcfacelo2
     &           ,idcalcfacehi0,idcalcfacehi1,idcalcfacehi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ifacevello0,ifacevello1,ifacevello2
      integer ifacevelhi0,ifacevelhi1,ifacevelhi2
      REAL_T facevel(
     &           ifacevello0:ifacevelhi0,
     &           ifacevello1:ifacevelhi1,
     &           ifacevello2:ifacevelhi2)
      integer ncellvelcomp
      integer icellvello0,icellvello1,icellvello2
      integer icellvelhi0,icellvelhi1,icellvelhi2
      REAL_T cellvel(
     &           icellvello0:icellvelhi0,
     &           icellvello1:icellvelhi1,
     &           icellvello2:icellvelhi2,
     &           0:ncellvelcomp-1)
      integer idir
      integer idcalcfacelo0,idcalcfacelo1,idcalcfacelo2
      integer idcalcfacehi0,idcalcfacehi1,idcalcfacehi2
      integer i0,i1,i2
      integer f2cLo0,f2cLo1,f2cLo2
      integer f2cHi0,f2cHi1,f2cHi2
      
      f2cLo0= -1*CHF_ID(0, idir)

      f2cLo1= -1*CHF_ID(1, idir)

      f2cLo2= -1*CHF_ID(2, idir)

      
      f2cHi0= 0*CHF_ID(0, idir)

      f2cHi1= 0*CHF_ID(1, idir)

      f2cHi2= 0*CHF_ID(2, idir)

      
      do i2 = idcalcfacelo2,idcalcfacehi2
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0

         facevel(i0,i1,i2) =  half * (
     &     cellvel(i0 +f2cLo0,i1 +f2cLo1,i2 +f2cLo2, idir) +
     &     cellvel(i0 +f2cHi0,i1 +f2cHi1,i2 +f2cHi2, idir) )
      
      enddo
      enddo
      enddo
      return
      end
      subroutine AVEFACETOCELL(
     &           cellvel
     &           ,icellvello0,icellvello1,icellvello2
     &           ,icellvelhi0,icellvelhi1,icellvelhi2
     &           ,ncellvelcomp
     &           ,facevel
     &           ,ifacevello0,ifacevello1,ifacevello2
     &           ,ifacevelhi0,ifacevelhi1,ifacevelhi2
     &           ,idir
     &           ,idcalccelllo0,idcalccelllo1,idcalccelllo2
     &           ,idcalccellhi0,idcalccellhi1,idcalccellhi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ncellvelcomp
      integer icellvello0,icellvello1,icellvello2
      integer icellvelhi0,icellvelhi1,icellvelhi2
      REAL_T cellvel(
     &           icellvello0:icellvelhi0,
     &           icellvello1:icellvelhi1,
     &           icellvello2:icellvelhi2,
     &           0:ncellvelcomp-1)
      integer ifacevello0,ifacevello1,ifacevello2
      integer ifacevelhi0,ifacevelhi1,ifacevelhi2
      REAL_T facevel(
     &           ifacevello0:ifacevelhi0,
     &           ifacevello1:ifacevelhi1,
     &           ifacevello2:ifacevelhi2)
      integer idir
      integer idcalccelllo0,idcalccelllo1,idcalccelllo2
      integer idcalccellhi0,idcalccellhi1,idcalccellhi2
      integer i0,i1,i2
      integer c2fLo0,c2fLo1,c2fLo2
      integer c2fHi0,c2fHi1,c2fHi2
      
      c2fLo0= 0*CHF_ID(0, idir)

      c2fLo1= 0*CHF_ID(1, idir)

      c2fLo2= 0*CHF_ID(2, idir)

      
      c2fHi0= 1*CHF_ID(0, idir)

      c2fHi1= 1*CHF_ID(1, idir)

      c2fHi2= 1*CHF_ID(2, idir)

      
      do i2 = idcalccelllo2,idcalccellhi2
      do i1 = idcalccelllo1,idcalccellhi1
      do i0 = idcalccelllo0,idcalccellhi0

         cellvel(i0,i1,i2, idir) =  half * (
     &     facevel(i0 +c2fLo0,i1 +c2fLo1,i2 +c2fLo2) +
     &     facevel(i0 +c2fHi0,i1 +c2fHi1,i2 +c2fHi2) )
      
      enddo
      enddo
      enddo
      return
      end
      subroutine MAGNITUDEF(
     &           magdata
     &           ,imagdatalo0,imagdatalo1,imagdatalo2
     &           ,imagdatahi0,imagdatahi1,imagdatahi2
     &           ,data
     &           ,idatalo0,idatalo1,idatalo2
     &           ,idatahi0,idatahi1,idatahi2
     &           ,ndatacomp
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer imagdatalo0,imagdatalo1,imagdatalo2
      integer imagdatahi0,imagdatahi1,imagdatahi2
      REAL_T magdata(
     &           imagdatalo0:imagdatahi0,
     &           imagdatalo1:imagdatahi1,
     &           imagdatalo2:imagdatahi2)
      integer ndatacomp
      integer idatalo0,idatalo1,idatalo2
      integer idatahi0,idatahi1,idatahi2
      REAL_T data(
     &           idatalo0:idatahi0,
     &           idatalo1:idatahi1,
     &           idatalo2:idatahi2,
     &           0:ndatacomp-1)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i0,i1,i2
      integer iv
      real_t cur,sum
      
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

         sum = zero
         do iv = 0, ndatacomp-1
            cur = data(i0,i1,i2, iv)
            sum = sum + cur*cur
         enddo
         magdata(i0,i1,i2) = sqrt(sum)
      
      enddo
      enddo
      enddo
      return
      end
        subroutine GETRELGRADF(
     &           du
     &           ,idulo0,idulo1,idulo2
     &           ,iduhi0,iduhi1,iduhi2
     &           ,u
     &           ,iulo0,iulo1,iulo2
     &           ,iuhi0,iuhi1,iuhi2
     &           ,idir
     &           ,iloBoxlo0,iloBoxlo1,iloBoxlo2
     &           ,iloBoxhi0,iloBoxhi1,iloBoxhi2
     &           ,hasLo
     &           ,ihiBoxlo0,ihiBoxlo1,ihiBoxlo2
     &           ,ihiBoxhi0,ihiBoxhi1,ihiBoxhi2
     &           ,hasHi
     &           ,icenterBoxlo0,icenterBoxlo1,icenterBoxlo2
     &           ,icenterBoxhi0,icenterBoxhi1,icenterBoxhi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer idulo0,idulo1,idulo2
      integer iduhi0,iduhi1,iduhi2
      REAL_T du(
     &           idulo0:iduhi0,
     &           idulo1:iduhi1,
     &           idulo2:iduhi2)
      integer iulo0,iulo1,iulo2
      integer iuhi0,iuhi1,iuhi2
      REAL_T u(
     &           iulo0:iuhi0,
     &           iulo1:iuhi1,
     &           iulo2:iuhi2)
      integer idir
      integer iloBoxlo0,iloBoxlo1,iloBoxlo2
      integer iloBoxhi0,iloBoxhi1,iloBoxhi2
      integer hasLo
      integer ihiBoxlo0,ihiBoxlo1,ihiBoxlo2
      integer ihiBoxhi0,ihiBoxhi1,ihiBoxhi2
      integer hasHi
      integer icenterBoxlo0,icenterBoxlo1,icenterBoxlo2
      integer icenterBoxhi0,icenterBoxhi1,icenterBoxhi2
        integer i0,i1,i2
        integer ioff0,ioff1,ioff2
        real_t diff,aver
        
      ioff0=CHF_ID(0, idir)

      ioff1=CHF_ID(1, idir)

      ioff2=CHF_ID(2, idir)

        
      do i2 = icenterBoxlo2,icenterBoxhi2
      do i1 = icenterBoxlo1,icenterBoxhi1
      do i0 = icenterBoxlo0,icenterBoxhi0

          diff = half*(u(i0 +ioff0,i1 +ioff1,i2 +ioff2)
     &                -u(i0 -ioff0,i1 -ioff1,i2 -ioff2))
          aver = half*(u(i0 +ioff0,i1 +ioff1,i2 +ioff2)
     &                +u(i0 -ioff0,i1 -ioff1,i2 -ioff2))
          du(i0,i1,i2) = diff / aver
        
      enddo
      enddo
      enddo
        if (hasLo .eq. 1) then
          
      do i2 = iloBoxlo2,iloBoxhi2
      do i1 = iloBoxlo1,iloBoxhi1
      do i0 = iloBoxlo0,iloBoxhi0

            diff = u(i0 +ioff0,i1 +ioff1,i2 +ioff2) - u(i0,i1,i2)
            aver = half*(u(i0 +ioff0,i1 +ioff1,i2 +ioff2) + u(i0,i1,i2))
            du(i0,i1,i2) = diff / aver
          
      enddo
      enddo
      enddo
        endif
        if (hasHi .eq. 1) then
          
      do i2 = ihiBoxlo2,ihiBoxhi2
      do i1 = ihiBoxlo1,ihiBoxhi1
      do i0 = ihiBoxlo0,ihiBoxhi0

            diff = u(i0,i1,i2) - u(i0 -ioff0,i1 -ioff1,i2 -ioff2)
            aver = half*(u(i0,i1,i2) + u(i0 -ioff0,i1 -ioff1,i2 -ioff2))
            du(i0,i1,i2) = diff / aver
          
      enddo
      enddo
      enddo
        endif
        return
        end
      subroutine POSTNORMALSOURCE(
     &           dWminus
     &           ,idWminuslo0,idWminuslo1,idWminuslo2
     &           ,idWminushi0,idWminushi1,idWminushi2
     &           ,ndWminuscomp
     &           ,dWplus
     &           ,idWpluslo0,idWpluslo1,idWpluslo2
     &           ,idWplushi0,idWplushi1,idWplushi2
     &           ,ndWpluscomp
     &           ,W
     &           ,iWlo0,iWlo1,iWlo2
     &           ,iWhi0,iWhi1,iWhi2
     &           ,nWcomp
     &           ,advVel
     &           ,iadvVello0,iadvVello1,iadvVello2
     &           ,iadvVelhi0,iadvVelhi1,iadvVelhi2
     &           ,dt
     &           ,dx
     &           ,idir
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ndWminuscomp
      integer idWminuslo0,idWminuslo1,idWminuslo2
      integer idWminushi0,idWminushi1,idWminushi2
      REAL_T dWminus(
     &           idWminuslo0:idWminushi0,
     &           idWminuslo1:idWminushi1,
     &           idWminuslo2:idWminushi2,
     &           0:ndWminuscomp-1)
      integer ndWpluscomp
      integer idWpluslo0,idWpluslo1,idWpluslo2
      integer idWplushi0,idWplushi1,idWplushi2
      REAL_T dWplus(
     &           idWpluslo0:idWplushi0,
     &           idWpluslo1:idWplushi1,
     &           idWpluslo2:idWplushi2,
     &           0:ndWpluscomp-1)
      integer nWcomp
      integer iWlo0,iWlo1,iWlo2
      integer iWhi0,iWhi1,iWhi2
      REAL_T W(
     &           iWlo0:iWhi0,
     &           iWlo1:iWhi1,
     &           iWlo2:iWhi2,
     &           0:nWcomp-1)
      integer iadvVello0,iadvVello1,iadvVello2
      integer iadvVelhi0,iadvVelhi1,iadvVelhi2
      REAL_T advVel(
     &           iadvVello0:iadvVelhi0,
     &           iadvVello1:iadvVelhi1,
     &           iadvVello2:iadvVelhi2)
      REAL_T dt
      REAL_T dx
      integer idir
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i0,i1,i2
      integer c2fLo0,c2fLo1,c2fLo2
      integer c2fHi0,c2fHi1,c2fHi2
      integer n
      real_t dudx
      
      c2fLo0= 0*CHF_ID(0, idir)

      c2fLo1= 0*CHF_ID(1, idir)

      c2fLo2= 0*CHF_ID(2, idir)

      
      c2fHi0= 1*CHF_ID(0, idir)

      c2fHi1= 1*CHF_ID(1, idir)

      c2fHi2= 1*CHF_ID(2, idir)

      
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

         dudx = advVel(i0 +c2fHi0,i1 +c2fHi1,i2 +c2fHi2)
     &     - advVel(i0 +c2fLo0,i1 +c2fLo1,i2 +c2fLo2)
         dudx = dt*dudx / dx * half
         do n=0, nWcomp-1
            dWplus(i0,i1,i2, n) =
     &           dWPlus(i0,i1,i2, n) - dudx * W(i0,i1,i2,n)
            dWMinus(i0,i1,i2, n) =
     &           dWMinus(i0,i1,i2, n) - dudx * W(i0,i1,i2,n)
         end do
      
      enddo
      enddo
      enddo
      return
      end
      subroutine RIEMANNF(
     &           Wgdnv
     &           ,iWgdnvlo0,iWgdnvlo1,iWgdnvlo2
     &           ,iWgdnvhi0,iWgdnvhi1,iWgdnvhi2
     &           ,nWgdnvcomp
     &           ,WLeft
     &           ,iWLeftlo0,iWLeftlo1,iWLeftlo2
     &           ,iWLefthi0,iWLefthi1,iWLefthi2
     &           ,nWLeftcomp
     &           ,WRight
     &           ,iWRightlo0,iWRightlo1,iWRightlo2
     &           ,iWRighthi0,iWRighthi1,iWRighthi2
     &           ,nWRightcomp
     &           ,advVel
     &           ,iadvVello0,iadvVello1,iadvVello2
     &           ,iadvVelhi0,iadvVelhi1,iadvVelhi2
     &           ,idir
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer nWgdnvcomp
      integer iWgdnvlo0,iWgdnvlo1,iWgdnvlo2
      integer iWgdnvhi0,iWgdnvhi1,iWgdnvhi2
      REAL_T Wgdnv(
     &           iWgdnvlo0:iWgdnvhi0,
     &           iWgdnvlo1:iWgdnvhi1,
     &           iWgdnvlo2:iWgdnvhi2,
     &           0:nWgdnvcomp-1)
      integer nWLeftcomp
      integer iWLeftlo0,iWLeftlo1,iWLeftlo2
      integer iWLefthi0,iWLefthi1,iWLefthi2
      REAL_T WLeft(
     &           iWLeftlo0:iWLefthi0,
     &           iWLeftlo1:iWLefthi1,
     &           iWLeftlo2:iWLefthi2,
     &           0:nWLeftcomp-1)
      integer nWRightcomp
      integer iWRightlo0,iWRightlo1,iWRightlo2
      integer iWRighthi0,iWRighthi1,iWRighthi2
      REAL_T WRight(
     &           iWRightlo0:iWRighthi0,
     &           iWRightlo1:iWRighthi1,
     &           iWRightlo2:iWRighthi2,
     &           0:nWRightcomp-1)
      integer iadvVello0,iadvVello1,iadvVello2
      integer iadvVelhi0,iadvVelhi1,iadvVelhi2
      REAL_T advVel(
     &           iadvVello0:iadvVelhi0,
     &           iadvVello1:iadvVelhi1,
     &           iadvVello2:iadvVelhi2)
      integer idir
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i0,i1,i2
      integer n
      real_t sl,sr
      real_t so
      real_t ustar
      do n=0, nWLeftcomp-1
         
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

            sl =  WLeft(i0,i1,i2, n)
            sr = WRight(i0,i1,i2, n)
            ustar = advVel(i0,i1,i2)
            if (ustar .gt. zero) then
               so = sl
            else
               so = sr
            endif
            if (abs(ustar) .lt. SMALLNUMBER) then
               so = half*(sl+sr)
            endif
            Wgdnv(i0,i1,i2, n) = so
         
      enddo
      enddo
      enddo
      end do
      return
      end
      subroutine QUASILINEARUPDATE(
     &           AdWdx
     &           ,iAdWdxlo0,iAdWdxlo1,iAdWdxlo2
     &           ,iAdWdxhi0,iAdWdxhi1,iAdWdxhi2
     &           ,nAdWdxcomp
     &           ,WHalf
     &           ,iWHalflo0,iWHalflo1,iWHalflo2
     &           ,iWHalfhi0,iWHalfhi1,iWHalfhi2
     &           ,nWHalfcomp
     &           ,cellVel
     &           ,icellVello0,icellVello1,icellVello2
     &           ,icellVelhi0,icellVelhi1,icellVelhi2
     &           ,scale
     &           ,idir
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nAdWdxcomp
      integer iAdWdxlo0,iAdWdxlo1,iAdWdxlo2
      integer iAdWdxhi0,iAdWdxhi1,iAdWdxhi2
      REAL_T AdWdx(
     &           iAdWdxlo0:iAdWdxhi0,
     &           iAdWdxlo1:iAdWdxhi1,
     &           iAdWdxlo2:iAdWdxhi2,
     &           0:nAdWdxcomp-1)
      integer nWHalfcomp
      integer iWHalflo0,iWHalflo1,iWHalflo2
      integer iWHalfhi0,iWHalfhi1,iWHalfhi2
      REAL_T WHalf(
     &           iWHalflo0:iWHalfhi0,
     &           iWHalflo1:iWHalfhi1,
     &           iWHalflo2:iWHalfhi2,
     &           0:nWHalfcomp-1)
      integer icellVello0,icellVello1,icellVello2
      integer icellVelhi0,icellVelhi1,icellVelhi2
      REAL_T cellVel(
     &           icellVello0:icellVelhi0,
     &           icellVello1:icellVelhi1,
     &           icellVello2:icellVelhi2)
      REAL_T scale
      integer idir
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i0,i1,i2
      integer c2fLo0,c2fLo1,c2fLo2
      integer c2fHi0,c2fHi1,c2fHi2
      integer n
      
      c2fLo0= 0*CHF_ID(0, idir)

      c2fLo1= 0*CHF_ID(1, idir)

      c2fLo2= 0*CHF_ID(2, idir)

      
      c2fHi0= 1*CHF_ID(0, idir)

      c2fHi1= 1*CHF_ID(1, idir)

      c2fHi2= 1*CHF_ID(2, idir)

      do n=0, nAdWdxcomp-1
         
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

            AdWdx(i0,i1,i2, n) =
     &        scale * cellVel(i0,i1,i2) *
     &        (WHalf(i0 +c2fHi0,i1 +c2fHi1,i2 +c2fHi2, n) -
     &         WHalf(i0 +c2fLo0,i1 +c2fLo1,i2 +c2fLo2, n))
         
      enddo
      enddo
      enddo
      end do
      return
      end
      subroutine AVECELLVECTOFACEVEC(
     &           facevec
     &           ,ifaceveclo0,ifaceveclo1,ifaceveclo2
     &           ,ifacevechi0,ifacevechi1,ifacevechi2
     &           ,nfaceveccomp
     &           ,cellvec
     &           ,icellveclo0,icellveclo1,icellveclo2
     &           ,icellvechi0,icellvechi1,icellvechi2
     &           ,ncellveccomp
     &           ,facedir
     &           ,idcalcfacelo0,idcalcfacelo1,idcalcfacelo2
     &           ,idcalcfacehi0,idcalcfacehi1,idcalcfacehi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfaceveccomp
      integer ifaceveclo0,ifaceveclo1,ifaceveclo2
      integer ifacevechi0,ifacevechi1,ifacevechi2
      REAL_T facevec(
     &           ifaceveclo0:ifacevechi0,
     &           ifaceveclo1:ifacevechi1,
     &           ifaceveclo2:ifacevechi2,
     &           0:nfaceveccomp-1)
      integer ncellveccomp
      integer icellveclo0,icellveclo1,icellveclo2
      integer icellvechi0,icellvechi1,icellvechi2
      REAL_T cellvec(
     &           icellveclo0:icellvechi0,
     &           icellveclo1:icellvechi1,
     &           icellveclo2:icellvechi2,
     &           0:ncellveccomp-1)
      integer facedir
      integer idcalcfacelo0,idcalcfacelo1,idcalcfacelo2
      integer idcalcfacehi0,idcalcfacehi1,idcalcfacehi2
      integer i0,i1,i2
      integer f2cLo0,f2cLo1,f2cLo2
      integer f2cHi0,f2cHi1,f2cHi2
      integer idir
      
      f2cLo0= -1*CHF_ID(0, facedir)

      f2cLo1= -1*CHF_ID(1, facedir)

      f2cLo2= -1*CHF_ID(2, facedir)

      
      f2cHi0= 0*CHF_ID(0, facedir)

      f2cHi1= 0*CHF_ID(1, facedir)

      f2cHi2= 0*CHF_ID(2, facedir)

      do idir = 0, CH_SPACEDIM-1
         
      do i2 = idcalcfacelo2,idcalcfacehi2
      do i1 = idcalcfacelo1,idcalcfacehi1
      do i0 = idcalcfacelo0,idcalcfacehi0

            facevec(i0,i1,i2, idir) =  half * (
     &           cellvec(i0 +f2cLo0,i1 +f2cLo1,i2 +f2cLo2, idir) +
     &           cellvec(i0 +f2cHi0,i1 +f2cHi1,i2 +f2cHi2, idir) )
         
      enddo
      enddo
      enddo
      enddo
      return
      end
