#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine PROLONGROP(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,m
     &           )

      implicit none
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           icoarselo2:icoarsehi2,
     &           0:ncoarsecomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer m
      INTEGER ncomp, n
      integer i,j,k
      integer ii,jj,kk
      ncomp = nphicomp
      do n = 0, ncomp-1
          
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          
          ii = (i-abs(mod(i,m)))/m
          jj = (j-abs(mod(j,m)))/m
          kk = (k-abs(mod(k,m)))/m
          phi(i,j,k,n) =  phi(i,j,k,n) +
     &        coarse(ii,jj,kk,n)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine RESTRICTRESROP(
     &           res
     &           ,ireslo0,ireslo1,ireslo2
     &           ,ireshi0,ireshi1,ireshi2
     &           ,nrescomp
     &           ,resFine
     &           ,iresFinelo0,iresFinelo1,iresFinelo2
     &           ,iresFinehi0,iresFinehi1,iresFinehi2
     &           ,nresFinecomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           )

      implicit none
      integer nrescomp
      integer ireslo0,ireslo1,ireslo2
      integer ireshi0,ireshi1,ireshi2
      REAL_T res(
     &           ireslo0:ireshi0,
     &           ireslo1:ireshi1,
     &           ireslo2:ireshi2,
     &           0:nrescomp-1)
      integer nresFinecomp
      integer iresFinelo0,iresFinelo1,iresFinelo2
      integer iresFinehi0,iresFinehi1,iresFinehi2
      REAL_T resFine(
     &           iresFinelo0:iresFinehi0,
     &           iresFinelo1:iresFinehi1,
     &           iresFinelo2:iresFinehi2,
     &           0:nresFinecomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T denom
      integer n
      integer i,j,k
      integer ii,jj,kk
      integer ncomp
      ncomp = nrescomp
      denom = two *two *two
      do n = 0, ncomp-1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         
         ii = (i-abs(mod(i,2)))/2
         jj = (j-abs(mod(j,2)))/2
         kk = (k-abs(mod(k,2)))/2
         res(ii,jj,kk,n) = res(ii,jj,kk,n) + resFine(i,j,k,n)/denom
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine CELLGRADROP(
     &           grad
     &           ,igradlo0,igradlo1,igradlo2
     &           ,igradhi0,igradhi1,igradhi2
     &           ,vel
     &           ,ivello0,ivello1,ivello2
     &           ,ivelhi0,ivelhi1,ivelhi2
     &           ,igridlo0,igridlo1,igridlo2
     &           ,igridhi0,igridhi1,igridhi2
     &           ,dx
     &           ,divdir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradlo0,igradlo1,igradlo2
      integer igradhi0,igradhi1,igradhi2
      REAL_T grad(
     &           igradlo0:igradhi0,
     &           igradlo1:igradhi1,
     &           igradlo2:igradhi2)
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1,
     &           ivello2:ivelhi2)
      integer igridlo0,igridlo1,igridlo2
      integer igridhi0,igridhi1,igridhi2
      REAL_T dx
      integer divdir
      integer ii,i,jj,j,kk,k
      
      ii = chf_id(divdir, 0)
      jj = chf_id(divdir, 1)
      kk = chf_id(divdir, 2)
      
      do k = igridlo2,igridhi2
      do j = igridlo1,igridhi1
      do i = igridlo0,igridhi0

      grad(i,j,k) =
     $     (    vel(i+ii,j+jj,k+kk)
     $     -    vel(i-ii,j-jj,k-kk) )/(two*dx)
      
      enddo
      enddo
      enddo
      return
      end
      subroutine ADDGRADTOFLUXROP(
     &           flux
     &           ,ifluxlo0,ifluxlo1,ifluxlo2
     &           ,ifluxhi0,ifluxhi1,ifluxhi2
     &           ,nfluxcomp
     &           ,fluxComp
     &           ,grad
     &           ,igradlo0,igradlo1,igradlo2
     &           ,igradhi0,igradhi1,igradhi2
     &           ,ngradcomp
     &           ,gradComp
     &           ,sign
     &           ,iregionfacelo0,iregionfacelo1,iregionfacelo2
     &           ,iregionfacehi0,iregionfacehi1,iregionfacehi2
     &           )

      implicit none
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           ifluxlo1:ifluxhi1,
     &           ifluxlo2:ifluxhi2,
     &           0:nfluxcomp-1)
      integer fluxComp
      integer ngradcomp
      integer igradlo0,igradlo1,igradlo2
      integer igradhi0,igradhi1,igradhi2
      REAL_T grad(
     &           igradlo0:igradhi0,
     &           igradlo1:igradhi1,
     &           igradlo2:igradhi2,
     &           0:ngradcomp-1)
      integer gradComp
      REAL_T sign
      integer iregionfacelo0,iregionfacelo1,iregionfacelo2
      integer iregionfacehi0,iregionfacehi1,iregionfacehi2
      integer i,j,k
      
      do k = iregionfacelo2,iregionfacehi2
      do j = iregionfacelo1,iregionfacehi1
      do i = iregionfacelo0,iregionfacehi0

      flux(i,j,k, fluxcomp) =  flux(i,j,k, fluxcomp) +
     $     ( sign*grad(i,j,k, gradcomp))
      
      enddo
      enddo
      enddo
      return
      end
      subroutine GETFACEGRADROP(
     &           gradvelface
     &           ,igradvelfacelo0,igradvelfacelo1,igradvelfacelo2
     &           ,igradvelfacehi0,igradvelfacehi1,igradvelfacehi2
     &           ,gradvelcell
     &           ,igradvelcelllo0,igradvelcelllo1,igradvelcelllo2
     &           ,igradvelcellhi0,igradvelcellhi1,igradvelcellhi2
     &           ,velcomp
     &           ,ivelcomplo0,ivelcomplo1,ivelcomplo2
     &           ,ivelcomphi0,ivelcomphi1,ivelcomphi2
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     &           ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     &           ,iloboxlo0,iloboxlo1,iloboxlo2
     &           ,iloboxhi0,iloboxhi1,iloboxhi2
     &           ,haslo
     &           ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     &           ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     &           ,hashi
     &           ,dx
     &           ,facedir
     &           ,divdir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer igradvelfacelo0,igradvelfacelo1,igradvelfacelo2
      integer igradvelfacehi0,igradvelfacehi1,igradvelfacehi2
      REAL_T gradvelface(
     &           igradvelfacelo0:igradvelfacehi0,
     &           igradvelfacelo1:igradvelfacehi1,
     &           igradvelfacelo2:igradvelfacehi2)
      integer igradvelcelllo0,igradvelcelllo1,igradvelcelllo2
      integer igradvelcellhi0,igradvelcellhi1,igradvelcellhi2
      REAL_T gradvelcell(
     &           igradvelcelllo0:igradvelcellhi0,
     &           igradvelcelllo1:igradvelcellhi1,
     &           igradvelcelllo2:igradvelcellhi2)
      integer ivelcomplo0,ivelcomplo1,ivelcomplo2
      integer ivelcomphi0,ivelcomphi1,ivelcomphi2
      REAL_T velcomp(
     &           ivelcomplo0:ivelcomphi0,
     &           ivelcomplo1:ivelcomphi1,
     &           ivelcomplo2:ivelcomphi2)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      REAL_T dx
      integer facedir
      integer divdir
      integer ii,i,jj,j,kk,k
      
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      kk = chf_id(facedir, 2)
      if (facedir .eq. divdir) then
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         gradvelface(i,j,k) =
     $        ( velcomp(i   ,j   ,k   )
     $        - velcomp(i-ii,j-jj,k-kk) )/dx
         
      enddo
      enddo
      enddo
      else
         
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

         gradvelface(i,j,k) =
     $        ( gradvelcell(i   ,j   ,k   )
     $        + gradvelcell(i-ii,j-jj,k-kk) )/two
         
      enddo
      enddo
      enddo
         if(haslo .eq. 1) then
            
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

            gradvelface(i,j,k) =
     $           (three*gradvelcell(i   ,j   ,k   )
     $           -      gradvelcell(i+ii,j+jj,k+kk))/two
            
      enddo
      enddo
      enddo
         endif
         if(hashi .eq. 1) then
            
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

            gradvelface(i,j,k) =
     $           (three*gradvelcell(i-  ii,j-  jj,k-  kk)
     $           -      gradvelcell(i-2*ii,j-2*jj,k-2*kk))/two
            
      enddo
      enddo
      enddo
         endif
      endif
      return
      end
      subroutine CELLDIVINCRROP(
     &           divvel
     &           ,idivvello0,idivvello1,idivvello2
     &           ,idivvelhi0,idivvelhi1,idivvelhi2
     &           ,vel
     &           ,ivello0,ivello1,ivello2
     &           ,ivelhi0,ivelhi1,ivelhi2
     &           ,nvelcomp
     &           ,dx
     &           ,divdir
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer idivvello0,idivvello1,idivvello2
      integer idivvelhi0,idivvelhi1,idivvelhi2
      REAL_T divvel(
     &           idivvello0:idivvelhi0,
     &           idivvello1:idivvelhi1,
     &           idivvello2:idivvelhi2)
      integer nvelcomp
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1,
     &           ivello2:ivelhi2,
     &           0:nvelcomp-1)
      REAL_T dx
      integer divdir
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer ii,i,jj,j,kk,k
      
      ii = chf_id(divdir, 0)
      jj = chf_id(divdir, 1)
      kk = chf_id(divdir, 2)
      
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

      divvel(i,j,k) = divvel(i,j,k) +
     $     (    vel(i+ii,j+jj,k+kk,divdir)
     $     -    vel(i-ii,j-jj,k-kk,divdir) )/(two*dx)
      
      enddo
      enddo
      enddo
      return
      end
      subroutine FACEDIVINCRROP(
     &           divvel
     &           ,idivvello0,idivvello1,idivvello2
     &           ,idivvelhi0,idivvelhi1,idivvelhi2
     &           ,vel
     &           ,ivello0,ivello1,ivello2
     &           ,ivelhi0,ivelhi1,ivelhi2
     &           ,nvelcomp
     &           ,gradvel
     &           ,igradvello0,igradvello1,igradvello2
     &           ,igradvelhi0,igradvelhi1,igradvelhi2
     &           ,ngradvelcomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,icenterboxlo0,icenterboxlo1,icenterboxlo2
     &           ,icenterboxhi0,icenterboxhi1,icenterboxhi2
     &           ,iloboxlo0,iloboxlo1,iloboxlo2
     &           ,iloboxhi0,iloboxhi1,iloboxhi2
     &           ,haslo
     &           ,ihiboxlo0,ihiboxlo1,ihiboxlo2
     &           ,ihiboxhi0,ihiboxhi1,ihiboxhi2
     &           ,hashi
     &           ,dx
     &           ,facedir
     &           ,divdir
     &           ,gradcomp
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer idivvello0,idivvello1,idivvello2
      integer idivvelhi0,idivvelhi1,idivvelhi2
      REAL_T divvel(
     &           idivvello0:idivvelhi0,
     &           idivvello1:idivvelhi1,
     &           idivvello2:idivvelhi2)
      integer nvelcomp
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL_T vel(
     &           ivello0:ivelhi0,
     &           ivello1:ivelhi1,
     &           ivello2:ivelhi2,
     &           0:nvelcomp-1)
      integer ngradvelcomp
      integer igradvello0,igradvello1,igradvello2
      integer igradvelhi0,igradvelhi1,igradvelhi2
      REAL_T gradvel(
     &           igradvello0:igradvelhi0,
     &           igradvello1:igradvelhi1,
     &           igradvello2:igradvelhi2,
     &           0:ngradvelcomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer icenterboxlo0,icenterboxlo1,icenterboxlo2
      integer icenterboxhi0,icenterboxhi1,icenterboxhi2
      integer iloboxlo0,iloboxlo1,iloboxlo2
      integer iloboxhi0,iloboxhi1,iloboxhi2
      integer haslo
      integer ihiboxlo0,ihiboxlo1,ihiboxlo2
      integer ihiboxhi0,ihiboxhi1,ihiboxhi2
      integer hashi
      REAL_T dx
      integer facedir
      integer divdir
      integer gradcomp
      integer ii,i,jj,j,kk,k
      
      ii = chf_id(facedir, 0)
      jj = chf_id(facedir, 1)
      kk = chf_id(facedir, 2)
      if (facedir .eq. divdir) then
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         divvel(i,j,k) = divvel(i,j,k) +
     $        (    vel(i   ,j   ,k   ,facedir)
     $        -    vel(i-ii,j-jj,k-kk,facedir) )/dx
         
      enddo
      enddo
      enddo
      else
         
      do k = icenterboxlo2,icenterboxhi2
      do j = icenterboxlo1,icenterboxhi1
      do i = icenterboxlo0,icenterboxhi0

         divvel(i,j,k) = divvel(i,j,k) +
     $        ( gradvel(i   ,j   ,k   , gradcomp)
     $        + gradvel(i-ii,j-jj,k-kk, gradcomp) )/two
         
      enddo
      enddo
      enddo
         if(haslo .eq. 1) then
            
      do k = iloboxlo2,iloboxhi2
      do j = iloboxlo1,iloboxhi1
      do i = iloboxlo0,iloboxhi0

            divvel(i,j,k) = divvel(i,j,k) +
     $           (three*gradvel(i   ,j   ,k   , gradcomp)
     $           -      gradvel(i+ii,j+jj,k+kk, gradcomp))/two
            
      enddo
      enddo
      enddo
         endif
         if(hashi .eq. 1) then
            
      do k = ihiboxlo2,ihiboxhi2
      do j = ihiboxlo1,ihiboxhi1
      do i = ihiboxlo0,ihiboxhi0

            divvel(i,j,k) = divvel(i,j,k) +
     $           (three*gradvel(i-  ii,j-  jj,k-  kk, gradcomp)
     $           -      gradvel(i-2*ii,j-2*jj,k-2*kk, gradcomp))/two
            
      enddo
      enddo
      enddo
         endif
      endif
      return
      end
      subroutine DECRINVLAMBDAROP(
     &           lambda
     &           ,ilambdalo0,ilambdalo1,ilambdalo2
     &           ,ilambdahi0,ilambdahi1,ilambdahi2
     &           ,nlambdacomp
     &           ,eta
     &           ,ietalo0,ietalo1,ietalo2
     &           ,ietahi0,ietahi1,ietahi2
     &           ,netacomp
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,beta
     &           ,dx
     &           ,idir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlambdacomp
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1,
     &           ilambdalo2:ilambdahi2,
     &           0:nlambdacomp-1)
      integer netacomp
      integer ietalo0,ietalo1,ietalo2
      integer ietahi0,ietahi1,ietahi2
      REAL_T eta(
     &           ietalo0:ietahi0,
     &           ietalo1:ietahi1,
     &           ietalo2:ietahi2,
     &           0:netacomp-1)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      REAL_T beta
      REAL_T dx
      integer idir
      integer ii,jj,kk
      integer i,j,k
      integer icomp
      REAL_T etahi, etalo
      integer ncomp
      ncomp = nlambdacomp
      
      ii = CHF_ID(idir, 0)
      jj = CHF_ID(idir, 1)
      kk = CHF_ID(idir, 2)
      do icomp = 0, ncomp-1
         
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

         etahi =  eta(i+ii,j+jj,k+kk,0)
         etalo =  eta(i   ,j   ,k   ,0)
         lambda(i,j,k, icomp) = lambda(i,j,k,icomp)
     $        - beta*(
     $        eta(i+ii,j+jj,k+kk,0) +
     $        eta(i   ,j   ,k   ,0))/(dx*dx)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine INVERTLAMBDAROP(
     &           lambda
     &           ,ilambdalo0,ilambdalo1,ilambdalo2
     &           ,ilambdahi0,ilambdahi1,ilambdahi2
     &           ,nlambdacomp
     &           ,safety
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1,
     &           ilambdalo2:ilambdahi2,
     &           0:nlambdacomp-1)
      REAL_T safety
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i,j,k
      integer icomp
      integer ncomp
      ncomp = nlambdacomp
      do icomp = 0, ncomp-1
         
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

         lambda(i,j,k, icomp) =
     $        safety/(lambda(i,j,k,icomp))
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBROP(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,lphi
     &           ,ilphilo0,ilphilo1,ilphilo2
     &           ,ilphihi0,ilphihi1,ilphihi2
     &           ,nlphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,nrhscomp
     &           ,lambda
     &           ,ilambdalo0,ilambdalo1,ilambdalo2
     &           ,ilambdahi0,ilambdahi1,ilambdahi2
     &           ,nlambdacomp
     &           ,icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
     &           ,icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
     &           )

      implicit none
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer nlphicomp
      integer ilphilo0,ilphilo1,ilphilo2
      integer ilphihi0,ilphihi1,ilphihi2
      REAL_T lphi(
     &           ilphilo0:ilphihi0,
     &           ilphilo1:ilphihi1,
     &           ilphilo2:ilphihi2,
     &           0:nlphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           0:nrhscomp-1)
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL_T lambda(
     &           ilambdalo0:ilambdahi0,
     &           ilambdalo1:ilambdahi1,
     &           ilambdalo2:ilambdahi2,
     &           0:nlambdacomp-1)
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
      integer i,j,k
      integer icomp
      integer ncomp
      ncomp = nphicomp
      do icomp = 0, ncomp-1
         
      do k = icoloredBoxlo2,icoloredBoxhi2,2
      do j = icoloredBoxlo1,icoloredBoxhi1,2
      do i = icoloredBoxlo0,icoloredBoxhi0,2

         phi(i,j,k, icomp) =
     $        phi(   i,j,k,icomp) +
     &        lambda(i,j,k,icomp)*(
     $        rhs(   i,j,k,icomp) -
     $        lphi(  i,j,k,icomp))
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine ADDDIVFLUXDIRROP(
     &           lhs
     &           ,ilhslo0,ilhslo1,ilhslo2
     &           ,ilhshi0,ilhshi1,ilhshi2
     &           ,nlhscomp
     &           ,flux
     &           ,ifluxlo0,ifluxlo1,ifluxlo2
     &           ,ifluxhi0,ifluxhi1,ifluxhi2
     &           ,nfluxcomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,dx
     &           ,facedir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nlhscomp
      integer ilhslo0,ilhslo1,ilhslo2
      integer ilhshi0,ilhshi1,ilhshi2
      REAL_T lhs(
     &           ilhslo0:ilhshi0,
     &           ilhslo1:ilhshi1,
     &           ilhslo2:ilhshi2,
     &           0:nlhscomp-1)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           ifluxlo1:ifluxhi1,
     &           ifluxlo2:ifluxhi2,
     &           0:nfluxcomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T dx
      integer facedir
      integer ii,i,jj,j,kk,k
      integer icomp
      integer ncomp
      ncomp = nlhscomp
      
      ii = CHF_ID(facedir, 0)
      jj = CHF_ID(facedir, 1)
      kk = CHF_ID(facedir, 2)
      do icomp = 0, ncomp-1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         lhs(i,j,k, icomp) = lhs(i,j,k, icomp) +
     $        (flux(i+ii,j+jj,k+kk, icomp)
     $        -flux(i   ,j   ,k   , icomp))/dx
         
      enddo
      enddo
      enddo
      enddo
      return
      end
