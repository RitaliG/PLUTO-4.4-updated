      subroutine AVEFACESCALTOFACEVECT(
     & facevect
     & ,ifacevectlo0
     & ,ifacevecthi0
     & ,nfacevectcomp
     & ,facescal
     & ,ifacescallo0
     & ,ifacescalhi0
     & ,facedir
     & ,vectdir
     & ,idcalcfacelo0
     & ,idcalcfacehi0
     & ,ioffboxlo0
     & ,ioffboxhi0
     & )
      implicit none
      integer nfacevectcomp
      integer ifacevectlo0
      integer ifacevecthi0
      REAL*8 facevect(
     & ifacevectlo0:ifacevecthi0,
     & 0:nfacevectcomp-1)
      integer ifacescallo0
      integer ifacescalhi0
      REAL*8 facescal(
     & ifacescallo0:ifacescalhi0)
      integer facedir
      integer vectdir
      integer idcalcfacelo0
      integer idcalcfacehi0
      integer ioffboxlo0
      integer ioffboxhi0
      integer i0
      integer ioff0
      integer numpts
      REAL*8 wt, tot
      if (facedir .eq. vectdir) then
      do i0 = idcalcfacelo0,idcalcfacehi0
            facevect(i0, vectdir) = facescal(i0)
      enddo
      else
         numpts = 2**1
         wt = (1.0d0) / (numpts * (1.0d0))
      do i0 = idcalcfacelo0,idcalcfacehi0
            tot = (0.0d0)
      do ioff0 = ioffboxlo0,ioffboxhi0
               tot = tot + facescal(i0 +ioff0)
      enddo
            facevect(i0, vectdir) = wt * tot
      enddo
      endif
      return
      end
      subroutine AVESCALTOFACE(
     & facescal
     & ,ifacescallo0
     & ,ifacescalhi0
     & ,cellscal
     & ,icellscallo0
     & ,icellscalhi0
     & ,idir
     & ,idcalcfacelo0
     & ,idcalcfacehi0
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacescallo0
      integer ifacescalhi0
      REAL*8 facescal(
     & ifacescallo0:ifacescalhi0)
      integer icellscallo0
      integer icellscalhi0
      REAL*8 cellscal(
     & icellscallo0:icellscalhi0)
      integer idir
      integer idcalcfacelo0
      integer idcalcfacehi0
      integer i0
      integer f2cLo0
      integer f2cHi0
      f2cLo0= -1*CHF_ID(0, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      do i0 = idcalcfacelo0,idcalcfacehi0
         facescal(i0) = (0.500d0) * (
     & cellscal(i0 +f2cLo0) +
     & cellscal(i0 +f2cHi0) )
      enddo
      return
      end
      subroutine AVECELLTOFACE(
     & facevel
     & ,ifacevello0
     & ,ifacevelhi0
     & ,cellvel
     & ,icellvello0
     & ,icellvelhi0
     & ,ncellvelcomp
     & ,idir
     & ,idcalcfacelo0
     & ,idcalcfacehi0
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacevello0
      integer ifacevelhi0
      REAL*8 facevel(
     & ifacevello0:ifacevelhi0)
      integer ncellvelcomp
      integer icellvello0
      integer icellvelhi0
      REAL*8 cellvel(
     & icellvello0:icellvelhi0,
     & 0:ncellvelcomp-1)
      integer idir
      integer idcalcfacelo0
      integer idcalcfacehi0
      integer i0
      integer f2cLo0
      integer f2cHi0
      f2cLo0= -1*CHF_ID(0, idir)
      f2cHi0= 0*CHF_ID(0, idir)
      do i0 = idcalcfacelo0,idcalcfacehi0
         facevel(i0) = (0.500d0) * (
     & cellvel(i0 +f2cLo0, idir) +
     & cellvel(i0 +f2cHi0, idir) )
      enddo
      return
      end
      subroutine AVEFACETOCELL(
     & cellvel
     & ,icellvello0
     & ,icellvelhi0
     & ,ncellvelcomp
     & ,facevel
     & ,ifacevello0
     & ,ifacevelhi0
     & ,idir
     & ,idcalccelllo0
     & ,idcalccellhi0
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ncellvelcomp
      integer icellvello0
      integer icellvelhi0
      REAL*8 cellvel(
     & icellvello0:icellvelhi0,
     & 0:ncellvelcomp-1)
      integer ifacevello0
      integer ifacevelhi0
      REAL*8 facevel(
     & ifacevello0:ifacevelhi0)
      integer idir
      integer idcalccelllo0
      integer idcalccellhi0
      integer i0
      integer c2fLo0
      integer c2fHi0
      c2fLo0= 0*CHF_ID(0, idir)
      c2fHi0= 1*CHF_ID(0, idir)
      do i0 = idcalccelllo0,idcalccellhi0
         cellvel(i0, idir) = (0.500d0) * (
     & facevel(i0 +c2fLo0) +
     & facevel(i0 +c2fHi0) )
      enddo
      return
      end
      subroutine MAGNITUDEF(
     & magdata
     & ,imagdatalo0
     & ,imagdatahi0
     & ,data
     & ,idatalo0
     & ,idatahi0
     & ,ndatacomp
     & ,iboxlo0
     & ,iboxhi0
     & )
      implicit none
      integer imagdatalo0
      integer imagdatahi0
      REAL*8 magdata(
     & imagdatalo0:imagdatahi0)
      integer ndatacomp
      integer idatalo0
      integer idatahi0
      REAL*8 data(
     & idatalo0:idatahi0,
     & 0:ndatacomp-1)
      integer iboxlo0
      integer iboxhi0
      integer i0
      integer iv
      REAL*8 cur,sum
      do i0 = iboxlo0,iboxhi0
         sum = (0.0d0)
         do iv = 0, ndatacomp-1
            cur = data(i0, iv)
            sum = sum + cur*cur
         enddo
         magdata(i0) = sqrt(sum)
      enddo
      return
      end
        subroutine GETRELGRADF(
     & du
     & ,idulo0
     & ,iduhi0
     & ,u
     & ,iulo0
     & ,iuhi0
     & ,idir
     & ,iloBoxlo0
     & ,iloBoxhi0
     & ,hasLo
     & ,ihiBoxlo0
     & ,ihiBoxhi0
     & ,hasHi
     & ,icenterBoxlo0
     & ,icenterBoxhi0
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idulo0
      integer iduhi0
      REAL*8 du(
     & idulo0:iduhi0)
      integer iulo0
      integer iuhi0
      REAL*8 u(
     & iulo0:iuhi0)
      integer idir
      integer iloBoxlo0
      integer iloBoxhi0
      integer hasLo
      integer ihiBoxlo0
      integer ihiBoxhi0
      integer hasHi
      integer icenterBoxlo0
      integer icenterBoxhi0
        integer i0
        integer ioff0
        REAL*8 diff,aver
      ioff0=CHF_ID(0, idir)
      do i0 = icenterBoxlo0,icenterBoxhi0
          diff = (0.500d0)*(u(i0 +ioff0)
     & -u(i0 -ioff0))
          aver = (0.500d0)*(u(i0 +ioff0)
     & +u(i0 -ioff0))
          du(i0) = diff / aver
      enddo
        if (hasLo .eq. 1) then
      do i0 = iloBoxlo0,iloBoxhi0
            diff = u(i0 +ioff0) - u(i0)
            aver = (0.500d0)*(u(i0 +ioff0) + u(i0))
            du(i0) = diff / aver
      enddo
        endif
        if (hasHi .eq. 1) then
      do i0 = ihiBoxlo0,ihiBoxhi0
            diff = u(i0) - u(i0 -ioff0)
            aver = (0.500d0)*(u(i0) + u(i0 -ioff0))
            du(i0) = diff / aver
      enddo
        endif
        return
        end
      subroutine POSTNORMALSOURCE(
     & dWminus
     & ,idWminuslo0
     & ,idWminushi0
     & ,ndWminuscomp
     & ,dWplus
     & ,idWpluslo0
     & ,idWplushi0
     & ,ndWpluscomp
     & ,W
     & ,iWlo0
     & ,iWhi0
     & ,nWcomp
     & ,advVel
     & ,iadvVello0
     & ,iadvVelhi0
     & ,dt
     & ,dx
     & ,idir
     & ,iboxlo0
     & ,iboxhi0
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ndWminuscomp
      integer idWminuslo0
      integer idWminushi0
      REAL*8 dWminus(
     & idWminuslo0:idWminushi0,
     & 0:ndWminuscomp-1)
      integer ndWpluscomp
      integer idWpluslo0
      integer idWplushi0
      REAL*8 dWplus(
     & idWpluslo0:idWplushi0,
     & 0:ndWpluscomp-1)
      integer nWcomp
      integer iWlo0
      integer iWhi0
      REAL*8 W(
     & iWlo0:iWhi0,
     & 0:nWcomp-1)
      integer iadvVello0
      integer iadvVelhi0
      REAL*8 advVel(
     & iadvVello0:iadvVelhi0)
      REAL*8 dt
      REAL*8 dx
      integer idir
      integer iboxlo0
      integer iboxhi0
      integer i0
      integer c2fLo0
      integer c2fHi0
      integer n
      REAL*8 dudx
      c2fLo0= 0*CHF_ID(0, idir)
      c2fHi0= 1*CHF_ID(0, idir)
      do i0 = iboxlo0,iboxhi0
         dudx = advVel(i0 +c2fHi0)
     & - advVel(i0 +c2fLo0)
         dudx = dt*dudx / dx * (0.500d0)
         do n=0, nWcomp-1
            dWplus(i0, n) =
     & dWPlus(i0, n) - dudx * W(i0,n)
            dWMinus(i0, n) =
     & dWMinus(i0, n) - dudx * W(i0,n)
         end do
      enddo
      return
      end
      subroutine RIEMANNF(
     & Wgdnv
     & ,iWgdnvlo0
     & ,iWgdnvhi0
     & ,nWgdnvcomp
     & ,WLeft
     & ,iWLeftlo0
     & ,iWLefthi0
     & ,nWLeftcomp
     & ,WRight
     & ,iWRightlo0
     & ,iWRighthi0
     & ,nWRightcomp
     & ,advVel
     & ,iadvVello0
     & ,iadvVelhi0
     & ,idir
     & ,iboxlo0
     & ,iboxhi0
     & )
      implicit none
      integer nWgdnvcomp
      integer iWgdnvlo0
      integer iWgdnvhi0
      REAL*8 Wgdnv(
     & iWgdnvlo0:iWgdnvhi0,
     & 0:nWgdnvcomp-1)
      integer nWLeftcomp
      integer iWLeftlo0
      integer iWLefthi0
      REAL*8 WLeft(
     & iWLeftlo0:iWLefthi0,
     & 0:nWLeftcomp-1)
      integer nWRightcomp
      integer iWRightlo0
      integer iWRighthi0
      REAL*8 WRight(
     & iWRightlo0:iWRighthi0,
     & 0:nWRightcomp-1)
      integer iadvVello0
      integer iadvVelhi0
      REAL*8 advVel(
     & iadvVello0:iadvVelhi0)
      integer idir
      integer iboxlo0
      integer iboxhi0
      integer i0
      integer n
      REAL*8 sl,sr
      REAL*8 so
      REAL*8 ustar
      do n=0, nWLeftcomp-1
      do i0 = iboxlo0,iboxhi0
            sl = WLeft(i0, n)
            sr = WRight(i0, n)
            ustar = advVel(i0)
            if (ustar .gt. (0.0d0)) then
               so = sl
            else
               so = sr
            endif
            if (abs(ustar) .lt. 1.0d-9) then
               so = (0.500d0)*(sl+sr)
            endif
            Wgdnv(i0, n) = so
      enddo
      end do
      return
      end
      subroutine QUASILINEARUPDATE(
     & AdWdx
     & ,iAdWdxlo0
     & ,iAdWdxhi0
     & ,nAdWdxcomp
     & ,WHalf
     & ,iWHalflo0
     & ,iWHalfhi0
     & ,nWHalfcomp
     & ,cellVel
     & ,icellVello0
     & ,icellVelhi0
     & ,scale
     & ,idir
     & ,iboxlo0
     & ,iboxhi0
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nAdWdxcomp
      integer iAdWdxlo0
      integer iAdWdxhi0
      REAL*8 AdWdx(
     & iAdWdxlo0:iAdWdxhi0,
     & 0:nAdWdxcomp-1)
      integer nWHalfcomp
      integer iWHalflo0
      integer iWHalfhi0
      REAL*8 WHalf(
     & iWHalflo0:iWHalfhi0,
     & 0:nWHalfcomp-1)
      integer icellVello0
      integer icellVelhi0
      REAL*8 cellVel(
     & icellVello0:icellVelhi0)
      REAL*8 scale
      integer idir
      integer iboxlo0
      integer iboxhi0
      integer i0
      integer c2fLo0
      integer c2fHi0
      integer n
      c2fLo0= 0*CHF_ID(0, idir)
      c2fHi0= 1*CHF_ID(0, idir)
      do n=0, nAdWdxcomp-1
      do i0 = iboxlo0,iboxhi0
            AdWdx(i0, n) =
     & scale * cellVel(i0) *
     & (WHalf(i0 +c2fHi0, n) -
     & WHalf(i0 +c2fLo0, n))
      enddo
      end do
      return
      end
      subroutine AVECELLVECTOFACEVEC(
     & facevec
     & ,ifaceveclo0
     & ,ifacevechi0
     & ,nfaceveccomp
     & ,cellvec
     & ,icellveclo0
     & ,icellvechi0
     & ,ncellveccomp
     & ,facedir
     & ,idcalcfacelo0
     & ,idcalcfacehi0
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfaceveccomp
      integer ifaceveclo0
      integer ifacevechi0
      REAL*8 facevec(
     & ifaceveclo0:ifacevechi0,
     & 0:nfaceveccomp-1)
      integer ncellveccomp
      integer icellveclo0
      integer icellvechi0
      REAL*8 cellvec(
     & icellveclo0:icellvechi0,
     & 0:ncellveccomp-1)
      integer facedir
      integer idcalcfacelo0
      integer idcalcfacehi0
      integer i0
      integer f2cLo0
      integer f2cHi0
      integer idir
      f2cLo0= -1*CHF_ID(0, facedir)
      f2cHi0= 0*CHF_ID(0, facedir)
      do idir = 0, 1 -1
      do i0 = idcalcfacelo0,idcalcfacehi0
            facevec(i0, idir) = (0.500d0) * (
     & cellvec(i0 +f2cLo0, idir) +
     & cellvec(i0 +f2cHi0, idir) )
      enddo
      enddo
      return
      end
