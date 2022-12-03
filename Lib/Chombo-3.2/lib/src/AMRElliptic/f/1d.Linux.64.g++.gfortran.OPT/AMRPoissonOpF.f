      subroutine GSMCAMRPOP(
     & phi
     & ,iphilo0
     & ,iphihi0
     & ,rhs
     & ,irhslo0
     & ,irhshi0
     & ,icoloredboxlo0
     & ,icoloredboxhi0
     & ,dx
     & ,alpha
     & ,beta
     & )
      implicit none
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0)
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     & irhslo0:irhshi0)
      integer icoloredboxlo0
      integer icoloredboxhi0
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      REAL*8 lambda, dxinv, sum_b, lphi
      integer i
      integer idir
      dxinv = (1.0d0)/(dx*dx)
      sum_b = 0.0
      do idir = 0, 1 -1
         sum_b = sum_b + (2.0d0)*dxinv
      enddo
      lambda = (1.0d0)/(alpha - beta*sum_b)
      do i = icoloredBoxlo0,icoloredBoxhi0,2
        lphi =
     & ( phi(i+1)
     & + phi(i-1)
     $ -(2.0d0)*phi(i ))
        lphi = lphi*dxinv
        phi(i) =
     $ phi( i) +
     & lambda*( rhs( i) - lphi)
      enddo
      return
      end
      subroutine GSRBHELMHOLTZ(
     & phi
     & ,iphilo0
     & ,iphihi0
     & ,nphicomp
     & ,rhs
     & ,irhslo0
     & ,irhshi0
     & ,nrhscomp
     & ,iregionlo0
     & ,iregionhi0
     & ,dx
     & ,alpha
     & ,beta
     & ,redBlack
     & )
      implicit none
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0,
     & 0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & 0:nrhscomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      integer redBlack
      REAL*8 lambda, dxinv, sum_b, lphi, helmop
      integer i
      integer n,ncomp,idir,indtot,imin,imax
      dxinv = (1.0d0)/(dx*dx)
      sum_b = 0.0
      do idir = 0, 1 -1
         sum_b = sum_b + (2.0d0)*dxinv
      enddo
      lambda = -(1.0d0)/(alpha - beta*sum_b)
      ncomp = nphicomp
      if (ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
              lphi = (phi(i+1,n)
     & + phi(i-1,n)
     & -(2.0d0)*1*phi(i,n))*dxinv
              helmop = alpha*phi(i,n) + beta*lphi
              phi(i,n) = phi(i,n) +
     & lambda*(helmop - rhs(i,n))
            enddo
      enddo
      return
      end
      subroutine GSRBLAPLACIAN(
     & phi
     & ,iphilo0
     & ,iphihi0
     & ,nphicomp
     & ,rhs
     & ,irhslo0
     & ,irhshi0
     & ,nrhscomp
     & ,iregionlo0
     & ,iregionhi0
     & ,dx
     & ,redBlack
     & )
      implicit none
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0,
     & 0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & 0:nrhscomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      integer redBlack
      REAL*8 lambda, dxinv, sum_b, lphi, lap
      integer i
      integer n,ncomp,idir,indtot,imin,imax
      dxinv = (1.0d0)/(dx*dx)
      sum_b = 0.0
      do idir = 0, 1 -1
         sum_b = sum_b + (2.0d0)*dxinv
      enddo
      lambda = -(1.0d0)/sum_b
      ncomp = nphicomp
      if (ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
            imin = iregionlo0
            indtot = imin
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
              lphi = (
     & phi(i+1,n)
     & + phi(i-1,n)
     & ) * dxinv
              phi(i,n) = lambda*(rhs(i,n)-lphi)
            enddo
      enddo
      return
      end
      subroutine GSRBLAZY(
     & phi
     & ,iphilo0
     & ,iphihi0
     & ,lphi
     & ,ilphilo0
     & ,ilphihi0
     & ,rhs
     & ,irhslo0
     & ,irhshi0
     & ,icoloredboxlo0
     & ,icoloredboxhi0
     & ,alpha
     & ,beta
     & ,dx
     & )
      implicit none
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0)
      integer ilphilo0
      integer ilphihi0
      REAL*8 lphi(
     & ilphilo0:ilphihi0)
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     & irhslo0:irhshi0)
      integer icoloredboxlo0
      integer icoloredboxhi0
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dx
      integer i, idir
      REAL*8 dxinv, sum_b, lambda
      dxinv = (1.0d0)/(dx*dx)
      sum_b = 0.0
      do idir = 0, 1 -1
         sum_b = sum_b + (2.0d0)*dxinv
      enddo
      lambda = -(1.0d0)/(alpha - beta*sum_b)
      do i = icoloredBoxlo0,icoloredBoxhi0,2
      phi(i) =
     $ phi( i) -
     & lambda*(
     $ rhs( i) -
     $ lphi( i))
      enddo
      return
      end
      subroutine AMRPMULTICOLOR(
     & phi
     & ,iphilo0
     & ,iphihi0
     & ,rhs
     & ,irhslo0
     & ,irhshi0
     & ,weight
     & ,alpha
     & ,beta
     & ,dx
     & ,icoloredboxlo0
     & ,icoloredboxhi0
     & )
      implicit none
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0)
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     & irhslo0:irhshi0)
      REAL*8 weight
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dx(0:0)
      integer icoloredboxlo0
      integer icoloredboxhi0
      integer i
      REAL*8 laplphi, dx0
      dx0 = beta/(dx(0) * dx(0))
      do i = icoloredBoxlo0,icoloredBoxhi0,2
        laplphi =
     & ( phi(i+1)
     & + phi(i-1)
     $ -(2.0d0)*phi(i ))*dx0
        laplphi = laplphi + alpha * phi(i)
        phi(i) = phi(i) +
     & weight*(rhs(i) - laplphi)
      enddo
      return
      end
      subroutine OPERATORLAP(
     & lofphi
     & ,ilofphilo0
     & ,ilofphihi0
     & ,nlofphicomp
     & ,phi
     & ,iphilo0
     & ,iphihi0
     & ,nphicomp
     & ,iregionlo0
     & ,iregionhi0
     & ,dx
     & ,alpha
     & ,beta
     & )
      implicit none
      integer nlofphicomp
      integer ilofphilo0
      integer ilofphihi0
      REAL*8 lofphi(
     & ilofphilo0:ilofphihi0,
     & 0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0,
     & 0:nphicomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dxinv, lap
      integer n,ncomp
      integer i
      ncomp = nphicomp
      if (ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          lap = ( phi(i+1,n)
     & + phi(i-1,n)
     & -(2*1)*phi(i,n) )
     & * dxinv
          lofphi(i,n) = alpha*phi(i,n)+beta*lap
      enddo
      enddo
      return
      end
      subroutine OPERATORLAPRES(
     & r
     & ,irlo0
     & ,irhi0
     & ,nrcomp
     & ,phi
     & ,iphilo0
     & ,iphihi0
     & ,nphicomp
     & ,rhs
     & ,irhslo0
     & ,irhshi0
     & ,nrhscomp
     & ,iregionlo0
     & ,iregionhi0
     & ,dx
     & ,alpha
     & ,beta
     & )
      implicit none
      integer nrcomp
      integer irlo0
      integer irhi0
      REAL*8 r(
     & irlo0:irhi0,
     & 0:nrcomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0,
     & 0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & 0:nrhscomp-1)
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 alpha
      REAL*8 beta
      REAL*8 dxinv, lap
      integer n,ncomp
      integer i
      ncomp = nphicomp
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          lap = ( phi(i+1,n)
     & + phi(i-1,n)
     & -(2*1)*phi(i,n) )
     & * dxinv
         r(i,n) = -alpha*phi(i,n) -beta*lap +
     & rhs(i,n)
      enddo
      enddo
      return
      end
      subroutine RESTRICTRES(
     & res
     & ,ireslo0
     & ,ireshi0
     & ,nrescomp
     & ,phi
     & ,iphilo0
     & ,iphihi0
     & ,nphicomp
     & ,rhs
     & ,irhslo0
     & ,irhshi0
     & ,nrhscomp
     & ,alpha
     & ,beta
     & ,iregionlo0
     & ,iregionhi0
     & ,dx
     & )
      implicit none
      integer nrescomp
      integer ireslo0
      integer ireshi0
      REAL*8 res(
     & ireslo0:ireshi0,
     & 0:nrescomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0,
     & 0:nphicomp-1)
      integer nrhscomp
      integer irhslo0
      integer irhshi0
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & 0:nrhscomp-1)
      REAL*8 alpha
      REAL*8 beta
      integer iregionlo0
      integer iregionhi0
      REAL*8 dx
      REAL*8 denom,dxinv,lofphi
      integer n,ncomp
      integer i
      integer ii
      ncomp = nphicomp
      dxinv = (1.0d0) / (dx*dx)
      denom = 2
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          ii = i/2
          lofphi = alpha * phi(i,n)
     & + beta *
     & ( phi(i+1,n)
     & + phi(i-1,n)
     & - phi(i ,n) * 2 * 1
     & ) * dxinv
          res(ii,n) = res(ii,n)
     & + (rhs(i,n) - lofphi) / denom
      enddo
      enddo
      return
      end
      subroutine PROLONG(
     & phi
     & ,iphilo0
     & ,iphihi0
     & ,nphicomp
     & ,coarse
     & ,icoarselo0
     & ,icoarsehi0
     & ,ncoarsecomp
     & ,iregionlo0
     & ,iregionhi0
     & ,m
     & )
      implicit none
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0,
     & 0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & 0:ncoarsecomp-1)
      integer iregionlo0
      integer iregionhi0
      integer m
      INTEGER ncomp, n
      integer i
      integer ii
      ncomp = nphicomp
      do n = 0, ncomp-1
      do i = iregionlo0,iregionhi0
          ii = i/m
          phi(i,n) = phi(i,n) +
     & coarse(ii,n)
      enddo
      enddo
      return
      end
      subroutine PROLONG_2(
     & phi
     & ,iphilo0
     & ,iphihi0
     & ,nphicomp
     & ,coarse
     & ,icoarselo0
     & ,icoarsehi0
     & ,ncoarsecomp
     & ,iregionlo0
     & ,iregionhi0
     & ,m
     & )
      implicit none
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0,
     & 0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0
      integer icoarsehi0
      REAL*8 coarse(
     & icoarselo0:icoarsehi0,
     & 0:ncoarsecomp-1)
      integer iregionlo0
      integer iregionhi0
      integer m
      INTEGER ncomp, n
      integer i
      integer offs(1)
      integer ic
      REAL*8 f0, den, fx(1)
      den = (1.0d0)/(4**1)
      fx(1) = (3.0d0)*den
      f0 = (1.0d0)*den
      ncomp = nphicomp
      do i = iregionlo0,iregionhi0
        ic = i/m
        offs(1) = 2*mod(i,2) - 1
        do n = 0, ncomp-1
          phi(i,n) = phi(i,n)
     $ + fx(1)*
     $ coarse(ic,n)
     $ + f0*coarse(ic+offs(1),n)
        enddo
      enddo
      return
      end
      subroutine NEWGETFLUX(
     & flux
     & ,ifluxlo0
     & ,ifluxhi0
     & ,nfluxcomp
     & ,phi
     & ,iphilo0
     & ,iphihi0
     & ,nphicomp
     & ,iboxlo0
     & ,iboxhi0
     & ,beta_dx
     & ,a_idir
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfluxcomp
      integer ifluxlo0
      integer ifluxhi0
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & 0:nfluxcomp-1)
      integer nphicomp
      integer iphilo0
      integer iphihi0
      REAL*8 phi(
     & iphilo0:iphihi0,
     & 0:nphicomp-1)
      integer iboxlo0
      integer iboxhi0
      REAL*8 beta_dx
      integer a_idir
      INTEGER ncomp,n
      integer ii
      integer i
      ncomp = nphicomp
      ii = CHF_ID(a_idir, 0)
      do n = 0, ncomp-1
      do i = iboxlo0,iboxhi0
          flux(i,n) =
     & (phi(i,n)-
     & phi(i-ii,n))*beta_dx
      enddo
      enddo
      return
      end
