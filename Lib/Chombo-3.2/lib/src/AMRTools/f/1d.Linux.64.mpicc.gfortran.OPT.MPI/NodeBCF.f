      subroutine FACENODEBC(
     & state
     & ,istatelo0
     & ,istatehi0
     & ,nstatecomp
     & ,neumfac
     & ,ineumfaclo0
     & ,ineumfachi0
     & ,nneumfaccomp
     & ,dircfac
     & ,idircfaclo0
     & ,idircfachi0
     & ,ndircfaccomp
     & ,inhmval
     & ,iinhmvallo0
     & ,iinhmvalhi0
     & ,ninhmvalcomp
     & ,ifaceboxlo0
     & ,ifaceboxhi0
     & ,idir
     & ,side
     & ,dx
     & ,startcomp
     & ,endcomp
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nstatecomp
      integer istatelo0
      integer istatehi0
      REAL*8 state(
     & istatelo0:istatehi0,
     & 0:nstatecomp-1)
      integer nneumfaccomp
      integer ineumfaclo0
      integer ineumfachi0
      REAL*8 neumfac(
     & ineumfaclo0:ineumfachi0,
     & 0:nneumfaccomp-1)
      integer ndircfaccomp
      integer idircfaclo0
      integer idircfachi0
      REAL*8 dircfac(
     & idircfaclo0:idircfachi0,
     & 0:ndircfaccomp-1)
      integer ninhmvalcomp
      integer iinhmvallo0
      integer iinhmvalhi0
      REAL*8 inhmval(
     & iinhmvallo0:iinhmvalhi0,
     & 0:ninhmvalcomp-1)
      integer ifaceboxlo0
      integer ifaceboxhi0
      integer idir
      integer side
      REAL*8 dx
      integer startcomp
      integer endcomp
      REAL*8 nfac, dfac, ival, sval,denom,numer
      integer ncomp,nc
      integer i0, ii0
      ncomp = nstatecomp
      if(ncomp .ne. nneumfaccomp) then
          call MAYDAY_ERROR()
      endif
      if(ncomp .ne. ndircfaccomp) then
          call MAYDAY_ERROR()
      endif
      if(ncomp .ne. ninhmvalcomp) then
          call MAYDAY_ERROR()
      endif
      if ((side .ne. -1) .and. (side .ne. 1)) then
          call MAYDAY_ERROR()
      endif
      ii0= side*CHF_ID(0, idir)
      do nc = startcomp, endcomp
      do i0 = ifaceboxlo0,ifaceboxhi0
              nfac = neumfac(i0, nc)
              dfac = dircfac(i0, nc)
              ival = inhmval(i0, nc)
              sval = state(i0-ii0, nc)
              denom = dfac + side*(nfac/dx)
              numer = ival + side*(nfac/dx)*sval
              if (abs(denom) .lt. 1.0e-9) then
                  call MAYDAY_ERROR()
              endif
              state(i0, nc) = numer/denom
      enddo
      enddo
      return
      end
