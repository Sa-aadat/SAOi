! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: objective and constraint functions
! SAOi:

      subroutine SAOi_funcs (n, m, ni, ne, x, f, c, iuser, luser, cuser,
     &                       ruser, eqn, lin, ictrl, lctrl, rctrl,
     &                       cctrl)
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the objective function f and the inequality constraint      !
!  functions c(j), j=1,ni                                              !
!                                                                      !
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          i, j, n, m, ni, ne, iw, npt
      double precision f, x(n), c(ni)
      double precision fsum, y(10,10)
!
      do j=1,n
        y(1,j)=1.0d0
        y(2,j)=2.0d0*x(j)-1.0d0
      enddo
      do i=2,n
        do j=1,n
          y(i+1,j)=2.0d0*y(2,j)*y(i,j)-y(i-1,j)
        enddo
      enddo
      f=0.0d0
      npt=n+1
      iw=1
      do i=1,npt
        fsum=0.0d0
        do j=1,n
          fsum=fsum+y(i,j)
        enddo
        fsum=fsum/dfloat(n)
        if (iw.gt.0) fsum=fsum+1.0d0/dfloat(i*i-2*i)
        iw=-iw
        f=f+fsum*fsum
      enddo
!
      return
      end subroutine SAOi_funcs
