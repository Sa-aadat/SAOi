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
      integer          i, j, n, m, ni, ne
      double precision f, x(n), c(ni)
      double precision pi, rads
!
      pi = 4.0d0*datan(1.0d0)
!
      rads=pi/180.d0
      f=0.d0
      do i=1,10
        f=f+x(i)*x(i+1)*dsin(x(i+11)*rads)
      end do
      f=-0.5d0*f
      c(1)=0.d0
      do i=1,10
        c(1)=c(1)+dsqrt(x(i)**2+x(i+1)**2
     &  -2.d0*x(i)*x(i+1)*dcos(x(i+11)*rads))
      end do
      c(1)=c(1)+x(1)+x(11)-60.d0
!
      return
      end subroutine SAOi_funcs
