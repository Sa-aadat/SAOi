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
      double precision pi, x1, x2, t1, t2
!
      pi = 4.0d0*datan(1.0d0)
!
      x1 = x(1)
      x2 = x(2)
      t1 = x2-(5.1d0/(4.d0*pi*pi))*x1*x1 + 5.d0/pi*x1 - 6.d0
      t2 = 10.d0*(1.d0-(1.d0/(8.d0*pi)))
      f  = t1*t1 + t2*dcos(x1) + 10.d0
!
      return
      end subroutine SAOi_funcs
