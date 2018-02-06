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
      double precision f, x(n), c(ni), x1, x2, x3, pi, theta
!
      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      pi = 4.0d0*datan(1.0d0)
!
      if (x1.gt.0.d0) then
        theta = (1.d0/(2.d0*pi))*datan(x2/x1)
      elseif (x1.lt.0.d0) then
        theta = (1.d0/(2.d0*pi))*(pi + datan(x2/x1))
      else
        theta = 0.25;
      endif
!
      f = 100.d0*((x3 - 10.d0*theta)**2 
     &  + (dsqrt(x1**2 + x2**2) - 1.d0)**2) + x3**2
!
      return
      end subroutine SAOi_funcs
