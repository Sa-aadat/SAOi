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
      double precision f, x(n), c(ni), x1, x2
!
      x1 = x(1)
      x2 = x(2)
!
      f = (1.5d0 - x1*(1.d0 - x2))**2 + (2.25d0 - x1*(1.d0 - x2**2))**2
     &  + (2.625d0 - x1*(1.d0 - x2**3))**2
!
      return
      end subroutine SAOi_funcs
