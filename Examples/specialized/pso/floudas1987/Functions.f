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
!
! the objective function
      f    = (x(1) - 10.d0)**3 + (x(2) - 20.d0)**3
!
! the first constraint
      c(1) =   -(x(1) - 5.d0)**2 - (x(2) - 5.d0)**2 + 100.d0 
!
! the second constraint
      c(2) =   (x(1) - 6.d0)**2 + (x(2) - 5.d0)**2 - 82.81d0
!
      return
      end subroutine SAOi_funcs
