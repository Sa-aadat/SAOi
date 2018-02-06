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
      double precision temp0, temp1, temp2, temp3
!
! some often occurring terms
      temp0 = 0.124d0
      temp1 = dsqrt(1.d0 + x(2)**2)
      temp2 = 8.d0/x(1) + 1.d0/(x(1)*x(2))
      temp3 = 8.d0/x(1) - 1.d0/(x(1)*x(2))
!
! the objective function
      f    = x(1)*temp1
!
! the first constraint
      c(1) = temp0*temp1*temp2 - 1.d0
!
! the second constraint
      c(2) = temp0*temp1*temp3 - 1.d0
!
      return
      end subroutine SAOi_funcs
