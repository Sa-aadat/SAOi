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
      double precision y, one, x1, x2, x3
      data             one /1.d0/
!
      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
!
      f = dexp(x1) + dexp(x2) + dexp(x3)
!
      do j=1,ni
        y = dble(j-1)/dble(ni-1)
        c(j) = - x1 - x2*y - x3*y**2 + one/(one+y**2)
      enddo
!
      return
      end subroutine SAOi_funcs
