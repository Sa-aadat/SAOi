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
      double precision c1, c2
!
      c1 = 0.0624d0
      f  = c1*(x(1)+x(2)+x(3)+x(4)+x(5))
      c2 = 1.d0
      c(1) = 61.d0/x(1)**3+37.d0/x(2)**3+19/x(3)**3+7.d0/x(4)**3
     &     +  1.d0/x(5)**3-c2
!
      return
      end subroutine SAOi_funcs
