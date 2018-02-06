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
      double precision f, x(n), c(ni), x1, x2, x3, x4
!
      x1 = x(1)
      x2 = x(2)
      x3 = x(3)
      x4 = x(4)
!      
      f = (x1 + 10.d0*x2)**2 +  5.d0*(x3 - x4)**2 
     &  + (x2 -  2.d0*x3)**4 + 10.d0*(x1 - x4)**4
!
      return
      end subroutine SAOi_funcs
