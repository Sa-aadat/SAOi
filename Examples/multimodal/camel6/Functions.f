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
      double precision x2, x22, x23, x24, x1, x12, x13, x14, x15, x16
!
      x1=x(1)
      x12=x1*x1
      x13=x12*x1
      x14=x12*x12
      x15=x12*x13
      x16=x13*x13
      x2=x(2)
      x22=x2*x2
      x23=x22*x2
      x24=x22*x22
      f=4.d0*x12-2.1d0*x14+x16/3.d0+x1*x2-4.d0*x22+4.d0*x24
!
      return
      end subroutine SAOi_funcs
