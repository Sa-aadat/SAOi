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
      double precision x1, x2, t1, t2, t3, t4, t5, t6, t7
!
      x1=x(1)
      x2=x(2)
      t1=x1+x2+1.d0
      t2=19.d0-14.d0*x1+3.d0*x1*x1-14.d0*x2+6.d0*x1*x2+3.d0*x2*x2
      t3=2.d0*x1-3.d0*x2
      t4=18.d0-32.d0*x1+12.d0*x1*x1+48.d0*x2-36.d0*x1*x2+27.d0*x2*x2
      t5=-14.d0+6.d0*x1+6.d0*x2
      t6=-32.d0+24.d0*x1-36.d0*x2
      t7=48.d0-36.d0*x1+54.d0*x2
!
      f=(1.d0+(t1*t1)*t2)*(30.d0+(t3*t3)*t4)
!
      return
      end subroutine SAOi_funcs
