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
      f=5.3578547*x(3)**2+0.8356891*x(1)*x(5)+37.2932239*x(1)
     &  -40792.141
!
! the constraints
      C(1)=-(85.334407+0.0056858*x(2)*x(5)+0.00026*x(1)*x(4)
     &     -0.0022053*x(3)*x(5))
      C(2)=85.334407+0.0056858*x(2)*x(5)+0.00026*x(1)*x(4)
     &    -0.0022053*x(3)*x(5)-92
      C(3)=90-(80.51249+0.0071317*x(2)*x(5)+0.0029955*x(1)*x(2)
     &    +0.0021813*x(3))
      C(4)=80.51249+0.0071317*x(2)*x(5)+0.0029955*x(1)*x(2)
     &    +0.0021813*x(3)-110
      C(5)=20-(9.300961+0.0047026*x(3)*x(5)+0.0012547*x(1)*x(3)
     &    +0.0019085*x(3)*x(4))
      C(6)=9.300961+0.0047026*x(3)*x(5)+0.0012547*x(1)*x(3)
     &    +0.0019085*x(3)*x(4)-25
!
      return
      end subroutine SAOi_funcs
