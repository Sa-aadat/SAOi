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

! the objective function
      f = 5.d0*x(1) - 80.d0*x(2) + 160.d0*x(3) + 15.d0*x(4)
  
! the constraints
      c(1) = 0.05d0*x(1) + 0.05d0*x(2) + 0.1d0*x(3) + 0.15d0*x(4) 
      c(2) = 0.10d0*x(1) + 0.15d0*x(2) + 0.2d0*x(3) + 0.05d0*x(4)  
      c(3) = 0.05d0*x(1) + 0.10d0*x(2) + 0.1d0*x(3) + 0.15d0*x(4)  

! the constraints in negative-null form
      c(1) = c(1) - 1000.d0 
      c(2) = c(2) - 2000.d0 
      c(3) = c(3) - 1500.d0 
!
      return
      end subroutine SAOi_funcs
