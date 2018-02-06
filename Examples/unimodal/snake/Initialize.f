! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: initialization
! SAOi:

      subroutine SAOi_init (n, ni, ne, x, x_lower, x_upper, ictrl,
     &                      lctrl, rctrl, cctrl, iuser, luser, cuser,
     &                      ruser, nnz, nnzh, eqn, lin, shift)
!----------------------------------------------------------------------!
!                                                                      !
!  Initialization of the SAOi algorithm; please see the users manual   !
!  for type declarations and comments                                  !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          i, j, n, m, ni, ne, nnz, nnzh, l2030
      double precision x(*), x_lower(*), x_upper(*), shift(*)
      double precision ai
      double precision pi
!
      pi = 4.0d0*datan(1.0d0)
!
      include          'ctrl_get.inc'
!
      l2030          = 10            !   problem specific parameter
      n              = l2030*3       !  the number of design variables
      ni             = l2030*4+1     !  the number of inequality constraints

!  the starting point
      do i=1,l2030
        ai           = ((3.d0*dble(i)-2.d0*dble(l2030))*pi)/(6.d0*l2030)
        x(i)         = dcos(ai+pi/12.d0)
        x(i+l2030)   = dsin(ai+pi/12.d0)
        x(i+2*l2030) = dsin(2.d0*ai+pi/6.d0)
      enddo
!
      do i=1,n
        x_lower(i)   = -2.d0          !  the lower bounds
      end do
!
      do i=1,n
        x_upper(i)   =  2.d0          !  the upper bounds
      end do
!
      cname1 = 'SNAKE'
!
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
!
!     finite_diff    =  .false.  ! use finite differences?
      force_converge =  2        ! enforce convergence? 0 = No; 1 through 3 = Yes (see manual)
!     approx_f       =  1        ! select approx_f (see manual)
!     approx_c       =  1        ! select approx_c (see manual)
!
!
!  specify a few more OPTIONAL parameters - normally not necessary
!
!
!     check_grad     = .true.    ! check the user supplied gradients
!     subsolver      =  1        ! specify the solver for the subproblems
      xtol           =  1.d-5    ! specify the x-tolerance (Euclidian norm)
      xtol_inf       =  1.d-7    ! specify the x-tolerance (infinity norm)
!     ftol           = -1.d-10   ! specify the f-tolerance. NOT recommended in general. Values < 0 disable this
!     outermax       =  500      ! the number of outer loops allowed
!     innermax       =  50       ! the number of inner loops allowed
!     dml_infinity   =  0.2d0    ! the infinity move limit
      DualTrustRadius=  1.d10    ! the allowed dual step size
      biglam         =  1.d10    ! upper bound on the dual variables
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
