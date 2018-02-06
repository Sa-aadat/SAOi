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
      integer          i, j, n, m, ni, ne, nnz, nseg, nnzh
      double precision x(*), x_lower(*), x_upper(*), shift(*)
      include          'ctrl_get.inc'
!
      nseg             = 100              ! problem specific parameter
!
      if (nseg.lt.1)     stop ' nseg < 1 '
!
      if (n.gt.nmax) stop ' n > nmax - increase nmax in size.h'
!
      n                = nseg*2           !  the number of design variables
      ni               = n                !  the number of inequality constraints
!
      do i=1,nseg
        x(2*i-1)       = 5.d0             !  the starting point
        x(2*i  )       = 60.d0
      enddo
!
      do i=1,nseg
        x_lower(2*i-1) = 1.d-9            !  the lower bounds
        x_lower(2*i  ) = 1.d-9
      end do
!
      do i=1,nseg
        x_upper(2*i-1) = 80.d0            !  the upper bounds
        x_upper(2*i  ) = 80.d0
      end do
!
      cname1 = 'VDPLAATS-CANT-NOTHING'
!
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
!
!     finite_diff    =  .false.  ! use finite differences?
!     force_converge =  0        ! enforce convergence? 0 = No; 1 through 3 = Yes (see manual)
      approx_f       =  4        ! select approx_f (see manual)
      approx_c       =  4        ! select approx_c (see manual)
!
!
!  specify a few more OPTIONAL parameters - normally not necessary
!
!
!     check_grad     = .true.    ! check the user supplied gradients
!     subsolver      =  1        ! specify the solver for the subproblems
!     xtol           =  1.d-4    ! specify the x-tolerance (Euclidian norm)
!     xtol_inf       =  1.d-6    ! specify the x-tolerance (infinity norm)
!     ftol           = -1.d-10   ! specify the f-tolerance. NOT recommended in general. Values < 0 disable this
!     outermax       =  500      ! the number of outer loops allowed
!     innermax       =  50       ! the number of inner loops allowed
!     dml_infinity   =  0.2d0    ! the infinity move limit
!     biglam         =  1.d8     ! upper bound on the dual variables
      structure      =  2        ! specify the storage scheme
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init