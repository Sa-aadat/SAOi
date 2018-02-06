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
      integer          i, j, n, m, ni, ne, nnz, nnzh
      double precision x(*), x_lower(*), x_upper(*), shift(*)
      double precision eps, lo, zero
      data             zero /0.d0/
      include          'ctrl_get.inc'
!
      n              =  2        !  the number of design variables
      ni             =  2        !  the number of inequality constraints
!
      x(1)           =  3.0d0    !  the starting point
      x(2)           =  3.0d0
!
      lo             =  1.d-9    !  or zero
      eps            =  sqrt(lo) !  actually, we merely require eps > lo,
                                 !      although it is unclear how much is
                                 !      enough to `open up' the degenerated
                                 !      domain 
      ruser(1)       =  eps
!
      x_lower(1)     =  lo       !  the lower bounds
      x_lower(2)     =  lo
!
      x_upper(1)     =  8.0d0    !  the upper bounds
      x_upper(2)     =  8.0d0
!
      cname1 = 'RELAXED-4-BAR'   !  problem name; max 24 characters
!
      multimodal     =  .true.
      ptarget        =  0.9999d0
      force_converge =  1        !  enforce local convergence? 0 = No; 1 through 3 = Yes (see manual)
!
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
!
!     finite_diff    =  .true.  ! use finite differences?
!     force_converge =  0        ! enforce local convergence? 0 = No; 1 through 3 = Yes (see manual)
      approx_f       =  1        ! select approx_f (see manual)
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
!     random_start   =  .false.  ! use random starting point?
!     biglam         =  1.d8     ! upper bound on the dual variables
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
