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
      include          'ctrl_get.inc'
!
      n              =  8                       !  the number of design variables (possible: 2,4,6,8)
      ni             =  0                       !  the number of inequality constraints
!
      do i = 1,n
        x(i)         =  dfloat(i)/dfloat(n+1)   !  the starting point
      enddo
!
!     x(1)           =  4.315266D-02            !  the solution for n=8
!     x(2)           =  1.930906D-01
!     x(3)           =  2.663288D-01
!     x(4)           =  4.999998D-01
!     x(5)           =  5.000001D-01
!     x(6)           =  7.336713D-01
!     x(7)           =  8.069091D-01
!     x(8)           =  9.568472D-01
!     fapriori       =  0.35168743E-02          !  the theoretical solution (for printing only, may be ignored)
!
      do i = 1,n
        x_lower(i)   = -10.d0                   !  the lower bounds (here selected arbitrarily)
      enddo
!
      do i = 1,n
        x_upper(i)   =  10.d0                   !  the upper bounds (here selected arbitrarily)
      enddo
!
      cname1 = 'NEWUOA-01'
!
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
      finite_diff      =  .true.                ! use finite differences?
      deltx            =  1.d-6                 ! finite difference step
      force_converge   =  1                     ! enforce convergence? 0 = No; 1 through 3 = Yes (see manual)
!     approx_f         =  1                     ! select approx_f (see manual)
!     approx_c         =  1                     ! select approx_c (see manual)
!     subsolver        =  20                    ! specify the solver for the subproblems
!     check_grad       = .true.                 ! check the user supplied gradients
!
!
!  specify a few more OPTIONAL parameters - normally not necessary
!
      xtol             =  1.d-5                 ! specify the x-tolerance (Euclidian norm)
      xtol_inf         =  1.d-7                 ! specify the x-tolerance (infinity norm)
!     ftol             = -1.d-10                ! specify the f-tolerance. NOT recommended in general. Values < 0 disable this
!     outermax         =  500                   ! the number of outer loops allowed
!     innermax         =  50                    ! the number of inner loops allowed
!     dml_infinity     =  0.2d0                 ! the infinity move limit
!     biglam           =  1.d8                  ! upper bound on the dual variables
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
