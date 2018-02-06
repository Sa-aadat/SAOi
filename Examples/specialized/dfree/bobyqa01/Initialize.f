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
      double precision twopi, temp
      include          'ctrl_get.inc'
!
      twopi          =  8.0d0*datan(1.0d0)
      m              =  10                      !  problem specific parameter (5 or 10)
!
      n              =  2*m                     !  the number of design variables
      ni             =  0                       !  the number of inequality constraints
!
      do j=1,m                                  !  the starting point
        temp         =  dfloat(j)*twopi/dfloat(m)
        x(2*j-1)     =  dcos(temp)
        x(2*j)       =  dsin(temp)
      enddo
!
!     x(1)           =  1.000000D+00            !  the solution for m=10 (and hence n=20)
!     x(2)           =  1.000000D+00
!     x(3)           =  3.616079D-01
!     x(4)           =  1.000000D+00
!     x(5)           = -3.616080D-01
!     x(6)           =  1.000000D+00
!     x(7)           = -1.000000D+00
!     x(8)           =  1.000000D+00
!     x(9)           = -1.000000D+00
!     x(10)          = -7.834429D-08
!     x(11)          = -1.000000D+00
!     x(12)          = -1.000000D+00
!     x(13)          = -3.616080D-01
!     x(14)          = -1.000000D+00
!     x(15)          =  3.616079D-01
!     x(16)          = -1.000000D+00
!     x(17)          =  1.000000D+00
!     x(18)          = -1.000000D+00
!     x(19)          =  1.000000D+00
!     x(20)          = -4.504914D-08
!     fapriori       =  3.220305336883035D+01   !  the theoretical solution (for printing only, may be ignored)
!
      do i = 1,n
        x_lower(i)   = -1.d0                    !  the lower bounds
      enddo
!
      do i = 1,n
        x_upper(i)   =  1.d0                    !  the upper bounds
      enddo
!
      cname1 = 'BOBYQA-01'
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
