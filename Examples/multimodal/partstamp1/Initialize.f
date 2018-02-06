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
      integer          i, j, n, m, ni, ne, nnz, nnzh, ndisk, isum
      double precision x(*), x_lower(*), x_upper(*), shift(*)
      include          'ctrl_get.inc'
!
      ndisk            =  5
!
      n                =  2*ndisk+2         !  the number of design variables
!
      isum = 0
      do i=1,ndisk-1
        isum = isum+i
      enddo
!
      ni               =  2*ndisk + isum    !  the number of inequality constraints
!
      do i=1,ndisk
        x(i)           =  dble(i)     !  the starting point
        x(i+ndisk)     =  dble(i)
      enddo
      x(n-1)           = 2.d0*dble(ndisk)
      x(n)             = 2.d0*dble(ndisk)
!
      do i=1,n
        x_lower(i)     =  1.d0     !  lower bounds
      enddo
!
      do i=1,n
        x_upper(i)     =  2.d0*dble(ndisk)      !  upper bounds
      enddo
!
      cname1 = 'PARTSTAMP-1'
!
!
!  specify REQUIRED global optimization parameters
!
!
      multimodal       = .true.    !  force a global search
!
!
!  specify a few OPTIONAL global optimization parameters - often desirable
!
      force_converge   = 1                      !  enforce local convergence? 0 = No; 1 through 3 = Yes (see manual)
      xtol             = 1.d-4                  !  tighten the convergence tolerance
!     fapriori         = .........              !  the theoretical solution (for printing only, may be ignored)
!     itglobalmax      = 50                     !  max number of random restarts
      ptarget          = 0.99d0                 !  the desired probability of convergence to the global optimum
!     tol_bayes        = 1.d-4                  !  assume local minima within tol_bayes to be equal
!
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
!     finite_diff      =  .true.   ! use finite differences?
!     check_grad       =  .true.   ! check the user supplied gradients
!     approx_f         =  1        ! select approx_f (see manual)
!     approx_c         =  1        ! select approx_c (see manual)
!     subsolver        =  1        ! specify the solver for the subproblems
!
!
!  specify a few more OPTIONAL parameters - normally not necessary
!
!     xtol_inf         =  1.d-6    ! specify the x-tolerance (infinity norm)
!     ftol             = -1.d-10   ! specify the f-tolerance. NOT recommended in general. Values < 0 disable this
!     outermax         =  500      ! the number of outer loops allowed
!     innermax         =  50       ! the number of inner loops allowed
!     dml_infinity     =  0.2d0    ! the infinity move limit
      biglam           =  1.d1     ! upper bound on the dual variables
                                   ! for this problem: very badly scaled, dual vars less than 1...
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
