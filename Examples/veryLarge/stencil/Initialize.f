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
      n                =  0           !  the number of design variables
      ni               =  0           !  the number of inequality constraints
      ne               =  0           !  the number of equality constraints 
!
      if (n.gt.nmax) stop ' n > nmax - increase nmax in size.h'
!
      do i=1,n
        x(i)           =  0.d0        !  the starting point
      enddo
!
      do i=1,n
        x_lower(i)     =  0.d0        !  the lower bounds
      end do
!
      do i=1,n
        x_upper(i)     =  0.d0        !  the upper bounds
      end do
!
      nnz              =  0           !  the number of non-zero terms in gc
      structure        =  3           !  invoke the sparse compressed sparse row (CSR) storage scheme
!
      cname1 = 'STENCIL'
!
!
!  specify a few OPTIONAL parameters - sometimes desirable (not all are available in sparse)
!
!
      finite_diff      = .true.      !  use finite differences?
!     force_converge   =  0           !  enforce convergence? 0 = No; 1 through 3 = Yes (see manual)
!     approx_f         =  4           !  select approx_f (see manual)
!     approx_c         =  4           !  select approx_c (see manual)
!
!
!  specify a few more OPTIONAL parameters - normally not necessary (not all are available in sparse)
!
!
!     check_grad       = .true.       !  check the user supplied gradients
!     subsolver        =  1           !  specify the subproblem solver
!     xtol             =  1.d-3       !  specify the x-tolerance (Euclidian norm)
!     xtol_inf         =  1.d-6       !  specify the x-tolerance (infinity norm)
!     ftol             = -1.d-10      !  specify the f-tolerance. NOT recommended in general. Values < 0 disable this
!     outermax         =  500         !  the number of outer loops allowed
!     innermax         =  50          !  the number of inner loops allowed
!     dml_infinity     =  0.2d0       !  the infinity move limit
!     random_start     =  .false.     !  use random starting point?
!     biglam           =  1.d8        !  upper bound on the dual variables
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
