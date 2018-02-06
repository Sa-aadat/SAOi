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
      double precision pi, d_theta, R_v, R_min, R_max, alpha_c
      include          'ctrl_get.inc'
!
      n                =  15          !  the number of design variables
      ni               =  3*n+4       !  the number of inequality constraints
!
      pi               =  4.0d0*datan(1.0d0)           !  tproblem specific stuff
      d_theta          =  2.d0*pi/(5.d0*dble(n+1))
      R_v              =  1.0d0
      R_min            =  1.0d0
      R_max            =  2.0d0
      alpha_c          =  1.5d0
!
      if (n.gt.nmax) stop ' n > nmax - increase nmax in size.h'
!
      do i=1,n
        x(i)           =  (R_min + R_max)/2.0d0    !  the starting point
      enddo
!
      do i=1,n
        x_lower(i)     =  R_min       !  the lower bounds
      end do
!
      do i=1,n
        x_upper(i)     =  R_max       !  the upper bounds
      end do
!
      cname1 = 'CAM-DESIGN-1'
!
      structure        =  2           !  invoke the sparse compressed sparse row (CSR) storage scheme
!
!
!  specify a few OPTIONAL parameters - sometimes desirable (not all are available in sparse)
!
!
!     finite_diff      = .false.      !  use finite differences?
!     force_converge   =  0           !  enforce convergence? 0 = No; 1 through 3 = Yes (see manual)
!     approx_f         =  4           !  select approx_f (see manual)
!     approx_c         =  4           !  select approx_c (see manual)
!
!
!  specify a few more OPTIONAL parameters - normally not necessary (not all are available in sparse)
!
!
!     check_grad       = .true.       !  check the user supplied gradients
      subsolver        =  25          !  specify the subproblem solver - beware of dual statements for n large ! 
!     xtol             =  1.d-5       !  specify the x-tolerance (Euclidian norm)
!     xtol_inf         =  1.d-7       !  specify the x-tolerance (infinity norm)
!     ftol             = -1.d-10      !  specify the f-tolerance. NOT recommended in general. Values < 0 disable this
!     outermax         =  500         !  the number of outer loops allowed
!     innermax         =  50          !  the number of inner loops allowed
!     dml_infinity     =  0.2d0       !  the infinity move limit
!     random_start     =  .false.     !  use random starting point?
!     biglam           =  1.d8        !  upper bound on the dual variables
!     atol1            =  0.d-6
!     btol1            =  0.d-6
!     pen1             =  1.d3
!
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
