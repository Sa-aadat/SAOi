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
      integer          i, j, n, m, ni, ne, nnz, nnzh, nseg
      double precision x(*), x_lower(*), x_upper(*), shift(*), tmp
      include          'ctrl_get.inc'
!
      nseg             =  1280       !  number of segments - problem specific parameter
      n                =  nseg*2     !  the number of design variables
      ni               =  n          !  the number of inequality constraints
!
      if (n.gt.nmax) stop ' n > nmax - increase nmax in size.h'
!
      do i=1,nseg
        x(2*i-1)       =  5.d0       !  the starting point
        x(2*i  )       =  60.d0
      enddo
!
      do i=1,nseg
        x_lower(2*i-1) = 1.d0        !  the lower bounds
        x_lower(2*i  ) = 5.d0
      end do
!
      do i=1,nseg
        x_upper(2*i-1) =  80.d0      !  the upper bounds
        x_upper(2*i  ) =  80.d0
      end do
!
!
      nnz              =  n*2        !  the number of non-zero terms in gc
      structure        =  3          !  invoke the sparse compressed sparse row (CSR) storage scheme
!
      cname1 = 'VDPLTS-CANT-NO-DISP-VLSO'
!

!  specify a few OPTIONAL parameters - sometimes desirable
!
!
!     finite_diff      = .false.     !   use finite differences?
!     force_converge   =  0          !   enforce convergence? 0 = No; 1 through 3 = Yes (see manual)
      approx_f         =  4          !   select approx_f (see manual)
      approx_c         =  4          !   select approx_c (see manual)
!
!
!  specify a few more OPTIONAL parameters - normally not necessary
!
!
!     check_grad       = .true.      !   check the user supplied gradients
!     subsolver        =  1
      xtol             =  1.d-3      !   specify the x-tolerance (Euclidian norm)
!     xtol_inf         =  1.d-6      !   specify the x-tolerance (infinity norm)
!     ftol             = -1.d-10     !   specify the f-tolerance. NOT recommended in general. Values < 0 disable this
!     outermax         =  500        !   the number of outer loops allowed
!     innermax         =  50         !   the number of inner loops allowed
!     dml_infinity     =  0.2d0      !   the infinity move limit
!     random_start     =  .false.    !   use random starting point?
!     biglam           =  1.d8       !   upper bound on the dual variables
!     tlimit           =  1000.d0    !   a time limit for the dual
      iprint           =  2
      atol1            =  0.d0
      btol1            =  0.d0
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
