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
      integer          i, j, n, m, ni, ne, nnz, nnzh, nmat  
      integer          nstate, ground, prob
      double precision x(*), x_lower(*), x_upper(*), shift(*), du, p
      include          'ctrl_get.inc'
!
      call globalD1(1, n, nmat, nstate, m, nnz, du, ground, prob)
!
!     n                =             !  the number of design variables
      ni               = nmat        !  the number of equality constraints
      ne               = m-nmat      !  the number of inequality constraints
!
      do i=1,n
        if(i.le.nstate) then              ! state variables
          x(i)         =  0.0d0
          x_lower(i)   = -du
          x_upper(i)   =  du
        else
          if(ground.ne.3) x(i) = 0.5d0    ! design (density) variables 
          if(ground.eq.3) x(i) = 0.25d0
          x_lower(i)           = 0.d0
          x_upper(i)           = 1.d0
        endif
      enddo
!
      p = 3.d0
      ruser(1) = du
      ruser(2) = p      ! SIMP penalty
      iuser(1) = nnz
      iuser(2) = nmat
      iuser(3) = nstate
!
      cname1 = 'seasand'   ! problem name; max 24 characters
!      
      structure = 3
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
!
!     multimodal     =  .true.   ! assume multimodality?
!     finite_diff    =  .false.  ! use finite differences?
      force_converge =  1        ! enforce local convergence? 0 = No; 1 through 3 = Yes (see manual)
!     approx_f       =  1        ! select approx_f (see manual)
!     approx_c       =  1        ! select approx_c (see manual)
!
!
!  specify a few more OPTIONAL parameters - normally not necessary
!
!
!     check_grad     = .true.    ! check the user supplied gradients
!     subsolver      =  21        ! specify the solver for the subproblems
!     xtol           =  1.d-4    ! specify the x-tolerance (Euclidian norm)
!     xtol_inf       =  1.d-6    ! specify the x-tolerance (infinity norm)
!     ftol           = -1.d-10   ! specify the f-tolerance. NOT recommended. Values < 0 disable this
!     outermax       =  500      ! the number of outer loops allowed
!     innermax       =  50       ! the number of inner loops allowed
!     dml_infinity   =  0.2d0    ! the infinity move limit
!     random_start   =  .false.  ! use random starting point?
!     biglam         =  1.d8     ! upper bound on the dual variables
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
