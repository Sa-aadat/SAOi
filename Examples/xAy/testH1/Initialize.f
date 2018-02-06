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
      integer          i, j, k, n, m, p, ni, ne, nnz, nnzh
      integer          iscale1, iscale2
      double precision x(*), x_lower(*), x_upper(*), shift(*)
      include          'ctrl_get.inc'
!
      k              =  3         !  problem specific - number of requested eigen values
      p              =  5         !  problem specific - length of each eigen vector column
!
      n              =  k*p       !  the number of design variables
      ni             =  0         !  the number of inequality constraints
      ne             =  k*(k+1)/2 !  the number of equality constraints, i.e. the number of entries
!                                 !      in the upper triangle + the diagonal, or k*(k-1)/2 + n 
!
!    
!     structure      = 3                   ! sparse storage structure
!     nnz            = k*(k-1)/2*p*2 + k*p ! number of off-diagonal and diagonal constraint entries
!                                          ! sparsity = 1-nnz/(ne*n)
!
!     mode length       sparsity fraction of Jacobian
!
!           1                       0.000
!           2                       0.500
!           5                       0.800
!          10                       0.900 
!          20                       0.950 
!          30                       0.967 
!          50                       0.980 
!         100                       0.990 
!        1000                       0.999 
! 
!  The sparsity fraction remains constant irrespective of the number of modes k
!     requested for a given mode length p, since the Jacobian = [m x n], while 
!     the number of non-zero entries nnz = k*(k-1)/2*p*2 + k*p  
!
!
!
!
!
!
!
!
!
!
!
!
!  The number of non-zero lower-triangular Hessian entries due to the constraints is:
!
      QPform_H = 0
      nnzh = n*(n+1)/2
!     nnzh = k*(k-1)/2*p + k*p
!
!  This includes diagonal entries of the objective. Off diagonal entries of the objective
!     may be added if a diagonal approximation is not used for the objective. An approximation
!     for the objective seems attractive if the objective is not convex.
!
!
!
!
!
!
!
!
!
      iscale1        = 1   
      iscale2        = 1 
!
      do i = 1,ne
        eqn(i)=.true.
      enddo
!
      iuser(1)       =  k
      iuser(2)       =  p
      iuser(3)       =  iscale1
      iuser(4)       =  iscale2
!
      if (mayReturn) then
!      
        do i = 1,n
          x(i) = 0.5d0*rand()       !  the starting point
        enddo
!       
        do i = 1,n 
          x_lower(i)   =  -1.0d0    !  the lower bounds
          x_upper(i)   =   1.0d0    !  the upper bounds
        enddo
!        
      endif  
!
      mayReturn = .false.
!
      if (current_mode_num.eq.3) mayReturn = .true.
!
      cname1 = 'testH1'             ! problem name; max 24 characters
!
      GetFinDiffGradHessDense = .true.
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
!
!     multimodal     =  .true.     ! assume multimodality?
!     finite_diff    =  .true.     ! use finite differences?
      force_converge =  2          ! enforce local convergence? 0 = No; 1 through 3 = Yes (see manual)
!     approx_f       =  50         ! select approx_f (see manual)
!     approx_c       =  50         ! select approx_c (see manual)
!
!
!  specify a few more OPTIONAL parameters - normally not necessary
!
!
!     check_grad     = .true.     ! check the user supplied gradients
!     mayStop        = .false.
      subsolver      =  27        ! specify the solver for the subproblems
!     xtol           =  1.d-9     ! specify the x-tolerance (Euclidian norm)
!     xtol_inf       =  1.d-8     ! specify the x-tolerance (infinity norm)
!     ftol           = -1.d-10    ! specify the f-tolerance. NOT recommended in general. Values < 0 disable this
      outermax       =  500       ! the number of outer loops allowed
      innermax       =  50        ! the number of inner loops allowed
!     dml_infinity   =  0.2d0     ! the infinity move limit
!     random_start   =  .false.   ! use random starting point?
!     biglam         =  1.d8      ! upper bound on the dual variables
!
      include          'ctrl_set.inc'      
      return
      end subroutine SAOi_init
