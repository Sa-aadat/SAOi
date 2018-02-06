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
      implicit          none
      include           'ctrl.h'
      logical           eqn(*), lin(*)
      integer           i, j, k, n, m, p, ni, ne, n2, nnz, nnzh
      integer           iscale1, iscale2
      double precision  x(*), x_lower(*), x_upper(*), shift(*)
      !double precision, dimension(:), allocatable :: xtemp
      double precision  xtemp(nmax)
      
      include          'ctrl_get.inc'
!
      k              =  current_mode_num   !  problem specific - number of requested eigen values
      p              =  5                  !  problem specific - length of each eigen vector column
!
      n              =  k*p                !  the number of design variables
      ni             =  0                  !  the number of inequality constraints
      ne             =  k*(k+1)/2          !  the number of equality constraints, i.e. the number of
!                                          !      entries in the upper triangle + the diagonal, 
!                                          !      or k*(k-1)/2 + n 
!
!    
      structure      = 3                   ! sparse storage structure
      nnz            = k*(k-1)/2*p*2 + k*p ! number of off-diagonal and diagonal constraint entries
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
      iscale1        = 1   
      iscale2        = 1  
! 
      n2 = n*(n+1)/2
!
!
!  determine nnzh - coded in this file for each specific problem
!
      ResetDualVars = .true.
      FormHessSparse = .true.
      QPform_H = 4
!
      call Get_nnzh (n, ni+ne, ni, ne, n2, k, p, iscale1, iscale2,
     &               ictrl, lctrl, rctrl, cctrl, iuser, luser, 
     &               cuser, ruser, nnz, nnzh, eqn, lin)
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
      if ((k-1)*p.gt.0) then
!      
        !allocate ( xtemp(n) )
        do i = 1,n
          xtemp(i) = x(i) 
        enddo
!        
        do i=n,p+1,-1 
          x(i) = xtemp(i-p)
          x_lower(i) = -1.0d0
          x_upper(i) =  1.0d0  
        enddo
!
        do i = 1,p
          x(i) = 2.d0*rand() - 1.d0
          x_lower(i) = -1.0d0 
          x_upper(i) =  1.0d0   
        enddo
!
      else
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
      
      !stop
      
!
      !deallocate (xtemp)
!
      mayReturn = .false.
!
      if (k.eq.p) mayReturn = .true.
!
      cname1 = 'testH7'          ! problem name; max 24 characters
!
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
!
!     multimodal     =  .true.    ! assume multimodality?
!     finite_diff    =  .true.    ! use finite differences?
      force_converge =  2         ! enforce local convergence? 0 = No; 1 through 3 = Yes (see manual)
      approx_f       =  100       ! select approx_f (see manual)
      approx_c       =  100       ! select approx_c (see manual)
!
!
!  specify a few more OPTIONAL parameters - normally not necessary
!
!
!     check_grad     = .true.     ! check the user supplied gradients
!     mayStop        = .false.
      subsolver      =  27        ! specify the solver for the subproblems
      xtol           =  1.d-6     ! specify the x-tolerance (Euclidian norm)
      xtol_inf       =  1.d-7     ! specify the x-tolerance (infinity norm)
!     ftol           = -1.d-10    ! specify the f-tolerance. NOT recommended in general. Values < 0 disable this
      outermax       =  100       ! the number of outer loops allowed
!     innermax       =  50        ! the number of inner loops allowed
!     dml_infinity   =  0.2d0     ! the infinity move limit
!     random_start   =  .false.   ! use random starting point?
!     biglam         =  1.d8      ! upper bound on the dual variables
!
      include          'ctrl_set.inc'      
      return
      end subroutine SAOi_init
!---------------------------------------------------------------------------------------------!
!                                                                       
!  Determine the number of Hessian entries.                             
!                                                                       
!
!  The number of non-zero lower-triangular Hessian entries due to the constraints is:
!     
!      nnzh = k*(k-1)/2*p + k*p
!
!  This includes diagonal entries of the objective. Off diagonal entries of the objective
!     may be added, as is done below, if a diagonal approximation is not used for the 
!     objective. An approximation for the objective seems attractive if the objective 
!     is not convex.
!
!  Note that duplicate entries for the Galahad solvers are allowed; the are summed in the   
!     solvers qpa, qpb or qpc. (lsqp should be used if everything is approximated diagonal.
!
!---------------------------------------------------------------------------------------------!
      subroutine Get_nnzh (n, m, ni, ne, n2, k, p, iscale1, iscale2,
     &                     ictrl, lctrl, rctrl, cctrl, iuser, luser, 
     &                     cuser, ruser, nnz, nnzh, eqn, lin)
      implicit          none
      include           'ctrl.h'
      logical           eqn(*), lin(*), Hrowflag(n2), Hcolflag(n2)
      integer           i, j, k, n, m, n2, p, ni, ne, nnz, nnzh
      integer           a, d, i1, k1, l, allocatestatus 
      integer           iscale1, iscale2, nnzh_off, nnzh_dia
      double precision  temp
      double precision, dimension(:),    allocatable :: xs
      double precision, dimension(:, :), allocatable :: toeplitz
      include          'ctrl_get.inc'
!
      if (k.gt.p) stop ' k.gt.p in Get_nnzh: not physically possible'
!
      allocate ( xs(p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate xs"
!
      allocate ( toeplitz(p, p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate toeplitz"
!
      xs    = 0.d0
      xs(1) = 4.d0
      xs(2) = 1.d0
!
      do i = 1,p
        toeplitz(i,i:p) = xs(1:(p+1-i))
        toeplitz(i:p,i) = xs(1:(p+1-i))
      end do  
!
      nnzh = k*(k-1)/2*p + k*p
!
      nnzh_off = 0 ! not yet used
      nnzh_dia = 0 ! not yet used 
!
      i1 = 0
      do l=1,k
        do a=1,p 
          do j=1,a
            temp = toeplitz(a,j)
            if (temp.ne.0.d0) then
!            
              nnzh = nnzh + 1
              if (a+(l-1)*p.eq.j+(l-1)*p) then
                nnzh_dia = nnzh_dia + 1
              else
                nnzh_off = nnzh_off + 1
              endif  
!
            endif
          enddo
        enddo  
      enddo      
!      
      return
      end subroutine Get_nnzh
!----------------------------------------------------------------------!
