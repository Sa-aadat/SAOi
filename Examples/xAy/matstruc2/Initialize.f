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
      double precision x(*), x_lower(*), x_upper(*), shift(*)
      include          'ctrl_get.inc'
!
!       do i=1,1000
!         write(*,*)i,i*i,i*(i+1)/2,dble(i*(i+1)/2)/(0.5d0*dble(i*i))
!       enddo
!       stop
!
      k              =  5         !  problem specific - number of requested eigen values
      p              =  5         !  problem specific - length of each eigen vector column
!
      n              =  k*p       !  the number of design variables
      ni             =  0         !  the number of inequality constraints
      ne             =  k*(k+1)/2 !  the number of equality constraints, i.e. the number of entries
!                                 !      in the upper triangle + the diagonal, i.e. k*(k-1)/2 + n      
!
      do i = 1,ne
        eqn(i)=.true.
      enddo
!
      iuser(1)       =  k
      iuser(2)       =  p
!
      do i = 1,n
        x(i) = 0.51d0*rand()      !  the starting point
      enddo
!       
!
!        x(1)=    -0.2768           
!        x(2)=     0.5565           
!        x(3)=    -0.4769           
!        x(4)=     0.5565           
!        x(5)=    -0.2768           
!        x(6)=     0.4068           
!        x(7)=     0.4340           
!        x(8)=     0.5406           
!        x(9)=     0.4340           
!        x(10)=    0.4068           
!                                  
!

!       x(1)=    -0.2733 
!       x(2)=     0.5153 
!       x(3)=    -0.5654 
!       x(4)=     0.5153 
!       x(5)=    -0.2733 
!       x(6)=     0.5684 
!       x(7)=    -0.4206 
!       x(8)=    -0.0000 
!       x(9)=     0.4206 
!       x(10)=   -0.5684 
!       x(11)=   -0.5065 
!       x(12)=    0.1053 
!       x(13)=    0.6817 
!       x(14)=    0.1053 
!       x(15)=   -0.5065 
!       x(16)=   -0.4206 
!       x(17)=   -0.5684 
!       x(18)=    0.0000 
!       x(19)=    0.5684 
!       x(20)=    0.4206 
!       x(21)=   -0.4108
!       x(22)=   -0.4727
!       x(23)=   -0.4644
!       x(24)=   -0.4727
!       x(25)=   -0.4108


      do i = 1,n 
        x_lower(i)   =  -1.0d0    !  the lower bounds
        x_upper(i)   =   1.0d0    !  the upper bounds
      enddo 
!
      cname1 = 'matstruc2'        ! problem name; max 24 characters
!
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
!
!     multimodal     =  .true.    ! assume multimodality?
      finite_diff    =  .true.    ! use finite differences?
      force_converge =  2         ! enforce local convergence? 0 = No; 1 through 3 = Yes (see manual)
!     approx_f       =  50         ! select approx_f (see manual)
!     approx_c       =  50         ! select approx_c (see manual)
!
!
!  specify a few more OPTIONAL parameters - normally not necessary
!
!
!     check_grad     = .true.     ! check the user supplied gradients
!     mayStop        = .false.
      subsolver      =  25        ! specify the solver for the subproblems
!     xtol           =  1.d-5     ! specify the x-tolerance (Euclidian norm)
!     xtol_inf       =  1.d-7     ! specify the x-tolerance (infinity norm)
!     ftol           = -1.d-10    ! specify the f-tolerance. NOT recommended in general. Values < 0 disable this
!     outermax       =  500       ! the number of outer loops allowed
!     innermax       =  50        ! the number of inner loops allowed
!     dml_infinity   =  0.2d0     ! the infinity move limit
!     random_start   =  .false.   ! use random starting point?
!     biglam         =  1.d8      ! upper bound on the dual variables
!
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
