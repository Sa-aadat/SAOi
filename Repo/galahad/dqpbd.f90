! SAOi:
! SAOi: Wed Jul 14 14:26:25 SAST 2010, Albert Groenwold, Stellenbosch
! SAOi: Driver for qpb (double) to solve sparse QP-like DQ subproblems
! SAOi:


! THIS VERSION: GALAHAD 2.2 - 23/04/2008 AT 16:30 GMT.

! double precision version !

   subroutine GALAHAD_QPB (n, m, a_ne, f_a, c_l, c_u, x_l, x_u,  &
                           A_row, A_col, A_ptr, x, g, w, y, z,   &
                           A_val, A_string, ierr, ks, H_string,  &
                           H_row, H_col, H_ptr, H_val, h_ne, f)

   USE GALAHAD_QPB_double
   IMPLICIT NONE
   INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )              ! set precision 
   REAL ( KIND = wp ), PARAMETER :: infinity = 10.0_wp ** 20
   TYPE ( QPT_problem_type ) :: p
   TYPE ( QPB_data_type ) :: data
   TYPE ( QPB_control_type ) :: control
   TYPE ( QPB_inform_type ) :: inform
   INTEGER :: n, m, h_ne, a_ne, ks
   INTEGER :: i, s, ierr
   INTEGER :: A_row(a_ne), A_col(a_ne), A_ptr(m+1)
   INTEGER :: H_row(h_ne), H_col(h_ne), H_ptr(n+1)
   REAL ( KIND = wp ) :: f, f_a, c_l(m), c_u(m), x_l(n), x_u(n)
   REAL ( KIND = wp ) :: x(n), y(m), z(n), A_val(a_ne)
   REAL ( KIND = wp ) :: g(n), w(n), H_val(h_ne)
   CHARACTER(len=*) :: A_string, H_string

! Initialize control parameters
   CALL QPB_initialize( data, control, inform )

! start problem data
   ALLOCATE( p%G( n ), p%X_l( n ), p%X_u( n ) )
   ALLOCATE( p%C( m ), p%C_l( m ), p%C_u( m ) )
   ALLOCATE( p%X( n ), p%Y( m ), p%Z( n ) )
   p%new_problem_structure = .TRUE.                         ! new structure
   p%n = n                                                  ! dimensions & objective constant
   p%m = m   
   p%f = f                               
   p%G = g                                                  ! objective gradient
   p%C_l = c_l                                              ! constraint lower bounds
   p%C_u = c_u                                              ! constraint upper bounds
   p%X_l = x_l                                              ! variable lower bounds
   p%X_u = x_u                                              ! variable upper bounds
   !p%rho_g = 1.0_wp ; p%rho_b = 1.0_wp                      ! initial penalty parameters
   p%X = 0.0_wp                                             ! start from zero
   p%Y = Y
   p%Z = Z

   ! storage format
   CALL SMT_put( p%H%type, H_string, s )                    ! storage for H
   
   if (H_string.eq.'DIAGONAL') then
   
     ALLOCATE( p%H%val( n ) )
     p%H%val = H_val
     
   else
   
     ALLOCATE( p%H%val( h_ne ), p%H%row( h_ne ),   &
               p%H%col( h_ne ), p%H%ptr( n+1 ) )
               
     p%H%row = H_row
     p%H%col = H_col
     p%H%ptr = H_ptr
     p%H%ne  = h_ne
     p%H%val = H_val
     
   endif

   CALL SMT_put( p%A%type, A_string, s )                    ! storage for A
   
   ALLOCATE( p%A%val( a_ne ), p%A%row( a_ne ),   &
             p%A%col( a_ne ), p%A%ptr( m+1 ) )
             
   p%A%row = A_row
   p%A%col = A_col
   p%A%ptr = A_ptr
   p%A%ne  = a_ne
   p%A%val = A_val

! some accuracy stuff
   control%stop_c = 10.0_wp ** ( - 12 )
   control%itref_max = 2
   control%infinity = infinity
   control%restore_problem = 2
   control%stop_p = 10.0_wp ** ( - 12 )
   control%print_level = 0

   CALL QPB_solve( p,  data,  control, inform )              ! Solve problem

   X = p%X
   Y = p%Y
   Z = p%Z
   f_a = inform%obj
   ierr = inform%status
   ks = inform%iter

   CALL QPB_terminate( data, control, inform )                ! delete internal workspace

   END SUBROUTINE GALAHAD_QPB
