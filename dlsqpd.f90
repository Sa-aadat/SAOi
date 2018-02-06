! SAOi:
! SAOi: Wed Jul 14 14:26:25 SAST 2010, Salomon van Huyssteen, Stellenbosch
! SAOi: in collaboration with Albert Groenwold. Driver for
! SAOi: lsqp (double) to solve sparse QP-like DQ subproblems,
! SAOi: based on a Galahad example by the original authors
! SAOi:


! THIS VERSION: GALAHAD 2.2 - 23/04/2008 AT 16:30 GMT.

! double precision version !

   subroutine GALAHAD_LSQP (n, m, a_ne, f_a, c_l, c_u, x_l, x_u,  &
                            A_row, A_col, A_ptr, x, g, w, y, z,   &
                            A_val, string, ierr, ks, f)

   USE GALAHAD_LSQP_double
   IMPLICIT NONE
   INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )                   ! set precision
   REAL ( KIND = wp ), PARAMETER :: infinity = 10.0_wp ** 20
   TYPE ( QPT_problem_type ) :: p
   TYPE ( LSQP_data_type ) :: data
   TYPE ( LSQP_control_type ) :: control
   TYPE ( LSQP_inform_type ) :: inform
   INTEGER :: n, m, a_ne, ks
   REAL ( KIND = wp ) :: f, f_a, c_l(m), c_u(m), x_l(n), x_u(n)
   REAL ( KIND = wp ) :: x(n), y(m), z(n), A_val(a_ne)
   REAL ( KIND = wp ) :: g(n), w(n)
   INTEGER :: i, s, ierr
   INTEGER :: A_row(a_ne), A_col(a_ne), A_ptr(m+1)
   CHARACTER(len=*) :: string

! Initialize control parameters
   CALL LSQP_initialize( data, control, inform )                      

! start problem data
   ALLOCATE( p%X_l( n ), p%X_u( n ) )
   ALLOCATE( p%C( m ), p%C_l( m ), p%C_u( m ) )
   ALLOCATE( p%X( n ), p%Y( m ), p%Z( n ) )
   ALLOCATE( p%A%val( a_ne ), p%A%row( a_ne ), p%A%col( a_ne ),    &
             p%A%ptr(m+1))
             
   p%new_problem_structure = .TRUE.                           ! new structure 
   p%n = n                                                    ! primal unknowns
   p%m = m                                                    ! dual unknowns
   p%f = f                                                    ! function value
   p%C_l = C_l                                                ! constraint lower bounds
   p%C_u = C_u                                                ! constraint upper bounds
   p%X_l = X_l                                                ! variable lower bounds
   p%X_u = X_u                                                ! variable upper bounds

! sparse co-ordinate storage format: integer components
   CALL SMT_put( p%A%type, string, s )                        ! storage for A
   p%A%row = A_row
   p%A%col = A_col
   p%A%ptr = A_ptr
   p%A%ne  = a_ne
   p%A%val = A_val

! set initial primal and dual variables
   p%Y = Y
   p%Z = Z
   ALLOCATE( p%X0( n ) )
   p%X0 = 0.d0

! the gradient
   p%gradient_kind = 2
   ALLOCATE( p%g( n ) )
   p%g=g

! sqrt of the Hessian
   p%Hessian_kind = 2
   ALLOCATE( p%weight( n ) )
   p%weight = dsqrt(w)

! some accuracy stuff
   control%stop_c = 10.0_wp ** ( - 12 )
   control%itref_max = 2
   control%infinity = infinity
   control%restore_problem = 2                              
   control%stop_p = 10.0_wp ** ( - 12 )
   control%print_level = 0

! Solve problem
   CALL LSQP_solve( p, data, control, inform )

   X = p%X
   Y = p%Y
   Z = p%Z
   f_a = inform%obj
   ierr = inform%status
   ks = inform%iter

! delete internal workspace
   CALL LSQP_terminate( data, control, inform )

   END SUBROUTINE GALAHAD_LSQP
