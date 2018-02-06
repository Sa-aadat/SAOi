! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: initialization
! SAOi: Fri 02 Dec 2016 07:54:04 SAST, Sa-aadat Parker UCT

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
      integer          i,j,n,m,ni,ne,nnz,nnzh
      double precision x(*), x_lower(*), x_upper(*), shift(*)
      double precision gc(nnz)
      double precision pi, dfactor 
      include          'ctrl_get.inc'
      
      pi =  2.d0*dacos(0.d0)
!       n                the number of design variables
!       ni               the number of inequality constraints
!       ne               the number of equality constraints

      if (n.gt.nmax) stop ' n > nmax - increase nmax in size.h'
      
      CALL FemInit(n,x)!Problem parameters
      
! angle upper and lower bounds
      do i = 1,n-ne
        x_lower(i)   =  -2.d0*pi!/-    0.d0!0.d0!  the lower bounds
      enddo
!
      do i = 1,n-ne
        x_upper(i)   = 2.d0*pi!/ !1.0d0   !/2.d0  the upper bounds
      enddo
!      
!       do j = 1,ne
!         eqn(j)       = .true.  !uncomment for equality constraints, but note that
!       enddo                      !that the equalities may be given in any order

!       do j = ne+1,ne+ni
!         eqn(j)         = .false. !uncomment for equality constraints, but note that
!       enddo                      !that the equalities may be given in any order
!
      cname1         = 'Topology NAND'! problem name; max 24 characters
!
!  specify REQUIRED global optimization parameters
!
      multimodal     = .false.        !  force a global search
      
!       tol_bayes      =  5.d-2
      ptarget        =  0.99d0
!
!  specify a few OPTIONAL parameters - sometimes desirable
!
      finite_diff    = .false.      ! use finite differences?
      force_converge =  2           ! enforce local convergence?
!       0 = No; 1 through 3 = Yes (see manual)
      approx_f       =  1           ! select approx_f (see manual)
      approx_c       =  1           ! select approx_c (see manual)

      subsolver      =  21           ! specify the solver for the subproblems
      check_grad     = .false.       ! check the user supplied gradients
      structure       =  1           ! invoke the sparse CSR storage scheme

      !  specify a few more OPTIONAL parameters - normally not necessary
!
!       xtol           =  1.d-5      ! specify x-tolerance(Euclidian norm)
!       xtol_inf       =  1.d-10      ! specify x-tolerance(Infinity norm)
      kkt_min        =  1.d1

      outermax       =  5000     ! the number of outer loops allowed
      itglobalmax    =  5000
!       innermax       =  1         ! the number of inner loops allowed
      dml_infinity   =  1.0d-2       ! the infinity move limit
!       biglam         =  1.d8        ! upper bound on the dual variables
!       GetFinDiffGradHessDense = .true.
! !       deltxH =  1.d-6!default = 1.d-6
      include          'ctrl_set.inc'
      return
      end subroutine SAOi_init
