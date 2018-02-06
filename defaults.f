! SAOi: 
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Default settings and version specific restrictions
! SAOi:


      subroutine SAOi_defaults (n, m, ni, ne, ictrl, lctrl, rctrl, 
     &                          cctrl, nnz, nnzh, eqn, lin, loop, 
     &                          iseed, shift, iuser, luser, cuser, 
     &                          ruser)
!----------------------------------------------------------------------!
!                                                                      !
!  do not edit this routine; the values specified here may be          !
!  overridden in routine initialize.f                                  !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          i, j, n, m, ni, ne, nq, iseed, nnz, nnzh, loop
      double precision drandu1, shift(*)
      !include          'ctrl_get.inc'
!
      if (loop.le.1) then
        iseed                  =  1
        drand0                 =  drandu1(iseed)
      endif
      n                        =  0
      m                        =  0
      ni                       =  0
      ne                       =  0
      n_lbfgs                  =  1  ! 1 = inactive / 2 = the current point + the last, etc !
      nnz                      = -1
      nnzh                     = -1
      approx_f                 =  1
      approx_c                 =  1
      force_converge           =  0
      subsolver                =  1
      outermax                 =  500
      innermax                 =  100
      ksubmax                  =  100000
      itglobalmax              =  10000
      iprint                   =  3
      nprob                    =  loop
      iaggressive              =  0
      structure                =  1
      ifalk                    =  0
      iwolfe                   =  0
      cplex_method             =  4
      cplex_feas               =  1
      alg_unconstr             =  1
      sph_dia_hess             =  1
      current_mode_num         =  0
      QPform_A                 =  2 ! CSR
      QPform_H                 =  3 ! 3 = DIA / 4 = COO
!
      convex                   =  .true.
      random_start             =  .false.
      finite_diff              =  .false.
      check_grad               =  .false.
      conservative             =  .false.
      trust_region             =  .false.
      unconstrained            =  .false.
      use_active_set           =  .false.
      allow_f                  =  .true.
      allow_c                  =  .true.
      classic_conserv          =  .true.
      classic_trust            =  .true.
      slack                    =  .false.
      relax                    =  .false.
      multimodal               =  .false.
      debug                    =  .false.
      special                  =  .false.
      calc_kkt                 =  .true.
      strict_struct            =  .true.
      mayStop                  =  .false.
      mayReturn                =  .true.
      GetFinDiffGradHessDense  =  .false.
      FormHessSparse           =  .false.
      ResetDualVars            =  .false.

!
      big                      =  1.d20
      dml_infinity             =  0.2d0
      xtol                     =  1.d-4
      ftol                     = -1.d-10     ! not used - negative !
      feaslim                  =  1.d-5
      deltx                    =  1.d-6
      deltxH                   =  1.d-6
      ci_slack                 =  0.d0
      di_slack                 =  0.d0
      actlim                   =  0.25
      tlimit                   =  100000.d0
      xtol_inf                 =  1.d-6
      rho_l_min                =  1.00d-8
      rho_l_max                =  1.d0
      biglam                   =  1.d8
      kkt_tol                  = -1.d-8      ! not used - negative ! 
      sfracmin                 =  0.3d0
      atol1                    =  1.d-5
      btol1                    =  0.d00
      fapriori                 =  big
      atol2                    =  5.d-5
      btol2                    =  5.d-5
      filter_hi                =  1.d10
      pen1                     =  0.d0
      pen2                     =  1.d0
      dummy1                   =  0.d0
      ymax                     =  1.d5  ! not to be confused with filter_hi !
      mma_lo                   =  0.7d0
      mma_hi                   =  1.2d0
      ptarget                  =  0.99d0
      tol_bayes                =  1.d-3
      rangemax                 =  1.d5
      ztol                     =  1.d-6
      kkt_min                  =  1.d-1
      DualTrustRadius          =  1.d5
      shift_tol                =  1.d-12
!
      do j=1,mmax
        eqn(j) = .false.
      enddo 
!
      if (loop.eq.1) then
          priv1 = 1.d0  
      else
        priv1 =  rctrl(35) * 1.1d0
        if (priv1.ge.1.d3) stop
      endif
!
!       do i=1,nictrl
!         ictrl = -1
!       enddo
!       do i=1,nrctrl
!         rctrl = -1.d0
!       enddo
!       do i=1,nlctrl
!         lctrl = .false.
!       enddo
!
      cname1 = '                        '
      cname2 = '++++++++++++++++++++++++'
!
      ictrl(nictrl)        = 123456789
      rctrl(nrctrl)        = 123456789.d0
      cctrl(ncctrl)        = '0a1b2c3d4e5f6g7h8i9j10kX'
      lctrl(nlctrl)        = .true.
      lctrl(nlctrl-1)      = .false.
      lctrl(nlctrl-2)      = .true.
!      
      iuser(imax)          = 123456789
      ruser(rmax)          = 123456789.d0
      cuser(cmax)          = '0a1b2c3d4e5f6g7h8i9j10kX'
      luser(lmax)          = .true.
      luser(lmax-1)        = .false.
      luser(lmax-2)        = .true.
!
      include 'ctrl_set.inc'
!
      return
      end subroutine SAOi_defaults
!----------------------------------------------------------------------!
      subroutine SAOi_assoc (ictrl, lctrl, rctrl, cctrl, 
     &                       m, n, lambda, z, iuser, luser, cuser, 
     &                       ruser)
!----------------------------------------------------------------------!
!                                                                      !
!  some sensible associated settings                                   !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      integer          i, j, m, n
      double precision lambda(m), z(n)
!
      include 'ctrl_get.inc'
!
      if (multimodal) iprint = 0
!
      if (ResetDualVars) then
        do j=1,m
          lambda(j) = 0.d0
        enddo  
        do i=1,n
          z(i) = 0.d0
        enddo  
      endif

!
      include 'ctrl_set.inc'
!
      return
      end subroutine SAOi_assoc
!----------------------------------------------------------------------!
      subroutine SAOi_assocM (ictrl, lctrl, rctrl, cctrl, sloop,
     &                        iuser, luser, cuser, ruser)
      implicit         none
      include          'ctrl.h'
      integer          sloop
!
      include 'ctrl_get.inc'
!
      current_mode_num = sloop ! for mode extraction
!
      include 'ctrl_set.inc'
!
      return
      end subroutine SAOi_assocM
!----------------------------------------------------------------------!
      subroutine SAOi_check (n, m, ni, ne, x, x_l, x_u, ictrl, lctrl, 
     &                       rctrl, cctrl, nnz, nnzh, eqn, lin, ierr,
     &                       iseed, shift, sloop, lambda, z,
     &                       iuser, luser, cuser, ruser)
!----------------------------------------------------------------------!
!                                                                      !
! check validity of the data for the current version of the            !
! algorithm. Do not edit this routine, it is version-specific          !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          i, j, n, m, ni, ne, ifc, nnz, nnzh, ierr, iseed
      integer          sloop
      double precision x(*), x_l(*), x_u(*), rdm
      double precision drandu1, shift(*), lambda(*), z(*) ! m temporarily not yet defined 
      equivalence      (ifc, force_converge)

      include 'ctrl_get.inc'
!
      m = ni+ne
!
      if (sloop.gt.mmax) stop ' sloop > mmax in SAOi_check'
!
      if (cname1.eq.'                        ') 
     &    cname1  = '--- (Not Specified...)  '

!
!  initialize the dual variables 
      if (sloop.eq.1.or.ResetDualVars) then
        do j=1,m
          lambda(j) = 0.d0
        enddo  
        do i=1,n
          z(i) = 0.d0
        enddo  
      endif
!
      ierr            = 0
      unconstrained   = .false.
      conservative    = .false.
      trust_region    = .false.
      classic_conserv = .true.
      classic_trust   = .true.
      relax           = .false.
!
      if (subsolver.ge.20.and.ifc.ne.0) ifc = 2 ! SQP only accepts trust
                                                ! region at this stage
!
      if (ifc.lt.0.or.ifc.gt.3) ifc = 0
!
      if (ifc.eq.1) conservative = .true.
!
      if (ifc.eq.2.or.ifc.eq.3) trust_region = .true.
!
      if (ifc.eq.3) classic_trust = .false.
!
      if (ifc.eq.0) innermax = 0
!
      pen2 = max (1.d-8,pen2)
      ymax = max (1.d-8,ymax)
!
      if (pen1.gt.0.d0) relax = .true.
!
      if (n.gt.nmax)  stop ': n > nmax. Increase nmax in size.h'
      if (n.le.0)               stop ': n <= 0'
      if (ni.lt.0)              stop ': ni < 0'
      if (ne.lt.0)              stop ': ne < 0'
      if (.not.convex)          stop ': not convex'
      if (xtol.le.0.0d0)        stop ': xtol <= 0.0d0'
      if (outermax.lt.0)        stop ': outermax < 0'
      if (innermax.lt.0)        stop ': innermax < 0'
      if (dml_infinity.le.0.d0) stop ': dml_inf <= 0.d0'
      if (feaslim.lt.0.d0)      stop ': feaslim < 0.d0'
!
      do j=1,mmax
        if (eqn(j)) then
          ifalk = 1         ! only effective for subsolver.eq.1
          !atol1 = -1.d3
          btol1 = -1.d3
          exit
        endif  
      enddo
!      
      if (m.eq.0) then 
        unconstrained=.true.
        subsolver = 0
      endif  
!
      ! kkt_tol = max(0.d0,kkt_tol)
!
      if (random_start) then
        do i=1,n
          rdm=drandu1(iseed)
          x(i) = (x_u(i)-x_l(i))*(rdm)+x_l(i)
        enddo
      endif
!
      write(6,*) ' '
      write(9,*) ' '
      do i=1,n
        if (x_l(i).eq.x_u(i)) then
          !write(6,1000) i
          write(9,1000) i
        endif
         if (x_l(i).gt.x_u(i)) then
          !write(6,1001) i
          write(9,1001) i
          stop
        endif
        if (x(i).lt.x_l(i)) then
          !write(6,1002) i
          write(9,1002) i
!         x(i) = (x_l(i)+x_u(i))/2.d0
          x(i) = x_l(i)
        endif
          if (x(i).gt.x_u(i)) then
          !write(6,1003) i
          write(9,1003) i
!         x(i) = (x_l(i)+x_u(i))/2.d0
          x(i) = x_u(i)
        endif
      enddo
!
      if (approx_f.eq.4.or.approx_c.eq.4.or.
     &    approx_f.eq.5.or.approx_c.eq.5.or.
     &    approx_f.eq.7.or.approx_c.eq.7) then 
        do i=1,n
          if (x_l(i).le.shift_tol) then
            approx_f = 8
            approx_c = 8
            !write(6,1004) i
            write(9,1004) i
            exit
          endif
        enddo
      endif
!
      if (approx_f.eq.8.or.approx_c.eq.8) then
        do i=1,n
          if (x_l(i).lt.shift_tol) then
            shift(i) = dabs(x_l(i)) + shift_tol
          else
            shift(i) = 0.d0
          endif
        enddo
      endif
!
      if (unconstrained) structure = 1
!      
      if (structure.eq.1.or.structure.eq.2) nnz = m*n
!
      if (nnzh.lt.0) nnzh = n ! assume diagonal
!
      !if (GetFinDiffGradHessDense) nnzh = n*(n+1)/2
!
      if (structure.eq.4) then
        strict_struct = .false.
        structure = 3
      endif
!
      if (structure.eq.3.and.m.eq.0) nnz = 0
!
      if (nnz.lt.0) then
        write(*,*)  ' Illegal number of entries in gc - try',
     &              ' compiling with -fdefault-integer-8, '
        write(*,*)  '    or check that nnz is correctly specified -',
     &              ' currently, (nnz,m,n) = (',nnz,',',m,',',n,')'
        write(*,*)  ' '
        stop
      endif
!
      if (structure.lt.3.and.GetFinDiffGradHessDense) then
        QPform_H = 0
        nnzh = n*(n+1)/2
        subsolver = 27
      endif
!
      include 'ctrl_set.inc'
!
      return
!
 1000 format (' WARNING : variable ',i9,': x_lower = x_upper')
 1001 format (' TERMINAL: variable ',i9,': x_lower > x_upper')
 1002 format (' WARNING : variable ',i9,': moved x_initial to x_lower')
 1003 format (' WARNING : variable ',i9,': moved x_initial to x_upper')
 1004 format (' WARNING : lower bound ',i9,' < e-12; set approx_* = 1')
!
      end subroutine SAOi_check
!----------------------------------------------------------------------!
