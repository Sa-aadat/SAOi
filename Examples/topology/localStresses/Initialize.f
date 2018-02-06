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
      integer          i, j, n, ni, ne, nnz, nnzh
      double precision x(*), x_lower(*), x_upper(*), shift(*)
!
      logical          loaded_knee
      integer          nelx, nely, istrgrad, iadjoint, k, l
      integer          ind, icount,loop, ii, jj, meshmult, ifleury
      integer          isdense, jacobfilt, jblockwrite, mcompliance
      double precision penal, volfrac, rmin, p1, p2, penalq, ds
      double precision rho_min, sall, conlim, gcvol(n)
      double precision eps, eta, penal_bruns, penal2
      double precision gradlim, divval, theta, fscale
!
      include          'ctrl_get.inc'
!
      cname1 = 'TOPO-LOCAL-STRESS'   ! problem name; max 24 characters
!
      approx_f       = 4
      approx_c       = 4
!     
      force_converge = 0
!
      ifleury        = 0 
!
! ------------------------------------------------------------------------- 
! For ifleury:     0 = ignore and use approx_f and approx_c
!                  1 = CONLIN, approx_f and approx_c modified at bottom of page
!                  2 = Nonconvex reciprocal used for constraints. CONLIN for objective.
!                      approx_f and approx_c modified at bottom of page.
!     If a volume constraint is present, it is approximated as a direct linear function. 
! ------------------------------------------------------------------------- 
!
      xtol           =  5.d-2
      xtol_inf       =  1.d-3
      outermax       =  250
      innermax       =  50
      dml_infinity   =  0.2d0  
!
      subsolver      =  1 
!
      loaded_knee    =  .false.
      finite_diff    =  .false.        
!
! ------------------------------------------------------------------------- 
! normal compliance:     mcompliance    =  1 + istrgrad       =  0
! weight + stress:       mcompliance    =  0 + istrgrad       =  4
! compliance + stress:   mcompliance    =  1 + istrgrad       =  4
! ------------------------------------------------------------------------- 
!
      mcompliance    =  0              !  (0) Weight min topology, (1) min Compliance topology.
      sall           =  20.d0          !  Stress limit
      conlim         = -1.d0           !  Activation limit for relaxed stress constraint
      gradlim        =  1.d-4          !  User defined minimum absolute value for Jacobian terms
      divval         =  1.d3           !  Jacobian terms are saved if [JT.gt.(max|J_diag|/divval)]
      tlimit         =  150000000.d0
! 
! ------------------------------------------------------------------------- 
! Switches for FEM (see below):
! ------------------------------------------------------------------------- 
!
      istrgrad       =  4              !  Evaluate which stress and sensitivities  ?
      iadjoint       =  1              !  Calculate stress gradients using adjoint method (or not)
      isdense        =  0              !  Evaluate all constraints and return all Jacobian entries
      jacobfilt      =  0              !  Get rid of small entries in the jacobian (and method)
      jblockwrite    =  1              !  Blockwrite/blockread (1) - or not
! ------------------------------------------------------------------------- 
! Some protection and logic for switches:
! ------------------------------------------------------------------------- 
      if (istrgrad.gt.4) istrgrad = 4
!     if (istrgrad.eq.0) istrgrad = 1 ! remove when SAO properly updated ! 
      if (finite_diff) iadjoint = 0
      if ((.not.finite_diff).and.(istrgrad.eq.0).and.(mcompliance.eq.0)) 
     &     istrgrad = 4
      if (istrgrad.eq.0) iadjoint = 0
      if ((.not.finite_diff).and.(istrgrad.ne.0).and.(iadjoint.eq.0)) 
     &     iadjoint = 1
      if (iadjoint.gt.1) iadjoint = 1
      if (jacobfilt.ne.3) gradlim = 1.0d-16
! ------------------------------------------------------------------------------- !
! ----------- Note: The FEM will return the sensitivities of only the selected -- !
! -----------       stress measure (sigma_1, sigma_2, tau_12 or VonMises) ------- !
! -----------       Neither stresses nor stress sensitivities will be ----------- !
! -----------       calculated by the FEM if istrgrad = 0 ----------------------- !
! ------------------------------------------------------------------------------- !
!                                                                                 !
!      istrgrad  = ?    ! 1 = Sigma_1, 2 = Sigma_2, 3 = Shear, 4 = Von Mises.     !
!      iadjoint  = ?    ! 0 = Don't return stress gradients, 1 = Adjoint method.  !
!      isdense   = ?    ! 0 = Not dense. Select constraints according to 'conlim' !
!                       !     and record Jacobian terms according 'jacobfilt'.    !
!      jacobfilt = ?    ! 0 = Return all Jacobian entries for active constraints, !
!                       ! 1 = Filter, 'gradlim' calculated from iter(k-1) info,   !
!                       ! 2 = Filter, and 'gradlim' updates as you go,            !
!                       ! 3 = Filter using user defined 'gradlim'.                !  
!                                                                                 !
! ------------------------------------------------------------------------------- !
! ------------------------------------------------------------------------------- !
!
! the mbb beam                    !  problem related thingies 
!       meshmult=1                !               "
!       nelx=5*3*meshmult         !               "
!       nely=5*meshmult           !               "
!       rmin=1.d0*dble(meshmult)  !               "
!       volfrac=.5d0              !               "              
!
! the twobar problem              !  problem related thingies 
      meshmult=2                  !               "
      nelx=5*3*meshmult           !  NOTE: nelx must be integer devisible by 5 for auto meshgen
      nely=5*meshmult             !               "
      rmin=1.d0                   !               "
      volfrac=0.5d0               !               "
!
! the flat plate                  !  problem related thingies 
!      meshmult=4                 !               "
!      nelx=25*meshmult           !               "
!      nely=25*meshmult           !               "
!      rmin=2.d0*dble(meshmult)   !               "
!      volfrac=.5d0               !               "
! 
! the loaded knee general
!      loaded_knee=.true.         !               "
!      meshmult=1                 !               "
!      nelx=60*meshmult           !               "
!      nely=60*meshmult           !               "
!      rmin=2.d0*dble(meshmult)   !               "
!      volfrac=0.75d0*13.d0/27.d0 !               "
! 
! the loaded knee for Svanberg's 96x96 nhs example 
!      loaded_knee=.true.         !               "
!      meshmult=1                 !               "
!      nelx=96*meshmult           !               "
!      nely=96*meshmult           !               "
!      rmin=1.4d0*dble(meshmult)  !               "
!      volfrac=0.75d0*3064.0d0/(6912.0d0) ! = 3064/(nelx*nely)
!
      fscale           =  1.d0
      theta            =  1.0d0
      eps              =  0.1d0   ! 1.d-2
      rho_min          =  1.d-4   ! Lower bound on density
      penal            =  3.d0
      penal2           =  1.d0  
      penalq           =  1.d0  
      penal_bruns      =  1.d0  
      n                =  nelx*nely
      ni               =  n                 ! For stress constrained weight min problem
      if (mcompliance.eq.1) then
         if (istrgrad.eq.0) then
            ni = 1                          ! For min compliance problem with only a volume constraint 
         else
            ni = n + 1                      ! For stress constrained min compliance problem 
         endif
      endif
!
      if (ifleury.eq.1) then                ! To select CONLIN
         approx_f = 101
         approx_c = 101
      endif
      if (ifleury.eq.2) then                ! To select non-convex reciprocal
         approx_f = 102
         approx_c = 102
      endif
!
      DualTrustRadius  =  1.d6
      feaslim          =  1.d-5
      special          =  .true. 
      strict_struct    =  .false.
      nnz              =  ni*n     ! the max number of non-zero terms in gc
      structure        =  3        ! invoke the sparse compressed sparse row (CSR) storage scheme
      iprint           =  2        ! a bit less output
      atol1            =  1.d-6
      btol1            =  0.d0
      biglam           =  1.d8
      pen1             =  0.d4
!
      if (structure.le.2) isdense = 1
!
      if (mcompliance.eq.1.and.structure.ne.3) then
        stop ' Oeps, something went wrong, ask DWW/AG what''s up'
      endif
!
      call system('rm -f DesCons')
      call system('rm -f DesVars')
      open (22,file='DesCons',status='unknown',err=1005)
      write(22,*) n
      write(22,*) nelx
      write(22,*) nely
      close(22)
!
      random_start=.false.
      do i=1,n
        x(i)=volfrac
      enddo
!
      do i=1,n
        x_lower(i)= rho_min  
      end do
!
      do i=1,n
        x_upper(i)=1.d0
      end do
!
!       if (loaded_knee) then
!         if (mod(nelx,2).ne.0) stop 'mod(nelx,2).ne.0'
!         ns=(nelx*nely)/2
!         do j=ns,nelx*nely-nelx,nelx
!           do i=1,nely/2
!             x_upper(j+i)= rho_min + 1.d-10
!           enddo
!         enddo
!
!         ns=nelx*nely-(nely-1)
!         ns1=nely/(2*3)
!
!         do j=1,ns1
!           x_lower(ns + nely*2/3 + j-1)=1.d0-1.d-10 
!         enddo
!       endif
!
      write(*,*) '  '
      write(*,*) ' topological settings: '
      write(*,*) '    rho_min       : ',rho_min
      write(*,*) '    nelx          : ',nelx
      write(*,*) '    nely          : ',nely
      write(*,*) '    n             : ',n
      write(*,*) '    ni            : ',ni
      write(*,*) '    penal         : ',penal
      write(*,*) '    volfrac       : ',volfrac
      write(*,*) '    rmin          : ',rmin
      write(*,*) '    dml_infinity  : ',dml_infinity
      write(*,*) '    random_start  : ',random_start
      write(*,*) '    xtol          : ',xtol
      write(*,*) '    sall          : ',sall
      write(*,*) '    conlim        : ',conlim
      write(*,*) '    tlimit        : ',tlimit
      write(*,*) '  '
!
      ruser(1)  = rho_min
      ruser(2)  = penal
      ruser(3)  = volfrac
      ruser(4)  = rmin
      ruser(5)  = p1
      ruser(6)  = p2
      ruser(7)  = penalq
      ruser(8)  = penal_bruns
      ruser(9)  = penal2
      ruser(10) = sall
      ruser(11) = conlim
      ruser(12) = gradlim
      ruser(13) = divval
      ruser(14) = theta
      ruser(15) = eps
      ruser(16) = fscale
      iuser(1)  = nelx
      iuser(2)  = nely
      iuser(3)  = istrgrad
      iuser(4)  = iadjoint
      iuser(5)  = isdense
      iuser(6)  = jacobfilt
      iuser(7)  = jblockwrite
      iuser(8)  = mcompliance
      luser(1)  = loaded_knee
!
      include          'ctrl_set.inc'
!      
      return
!
 1005 stop ' fatal i/o error in Initialize.f. period. '
!
      end subroutine SAOi_init
