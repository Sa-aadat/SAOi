! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Include file to overcome the need for common blocks present
! SAOi: in earlier versions.
! SAOi:
!


      approx_f                   =  ictrl(1)
      approx_c                   =  ictrl(2)
      force_converge             =  ictrl(3)
      subsolver                  =  ictrl(4)
      outermax                   =  ictrl(5)
      innermax                   =  ictrl(6)
      ksubmax                    =  ictrl(7)
      itglobalmax                =  ictrl(8)
      iprint                     =  ictrl(9)
      iaggressive                =  ictrl(10)
      structure                  =  ictrl(11)
      nprob                      =  ictrl(12)
      ifalk                      =  ictrl(13)
      iwolfe                     =  ictrl(14)
      n_lbfgs                    =  ictrl(15)
      cplex_method               =  ictrl(16)
      cplex_feas                 =  ictrl(17)
      alg_unconstr               =  ictrl(18)
!     idum1                      =  ictrl(19)
!     nadd                       =  ictrl(20)
!     idum2                      =  ictrl(21) 
      sph_dia_hess               =  ictrl(22)
      current_mode_num           =  ictrl(23)
      QPform_A                   =  ictrl(24)
      QPform_H                   =  ictrl(25) 
!
      convex                     =  lctrl(1)
      random_start               =  lctrl(2)
      finite_diff                =  lctrl(3)
      check_grad                 =  lctrl(4)
      conservative               =  lctrl(5)
      trust_region               =  lctrl(6)
      unconstrained              =  lctrl(7)
      use_active_set             =  lctrl(8)
      allow_f                    =  lctrl(9)
      allow_c                    =  lctrl(10)
      classic_conserv            =  lctrl(11)
      classic_trust              =  lctrl(12)
      slack                      =  lctrl(13)
      relax                      =  lctrl(14)
      multimodal                 =  lctrl(15)
      debug                      =  lctrl(16)
      special                    =  lctrl(17)
      calc_kkt                   =  lctrl(18)
      strict_struct              =  lctrl(19)
      mayStop                    =  lctrl(20)
!     ldum1                      =  lctrl(21)
      mayReturn                  =  lctrl(22)
      GetFinDiffGradHessDense    =  lctrl(23)
      FormHessSparse             =  lctrl(24)
      ResetDualVars              =  lctrl(25)

!
      big                        =  rctrl(1)
      xtol                       =  rctrl(2)
      ftol                       =  rctrl(3)
      feaslim                    =  rctrl(4)
      deltx                      =  rctrl(5)
      ci_slack                   =  rctrl(6)
      di_slack                   =  rctrl(7)
      actlim                     =  rctrl(8)
      tlimit                     =  rctrl(9)
      xtol_inf                   =  rctrl(10)
      rho_l_min                  =  rctrl(11)
      rho_l_max                  =  rctrl(12)
      biglam                     =  rctrl(13)
      kkt_tol                    =  rctrl(14)
      sfracmin                   =  rctrl(15)
      atol1                      =  rctrl(16)
      btol1                      =  rctrl(17)
      fapriori                   =  rctrl(18)
      atol2                      =  rctrl(19)
      btol2                      =  rctrl(20)
      filter_hi                  =  rctrl(21)
      pen1                       =  rctrl(22)
      pen2                       =  rctrl(23)
      dummy1                     =  rctrl(24)
      ymax                       =  rctrl(25)
      mma_lo                     =  rctrl(26)
      mma_hi                     =  rctrl(27)
      ptarget                    =  rctrl(28)
      tol_bayes                  =  rctrl(29)
      dml_infinity               =  rctrl(30)
      rangemax                   =  rctrl(31)
      ztol                       =  rctrl(32)
      kkt_min                    =  rctrl(33)
      drand0                     =  rctrl(34)
      priv1                      =  rctrl(35)
      DualTrustRadius            =  rctrl(36)
!     dmatengine                 =  rctrl(37) 
!     rdum1                      =  rctrl(38)
      shift_tol                  =  rctrl(39) 
      deltxH                     =  rctrl(40)
!
      cname1                     =  cctrl(1)
      cname2                     =  cctrl(2)
!     func_f_name                =  cctrl(3)
!     grad_f_name                =  cctrl(4)
!     hess_f_name                =  cctrl(5)
!
! check if ictrl,lctrl,rctrl, cctrl are passed to the subroutine where the headers are invoked
!
      if (ictrl(nictrl).ne.123456789) then
        stop ': header files are called or passed incorrectly (i)'
      endif
      if (rctrl(nrctrl).ne.123456789.d0) then
        stop ': header files are called or passed incorrectly (r)'
      endif
      if (cctrl(ncctrl).ne.'0a1b2c3d4e5f6g7h8i9j10kX') then
        stop ': header files are called or passed incorrectly (c)'
      endif
      if (.not.lctrl(nlctrl)) then
        stop ': header files are called or passed incorrectly (l_a)'
      endif
      if (lctrl(nlctrl-1)) then
        stop ': header files are called or passed incorrectly (l_b)'
      endif
      if (.not.lctrl(nlctrl-2)) then
        stop ': header files are called or passed incorrectly (l_c)'
      endif
!
! check if the user routines are passed to the subroutine where these routines are used
!
      if (iuser(imax).ne.123456789) then
        stop ': user files are called or passed incorrectly (i)'
      endif
      if (ruser(rmax).ne.123456789.d0) then
        stop ': user files are called or passed incorrectly (r)'
      endif
      if (cuser(cmax).ne.'0a1b2c3d4e5f6g7h8i9j10kX') then
        stop ': user files are called or passed incorrectly (c)'
      endif
      if (.not.luser(lmax)) then
        stop ': user files are called or passed incorrectly (l_a)'
      endif
      if (luser(lmax-1)) then
        stop ': user files are called or passed incorrectly (l_b)'
      endif
      if (.not.luser(lmax-2)) then
        stop ': user files are called or passed incorrectly (l_c)'
      endif
!
