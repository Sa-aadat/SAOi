! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: Include file to overcome the need for common blocks present
! SAOi: in earlier versions.
! SAOi:
!


      ictrl(1)  =  approx_f
      ictrl(2)  =  approx_c
      ictrl(3)  =  force_converge
      ictrl(4)  =  subsolver
      ictrl(5)  =  outermax
      ictrl(6)  =  innermax
      ictrl(7)  =  ksubmax
      ictrl(8)  =  itglobalmax
      ictrl(9)  =  iprint
      ictrl(10) =  iaggressive
      ictrl(11) =  structure
      ictrl(12) =  nprob
      ictrl(13) =  ifalk
      ictrl(14) =  iwolfe
      ictrl(15) =  n_lbfgs
      ictrl(16) =  cplex_method
      ictrl(17) =  cplex_feas
      ictrl(18) =  alg_unconstr
!     ictrl(19) =  idum1
!     ictrl(20) =  nadd
!     ictrl(21) =  idum2
      ictrl(22) =  sph_dia_hess
      ictrl(23) =  current_mode_num
      ictrl(24) =  QPform_A 
      ictrl(25) =  QPform_H 
!
      lctrl(1)  =  convex
      lctrl(2)  =  random_start
      lctrl(3)  =  finite_diff
      lctrl(4)  =  check_grad
      lctrl(5)  =  conservative
      lctrl(6)  =  trust_region
      lctrl(7)  =  unconstrained
      lctrl(8)  =  use_active_set
      lctrl(9)  =  allow_f
      lctrl(10) =  allow_c
      lctrl(11) =  classic_conserv
      lctrl(12) =  classic_trust
      lctrl(13) =  slack
      lctrl(14) =  relax
      lctrl(15) =  multimodal
      lctrl(16) =  debug
      lctrl(17) =  special
      lctrl(18) =  calc_kkt
      lctrl(19) =  strict_struct
      lctrl(20) =  mayStop
!     lctrl(21) =  ldum1
      lctrl(22) =  mayReturn
      lctrl(23) =  GetFinDiffGradHessDense 
      lctrl(24) =  FormHessSparse
      lctrl(25) =  ResetDualVars
!
      rctrl(1)  =  big
      rctrl(2)  =  xtol
      rctrl(3)  =  ftol
      rctrl(4)  =  feaslim
      rctrl(5)  =  deltx
      rctrl(6)  =  ci_slack
      rctrl(7)  =  di_slack
      rctrl(8)  =  actlim
      rctrl(9)  =  tlimit
      rctrl(10) =  xtol_inf
      rctrl(11) =  rho_l_min
      rctrl(12) =  rho_l_max
      rctrl(13) =  biglam
      rctrl(14) =  kkt_tol
      rctrl(15) =  sfracmin
      rctrl(16) =  atol1
      rctrl(17) =  btol1
      rctrl(18) =  fapriori
      rctrl(19) =  atol2
      rctrl(20) =  btol2
      rctrl(21) =  filter_hi
      rctrl(22) =  pen1
      rctrl(23) =  pen2
      rctrl(24) =  dummy1
      rctrl(25) =  ymax
      rctrl(26) =  mma_lo
      rctrl(27) =  mma_hi
      rctrl(28) =  ptarget
      rctrl(29) =  tol_bayes
      rctrl(30) =  dml_infinity
      rctrl(31) =  rangemax
      rctrl(32) =  ztol
      rctrl(33) =  kkt_min
      rctrl(34) =  drand0
      rctrl(35) =  priv1
      rctrl(36) =  DualTrustRadius
!     rctrl(37) =  dmatengine
!     rctrl(38) =  rdum1           
      rctrl(39) =  shift_tol
      rctrl(40) =  deltxH
!
      cctrl(1)  =  cname1
      cctrl(2)  =  cname2
!     cctrl(3)  =  func_f_name
!     cctrl(4)  =  grad_f_name
!     cctrl(5)  =  hess_f_name
!
