! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: initialization for automatic testing
! SAOi:


      subroutine SAOi_force (n,m,ni,ne,ictrl,lctrl,rctrl,cctrl,
     &                       iuser,luser,cuser,ruser)
!----------------------------------------------------------------------!
!                                                                      !
!  do not edit this routine, it is only of interest during automatic   !
!  testing of the SAOi algorithm; the values specified here may        !
!  be changed in subroutine initialize.f, if necessary                 !
!                                                                      !
!  NOTE: the settings enforced here will not be checked for            !
!        consistency or validity...                                    !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      integer          n,m,ni,ne
      include          'ctrl.h'
!
      include 'ctrl_get.inc'
!
!     mayStop                                      =  .true.
!     GetFinDiffGradHessDense                      =  .true.
!     finite_diff                                  =  .true.
!     check_grad                                   =  .true.
!     deltx                                        =  1.0d-6
!     deltxH                                       =  1.0d-4
!     ifalk                                        =  1
!     force_converge                               =  2
!     approx_f                                     =  1
!     approx_c                                     =  1
!     if (subsolver.ne.27) subsolver               =  25 
!     subsolver                                    =  27
!     xtol                                         =  1.d-4
!     xtol_inf                                     =  1.d-6
!     kkt_min                                      =  1.d-3
!     kkt_tol                                      =  1.d-10
!     ftol                                         =  1.d-10
!     outermax                                     =  500
!     innermax                                     =  50
!     dml_infinity                                 =  1.d0
!     random_start                                 =  .true.
!     biglam                                       =  1.d8
!     atol1                                        =  1.d-5
!     btol1                                        = -1.3d0
!     pen1                                         =  1.d5
! !
!       if (ni+ne.eq.0) subsolver = 1
! !
      include 'ctrl_set.inc'
!
      return
      end subroutine SAOi_force
