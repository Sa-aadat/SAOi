! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi:
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!               This is a public SAOi file, edit at will               !
!                                                                      !
!                                                                      !
!  The user arrays below are available in Initialize.f, Gradients.f    !
!       Functions.f and Hessian.f, and may be used to conveniently     !
!       pass information around. Do not be tempted to add COMMON       !
!       blocks, which may limit very large scale functionality;        !
!       the user may specify arbitrary non-zero lengths for these      !
!       arrays.                                                        !
!                                                                      !
!  nmax is the maximum number of primal design variables, which has    !
!       to be set adequately large by the user.                        !
!                                                                      !
!  mmax is the maximum number of dual variables, which has to be set   !
!       adequately large by the user.                                  !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      integer           nmax             ! maximum number of primal vars
      parameter        (nmax=500000)
!
      integer           mmax             ! maximum number of dual vars
      parameter        (mmax=500000)
!
      integer           iuser(50)        ! integer user array
!
      logical           luser(50)        ! logical user array
!
      character*15      cuser(50)        ! character user array
!
      double precision  ruser(50)        ! real user array
!

