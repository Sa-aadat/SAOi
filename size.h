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
      parameter        (nmax=10000000)
!
      integer           mmax             ! maximum number of dual vars (i.e. the number of constraints,
      parameter        (mmax=10000000)     !    but not including bound constraints on primal or dual vars)
!
      integer           imax             ! size of integer user array iuser
      parameter        (imax=50)
!
      integer           lmax             ! size of logical user array luser 
      parameter        (lmax=50)
!
      integer           cmax             ! size of character user array cuser 
      parameter        (cmax=50)
!
      integer           rmax             ! size of real (dble precision) user array ruser 
      parameter        (rmax=50)
!
      integer           iuser(imax)      ! integer user array
!
      logical           luser(lmax)      ! logical user array
!
      character*24      cuser(cmax)      ! character user array
!
      double precision  ruser(rmax)      ! real user array
!

