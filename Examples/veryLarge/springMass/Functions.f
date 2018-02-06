! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: objective and constraint functions
! SAOi:

      subroutine SAOi_funcs (n, m, ni, ne, x, f, c, iuser, luser, cuser,
     &                       ruser, eqn, lin, ictrl, lctrl, rctrl,
     &                       cctrl)
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the objective function f and the inequality constraint      !
!  functions c(j), j=1,ni                                              !
!                                                                      !
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          i, j, n, m, ni, ne, nmass
      double precision f, x(n), c(ni)
      double precision al, ali, pe1, pe2, xi, yi, wi, xip1, yip1
      double precision aki, dli, dli1
!
      if (mod(n,2).ne.0) stop ' n is odd '
!
      nmass=n/2
!
      al=60.
      ali=al/(dble(nmass)+1.)
      pe1=0.
      pe2=0.
      xi=0.
      yi=0.
      do 20 i=1,nmass
         j=i-1
         wi=50.*dble(i)
         pe1=pe1+wi*x(i+nmass)
         xip1=x(i)
         yip1=x(i+nmass)
         aki=500.+200.*((dble(nmass)/3.-dble(i))**2)
         dli1=dsqrt((ali+xip1-xi)**2+(yip1-yi)**2)
         dli=dli1-ali
         pe2=pe2+0.5*aki*(dli**2)
         xi=xip1
         yi=yip1
20    continue

!  last spring
      aki=500.+200.*((dble(nmass)/3.-dble(nmass+1))**2)
      dli1=dsqrt((ali-xi)**2+yi**2)
      dli=dli1-ali
      pe2=pe2+0.5*aki*(dli**2)
      f=pe1+pe2
!
      return
      end subroutine SAOi_funcs
