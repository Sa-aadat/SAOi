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
      integer          i, j, n, m, ni, ne
      double precision f, x(n), c(ni) 
      double precision xij, xij2, prod
!
      xij=4000.d0
      xij2=xij/2.d0
!
      prod=1.0d0
      do i=1,n
        prod=prod*dcos((x(i)-0.d0)/i**0.5d0)
      enddo
      f=0.d0
      do i=1,n
        f=f+(x(i)-0.d0)**2
      enddo
      f=f/xij+1.d0-prod
!
      return
      end subroutine SAOi_funcs
