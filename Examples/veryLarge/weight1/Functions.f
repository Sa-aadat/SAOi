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
      logical          eqn(*), lin(*)
      include          'ctrl.h'
      integer          i, j, n, m, ni, ne, ns
      double precision f, x(n), c(ni)
!
      ruser(1) = 1.d0
!
      f = 0.d0
      do i=1,n
        f = f + x(i)
      enddo
!
      ns = iuser(1)
!
      c(1)= 0.d0
      do i=1,ns
        c(1) = c(1) + 1.d0/x(i)
      enddo
!
      c(2) = c(1)
      do i=ns+1,n
        c(1) = c(1) + 1.d-6/x(i)
        c(2) = c(2) - 1.d-6/x(i)
      enddo
!
      c(1) = c(1) - 1000.d0
      c(2) = c(2) - 900.d0
!
      if (ruser(1).ne.1.d0) c = c*ruser(1) ! scale
!
      return
      end subroutine SAOi_funcs
