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
      double precision f1, f2, t1, t2
!
      f1=0.d0
      f2=0.d0
!
      do i=1,5
        t1=i*dcos((i+1.d0)*x(1)+i)
        t2=i*dcos((i+1.d0)*x(2)+i)
        f1=f1+t1
        f2=f2+t2
      enddo
      f=f1*f2
!
      return
      end subroutine SAOi_funcs
