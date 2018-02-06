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
      double precision wp6, he6(4,6), hp6(4,6), hcd6(4), fr6(4)
!
      data he6 /1.0d1,5.0d-2,3.0d0,1.7d1,3.0d0,1.0d1,3.5d0,8.0d0,1.7d1,
     &          1.7d1,1.7d0,5.0d-2,3.5d0,1.0d-1,1.0d1,
     &          1.0d1,1.7d0,8.0d0,1.7d1,1.0d-1,8.0d0,1.4d1,8.0d0,1.4d1/
      data hcd6/1.0d0,1.2d0,3.0d0,3.2d0/
      data hp6 /1.312d-1,2.329d-1,2.348d-1,4.047d-1,1.696d-1,4.135d-1,
     &          1.451d-1,8.828d-1,5.569d-1,8.307d-1,3.522d-1,8.732d-1,
     &          1.24d-2,3.736d-1,2.883d-1,5.743d-1,8.283d-1,1.004d-1,
     &          3.047d-1,1.091d-1,5.886d-1,9.991d-1,6.65d-1,3.81d-2/
!
      f=0.0d0
      do i=1,4
        wp6=0.d0
        fr6(i)=0.d0
        do j=1,n
          wp6=wp6-he6(i,j)*(x(j)-hp6(i,j))**2
        enddo
        f=f-hcd6(i)*dexp(wp6)
      enddo
!
      return
      end subroutine SAOi_funcs
