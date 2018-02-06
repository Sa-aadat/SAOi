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
      double precision wp3, he3(4,3), hp3(4,3), hcd3(4), fr3(4)
!
      data he3/3.d0,1.0d-1,3.0d0,1.0d-1,4*1.0d1,3.0d1,3.5d1,3.0d1,3.5d1/
      data hcd3/1.0d0,1.2d0,3.0d0,3.2d0/
      data hp3/3.689d-1,4.699d-1,1.091d-1,3.815d-2,1.17d-1,4.387d-1,
     &         8.732d-1,5.743d-1,2.673d-1,7.47d-1,5.547d-1,8.828d-1/
!
      f=0.0d0
      do i=1,4
        wp3=0.d0
        fr3(i)=0.d0
        do j=1,n
          wp3=wp3-he3(i,j)*(x(j)-hp3(i,j))**2
        enddo
        f=f-hcd3(i)*dexp(wp3)
      enddo
!
      return
      end subroutine SAOi_funcs
