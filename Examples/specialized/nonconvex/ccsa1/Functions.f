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
      double precision temp1, alphaij, sij, pij, qij, pi
!
      pi = 4.0d0*datan(1.0d0)
!
      if (n.le.1) stop ' n <= 1 '
!
      f = 0.d0
      c(1) = dble(n)/2.d0
      c(2) = dble(n)/2.d0
!
      do i = 1,n ! note that Svanberg required that S, P, Q are not stored, but re-calculated every time !
        do j = 1,n
          alphaij = dble(i+j-2)/dble(2*n-2)  ! hence alpha \in [0 1] \forall i,j
          temp1 = (1.d0+dabs(dble(i-j)))*dlog(dble(n))
          sij = (2.d0+dsin(4.d0*pi*alphaij))/temp1
          pij = (1.d0+2.d0*alphaij)/temp1
          qij = (3.d0-2.d0*alphaij)/temp1
          f = f+sij*x(i)*x(j)
          c(1) = c(1)-pij*x(i)*x(j)
          c(2) = c(2)-qij*x(i)*x(j)
        enddo
      enddo
!
      return
      end subroutine SAOi_funcs
