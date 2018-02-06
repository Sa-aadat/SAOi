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
      double precision pi, d_theta, R_v, R_min, R_max, alpha_c
!
      pi               =  4.0d0*datan(1.0d0)
      d_theta          =  2.d0*pi/(5.d0*dble(n+1))
      R_v              =  1.0d0
      R_min            =  1.0d0
      R_max            =  2.0d0
      alpha_c          =  1.5d0
!
      f = 0.d0
      do i=1,n
        f = f + x(i)
      enddo
      f = - f*(pi*R_v)/dble(n)
!
      do j = 2,n-1
        c(j-1) = - x(j-1)*x(j) - x(j)*x(j+1)
     &           + 2.d0*x(j-1)*x(j+1)*dcos(d_theta)
      enddo
!
      c(n-1) = - R_min*x(1) - x(1)*x(2) + 2.d0*R_min*x(2)*dcos(d_theta)
      c(n)   = - R_min**2 - R_min*x(1) + 2.d0*R_min*x(1)*dcos(d_theta)
      c(n+1) = - x(n-1)*x(n) - x(n)*R_max
     &         + 2.d0*x(n-1)*R_max*dcos(d_theta)
      c(n+2) = - 2.d0*R_max*x(n) + 2*x(n)**2*dcos(d_theta)
!
      do j=1,n-1
        c(j+n+2) = -alpha_c*d_theta - (x(j+1) - x(j))
      enddo
!
      c(2*n+2) = -alpha_c*d_theta - (x(1) - R_min)
      c(2*n+3) = -alpha_c*d_theta - (R_max - x(n))
!
      do j=1,n-1
        c(j+2*n+3) =  (x(j+1) - x(j)) - alpha_c*d_theta
      enddo
!
      c(3*n+3) = (x(1) - R_min) - alpha_c*d_theta
      c(3*n+4) = (R_max - x(n)) - alpha_c*d_theta
!
      return
      end subroutine SAOi_funcs
