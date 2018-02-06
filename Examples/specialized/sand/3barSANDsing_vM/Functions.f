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
      double precision f, x(n), c(m) 
      double precision l(3), ss(3), a(3), p(3), u, v, k11, k12, k22
      double precision E, Sy, Fl, d, h 
!
      E  = 3.d6
      Sy = 3.d4
      Fl = 1.d4
      d  = 0.29d0
      h  = 1.d2
      p(1) = -100.d0                   ! the positions
      p(2) =    0.d0                   ! the positions
      p(3) =  100.d0                   ! the positions
!
      do i=1,3
        a(i) = x(i)                    ! the areas
        l(i) = dsqrt(p(i)**2 + h**2)   ! the lengths
      enddo
!
      u = x(4)
      v = x(5)
! 
      k11 = 0.d0
      k12 = 0.d0
      k22 = 0.d0
      do i=1,3
        k11 = k11 + a(i)*p(i)**2/l(i)**3
        k12 = k12 - h*a(i)*p(i)/l(i)**3
        k22 = k22 + Fl*a(i)/l(i)**3
      enddo
      k11 = k11*E
      k12 = k12*E
      k22 = k22*E
!
      do i=1,3
        ss(i) = E*(-u*p(i) + h*v)/l(i)**2
      enddo
!
      f = d*(a(1)*l(1) + a(2)*l(2) + a(3)*l(3))
!
      do i=1,3
        c(i)   =  dsqrt(ss(i)**2)/Sy - 1.d0     ! Von Mises stress constraints without vanishing
      enddo
      c(4) = k11*u + k12*v - Fl  ! equilibrium in x_1 
      c(5) = k12*u + k22*v       ! equilibrium in x_2 
!
      return
      end subroutine SAOi_funcs
