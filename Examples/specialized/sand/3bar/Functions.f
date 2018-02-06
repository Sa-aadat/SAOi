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
      double precision l(3), s(3), a(3), p(3), u, v, k11, k12, k22
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
      if (k12.ne.0.d0) then
        v = Fl/(k12 - k11*k22/k12)        ! not really :-) Need to solve for q = K^-1 P ... !
        u = -k22/k12*v
      else
        v = 0.d0
        u = Fl/k11
      endif
! 
      do i=1,3
        s(i) = E*(-u*p(i) + h*v)/l(i)**2  ! not really :-) Need to do post-processing, which gives div by 0 !
      enddo
!
      f = d*(a(1)*l(1) + a(2)*l(2) + a(3)*l(3))
!
      do i=1,3
        c(i)   =  s(i)/Sy - 1.d0     ! tension
        c(i+3) = -s(i)/Sy - 1.d0     ! compression
      enddo
!
      do i=1,3
        write(*,*) ' Stress: ',i,   Sy*(c(i) + 1.d0)   ! =  s(i)       ! tension
        write(*,*) ' Stress: ',i+3, Sy*(c(i+3) + 1.d0) ! =  -s(i)      ! compression      
      enddo
!      
      return
      end subroutine SAOi_funcs
