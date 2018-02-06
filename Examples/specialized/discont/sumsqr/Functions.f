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
      double precision factr, temp, temp1, factr1, di
!
      factr = 1.5d0
      temp = 0.d0
      do i = 1,n
        temp = temp + x(i)
      enddo
      temp1 = dsin(temp/10.d0)
!
      f = 0.d0
!
      do i = 1,n
        di = dble(i)
        if (temp1.gt.0.5d0) then
          f = f + di*x(i)**2/factr
!
        elseif (temp1.lt.-0.5d0) then
          f = f + di*x(i)**2*factr
!
        elseif (temp1.gt.0.d0) then
          f = f + di*x(i)**2
!
        else
          f = f + di*x(i)**2 + 1.d0/dble(n)
!
        endif
      enddo
!
      return
      end subroutine SAOi_funcs
