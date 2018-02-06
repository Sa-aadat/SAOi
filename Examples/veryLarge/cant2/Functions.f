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
      double precision Length, P, E, li, ymax1, t, Ro1, inerti
      double precision sumlj, temp1, temp2, ydii, y1ii, ydi, yi
!
      Length = 500.d0
      P      = 40.d0
      E      = 2.0d8
      li     = Length/dble(n)
      ymax1  = 0.5d0
      t      = 0.2d0
      Ro1    = 0.00078d0
!
      f=0.d0
      do i=1,n
        f=f+4.d0*x(i)*t*li*Ro1
      enddo
!
      ydi = 0.d0
      yi  = 0.d0
      do i=1,n
        inerti= 2.d0*(x(i)**3)*t/3.d0
        sumlj = li*i
        temp1 = (Length+li/2.d0-sumlj)
        temp2 = (Length+2.d0*li/3.d0-sumlj)
!
        ydii = P*li   /(     E*inerti)*temp1      + ydi
        y1ii = P*li*li/(2.d0*E*inerti)*temp2 + yi + ydi*li
        ydi = ydii
        yi  = y1ii
!
      enddo
      c(1)=y1ii/ymax1 - 1.d0 ! constraint on tip displacement
!
      return
      end subroutine SAOi_funcs
