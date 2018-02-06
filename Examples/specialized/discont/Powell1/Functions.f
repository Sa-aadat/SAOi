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
      double precision x1, x2
!
      x1 = x(1)
      x2 = x(2)
!      
      f = x1**2 + x2**2   
!
      if (dabs(x2).le.0.5d0) then
        c(1) = -3.d0/4.d0*x1**2 - 1.d0/4.d0*x1*x2 -  4.d0/ 3.d0*x2**2
     &         -1.d0/2.d0*x1    + 1.d0/3.d0*x2   + 29.d0/12.d0         
!        
      elseif (x2.lt.-0.5d0) then
        c(1) = -3.d0/4.d0*x1**2 - 1.d0/4.d0*x1*x2 - 1.d0/3.d0*x2**2
     &         -1.d0/2.d0*x1    + 4.d0/3.d0*x2   + 8.d0/3.d0           
!     
      elseif (x2.gt. 0.5d0) then
        c(1) = -3.d0/4.d0*x1**2 - 1.d0/4.d0*x1*x2 - x2**2
     &         -1.d0/2.d0*x1    + 5.d0/2.d0           
!
      endif
!
      return
      end subroutine SAOi_funcs
