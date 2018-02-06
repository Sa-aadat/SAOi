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
      integer          i, j, n, m, ni, ne, ndisk, l
      double precision f, x(n), c(ni) 
!
      ndisk = (n-2)/2
!
      f=x(n-1)*x(n)
!
      do i=1,ndisk
        c(i)       = x(i)       + 1.d0 - x(n-1)
        c(i+ndisk) = x(i+ndisk) + 1.d0 - x(n)
      enddo
!
      l=0
      do i=1,ndisk-1
        do j=i+1,ndisk
          l=l+1
          c(2*ndisk+l) = -(x(i)-x(j))**2 - (x(i+ndisk)-x(j+ndisk))**2
     &                 + 4.d0
        enddo
      enddo
!
      return
      end subroutine SAOi_funcs
