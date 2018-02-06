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
      integer          i, j, n, m, ni, ne, nmat, nsta, nnz 
      double precision f, x(n), c(m), du, p
!
      du   = ruser(1)
      p    = ruser(2) 
      nnz  = iuser(1)
      nmat = iuser(2)
      nsta = iuser(3)
!      
      call globalD2 (2, x, n, nmat, nsta, m, nnz, du, f, c, eqn, p) 
!
      return
      end subroutine SAOi_funcs
