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
      double precision temp0, temp1, temp2, temp3
!
!
! the objective function
      temp0 =0.d0
      do i=1,n
       temp1 =  0.1D0**(DFLOAT(N-I)/DFLOAT(N-1))
       temp2 =  2.0D0
       temp0 = temp0 + (temp2 + 0.5d0*temp1*x(i))*x(i)
      enddo
      f = temp0
!
! the constraints
      do j=1,ni
        temp3 = -0.1d0
        do i=1,n
          temp3 = temp3 + x(i)*(DABS(DFLOAT(I-J))/DFLOAT(I+J)-0.5D0)
        enddo
      c(j) = temp3
      enddo
      
      return
      end subroutine SAOi_funcs
