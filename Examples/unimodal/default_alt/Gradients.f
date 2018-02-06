! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: objective and constraint functions gradients
! SAOi:

      subroutine SAOi_grads (n, m, ni, ne, x, gf, gc, nnz, Acol, Aptr,
     &                       iuser, luser, cuser, ruser, eqn, lin,
     &                       ictrl, lctrl, rctrl, cctrl)
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the gradients gf(i) of the objective function f, and the    !
!  derivatives gc(k) of the inequality constraint functions c          !
!  w.r.t. the variables x(i)                                           !
!                                                                      !
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!                                                                      !
!  The default storage scheme is the dense storage scheme. In          !
!      this scheme, it is only required to specify the vector          !
!      gf and the matrix gc. The dimensions of gf and gc are           !
!      gf(n) and gc(ni,n) respectively. For details, please            !
!      see the users manual                                            !
!                                                                      !
!  The default storage scheme if the sparse implementation is          !
!      selected is the compressed sparse row (CSR) storage scheme.     !
!      The dimensions of Acol and Aptr in the CSR storage sheme        !
!      are Acol(nnz) and Aptr(ni+1) respectively. For details,         !
!      please see the users manual. In short, the derivative           !
!      of constraint j w.r.t. variable i is stored in location         !
!      k = (i-1)*ni + j, where ni is the total number of               !
!      constraints                                                     !
!                                                                      !
!  nnz is the number of non-zero entries in gc. If the standard        !
!      sparse storage structure is used, nnz should be declared        !
!      in Initialize.f, and its value should not be changed            !
!                                                                      !
!  Acol and Aptr are to be declared here if the sparse storage         !
!      structure is used. For details, please see the users manual     !
!                                                                      !
!  iuser, luser, cuser and ruser are user arrays, which may be used    !
!      at will to pass arbitrary data around between the user          !
!      routines                                                        !
!                                                                      !
!----------------------------------------------------------------------!
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          i, j, n, m, ni, ne, nnz, Acol(*), Aptr(*)
      double precision x(n), gf(n), gc(*)
      double precision temp0, temp1, temp2, temp3, temp4
!
! some often occurring terms
      temp0 = 0.124d0
      temp1 = dsqrt(1.d0 + x(2)**2)
      temp2 = 8.d0/x(1) + 1.d0/(x(1)*x(2))
      temp3 = 8.d0/x(1) - 1.d0/(x(1)*x(2))
      temp4 = 2.d0*x(2)

! derivatives of the objective function
!----------------------------------------------------------------------!
      gf(1) = temp1
      gf(2) = x(1)/(2.d0*temp1)*temp4

! derivative of constraint 1 w.r.t variable 1: k = (1-1)*2 + 1 = 1
!----------------------------------------------------------------------!
      gc(1) = -temp0*temp1*(8.d0/x(1)**2 + 1.d0/(x(1)**2*x(2)))

! derivative of constraint 2 w.r.t variable 1: k = (1-1)*2 + 2 = 2
!----------------------------------------------------------------------!
      gc(2) = -temp0*temp1*(8.d0/x(1)**2 - 1.d0/(x(1)**2*x(2)))

! derivative of constraint 1 w.r.t variable 2: k = (2-1)*2 + 1 = 3
!----------------------------------------------------------------------!
      gc(3) = temp0/(2.d0*temp1)*temp4*temp2
     &      - temp0*temp1/(x(1)*x(2)**2)

! derivative of constraint 2 w.r.t variable 2: k = (2-1)*2 + 2 = 4
!----------------------------------------------------------------------!
      gc(4) = temp0/(2.d0*temp1)*temp4*temp3
     &      + temp0*temp1/(x(1)*x(2)**2)
!
      return
      end subroutine SAOi_grads
