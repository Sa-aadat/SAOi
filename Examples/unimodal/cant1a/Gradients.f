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
!      please see the users manual                                     !
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
      double precision c1, c2
!
      c1=0.0624
      do i=1,n
        gf(i)=c1
      enddo
      gc(1)=-3.d0 * 61.d0 / x(1)**4
      gc(2)=-3.d0 * 37.d0 / x(2)**4
      gc(3)=-3.d0 * 19.d0 / x(3)**4
      gc(4)=-3.d0 *  7.d0 / x(4)**4
      gc(5)=-3.d0 *  1.d0 / x(5)**4
!
      return
      end subroutine SAOi_grads
