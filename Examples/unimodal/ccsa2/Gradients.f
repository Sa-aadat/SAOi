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
      integer          i, j, k, n, m, ni, ne, nnz, Acol(*), Aptr(*)
      double precision x(n), gf(n), gc(*)
      double precision sij, pij, qij, temp1, alphaij, pi
!
      pi = 4.0d0*datan(1.0d0)
!
      do i = 1,n ! note that Svanberg required that S, P, Q are not stored, but re-calculated every time !
        gf(i) = 0.d0
        k=(i-1)*ni+1
        gc(k) = 0.d0
        k=(i-1)*ni+2
        gc(k) = 0.d0
        do j = 1,n
          alphaij = dble(i+j-2)/dble(2*n-2)
          temp1 = (1.d0+dabs(dble(i-j)))*dlog(dble(n))
          sij = (2.d0+dsin(4.d0*pi*alphaij))/temp1
          pij = (1.d0+2.d0*alphaij)/temp1
          qij = (3.d0-2.d0*alphaij)/temp1
          gf(i) = gf(i)-sij*x(j)*2.d0
          k = (i-1)*ni+1
          gc(k) = gc(k)+pij*x(j)*2.d0
          k = (i-1)*ni+2
          gc(k) = gc(k)+qij*x(j)*2.d0
        enddo
      enddo
!
      return
      end subroutine SAOi_grads
