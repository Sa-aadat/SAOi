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
      integer          i, j, n, m, ni, ne, nnz, Acol(*), Aptr(*), k, l
      integer          ndisk, isum
      double precision x(n), gf(n), gc(*)
!
      ndisk = (n-2)/2
!
      do i=1,n-2
        gf(i) = 0.d0
      enddo
      gf(n-1) = x(n)
      gf(n)   = x(n-1)
!
      k=0
      do i=1,n
        do j=1,ni
          k = (i-1)*ni + j
          gc(k) = 0.d0
        enddo
      enddo
!
!
      do i=1,ndisk
        k = (i-1)*ni + i
        gc(k) = 1.d0
!
        k = (i+ndisk-1)*ni + i
        gc(k+ndisk) = 1.d0
      enddo
!
      do i=1,ndisk
        k = (n-1-1)*ni + i
        gc(k) = -1.d0
!
        k = (n-1)*ni + i
        gc(k+ndisk) = -1.d0
      enddo
!
      l=0
      do i=1,ndisk-1
        do j=i+1,ndisk
          l = l+1
!
          k = (i-1)*ni + l+ndisk
          gc(k+ndisk) = -2.d0*(x(i)-x(j))
!
          k = (j-1)*ni + l+ndisk
          gc(k+ndisk) =  2.d0*(x(i)-x(j))
!
          k = (i+ndisk-1)*ni + l+ndisk
          gc(k+ndisk) = -2.d0*(x(i+ndisk)-x(j+ndisk))
!
          k = (j+ndisk-1)*ni + l+ndisk
          gc(k+ndisk) =  2.d0*(x(i+ndisk)-x(j+ndisk))
         enddo
      enddo
!
      return
      end subroutine SAOi_grads
