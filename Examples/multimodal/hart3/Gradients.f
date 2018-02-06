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
      double precision wp3, he3(4,3), hp3(4,3), hcd3(4), fr3(4)
!
      data he3/3.d0,1.0d-1,3.0d0,1.0d-1,4*1.0d1,3.0d1,3.5d1,3.0d1,3.5d1/
      data hcd3/1.0d0,1.2d0,3.0d0,3.2d0/
      data hp3/3.689d-1,4.699d-1,1.091d-1,3.815d-2,1.17d-1,4.387d-1,
     &         8.732d-1,5.743d-1,2.673d-1,7.47d-1,5.547d-1,8.828d-1/
!
      do i=1,4
        wp3=0.d0
        fr3(i)=0.d0
        do j=1,n
          wp3=wp3-he3(i,j)*(x(j)-hp3(i,j))**2
        enddo
        fr3(i)=fr3(i)-hcd3(i)*dexp(wp3)
      enddo
!
      do j=1,n
        gf(j) = 0.d0
        do i=1,4
          gf(j) = gf(j) +  fr3(i) * (-2.d0*he3(i,j)*(x(j)-hp3(i,j)))
        enddo
      enddo
!
      return
      end subroutine SAOi_grads
