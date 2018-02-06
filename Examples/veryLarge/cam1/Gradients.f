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
      integer          i, j, k, n, m, ni, ne, nseg, ij
      integer          nnz, Acol(*), Aptr(*), itmp(ni)
      double precision x(n), gf(n), gc(*)
      double precision pi, d_theta, R_v, R_min, R_max, alpha_c, temp
!
      pi               =  4.0d0*datan(1.0d0)
      d_theta          =  2.d0*pi/(5.d0*dble(n+1))
      R_v              =  1.0d0
      R_min            =  1.0d0
      R_max            =  2.0d0
      alpha_c          =  1.5d0
!
      temp = -(pi*R_v)/dble(n)
      do i=1,n
        gf(i) = temp
      enddo
!
      do j=1,ni
        itmp(j) = 0
      enddo
!
      ij=0
      do j=2,n-1
        ij=ij+1
        Acol(ij)=j-1
        gc(ij) = 2.d0*x(j+1)*dcos(d_theta)-x(j)
        ij=ij+1
        Acol(ij)=j
        gc(ij)   =-x(j-1)-x(j+1)
        ij=ij+1
        Acol(ij)=j+1
        gc(ij) = 2.d0*x(j-1)*dcos(d_theta)-x(j)
        itmp(j-1) = itmp(j-1) + 3
      enddo
!
      ij=ij+1
      Acol(ij)=1
      gc(ij)   =-R_min-x(2)
      ij=ij+1
      Acol(ij)=2
      gc(ij)   = 2.d0*R_min*dcos(d_theta)-x(1)
      ij=ij+1
      Acol(ij)=1
      gc(ij)     = 2.d0*R_min*dcos(d_theta)-x(1)
      ij=ij+1
      Acol(ij)=n-1
      gc(ij) = 2.d0*R_max*dcos(d_theta)-x(n)
      ij=ij+1
      Acol(ij)=n
      gc(ij)   =-R_max-x(n-1)
      ij=ij+1
      Acol(ij)=n
      gc(ij)   = 4.d0*x(n)*dcos(d_theta)-2.d0*R_max
      itmp(n-1)   = itmp(n-1) + 2
      itmp(n  )   = itmp(n  ) + 1
      itmp(n+1)   = itmp(n+1) + 2
      itmp(n+2)   = itmp(n+2) + 1
!
      do j=1,n-1
        ij=ij+1
        Acol(ij)=j
        gc(ij)   = 1.d0
        ij=ij+1
        Acol(ij)=j+1
        gc(ij) =-1.d0
        itmp(j+n+2)   = itmp(j+n+2) + 2
      enddo
!
      ij=ij+1
      Acol(ij)=1
      gc(ij) =-1.d0
      ij=ij+1
      Acol(ij)=n
      gc(ij) = 1.d0
      itmp(2*n+2)   = itmp(2*n+2) + 1
      itmp(2*n+3)   = itmp(2*n+3) + 1
!
      do j=1,n-1
        ij=ij+1
        Acol(ij)=j
        gc(ij)   =-1.d0
        ij=ij+1
        Acol(ij)=j+1
        gc(ij) = 1.d0
        itmp(j+2*n+3)   = itmp(j+2*n+3) + 2
      enddo
!
      ij=ij+1
      Acol(ij)=1
      gc(ij) = 1.d0
      ij=ij+1
      Acol(ij)=n
      gc(ij) =-1.d0
      itmp(3*n+3)   = itmp(3*n+3) + 1
      itmp(3*n+4)   = itmp(3*n+4) + 1
!
      Aptr(1)=1
      do i=2,ni+1
        Aptr(i)=Aptr(i-1)+itmp(i-1)
      enddo
!
      return
      end subroutine SAOi_grads
