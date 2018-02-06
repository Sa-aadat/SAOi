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
      integer          i, j, n, m, ni, ne, k, nnz, Acol(*), Aptr(*)
      double precision x(n), gf(n), gc(*)
      double precision rads, pi
!
      pi = 4.0d0*datan(1.0d0)
!
      rads=pi/180.d0
!
      gf( 1)=-x( 2)*dsin(x(12)*rads)/2.d0
      gf(11)=-x(10)*dsin(x(21)*rads)/2.d0
      do i=2,10
        gf(i)=-x(i-1)*dsin(x(i+10)*rads)/2.d0
     &        -x(i+1)*dsin(x(i+11)*rads)/2.d0
      enddo
!
      do i=12,21
        gf(i)=-x(i-11)*x(i-10)*dcos(x(i)*rads)*rads/2.d0
      enddo
!
      i=1
      gc(i)=0.5d0*(x(i)**2 + x(i+1)**2
     &       -2.d0*x(i)*x(i+1)*dcos(x(i+11)*rads))**(-0.5d0)
     &       *(2.d0*x(i)-2.d0*x(i+1)*dcos(x(i+11)*rads))
     &       +1.d0
!
      i=11
      gc(i)=0.5d0*(x(i-1)**2 + x(i)**2
     &       -2.d0*x(i-1)*x(i)*dcos(x(i+10)*rads))**(-0.5d0)
     &       *(2.d0*x(i)-2.d0*x(i-1)*dcos(x(i+10)*rads))
     &       +1.d0
!
      do i=2,10
        gc(i)=0.5d0*(x(i)**2 + x(i+1)**2
     &         -2.d0*x(i)*x(i+1)*dcos(x(i+11)*rads))**(-0.5d0)
     &         *(2.d0*x(i)-2.d0*x(i+1)*dcos(x(i+11)*rads))
     &         +0.5d0*(x(i-1)**2 + x(i)**2
     &         -2.d0*x(i-1)*x(i)*dcos(x(i+10)*rads))**(-0.5d0)
     &         *(2.d0*x(i)-2.d0*x(i-1)*dcos(x(i+10)*rads))
      enddo
!
      do i=12,21
        gc(i)=0.5d0*(x(i-11)**2 + x(i-10)**2
     &         -2.d0*x(i-11)*x(i-10)*dcos(x(i)*rads))**(-0.5d0)
     &         *2.d0*x(i-11)*x(i-10)*dsin(x(i)*rads)*rads
      enddo
!
      return
      end subroutine SAOi_grads
