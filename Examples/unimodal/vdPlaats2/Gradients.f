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
      integer          i, j, k, n, m, ni, ne, nseg, i1, j1
      integer          nnz, Acol(*), Aptr(*)
      double precision x(n), gf(n), gc(*)
      double precision pload, young, xlen, xleni
      double precision ybar, a1, b1, xinert, xmomi, sigi
      double precision xidb, xidh, sbar,sigmabar
!
      do k=1,n*ni
        gc(k)=0.d0
      enddo
      nseg=n/2
      pload=50000.d0
      young=2.d7
      xlen=500.d0
      xleni=xlen/dble(nseg)
      sigmabar=14000.d0
      do i=1,nseg
        xmomi=pload*(xlen+xleni-dble(i)*xleni)
!
        gf(2*i-1) =          x(2*i)*xleni
        gf(2*i  ) = x(2*i-1)*       xleni
!
        k = (2*i-1-1)*ni + 2*i-1
        gc(k)= -6.d0*xmomi/(sigmabar*x(2*i-1)**2*x(2*i)**2)
!
        k = (2*i-1)*ni + 2*i-1
        gc(k)= -12.d0*xmomi/(sigmabar*x(2*i-1)*x(2*i)**3)
!
        k = (2*i-1-1)*ni + 2*i
        gc(k) = -20.d0
!
        k = (2*i-1)*ni + 2*i
        gc(k) =  1.d0
      enddo
!
      return
      end subroutine SAOi_grads
