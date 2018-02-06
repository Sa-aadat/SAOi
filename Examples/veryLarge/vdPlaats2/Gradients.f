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
      integer          i, j, k, n, m, ni, ne, nseg, ijnum
      integer          nnz, Acol(*), Aptr(*)
      double precision x(n), gf(n), gc(*)
      double precision pload, young, xlen, xleni
      double precision a1, b1, xinert, xmomi, sigi
      double precision xidb, xidh, sbar, gcji
      double precision aterm, bterm, sigmabar, xinertdb
      double precision xinertdh
!
      nseg     = n/2
      pload    = 50000.d0
      young    = 2.d7
      xlen     = 500.d0
      xleni    = xlen/dble(nseg)
      sigmabar = 14000.d0
!
      Aptr(1)=1
      ijnum=0
      do i=1,nseg
!
        gf(2*i-1) =          x(2*i)*xleni
        gf(2*i  ) = x(2*i-1)*       xleni
!
        Aptr(i+1)=Aptr(i)+2
        xmomi=pload*(xlen+xleni-dble(i)*xleni)
!
        ijnum=ijnum+1
        Acol(ijnum)=2*i-1
        gc(ijnum)= -6.d0*xmomi/
     &                   (sigmabar*x(2*i-1)**2*x(2*i  )**2)
!
        ijnum=ijnum+1
        Acol(ijnum)=2*i
        gc(ijnum)= -12.d0*xmomi/
     &                   (sigmabar*x(2*i-1)*x(2*i  )**3)
!
      enddo
!
      do i=1,nseg
        Aptr(nseg+i+1)=Aptr(nseg+i)+2
!
        ijnum=ijnum+1
        Acol(ijnum)=2*i-1
        gc(ijnum) =-20.d0
!
        ijnum=ijnum+1
        Acol(ijnum)=2*i
        gc(ijnum) =  1.d0
      enddo
!
      return
      end subroutine SAOi_grads
