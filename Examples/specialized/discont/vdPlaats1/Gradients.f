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
      double precision pload, young, xlen, slen(n/2), slenc
      double precision ybar, a1, b1, xinert, xmomi, sigi
      double precision xidb, xidh, sbar, fact
!
      do k=1,ni*n
        gc(k) = 0.d0
      enddo
!
      nseg    = n/2
      pload   = 50000.d0
      young   = 2.d7
      xlen    = 500.d0
      slenc   = xlen/dble(nseg)
      sbar    = 14000.d0
      ybar    = 2.5d0
      fact    = ruser(1)
!
      do i=1,nseg-1,2
        slen(i) = slenc+slenc*rctrl(34)*fact
      enddo
!
      do i=2,nseg,2
        slen(i) = slenc-slenc*rctrl(34)*fact
      enddo
!
      do i=1,nseg
        xmomi=pload*(xlen+slen(i)-dble(i)*slen(i))
        gf(2*i-1)=x(2*i)*slen(i)
        gf(2*i  )=x(2*i-1)*slen(i)
        j1=2*i-1
        i1=2*i-1
        k=(i1-1)*ni+j1
        gc(k)=-6.d0*xmomi/(sbar*x(2*i-1)**2*x(2*i)**2)
        i1=2*i
        k=(i1-1)*ni+j1
        gc(k)=-12.d0*xmomi/(sbar*x(2*i-1)*x(2*i)**3)
        j1=2*i
        i1=2*i-1
        k=(i1-1)*ni+j1
        gc(k)=-20.d0
        i1=2*i
        k=(i1-1)*ni+j1
        gc(k)=1.d0
      enddo
!
      do i=1,nseg
        xidb=-(x(2*i-1)**2*x(2*i)**3)
        xidh=-(x(2*i-1)*x(2*i)**4)/3.d0
        a1=(xlen-dble(i)*slen(i)+2.d0/3.d0*slen(i))
        j1=ni
        i1=2*i-1
        k=(i1-1)*ni+j1
        gc(k)=(12.d0*pload*slen(i)**2)/(2.d0*young*xidb)*a1/ybar
        i1=2*i
        k=(i1-1)*ni+j1
        gc(k)=(12.d0*pload*slen(i)**2)/(2.d0*young*xidh)*a1/ybar
      enddo
!
      do i=1,nseg-1
        xidb=-(x(2*i-1)**2*x(2*i)**3)
        xidh=-(x(2*i-1)*x(2*i)**4)/3.d0
        b1=(xlen+slen(i)/2.d0-dble(i)*slen(i))
        j1=ni
        i1=2*i-1
        k=(i1-1)*ni+j1
        gc(k)=gc(k)+(12.d0*pload*slen(i)**2)/(young*xidb)*
     &                 b1*(nseg-i)/ybar
        i1=2*i
        k=(i1-1)*ni+j1
        gc(k)=gc(k)+(12.d0*pload*slen(i)**2)/(young*xidh)*
     &                 b1*(nseg-i)/ybar
      enddo
!
      return
      end subroutine SAOi_grads
