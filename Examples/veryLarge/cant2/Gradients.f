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
      double precision Length, P, E, li, ymax1, t, Ro1, dinerti
      double precision sumlj, temp1, temp2, ydi, y1i
!
      Length = 500.d0
      P      = 40.d0
      E      = 2.0d8
      li     = Length/dble(n)
      ymax1  = 0.5d0
      t      = 0.2d0
      Ro1    = 0.00078d0
!
      do i = 1,n
        gf(i) = 4.d0*t*li*Ro1
      enddo
!
      do i = 1,n
        dinerti = -2.d0*(x(i)**4)*t/9.d0
        sumlj = li*i
        temp1 = (Length+li/2.d0-sumlj)
        temp2 = (Length+2.d0*li/3.d0-sumlj)
!
        ydi   = P*li/(E*dinerti)*temp1
        y1i   = P*li*li/(2.d0*E*dinerti)*temp2
!
        gc(i) = y1i/ymax1 + (n-i)*ydi*li/ymax1
      enddo
!
      return
      end subroutine SAOi_grads
