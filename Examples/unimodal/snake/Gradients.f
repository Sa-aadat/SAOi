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
      integer          i, j, n, m, ni, ne, l2030, nnz, Acol(*), Aptr(*)
      integer          ix, iy, iz, i1, i2, i3, i4, k
      integer          i1jx, i1jy, i1jz, i2jx, i2jy, i2jz
      integer          i3jx, i3jy, i3jz, i4jx, i4jy, i4jz
      integer          ij, jx, jy, jz, m1, mjx, mjy
      double precision x(n), gf(n), gc(*)
      double precision pi, zzk, pii2, zzi, alfaj, sumkva, xy1, zxy, aj
      double precision zzj
!
      pi = 4.0d0*datan(1.0d0)
!
      l2030=10
      k=l2030
      zzk=dble(k)
      pii2=pi/2.d0
      do j=1,k
        zzj=dble(j)
        aj=pii2*zzj/zzk - 2.d0*pii2/3.d0
        jx=j
        jy=k+j
        jz=2*k+j
        gf(jx)=dcos(aj)
        gf(jy)=dsin(aj)
        gf(jz)=-0.1d0
      enddo
!
      do ij=1,ni*n
        gc(ij)=0.d0
      enddo
!
      m1=ni
      do i=1,k
        jx=i
        jy=k+i
        jz=2*k+i
        i1=i
        i2=k+i
        i3=2*k+i
        i4=3*k+i
        xy1 = x(jx)**2 + x(jy)**2 - 1.d0
        zxy = x(jz) - 2.d0*x(jx)*x(jy)
        i1jx=m1*(jx-1)+i1
        i1jy=m1*(jy-1)+i1
        i1jz=m1*(jz-1)+i1
        i2jx=m1*(jx-1)+i2
        i2jy=m1*(jy-1)+i2
        i2jz=m1*(jz-1)+i2
        i3jx=m1*(jx-1)+i3
        i3jy=m1*(jy-1)+i3
        i3jz=m1*(jz-1)+i3
        i4jx=m1*(jx-1)+i4
        i4jy=m1*(jy-1)+i4
        i4jz=m1*(jz-1)+i4
        mjx=m1*jx
        mjy=m1*jy
        gc(i1jx)= 2.d0*x(jx)*(10.d0+70.d0*(10.d0*xy1)**6)
        gc(i1jy)= 2.d0*x(jy)*(10.d0+70.d0*(10.d0*xy1)**6)
        gc(i1jz)= 0.d0
        gc(i2jx)=-2.d0*x(jy)*(10.d0+70.d0*(10.d0*zxy)**6)
        gc(i2jy)=-2.d0*x(jx)*(10.d0+70.d0*(10.d0*zxy)**6)
        gc(i2jz)= 10.d0+70.d0*(10.d0*zxy)**6
        gc(i3jx)=-gc(i1jx)
        gc(i3jy)=-gc(i1jy)
        gc(i3jz)=-gc(i1jz)
        gc(i4jx)=-gc(i2jx)
        gc(i4jy)=-gc(i2jy)
        gc(i4jz)=-gc(i2jz)
        gc(mjx)= 2.d0*x(jx)
        gc(mjy)= 2.d0*x(jy)
      enddo
!
      return
      end subroutine SAOi_grads
