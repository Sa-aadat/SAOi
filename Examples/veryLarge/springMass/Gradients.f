! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, albert Groenwold, Stellenbosch
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
      integer          nmass
      double precision x(n), gf(n), gc(*)
      double precision al, ali, pe1, pe2, xi, yi, wi, xip1, yip1
      double precision aki, dli, dli1
!
      if (mod(n,2).ne.0) stop ' n is odd '
!
      NMASS=n/2
!
      do i=1,n
        gf(i)=0.d0
      enddo
!
      al=60.
      ali=al/(dble(nmass)+1.)
      pe1=0.
      pe2=0.
      xi=0.
      yi=0.
      do 20 i=1,nmass
         j=i-1
         wi=50.*dble(i)
         pe1=pe1+wi*x(i+nmass)
         gf(i+nmass)=gf(i+nmass)+wi
         xip1=x(i)
         yip1=x(i+nmass)
         aki=500.+200.*((dble(nmass)/3.-dble(i))**2)
         dli1=dsqrt((ali+xip1-xi)**2+(yip1-yi)**2)
         dli=dli1-ali
         gf(i)=gf(i)+aki*dli*(ali+xip1-xi)/dli1
         gf(i+nmass)=gf(i+nmass)+aki*dli*(yip1-yi)/dli1
         if(j.gt.0) then
           gf(j)=gf(j)-aki*dli*(ali+xip1-xi)/dli1
           gf(j+nmass)=gf(j+nmass)-aki*dli*(yip1-yi)/dli1
         endif
         pe2=pe2+0.5*aki*(dli**2)
         xi=xip1
         yi=yip1
20    continue

!  last spring
      aki=500.+200.*((dble(nmass)/3.-dble(nmass+1))**2)
      dli1=dsqrt((ali-xi)**2+yi**2)
      dli=dli1-ali
      pe2=pe2+0.5*aki*(dli**2)
      gf(nmass)=gf(nmass)-aki*dli*(ali-xi)/dli1
      gf(n)=gf(n)+aki*dli*yi/dli1
!
      return
      end subroutine SAOi_grads
