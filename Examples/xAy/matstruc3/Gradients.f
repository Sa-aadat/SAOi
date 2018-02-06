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
!  derivatives gc(k) of the (in)equality constraint functions c        !
!  w.r.t. the variables x(i)                                           !
!                                                                      !
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!                                                                      !
!  The default storage scheme is the dense storage scheme. In          !
!      this scheme, it is only required to specify the vector          !
!      gf and the matrix gc. The dimensions of gf and gc are           !
!      gf(n) and gc(m,n) respectively. For details, please             !
!      see the users manual                                            !
!                                                                      !
!  The default storage scheme if the sparse implementation is          !
!      selected is the compressed sparse row (CSR) storage scheme.     !
!      The dimensions of Acol and Aptr in the CSR storage sheme        !
!      are Acol(nnz) and Aptr(m+1) respectively. For details,          !
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
      implicit          none
      include           'ctrl.h'
      logical           eqn(*), lin(*)
      integer           i, j, n, m, ni, ne, nnz, Acol(*), Aptr(*)
      integer           allocatestatus, l, k, p, d, k1, a
      double precision  x(*), gf(*), gc(m,n)
      double precision, dimension(:),    allocatable :: xs, mode, temp
      double precision, dimension(:, :), allocatable :: toeplitz, xx
!
      k = iuser(1)  !  problem specific - number of requested eigen values
      p = iuser(2)  !  problem specific - length of each eigen vector column  
!
      allocate ( xs(p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate xs"
!
      allocate ( temp(p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate temp"
!
      allocate ( mode(p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate mode"
!
      allocate ( xx(p, k), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate xx"
!
      allocate ( toeplitz(p, p), stat = allocatestatus)
      if (allocatestatus /= 0) 
     &          STOP " not enough memory to allocate toeplitz"
!
      xs = (/ 4.d0, 1.d0, 0.d0, 0.d0, 0.d0 /)
!
      do i = 1,p
        toeplitz(i,i:p) = xs(1:(p+1-i))
        toeplitz(i:p,i) = xs(1:(p+1-i))
      end do  
!      
!       do i=1,p 
!         write(*,*) i,(toeplitz(i,j),j=1,p)
!       enddo
!       stop
!
!  Construct k eigenvectors in xx from the single design variable array x
!
      do l = 1,k
        xx(:,l) = x((l-1)*p+1:l*p)
      enddo
!
!  Construct objective gradients
!
      i = 0
      do l = 1,k
        mode(:) = xx(:,l) 
        temp = matmul(toeplitz,mode)           ! toeplitz * mode 
        do j = 1,p
          i = i + 1
          gf(i) = -2.d0*temp(j)*dble(k+100*l-100)
        enddo  
      enddo 
!
!  Construct constraint gradients
!
      gc = 0.d0
!
!  Enforce off-diagonal terms = 0
!       
      i=0 
      do l=1,k-1
        mode(:) = xx(:,l)
        do d = l+1,k
          i = i+1
          temp(:) = xx(:,d)
          do j = 1,p
            gc(i,(l-1)*p+j) = temp(j)
            gc(i,(d-1)*p+j) = mode(j)
          enddo
        enddo  
      enddo  
!
!  Enforce diagonal terms - 1 = 0
!
      k1=k*(k-1)/2
      do i=1,k
        mode(:) = xx(:,i)      
        do j=1,p
          gc(i+k1,(i-1)*p+j) = 2.d0*mode(j)
        enddo
      enddo      
!            
!  Deallocate some stuff here ...
!
      return
      end subroutine SAOi_grads
