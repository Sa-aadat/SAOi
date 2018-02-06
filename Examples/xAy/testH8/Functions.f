! SAOi:
! SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
! SAOi: objective and constraint functions
! SAOi:

      subroutine SAOi_funcs (n, m, ni, ne, x, f, c, iuser, luser, cuser,
     &                       ruser, eqn, lin, ictrl, lctrl, rctrl,
     &                       cctrl)
!----------------------------------------------------------------------!
!                                                                      !
!  Compute the objective function f and the inequality constraint      !
!  functions c(j), j=1,m                                               !
!                                                                      !
!  Please see the users manual for type declarations and comments      !
!                                                                      !
!----------------------------------------------------------------------!
      implicit          none
      include           'ctrl.h'
      logical           eqn(*), lin(*)
      integer           i, j, n, m, ni, ne
      integer           allocatestatus, l, k, p, d, k1
      integer           iscale1, iscale2
      double precision  f, x(n), c(m)
      double precision, dimension(:),    allocatable :: xs, mode, temp
      double precision, dimension(:, :), allocatable :: toeplitz, xx
!
      k       = iuser(1)  !  problem specific - number of requested eigen values
      p       = iuser(2)  !  problem specific - length of each eigen vector column  
      iscale1 = iuser(3)  !  problem specific 
      iscale2 = iuser(4)  !  problem specific 
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
      xs    = 0.d0
      xs(1) = 4.d0
      xs(2) = 1.d0
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
!  Construct objective
!
      f = 0.d0
!
      do l = 1,k
        mode(:) = xx(:,l) 
        !temp = matmul(toeplitz,mode)           ! toeplitz * mode 
        temp = matmul(toeplitz,xx(:,l))
        f = f - dot_product(temp,mode)*dble(k+iscale1*l-iscale2)  ! maximize weighted objectve: mode' * toeplitz * mode * weight 
      enddo 
!
!  Construct constraints
!
!  Enforce off-diagonal terms = 0
!       
      i=0 
      do l=1,k-1
        mode(:) = xx(:,l)
        do d = l+1,k
          i = i+1
          temp(:) = xx(:,d)
          !c(i) = dot_product(mode,temp)       !   either this one line...
          c(i) = 0.d0                          ! }
          do j = 1,p                           ! } or all of this 
            c(i) = c(i) + mode(j)*temp(j)      ! }
          enddo                                ! }
        enddo  
      enddo  
!
!  Enforce diagonal terms - 1 = 0
!
      k1=k*(k-1)/2
      do i=1,k
        mode(:) = xx(:,i)      
        !c(i+k1) = dot_product(mode,mode) - 1.d0   !   either this one line...
        c(i+k1) = -1.d0                            ! }
        do j = 1,p                                 ! } or all of this 
          c(i+k1) = c(i+k1) + mode(j)*mode(j)      ! }
        enddo                                      ! }
      enddo      
!            
!  Deallocate some stuff here ...
!

!
!  Write in mode form if desired - move to subroutine Specialized.f 
!
!      write(*,*) ' '
!      do i=1,p
!        write(*,1000) (xx(i,l),l=1,k)
!      enddo
! 1000 format(1000f10.4)
!      write(*,*) ' '
!
!
!
      return
      end subroutine SAOi_funcs
!----------------------------------------------------------------------!
