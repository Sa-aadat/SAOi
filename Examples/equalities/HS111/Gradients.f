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
      integer          i, j, n, m, ni, ne, nnz, acol(*), aptr(*)
      double precision x(n), gf(n), gc(m,*)
      double precision y(10), t
!
      y(1)=-6.089d0     
      y(2)=-17.164d0    
      y(3)=-34.054d0    
      y(4)=-5.914d0     
      y(5)=-24.721d0    
      y(6)=-14.986d0    
      y(7)=-24.1d0      
      y(8)=-10.708d0    
      y(9)=-26.662d0    
      y(10)=-22.179d0   
!
      t=0.d0    
      do i=1,10      
        t=t+dexp(x(i))  
      enddo
!      
      do j=1,10            
        do i=1,3       
          gc(i,j)=0.d0  
        enddo
      enddo
!
      do i=1,10      
        gf(i)=dexp(x(i))*(y(i)+x(i)-dlog(t)) 
      enddo     
!
      gc(1,1) = dexp(x(1))
      gc(1,2) = 2.d0*dexp(x(2))   
      gc(1,3) = 2.d0*dexp(x(3))   
      gc(1,6) = dexp(x(6))
      gc(1,10)= dexp(x(10)) 
! 
      gc(2,4) = dexp(x(4))
      gc(2,5) = 2.d0*dexp(x(5))   
      gc(2,6) = dexp(x(6))
      gc(2,7) = dexp(x(7))
!                
      gc(3,3) = dexp(x(3))
      gc(3,7) = dexp(x(7))
      gc(3,8) = dexp(x(8))
      gc(3,9) = 2.d0*dexp(x(9))   
      gc(3,10)= dexp(x(10))  
!
      return
      end subroutine SAOi_grads
