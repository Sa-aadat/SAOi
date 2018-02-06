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
      implicit         none
      include          'ctrl.h'
      logical          eqn(*), lin(*)
      integer          i, j, n, m, ni, ne, nnz, acol(*), aptr(*)
      double precision x(n), gf(n), gc(m,*) 
      double precision v1, v2, c, d, v3, v4, v5, v6, v7, v8, v9, v10
      double precision v13, v14, v15, y1, y2, y3, y4, y5, y6, v11, v12
!
      gf(1)=3.d+3+3.d+3*x(1)**2 
      gf(2)=2.d+3+2.000001d+3*x(2)**2  
!
      do i=3,9,1
        gf(i) = 0.d0
      enddo
!
      y1=dsin(x(8))     
      y2=dcos(x(8))     
      y3=dsin(x(9))     
      y4=dcos(x(9))     
      y5=dsin(x(8)-x(9))
      y6=dcos(x(8)-x(9))
!
      v1=48.4d0/50.176d0
      c=v1*dsin(0.25d0) 
      d=v1*dcos(0.25d0) 
!
      do i=1,6       
        do j=1,9       
          gc(i,j)=0.d0  
        enddo
      enddo
!
      gc(1,1)=-1.d0     
      gc(2,2)=-1.d0     
      gc(4,3)=-1.d0     
      gc(5,4)=-1.d0     
!
      v1=-d*y1-c*y2     
      v2=-d*y3-c*y4     
      gc(1,5)=4.d0*c*x(5)+x(6)*v1+x(7)*v2       
      gc(1,6)=x(5)*v1   
      gc(1,7)=x(5)*v2   
      gc(1,8)=x(5)*x(6)*(-d*y2+c*y1)    
      gc(1,9)=x(5)*x(7)*(-d*y4+c*y3)    
!
      v2=d*y1-c*y2      
      v3=d*y6+c*y5      
      v4=d*y5-c*y6      
      gc(2,5)=x(6)*v2   
      gc(2,6)=4.d0*c*x(6)+x(5)*v2+x(7)*v4       
      gc(2,7)=x(6)*v4   
      gc(2,8)=x(5)*x(6)*(d*y2+c*y1)+x(6)*x(7)*v3
      gc(2,9)=-x(6)*x(7)*v3     
!
      v5=d*y3-c*y4      
      v6=-d*y5-c*y6     
      v7=-d*y6+c*y5     
      gc(3,5)=x(7)*v5   
      gc(3,6)=x(7)*v6   
      gc(3,7)=4.d0*c*x(7)+x(5)*v5+x(6)*v6       
      gc(3,8)=x(6)*x(7)*v7      
      gc(3,9)=x(5)*x(7)*(d*y4+c*y3)-x(6)*x(7)*v7
!    
      v8=-c*y1+d*y2     
      v9=-c*y3+d*y4     
      gc(4,5)=4.d0*d*x(5)-x(6)*v8-x(7)*v9       
      gc(4,6)=-x(5)*v8  
      gc(4,7)=-x(5)*v9  
      gc(4,8)=x(5)*x(6)*(c*y2+d*y1)     
      gc(4,9)=x(5)*x(7)*(c*y4+d*y3)     
!      
      v10=c*y1+d*y2     
      v11=c*y5+d*y6     
      v12=(c*y6-d*y5)*x(6)      
      gc(5,5)=-x(6)*v10 
      gc(5,6)=4.d0*d*x(6)-x(5)*v10-x(7)*v11     
      gc(5,7)=-x(6)*v11 
      gc(5,8)=-x(5)*x(6)*(c*y2-d*y1)-x(7)*v12   
      gc(5,9)=x(7)*v12  
!     
      v13=c*y3+d*y4     
      v14=-c*y5+d*y6    
      v15=(c*y6+d*y5)*x(6)*x(7) 
      gc(6,5)=-x(7)*v13 
      gc(6,6)=-x(7)*v14 
      gc(6,7)=4.d0*d*x(7)-x(5)*v13-x(6)*v14     
      gc(6,8)=v15       
      gc(6,9)=-x(5)*x(7)*(c*y4-d*y3)-v15
!
      return
      end subroutine SAOi_grads
