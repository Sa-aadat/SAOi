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
      double precision x(n), gf(n), gc(m,n) 
      double precision l(3), s(3), a(3), p(3), u, v, k11, k12, k22
      double precision E, Sy, Fl, d, h 
!
      E  = 3.d6
      Sy = 3.d4
      Fl = 1.d4
      d  = 0.29d0
      h  = 1.d2
      p(1) = -100.d0                   ! the positions
      p(2) =    0.d0                   ! the positions
      p(3) =  100.d0                   ! the positions
!
      do i=1,3
        a(i) = x(i)                    ! the areas
        l(i) = dsqrt(p(i)**2 + h**2)   ! the lengths
      enddo
!
      u = x(4)
      v = x(5)
! 
      k11 = 0.d0
      k12 = 0.d0
      k22 = 0.d0
      do i=1,3
        k11 = k11 + a(i)*p(i)**2/l(i)**3
        k12 = k12 - h*a(i)*p(i)/l(i)**3
        k22 = k22 + Fl*a(i)/l(i)**3
      enddo
      k11 = k11*E
      k12 = k12*E
      k22 = k22*E
!
      do i=1,3
        s(i) = E*(-u*p(i) + h*v)/l(i)**2
      enddo
!
!      f = d*(a(1)*l(1) + a(2)*l(2) + a(3)*l(3)) ! objective
!
!  add gradients for Marlize - gf 
!
      gf(1) = d*l(1) ! objective w.r.t. design variables
      gf(2) = d*l(2)
      gf(3) = d*l(3)
      gf(4) = 0.d0   ! objective w.r.t. state variables
      gf(5) = 0.d0
!      
!  end add gradients for Marlize - gf     
!
!      do i=1,3
!        c(i)   =  s(i)/Sy - 1.d0     ! tension constraints
!        c(i+3) = -s(i)/Sy - 1.d0     ! compression constraints
!      enddo
!      
!      c(7) = k11*u + k12*v - Fl  ! equilibrium constraint in x_1 
!      c(8) = k12*u + k22*v       ! equilibrium constraint in x_2 
!
!  add gradients for Marlize - gc  (primal variables = 1-3, state variables = 4-5) 
!
      do i=1,n
        do j=1,m
          gc(j,i) = 0.d0
        enddo
      enddo

      gc(1,4) = (E*(-p(1))/l(1)**2)/Sy ! tension constraints w.r.t. state vars  
      gc(1,5) = (E*(h)/l(1)**2)/Sy     
      gc(2,4) = (E*(-p(2))/l(2)**2)/Sy      
      gc(2,5) = (E*(h)/l(2)**2)/Sy     
      gc(3,4) = (E*(-p(3))/l(3)**2)/Sy        
      gc(3,5) = (E*(h)/l(3)**2)/Sy
      
      gc(4,4) = -gc(1,4) ! compression constraints w.r.t. state vars 
      gc(4,5) = -gc(1,5)
      gc(5,4) = -gc(2,4)
      gc(5,5) = -gc(2,5)
      gc(6,4) = -gc(3,4)
      gc(6,5) = -gc(3,5)
      
      gc(7,1) = (p(1)**2/l(1)**3*u - h*p(1)/l(1)**3*v)*E !equilibrium constraints w.r.t. design vars
      gc(7,2) = (p(2)**2/l(2)**3*u - h*p(2)/l(2)**3*v)*E
      gc(7,3) = (p(3)**2/l(3)**3*u - h*p(3)/l(3)**3*v)*E
      gc(8,1) = (- h*p(1)/l(1)**3*u + Fl/l(1)**3*v)*E
      gc(8,2) = (- h*p(2)/l(2)**3*u + Fl/l(2)**3*v)*E
      gc(8,3) = (- h*p(3)/l(3)**3*u + Fl/l(3)**3*v)*E
      
      gc(7,4) = k11 !equilibrium constraints w.r.t. state vars 
      gc(7,5) = k12  
      gc(8,4) = k12  
      gc(8,5) = k22      
!      
!  end add gradients for Marlize - gc
!
      return
      end subroutine SAOi_grads
