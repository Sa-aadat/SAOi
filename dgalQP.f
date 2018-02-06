! SAOi:
! SAOi: Wed Jul 14 14:26:25 SAST 2010, Salomon van Huyssteen
! SAOi: and Albert Groenwold, Stellenbosch. Driver for
! SAOi: the Galahad QP solvers (double) to solve sparse QP-like
! SAOi: DQ subproblems.
! SAOi:


      subroutine drive_galQP (n,x,m,nq,f_a,c_a,h_a,xlam,x_h,
     &                        acurvn,bcurvn,ccurvn,
     &                        x_l,x_u,f,c,h,gf,gc,gh,ksubiter,
     &                        xkkt,flam,cstage0,ostring,xlam_h,
     &                        nsx,ictrl,lctrl,rctrl,cctrl,message,yhi,
     &                        outerloop,n2,ierr1,eqn,lin,
     &                        iuser,luser,cuser,ruser,sloop,z)
      implicit none
      include          'ctrl.h'
      character*16     A_string,H_string
      character*90     ostring
      logical          kkt_local,cstage0
      logical          eqn(*),lin(*)
      integer          m,nq,n,n2,n3,ierr,a_ne,h_ne,ksubiter,outerloop
      integer          i,j,k,sdr,relaxed,nsx,message
      integer          A_row(n*m),A_col(n*m),A_ptr(m+1),ierr1
      integer          H_row(n2),H_col(n2),H_ptr(n+1),sloop
      double precision w(n),xlamj,bji,gh(nq,n)
      double precision x(n),xi(n),x_li(n),x_ui(n)
      double precision x_h(n),x_l(n),x_u(n),xlam(m)
      double precision acurvn(n),bcurvn(m,n),gc_h(m,n),gf_h(n)
      double precision ccurvn(nq,n),f,c(m),h(nq),gf(n),gc(m,n)
      double precision f_a,c_a(m),h_a(nq),xkkt,gf_a(n),gc_a(m,n)
      double precision y(m),z(n),xlam_h(m),yhi,htol
      double precision flam,c_l(m),c_u(m),A_val(n*m),g(n)
      double precision H_val(n2)
      include          'ctrl_get.inc'
!
      ierr = 0
      QPform_A = 2 ! 0 = dense / 1 = coordinate / 2 = sparse_by_rows                   ! 2 = default 
      QPform_H = 3 ! 0 = dense / 1 = coordinate / 2 = sparse_by_rows / 3 = diagonal    ! 3 = default 
!
!     QPform_H = 0 will only be retained if subsolver .ne. 25; else the following happens
!
      if (subsolver.eq.25) QPform_H = 3 ! Galahad lsqp solves diagonal QP problems only
!
      if (GetFinDiffGradHessDense) then
        QPform_H = 0
        if (subsolver.eq.25) then 
          write(*,*) ' STOP: diagonal subsolver 25 not',
     &               ' possible with GetFinDiffGradHessDense' 
          stop
        endif  
      endif
!
      htol = 1.d-6
!
      A_row=0
!      
      H_ptr=0
!
      if (QPform_A.eq.0.and.outerloop.eq.1.and.sloop.eq.1) then
        write(*,*) ' '
        write(*,*) 'Warning: QPform_A.eq.0 (dense) in drive_galQP'
        write(*,*) ' '
      endif  
      if (QPform_H.eq.0.and.outerloop.eq.1.and.sloop.eq.1) then
        write(*,*) ' '
        write(*,*) 'Warning: QPform_H.eq.0 (dense) in drive_galQP'
        write(*,*) ' '
      endif  
!      
      call QP_set (n,m,n2,a_ne,h_ne,gc,QPform_A,QPform_H,
     &             A_string,H_string,acurvn,bcurvn,
     &             ictrl,lctrl,rctrl,cctrl,
     &             iuser,luser,cuser,ruser)
!
      call QP_pre (n,m,a_ne,f,c,acurvn,
     &             bcurvn,xlam,gc,x,xi,x_l,x_u,x_li,x_ui,
     &             z,c_l,c_u,w,A_row,A_col,A_val,A_ptr,
     &             A_string,QPform_A,QPform_H,H_string,
     &             H_row,H_col,H_ptr,H_val,h_ne,ictrl,
     &             lctrl,rctrl,cctrl,outerloop,eqn,lin,htol,
     &             iuser,luser,cuser,ruser,n2,xlam_h)

! call the Galahad QP solvers (double)
      if     (subsolver.eq.25) then
        call galahad_lsqp (n,m,a_ne,f_a,c_l,c_u,x_li,
     &                     x_ui,A_row,A_col,A_ptr,xi,gf,
     &                     w,xlam,z,A_val,A_string,ierr,
     &                     ksubiter,f)

      elseif (subsolver.eq.26) then
        call galahad_qpa  (n,m,a_ne,f_a,c_l,c_u,x_li,
     &                     x_ui,A_row,A_col,A_ptr,xi,gf,
     &                     w,xlam,z,A_val,A_string,ierr,
     &                     ksubiter,H_string,H_row,H_col,
     &                     H_ptr,H_val,h_ne,f)

      elseif (subsolver.eq.27) then
        call galahad_qpb  (n,m,a_ne,f_a,c_l,c_u,x_li,
     &                     x_ui,A_row,A_col,A_ptr,xi,gf,
     &                     w,xlam,z,A_val,A_string,ierr,
     &                     ksubiter,H_string,H_row,H_col,
     &                     H_ptr,H_val,h_ne,f)

      elseif (subsolver.eq.28) then
        call galahad_qpc  (n,m,a_ne,f_a,c_l,c_u,x_li,
     &                     x_ui,A_row,A_col,A_ptr,xi,gf,
     &                     w,xlam,z,A_val,A_string,ierr,
     &                     ksubiter,H_string,H_row,H_col,
     &                     H_ptr,H_val,h_ne,f)
      else
!
        stop ' subsolver invalid in drive_galQP'
!
      endif

! update multipliers
      do j=1,m
        xlam(j)= max(-biglam,min(-xlam(j),biglam))
      enddo
      do i=1,n
        z(i)= max(-biglam,min(z(i),biglam))
      enddo

! update x in current point
      do i=1,n
        x(i)=min(x_u(i),max(x_l(i),x(i)+xi(i)))
      enddo
!
      message = ierr
!
      return
      end subroutine drive_galQP
!----------------------------------------------------------------------!
      subroutine QP_pre (n,m,a_ne,f,c,acurvn,
     &                   bcurvn,xlam,gc,x,xi,x_l,x_u,x_li,x_ui,
     &                   z,c_l,c_u,w,A_row,A_col,A_val,A_ptr,
     &                   A_string,QPform_A,QPform_H,H_string,
     &                   H_row,H_col,H_ptr,H_val,h_ne,
     &                   ictrl,lctrl,rctrl,cctrl,outerloop,
     &                   eqn,lin,htol,
     &                   iuser,luser,cuser,ruser,n2,xlam_h)
      implicit none
      include          'ctrl.h'
      integer          n,m,a_ne,h_ne,ierr,i,j,k,l1,k1
      integer          A_row(*),A_col(*),A_ptr(*),q,test,n2
      integer          H_row(*),H_col(*),H_ptr(*),outerloop
      logical          eqn(*),lin(*)
      double precision w(*),xlamj,bji,x(*),xi(*),x_li(*),x_ui(*)
      double precision x_l(*),x_u(*),xlam(m),gc(m,*),z(*),c(*),f
      double precision c_l(*),c_u(*),A_val(*),acurvn(*),bcurvn(m,*)
      double precision H_val(*),ddot,htol,xlam_h(m)

!     double precision l(3), p(3)
!     double precision E, Sy, Fl, d, h 

      character*16     A_string
      character*16     H_string
      include          'ctrl_get.inc'

!  Initialize n side constraint multipliers
!       do i=1,n
!         z(i) = 0.d0
!       enddo

!  m constraint lower bounds
      do j=1,m
        if (eqn(j)) then
          c_l(j) = -c(j)
        else
          c_l(j) = -1.d20
        endif
      enddo

!  m constraint upper bounds
      do j=1,m
        c_u(j) = -c(j)
      enddo

!  Determine approximate dual variables for the first iteration
!       if (outerloop.eq.1.and.eqn(m)) then
!         call drive_lbfgsb24u2 (m,xlam,f,c,eqn,lin,
!       &                        ictrl,lctrl,rctrl,cctrl,
!       &                        iuser,luser,cuser,ruser)
!         xlam_h = xlam
!       endif  
 
!  Construct the diagonal Hessian for SAOi - we work with acurvn, bcurvn, but approx_f and approx_c are assumed to both equal 1 and NOT 3 (nabla), i.e. acurv and bcurv could instead be used.
      if (approx_f.lt.100.and.approx_c.lt.100) then
        do i=1,n
          w(i) = acurvn(i)
          do j=1,m
            bji   = bcurvn(j,i)
            xlamj = xlam(j)
!            
            if (sph_dia_hess.eq.0) then                        ! all approx incl. SQ form 0 - every term positive
              w(i) = max(htol,acurvn(i))
              if (bji.ge.0.d0.and.xlamj.gt.0.d0) then
                w(i) =  w(i) + bji*xlamj
              endif  
!
            elseif (sph_dia_hess.eq.1.or.sph_dia_hess.eq.2) then
              if (xlamj.ne.0.d0) then
                w(i) =  w(i) + bji*xlamj
              endif
!
            else
              stop ' sph_dia_hess not defined in dgalQP.f'
            endif
          enddo
!          
          if (sph_dia_hess.eq.1) w(i) = max(htol,w(i))        ! all approx incl. SQ form 1 - final sum positive
          if (sph_dia_hess.eq.2) w(i) = max(htol,dabs(w(i)))  ! all approx incl. SQ form 2 - with abs operator 
        enddo 
!        
      elseif (approx_f.eq.100.and.approx_c.eq.100) then
        do i=1,n 
          w(i) = acurvn(i)
          do j=1,m
            bji   = bcurvn(j,i)
            xlamj = xlam(j)
            if (xlamj.ne.0.d0) then
              w(i) =  w(i) + bji*xlamj
            endif
          enddo
!          l
          !if (sph_dia_hess.eq.1) w(i) = max(htol,w(i)) 
          !if (sph_dia_hess.eq.2) w(i) = max(htol,dabs(w(i))) 
        enddo
      endif  

! Transform the Hessian to one of four QP forms. 0, 1, 2 are not sensible for diagonal; only given below for illustration
      if (QPform_H.eq.1) then
        H_string = 'COORDINATE'
        h_ne = n
        do i=1,n
          H_val(i) = w(i)
          H_row(i) = i
          H_col(i) = i
        enddo
!
      else if (QPform_H.eq.2) then
        H_string = 'SPARSE_BY_ROWS'
        h_ne = n
        do i=1,n
          H_val(i) = w(i)
          H_col(i) = i
          H_ptr(i) = i
        enddo
        H_ptr(n+1) = n+1
!
      else if (QPform_H.eq.3) then
        H_string = 'DIAGONAL'
!         h_ne = n
!         do i=1,n
!           H_val(i) = w(i)
!         enddo
!
      else ! 0 or anything else
        H_string = 'DENSE'
        h_ne = n2 ! = (n*n+n)/2
!
!
!  Now fill lower triangular part of symmetric Hessian. 
!
!      First example here is for the exact (indefinite) Hessian of 3barSAND and 3barSANDsing. 
!
!      Second example here is for the Hessian of SVD, using finite differences for the Hessian in ruser.
!
!      For other examples, iuser and ruser may (for example) be used.
!
!     
!  Begin first example 
!
!         do j=1,h_ne
!           H_val(j) = 0.d0
!         enddo
! 
!         E  = 3.d6
!         Sy = 3.d4
!         Fl = 1.d4
!         d  = 0.29d0
!         h  = 1.d2
!         p(1) = -100.d0                   ! the positions
!         p(2) =    0.d0                   ! the positions
!         p(3) =  100.d0                   ! the positions
! 
!         do i=1,3
!           l(i) = dsqrt(p(i)**2 + h**2)   ! the lengths
!         enddo
! 
!         H_val(4)=xlam(7)*(E*p(1)**2/l(1)**3)-xlam(8)*(E*h*p(1)/l(1)**3)
!         H_val(5)=-xlam(7)*(E*h*p(1)/l(1)**3)+ xlam(8)*(E*Fl/l(1)**3)
!         H_val(8)=xlam(7)*(E*p(2)**2/l(2)**3)-xlam(8)*(E*h*p(2)/l(2)**3)
!         H_val(9)=-xlam(7)*(E*h*p(2)/l(2)**3)+xlam(8)*(E*Fl/l(2)**3)
!         H_val(11)=xlam(7)*(E*p(3)**2/l(3)**3)-xlam(8)*(E*h*p(3)/l(3)**3)
!         H_val(12)=-xlam(7)*(E*h*p(3)/l(3)**3)+xlam(8)*(E*Fl/l(3)**3)
!
!  End first example 
!  
!
!
!
!
!  Begin second example
!
        k = 0
        k1 = 0
        do i=1,n
          do l1=1,i
            k=k+1
            k1=k1+1
            H_val(k1) = ruser(k)
          enddo
        enddo
!
        do j=1,m 
          k1 = 0
          xlamj = xlam(j)      
          do i=1,n
            do l1=1,i
              k=k+1
              k1=k1+1
              H_val(k1) = H_val(k1) + xlamj*ruser(k)
            enddo
          enddo
        enddo
!
!  End second example 
!  
      endif      
      
! Transform the Jacobian to one of three QP forms
      if (QPform_A.eq.1) then
        A_string = 'COORDINATE'
        call co_ordinate_A(n,m,gc,A_row,A_col,A_val,a_ne)
!
      else if (QPform_A.eq.2) then
        A_string = 'SPARSE_BY_ROWS'
        call sparse_by_rows_A(n,m,gc,A_ptr,A_col,A_val,a_ne)
!
      else ! 0 or anything else
        A_string = 'DENSE'
        do j=1,m
          do i=1,n
            A_val(n*(j-1)+i) = gc(j,i)
          enddo
        enddo
        a_ne = n*m
      endif

! Initialise move limits
      do i=1,n
        x_li(i)=x_l(i)-x(i)
        x_ui(i)=x_u(i)-x(i)
      enddo

! create copy of x
      call dcopy(n,x,1,xi,1)

      return
      end subroutine QP_pre
!----------------------------------------------------------------------!
      subroutine co_ordinate_A(n,m,gc,A_row,A_col,A_val,a_ne)
      implicit none
      integer          i,j,n,m,a_ne
      integer          A_row(n*m),A_col(n*m)
      double precision A_val(n*m),gc(m,n)
!
      a_ne = 0
      do j=1,m
        do i=1,n
          if (gc(j,i).ne.0.d0) then
             a_ne = a_ne+1
             A_val(a_ne) = gc(j,i)
             A_row(a_ne) = j
             A_col(a_ne) = i
          endif
        enddo
      enddo
!
      return
      end subroutine co_ordinate_A
!----------------------------------------------------------------------!
      subroutine sparse_by_rows_A(n,m,gc,A_ptr,A_col,A_val,a_ne)
      implicit none
      integer          i,j,n,m,a_ne,ptr
      integer          A_ptr(m+1),A_col(n*m)
      double precision A_val(n*m),gc(m,n)
!
!  ne  : number of entries
!  ptr : first entry of each row in A_val
!
      a_ne = 0
      do j=1,m   ! loop through rows
        A_ptr(j) = a_ne + 1
        do i=1,n   ! loop through columns
          if (gc(j,i).ne.0.d0) then
             a_ne = a_ne+1
             A_val(a_ne) = gc(j,i)
             A_col(a_ne) = i
          endif
        enddo
      enddo

! A%ptr(m+1) holds the total number of entries plus one
      A_ptr(m+1) = a_ne+1
      return
      end subroutine sparse_by_rows_A
!----------------------------------------------------------------------!
      subroutine QP_set (n,m,n2,a_ne,h_ne,gc,QPform_A,QPform_H,
     &                   A_string,H_string,acurvn,bcurvn,
     &                   ictrl,lctrl,rctrl,cctrl,
     &                   iuser,luser,cuser,ruser)
      implicit none
      include          'ctrl.h'
      character*16     A_string,H_string
      integer          n,m,n2,a_ne,h_ne
      double precision gc(m,n),acurvn(n),bcurvn(n,m) 
      include          'ctrl_get.inc'
!
      if (QPform_A.eq.0) then      ! 0 = dense 
        A_string = 'DENSE'
        a_ne = n*m
        
      elseif (QPform_A.eq.1) then  ! 1 = coordinate 
        A_string = 'COORDINATE' 
        
      elseif (QPform_A.eq.2) then  ! 2 = sparse by rows
        A_string = 'SPARSE_BY_ROWS'
        
      else
        stop ' QPform_A invalid in QP_set' 
        
      endif
!      
      if (QPform_H.eq.0) then      ! 0 = dense 
        H_string = 'DENSE'
        h_ne = n2
      
      elseif (QPform_H.eq.1) then  ! 1 = coordinate 
        H_string = 'COORDINATE' 
      
      elseif (QPform_H.eq.2) then  ! 2 = sparse by rows
        H_string = 'SPARSE_BY_ROWS'
      
      elseif (QPform_H.eq.3) then  ! 3 = diagonal
        H_string = 'DIAGONAL'

      else
        stop ' QPform_A invalid in QP_set' 
        
      endif
!
      return
      end subroutine QP_set
!----------------------------------------------------------------------!


